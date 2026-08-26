const _SPEEDY_TERRARIUM_EXTENSION =
    Base.get_extension(SpeedyWeather, :SpeedyWeatherTerrariumExt)

isnothing(_SPEEDY_TERRARIUM_EXTENSION) && error(
    "SpeedyWeatherTerrariumExt was not loaded after importing SpeedyWeather and Terrarium",
)

const _ReadyESMTerrariumLand = _SPEEDY_TERRARIUM_EXTENSION.TerrariumLand
const _ReadyESMAbstractTerrariumLandModel =
    _SPEEDY_TERRARIUM_EXTENSION.AbstractTerrariumLandModel
const _ReadyESMTerrariumVars = _SPEEDY_TERRARIUM_EXTENSION.TerrariumVars

# SpeedyWeather's current Terrarium extension predates Terrarium 0.1.6 and
# expects `Terrarium.initialize(model)` to return StateVariables. In 0.1.6 that
# public function returns a complete ModelIntegrator. Allocate only the state
# here; SpeedyWeather.initialize! below constructs and initializes the
# integrator after its atmospheric clock has been populated.
function SpeedyWeather.allocate(
    ::SpeedyWeather.AbstractVariable{_ReadyESMTerrariumVars},
    model::SpeedyWeather.PrimitiveWetModel,
)
    land = model.land
    return Terrarium.StateVariables(
        land.model;
        clock = Terrarium.Clock(time = SpeedyWeather.DEFAULT_DATE),
        boundary_conditions = land.boundary_conditions,
        input_variables = land.input_variables,
        fields = land.fields,
    )
end

# SpeedyWeather's Variables tree can contain Terrarium's complete land state,
# but the two packages deliberately use different architecture type systems.
# Teach SpeedyWeather's host-transfer recursion how to cross that integration
# boundary. Terrarium's Adapt implementation reconstructs the state while
# moving every nested Oceananigans field and array to ordinary host storage.
function SpeedyWeather.Architectures.on_architecture(
    ::SpeedyWeather.Architectures.CPU,
    state::Terrarium.StateVariables,
)
    return Adapt.adapt(Array, state)
end

# Adapt intentionally reduces an Oceananigans Field to its full data array.
# That is desirable for a portable checkpoint (it omits grids, boundary
# functions, and model construction state), but it means SpeedyWeather's
# same-concrete-type Variables copier cannot restore the result by itself.
# Restore every Terrarium buffer, including halos, caches, and nested
# namespaces, into the freshly reconstructed land state.
const _ReadyESMOffsetArray = Oceananigans.Fields.OffsetArray

function _copy_terrarium_restart_array!(
    destination::_ReadyESMOffsetArray,
    source::_ReadyESMOffsetArray,
    label,
)
    axes(destination) == axes(source) || throw(DimensionMismatch(
        "$label checkpoint axes $(axes(source)) do not match destination " *
        "$(axes(destination))",
    ))
    # OffsetArrays' generic copyto! iterates logical indices. With a CuArray
    # parent that becomes forbidden scalar GPU indexing. The offsets are
    # already validated above, so transfer the identical contiguous parents
    # through CUDA's supported host-to-device copy path.
    copyto!(parent(destination), parent(source))
    return destination
end

function _copy_terrarium_restart_array!(
    destination::AbstractArray,
    source::AbstractArray,
    label,
)
    axes(destination) == axes(source) || throw(DimensionMismatch(
        "$label checkpoint axes $(axes(source)) do not match destination " *
        "$(axes(destination))",
    ))
    copyto!(destination, source)
    return destination
end

function _copy_terrarium_restart_entry!(
    destination::Oceananigans.Fields.Field,
    source::AbstractArray,
    label,
)
    source_data = source isa Oceananigans.Fields.Field ? source.data : source
    size(destination.data) == size(source_data) || throw(DimensionMismatch(
        "$label checkpoint size $(size(source_data)) does not match destination " *
        "$(size(destination.data))",
    ))
    _copy_terrarium_restart_array!(destination.data, source_data, label)
    return destination
end

function _copy_terrarium_restart_entry!(
    destination::AbstractArray,
    source::AbstractArray,
    label,
)
    source_data = source isa Oceananigans.Fields.Field ? source.data : source
    size(destination) == size(source_data) || throw(DimensionMismatch(
        "$label checkpoint size $(size(source_data)) does not match destination " *
        "$(size(destination))",
    ))
    _copy_terrarium_restart_array!(destination, source_data, label)
    return destination
end

function _copy_terrarium_restart_entry!(
    destination::Base.RefValue,
    source::Base.RefValue,
    label,
)
    destination[] = source[]
    return destination
end

function _copy_terrarium_restart_entry!(
    destination::NamedTuple,
    source::NamedTuple,
    label,
)
    keys(destination) == keys(source) || error(
        "$label checkpoint keys $(keys(source)) do not match destination " *
        "$(keys(destination))",
    )
    for name in keys(destination)
        _copy_terrarium_restart_entry!(
            getproperty(destination, name),
            getproperty(source, name),
            "$label.$name",
        )
    end
    return destination
end


function _copy_terrarium_restart_entry!(
    destination::Terrarium.StateVariables,
    source::Terrarium.StateVariables,
    label,
)
    for group in (
        :prognostic,
        :tendencies,
        :auxiliary,
        :inputs,
        :namespaces,
        :timestepper_cache,
    )
        _copy_terrarium_restart_entry!(
            getproperty(destination, group),
            getproperty(source, group),
            "$label.$group",
        )
    end
    Oceananigans.restore_prognostic_state!(destination.clock, source.clock)
    return destination
end

function _copy_terrarium_restart_entry!(destination, source, label)
    isequal(destination, source) || error(
        "$label contains unsupported non-array restart state of types " *
        "$(typeof(destination)) and $(typeof(source))",
    )
    return destination
end

function SpeedyWeather._copy_entry!(
    destination::Terrarium.StateVariables,
    source::Terrarium.StateVariables,
)
    return _copy_terrarium_restart_entry!(destination, source, "Terrarium")
end

KernelAbstractions.@kernel function _accumulate_land_water_components_kernel!(
    cumulative_precipitation,
    cumulative_evapotranspiration,
    cumulative_runoff,
    precipitation,
    ground_evaporation,
    transpiration,
    runoff,
    scale,
)
    column = @index(Global, Linear)
    cumulative_precipitation[column] += scale * precipitation[column]
    cumulative_evapotranspiration[column] += scale * (
        ground_evaporation[column] + transpiration[column]
    )
    cumulative_runoff[column] += scale * runoff[column]
end

KernelAbstractions.@kernel function _accumulate_land_water_latent_kernel!(
    cumulative_precipitation,
    cumulative_evapotranspiration,
    cumulative_runoff,
    precipitation,
    latent_heat_flux,
    runoff,
    scale,
    inverse_latent_heat_density,
)
    column = @index(Global, Linear)
    cumulative_precipitation[column] += scale * precipitation[column]
    cumulative_evapotranspiration[column] +=
        scale * latent_heat_flux[column] * inverse_latent_heat_density
    cumulative_runoff[column] += scale * runoff[column]
end

function _accumulate_terrarium_land_water_budget!(
    callback::GlobalAtmosphereDiagnosticsCallback,
    state,
    model,
    step_seconds,
)
    water_density = model.constants.material.density_water
    precipitation = Oceananigans.interior(state.inputs.rainfall)
    runoff = Oceananigans.interior(state.surface_runoff)
    backend = KernelAbstractions.get_backend(
        callback.column_cumulative_land_precipitation,
    )
    scale = water_density * step_seconds
    if hasproperty(state, :ground_water_flux) &&
       hasproperty(state, :transpiration_water_flux)
        _accumulate_land_water_components_kernel!(backend)(
            callback.column_cumulative_land_precipitation,
            callback.column_cumulative_land_evapotranspiration,
            callback.column_cumulative_land_surface_runoff,
            precipitation,
            Oceananigans.interior(state.ground_water_flux),
            Oceananigans.interior(state.transpiration_water_flux),
            runoff,
            scale;
            ndrange = length(callback.column_cumulative_land_precipitation),
        )
    else
        latent_heat = model.constants.thermodynamics.latent_heat_vaporization
        _accumulate_land_water_latent_kernel!(backend)(
            callback.column_cumulative_land_precipitation,
            callback.column_cumulative_land_evapotranspiration,
            callback.column_cumulative_land_surface_runoff,
            precipitation,
            Oceananigans.interior(state.latent_heat_flux),
            runoff,
            scale,
            inv(latent_heat * water_density);
            ndrange = length(callback.column_cumulative_land_precipitation),
        )
    end
    return nothing
end

# CPU half of ReadyESM's concrete SpeedyWeather/Terrarium adapter. The 0.1.6
# compatibility baseline explicitly disables its new prognostic snow process,
# so the selected no-canopy/direct-runoff hydrology does not consume the
# otherwise-populated `snowfall` input. Apply snowfall as immediate liquid-water
# equivalent rather than deleting its mass. `_run_terrarium_land!` also
# integrates the actually applied water fluxes at the native land timestep.
function _timestep_terrarium_cpu!(
    vars::SpeedyWeather.Variables,
    land::_ReadyESMTerrariumLand,
    model::SpeedyWeather.PrimitiveWetModel,
)
    state = vars.prognostic.land.terrarium
    terrarium_model = land.model
    constants = terrarium_model.constants
    NF = eltype(state)
    mask = land.model.grid.mask.data

    air_temperature = vars.grid.temperature[mask, end]
    humidity = vars.grid.humidity[mask, end]
    pressure = vars.grid.pressure[mask]
    wind = vars.parameterizations.surface_wind_speed[mask]
    rain = vars.parameterizations.rain_rate[mask]
    snow = vars.parameterizations.snow_rate[mask]
    shortwave_down = vars.parameterizations.surface_shortwave_down[mask]
    longwave_down = vars.parameterizations.surface_longwave_down[mask]

    inputs = state.inputs
    Terrarium.set!(inputs.air_temperature, air_temperature)
    Terrarium.set!(inputs.air_temperature, inputs.air_temperature - NF(273.15))
    Terrarium.set!(inputs.air_pressure, pressure)
    Terrarium.set!(inputs.air_pressure, exp(inputs.air_pressure))
    Terrarium.set!(inputs.specific_humidity, humidity)
    Terrarium.set!(inputs.rainfall, rain + snow)
    Terrarium.set!(inputs.snowfall, snow)
    Terrarium.set!(inputs.windspeed, wind)
    Terrarium.set!(inputs.surface_shortwave_down, shortwave_down)
    Terrarium.set!(inputs.surface_longwave_down, longwave_down)

    integrator = Terrarium.ModelIntegrator(
        state.clock,
        terrarium_model,
        Terrarium.InputSources(NF),
        state,
        land.initializers,
    )
    diagnostics = haskey(model.callbacks, :global_atmosphere_diagnostics) ?
        model.callbacks[:global_atmosphere_diagnostics] : nothing
    coupling_period = _terrarium_coupling_period(vars, model)
    _run_terrarium_land!(
        integrator,
        coupling_period,
        land.Δt;
        diagnostics,
    )

    vars.prognostic.land.soil_temperature[mask] .=
        Oceananigans.interior(state.skin_temperature) .+ NF(273.15)
    vars.prognostic.land.soil_moisture[mask] .=
        @view Oceananigans.interior(state.saturation_water_ice)[:, 1, end]
    if haskey(vars.prognostic.land, :sensible_heat_flux)
        vars.prognostic.land.sensible_heat_flux[mask] .=
            Oceananigans.interior(state.sensible_heat_flux)
    end
    if haskey(vars.prognostic.land, :surface_humidity_flux)
        vars.prognostic.land.surface_humidity_flux[mask] .=
            Oceananigans.interior(state.latent_heat_flux) ./
            constants.thermodynamics.latent_heat_vaporization
    end
    if haskey(vars.parameterizations, :surface_longwave_up)
        vars.parameterizations.surface_longwave_up[mask] .=
            Oceananigans.interior(state.surface_longwave_up)
    end
    if haskey(vars.parameterizations, :surface_shortwave_up)
        vars.parameterizations.surface_shortwave_up[mask] .=
            Oceananigans.interior(state.surface_shortwave_up)
    end
    return nothing
end

# Terrarium applies infiltration (m s⁻¹ of liquid-water depth) as the top flux
# boundary condition of relative saturation. The saturation equation requires
# that flux divided by soil porosity; applying it directly stores only the
# porosity fraction of infiltrated water. Specialize the exact vegetation-free,
# snow-free, zero-organic-carbon Richards land model used by ReadyESM and
# construct its 0.1.6 StateVariables with the required conversion. The narrow
# type signature leaves every other Terrarium model on its upstream path.
function Terrarium.StateVariables(
    model::Terrarium.LandModel{
        NF,
        GridType,
        Nothing,
        Terrarium.SoilEnergyWaterCarbon{
            NF,
            Stratigraphy,
            Energy,
            Terrarium.SoilHydrology{
                NF,
                Terrarium.RichardsEq,
                HydrologyClosure,
                HydraulicProperties,
                VWCForcing,
            },
            Terrarium.ConstantSoilCarbonDensity{NF},
        },
        Nothing,
        SurfaceEnergy,
        Terrarium.SurfaceHydrology{
            NF,
            CanopyInterception,
            Evapotranspiration,
            Terrarium.DirectSurfaceRunoff{NF},
        },
        Atmosphere,
        Initializer,
        Timestepper,
    };
    clock = Terrarium.Clock(time = zero(NF)),
    boundary_conditions = (;),
    fields = (;),
    input_variables = (),
) where {
    NF,
    GridType <: Terrarium.AbstractLandGrid{NF},
    Stratigraphy <: Terrarium.SoilStratigraphy{NF, 1},
    Energy <: Terrarium.AbstractSoilThermodynamics{NF},
    HydrologyClosure,
    HydraulicProperties,
    VWCForcing,
    SurfaceEnergy <: Terrarium.AbstractSurfaceEnergyBalance{NF},
    CanopyInterception,
    Evapotranspiration,
    Atmosphere <: Terrarium.AbstractAtmosphere{NF},
    Initializer <: Terrarium.AbstractInitializer{NF},
    Timestepper <: Terrarium.AbstractTimeStepper{NF},
}
    iszero(model.soil.biogeochem.ρ_soc) || error(
        "ReadyESM's porosity-scaled infiltration boundary requires zero " *
        "prognostic soil organic carbon",
    )
    horizon = only(model.soil.strat.horizons)
    horizon isa Terrarium.ConstantSoilHorizon || error(
        "ReadyESM's porosity-scaled infiltration boundary requires one " *
        "constant soil horizon",
    )
    horizon.porosity isa Terrarium.ConstantSoilPorosity || error(
        "ReadyESM's porosity-scaled infiltration boundary requires constant porosity",
    )
    porosity = horizon.porosity.mineral_porosity
    zero(NF) < porosity <= one(NF) || error(
        "ReadyESM soil porosity must be in (0, 1], got $porosity",
    )

    grid = Terrarium.get_grid(model)
    interface_variables = Terrarium.interface_variables(model)
    variables = Terrarium.Variables(
        Terrarium.variables(model)...,
        interface_variables...,
        input_variables...,
    )
    ground_heat_flux = Terrarium.initialize(
        variables.ground_heat_flux,
        grid,
        clock,
        fields,
        boundary_conditions,
    )
    infiltration = Terrarium.initialize(
        variables.infiltration,
        grid,
        clock,
        fields,
        boundary_conditions,
    )
    soil_heat_flux = ground_heat_flux
    soil_heat_flux_boundary = Terrarium.SoilHeatFlux(soil_heat_flux)
    saturation_infiltration = -(infiltration * inv(porosity))
    infiltration_boundary = Terrarium.InfiltrationFlux(
        saturation_infiltration,
    )
    boundaries = Terrarium.merge_boundary_conditions(
        boundary_conditions,
        soil_heat_flux_boundary,
        infiltration_boundary,
    )
    initialized_fields = merge(
        (; ground_heat_flux, infiltration, soil_heat_flux),
        fields,
    )
    return Terrarium.StateVariables(
        variables,
        grid;
        clock,
        timestepper = Terrarium.get_timestepper(model),
        model,
        boundary_conditions = boundaries,
        fields = initialized_fields,
    )
end

# Upstream Terrarium diagnoses DirectSurfaceRunoff drainage as +S / τ_r. If
# the runoff scheme reaches the Richards tendency, the package adds that positive
# value to the prognostic surface-excess-water reservoir while simultaneously
# exporting it as runoff. Drainage is a loss from the reservoir, so apply the
# physically consistent negative tendency for the exact RichardsEq +
# DirectSurfaceRunoff combination used by ReadyESM. Keeping this method
# narrowly typed avoids changing any other Terrarium hydrology scheme.
@inline function _surface_excess_drainage_tendency(surface_excess_water, runoff)
    return -Terrarium.compute_surface_drainage(runoff, surface_excess_water)
end

Base.@propagate_inbounds function Terrarium.compute_surface_excess_water_tendency(
    i,
    j,
    k,
    grid,
    clock,
    fields,
    hydrology::Terrarium.SoilHydrology{NF, Terrarium.RichardsEq},
    runoff::Terrarium.DirectSurfaceRunoff,
) where {NF}
    surface_excess_water = fields.surface_excess_water[i, j, k]
    top_layer = Terrarium.get_field_grid(grid).Nz
    tendency = _surface_excess_drainage_tendency(surface_excess_water, runoff)
    return ifelse(k == top_layer, tendency, zero(NF))
end

# Terrarium's upstream Richards kernel is launched over XYZ for the 3-D soil
# saturation tendency, but it also updates the 2-D surface-excess-water
# tendency inside every vertical thread. Those threads race on the same XY
# element on a GPU. In particular, the single negative drainage update can be
# overwritten by one of the other layers' zero updates. Run the two prognostic
# fields on their natural domains instead: XYZ for soil water and XY exactly
# once for surface water.
KernelAbstractions.@kernel function _richards_saturation_tendency_kernel!(
    saturation_tendency,
    grid,
    clock,
    fields,
    hydrology,
    stratigraphy,
    biogeochemistry,
    constants,
    evapotranspiration,
)
    i, j, k = @index(Global, NTuple)
    Terrarium.compute_saturation_tendency!(
        saturation_tendency,
        i,
        j,
        k,
        grid,
        clock,
        fields,
        hydrology,
        stratigraphy,
        biogeochemistry,
        constants,
        evapotranspiration,
    )
end

KernelAbstractions.@kernel function _direct_surface_drainage_tendency_kernel!(
    surface_water_tendency,
    grid,
    fields,
    runoff,
)
    i, j = @index(Global, NTuple)
    @inbounds begin
        surface_water = fields.surface_excess_water[i, j, 1]
        # DirectSurfaceRunoff is the complete tendency for this prognostic
        # reservoir. Assign it rather than accumulating into an old tendency
        # value; this also makes the intended split-kernel contract explicit.
        surface_water_tendency[i, j, 1] =
            _surface_excess_drainage_tendency(surface_water, runoff)
    end
end

function _compute_split_richards_tendencies!(
    state,
    grid,
    hydrology,
    soil,
    constants,
    evapotranspiration,
    runoff,
    fields,
)
    stratigraphy = Terrarium.get_stratigraphy(soil)
    biogeochemistry = Terrarium.get_biogeochemistry(soil)
    tendencies = Terrarium.tendency_fields(state, hydrology)
    Terrarium.launch!(
        grid,
        Terrarium.XYZ,
        _richards_saturation_tendency_kernel!,
        tendencies.saturation_water_ice,
        state.clock,
        fields,
        hydrology,
        stratigraphy,
        biogeochemistry,
        constants,
        evapotranspiration,
    )
    Terrarium.launch!(
        grid,
        Terrarium.XY,
        _direct_surface_drainage_tendency_kernel!,
        tendencies.surface_excess_water,
        fields,
        runoff,
    )
    return nothing
end

# Apply the split kernel to every ReadyESM Richards + direct-runoff land stack.
# Terrarium 0.1.6 stores both bare-ground and prescribed-vegetation
# evapotranspiration as physical liquid-water-depth fluxes.
function Terrarium.compute_tendencies!(
    state,
    grid,
    hydrology::Terrarium.SoilHydrology{NF, Terrarium.RichardsEq},
    soil::Terrarium.AbstractSoil,
    constants::Terrarium.PhysicalConstants,
    evapotranspiration::Terrarium.AbstractEvapotranspiration,
    runoff::Terrarium.DirectSurfaceRunoff,
    args...,
) where {NF}
    stratigraphy = Terrarium.get_stratigraphy(soil)
    biogeochemistry = Terrarium.get_biogeochemistry(soil)
    fields = Terrarium.get_fields(
        state,
        hydrology,
        stratigraphy,
        biogeochemistry,
        evapotranspiration,
    )
    return _compute_split_richards_tendencies!(
        state,
        grid,
        hydrology,
        soil,
        constants,
        evapotranspiration,
        runoff,
        fields,
    )
end

# In LandModel, infiltration is applied explicitly through the top flux
# boundary condition of `saturation_water_ice`. Terrarium's Richards divergence
# nevertheless evaluates an additional Darcy flux at the top face from the
# auxiliary pressure-head halo. That halo is not part of the closure and
# defaults to zero, producing a large artificial influx into an otherwise
# constant-head column. Suppress only this duplicate top-face flux on the
# bounded vertical column grid; all interior and lower faces retain the package
# implementation.
Base.@propagate_inbounds function Terrarium.darcy_flux(
    i,
    j,
    k,
    grid::Oceananigans.Grids.RectilinearGrid{
        FT,
        TX,
        TY,
        Oceananigans.Grids.Bounded,
    },
    pressure_head,
    hydraulic_conductivity,
) where {FT, TX, TY}
    k > grid.Nz && return zero(FT)
    return invoke(
        Terrarium.darcy_flux,
        Tuple{Any, Any, Any, Any, Any, Any},
        i,
        j,
        k,
        grid,
        pressure_head,
        hydraulic_conductivity,
    )
end

# Terrarium's LandModel tendency driver calls the coupled soil tendency without
# forwarding the surface hydrology's evapotranspiration or runoff processes.
# The Richards implementation therefore receives `nothing` for both optional
# arguments: bare-ground evaporation is omitted from soil water and, crucially,
# surface drainage is never applied to `surface_excess_water` even though the
# auxiliary runoff field exports that water. Thread both processes into the
# prognostic soil hydrology for ReadyESM's exact bare-ground Richards + direct
# runoff stack. Other LandModel combinations retain the upstream method.
function Terrarium.compute_tendencies!(
    state,
    model::Terrarium.LandModel{NF, GridType, Nothing, Soil, Nothing},
) where {NF, GridType <: Terrarium.AbstractLandGrid{NF}, Soil}
    grid = Terrarium.get_grid(model)
    surface_hydrology = model.surface_hydrology
    soil = model.soil

    if !(
        soil.hydrology.vertical_flow isa Terrarium.RichardsEq &&
        surface_hydrology.surface_runoff isa Terrarium.DirectSurfaceRunoff
    )
        return invoke(
            Terrarium.compute_tendencies!,
            Tuple{Any, Terrarium.LandModel},
            state,
            model,
        )
    end

    Terrarium.compute_tendencies!(state, grid, surface_hydrology)
    Terrarium.compute_tendencies!(
        state,
        grid,
        soil.hydrology,
        soil,
        model.constants,
        surface_hydrology.evapotranspiration,
        surface_hydrology.surface_runoff,
    )
    Terrarium.compute_tendencies!(
        state,
        grid,
        soil.biogeochem,
        soil,
        model.constants,
    )
    Terrarium.compute_tendencies!(
        state,
        grid,
        soil.energy,
        soil,
        model.constants,
    )
    Terrarium.compute_tendencies!(state, grid, model.vegetation)
    return nothing
end

# Terrarium expresses its top soil-water boundary condition as the vertically
# reduced unary operation `-infiltration`. Oceananigans' generic AbstractArray
# boundary accessor supplies only two indices, but UnaryOperation defines its
# fast/device-safe indexing with three. The missing reduced-operation method
# otherwise falls through dynamic axes machinery and fails CUDA compilation.
Base.@propagate_inbounds function Oceananigans.BoundaryConditions.getbc(
    condition::Oceananigans.AbstractOperations.UnaryOperation{LX, LY, Nothing},
    i::Integer,
    j::Integer,
    grid::Oceananigans.Grids.AbstractGrid,
    args...,
) where {LX, LY}
    return condition[i, j, 1]
end

# Terrarium 0.1.2 interpolates the values of a failed assertion into a String in
# `saturation_to_pressure!`.  CUDA must compile both branches and therefore
# rejects the string allocation even when every value is finite.  The adapted
# Both adapted ColumnGrid and ColumnRingGrid have `Nothing` as their architecture
# parameter, so the two wrappers below are selected only in device code; CPU
# execution retains Terrarium's diagnostic assertion.  Non-finite device results
# are still caught by Terrarium's field debug hooks and ReadyESM's output
# validation.  Keeping the calculation in one helper ensures the focused
# standalone-column conservation gate exercises the same implementation as the
# production ring grid.
Base.@propagate_inbounds function _terrarium_device_saturation_to_pressure!(
    pressure_head,
    i,
    j,
    k,
    grid,
    fields,
    closure::Terrarium.SoilSaturationPressureClosure,
    hydrology::Terrarium.SoilHydrology{NF, Terrarium.RichardsEq},
    strat::Terrarium.AbstractStratigraphy,
    bgc::Terrarium.AbstractSoilBiogeochemistry,
) where {NF}
    fgrid = Terrarium.get_field_grid(grid)
    sat = fields.saturation_water_ice[i, j, k]
    z = Terrarium.znode(i, j, k, fgrid, Terrarium.Center(), Terrarium.Center(), Terrarium.Center())
    z_ref = Terrarium.znode(
        i,
        j,
        fgrid.Nz + 1,
        fgrid,
        Terrarium.Center(),
        Terrarium.Center(),
        Terrarium.Face(),
    )
    inv_swrc = inv(Terrarium.get_swrc(hydrology))
    por = Terrarium.porosity(i, j, k, grid, fields, strat, bgc)
    sat_res = Terrarium.residual_saturation(Terrarium.get_hydraulic_properties(hydrology))
    psi_m = inv_swrc(max(sat * por, sat_res); θsat = por)
    psi_z = z - z_ref
    z_water_table = fields.water_table[i, j, 1]
    psi_h = max(zero(z_water_table), z_water_table - z)
    pressure_head[i, j, k] = psi_h + psi_m + psi_z
    return nothing
end

Base.@propagate_inbounds function Terrarium.saturation_to_pressure!(
    pressure_head,
    i,
    j,
    k,
    grid::Terrarium.ColumnGrid{FT, Nothing},
    fields,
    closure::Terrarium.SoilSaturationPressureClosure,
    hydrology::Terrarium.SoilHydrology{NF, Terrarium.RichardsEq},
    strat::Terrarium.AbstractStratigraphy,
    bgc::Terrarium.AbstractSoilBiogeochemistry,
) where {FT, NF}
    return _terrarium_device_saturation_to_pressure!(
        pressure_head, i, j, k, grid, fields, closure, hydrology, strat, bgc,
    )
end

Base.@propagate_inbounds function Terrarium.saturation_to_pressure!(
    pressure_head,
    i,
    j,
    k,
    grid::Terrarium.ColumnRingGrid{FT, Nothing},
    fields,
    closure::Terrarium.SoilSaturationPressureClosure,
    hydrology::Terrarium.SoilHydrology{NF, Terrarium.RichardsEq},
    strat::Terrarium.AbstractStratigraphy,
    bgc::Terrarium.AbstractSoilBiogeochemistry,
) where {FT, NF}
    return _terrarium_device_saturation_to_pressure!(
        pressure_head, i, j, k, grid, fields, closure, hydrology, strat, bgc,
    )
end

# Terrarium's ColumnRingGrid stores only a Boolean mask, while CUDA.jl deliberately
# does not support logical indexing. Cache the equivalent integer indices once per
# model. All subsequent atmosphere/land exchange remains on the selected device.
const _TERRARIUM_LAND_INDEX_CACHE = IdDict{Any, Any}()

_is_gpu_terrarium(land) = Terrarium.architecture(land.model.grid) isa Terrarium.GPU

function _terrarium_land_indices(land)
    mask = land.model.grid.mask.data
    return get!(_TERRARIUM_LAND_INDEX_CACHE, mask) do
        host_indices = findall(Array(mask))
        Terrarium.on_architecture(Terrarium.architecture(land.model.grid), host_indices)
    end
end

KernelAbstractions.@kernel function _initialize_terrarium_land_kernel!(
    soil_temperature,
    soil_moisture,
    atmospheric_temperature,
    terrarium_temperature,
    terrarium_skin_temperature,
    saturation,
    land_indices,
    kelvin_offset,
    lowest_atmosphere_level,
    top_soil_level,
)
    column = @index(Global, Linear)
    grid_point = land_indices[column]
    initial_temperature_kelvin = atmospheric_temperature[
        grid_point,
        lowest_atmosphere_level,
    ]
    initial_temperature_celsius = initial_temperature_kelvin - kelvin_offset
    for level in 1:top_soil_level
        terrarium_temperature[column, 1, level] = initial_temperature_celsius
    end
    terrarium_skin_temperature[column, 1, 1] = initial_temperature_celsius
    soil_temperature[grid_point] = initial_temperature_kelvin
    soil_moisture[grid_point] = saturation[column, 1, top_soil_level]
end

KernelAbstractions.@kernel function _gather_terrarium_forcing_kernel!(
    air_temperature,
    air_pressure,
    specific_humidity,
    rainfall,
    snowfall,
    windspeed,
    shortwave_down,
    longwave_down,
    atmospheric_temperature,
    atmospheric_humidity,
    log_surface_pressure,
    surface_wind_speed,
    rain_rate,
    snow_rate,
    atmospheric_shortwave_down,
    atmospheric_longwave_down,
    land_indices,
    kelvin_offset,
    lowest_level,
)
    column = @index(Global, Linear)
    grid_point = land_indices[column]
    air_temperature[column] = atmospheric_temperature[grid_point, lowest_level] - kelvin_offset
    air_pressure[column] = exp(log_surface_pressure[grid_point])
    specific_humidity[column] = atmospheric_humidity[grid_point, lowest_level]
    snowfall_water_equivalent = snow_rate[grid_point]
    rainfall[column] = rain_rate[grid_point] + snowfall_water_equivalent
    snowfall[column] = snowfall_water_equivalent
    windspeed[column] = surface_wind_speed[grid_point]
    shortwave_down[column] = atmospheric_shortwave_down[grid_point]
    longwave_down[column] = atmospheric_longwave_down[grid_point]
end

KernelAbstractions.@kernel function _scatter_terrarium_surface_kernel!(
    destination,
    source,
    land_indices,
    multiplier,
    offset,
)
    column = @index(Global, Linear)
    grid_point = land_indices[column]
    destination[grid_point] = multiplier * source[column] + offset
end

KernelAbstractions.@kernel function _scatter_terrarium_bottom_kernel!(
    destination,
    source,
    land_indices,
    multiplier,
    offset,
    bottom_level,
)
    column = @index(Global, Linear)
    grid_point = land_indices[column]
    destination[grid_point] = multiplier * source[column, 1, bottom_level] + offset
end

function _scatter_terrarium_surface!(
    backend,
    destination,
    source,
    land_indices;
    multiplier = one(eltype(source)),
    offset = zero(eltype(source)),
)
    _scatter_terrarium_surface_kernel!(backend)(
        destination,
        Oceananigans.interior(source),
        land_indices,
        multiplier,
        offset;
        ndrange = length(land_indices),
    )
    return nothing
end


"""
Terrarium 0.1.6-compatible initialization for SpeedyWeather's wet primitive model.

The upstream extension targets the older state-returning initializer and, on a
GPU, also uses a `CuArray{Bool}` as an array index. This concrete-land method is
more specific than the upstream abstract method on both CPU and GPU: it builds
the new model-owned integrator and uses an integer-index kernel for the device
scatter.
"""
function SpeedyWeather.initialize!(
    vars::SpeedyWeather.Variables,
    land::_ReadyESMTerrariumLand,
    model::SpeedyWeather.PrimitiveWetModel,
)
    state = vars.prognostic.land.terrarium
    NF = eltype(vars.prognostic.land.soil_temperature)
    state.clock.time = vars.prognostic.clock.time
    integrator = Terrarium.ModelIntegrator(
        state.clock,
        land.model,
        Terrarium.InputSources(NF),
        state,
        land.initializers,
    )
    Terrarium.initialize!(integrator)

    land_indices = _terrarium_land_indices(land)
    temperature = Oceananigans.interior(state.temperature)
    saturation = Oceananigans.interior(state.saturation_water_ice)
    @assert length(land.model.grid.mask.data) ==
        length(vars.prognostic.land.soil_temperature)
    @assert length(land_indices) == size(temperature, 1)

    backend = KernelAbstractions.get_backend(land_indices)
    _scatter_terrarium_bottom_kernel!(backend)(
        vars.prognostic.land.soil_temperature,
        temperature,
        land_indices,
        one(NF),
        NF(273.15),
        size(temperature, 3);
        ndrange = length(land_indices),
    )
    _scatter_terrarium_bottom_kernel!(backend)(
        vars.prognostic.land.soil_moisture,
        saturation,
        land_indices,
        one(NF),
        zero(NF),
        size(saturation, 3);
        ndrange = length(land_indices),
    )
    KernelAbstractions.synchronize(backend)
    return nothing
end

"""
Initialize Terrarium soil and skin temperature from SpeedyWeather's resolved
lowest atmospheric level after the spectral initial conditions have been
transformed to grid space.

SpeedyWeather initializes land before atmospheric initial conditions, so its
ordinary land hook cannot read `vars.grid.temperature`: that array is still
zero at that point. ReadyESM invokes this second-stage initializer from the
atmospheric diagnostics callback, whose initialization occurs after the
simulation-level initial transform and before the first timestep.
"""
function _initialize_terrarium_temperature_from_atmosphere!(vars, model)
    haskey(vars.prognostic.land, :terrarium) || return nothing
    land = model.land
    land isa _ReadyESMAbstractTerrariumLandModel || return nothing

    state = vars.prognostic.land.terrarium
    NF = eltype(state)
    land_indices = _terrarium_land_indices(land)
    temperature = Oceananigans.interior(state.temperature)
    skin_temperature = Oceananigans.interior(state.skin_temperature)
    saturation = Oceananigans.interior(state.saturation_water_ice)
    backend = KernelAbstractions.get_backend(land_indices)
    prescribed_soil_state = all(
        name -> haskey(land.initializers, name),
        (:temperature, :saturation_water_ice, :skin_temperature),
    )

    if prescribed_soil_state && _is_gpu_terrarium(land)
        # An exact land reanalysis has already initialized the Terrarium
        # fields. Preserve that vertical profile and only refresh the two
        # single-layer mirrors consumed by SpeedyWeather.
        _scatter_terrarium_bottom_kernel!(backend)(
            vars.prognostic.land.soil_temperature,
            temperature,
            land_indices,
            one(NF),
            NF(273.15),
            size(temperature, 3);
            ndrange = length(land_indices),
        )
        _scatter_terrarium_bottom_kernel!(backend)(
            vars.prognostic.land.soil_moisture,
            saturation,
            land_indices,
            one(NF),
            zero(NF),
            size(saturation, 3);
            ndrange = length(land_indices),
        )
        KernelAbstractions.synchronize(backend)
    elseif prescribed_soil_state
        bottom_soil_level = size(temperature, 3)
        for (column, grid_point) in enumerate(land_indices)
            vars.prognostic.land.soil_temperature[grid_point] =
                temperature[column, 1, bottom_soil_level] + NF(273.15)
            vars.prognostic.land.soil_moisture[grid_point] =
                saturation[column, 1, bottom_soil_level]
        end
    elseif _is_gpu_terrarium(land)
        _initialize_terrarium_land_kernel!(backend)(
            vars.prognostic.land.soil_temperature,
            vars.prognostic.land.soil_moisture,
            vars.grid.temperature,
            temperature,
            skin_temperature,
            saturation,
            land_indices,
            NF(273.15),
            size(vars.grid.temperature, 2),
            size(temperature, 3);
            ndrange = length(land_indices),
        )
        KernelAbstractions.synchronize(backend)
    else
        lowest_atmosphere_level = size(vars.grid.temperature, 2)
        top_soil_level = size(temperature, 3)
        for (column, grid_point) in enumerate(land_indices)
            initial_temperature_kelvin =
                vars.grid.temperature[grid_point, lowest_atmosphere_level]
            initial_temperature_celsius = initial_temperature_kelvin - NF(273.15)
            temperature[column, 1, :] .= initial_temperature_celsius
            skin_temperature[column, 1, 1] = initial_temperature_celsius
            vars.prognostic.land.soil_temperature[grid_point] =
                initial_temperature_kelvin
            vars.prognostic.land.soil_moisture[grid_point] =
                saturation[column, 1, top_soil_level]
        end
    end

    soil = land.model.soil
    energy = soil.energy
    Terrarium.invclosure!(
        state,
        land.model.grid,
        energy.closure,
        energy,
        soil,
        land.model.constants,
    )
    return nothing
end


function _trace_terrarium_substeps(state)
    enabled = lowercase(
        get(ENV, "READYESM_TRACE_TERRARIUM_SUBSTEPS", "false"),
    ) in ("1", "true", "yes", "on")
    enabled || return false
    first_iteration = parse(
        Int,
        get(ENV, "READYESM_TRACE_TERRARIUM_AFTER_ITERATION", "0"),
    )
    first_iteration >= 0 || throw(
        ArgumentError("READYESM_TRACE_TERRARIUM_AFTER_ITERATION must be non-negative"),
    )
    return state.clock.iteration >= first_iteration
end

function _trace_speedy_phases(vars)
    enabled = lowercase(
        get(ENV, "READYESM_TRACE_SPEEDY_PHASES", "false"),
    ) in ("1", "true", "yes", "on")
    enabled || return false
    first_timestep = parse(
        Int,
        get(ENV, "READYESM_TRACE_SPEEDY_AFTER_TIMESTEP", "0"),
    )
    first_timestep >= 0 || throw(
        ArgumentError("READYESM_TRACE_SPEEDY_AFTER_TIMESTEP must be non-negative"),
    )
    return vars.prognostic.clock.timestep_counter >= first_timestep
end

function _synchronize_speedy_handoff()
    return lowercase(
        get(ENV, "READYESM_SYNC_SPEEDY_HANDOFF", "false"),
    ) in ("1", "true", "yes", "on")
end

const COUPLED_GPU_COMPONENT_HANDOFF_PROVENANCE =
    "speedyweather_completion_before_land_ice_ocean_v1"

@inline _coupled_gpu_handoff_required(has_terrarium, cuda_functional) =
    has_terrarium && cuda_functional

"""
Complete the SpeedyWeather GPU component before NumericalEarth advances land,
sea ice and ocean. Those frameworks launch work asynchronously and
NumericalEarth's sequential component loop does not itself establish a
cross-framework completion contract. Without this boundary the exact coupled
T31 graph can enter ClimaSeaIce momentum while prior atmospheric device work is
still outstanding; isolated and sanitizer runs are clean, while a completion
at this one physical handoff removes the first-step stall.
"""
function _complete_speedy_component_handoff!(has_terrarium)
    cuda_functional = CUDA.functional()
    required = _coupled_gpu_handoff_required(has_terrarium, cuda_functional)
    requested = _synchronize_speedy_handoff() && cuda_functional
    (required || requested) && CUDA.synchronize()
    return required || requested
end

function _synchronize_speedy_internal_phases()
    return lowercase(
        get(ENV, "READYESM_SYNC_SPEEDY_INTERNAL_PHASES", "false"),
    ) in ("1", "true", "yes", "on")
end

function _terrarium_coupling_period(vars, model)
    clock = vars.prognostic.clock
    full_period = clock.Δt
    time_stepping = model.time_stepping
    startup = time_stepping.first_step_euler &&
              clock.timestep_counter == 0
    startup || return full_period

    full_milliseconds = Millisecond(full_period).value
    first_half = Millisecond(full_milliseconds ÷ 2)
    return clock.time == clock.start ?
        first_half : Millisecond(full_milliseconds - first_half.value)
end

function _terrarium_substep_schedule(period, land_timestep)
    period_seconds = Millisecond(period).value / 1000
    step_seconds = Float64(land_timestep)
    nfull = floor(Int, period_seconds / step_seconds)
    remainder = period_seconds - nfull * step_seconds
    abs(remainder) <= 8eps(period_seconds) && (remainder = 0.0)
    return (; nfull, remainder)
end

function _terrarium_host_extrema(field)
    values = Float64.(Array(Oceananigans.interior(field)))
    limits = extrema(values)
    maximum_value, maximum_index = findmax(values)
    minimum_value, minimum_index = findmin(values)
    return (; limits, minimum_value, minimum_index, maximum_value, maximum_index)
end

function _speedy_host_summary(field)
    # CUDA arrays expose an internal `data` property that is a GPUArrays.DataRef
    # rather than a copyable array. Conversely, SpeedyWeather's triangular
    # spectral wrappers are AbstractArrays whose `.data` is the copyable
    # CuArray; calling `Array(wrapper)` otherwise takes a scalar-indexing path.
    # Unwrap only when `.data` itself implements AbstractArray.
    data = hasproperty(field, :data) ? getproperty(field, :data) : nothing
    source = data isa AbstractArray ? data : field
    values = Array(source)
    finite = isfinite.(values)
    nonfinite_indices = findall(!, finite)
    finite_values = values[finite]
    absolute_limits = isempty(finite_values) ?
        (NaN, NaN) : extrema(Float64.(abs.(finite_values)))
    value_limits = if isempty(finite_values)
        (NaN, NaN)
    elseif eltype(values) <: Real
        extrema(Float64.(finite_values))
    else
        absolute_limits
    end
    return (;
        value_limits,
        absolute_limits,
        nonfinite_count = length(nonfinite_indices),
        first_nonfinite = if isempty(nonfinite_indices)
            nothing
        else
            index = first(nonfinite_indices)
            index isa CartesianIndex ? Tuple(index) : index
        end,
    )
end

@inline function _speedy_optional_summary(fields, name::Symbol)
    return hasproperty(fields, name) ?
        _speedy_host_summary(getproperty(fields, name)) : nothing
end

function _trace_speedy_step_phase!(vars, phase)
    CUDA.functional() && CUDA.synchronize()
    summaries = (
        grid_vorticity = _speedy_host_summary(vars.grid.vorticity),
        grid_divergence = _speedy_host_summary(vars.grid.divergence),
        grid_temperature = _speedy_host_summary(vars.grid.temperature),
        grid_humidity = _speedy_host_summary(vars.grid.humidity),
        grid_log_surface_pressure = _speedy_host_summary(vars.grid.pressure),
        grid_u = _speedy_host_summary(vars.grid.u),
        grid_v = _speedy_host_summary(vars.grid.v),
        grid_vorticity_tendency =
            _speedy_host_summary(vars.tendencies.grid.vorticity),
        grid_divergence_tendency =
            _speedy_host_summary(vars.tendencies.grid.divergence),
        grid_temperature_tendency =
            _speedy_host_summary(vars.tendencies.grid.temperature),
        grid_humidity_tendency =
            _speedy_host_summary(vars.tendencies.grid.humidity),
        grid_pressure_tendency =
            _speedy_host_summary(vars.tendencies.grid.pressure),
        grid_u_tendency = _speedy_host_summary(vars.tendencies.grid.u),
        grid_v_tendency = _speedy_host_summary(vars.tendencies.grid.v),
        spectral_vorticity =
            _speedy_host_summary(vars.prognostic.vorticity),
        spectral_divergence =
            _speedy_host_summary(vars.prognostic.divergence),
        spectral_temperature =
            _speedy_host_summary(vars.prognostic.temperature),
        spectral_humidity = _speedy_host_summary(vars.prognostic.humidity),
        spectral_log_surface_pressure =
            _speedy_host_summary(vars.prognostic.pressure),
        spectral_vorticity_tendency =
            _speedy_host_summary(vars.tendencies.vorticity),
        spectral_divergence_tendency =
            _speedy_host_summary(vars.tendencies.divergence),
        spectral_temperature_tendency =
            _speedy_host_summary(vars.tendencies.temperature),
        spectral_humidity_tendency =
            _speedy_host_summary(vars.tendencies.humidity),
        spectral_pressure_tendency =
            _speedy_host_summary(vars.tendencies.pressure),
        sea_surface_temperature = _speedy_host_summary(
            vars.prognostic.ocean.sea_surface_temperature,
        ),
        sea_ice_concentration = _speedy_host_summary(
            vars.prognostic.ocean.sea_ice_concentration,
        ),
        prescribed_ocean_sensible_heat = _speedy_optional_summary(
            vars.prognostic.ocean,
            :sensible_heat_flux,
        ),
        prescribed_ocean_humidity_flux = _speedy_optional_summary(
            vars.prognostic.ocean,
            :surface_humidity_flux,
        ),
        surface_air_temperature = _speedy_host_summary(
            vars.parameterizations.surface_air_temperature,
        ),
        surface_sensible_heat = _speedy_host_summary(
            vars.parameterizations.sensible_heat_flux,
        ),
        surface_humidity_flux = _speedy_host_summary(
            vars.parameterizations.surface_humidity_flux,
        ),
    )
    println(
        "SPEEDY_STEP phase=$phase atmospheric_time=$(vars.prognostic.clock.time) " *
        "timestep_counter=$(vars.prognostic.clock.timestep_counter) summaries=$summaries",
    )
    flush(stdout)
    return all(
        isnothing(summary) || summary.nonfinite_count == 0 for
        summary in values(summaries)
    )
end

function SpeedyWeather.timestep!(
    vars::SpeedyWeather.Variables,
    dt::Real,
    model::SpeedyWeather.PrimitiveWetModel,
    lf1::Integer = 2,
    lf2::Integer = 2,
)
    has_terrarium = haskey(vars.prognostic.land, :terrarium)
    trace = _trace_speedy_phases(vars) ||
            (has_terrarium &&
             _trace_terrarium_substeps(vars.prognostic.land.terrarium))
    synchronize_internal_phases = _synchronize_speedy_internal_phases()
    if !trace && !synchronize_internal_phases
        invoke(
            SpeedyWeather.timestep!,
            Tuple{
                SpeedyWeather.Variables,
                Real,
                SpeedyWeather.PrimitiveEquation,
                Integer,
                Integer,
            },
            vars,
            dt,
            model,
            lf1,
            lf2,
        )
        _complete_speedy_component_handoff!(has_terrarium)
        return nothing
    end

    function finish_phase!(phase)
        if trace
            finite = _trace_speedy_step_phase!(vars, phase)
            finite || error(
                "SpeedyWeather first crossed a finite trace bound during $phase",
            )
        elseif synchronize_internal_phases && CUDA.functional()
            CUDA.synchronize()
        end
        return nothing
    end

    !isnothing(model.feedback) && model.feedback.nans_detected && return nothing
    finish_phase!(:entry)
    SpeedyWeather.reset_tendencies!(vars)
    finish_phase!(:after_reset_tendencies)

    if !model.dynamics_only
        SpeedyWeather.greenhouse_gases_time_step!(vars, model)
        finish_phase!(:after_greenhouse_gases)
        SpeedyWeather.parameterization_tendencies!(vars, model)
        finish_phase!(:after_parameterizations)
        SpeedyWeather.ocean_timestep!(vars, model)
        finish_phase!(:after_speedy_ocean)
        SpeedyWeather.sea_ice_timestep!(vars, model)
        finish_phase!(:after_speedy_sea_ice)
        SpeedyWeather.land_timestep!(vars, model)
        finish_phase!(:after_land)
    end

    if model.dynamics
        SpeedyWeather.dynamics_tendencies!(vars, lf2, model)
        finish_phase!(:after_dynamics_tendencies)
        SpeedyWeather.implicit_correction!(vars, model.implicit, model)
        finish_phase!(:after_implicit_correction)
    else
        SpeedyWeather.parameterization_tendencies_only!(vars, model)
        finish_phase!(:after_parameterization_transform)
    end

    SpeedyWeather.horizontal_diffusion!(
        vars,
        model.horizontal_diffusion,
        model,
    )
    finish_phase!(:after_horizontal_diffusion)
    SpeedyWeather.leapfrog!(vars, dt, lf1, model)
    finish_phase!(:after_leapfrog)
    SpeedyWeather.transform!(vars, lf2, model)
    finish_phase!(:after_transform)

    lf2 == 2 && SpeedyWeather.particle_advection!(vars, model)
    finish_phase!(:complete)
    _complete_speedy_component_handoff!(has_terrarium)
    return nothing
end

function _trace_terrarium_substep!(state, phase)
    air_temperature = _terrarium_host_extrema(state.inputs.air_temperature)
    air_pressure = _terrarium_host_extrema(state.inputs.air_pressure)
    specific_humidity = _terrarium_host_extrema(state.inputs.specific_humidity)
    windspeed = _terrarium_host_extrema(state.inputs.windspeed)
    shortwave_down = _terrarium_host_extrema(
        state.inputs.surface_shortwave_down,
    )
    longwave_down = _terrarium_host_extrema(
        state.inputs.surface_longwave_down,
    )
    temperature = _terrarium_host_extrema(state.temperature)
    skin_temperature = _terrarium_host_extrema(state.skin_temperature)
    saturation = _terrarium_host_extrema(state.saturation_water_ice)
    liquid_water_fraction = _terrarium_host_extrema(
        state.liquid_water_fraction,
    )
    saturation_tendency = _terrarium_host_extrema(
        state.tendencies.saturation_water_ice,
    )
    internal_energy = _terrarium_host_extrema(state.internal_energy)
    energy_tendency = _terrarium_host_extrema(state.tendencies.internal_energy)
    pressure_head = _terrarium_host_extrema(state.pressure_head)
    conductivity = _terrarium_host_extrema(state.hydraulic_conductivity)
    surface_water = _terrarium_host_extrema(state.surface_excess_water)
    rainfall = _terrarium_host_extrema(state.inputs.rainfall)
    snowfall = _terrarium_host_extrema(state.inputs.snowfall)
    infiltration = _terrarium_host_extrema(state.infiltration)
    runoff = _terrarium_host_extrema(state.surface_runoff)
    ground_heat = _terrarium_host_extrema(state.ground_heat_flux)
    latent_heat = _terrarium_host_extrema(state.latent_heat_flux)

    println(
        "TERRARIUM_SUBSTEP phase=$phase iteration=$(state.clock.iteration) " *
        "time=$(state.clock.time) " *
        "air_temperature=$(air_temperature.limits) " *
        "air_pressure=$(air_pressure.limits) " *
        "specific_humidity=$(specific_humidity.limits) " *
        "windspeed=$(windspeed.limits) " *
        "shortwave_down=$(shortwave_down.limits) " *
        "longwave_down=$(longwave_down.limits) " *
        "temperature=$(temperature.limits) Tmax_index=$(temperature.maximum_index) " *
        "skin=$(skin_temperature.limits) saturation=$(saturation.limits) " *
        "liquid_water_fraction=$(liquid_water_fraction.limits) " *
        "dSdt=$(saturation_tendency.limits) " *
        "dSdt_max_index=$(saturation_tendency.maximum_index) " *
        "energy=$(internal_energy.limits) dUdt=$(energy_tendency.limits) " *
        "pressure=$(pressure_head.limits) conductivity=$(conductivity.limits) " *
        "surface_water=$(surface_water.limits) " *
        "surface_water_max_index=$(surface_water.maximum_index) " *
        "rainfall=$(rainfall.limits) snowfall=$(snowfall.limits) " *
        "infiltration=$(infiltration.limits) " *
        "runoff=$(runoff.limits) ground_heat=$(ground_heat.limits) " *
        "latent_heat=$(latent_heat.limits)",
    )
    flush(stdout)

    all_finite = all(
        isfinite,
        (
            air_temperature.limits...,
            air_pressure.limits...,
            specific_humidity.limits...,
            windspeed.limits...,
            shortwave_down.limits...,
            longwave_down.limits...,
            temperature.limits...,
            skin_temperature.limits...,
            saturation.limits...,
            liquid_water_fraction.limits...,
            saturation_tendency.limits...,
            internal_energy.limits...,
            energy_tendency.limits...,
            pressure_head.limits...,
            conductivity.limits...,
            surface_water.limits...,
            rainfall.limits...,
            snowfall.limits...,
            infiltration.limits...,
            runoff.limits...,
            ground_heat.limits...,
            latent_heat.limits...,
        ),
    )
    physically_bounded = temperature.limits[1] >= -100 &&
                         temperature.limits[2] <= 80 &&
                         skin_temperature.limits[1] >= -100 &&
                         skin_temperature.limits[2] <= 80 &&
                         saturation.limits[1] >= -1e-6 &&
                         saturation.limits[2] <= 1 + 1e-6 &&
                         liquid_water_fraction.limits[1] >= -1e-6 &&
                         liquid_water_fraction.limits[2] <= 1 + 1e-6 &&
                         surface_water.limits[1] >= -1e-8 &&
                         surface_water.limits[2] <= 1 &&
                         air_temperature.limits[1] >= -150 &&
                         air_temperature.limits[2] <= 100 &&
                         air_pressure.limits[1] > 0 &&
                         air_pressure.limits[2] <= 200_000 &&
                         specific_humidity.limits[1] >= -1e-6 &&
                         specific_humidity.limits[2] <= 0.1 &&
                         windspeed.limits[1] >= 0 &&
                         windspeed.limits[2] <= 500 &&
                         shortwave_down.limits[1] >= -1e-3 &&
                         shortwave_down.limits[2] <= 2_000 &&
                         longwave_down.limits[1] >= -1e-3 &&
                         longwave_down.limits[2] <= 1_000
    return all_finite && physically_bounded
end

function _terrarium_trace_sync!(phase)
    println("TERRARIUM_SYNC phase=$phase status=waiting")
    flush(stdout)
    CUDA.functional() && CUDA.synchronize()
    println("TERRARIUM_SYNC phase=$phase status=complete")
    flush(stdout)
    return nothing
end

function _terrarium_trace_state_sync!(state, phase)
    _terrarium_trace_sync!(phase)
    summaries = (
        saturation = _speedy_host_summary(
            Oceananigans.interior(state.saturation_water_ice),
        ),
        liquid_fraction = _speedy_host_summary(
            Oceananigans.interior(state.liquid_water_fraction),
        ),
        internal_energy = _speedy_host_summary(
            Oceananigans.interior(state.internal_energy),
        ),
        temperature = _speedy_host_summary(
            Oceananigans.interior(state.temperature),
        ),
        pressure_head = _speedy_host_summary(
            Oceananigans.interior(state.pressure_head),
        ),
        hydraulic_conductivity = _speedy_host_summary(
            Oceananigans.interior(state.hydraulic_conductivity),
        ),
        surface_excess_water = _speedy_host_summary(
            Oceananigans.interior(state.surface_excess_water),
        ),
        saturation_tendency = _speedy_host_summary(
            Oceananigans.interior(
                state.tendencies.saturation_water_ice,
            ),
        ),
        energy_tendency = _speedy_host_summary(
            Oceananigans.interior(state.tendencies.internal_energy),
        ),
    )
    println(
        "TERRARIUM_STATE phase=$phase iteration=$(state.clock.iteration) " *
        "summaries=$summaries",
    )
    flush(stdout)
    all_finite = all(
        summary.nonfinite_count == 0 for summary in values(summaries)
    )
    all_finite || error(
        "Terrarium first generated a non-finite traced state during $phase " *
        "at land iteration $(state.clock.iteration)",
    )
    return nothing
end

"""
Run Terrarium's ordinary final auxiliary sequence with a global CUDA barrier
after every independently launched stage. This path is enabled only by the
substep-trace environment switch and is deliberately more synchronized than
production. It identifies the exact asynchronous kernel group that launches a
device assertion without changing the operations or their order.
"""
function _compute_terrarium_auxiliary_traced!(
    state,
    model::Terrarium.LandModel,
)
    grid = Terrarium.get_grid(model)
    soil = model.soil
    surface_hydrology = model.surface_hydrology
    surface_energy_balance = model.surface_energy_balance
    constants = model.constants
    atmosphere = model.atmosphere
    vegetation = model.vegetation
    snow = model.snow

    _terrarium_trace_state_sync!(state, :before_final_auxiliary)

    Terrarium.compute_auxiliary!(state, grid, atmosphere)
    _terrarium_trace_state_sync!(state, :atmosphere_auxiliary)

    Terrarium.compute_auxiliary!(
        state,
        grid,
        soil.hydrology,
        soil,
        constants,
    )
    _terrarium_trace_state_sync!(state, :soil_hydrology_auxiliary)
    Terrarium.compute_auxiliary!(
        state,
        grid,
        soil.biogeochem,
        soil,
        constants,
    )
    _terrarium_trace_state_sync!(state, :soil_biogeochemistry_auxiliary)
    Terrarium.compute_auxiliary!(state, grid, soil.energy, soil, constants)
    _terrarium_trace_state_sync!(state, :soil_energy_auxiliary)

    Terrarium.compute_auxiliary!(state, grid, snow, constants)
    _terrarium_trace_state_sync!(state, :snow_auxiliary)

    Terrarium.compute_auxiliary!(
        state,
        grid,
        vegetation,
        constants,
        atmosphere,
        soil,
    )
    _terrarium_trace_state_sync!(state, :vegetation_auxiliary)

    Terrarium.compute_auxiliary!(
        state,
        grid,
        surface_hydrology.canopy_interception,
        atmosphere,
    )
    _terrarium_trace_state_sync!(state, :canopy_interception_auxiliary)
    Terrarium.compute_auxiliary!(
        state,
        grid,
        surface_hydrology.evapotranspiration,
        surface_hydrology.canopy_interception,
        constants,
        atmosphere,
        soil,
        vegetation,
        snow,
    )
    _terrarium_trace_state_sync!(state, :evapotranspiration_auxiliary)
    Terrarium.compute_auxiliary!(
        state,
        grid,
        surface_hydrology.surface_runoff,
        surface_hydrology.canopy_interception,
        soil,
        snow,
    )
    _terrarium_trace_state_sync!(state, :surface_runoff_auxiliary)

    Terrarium.compute_auxiliary!(
        state,
        grid,
        surface_energy_balance,
        vegetation,
        snow,
    )
    _terrarium_trace_state_sync!(state, :surface_energy_auxiliary)
    return nothing
end

"""
Run the exact Forward-Euler Terrarium substep with a CUDA barrier after each
independently launched process group. This is an opt-in failure-localization
path only: the untraced production branch continues to call
`Terrarium.timestep!` directly.

The extra barriers identify which asynchronous 2-D or 3-D land kernel launches
a device exception. Operations and their order match Terrarium's
`update_state! -> boundary conditions -> explicit_step! -> timestep! ->
closure! -> tick!` sequence,
including ReadyESM's narrowly specialized bare-ground Richards tendencies.
"""
function _terrarium_forward_euler_traced!(
    integrator,
    timestepper::Terrarium.ForwardEuler,
    Δt,
)
    state = integrator.state
    model = integrator.model
    grid = Terrarium.get_grid(model)
    soil = model.soil
    surface_hydrology = model.surface_hydrology

    Terrarium.reset_tendencies!(state)
    _terrarium_trace_state_sync!(state, :reset_tendencies)
    Terrarium.update_inputs!(state, grid, integrator.inputs)
    _terrarium_trace_state_sync!(state, :update_inputs)

    _compute_terrarium_auxiliary_traced!(state, model)
    Terrarium.compute_boundary_conditions!(state, model)
    _terrarium_trace_state_sync!(state, :boundary_conditions)

    Terrarium.compute_tendencies!(state, grid, surface_hydrology)
    _terrarium_trace_state_sync!(state, :surface_hydrology_tendencies)
    Terrarium.compute_tendencies!(
        state,
        grid,
        soil.hydrology,
        soil,
        model.constants,
        surface_hydrology.evapotranspiration,
        surface_hydrology.surface_runoff,
    )
    _terrarium_trace_state_sync!(state, :soil_hydrology_tendencies)
    Terrarium.compute_tendencies!(
        state,
        grid,
        soil.biogeochem,
        soil,
        model.constants,
    )
    _terrarium_trace_state_sync!(state, :soil_biogeochemistry_tendencies)
    Terrarium.compute_tendencies!(
        state,
        grid,
        soil.energy,
        soil,
        model.constants,
    )
    _terrarium_trace_state_sync!(state, :soil_energy_tendencies)
    Terrarium.compute_tendencies!(state, grid, model.vegetation)
    _terrarium_trace_state_sync!(state, :vegetation_tendencies)

    prognostic_names = Terrarium.prognostic_names(state)
    Terrarium.explicit_step!(state, grid, timestepper, Δt, prognostic_names)
    _terrarium_trace_state_sync!(state, :explicit_step)
    Terrarium.timestep!(state, model, timestepper, Δt)
    _terrarium_trace_state_sync!(state, :model_timestep_correction)

    Terrarium.closure!(state, model)
    _terrarium_trace_state_sync!(state, :closure)

    Terrarium.tick!(state.clock, Δt)
    return nothing
end

function _run_terrarium_land!(
    integrator,
    period,
    land_timestep;
    diagnostics = nothing,
)
    trace = _trace_terrarium_substeps(integrator.state)
    schedule = _terrarium_substep_schedule(period, land_timestep)
    timestepper = Terrarium.get_timestepper(integrator.model)
    if !trace && isnothing(diagnostics) && iszero(schedule.remainder)
        Terrarium.run!(integrator; period, Δt = land_timestep)
        return nothing
    end

    if trace
        physical = _trace_terrarium_substep!(integrator.state, :before)
        physical || error(
            "Terrarium received non-finite or nonphysical forcing before " *
            "land iteration $(integrator.state.clock.iteration)",
        )
    end
    substep = 0
    for _ in 1:schedule.nfull
        substep += 1
        if trace && integrator.model isa Terrarium.LandModel &&
           timestepper isa Terrarium.ForwardEuler
            _terrarium_forward_euler_traced!(
                integrator,
                timestepper,
                land_timestep,
            )
        else
            Terrarium.timestep!(integrator, land_timestep; finalize = false)
        end
        isnothing(diagnostics) || _accumulate_terrarium_land_water_budget!(
            diagnostics,
            integrator.state,
            integrator.model,
            land_timestep,
        )
        if trace
            physical = _trace_terrarium_substep!(
                integrator.state,
                Symbol("after_$substep"),
            )
            physical || error(
                "Terrarium first crossed its finite/physical trace bounds at " *
                "land iteration $(integrator.state.clock.iteration)",
            )
        end
    end
    if !iszero(schedule.remainder)
        substep += 1
        if trace && integrator.model isa Terrarium.LandModel &&
           timestepper isa Terrarium.ForwardEuler
            _terrarium_forward_euler_traced!(
                integrator,
                timestepper,
                schedule.remainder,
            )
        else
            Terrarium.timestep!(
                integrator,
                schedule.remainder;
                finalize = false,
            )
        end
        isnothing(diagnostics) || _accumulate_terrarium_land_water_budget!(
            diagnostics,
            integrator.state,
            integrator.model,
            schedule.remainder,
        )
        if trace
            physical = _trace_terrarium_substep!(
                integrator.state,
                Symbol("after_$(substep)_remainder"),
            )
            physical || error(
                "Terrarium first crossed its finite/physical trace bounds at " *
                "land iteration $(integrator.state.clock.iteration)",
            )
        end
    end
    if trace && integrator.model isa Terrarium.LandModel
        _compute_terrarium_auxiliary_traced!(
            integrator.state,
            integrator.model,
        )
        physical = _trace_terrarium_substep!(
            integrator.state,
            :after_final_auxiliary,
        )
        physical || error(
            "Terrarium final auxiliary evaluation crossed its finite/physical " *
            "trace bounds at land iteration $(integrator.state.clock.iteration)",
        )
    else
        Terrarium.compute_auxiliary!(integrator.state, integrator.model)
    end
    return nothing
end


"""
GPU-specialized SpeedyWeather/Terrarium timestep exchange.

Atmospheric forcing is gathered into Terrarium with one kernel. Updated land
state and fluxes are then scattered into full SpeedyWeather ring fields without
host copies or scalar/Boolean device indexing.
"""
function SpeedyWeather.timestep!(
    vars::SpeedyWeather.Variables,
    land::_ReadyESMTerrariumLand,
    model::SpeedyWeather.PrimitiveWetModel,
)
    if !_is_gpu_terrarium(land)
        return _timestep_terrarium_cpu!(vars, land, model)
    end

    state = vars.prognostic.land.terrarium
    tmodel = land.model
    NF = eltype(state)
    inputs = state.inputs
    land_indices = _terrarium_land_indices(land)
    backend = KernelAbstractions.get_backend(land_indices)

    _gather_terrarium_forcing_kernel!(backend)(
        Oceananigans.interior(inputs.air_temperature),
        Oceananigans.interior(inputs.air_pressure),
        Oceananigans.interior(inputs.specific_humidity),
        Oceananigans.interior(inputs.rainfall),
        Oceananigans.interior(inputs.snowfall),
        Oceananigans.interior(inputs.windspeed),
        Oceananigans.interior(inputs.surface_shortwave_down),
        Oceananigans.interior(inputs.surface_longwave_down),
        vars.grid.temperature,
        vars.grid.humidity,
        vars.grid.pressure,
        vars.parameterizations.surface_wind_speed,
        vars.parameterizations.rain_rate,
        vars.parameterizations.snow_rate,
        vars.parameterizations.surface_shortwave_down,
        vars.parameterizations.surface_longwave_down,
        land_indices,
        NF(273.15),
        size(vars.grid.temperature, 2);
        ndrange = length(land_indices),
    )
    KernelAbstractions.synchronize(backend)

    integrator = Terrarium.ModelIntegrator(
        state.clock,
        tmodel,
        Terrarium.InputSources(NF),
        state,
        land.initializers,
    )
    coupling_period = _terrarium_coupling_period(vars, model)
    diagnostics = haskey(model.callbacks, :global_atmosphere_diagnostics) ?
        model.callbacks[:global_atmosphere_diagnostics] : nothing
    _run_terrarium_land!(
        integrator,
        coupling_period,
        land.Δt;
        diagnostics,
    )
    trace = _trace_terrarium_substeps(state)
    trace && _terrarium_trace_sync!(:land_integrator_complete)

    one_NF = one(NF)
    zero_NF = zero(NF)
    _scatter_terrarium_surface!(
        backend,
        vars.prognostic.land.soil_temperature,
        state.skin_temperature,
        land_indices;
        offset = NF(273.15),
    )
    trace && _terrarium_trace_sync!(:soil_temperature_scatter)
    saturation = Oceananigans.interior(state.saturation_water_ice)
    _scatter_terrarium_bottom_kernel!(backend)(
        vars.prognostic.land.soil_moisture,
        saturation,
        land_indices,
        one_NF,
        zero_NF,
        size(saturation, 3);
        ndrange = length(land_indices),
    )
    trace && _terrarium_trace_sync!(:soil_moisture_scatter)
    if haskey(vars.prognostic.land, :sensible_heat_flux)
        _scatter_terrarium_surface!(
            backend,
            vars.prognostic.land.sensible_heat_flux,
            state.sensible_heat_flux,
            land_indices,
        )
        trace && _terrarium_trace_sync!(:sensible_heat_scatter)
    end
    if haskey(vars.prognostic.land, :surface_humidity_flux)
        latent_heat = tmodel.constants.thermodynamics.latent_heat_vaporization
        _scatter_terrarium_surface!(
            backend,
            vars.prognostic.land.surface_humidity_flux,
            state.latent_heat_flux,
            land_indices;
            multiplier = inv(latent_heat),
        )
        trace && _terrarium_trace_sync!(:humidity_flux_scatter)
    end
    if haskey(vars.parameterizations, :surface_longwave_up)
        _scatter_terrarium_surface!(
            backend,
            vars.parameterizations.surface_longwave_up,
            state.surface_longwave_up,
            land_indices,
        )
        trace && _terrarium_trace_sync!(:longwave_up_scatter)
    end
    if haskey(vars.parameterizations, :surface_shortwave_up)
        _scatter_terrarium_surface!(
            backend,
            vars.parameterizations.surface_shortwave_up,
            state.surface_shortwave_up,
            land_indices,
        )
        trace && _terrarium_trace_sync!(:shortwave_up_scatter)
    end
    KernelAbstractions.synchronize(backend)
    return nothing
end
