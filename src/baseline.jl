const RG = SpeedyWeather.RingGrids

# A coupled radiative boundary must be single-valued. RRTMGP and Terrarium
# separately compute the upward surface flux from the same downward flux; if
# their albedo or emissivity differs, the atmosphere and land budgets create or
# destroy the difference. These values match SpeedyWeather's current bare-land
# albedo and the RRTMGP/ocean emissivity used elsewhere in the coupled stack.
const COUPLED_LAND_SURFACE_ALBEDO = 0.4f0
const COUPLED_SURFACE_EMISSIVITY = 0.98f0

function _coupled_land_surface_energy_balance(::Type{FT}) where {FT}
    surface_properties = Terrarium.ConstantAlbedo(
        FT;
        albedo = FT(COUPLED_LAND_SURFACE_ALBEDO),
        emissivity = FT(COUPLED_SURFACE_EMISSIVITY),
    )
    return Terrarium.SurfaceEnergyBalance(
        FT;
        albedo = surface_properties,
        # Terrarium 0.1.6 changed this default from 2 to 1 W m⁻¹ K⁻¹.
        # Preserve the already-qualified ReadyESM surface-energy identity.
        skin_temperature = Terrarium.ImplicitSkinTemperature(FT; κₛ = FT(2)),
    )
end

_coupled_atmosphere_albedo(spectral_grid) = SpeedyWeather.OceanLandAlbedo(
    spectral_grid;
    land = SpeedyWeather.GlobalConstantAlbedo(
        spectral_grid;
        albedo = spectral_grid.NF(COUPLED_LAND_SURFACE_ALBEDO),
    ),
)

_speedy_architecture(config::ExperimentConfig) =
    config.device == :gpu ? SpeedyWeather.GPU() : SpeedyWeather.CPU()

_ocean_architecture(config::ExperimentConfig) =
    config.device == :gpu ? Oceananigans.GPU() : Oceananigans.CPU()

_terrarium_architecture(config::ExperimentConfig) =
    config.device == :gpu ? Terrarium.GPU() : Terrarium.CPU()

function _speedy_cuda_graphs()
    value = lowercase(get(ENV, "READYESM_SPEEDY_CUDA_GRAPHS", "true"))
    value in ("1", "true", "yes", "on") && return true
    value in ("0", "false", "no", "off") && return false
    throw(
        ArgumentError(
            "READYESM_SPEEDY_CUDA_GRAPHS must be a boolean, got $(repr(value))",
        ),
    )
end

_speedy_spectral_transform(spectral_grid) =
    SpeedyWeather.SpeedyTransforms.SpectralTransform(
        spectral_grid;
        cuda_graphs = _speedy_cuda_graphs(),
    )

"""Scale-selective atmospheric damping configured in physical time units."""
_atmosphere_horizontal_diffusion(config::ExperimentConfig, spectral_grid) =
    HyperDiffusion(
        spectral_grid;
        time_scale = Second(round(Int, 3600 * config.atmosphere_hyperdiffusion_hours)),
        time_scale_div = Second(round(
            Int,
            3600 * config.atmosphere_divergence_hyperdiffusion_hours,
        )),
    )

"""Quadratic protection against unresolved all-level and surface winds."""
_atmosphere_speed_limit_drag(config::ExperimentConfig, spectral_grid) =
    CoupledSurfaceSpeedLimitDrag(
        spectral_grid;
        speed_limit = config.atmosphere_speed_limit_ms,
        drag = config.atmosphere_speed_limit_drag_per_m,
        surface_speed_limit = config.atmosphere_surface_speed_limit_ms,
        surface_drag = config.atmosphere_surface_speed_limit_drag_per_m,
    )

# Terrarium's `SoilHydrology(FT)` default is `NoFlow`: saturation is an
# auxiliary field and its tendency method is intentionally empty. A coupled
# climate land surface needs prognostic soil water, supplied by RichardsEq.
# Use the Van Genuchten retention/conductivity pair exercised by Terrarium's
# coupled Richards tests. The package's generic Richards defaults combine a
# Brooks--Corey curve with conductivity linear in saturation. At the 5 cm
# production surface layer that pair has an explicit stability limit of only a
# few seconds in dry soil and catastrophically overshoots with our 60 s step.
function _prognostic_soil_hydrology(::Type{FT}) where {FT}
    swrc = Terrarium.VanGenuchten(α = FT(2), n = FT(2))
    hydraulic_properties = Terrarium.ConstantSoilHydraulics(
        FT;
        swrc,
        unsat_hydraulic_cond = Terrarium.UnsatKVanGenuchten(FT),
    )
    return Terrarium.SoilHydrology(
        FT,
        Terrarium.RichardsEq();
        hydraulic_properties,
    )
end

"""Hydrostatic soil saturation above and below a prescribed water table."""
struct HydrostaticSoilSaturation{FT, SW} <: Terrarium.AbstractInitializer{FT}
    water_table_elevation::FT
    porosity::FT
    swrc::SW
end

@inline function (profile::HydrostaticSoilSaturation)(x, z)
    z <= profile.water_table_elevation && return one(profile.porosity)
    matric_head = profile.water_table_elevation - z
    volumetric_water = profile.swrc(matric_head; θsat = profile.porosity)
    return clamp(volumetric_water / profile.porosity, zero(profile.porosity), one(profile.porosity))
end

# Oceananigans dispatches coordinate-dependent `set!` only for values whose
# type is a subtype of `Function`; a callable struct otherwise reaches its
# broadcast-data fallback. Wrap the isbits profile in an ordinary closure at
# the Terrarium field-initializer boundary.
_hydrostatic_soil_field_initializer(profile::HydrostaticSoilSaturation) =
    (x, z) -> profile(x, z)

function Terrarium.initialize!(
    state,
    ::Terrarium.AbstractModel,
    profile::HydrostaticSoilSaturation,
)
    Terrarium.set!(
        state.saturation_water_ice,
        _hydrostatic_soil_field_initializer(profile),
    )
    return nothing
end

function _aligned_water_table_elevation(
    ::Type{FT},
    spacing::Terrarium.AbstractVerticalSpacing,
    target_depth,
) where {FT}
    layer_thickness = Terrarium.get_spacing(spacing)
    faces = FT.(vcat(-reverse(cumsum(layer_thickness)), zero(eltype(layer_thickness))))
    target_elevation = -FT(target_depth)
    return faces[argmin(abs.(faces .- target_elevation))]
end

function _atmosphere_large_scale_condensation(
    config::ExperimentConfig,
    spectral_grid,
)
    if config.atmosphere_large_scale_precipitation == :upstream_observed
        return ObservedImplicitCondensation(spectral_grid)
    elseif config.atmosphere_large_scale_precipitation ==
           :net_column_conservative
        return NetColumnImplicitCondensation(spectral_grid)
    else
        error(
            "unsupported large-scale precipitation scheme " *
            "$(config.atmosphere_large_scale_precipitation)",
        )
    end
end

"""
Construct the first ReadyESM baseline: SpeedyWeather's moist primitive-equation
atmosphere, a prognostic slab ocean, and concentration-only thermodynamic sea ice.

The model deliberately uses SpeedyWeather's simplified radiation. CO₂ and aerosol
values are retained in configuration/provenance but do not affect this baseline;
the RRTMGP experiment rejects this radiation mode distinction explicitly.
"""
function build_baseline(config::ExperimentConfig)
    validate(config)
    config.forcing.radiation == :speedy_simplified || throw(
        ArgumentError(
            "build_baseline only supports radiation=speedy_simplified; " *
                "use the RRTMGP experiment for concentration-dependent forcing",
        ),
    )

    spectral_grid = SpectralGrid(
        ; NF = Float32,
        trunc = config.truncation,
        Grid = RG.FullGaussianGrid,
        nlayers = config.nlayers,
        architecture = _speedy_architecture(config),
    )
    time_stepping = Leapfrog(
        spectral_grid;
        Δt_at_T31 = Second(round(Int, 60 * config.atmosphere_timestep_at_t31_minutes)),
    )
    horizontal_diffusion = _atmosphere_horizontal_diffusion(config, spectral_grid)
    vertical_diffusion = _atmosphere_vertical_diffusion(config, spectral_grid)
    drag = _atmosphere_speed_limit_drag(config, spectral_grid)
    spectral_transform = _speedy_spectral_transform(spectral_grid)
    hole_filling = LayerConservativeHumidityHoleFilling(spectral_grid)
    large_scale_condensation =
        _atmosphere_large_scale_condensation(config, spectral_grid)
    convection = _atmosphere_convection(config, spectral_grid)
    ocean = SlabOcean(
        spectral_grid;
        mixed_layer_depth = Float32(config.mixed_layer_depth_m),
    )
    sea_ice = ThermodynamicSeaIce(spectral_grid)
    model = PrimitiveWetModel(
        spectral_grid;
        ocean,
        sea_ice,
        time_stepping,
        horizontal_diffusion,
        vertical_diffusion,
        drag,
        spectral_transform,
        hole_filling,
        large_scale_condensation,
        convection,
    )
    add!(model.callbacks, :global_surface_temperature => GlobalSurfaceTemperatureCallback(spectral_grid))
    return initialize!(model)
end

"""Run a configured baseline and return the simulation plus collected diagnostics."""
function run_baseline!(config::ExperimentConfig; steps::Union{Nothing, Int} = nothing)
    simulation = build_baseline(config)
    if isnothing(steps)
        milliseconds = round(Int, config.duration_days * 86_400_000)
        SpeedyWeather.run!(simulation; period = Millisecond(milliseconds))
    else
        steps > 0 || throw(ArgumentError("steps must be positive"))
        SpeedyWeather.run!(simulation; steps)
    end
    return simulation, collect_diagnostics(simulation, config)
end

"""Build an RRTMGP atmosphere model, optionally exposing prescribed ocean flux hooks."""
function _rrtmgp_atmosphere_model(config::ExperimentConfig; coupled_ocean = false)
    validate(config)
    config.forcing.radiation in (:rrtmgp_clear_sky, :rrtmgp_all_sky) || throw(
        ArgumentError("build_rrtmgp requires radiation=rrtmgp_clear_sky or rrtmgp_all_sky"),
    )
    spectral_grid = SpectralGrid(
        ; NF = Float32,
        trunc = config.truncation,
        Grid = RG.FullGaussianGrid,
        nlayers = config.nlayers,
        architecture = _speedy_architecture(config),
    )
    time_stepping = Leapfrog(
        spectral_grid;
        Δt_at_T31 = Second(round(Int, 60 * config.atmosphere_timestep_at_t31_minutes)),
    )
    horizontal_diffusion = _atmosphere_horizontal_diffusion(config, spectral_grid)
    vertical_diffusion = _atmosphere_vertical_diffusion(config, spectral_grid)
    drag = _atmosphere_speed_limit_drag(config, spectral_grid)
    spectral_transform = _speedy_spectral_transform(spectral_grid)
    hole_filling = LayerConservativeHumidityHoleFilling(spectral_grid)
    large_scale_condensation =
        _atmosphere_large_scale_condensation(config, spectral_grid)
    prognostic_cloud_condensate =
        config.forcing.cloud_scheme == :prognostic_condensate ?
        PrognosticCloudCondensateState(
            spectral_grid;
            retention_fraction =
                config.forcing.cloud_condensate_retention_fraction,
            residence_time_hours =
                config.forcing.cloud_condensate_residence_time_hours,
        ) : nothing
    base_convection = _atmosphere_convection(config, spectral_grid)
    convection = isnothing(prognostic_cloud_condensate) ? base_convection :
        PrognosticCondensateConvection(
            base_convection,
            prognostic_cloud_condensate,
        )
    ocean = if coupled_ocean
        PrescribedOcean(spectral_grid)
    else
        SlabOcean(
            spectral_grid;
            mixed_layer_depth = Float32(config.mixed_layer_depth_m),
        )
    end
    # The externally coupled ClimaSeaIce state must still be represented inside
    # SpeedyWeather.  In particular, OceanSeaIceAlbedo reads this prognostic
    # concentration when constructing RRTMGP's lower-boundary albedo.  Using
    # `nothing` silently makes RRTMGP see open ocean everywhere even while the
    # ocean/ice component applies a high ice albedo at the same surface.
    sea_ice = coupled_ocean ? PrescribedSeaIce(spectral_grid) :
        ThermodynamicSeaIce(spectral_grid)
    initial_conditions = _atmosphere_initial_conditions(config, spectral_grid)
    shortwave_radiation = RRTMGPNoShortwave()
    longwave_radiation = RRTMGPClearSkyRadiation(
        spectral_grid;
        co2_ppm = config.forcing.co2_ppm,
        aerosol_optical_depth_550nm = config.forcing.aerosol_optical_depth_550nm,
        surface_emissivity = COUPLED_SURFACE_EMISSIVITY,
        solve_every_n_steps = config.forcing.radiation_every_n_steps,
        cloud_scheme = config.forcing.cloud_scheme,
        cloud_humidity_search_min_sigma =
            config.forcing.cloud_humidity_search_min_sigma,
        cloud_liquid_water_path_gm2 = config.forcing.cloud_liquid_water_path_gm2,
        cloud_ice_water_path_gm2 = config.forcing.cloud_ice_water_path_gm2,
        prognostic_cloud_condensate,
    )
    # Diagnostics and any embedded land model inspect this mask during model
    # construction, before SpeedyWeather's later initialize! pass.
    land_sea_mask = EarthLandSeaMask(spectral_grid)
    SpeedyWeather.load_mask!(land_sea_mask)
    model = if config.land_model == :terrarium
        land_mask = land_sea_mask.mask .> 0
        soil_spacing = Terrarium.ExponentialSpacing(
            ; N = config.terrarium_soil_layers,
            Δz_min = 0.05,
            Δz_max = config.terrarium_max_layer_thickness_m,
        )
        column_grid = Terrarium.ColumnRingGrid(
            _terrarium_architecture(config),
            Float32,
            soil_spacing,
            spectral_grid.grid,
            land_mask,
        )
        FT = eltype(column_grid)
        soil_hydrology = _prognostic_soil_hydrology(FT)
        # The package default jumps from 0.75 to saturation at 5 m, causing a
        # large artificial Richards redistribution at startup. Align the water
        # table to the nearest grid face, then initialize a continuous
        # zero-head-gradient profile against that exact face and the default
        # homogeneous 0.49 porosity. This must be the model hydrology
        # initializer: Terrarium applies user field initializers first and its
        # model initializer afterwards.
        water_table_elevation = _aligned_water_table_elevation(
            FT,
            soil_spacing,
            5,
        )
        initial_saturation = HydrostaticSoilSaturation(
            water_table_elevation,
            FT(0.49),
            Terrarium.get_swrc(soil_hydrology),
        )
        soil_initializer = Terrarium.SoilInitializer(
            FT;
            hydrology = initial_saturation,
        )
        terrarium_field_initializers = (;)
        if config.terrarium_initial_conditions == :era5_instantaneous
            # ERA5 state arrays are applied as user field initializers. Use a
            # no-op model initializer so Terrarium does not overwrite them
            # after that pass; process initialization still diagnoses all
            # hydrological and thermal auxiliary fields from the supplied
            # prognostic state.
            era5_land = _era5_land_initializers(
                config,
                spectral_grid,
                column_grid,
                water_table_elevation,
            )
            terrarium_field_initializers = era5_land.initializers
            soil_initializer = Terrarium.DefaultInitializer(FT)
            @info "ERA5 Terrarium initialization prepared" era5_land.audit
        end
        soil = Terrarium.SoilEnergyWaterCarbon(
            FT;
            hydrology = soil_hydrology,
        )
        surface_hydrology = Terrarium.SurfaceHydrology(
            FT;
            canopy_interception = Terrarium.NoCanopyInterception(FT),
            evapotranspiration = _terrarium_evapotranspiration(config, FT),
            surface_runoff = Terrarium.DirectSurfaceRunoff(FT),
        )
        surface_energy_balance =
            _coupled_land_surface_energy_balance(FT)
        terrarium_model = Terrarium.LandModel(
            column_grid;
            initializer = soil_initializer,
            vegetation = nothing,
            # Land snow was not prognostic in the qualified Terrarium 0.1.2
            # stack. Enable the new 0.1.6 snow model only in a separately
            # attributable experiment with restart and storage-ledger support.
            snow = nothing,
            soil,
            surface_energy_balance,
            surface_hydrology,
            timestepper = Terrarium.ForwardEuler(FT),
        )
        land = SpeedyWeather.LandModel(
            spectral_grid,
            terrarium_model;
            timestepper = Terrarium.get_timestepper(terrarium_model),
            initializers = terrarium_field_initializers,
            Δt = config.terrarium_timestep_seconds,
        )
        ocean_heat_flux = coupled_ocean ?
            PrescribedOceanHeatFlux(spectral_grid) : SurfaceOceanHeatFlux(spectral_grid)
        ocean_humidity_flux = coupled_ocean ?
            PrescribedOceanHumidityFlux(spectral_grid) : SurfaceOceanHumidityFlux(spectral_grid)
        surface_heat_flux = SurfaceHeatFlux(
            spectral_grid;
            ocean = ocean_heat_flux,
            land = PrescribedLandHeatFlux(),
        )
        surface_humidity_flux = SurfaceHumidityFlux(
            spectral_grid;
            ocean = ocean_humidity_flux,
            land = PrescribedLandHumidityFlux(),
        )
        albedo = _coupled_atmosphere_albedo(spectral_grid)
        PrimitiveWetModel(
            spectral_grid;
            ocean,
            sea_ice,
            initial_conditions,
            land,
            land_sea_mask,
            albedo,
            surface_heat_flux,
            surface_humidity_flux,
            shortwave_radiation,
            longwave_radiation,
            time_stepping,
            horizontal_diffusion,
            vertical_diffusion,
            drag,
            spectral_transform,
            hole_filling,
            large_scale_condensation,
            convection,
        )
    elseif coupled_ocean
        surface_heat_flux = SurfaceHeatFlux(
            spectral_grid;
            ocean = PrescribedOceanHeatFlux(spectral_grid),
            land = SurfaceLandHeatFlux(spectral_grid),
        )
        surface_humidity_flux = SurfaceHumidityFlux(
            spectral_grid;
            ocean = PrescribedOceanHumidityFlux(spectral_grid),
            land = SurfaceLandHumidityFlux(spectral_grid),
        )
        PrimitiveWetModel(
            spectral_grid;
            ocean,
            sea_ice,
            initial_conditions,
            land_sea_mask,
            surface_heat_flux,
            surface_humidity_flux,
            shortwave_radiation,
            longwave_radiation,
            time_stepping,
            horizontal_diffusion,
            vertical_diffusion,
            drag,
            spectral_transform,
            hole_filling,
            large_scale_condensation,
            convection,
        )
    else
        PrimitiveWetModel(
            spectral_grid;
            ocean,
            sea_ice,
            initial_conditions,
            land_sea_mask,
            shortwave_radiation,
            longwave_radiation,
            time_stepping,
            horizontal_diffusion,
            vertical_diffusion,
            drag,
            spectral_transform,
            hole_filling,
            large_scale_condensation,
            convection,
        )
    end
    add!(
        model.callbacks,
        :rrtmgp_postphysics_update =>
            RRTMGPPostPhysicsUpdateCallback(; radiation = longwave_radiation),
        :global_surface_temperature => GlobalSurfaceTemperatureCallback(spectral_grid),
        :global_radiation_budget => GlobalRadiationBudgetCallback(spectral_grid),
        :global_atmosphere_diagnostics =>
            GlobalAtmosphereDiagnosticsCallback(model),
    )
    return model
end

"""Build the concentration-dependent RRTMGP atmosphere with slab ocean/ice."""
function build_rrtmgp(config::ExperimentConfig)
    config.ocean_model == :slab || throw(
        ArgumentError("build_rrtmgp requires ocean_model=slab; use build_dynamic_esm otherwise"),
    )
    return initialize!(_rrtmgp_atmosphere_model(config))
end

function run_rrtmgp!(config::ExperimentConfig; steps::Union{Nothing, Int} = nothing)
    simulation = build_rrtmgp(config)
    if isnothing(steps)
        milliseconds = round(Int, config.duration_days * 86_400_000)
        SpeedyWeather.run!(simulation; period = Millisecond(milliseconds))
    else
        steps > 0 || throw(ArgumentError("steps must be positive"))
        SpeedyWeather.run!(simulation; steps)
    end
    return simulation, collect_diagnostics(simulation, config)
end
