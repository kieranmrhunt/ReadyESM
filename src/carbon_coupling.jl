const _CARBON_TRACER_NAMES = (:T, :S, :DIC, :TA)
const _CARBON_PGC_PER_PPM = 2.124
const _CARBON_DIC_STORAGE_OFFSET_UMOL_KG = 2000.0
const _CARBON_TA_STORAGE_OFFSET_UMOL_KG = 2300.0

@inline _carbon_storage_offset(::Symbol, ::Val, ::Type{FT}) where FT = zero(FT)

@inline function _carbon_storage_offset(
    mode::Symbol,
    ::Val{:DIC},
    ::Type{FT},
) where FT
    return mode in (:prescribed, :emissions_driven) ? FT(_CARBON_DIC_STORAGE_OFFSET_UMOL_KG) :
        zero(FT)
end

@inline function _carbon_storage_offset(
    mode::Symbol,
    ::Val{:TA},
    ::Type{FT},
) where FT
    return mode in (:prescribed, :emissions_driven) ? FT(_CARBON_TA_STORAGE_OFFSET_UMOL_KG) :
        zero(FT)
end

@inline _carbon_absolute_value(stored, offset) = stored + offset

"""Return active-only ocean constructor keywords for the carbon path."""
function _carbon_ocean_kwargs(config::ExperimentConfig)
    config.carbon.mode == :off && return NamedTuple()
    config.carbon.mode in (:tracer_only, :prescribed, :emissions_driven) || error(
        "unsupported carbon mode $(config.carbon.mode)",
    )
    return (; tracers = _CARBON_TRACER_NAMES)
end

"""
Device-adaptable freshwater exchange for an active carbon tracer.

NumericalEarth gives salinity an explicit carried-tracer flux that cancels
z-star's ambient carry when freshwater with zero salinity enters the ocean.
DIC and TA need the same correction. The carrying field is rebound to the
constructor-owned salinity freshwater field immediately after construction.
The carrying field and storage offset follow the model device through Adapt.
"""
mutable struct _CarbonFreshwaterExchange{name, F, O} <: Function
    carrying_flux::F
    storage_offset::O
end

"""Immutable kernel view of the host-rebindable carbon boundary."""
struct _CarbonFreshwaterKernelExchange{name, F, O} <: Function
    carrying_flux::F
    storage_offset::O
end

function Adapt.adapt_structure(
    to,
    exchange::Union{_CarbonFreshwaterExchange{name},
                    _CarbonFreshwaterKernelExchange{name}},
) where name
    carrying_flux = Adapt.adapt(to, exchange.carrying_flux)
    storage_offset = Adapt.adapt(to, exchange.storage_offset)
    # Keep host fields rebindable. A device Field is an isbits view, so its
    # enclosing callable must also be immutable when passed to a kernel.
    Exchange = isbitstype(typeof(carrying_flux)) &&
               isbitstype(typeof(storage_offset)) ?
        _CarbonFreshwaterKernelExchange : _CarbonFreshwaterExchange
    return Exchange{
        name,
        typeof(carrying_flux),
        typeof(storage_offset),
    }(
        carrying_flux,
        storage_offset,
    )
end

@inline function (exchange::Union{_CarbonFreshwaterExchange{name},
                                  _CarbonFreshwaterKernelExchange{name}})(
    i,
    j,
    grid,
    clock,
    fields,
) where name
    stored = @inbounds getproperty(fields, name)[i, j, grid.Nz]
    concentration = _carbon_absolute_value(stored, exchange.storage_offset)
    freshwater_flux = @inbounds exchange.carrying_flux[i, j, 1]
    return concentration * freshwater_flux
end

function _carbon_surface_boundary_condition(
    grid,
    eta_grid,
    name::Symbol,
    storage_offset,
)
    FT = eltype(grid)
    net_flux = Oceananigans.Field{
        Oceananigans.Center,
        Oceananigans.Center,
        Nothing,
    }(grid)
    dummy_freshwater = Oceananigans.Field{
        Oceananigans.Center,
        Oceananigans.Center,
        Nothing,
    }(eta_grid)
    fill!(net_flux, zero(FT))
    fill!(dummy_freshwater, zero(FT))
    exchange = _CarbonFreshwaterExchange{
        name,
        typeof(dummy_freshwater),
        typeof(storage_offset),
    }(dummy_freshwater, storage_offset)
    multiple = NumericalEarth.Oceans.MultipleFluxes(net_flux, exchange)
    top = Oceananigans.FluxBoundaryCondition(
        multiple;
        discrete_form = true,
    )
    return Oceananigans.FieldBoundaryConditions(; top)
end

function _carbon_surface_ocean_kwargs(
    config::ExperimentConfig,
    grid,
    supplied_free_surface = nothing,
)
    config.carbon.mode == :off && return NamedTuple()
    free_surface = isnothing(supplied_free_surface) ?
        NumericalEarth.Oceans.default_free_surface(grid) :
        supplied_free_surface
    eta_grid = if free_surface isa Oceananigans.SplitExplicitFreeSurface
        TX, TY, _ = Oceananigans.Grids.topology(grid)
        NumericalEarth.Oceans.maybe_extend_halos(
            TX,
            TY,
            grid,
            free_surface.substepping,
        )
    else
        grid
    end
    FT = eltype(grid)
    boundary_conditions = (
        DIC = _carbon_surface_boundary_condition(
            grid,
            eta_grid,
            :DIC,
            _carbon_storage_offset(config.carbon.mode, Val(:DIC), FT),
        ),
        TA = _carbon_surface_boundary_condition(
            grid,
            eta_grid,
            :TA,
            _carbon_storage_offset(config.carbon.mode, Val(:TA), FT),
        ),
    )
    kwargs = (; boundary_conditions)
    return isnothing(supplied_free_surface) ?
        merge((; free_surface), kwargs) : kwargs
end

function _carbon_freshwater_exchange(field, name)
    condition = field.boundary_conditions.top.condition
    hasproperty(condition, :func) || error(
        "$name carbon boundary is not a discrete function",
    )
    multiple = condition.func
    multiple isa NumericalEarth.Oceans.MultipleFluxes || error(
        "$name carbon boundary is not a multiple flux",
    )
    exchange = multiple.additional_fluxes
    exchange isa _CarbonFreshwaterExchange || error(
        "$name carbon boundary lacks its freshwater exchange",
    )
    return exchange
end

function _bind_carbon_freshwater_exchange!(ocean, config::ExperimentConfig)
    config.carbon.mode == :off && return nothing
    tracers = ocean.model.tracers
    freshwater = NumericalEarth.Oceans.extract_freshwater_flux(
        tracers.S.boundary_conditions.top.condition,
    )
    for name in (:DIC, :TA)
        exchange = _carbon_freshwater_exchange(getproperty(tracers, name), name)
        typeof(exchange.carrying_flux) == typeof(freshwater) || error(
            "$name and salinity freshwater fields have different types",
        )
        exchange.carrying_flux = freshwater
    end
    return freshwater
end

"""Carbonate constants from the BS10 v1/v2 families with corrected KW."""
function _carbon_carbonate_constants(temperature_kelvin::Real, salinity::Real)
    T = Float64(temperature_kelvin)
    S = Float64(salinity)
    isfinite(T) && T > 0 || throw(ArgumentError(
        "carbonate temperature must be finite and positive",
    ))
    isfinite(S) && S >= 0 || throw(ArgumentError(
        "carbonate salinity must be finite and nonnegative",
    ))

    T100 = T / 100
    lnK0 = -60.2409 + 93.4517 / T100 + 23.3585 * log(T100) +
           S * (0.023517 - 0.023656 * T100 + 0.0047036 * T100^2)
    pK1 = 3633.86 / T - 61.2172 + 9.6777 * log(T) -
          0.011555 * S + 0.0001152 * S^2
    pK2 = 471.78 / T + 25.9290 - 3.16967 * log(T) -
          0.01781 * S + 0.0001122 * S^2
    sqrt_salinity = sqrt(S)
    lnKB = (-8966.90 - 2890.53 * sqrt_salinity - 77.942 * S +
            1.728 * S * sqrt_salinity - 0.0996 * S^2) / T +
           148.0248 + 137.1942 * sqrt_salinity + 1.62142 * S -
           (24.4344 + 25.085 * sqrt_salinity + 0.2474 * S) * log(T) +
           0.053105 * sqrt_salinity * T
    # Dickson et al. (2007), chapter 5, applies the published Millero (1995)
    # erratum by subtracting 0.015 from the constant term.
    lnKw = 148.9652 - 13847.26 / T - 23.6521 * log(T) +
           (-5.977 + 118.67 / T + 1.0495 * log(T)) * sqrt_salinity -
           0.01615 * S
    total_boron = 4.157e-4 * S / 35
    return (
        K0 = exp(lnK0),
        K1 = 10.0^(-pK1),
        K2 = 10.0^(-pK2),
        KB = exp(lnKB),
        Kw = exp(lnKw),
        total_boron,
    )
end

function _carbon_alkalinity_residual(
    hydrogen,
    dissolved_co2,
    total_alkalinity,
    constants,
)
    carbonate = dissolved_co2 * (
        constants.K1 / hydrogen +
        2 * constants.K1 * constants.K2 / hydrogen^2
    )
    borate = constants.total_boron * constants.KB /
             (constants.KB + hydrogen)
    water = constants.Kw / hydrogen - hydrogen
    return carbonate + borate + water - total_alkalinity
end

"""
Return DIC in micromoles per kilogram at equilibrium with a specified pCO2.

At fixed pCO2 the dissolved CO2 concentration follows directly from Henry's
law. A single hydrogen-ion bisection then satisfies total alkalinity. This is
the algebraic inverse of the passed BS10 v1 pCO2 calculation and avoids an
outer DIC bisection during full-field initialization.
"""
function _carbon_equilibrium_dic_umol_kg(
    co2_ppm::Real,
    alkalinity_umol_kg::Real,
    temperature_kelvin::Real,
    salinity::Real,
)
    ppm = Float64(co2_ppm)
    alkalinity = Float64(alkalinity_umol_kg) * 1e-6
    ppm > 0 && isfinite(ppm) || throw(ArgumentError(
        "equilibrium pCO2 must be finite and positive",
    ))
    alkalinity > 0 && isfinite(alkalinity) || throw(ArgumentError(
        "equilibrium alkalinity must be finite and positive",
    ))
    constants = _carbon_carbonate_constants(temperature_kelvin, salinity)
    dissolved_co2 = constants.K0 * ppm * 1e-6
    residual(h) = _carbon_alkalinity_residual(
        h,
        dissolved_co2,
        alkalinity,
        constants,
    )
    low = 1e-10
    high = 1e-6
    residual(low) > 0 && residual(high) < 0 || throw(ArgumentError(
        "equilibrium carbonate system is not bracketed over pH 6 to 10",
    ))
    for _ in 1:96
        midpoint = sqrt(low * high)
        if residual(midpoint) > 0
            low = midpoint
        else
            high = midpoint
        end
    end
    hydrogen = sqrt(low * high)
    abs(residual(hydrogen)) <= 1e-9 * alkalinity || error(
        "equilibrium carbonate solve did not converge",
    )
    dic = dissolved_co2 * (
        1 + constants.K1 / hydrogen +
        constants.K1 * constants.K2 / hydrogen^2
    ) * 1e6
    isfinite(dic) && dic > 0 || error(
        "equilibrium carbonate solve returned invalid DIC",
    )
    return dic
end

"""Return seawater pCO2 in microatmospheres from DIC and alkalinity."""
function _carbon_pco2_uatm(
    dic_umol_kg::Real,
    alkalinity_umol_kg::Real,
    temperature_kelvin::Real,
    salinity::Real,
)
    dic = Float64(dic_umol_kg) * 1e-6
    alkalinity = Float64(alkalinity_umol_kg) * 1e-6
    dic > 0 && isfinite(dic) || throw(ArgumentError(
        "DIC must be finite and positive",
    ))
    alkalinity > 0 && isfinite(alkalinity) || throw(ArgumentError(
        "alkalinity must be finite and positive",
    ))
    constants = _carbon_carbonate_constants(temperature_kelvin, salinity)
    residual(h) = dic * (
        h * constants.K1 + 2 * constants.K1 * constants.K2
    ) / (
        h^2 + h * constants.K1 + constants.K1 * constants.K2
    ) + constants.total_boron * constants.KB / (constants.KB + h) +
        constants.Kw / h - h - alkalinity
    low = 1e-10
    high = 1e-6
    residual(low) > 0 && residual(high) < 0 || throw(ArgumentError(
        "carbonate system is not bracketed over pH 6 to 10",
    ))
    for _ in 1:96
        midpoint = sqrt(low * high)
        if residual(midpoint) > 0
            low = midpoint
        else
            high = midpoint
        end
    end
    hydrogen = sqrt(low * high)
    co2 = dic * hydrogen^2 / (
        hydrogen^2 + hydrogen * constants.K1 + constants.K1 * constants.K2
    )
    pco2 = co2 / constants.K0 * 1e6
    isfinite(pco2) && pco2 > 0 || error(
        "carbonate solve returned invalid pCO2",
    )
    return pco2
end

function _carbon_active_cell_mask(grid)
    cpu_grid = Oceananigans.on_architecture(Oceananigans.CPU(), grid)
    Nx, Ny, Nz = size(cpu_grid)
    active = falses(Nx, Ny, Nz)
    center = Oceananigans.Center()
    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        active[i, j, k] = !Oceananigans.Grids.inactive_node(
            i,
            j,
            k,
            cpu_grid,
            center,
            center,
            center,
        )
    end
    return active
end

"""Initialize native DIC and TA tracers after temperature and salinity."""
function _initialize_carbon_tracers!(ocean, config::ExperimentConfig)
    config.carbon.mode == :off && return nothing
    tracers = ocean.model.tracers
    names = propertynames(tracers)
    length(names) >= 4 && names[1:4] == _CARBON_TRACER_NAMES || error(
        "carbon tracer order differs from T, S, DIC, TA: $names",
    )

    temperature = Float64.(Array(Oceananigans.interior(tracers.T)))
    salinity = Float64.(Array(Oceananigans.interior(tracers.S)))
    active = _carbon_active_cell_mask(ocean.model.grid)
    size(temperature) == size(salinity) == size(active) || throw(
        DimensionMismatch("carbon initializer fields do not share one grid"),
    )
    FT = eltype(ocean.model.grid)
    dic = zeros(FT, size(active))
    alkalinity = zeros(FT, size(active))
    floor_salinity = config.carbon.chemistry_salinity_floor
    alkalinity_at_35 = config.carbon.alkalinity_at_salinity_35_umol_kg
    equilibrium_ppm = config.carbon.ocean_equilibrium_co2_ppm
    dic_offset = _carbon_storage_offset(
        config.carbon.mode,
        Val(:DIC),
        FT,
    )
    ta_offset = _carbon_storage_offset(
        config.carbon.mode,
        Val(:TA),
        FT,
    )
    for index in CartesianIndices(active)
        active[index] || continue
        temperature_celsius = temperature[index]
        salinity_value = salinity[index]
        isfinite(temperature_celsius) || error(
            "active ocean carbon cell has nonfinite temperature",
        )
        isfinite(salinity_value) || error(
            "active ocean carbon cell has nonfinite salinity",
        )
        chemistry_salinity = max(salinity_value, floor_salinity)
        ta = alkalinity_at_35 * chemistry_salinity / 35
        value = _carbon_equilibrium_dic_umol_kg(
            equilibrium_ppm,
            ta,
            temperature_celsius + 273.15,
            chemistry_salinity,
        )
        alkalinity[index] = FT(ta) - ta_offset
        dic[index] = FT(value) - dic_offset
    end
    absolute_dic = dic[active] .+ dic_offset
    absolute_alkalinity = alkalinity[active] .+ ta_offset
    all(isfinite, absolute_dic) && all(>(0), absolute_dic) || error(
        "active-cell DIC initialization is invalid",
    )
    all(isfinite, absolute_alkalinity) &&
        all(>(0), absolute_alkalinity) || error(
        "active-cell alkalinity initialization is invalid",
    )
    all(iszero, dic[.!active]) && all(iszero, alkalinity[.!active]) || error(
        "inactive carbon cells must have zero fill",
    )
    Oceananigans.set!(tracers.DIC, dic)
    Oceananigans.set!(tracers.TA, alkalinity)
    Oceananigans.fill_halo_regions!(tracers.DIC)
    Oceananigans.fill_halo_regions!(tracers.TA)
    Oceananigans.Architectures.synchronize(
        Oceananigans.Architectures.architecture(ocean.model.grid),
    )
    return nothing
end

function _carbon_ocean_inventory_pgc(model, mode::Symbol)
    ocean = model.ocean.model
    # on_architecture may reconstruct a grid from reference coordinates and
    # reset its evolving z-star metrics. Preserve the current coordinate state
    # when copying to the host for this inventory reduction.
    grid = Adapt.adapt(Array, ocean.grid)
    values = Float64.(Array(Oceananigans.interior(ocean.tracers.DIC)))
    size(values) == size(grid) || throw(DimensionMismatch(
        "carbon inventory field and ocean grid differ",
    ))
    density = Float64(model.interfaces.ocean_properties.reference_density)
    offset = Float64(_carbon_storage_offset(mode, Val(:DIC), Float64))
    center = Oceananigans.Center()
    total = 0.0
    compensation = 0.0
    for k in axes(values, 3), j in axes(values, 2), i in axes(values, 1)
        Oceananigans.Grids.inactive_node(
            i,
            j,
            k,
            grid,
            center,
            center,
            center,
        ) && continue
        concentration = values[i, j, k] + offset
        isfinite(concentration) && concentration > 0 || error(
            "active ocean DIC inventory value is invalid",
        )
        volume = Float64(Oceananigans.volume(
            i,
            j,
            k,
            grid,
            center,
            center,
            center,
        ))
        term = density * concentration * 1e-6 * volume *
               _CARBON_MOLAR_MASS_G_MOL * _CARBON_PG_PER_G
        next_total = total + term
        if abs(total) >= abs(term)
            compensation += (total - next_total) + term
        else
            compensation += (term - next_total) + total
        end
        total = next_total
    end
    inventory = total + compensation
    isfinite(inventory) && inventory > 0 || error(
        "initial ocean carbon inventory is invalid",
    )
    return inventory
end

mutable struct CarbonController{P}
    mode::Symbol
    atmospheric_carbon_pgc::Float64
    atmospheric_co2_ppm::Float64
    ocean_carbon_pgc::Float64
    cumulative_emissions_pgc::Float64
    cumulative_air_sea_uptake_pgc::Float64
    cumulative_implied_emissions_pgc::Float64
    previous_air_sea_uptake_rate_pgc_s::Float64
    ledger_time_seconds::Float64
    ledger_iteration::Int
    time_step_seconds::Float64
    provenance::P
end

const _CARBON_CONTROLLERS = WeakKeyDict{Any, Any}()
const _CARBON_CONTROLLERS_LOCK = ReentrantLock()

function _carbon_controller(model)
    return lock(_CARBON_CONTROLLERS_LOCK) do
        get(_CARBON_CONTROLLERS, model, nothing)
    end
end

function _register_carbon_controller!(model, controller)
    lock(_CARBON_CONTROLLERS_LOCK) do
        haskey(_CARBON_CONTROLLERS, model) && error(
            "a carbon controller is already installed on this model",
        )
        _CARBON_CONTROLLERS[model] = controller
    end
    return controller
end

function _carbon_provenance(config::ExperimentConfig)
    exchange = config.carbon.mode in (:prescribed, :emissions_driven)
    provenance = (
        schema = exchange ? "readiesm_carbon_stage0c_v1" :
            "readiesm_carbon_stage0a_v2",
        mode = String(config.carbon.mode),
        tracer_names = _CARBON_TRACER_NAMES,
        atmospheric_co2_ppm = config.carbon.atmospheric_co2_ppm,
        ocean_equilibrium_co2_ppm =
            config.carbon.ocean_equilibrium_co2_ppm,
        alkalinity_at_salinity_35_umol_kg =
            config.carbon.alkalinity_at_salinity_35_umol_kg,
        chemistry_salinity_floor = config.carbon.chemistry_salinity_floor,
        dic_storage_offset_umol_kg = exchange ?
            _CARBON_DIC_STORAGE_OFFSET_UMOL_KG : 0.0,
        ta_storage_offset_umol_kg = exchange ?
            _CARBON_TA_STORAGE_OFFSET_UMOL_KG : 0.0,
    )
    config.carbon.mode == :emissions_driven || return provenance
    return merge(provenance, (
        schema = "readiesm_carbon_emissions_v1",
        emissions_pgc_per_year = config.carbon.emissions_pgc_per_year,
        seconds_per_year = _CARBON_EMISSIONS_SECONDS_PER_YEAR,
    ))
end

function _apply_carbon_radiation_link!(model, controller::CarbonController)
    radiation = model.atmosphere.model.longwave_radiation
    radiation isa RRTMGPRadiation || error(
        "carbon controller requires RRTMGP radiation",
    )
    expected_carbon = controller.atmospheric_co2_ppm * _CARBON_PGC_PER_PPM
    isapprox(
        controller.atmospheric_carbon_pgc,
        expected_carbon;
        rtol = 8eps(Float64),
        atol = 0,
    ) || error(
        "carbon controller PgC and ppm differ",
    )
    radiation.co2_ppm = convert(
        typeof(radiation.co2_ppm),
        controller.atmospheric_co2_ppm,
    )
    return radiation.co2_ppm
end

function _install_carbon_controller!(
    model,
    config::ExperimentConfig,
    time_step_seconds,
)
    config.carbon.mode == :off && return nothing
    ppm = config.carbon.atmospheric_co2_ppm
    ocean_carbon_pgc = _carbon_ocean_inventory_pgc(
        model,
        config.carbon.mode,
    )
    controller = CarbonController(
        config.carbon.mode,
        ppm * _CARBON_PGC_PER_PPM,
        ppm,
        ocean_carbon_pgc,
        0.0,
        0.0,
        0.0,
        0.0,
        NaN,
        -1,
        Float64(time_step_seconds),
        _carbon_provenance(config),
    )
    controller.time_step_seconds > 0 || throw(ArgumentError(
        "carbon controller timestep must be positive",
    ))
    _register_carbon_controller!(model, controller)
    _initialize_carbon_diagnostic_history!(controller)
    _apply_carbon_radiation_link!(model, controller)
    return controller
end

function Oceananigans.prognostic_state(controller::CarbonController)
    state = (
        mode = String(controller.mode),
        atmospheric_carbon_pgc = controller.atmospheric_carbon_pgc,
        atmospheric_co2_ppm = controller.atmospheric_co2_ppm,
        ocean_carbon_pgc = controller.ocean_carbon_pgc,
        cumulative_emissions_pgc = controller.cumulative_emissions_pgc,
        cumulative_air_sea_uptake_pgc =
            controller.cumulative_air_sea_uptake_pgc,
        cumulative_implied_emissions_pgc =
            controller.cumulative_implied_emissions_pgc,
        previous_air_sea_uptake_rate_pgc_s =
            controller.previous_air_sea_uptake_rate_pgc_s,
        ledger_time_seconds = controller.ledger_time_seconds,
        ledger_iteration = controller.ledger_iteration,
        time_step_seconds = controller.time_step_seconds,
        provenance = deepcopy(controller.provenance),
    )
    return merge(state, _carbon_history_restart_extension(controller))
end

function _validated_carbon_checkpoint_values(
    controller::CarbonController,
    state,
)
    Symbol(state.mode) == controller.mode || error(
        "checkpoint carbon mode differs from the reconstructed controller",
    )
    isequal(state.provenance, controller.provenance) || error(
        "checkpoint carbon provenance differs from the reconstructed controller",
    )
    values = (
        atmospheric_carbon_pgc = Float64(state.atmospheric_carbon_pgc),
        atmospheric_co2_ppm = Float64(state.atmospheric_co2_ppm),
        ocean_carbon_pgc = Float64(state.ocean_carbon_pgc),
        cumulative_emissions_pgc = Float64(state.cumulative_emissions_pgc),
        cumulative_air_sea_uptake_pgc =
            Float64(state.cumulative_air_sea_uptake_pgc),
        cumulative_implied_emissions_pgc =
            Float64(state.cumulative_implied_emissions_pgc),
        previous_air_sea_uptake_rate_pgc_s =
            Float64(state.previous_air_sea_uptake_rate_pgc_s),
        ledger_time_seconds = Float64(state.ledger_time_seconds),
        ledger_iteration = Int(state.ledger_iteration),
        time_step_seconds = Float64(state.time_step_seconds),
    )
    finite_values = (
        values.atmospheric_carbon_pgc,
        values.atmospheric_co2_ppm,
        values.ocean_carbon_pgc,
        values.cumulative_emissions_pgc,
        values.cumulative_air_sea_uptake_pgc,
        values.cumulative_implied_emissions_pgc,
        values.previous_air_sea_uptake_rate_pgc_s,
        values.time_step_seconds,
    )
    all(isfinite, finite_values) || error(
        "checkpoint carbon state contains a nonfinite value",
    )
    values.atmospheric_carbon_pgc > 0 &&
        values.atmospheric_co2_ppm > 0 &&
        values.ocean_carbon_pgc > 0 || error(
        "checkpoint atmospheric carbon must be positive",
    )
    isapprox(
        values.atmospheric_carbon_pgc,
        values.atmospheric_co2_ppm * _CARBON_PGC_PER_PPM;
        rtol = 8eps(Float64),
        atol = 0,
    ) || error("checkpoint carbon PgC and ppm differ")
    values.time_step_seconds == controller.time_step_seconds || error(
        "checkpoint carbon timestep differs from the reconstructed controller",
    )
    if controller.mode == :emissions_driven
        values.cumulative_emissions_pgc >= 0 &&
            iszero(values.cumulative_implied_emissions_pgc) || error(
            "invalid emissions-driven carbon accounting in checkpoint",
        )
    end
    uninitialized = values.ledger_iteration == -1 &&
        isnan(values.ledger_time_seconds) &&
        iszero(values.previous_air_sea_uptake_rate_pgc_s)
    initialized = values.ledger_iteration >= 0 &&
        isfinite(values.ledger_time_seconds)
    uninitialized || initialized || error(
        "checkpoint carbon ledger position is invalid",
    )
    _validate_carbon_history_checkpoint(controller, state, values)
    return values
end

function Oceananigans.restore_prognostic_state!(
    controller::CarbonController,
    state,
)
    values = _validated_carbon_checkpoint_values(controller, state)
    for name in keys(values)
        setproperty!(controller, name, getproperty(values, name))
    end
    _restore_carbon_diagnostic_history!(controller, state)
    return controller
end

function _validate_carbon_restart_compatibility(model, state)
    controller = _carbon_controller(model)
    if isnothing(state)
        isnothing(controller) || error(
            "installed carbon controller is absent from the checkpoint",
        )
        return nothing
    end
    isnothing(controller) && error(
        "checkpoint contains carbon state but no controller is installed",
    )
    _validated_carbon_checkpoint_values(controller, state)
    return nothing
end

function _carbon_restart_state(model)
    controller = _carbon_controller(model)
    return isnothing(controller) ? nothing :
        Oceananigans.prognostic_state(controller)
end

function _restore_carbon_restart_state!(model, state)
    _validate_carbon_restart_compatibility(model, state)
    isnothing(state) && return nothing
    controller = _carbon_controller(model)
    Oceananigans.restore_prognostic_state!(controller, state)
    _apply_carbon_radiation_link!(model, controller)
    return controller
end

@kernel function _update_carbon_air_sea_boundary_kernel!(
    dic_boundary,
    dic,
    alkalinity,
    temperature,
    salinity,
    atmosphere_u,
    atmosphere_v,
    surface_pressure,
    ice_concentration,
    grid,
    atmospheric_co2_ppm,
    seawater_density,
    time_step_seconds,
    chemistry_salinity_floor,
    dic_storage_offset,
    ta_storage_offset,
)
    i, j = @index(Global, NTuple)
    k = size(grid, 3)
    inactive = Oceananigans.Grids.inactive_node(
        i,
        j,
        k,
        grid,
        Oceananigans.Center(),
        Oceananigans.Center(),
        Oceananigans.Center(),
    )
    @inbounds begin
        absolute_dic = _carbon_absolute_value(
            dic[i, j, k],
            dic_storage_offset,
        )
        absolute_ta = _carbon_absolute_value(
            alkalinity[i, j, k],
            ta_storage_offset,
        )
        temperature_kelvin = temperature[i, j, k] +
                             typeof(temperature[i, j, k])(273.15)
        chemistry_salinity = max(
            salinity[i, j, k],
            chemistry_salinity_floor,
        )
        wind_speed = sqrt(
            atmosphere_u[i, j, 1]^2 + atmosphere_v[i, j, 1]^2,
        )
        ice_fraction = ice_concentration[i, j, 1]
        pressure = surface_pressure[i, j, 1]
    end
    raw = _carbon_air_sea_flux_mol_m2_s(
        absolute_dic,
        absolute_ta,
        temperature_kelvin,
        chemistry_salinity,
        wind_speed,
        ice_fraction,
        atmospheric_co2_ppm,
        pressure,
        seawater_density,
    )
    horizontal_area = Oceananigans.Operators.Azᶜᶜᶜ(i, j, k, grid)
    cell_volume = Oceananigans.volume(
        i,
        j,
        k,
        grid,
        Oceananigans.Center(),
        Oceananigans.Center(),
        Oceananigans.Center(),
    )
    top_cell_thickness = ifelse(
        inactive,
        one(horizontal_area),
        cell_volume / horizontal_area,
    )
    limited = _carbon_apply_surface_cell_flux(
        absolute_dic,
        raw.flux,
        time_step_seconds,
        one(horizontal_area),
        top_cell_thickness,
        seawater_density,
    )
    valid = !inactive && raw.valid && limited.valid
    boundary_flux = -limited.flux_mol_m2_s *
                    typeof(limited.flux_mol_m2_s)(1e6) /
                    seawater_density
    @inbounds dic_boundary[i, j, 1] = ifelse(
        valid,
        boundary_flux,
        zero(boundary_flux),
    )
end

function _update_carbon_air_sea_boundary!(model)
    controller = _carbon_controller(model)
    isnothing(controller) && return nothing
    controller.mode in (:prescribed, :emissions_driven) || return nothing

    ocean = model.ocean.model
    grid = ocean.grid
    FT = eltype(grid)
    architecture = Oceananigans.Architectures.architecture(grid)
    launch_grid = grid isa
        Oceananigans.ImmersedBoundaries.ImmersedBoundaryGrid ?
        grid.underlying_grid : grid
    net_ocean = model.interfaces.net_fluxes.ocean
    atmosphere = model.interfaces.exchanger.atmosphere.state
    ice_concentration = model.sea_ice.model.ice_concentration
    provenance = controller.provenance
    Oceananigans.Utils.launch!(
        architecture,
        launch_grid,
        :xy,
        _update_carbon_air_sea_boundary_kernel!,
        net_ocean.DIC,
        ocean.tracers.DIC,
        ocean.tracers.TA,
        ocean.tracers.T,
        ocean.tracers.S,
        atmosphere.u,
        atmosphere.v,
        atmosphere.p,
        ice_concentration,
        grid,
        FT(controller.atmospheric_co2_ppm),
        FT(model.interfaces.ocean_properties.reference_density),
        FT(controller.time_step_seconds),
        FT(provenance.chemistry_salinity_floor),
        FT(provenance.dic_storage_offset_umol_kg),
        FT(provenance.ta_storage_offset_umol_kg),
    )
    fill!(net_ocean.TA, zero(FT))
    Oceananigans.fill_halo_regions!((net_ocean.DIC, net_ocean.TA))
    return net_ocean.DIC
end

struct CarbonLedgerCallback{F}
    boundary_integral::F
    reference_density::Float64
end

function CarbonLedgerCallback(model)
    controller = _carbon_controller(model)
    isnothing(controller) && error(
        "carbon ledger callback requires an installed controller",
    )
    controller.mode in (:prescribed, :emissions_driven) || error(
        "carbon ledger callback requires an air-sea exchange mode",
    )
    boundary = model.interfaces.net_fluxes.ocean.DIC
    return CarbonLedgerCallback(
        _budget_integral(boundary),
        Float64(model.interfaces.ocean_properties.reference_density),
    )
end

Oceananigans.initialize!(::CarbonLedgerCallback, simulation) = nothing

function _carbon_boundary_uptake_rate_pgc_s(callback::CarbonLedgerCallback)
    boundary_integral = _budget_scalar(callback.boundary_integral)
    return -boundary_integral * callback.reference_density * 1e-6 *
           _CARBON_MOLAR_MASS_G_MOL * _CARBON_PG_PER_G
end

function _advance_carbon_ledger!(
    controller::CarbonController,
    iteration::Int,
    time_seconds::Float64,
    current_rate::Float64,
)
    if controller.mode == :emissions_driven
        return _advance_emissions_carbon_ledger!(controller, iteration, time_seconds, current_rate)
    end
    controller.mode == :prescribed || error(
        "carbon ledger callback received mode $(controller.mode)",
    )
    isfinite(current_rate) || error(
        "carbon boundary integral is nonfinite",
    )

    if controller.ledger_iteration == -1
        controller.previous_air_sea_uptake_rate_pgc_s = current_rate
        controller.ledger_time_seconds = time_seconds
        controller.ledger_iteration = iteration
        return nothing
    end

    if iteration == controller.ledger_iteration
        time_seconds == controller.ledger_time_seconds || error(
            "carbon ledger time changed without an iteration change",
        )
        return nothing
    end

    iteration == controller.ledger_iteration + 1 || error(
        "carbon ledger callback skipped an iteration",
    )
    elapsed = time_seconds - controller.ledger_time_seconds
    elapsed > 0 || error("carbon ledger time did not advance")
    elapsed == controller.time_step_seconds || error(
        "carbon ledger timestep differs from the configured timestep",
    )
    uptake = controller.previous_air_sea_uptake_rate_pgc_s * elapsed
    previous_invariant = controller.atmospheric_carbon_pgc +
        controller.ocean_carbon_pgc - controller.cumulative_emissions_pgc -
        controller.cumulative_implied_emissions_pgc
    controller.ocean_carbon_pgc += uptake
    controller.cumulative_air_sea_uptake_pgc += uptake
    controller.cumulative_implied_emissions_pgc += uptake
    controller.ocean_carbon_pgc > 0 || error(
        "air-sea carbon flux exhausted the ocean inventory",
    )
    updated_invariant = controller.atmospheric_carbon_pgc +
        controller.ocean_carbon_pgc - controller.cumulative_emissions_pgc -
        controller.cumulative_implied_emissions_pgc
    isapprox(
        updated_invariant,
        previous_invariant;
        rtol = 64eps(Float64),
        atol = 0,
    ) || error("prescribed carbon ledger does not close")
    controller.previous_air_sea_uptake_rate_pgc_s = current_rate
    controller.ledger_time_seconds = time_seconds
    controller.ledger_iteration = iteration
    return nothing
end

function (callback::CarbonLedgerCallback)(simulation)
    model = simulation.model
    controller = _carbon_controller(model)
    isnothing(controller) && error("carbon controller disappeared")
    iteration = Int(model.clock.iteration)
    time_seconds = Float64(model.clock.time)
    current_rate = _carbon_boundary_uptake_rate_pgc_s(callback)
    previous_iteration = controller.ledger_iteration
    _advance_carbon_ledger!(
        controller,
        iteration,
        time_seconds,
        current_rate,
    )
    if controller.mode == :emissions_driven && iteration > 0 &&
        iteration != previous_iteration
        # The ocean just consumed the previous boundary. Use the updated
        # atmosphere for both the next gas-exchange boundary and radiation solve.
        _apply_carbon_radiation_link!(model, controller)
        _update_carbon_air_sea_boundary!(model)
        next_rate = _carbon_boundary_uptake_rate_pgc_s(callback)
        isfinite(next_rate) || error("updated carbon boundary integral is nonfinite")
        controller.previous_air_sea_uptake_rate_pgc_s = next_rate
    end
    _record_carbon_diagnostics!(controller, simulation)
    return nothing
end
