# BS08 runoff-enthalpy diagnostics for normal runs (option on only).
#
# The history belongs to the BS08 exchange. It is held on the host, keyed by the exchange's mutable receiver guard,
# so it never enters a device kernel and the exact `_BS08Exchange` layout and dispatch stay unchanged. It is saved in
# the exchange checkpoint (schema v3). Sampling follows the coupled-budget cadence and endpoint rule. Collection is
# read-only.
#
# Timing. In coupled step n, the land step consumes pending runoff into the held enthalpy rate. The coupled state
# update at the end of step n scatters that rate to the ocean receivers, and the flux assembly applies the receiver
# correction. The ocean uses that correction during step n + 1. So:
#   released(n)  = sum over completed land intervals k <= n of  source enthalpy rate after step k  * dt_k
#   received(n)  = sum over completed ocean intervals k <= n of receiver rate after step k - 1     * dt_k
#   in transit   = receiver rate after step n (applied during step n + 1)
# A repeated assembly at the same clock leaves the receiver field bit-identical, and the callback integrates only
# when the iteration advances, so each interval is counted once.

const _BS08_HISTORY_SCHEMA = "readiesm_bs08_history_v1"
const _BS08_HISTORY_ROW_KEYS = (
    :iteration, :time_seconds,
    :source_water_kgs, :receiver_water_kgs, :source_enthalpy_w, :receiver_enthalpy_w,
    :native_sst_heat_w, :bs08_minus_native_heat_w, :zero_flow,
    :cumulative_released_enthalpy_j, :cumulative_received_enthalpy_j,
    :surface_excess_enthalpy_j, :pending_runoff_enthalpy_j, :cumulative_clamp_enthalpy_j,
    :partition_inconsistent_columns,
)

"Relative closure limit for routed water and enthalpy (the frozen routing gates use the same bound)."
const _BS08_ROUTING_CLOSURE_LIMIT = 2e-6

mutable struct BS08DiagnosticHistory
    source_area_m2::Vector{Float64}
    sample_every_n_steps::Int
    last_iteration::Int
    last_time_seconds::Float64
    in_effect_receiver_enthalpy_w::Float64
    cumulative_released_enthalpy_j::Float64
    cumulative_received_enthalpy_j::Float64
    rows::Vector{NamedTuple}
end

const _BS08_DIAGNOSTIC_HISTORIES = WeakKeyDict{BS08ReceiverGuard, BS08DiagnosticHistory}()
const _BS08_DIAGNOSTIC_LOCK = ReentrantLock()

_bs08_diagnostic_history(land::_BS08RunoffLand) =
    lock(() -> get(_BS08_DIAGNOSTIC_HISTORIES, land.runoff_exchange.guard, nothing), _BS08_DIAGNOSTIC_LOCK)
_bs08_diagnostic_history(::Any) = nothing

"Land-fraction-weighted source areas of the Terrarium columns, in exchange column order (m2)."
function _bs08_source_areas(atmosphere)
    spectral_grid = atmosphere.model.spectral_grid
    land_fraction = Float64.(Array(atmosphere.model.land_sea_mask.mask.data))
    points = findall(land_fraction .> 0)
    weights = Float64.(Array(_global_point_weights(spectral_grid)))
    return weights[points] .* (4π * Float64(atmosphere.model.planet.radius)^2) .* land_fraction[points]
end

function _install_bs08_diagnostic_history!(earth, sample_every_n_steps::Integer)
    land = earth.land
    land isa _BS08RunoffLand || return nothing
    sample_every_n_steps > 0 || throw(ArgumentError("BS08 diagnostic cadence must be positive"))
    area = _bs08_source_areas(earth.atmosphere)
    length(area) == length(land.runoff_exchange.held_enthalpy_rate) ||
        error("BS08 source areas and exchange columns differ")
    lock(_BS08_DIAGNOSTIC_LOCK) do
        haskey(_BS08_DIAGNOSTIC_HISTORIES, land.runoff_exchange.guard) &&
            error("BS08 diagnostic history is already installed")
        _BS08_DIAGNOSTIC_HISTORIES[land.runoff_exchange.guard] = BS08DiagnosticHistory(
            area, Int(sample_every_n_steps), -1, NaN, 0.0, 0.0, 0.0, NamedTuple[])
    end
    return nothing
end

_bs08_column_values(field) = Float64.(vec(Array(Oceananigans.interior(field))))

"""Routed water and enthalpy now: source totals against receiver integrals and the heat the correction adds.
Zero flow gives a missing (NaN) routed-water temperature and `zero_flow = true`."""
function _bs08_routing_budget(earth, area)
    land = earth.land
    ex = land.runoff_exchange
    ocean = earth.ocean.model
    J = earth.interfaces.exchanger.land.state.freshwater_flux
    H = ex.receiver_enthalpy_flux
    Nz = size(ocean.grid, 3)
    TJ = Oceananigans.Field{Oceananigans.Center, Oceananigans.Center, Nothing}(ocean.grid)
    Oceananigans.interior(TJ) .= Oceananigans.interior(J) .* view(Oceananigans.interior(ocean.tracers.T), :, :, Nz:Nz)
    cp = Float64(earth.interfaces.ocean_properties.heat_capacity)
    source_water = BS08_ROUTING_FRESHWATER_DENSITY * sum(area .* Float64.(vec(Array(land.surface_runoff))))
    source_enthalpy = sum(area .* Float64.(vec(Array(ex.held_enthalpy_rate))))
    receiver_water = _budget_scalar(_budget_integral(J))
    receiver_enthalpy = _budget_scalar(_budget_integral(H))
    native_sst_heat = cp * _budget_scalar(_budget_integral(TJ))
    zero_flow = iszero(source_water)
    return (; source_water_kgs = source_water, receiver_water_kgs = receiver_water,
            source_enthalpy_w = source_enthalpy, receiver_enthalpy_w = receiver_enthalpy,
            native_sst_heat_w = native_sst_heat, bs08_minus_native_heat_w = receiver_enthalpy - native_sst_heat,
            routed_water_mean_temperature_c = zero_flow ? NaN :
                source_enthalpy / (BS08_SURFACE_CV / 1000 * source_water),
            zero_flow)
end

"Routed-water mean temperature (degC) from a history row: missing (NaN) at zero flow, never 0 degC."
_bs08_routed_temperature(row) = row.zero_flow ? NaN : row.source_enthalpy_w / (BS08_SURFACE_CV / 1000 * row.source_water_kgs)

"Closure of a source/receiver pair: zero when both are zero, NaN (rejected) for an unmatched nonzero side."
function _bs08_relative_closure(source, receiver)
    iszero(source) && iszero(receiver) && return 0.0
    (iszero(source) || iszero(receiver)) && return NaN
    return abs(receiver - source) / abs(source)
end

"Validate one history row: field set, finite values, zero-flow flag and routed water and enthalpy closure."
function _validate_bs08_row(row)
    keys(row) == _BS08_HISTORY_ROW_KEYS || error("BS08 history fields differ")
    for key in _BS08_HISTORY_ROW_KEYS
        value = getproperty(row, key)
        value isa Bool || isfinite(value) || error("nonfinite BS08 history value $key")
    end
    row.iteration isa Integer && row.iteration >= 0 && row.time_seconds >= 0 || error("invalid BS08 history clock")
    row.zero_flow == iszero(row.source_water_kgs) || error("BS08 zero-flow flag disagrees with routed water")
    for (source, receiver, name) in ((row.source_water_kgs, row.receiver_water_kgs, "water"),
                                     (row.source_enthalpy_w, row.receiver_enthalpy_w, "enthalpy"))
        closure = _bs08_relative_closure(source, receiver)
        isnan(closure) && error("BS08 history has an unmatched nonzero $name side")
        closure <= _BS08_ROUTING_CLOSURE_LIMIT || error("BS08 routed $name closure $closure exceeds the limit")
    end
    return nothing
end

"Validate a row sequence: launch row first, strictly increasing iterations and times."
function _validate_bs08_rows(rows)
    isempty(rows) && error("BS08 history has no launch row")
    first(rows).iteration == 0 && first(rows).time_seconds == 0 || error("BS08 history lacks its launch row")
    for (index, row) in enumerate(rows)
        _validate_bs08_row(row)
        if index > 1
            previous = rows[index - 1]
            row.iteration > previous.iteration && row.time_seconds > previous.time_seconds ||
                error("BS08 history rows are not increasing")
        end
    end
    return nothing
end

"""
Columns whose last land step did not use the native runoff partition. The land model evaluates the partition in
each land step, so the launch row (iteration 0) comes before the first evaluation and reports 0.
"""
_bs08_partition_inconsistent(iteration, flag) = iteration == 0 ? 0 : count(!=(1.0), flag)

function _bs08_diagnostic_row(earth, history, iteration, time_seconds)
    budget = _bs08_routing_budget(earth, history.source_area_m2)
    terrarium = earth.atmosphere.variables.prognostic.land.terrarium
    area = history.source_area_m2
    surface = _bs08_column_values(terrarium.surface_excess_enthalpy)
    pending = _bs08_column_values(terrarium.runoff_pending_enthalpy)
    clamp = _bs08_column_values(terrarium.bs08_cumulative_clamp_enthalpy)
    flag = _bs08_column_values(terrarium.bs08_partition_consistent)
    length(surface) == length(pending) == length(clamp) == length(area) || error("BS08 column fields differ")
    return (; iteration, time_seconds, budget.source_water_kgs, budget.receiver_water_kgs,
            budget.source_enthalpy_w, budget.receiver_enthalpy_w, budget.native_sst_heat_w,
            budget.bs08_minus_native_heat_w, budget.zero_flow,
            cumulative_released_enthalpy_j = history.cumulative_released_enthalpy_j,
            cumulative_received_enthalpy_j = history.cumulative_received_enthalpy_j,
            surface_excess_enthalpy_j = sum(area .* surface), pending_runoff_enthalpy_j = sum(area .* pending),
            cumulative_clamp_enthalpy_j = sum(area .* clamp),
            partition_inconsistent_columns = _bs08_partition_inconsistent(iteration, flag))
end

"""Record one sample. Integrates the completed interval once per iteration advance, then appends a row at launch, on
the coupled-budget cadence and at the endpoint. A repeated call at the same iteration changes nothing."""
function _record_bs08_diagnostics!(simulation)
    earth = simulation.model
    history = _bs08_diagnostic_history(earth.land)
    if isnothing(history)
        earth.land isa _BS08RunoffLand && error("BS08 option is on but its diagnostic history is missing")
        return nothing
    end
    iteration = earth.clock.iteration
    time_seconds = Float64(earth.clock.time)
    iteration == history.last_iteration && return nothing
    if history.last_iteration < 0
        isempty(history.rows) || error("BS08 diagnostic history lost its clock")
        iteration == 0 || error("BS08 diagnostic history must start at launch; a checkpoint without BS08 history " *
                                "can be restored only at launch")
    else
        iteration > history.last_iteration && time_seconds > history.last_time_seconds ||
            error("BS08 diagnostic clock moved backwards")
        dt = time_seconds - history.last_time_seconds
        source_now = sum(history.source_area_m2 .* Float64.(vec(Array(earth.land.runoff_exchange.held_enthalpy_rate))))
        history.cumulative_released_enthalpy_j += source_now * dt
        history.cumulative_received_enthalpy_j += history.in_effect_receiver_enthalpy_w * dt
    end
    budget_callback = simulation.callbacks[:global_coupled_budgets].func
    budget_callback.sample_every_n_steps == history.sample_every_n_steps ||
        error("BS08 diagnostic cadence differs from the coupled-budget cadence")
    at_end = iteration >= simulation.stop_iteration || time_seconds >= simulation.stop_time
    # Every step needs only the receiver rate now in effect (one reduction). The full row (water integrals, the
    # native-SST term and the store integrals) is built only on sample steps.
    if iteration == 0 || iteration % history.sample_every_n_steps == 0 || at_end
        row = _bs08_diagnostic_row(earth, history, iteration, time_seconds)
        history.in_effect_receiver_enthalpy_w = row.receiver_enthalpy_w
        push!(history.rows, row)
    else
        history.in_effect_receiver_enthalpy_w =
            _budget_scalar(_budget_integral(earth.land.runoff_exchange.receiver_enthalpy_flux))
    end
    history.last_iteration = iteration
    history.last_time_seconds = time_seconds
    return nothing
end

"Coupled-simulation callback for BS08 diagnostics, added only when the option is on."
struct BS08DiagnosticsCallback end
(callback::BS08DiagnosticsCallback)(simulation) = _record_bs08_diagnostics!(simulation)
Oceananigans.initialize!(::BS08DiagnosticsCallback, simulation) = nothing

# ------------------------------------------------------------------------------------------------------------------
# Checkpoint extension (exchange schema v3 = v2 + diagnostic history)
# ------------------------------------------------------------------------------------------------------------------
function _bs08_history_checkpoint(land)
    history = _bs08_diagnostic_history(land)
    isnothing(history) && return nothing
    return (; schema = _BS08_HISTORY_SCHEMA, sample_every_n_steps = history.sample_every_n_steps,
            last_iteration = history.last_iteration, last_time_seconds = history.last_time_seconds,
            in_effect_receiver_enthalpy_w = history.in_effect_receiver_enthalpy_w,
            cumulative_released_enthalpy_j = history.cumulative_released_enthalpy_j,
            cumulative_received_enthalpy_j = history.cumulative_received_enthalpy_j,
            rows = deepcopy(history.rows))
end

"Validate a saved history completely before any destination changes."
function _validate_bs08_history_checkpoint(history, saved)
    saved.schema == _BS08_HISTORY_SCHEMA || error("unknown BS08 history schema")
    saved.sample_every_n_steps == history.sample_every_n_steps || error("BS08 history cadence differs")
    saved.last_iteration >= 0 || error("BS08 history checkpoint has no clock")
    all(isfinite, (saved.last_time_seconds, saved.in_effect_receiver_enthalpy_w,
                   saved.cumulative_released_enthalpy_j, saved.cumulative_received_enthalpy_j)) ||
        error("nonfinite BS08 history ledger")
    _validate_bs08_rows(saved.rows)
    endpoint = last(saved.rows)
    endpoint.iteration <= saved.last_iteration && endpoint.time_seconds <= saved.last_time_seconds ||
        error("BS08 history extends beyond its clock")
    if endpoint.iteration == saved.last_iteration
        endpoint.time_seconds == saved.last_time_seconds &&
            endpoint.cumulative_released_enthalpy_j == saved.cumulative_released_enthalpy_j &&
            endpoint.cumulative_received_enthalpy_j == saved.cumulative_received_enthalpy_j &&
            endpoint.receiver_enthalpy_w == saved.in_effect_receiver_enthalpy_w ||
            error("BS08 history endpoint row disagrees with the saved ledger")
    end
    return nothing
end

"Return an installed history to its launch state (used by restores without a saved history)."
function _reset_bs08_history!(history)
    history.last_iteration = -1
    history.last_time_seconds = NaN
    history.in_effect_receiver_enthalpy_w = 0.0
    history.cumulative_released_enthalpy_j = 0.0
    history.cumulative_received_enthalpy_j = 0.0
    empty!(history.rows)
    return history
end

function _restore_bs08_history!(history, saved)
    history.last_iteration = saved.last_iteration
    history.last_time_seconds = saved.last_time_seconds
    history.in_effect_receiver_enthalpy_w = saved.in_effect_receiver_enthalpy_w
    history.cumulative_released_enthalpy_j = saved.cumulative_released_enthalpy_j
    history.cumulative_received_enthalpy_j = saved.cumulative_received_enthalpy_j
    empty!(history.rows)
    append!(history.rows, deepcopy(saved.rows))
    return history
end

# ------------------------------------------------------------------------------------------------------------------
# Output (read-only collection) and validation
# ------------------------------------------------------------------------------------------------------------------
function _bs08_diagnostic_fields(earth)
    land = earth.land
    land isa _BS08RunoffLand || return NamedTuple()
    history = _bs08_diagnostic_history(land)
    isnothing(history) && error("BS08 option is on but its diagnostic history is missing")
    history.last_iteration == earth.clock.iteration || error("BS08 diagnostics and model clock differ")
    rows = deepcopy(history.rows)
    isempty(rows) && error("BS08 diagnostic history is empty")
    last(rows).iteration == earth.clock.iteration ||
        error("BS08 history lacks the current endpoint; sampling must follow the coupled-budget endpoint rule")
    terrarium = earth.atmosphere.variables.prognostic.land.terrarium
    return (; bs08_diagnostics = (
        schema = _BS08_HISTORY_SCHEMA,
        interval_seconds = history.sample_every_n_steps * Float64(earth.clock.last_Δt),
        rows,
        in_transit_enthalpy_j = history.in_effect_receiver_enthalpy_w * Float64(earth.clock.last_Δt),
        source_area_m2 = copy(history.source_area_m2),
        surface_excess_enthalpy_jm2 = _bs08_column_values(terrarium.surface_excess_enthalpy),
        pending_runoff_enthalpy_jm2 = _bs08_column_values(terrarium.runoff_pending_enthalpy),
        cumulative_clamp_enthalpy_jm2 = _bs08_column_values(terrarium.bs08_cumulative_clamp_enthalpy),
        receiver_enthalpy_flux_wm2 = Float64.(Array(Oceananigans.interior(land.runoff_exchange.receiver_enthalpy_flux,
                                                                          :, :, 1))),
    ))
end

function _validate_bs08_diagnostics(diagnostics)
    hasproperty(diagnostics, :bs08_diagnostics) || return nothing
    d = diagnostics.bs08_diagnostics
    d.schema == _BS08_HISTORY_SCHEMA || error("unknown BS08 output schema")
    _validate_bs08_rows(d.rows)
    all(row -> row.partition_inconsistent_columns == 0, d.rows) || error("BS08 host partition is not the native law")
    isfinite(d.in_transit_enthalpy_j) || error("nonfinite BS08 in-transit enthalpy")
    all(isfinite, d.surface_excess_enthalpy_jm2) && all(isfinite, d.pending_runoff_enthalpy_jm2) &&
        all(isfinite, d.cumulative_clamp_enthalpy_jm2) && all(isfinite, d.receiver_enthalpy_flux_wm2) ||
        error("nonfinite BS08 output map")
    return nothing
end

function _write_bs08_diagnostics!(dataset, diagnostics, ocean_dimensions, time_origin)
    hasproperty(diagnostics, :bs08_diagnostics) || return nothing
    _validate_bs08_diagnostics(diagnostics)
    d = diagnostics.bs08_diagnostics
    dataset.attrib["bs08_diagnostics_schema"] = d.schema
    dataset.attrib["bs08_history_interval_seconds"] = d.interval_seconds
    dataset.attrib["bs08_in_transit_enthalpy_j"] = d.in_transit_enthalpy_j
    dataset.attrib["bs08_enthalpy_reference"] = "liquid water at 0 degC"
    dataset.attrib["bs08_ledger_timing"] =
        "released: land interval of each completed step; received: ocean interval, one step later; " *
        "released minus received equals the in-transit enthalpy"
    defDim(dataset, "bs08_time", length(d.rows))
    time = defVar(dataset, "bs08_time", Float64, ("bs08_time",);
                  attrib = Dict("units" => "seconds since $time_origin", "calendar" => "proleptic_gregorian"))
    time[:] = [row.time_seconds for row in d.rows]
    iteration = defVar(dataset, "bs08_iteration", Int64, ("bs08_time",))
    iteration[:] = [row.iteration for row in d.rows]
    for (key, units, comment) in (
        (:source_water_kgs, "kg s-1", "routed water leaving land"),
        (:receiver_water_kgs, "kg s-1", "routed water reaching ocean receivers"),
        (:source_enthalpy_w, "W", "runoff enthalpy rate leaving land"),
        (:receiver_enthalpy_w, "W", "runoff enthalpy rate at ocean receivers, applied in the next ocean step"),
        (:native_sst_heat_w, "W", "native land-water heat term, referenced to sea surface temperature"),
        (:bs08_minus_native_heat_w, "W", "receiver enthalpy minus the native term"),
        (:cumulative_released_enthalpy_j, "J", "enthalpy released from land since launch"),
        (:cumulative_received_enthalpy_j, "J", "enthalpy applied to the ocean since launch"),
        (:surface_excess_enthalpy_j, "J", "area-integrated surface excess water enthalpy"),
        (:pending_runoff_enthalpy_j, "J", "area-integrated runoff enthalpy awaiting exchange"),
        (:cumulative_clamp_enthalpy_j, "J", "enthalpy given to water created by the bottom saturation clamp"),
    )
        name = startswith(String(key), "bs08_") ? String(key) : "bs08_" * String(key)
        variable = defVar(dataset, name, Float64, ("bs08_time",);
                          attrib = Dict("units" => units, "comment" => comment))
        variable[:] = [Float64(getproperty(row, key)) for row in d.rows]
    end
    temperature = defVar(dataset, "bs08_routed_water_mean_temperature_c", Float64, ("bs08_time",);
                         attrib = Dict("units" => "degC", "comment" => "missing (NaN) when bs08_zero_flow is 1"))
    temperature[:] = [_bs08_routed_temperature(row) for row in d.rows]
    zero_flow = defVar(dataset, "bs08_zero_flow", Int8, ("bs08_time",); attrib = Dict("units" => "1"))
    zero_flow[:] = [Int8(row.zero_flow) for row in d.rows]
    inconsistent = defVar(dataset, "bs08_partition_inconsistent_columns", Int64, ("bs08_time",))
    inconsistent[:] = [row.partition_inconsistent_columns for row in d.rows]
    for (name, values, units) in (
        ("bs08_source_area", d.source_area_m2, "m2"),
        ("bs08_surface_excess_enthalpy", d.surface_excess_enthalpy_jm2, "J m-2"),
        ("bs08_pending_runoff_enthalpy", d.pending_runoff_enthalpy_jm2, "J m-2"),
        ("bs08_cumulative_clamp_enthalpy", d.cumulative_clamp_enthalpy_jm2, "J m-2"),
    )
        variable = defVar(dataset, name, Float64, ("terrarium_land_column",); attrib = Dict("units" => units))
        variable[:] = values
    end
    flux = defVar(dataset, "bs08_receiver_enthalpy_flux", Float64, ocean_dimensions;
                  attrib = Dict("units" => "W m-2", "comment" => "final receiver enthalpy flux, positive into ocean"))
    flux[:, :] = d.receiver_enthalpy_flux_wm2
    return nothing
end
