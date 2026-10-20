const _CARBON_HISTORY_SCHEMA = "readiesm_carbon_history_v1"
const _CARBON_HISTORY_INTERVAL_SECONDS = 86400.0

mutable struct CarbonDiagnosticHistory
    initial_total_pgc::Float64
    sample_every_n_steps::Int
    rows::Vector{NamedTuple}
end

const _CARBON_DIAGNOSTIC_ROW_KEYS = (
    :iteration, :time_seconds, :atmospheric_co2_ppm, :atmospheric_carbon_pgc,
    :ocean_carbon_pgc, :cumulative_emissions_pgc, :cumulative_air_sea_uptake_pgc,
    :cumulative_implied_emissions_pgc, :next_air_sea_uptake_rate_pgc_s,
    :ledger_residual_pgc,
)

const _CARBON_DIAGNOSTIC_HISTORIES = WeakKeyDict{Any, CarbonDiagnosticHistory}()
const _CARBON_DIAGNOSTIC_LOCK = ReentrantLock()

function _carbon_diagnostic_history(controller)
    return lock(_CARBON_DIAGNOSTIC_LOCK) do
        get(_CARBON_DIAGNOSTIC_HISTORIES, controller, nothing)
    end
end

function _initialize_carbon_diagnostic_history!(controller)
    controller.mode in (:prescribed, :emissions_driven) || return nothing
    initial_total = controller.atmospheric_carbon_pgc + controller.ocean_carbon_pgc -
        controller.cumulative_emissions_pgc - controller.cumulative_implied_emissions_pgc
    isfinite(initial_total) && initial_total > 0 || error("invalid initial carbon inventory")
    return lock(_CARBON_DIAGNOSTIC_LOCK) do
        haskey(_CARBON_DIAGNOSTIC_HISTORIES, controller) &&
            error("carbon diagnostic history is already installed")
        _CARBON_DIAGNOSTIC_HISTORIES[controller] =
            CarbonDiagnosticHistory(initial_total,
                max(1, round(Int, _CARBON_HISTORY_INTERVAL_SECONDS /
                                 controller.time_step_seconds)), NamedTuple[])
    end
end

function _carbon_diagnostic_row(controller, history)
    controller.ledger_iteration >= 0 || error("carbon ledger has not initialized")
    return (
        iteration = controller.ledger_iteration,
        time_seconds = controller.ledger_time_seconds,
        atmospheric_co2_ppm = controller.atmospheric_co2_ppm,
        atmospheric_carbon_pgc = controller.atmospheric_carbon_pgc,
        ocean_carbon_pgc = controller.ocean_carbon_pgc,
        cumulative_emissions_pgc = controller.cumulative_emissions_pgc,
        cumulative_air_sea_uptake_pgc = controller.cumulative_air_sea_uptake_pgc,
        cumulative_implied_emissions_pgc = controller.cumulative_implied_emissions_pgc,
        next_air_sea_uptake_rate_pgc_s = controller.previous_air_sea_uptake_rate_pgc_s,
        ledger_residual_pgc = controller.atmospheric_carbon_pgc +
            controller.ocean_carbon_pgc - controller.cumulative_emissions_pgc -
            controller.cumulative_implied_emissions_pgc - history.initial_total_pgc,
    )
end

function _record_carbon_diagnostics!(controller; sample_every_n_steps=nothing, at_end=false)
    history = _carbon_diagnostic_history(controller)
    isnothing(history) && return nothing
    cadence = isnothing(sample_every_n_steps) ?
        history.sample_every_n_steps : sample_every_n_steps
    cadence isa Integer && cadence > 0 || error("invalid carbon diagnostic cadence")
    row = _carbon_diagnostic_row(controller, history)
    if isempty(history.rows)
        row.iteration == 0 && row.time_seconds == 0 ||
            error("carbon diagnostic history must start at launch")
        history.sample_every_n_steps = cadence
        push!(history.rows, row)
    else
        cadence == history.sample_every_n_steps ||
            error("carbon diagnostic cadence changed after launch")
        previous = last(history.rows)
        row.iteration >= previous.iteration && row.time_seconds >= previous.time_seconds ||
            error("carbon diagnostic clock moved backwards")
        if row.iteration > previous.iteration &&
            (row.iteration % cadence == 0 || at_end)
            push!(history.rows, row)
        end
    end
    return nothing
end

function _record_carbon_diagnostics!(controller, simulation)
    isnothing(_carbon_diagnostic_history(controller)) && return nothing
    budget = simulation.callbacks[:global_coupled_budgets].func
    at_end = simulation.model.clock.iteration >= simulation.stop_iteration ||
        simulation.model.clock.time >= simulation.stop_time
    return _record_carbon_diagnostics!(controller;
        sample_every_n_steps=budget.sample_every_n_steps, at_end)
end

function _carbon_history_restart_extension(controller)
    history = _carbon_diagnostic_history(controller)
    isnothing(history) && return NamedTuple()
    return (; diagnostic_history=(
        schema = _CARBON_HISTORY_SCHEMA,
        sample_every_n_steps = history.sample_every_n_steps,
        interval_seconds = history.sample_every_n_steps * controller.time_step_seconds,
        initial_total_pgc = history.initial_total_pgc,
        rows = deepcopy(history.rows),
    ))
end

function _validate_carbon_history_checkpoint(controller, state, restored_values)
    history = _carbon_diagnostic_history(controller)
    present = hasproperty(state, :diagnostic_history)
    isnothing(history) && !present && return nothing
    isnothing(history) && error("checkpoint carbon history has no destination")
    present || error("checkpoint lacks carbon diagnostic history")
    saved = state.diagnostic_history
    saved.schema == _CARBON_HISTORY_SCHEMA || error("unknown carbon history schema")
    saved.sample_every_n_steps isa Integer && saved.sample_every_n_steps > 0 ||
        error("invalid checkpoint carbon diagnostic cadence")
    saved.interval_seconds == saved.sample_every_n_steps * restored_values.time_step_seconds ||
        error("carbon diagnostic interval differs")
    isfinite(saved.initial_total_pgc) && saved.initial_total_pgc > 0 ||
        error("invalid carbon history initial inventory")
    for (index, row) in enumerate(saved.rows)
        keys(row) == _CARBON_DIAGNOSTIC_ROW_KEYS || error("carbon history fields differ")
        row.iteration isa Integer && row.iteration >= 0 ||
            error("invalid carbon history iteration")
        all(isfinite, values(row)) || error("nonfinite carbon history")
        row.time_seconds == row.iteration * restored_values.time_step_seconds ||
            error("carbon history time and iteration differ")
        row.iteration <= restored_values.ledger_iteration ||
            error("carbon history extends beyond checkpoint")
        row.atmospheric_co2_ppm > 0 && row.atmospheric_carbon_pgc > 0 &&
            row.ocean_carbon_pgc > 0 || error("nonpositive carbon history inventory")
        isapprox(row.atmospheric_carbon_pgc,
            row.atmospheric_co2_ppm * _CARBON_PGC_PER_PPM;
            rtol=8eps(Float64), atol=0) || error("carbon history PgC and ppm differ")
        residual = row.atmospheric_carbon_pgc + row.ocean_carbon_pgc -
            row.cumulative_emissions_pgc - row.cumulative_implied_emissions_pgc -
            saved.initial_total_pgc
        row.ledger_residual_pgc == residual || error("carbon history ledger residual differs")
        if index == 1
            row.iteration == 0 || error("carbon history lacks launch sample")
        else
            previous = saved.rows[index - 1]
            row.iteration > previous.iteration ||
                error("carbon history iterations are not increasing")
            row.time_seconds > previous.time_seconds ||
                error("carbon history times are not increasing")
        end
    end
    restored_values.ledger_iteration < 0 || !isempty(saved.rows) ||
        error("initialized carbon checkpoint has empty diagnostic history")
    return nothing
end

function _restore_carbon_diagnostic_history!(controller, state)
    history = _carbon_diagnostic_history(controller)
    isnothing(history) && return nothing
    saved = state.diagnostic_history
    history.initial_total_pgc = saved.initial_total_pgc
    history.sample_every_n_steps = saved.sample_every_n_steps
    empty!(history.rows)
    append!(history.rows, deepcopy(saved.rows))
    return nothing
end

"""Copy daily history and append the current endpoint without sampling the model."""
function _carbon_output_rows(controller)
    history = _carbon_diagnostic_history(controller)
    isnothing(history) && error("carbon diagnostic history is absent")
    rows = deepcopy(history.rows)
    endpoint = _carbon_diagnostic_row(controller, history)
    if isempty(rows) || last(rows).iteration != endpoint.iteration
        push!(rows, endpoint)
    else
        isequal(last(rows), endpoint) || error("carbon history endpoint differs from ledger")
    end
    return rows
end

function _physical_carbon_output(values, wet, offset)
    size(values) == size(wet) || error("carbon output grid differs")
    physical = Float64.(Array(values)) .+ Float64(offset)
    physical[.!wet] .= NaN
    return physical
end

function _carbon_diagnostic_fields(model, config)
    config.carbon.mode == :off && return NamedTuple()
    controller = _carbon_controller(model)
    isnothing(controller) && error("active carbon model has no controller")
    controller.mode == config.carbon.mode || error("carbon diagnostic mode differs")
    ocean = model.ocean.model
    wet = _carbon_active_cell_mask(ocean.grid)
    dic = _physical_carbon_output(Oceananigans.interior(ocean.tracers.DIC), wet,
        _carbon_storage_offset(controller.mode, Val(:DIC), Float64))
    alkalinity = _physical_carbon_output(Oceananigans.interior(ocean.tracers.TA), wet,
        _carbon_storage_offset(controller.mode, Val(:TA), Float64))
    inventory = _carbon_ocean_inventory_pgc(model, controller.mode)
    rows = NamedTuple[]
    initial_total = nothing
    inventory_relative_error = nothing
    if controller.mode in (:prescribed, :emissions_driven)
        controller.ledger_iteration == model.clock.iteration &&
            controller.ledger_time_seconds == model.clock.time ||
            error("carbon diagnostics and model clock differ")
        rows = _carbon_output_rows(controller)
        initial_total = _carbon_diagnostic_history(controller).initial_total_pgc
        inventory_relative_error = (inventory - controller.ocean_carbon_pgc) /
            controller.ocean_carbon_pgc
    end
    return (; carbon_diagnostics=(
        schema = _CARBON_HISTORY_SCHEMA,
        mode = String(controller.mode),
        endpoint_time_seconds = Float64(model.clock.time),
        interval_seconds = isempty(rows) ? nothing :
            _carbon_diagnostic_history(controller).sample_every_n_steps * controller.time_step_seconds,
        initial_total_pgc = initial_total,
        ocean_tracer_inventory_pgc = inventory,
        ocean_inventory_relative_error = inventory_relative_error,
        reference_density_kg_m3 = Float64(model.interfaces.ocean_properties.reference_density),
        dic_storage_offset_umol_kg = _carbon_storage_offset(controller.mode, Val(:DIC), Float64),
        alkalinity_storage_offset_umol_kg = _carbon_storage_offset(controller.mode, Val(:TA), Float64),
        rows,
        wet,
        dissolved_inorganic_carbon_umol_kg = dic,
        alkalinity_umol_kg = alkalinity,
    ))
end

function _validate_carbon_diagnostics(diagnostics)
    hasproperty(diagnostics, :carbon_diagnostics) || return nothing
    carbon = diagnostics.carbon_diagnostics
    carbon.schema == _CARBON_HISTORY_SCHEMA || error("unknown carbon output schema")
    carbon.mode in ("tracer_only", "prescribed", "emissions_driven") ||
        error("unknown carbon output mode")
    isfinite(carbon.endpoint_time_seconds) && carbon.endpoint_time_seconds >= 0 ||
        error("invalid carbon output time")
    isfinite(carbon.ocean_tracer_inventory_pgc) && carbon.ocean_tracer_inventory_pgc > 0 ||
        error("invalid output carbon inventory")
    isfinite(carbon.reference_density_kg_m3) && carbon.reference_density_kg_m3 > 0 ||
        error("invalid carbon output reference density")
    all(isfinite, (carbon.dic_storage_offset_umol_kg, carbon.alkalinity_storage_offset_umol_kg)) ||
        error("invalid carbon output storage offset")
    for field in (carbon.dissolved_inorganic_carbon_umol_kg, carbon.alkalinity_umol_kg)
        size(field) == size(carbon.wet) || error("carbon map shape differs")
        all(isfinite, field[carbon.wet]) && all(>=(0), field[carbon.wet]) ||
            error("invalid active-cell carbon output")
        all(isnan, field[.!carbon.wet]) || error("dry carbon output is not masked")
    end
    if carbon.mode == "tracer_only"
        isempty(carbon.rows) && isnothing(carbon.initial_total_pgc) &&
            isnothing(carbon.ocean_inventory_relative_error) && isnothing(carbon.interval_seconds) ||
            error("tracer-only output has an exchange ledger")
    else
        !isempty(carbon.rows) || error("carbon output lacks ledger samples")
        isfinite(carbon.interval_seconds) && carbon.interval_seconds > 0 ||
            error("invalid carbon output interval")
        isfinite(carbon.initial_total_pgc) && carbon.initial_total_pgc > 0 ||
            error("invalid carbon output initial inventory")
        isfinite(carbon.ocean_inventory_relative_error) ||
            error("nonfinite carbon inventory difference")
        all(row->keys(row) == _CARBON_DIAGNOSTIC_ROW_KEYS &&
            all(isfinite, values(row)), carbon.rows) ||
            error("nonfinite carbon output history")
        first(carbon.rows).time_seconds == 0 || error("carbon output lacks launch time")
        last(carbon.rows).time_seconds == carbon.endpoint_time_seconds ||
            error("carbon output lacks current endpoint")
        times = [row.time_seconds for row in carbon.rows]
        all(>(0), diff(times)) || error("carbon output times are not increasing")
    end
    return nothing
end

function _write_carbon_diagnostics!(dataset, diagnostics, ocean_dimensions, time_origin)
    hasproperty(diagnostics, :carbon_diagnostics) || return nothing
    _validate_carbon_diagnostics(diagnostics)
    carbon = diagnostics.carbon_diagnostics
    dataset.attrib["carbon_diagnostics_schema"] = carbon.schema
    dataset.attrib["carbon_mode"] = carbon.mode
    dataset.attrib["carbon_endpoint_time_seconds"] = carbon.endpoint_time_seconds
    dataset.attrib["carbon_ocean_tracer_inventory_pgc"] = carbon.ocean_tracer_inventory_pgc
    dataset.attrib["carbon_reference_density_kg_m3"] = carbon.reference_density_kg_m3
    dataset.attrib["carbon_dic_storage_offset_umol_kg"] = carbon.dic_storage_offset_umol_kg
    dataset.attrib["carbon_alkalinity_storage_offset_umol_kg"] = carbon.alkalinity_storage_offset_umol_kg
    dataset.attrib["carbon_output_storage_convention"] =
        "physical concentrations with numerical storage offsets removed"
    dataset.attrib["carbon_inventory_geometry"] =
        "live native cell volumes, Float64 compensated fixed-order sum"
    for (name, field) in (
        ("ocean_dissolved_inorganic_carbon", carbon.dissolved_inorganic_carbon_umol_kg),
        ("ocean_total_alkalinity", carbon.alkalinity_umol_kg),
    )
        variable = defVar(dataset, name, Float64, (ocean_dimensions..., "ocean_layer");
                          attrib=Dict("units"=>"umol kg-1"))
        variable.attrib["long_name"] = name == "ocean_total_alkalinity" ?
            "Final total alkalinity, with dry cells masked" :
            "Final dissolved inorganic carbon, with storage offset removed and dry cells masked"
        variable[:, :, :] = field
    end
    isempty(carbon.rows) && return nothing
    dataset.attrib["carbon_initial_total_pgc"] = carbon.initial_total_pgc
    dataset.attrib["carbon_ocean_inventory_relative_error"] = carbon.ocean_inventory_relative_error
    dataset.attrib["carbon_history_interval_seconds"] = carbon.interval_seconds
    dataset.attrib["carbon_emissions_seconds_per_year"] = _CARBON_EMISSIONS_SECONDS_PER_YEAR
    defDim(dataset, "carbon_time", length(carbon.rows))
    time = defVar(dataset, "carbon_time", Float64, ("carbon_time",);
                  attrib=Dict("units"=>"seconds since $time_origin", "calendar"=>"proleptic_gregorian"))
    time[:] = [row.time_seconds for row in carbon.rows]
    iteration = defVar(dataset, "carbon_iteration", Int64, ("carbon_time",))
    iteration[:] = [row.iteration for row in carbon.rows]
    for (name, field, units) in (
        ("atmospheric_co2_ppm", :atmospheric_co2_ppm, "ppm"),
        ("atmospheric_carbon_pgc", :atmospheric_carbon_pgc, "Pg C"),
        ("ocean_carbon_ledger_pgc", :ocean_carbon_pgc, "Pg C"),
        ("cumulative_carbon_emissions_pgc", :cumulative_emissions_pgc, "Pg C"),
        ("cumulative_air_sea_carbon_uptake_pgc", :cumulative_air_sea_uptake_pgc, "Pg C"),
        ("cumulative_implied_carbon_emissions_pgc", :cumulative_implied_emissions_pgc, "Pg C"),
        ("next_air_sea_carbon_uptake_rate_pgc_s", :next_air_sea_uptake_rate_pgc_s, "Pg C s-1"),
        ("carbon_ledger_residual_pgc", :ledger_residual_pgc, "Pg C"),
    )
        variable = defVar(dataset, name, Float64, ("carbon_time",);
                          attrib=Dict("units"=>units))
        variable[:] = [getproperty(row, field) for row in carbon.rows]
        occursin("uptake", name) && (variable.attrib["positive"] = "into ocean")
    end
    mean_rate = zeros(Float64, length(carbon.rows))
    for i in 2:length(mean_rate)
        previous, current = carbon.rows[i - 1], carbon.rows[i]
        mean_rate[i] = (current.cumulative_air_sea_uptake_pgc -
            previous.cumulative_air_sea_uptake_pgc) /
            (current.time_seconds - previous.time_seconds)
    end
    mean_variable = defVar(dataset, "mean_air_sea_carbon_uptake_rate_pgc_s",
        Float64, ("carbon_time",);
        attrib=Dict("units"=>"Pg C s-1", "positive"=>"into ocean",
                    "long_name"=>"Mean uptake over the previous output interval, zero at launch"))
    mean_variable[:] = mean_rate
    dataset["next_air_sea_carbon_uptake_rate_pgc_s"].attrib["long_name"] =
        "Boundary uptake rate for the next coupling step"
    return nothing
end
