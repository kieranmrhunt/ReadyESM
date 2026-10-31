module CheckpointSnapshot

"""Collect a current checkpoint diagnostic without adding to restart history.

The model's budget callback accumulates fluxes every step but samples states
less often. A segment checkpoint can fall between those samples. Append the
endpoint only to copied history vectors, then collect through a read-only
proxy. Ordinary diagnostic reduction caches may refresh during sampling.
The original callback, its histories and the model checkpoint are retained.
"""
function _collect_coupled_diagnostics(R, simulation, config)
    callback = simulation.callbacks[:global_coupled_budgets].func
    seconds = Float64(simulation.model.clock.time)
    day = seconds / 86_400
    isempty(callback.time_days) && error("Checkpoint has no initialized budget history")
    last(callback.time_days) <= day || error("Budget history is ahead of checkpoint")
    callback.last_time == seconds || error("Budget flux ledger is behind checkpoint")
    appended = last(callback.time_days) != day
    if appended
        C = typeof(callback)
        values = map(fieldnames(C)) do name
            value = getfield(callback, name)
            value isa Vector{<:Real} ? copy(value) : value
        end
        snapshot = C(values...)
        R._sample_coupled_budget_state!(snapshot, simulation)
        last(snapshot.time_days) == day || error("Snapshot sampler did not record checkpoint time")
        proxy = (; model=simulation.model,
                   callbacks=(; global_coupled_budgets=(; func=snapshot)))
        diagnostics = R.collect_dynamic_diagnostics(proxy, config)
    else
        diagnostics = R.collect_dynamic_diagnostics(simulation, config)
    end
    metadata = merge(diagnostics.metadata, (
        checkpoint_snapshot_time_days=day,
        checkpoint_budget_endpoint_sample=appended ?
            "output_only_copy_of_budget_history" : "existing_model_budget_sample",
    ))
    return merge(diagnostics, (; metadata))
end

"Collect with a private atmospheric endpoint, restoring original vector identities."
function collect_checkpoint_diagnostics(R, simulation, config)
    atmosphere = simulation.model.atmosphere
    model = atmosphere.model
    callback = model.callbacks[:global_atmosphere_diagnostics]
    seconds = Float64(simulation.model.clock.time)
    day = seconds / 86_400
    isempty(callback.time_days) && error("Checkpoint has no initialized atmosphere history")
    tolerance = 4 * eps(eltype(callback.time_days)) * max(1.0, abs(day))
    last(callback.time_days) <= day + tolerance || error("Atmosphere history is ahead of checkpoint")
    counter_seconds = callback.timestep_counter * Float64(model.time_stepping.Δt_sec)
    isapprox(counter_seconds, seconds; rtol=0, atol=1e-6) ||
        error("Atmosphere diagnostic counter differs from checkpoint time")
    appended = !isapprox(Float64(last(callback.time_days)), day; rtol=0, atol=tolerance)
    if appended
        C = typeof(callback)
        names = filter(name -> getfield(callback, name) isa Vector{<:Real}, fieldnames(C))
        originals = map(name -> getfield(callback, name), names)
        values = map(fieldnames(C)) do name
            value = getfield(callback, name)
            value isa Vector{<:Real} ? copy(value) : value
        end
        snapshot = C(values...)
        R._sample_global_atmosphere!(snapshot, atmosphere.variables, model)
        isapprox(Float64(last(snapshot.time_days)), day; rtol=0, atol=tolerance) ||
            error("Atmosphere sampler did not record checkpoint time")
        length(snapshot.time_days) == length(callback.time_days) + 1 ||
            error("Atmosphere sampler must append exactly one endpoint")
        # The collector accesses this concrete model's callback. Bind only the
        # private host vectors during synchronous collection, then restore their
        # original identities even if collection throws. Model/counter/ledger
        # fields and device accumulators are never replaced here.
        try
            for name in names
                setfield!(callback, name, getfield(snapshot, name))
            end
            diagnostics = _collect_coupled_diagnostics(R, simulation, config)
        finally
            for (name, value) in zip(names, originals)
                setfield!(callback, name, value)
            end
        end
    else
        diagnostics = _collect_coupled_diagnostics(R, simulation, config)
    end
    metadata = merge(diagnostics.metadata, (
        checkpoint_atmosphere_endpoint_sample=appended ?
            "output_only_copy_of_atmosphere_history" : "existing_model_atmosphere_sample",
    ))
    return merge(diagnostics, (; metadata))
end

end
