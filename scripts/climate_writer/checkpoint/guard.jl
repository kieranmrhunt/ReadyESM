module ManualStepGuard

export advance_checked!, callback_nan_detected

"""Read callback flags on the single-process CPU/GPU path without resetting them."""
callback_nan_detected(O, sim) = any(O.Diagnostics.nan_detected(cb.func) for cb in values(sim.callbacks))

"""Advance exactly one step, respecting stop requests and the full-run horizon.

An earlier segment boundary is chosen by the caller. It does not authorise
continuing after a native stop request or ignoring NaNs at the final step.
"""
function advance_checked!(step!, nan_detected, sim, iteration::Integer,
                          step_seconds::Integer, horizon_steps::Integer)
    step_seconds > 0 && 0 < iteration <= horizon_steps || error("invalid step request")
    sim.running || error("simulation already stopped before step $iteration")
    nan_detected(sim) && error("NaN flag set before step $iteration")
    sim.model.clock.iteration == iteration-1 || error("unexpected prior iteration")
    sim.model.clock.time == (iteration-1)*step_seconds || error("unexpected prior time")
    step!(sim)
    sim.model.clock.iteration == iteration || error("incomplete coupled step $iteration")
    sim.model.clock.time == iteration*step_seconds || error("unexpected completed-step time")
    nan_detected(sim) && error("NaN detected at step $iteration")
    sim.running || iteration == horizon_steps || error("simulation stopped early at step $iteration")
    return nothing
end

end
