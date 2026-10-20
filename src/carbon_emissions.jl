const _CARBON_EMISSIONS_SECONDS_PER_YEAR = 365.25 * 86_400

"""Advance atmosphere and ocean carbon using the boundary consumed this step."""
function _advance_emissions_carbon_ledger!(
    controller::CarbonController,
    iteration::Int,
    time_seconds::Float64,
    current_rate::Float64,
)
    controller.mode == :emissions_driven || error("expected emissions-driven controller")
    isfinite(current_rate) && isfinite(time_seconds) && iteration >= 0 ||
        error("invalid carbon ledger sample")
    rate = controller.provenance.emissions_pgc_per_year
    isfinite(rate) && rate >= 0 || error("invalid carbon emissions rate")
    if controller.ledger_iteration == -1
        iteration == 0 || error("emissions ledger must initialize at iteration zero")
        controller.previous_air_sea_uptake_rate_pgc_s = current_rate
        controller.ledger_time_seconds = time_seconds
        controller.ledger_iteration = iteration
        return nothing
    end
    if iteration == controller.ledger_iteration
        time_seconds == controller.ledger_time_seconds || error("carbon time changed without a step")
        return nothing
    end
    iteration == controller.ledger_iteration + 1 || error("carbon ledger skipped an iteration")
    elapsed = time_seconds - controller.ledger_time_seconds
    elapsed > 0 && elapsed == controller.time_step_seconds || error("carbon timestep changed")
    emissions = rate * (elapsed / _CARBON_EMISSIONS_SECONDS_PER_YEAR)
    uptake = controller.previous_air_sea_uptake_rate_pgc_s * elapsed
    before = _Stage0bCarbonLedger(controller.mode,
        controller.atmospheric_carbon_pgc, controller.ocean_carbon_pgc,
        controller.cumulative_emissions_pgc, controller.cumulative_air_sea_uptake_pgc,
        controller.cumulative_implied_emissions_pgc)
    iszero(before.cumulative_implied_emissions_pgc) || error("emissions mode has implied emissions")
    after = _carbon_ledger_step(before, uptake, emissions)
    ppm = after.atmospheric_carbon_pgc / _CARBON_PGC_PER_PPM
    all(isfinite, (after.atmospheric_carbon_pgc, after.ocean_carbon_pgc,
        after.cumulative_emissions_pgc, after.cumulative_air_sea_uptake_pgc, ppm)) ||
        error("nonfinite updated carbon state")
    isapprox(_carbon_ledger_invariant(after), _carbon_ledger_invariant(before);
        rtol=64eps(Float64), atol=0) || error("emissions-driven carbon ledger does not close")
    # Validate the complete update before changing the controller.
    controller.atmospheric_carbon_pgc = after.atmospheric_carbon_pgc
    controller.atmospheric_co2_ppm = ppm
    controller.ocean_carbon_pgc = after.ocean_carbon_pgc
    controller.cumulative_emissions_pgc = after.cumulative_emissions_pgc
    controller.cumulative_air_sea_uptake_pgc = after.cumulative_air_sea_uptake_pgc
    controller.previous_air_sea_uptake_rate_pgc_s = current_rate
    controller.ledger_time_seconds = time_seconds
    controller.ledger_iteration = iteration
    return nothing
end
