"""Apply the existing ocean T/S bounds to every saved coupled-budget sample."""
function _validate_ocean_extrema_history(time_days, temperature_minimum,
                                       temperature_maximum, salinity_minimum,
                                       salinity_maximum)
    isempty(time_days) && error("ocean extrema history is empty")
    series = (temperature_minimum, temperature_maximum,
              salinity_minimum, salinity_maximum)
    all(length(x) == length(time_days) for x in series) ||
        error("ocean extrema history and time dimensions differ")
    all(isfinite, time_days) && all(diff(time_days) .> 0) ||
        error("ocean extrema history times must be finite and strictly increasing")
    for i in eachindex(time_days)
        Tmin, Tmax, Smin, Smax = (x[i] for x in series)
        all(isfinite, (Tmin, Tmax, Smin, Smax)) ||
            error("ocean extrema history is non-finite at day $(time_days[i])")
        -5 <= Tmin <= Tmax <= 40 || error(
            "full-depth ocean temperature history is outside [-5, 40] degree_Celsius " *
            "at day $(time_days[i]): minimum=$Tmin, maximum=$Tmax",
        )
        0 <= Smin <= Smax <= 50 || error(
            "full-depth ocean salinity history is outside [0, 50] " *
            "at day $(time_days[i]): minimum=$Smin, maximum=$Smax",
        )
    end
    return nothing
end
