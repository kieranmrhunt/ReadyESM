# Rates are calculated and corrected in Float32, then exported in Float64
# kg m⁻² s⁻¹. Conversion does not recover precision lost in the native sums.
# Keep the old near-zero floor and allow four native relative roundoff units
# for the add/subtract corrections; this is not a water-budget tolerance.
function _precipitation_components_close(total, component_sum)
    axes(total) == axes(component_sum) || return false
    return all(eachindex(total, component_sum)) do i
        a, b = total[i], component_sum[i]
        isfinite(a) && isfinite(b) &&
            isapprox(a, b; atol = 1e-10, rtol = 4eps(Float32))
    end
end

function _validate_precipitation_phase_components(
    rain, convective_rain, large_scale_rain,
    snow, convective_snow, large_scale_snow,
)
    for component in (rain, convective_rain, large_scale_rain,
                      snow, convective_snow, large_scale_snow)
        axes(component) == axes(rain) || error("precipitation component shapes differ")
        all(x -> isfinite(x) && x >= 0, component) ||
            error("precipitation component is negative or nonfinite")
    end
    _precipitation_components_close(rain, convective_rain .+ large_scale_rain) ||
        error("rainfall flux does not equal convective plus large-scale components")
    _precipitation_components_close(snow, convective_snow .+ large_scale_snow) ||
        error("snowfall flux does not equal convective plus large-scale components")
    return nothing
end
