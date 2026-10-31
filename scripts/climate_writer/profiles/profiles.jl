module ClimateProfiles

# Reuse the pressure interpolation used by the released nine-field score.
include(joinpath(@__DIR__, "..", "shared", "monthly_core.jl"))

const TARGETS = (("temperature_200hpa", :temperature, 20000.0, "K"),
                 ("temperature_500hpa", :temperature, 50000.0, "K"),
                 ("temperature_850hpa", :temperature, 85000.0, "K"),
                 ("zonal_wind_200hpa", :wind, 20000.0, "m s-1"),
                 ("zonal_wind_850hpa", :wind, 85000.0, "m s-1"))
const SPECS = let specs = Dict(
    "column_water" => (:atmosphere, "kg m-2", "column water vapour from specific humidity and sigma-layer mass"),
    "surface_pressure" => (:atmosphere, "Pa", "surface air pressure"))
    for (name, _, pressure, units) in TARGETS
        specs[name * "_contribution"] = (:atmosphere, units,
            "$(name): valid log-pressure sample, zero if outside full-level range; mean requires valid_fraction = 1")
        specs[name * "_valid_fraction"] = (:atmosphere, "1",
            "fraction of samples with $(Int(pressure)) Pa inside the native full-level range")
    end
    specs
end

function register_metadata!(registry)
    isempty(intersect(keys(registry),keys(SPECS))) || error("climate profile fields already registered")
    merge!(registry,SPECS)
    return nothing
end

"""Finite contributions and coverage, suitable for the existing streaming accumulator.
No pressure extrapolation. Zero-filled missing samples must be paired with
their valid fractions. Use complete_mean for the release-score convention.
"""
function samples(temperature, humidity, u, ps, sigma, dsigma, gravity)
    n, nz = size(temperature)
    n > 0 && nz >= 2 || error("empty profile grid")
    size(humidity) == size(u) == (n,nz) && length(ps) == n ||
        throw(DimensionMismatch("profile field shapes"))
    length(sigma) == length(dsigma) == nz || throw(DimensionMismatch("sigma grid"))
    all(isfinite,sigma) && all(>(0),sigma) && all(<=(1),sigma) && all(>(0),diff(sigma)) ||
        error("sigma must increase toward the surface")
    all(isfinite,dsigma) && all(>(0),dsigma) || error("invalid sigma thickness")
    isapprox(sum(Float64.(dsigma)),1.0; atol=2e-6,rtol=0) || error("sigma layers must span the column")
    all(isfinite,temperature) && all(isfinite,humidity) && all(isfinite,u) &&
        all(isfinite,ps) && all(>(0),ps) || error("invalid profile state")
    isfinite(gravity) && gravity > 0 || error("invalid gravity")
    sig = Float64.(sigma)
    thick = Float64.(dsigma)
    result = Dict{String,Vector{Float64}}(
        "surface_pressure" => Float64.(ps),
        "column_water" => [sum(Float64(humidity[i,k])*thick[k] for k in 1:nz)*Float64(ps[i])/Float64(gravity) for i in 1:n])
    for (name, kind, pressure, _) in TARGETS
        field = kind == :temperature ? temperature : u
        contribution, valid = zeros(n), zeros(n)
        for i in 1:n
            value = MonthlyCore.pressure_temperature(view(field,i,:),sig,Float64(ps[i]),pressure)
            isnothing(value) && continue
            contribution[i] = value
            valid[i] = 1.0
        end
        result[name * "_contribution"] = contribution
        result[name * "_valid_fraction"] = valid
    end
    return result
end

"""Recover a monthly mean only where every sample covered the pressure level."""
function complete_mean(fields, name)
    name in first.(TARGETS) || throw(ArgumentError("unknown pressure field"))
    values, coverage = fields[name * "_contribution"], fields[name * "_valid_fraction"]
    size(values) == size(coverage) || throw(DimensionMismatch("coverage shape"))
    all(isfinite,values) && all(isfinite,coverage) && all(x -> 0 <= x <= 1,coverage) || error("invalid monthly profile")
    result = fill(NaN,size(values))
    valid = coverage .== 1.0
    result[valid] .= values[valid]
    return result
end

end
