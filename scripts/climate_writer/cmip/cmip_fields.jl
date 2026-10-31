module CMIPFields

# Atmosphere fields that ClimateEval reads, as pure functions of host arrays (Claude, 2026-10-29, DESIGN.md).
# The pressure interpolation is the one the existing profile sampler and the released score use.
include(joinpath(@__DIR__, "..", "shared", "monthly_core.jl"))

# CMIP plev19 without 10, 5 and 1 hPa, which lie above the model's top full level.
const LEVELS_HPA = (1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20)
const PROFILE_FIELDS = (("ta", "K", "air temperature"), ("ua", "m s-1", "eastward wind"),
                        ("va", "m s-1", "northward wind"), ("hus", "1", "specific humidity"),
                        ("zg", "m", "geopotential height"))
const SURFACE_FIELDS = (
    ("cmip_tas", "K", "near-surface air temperature: lowest-level temperature lowered to the surface with the dry " *
        "adiabatic lapse rate (SpeedyWeather surface_air_temperature), not a 2 m similarity diagnostic"),
    ("cmip_ts", "K", "radiative surface temperature: ((1-f) T_sea^4 + f T_land^4)^(1/4), as passed to RRTMGP"),
    ("cmip_clt", "%", "total cloud cover, maximum-random overlap of the latest RRTMGP cloud fraction"),
    ("cmip_clwvi", "kg m-2", "grid-box mean cloud liquid plus ice water path, latest RRTMGP solve"),
    ("cmip_clivi", "kg m-2", "grid-box mean cloud ice water path, latest RRTMGP solve"),
    ("cmip_tauu", "Pa", "surface downward eastward stress applied by SpeedyWeather SurfaceMomentumFlux"),
    ("cmip_tauv", "Pa", "surface downward northward stress applied by SpeedyWeather SurfaceMomentumFlux"),
)

profile_name(field, level) = "cmip_$(field)_$(level)hpa"

const SPECS = let specs = Dict{String,Tuple{Symbol,String,String}}()
    for (name, units, label) in SURFACE_FIELDS
        specs[name] = (:atmosphere, units, label)
    end
    for (field, units, label) in PROFILE_FIELDS, level in LEVELS_HPA
        name = profile_name(field, level)
        specs[name * "_contribution"] = (:atmosphere, units,
            "$(label) at $(level) hPa: log-pressure sample, zero outside the full-level range; " *
            "mean requires valid_fraction = 1")
        specs[name * "_valid_fraction"] = (:atmosphere, "1",
            "fraction of samples with $(level) hPa inside the native full-level range")
    end
    specs
end

function register_metadata!(registry)
    isempty(intersect(keys(registry), keys(SPECS))) || error("CMIP fields already registered")
    merge!(registry, SPECS)
    return nothing
end

"""Total cover of one column with maximum-random overlap (Geleyn and Hollingsworth 1979): adjacent cloudy layers
overlap maximally, and layers separated by a clear layer overlap randomly. The order of `fraction` is the vertical
order (top-down or bottom-up give the same result)."""
function max_random_cover(fraction)
    clear = 1.0
    previous = 0.0
    for c in fraction
        c = clamp(Float64(c), 0.0, 1.0)
        previous >= 1 && return 1.0
        clear *= (1 - max(c, previous)) / (1 - previous)
        previous = c
    end
    return clamp(1 - clear, 0.0, 1.0)
end

"""Grid-box mean path (kg m-2) from per-layer cloud fraction and in-cloud path (g m-2)."""
function grid_mean_path(fraction, path_gm2)
    length(fraction) == length(path_gm2) || throw(DimensionMismatch("cloud layers"))
    return sum(clamp(Float64(c), 0.0, 1.0) * Float64(p) for (c, p) in zip(fraction, path_gm2)) / 1000
end

"""The RRTMGP lower-boundary temperature: emitted-flux mean of the sea (ocean and ice) and land temperatures."""
radiative_surface_temperature(land_fraction, t_sea, t_land) =
    ((1 - Float64(land_fraction)) * Float64(t_sea)^4 + Float64(land_fraction) * Float64(t_land)^4)^(1 / 4)

"""
Surface downward stress (Pa) exactly as SpeedyWeather's SurfaceMomentumFlux applies it: its upward flux into the
lowest layer is -rho C_D V_0 f u, so the downward stress on the surface is +rho C_D V_0 f u (CMIP sign: positive
when eastward momentum goes into the surface). `u` and `v` are the lagged lowest-level winds that the scheme uses.
`drag` is the coefficient that the scheme selects (the boundary-layer drag or the land-ocean mix).
"""
function surface_stress(density, drag, wind_speed, slowdown, u, v)
    n = length(density)
    all(length(x) == n for x in (drag, wind_speed, u, v)) || throw(DimensionMismatch("stress field lengths"))
    for x in (density, drag, wind_speed, u, v)
        all(isfinite, x) || error("nonfinite stress input")
    end
    isfinite(slowdown) && slowdown > 0 || error("invalid wind slowdown")
    factor = Float64.(density) .* Float64.(drag) .* Float64.(wind_speed) .* Float64(slowdown)
    return factor .* Float64.(u), factor .* Float64.(v)
end

"""Drag coefficient as SurfaceMomentumFlux selects it."""
selected_drag(use_boundary_layer_drag, boundary_layer_drag, land_fraction, drag_land, drag_ocean) =
    use_boundary_layer_drag ? Float64.(boundary_layer_drag) :
        Float64.(land_fraction) .* Float64(drag_land) .+ (1 .- Float64.(land_fraction)) .* Float64(drag_ocean)

"""
Surface fields at every atmospheric point. Cloud arrays are RRTMGP's (layer, column) arrays. `land_temperature` is
the top land temperature that the radiation reads.
"""
function surface_samples(surface_air_temperature, land_fraction, sea_temperature, land_temperature,
                         cloud_fraction, liquid_path, ice_path, tauu, tauv)
    n = length(surface_air_temperature)
    length(land_fraction) == length(sea_temperature) == length(land_temperature) == n ||
        throw(DimensionMismatch("surface field lengths"))
    size(cloud_fraction) == size(liquid_path) == size(ice_path) && size(cloud_fraction, 2) == n ||
        throw(DimensionMismatch("cloud arrays"))
    all(isfinite, surface_air_temperature) && all(>(0), surface_air_temperature) || error("invalid surface air temperature")
    all(f -> 0 <= f <= 1, land_fraction) || error("land fraction outside [0, 1]")
    all(isfinite, cloud_fraction) && all(isfinite, liquid_path) && all(isfinite, ice_path) || error("nonfinite clouds")
    ts = [radiative_surface_temperature(land_fraction[i], sea_temperature[i], land_temperature[i]) for i in 1:n]
    all(isfinite, ts) && all(>(0), ts) || error("invalid radiative surface temperature")
    clt = [100 * max_random_cover(view(cloud_fraction, :, i)) for i in 1:n]
    clivi = [grid_mean_path(view(cloud_fraction, :, i), view(ice_path, :, i)) for i in 1:n]
    clwvi = [grid_mean_path(view(cloud_fraction, :, i), view(liquid_path, :, i)) for i in 1:n] .+ clivi
    length(tauu) == length(tauv) == n && all(isfinite, tauu) && all(isfinite, tauv) || error("invalid surface stress")
    return Dict{String,Vector{Float64}}("cmip_tas" => Float64.(surface_air_temperature), "cmip_ts" => ts,
        "cmip_clt" => clt, "cmip_clwvi" => clwvi, "cmip_clivi" => clivi,
        "cmip_tauu" => Float64.(tauu), "cmip_tauv" => Float64.(tauv))
end

"""
Pressure-level contributions and coverage for ta, ua, va, hus and zg, with the same method and conventions as the
existing ClimateProfiles sampler. `geopotential` is in m2 s-2 on full levels and includes the surface geopotential.
"""
function profile_samples(temperature, u, v, humidity, geopotential, ps, sigma, gravity)
    n, nz = size(temperature)
    n > 0 && nz >= 2 || error("empty profile grid")
    all(size(x) == (n, nz) for x in (u, v, humidity, geopotential)) && length(ps) == n ||
        throw(DimensionMismatch("profile field shapes"))
    length(sigma) == nz || throw(DimensionMismatch("sigma grid"))
    all(isfinite, sigma) && all(>(0), sigma) && all(<=(1), sigma) && all(>(0), diff(sigma)) ||
        error("sigma must increase toward the surface")
    for x in (temperature, u, v, humidity, geopotential, ps)
        all(isfinite, x) || error("invalid profile state")
    end
    all(>(0), ps) || error("invalid surface pressure")
    isfinite(gravity) && gravity > 0 || error("invalid gravity")
    sig = Float64.(sigma)
    sources = Dict("ta" => temperature, "ua" => u, "va" => v, "hus" => humidity, "zg" => geopotential)
    scale = Dict("ta" => 1.0, "ua" => 1.0, "va" => 1.0, "hus" => 1.0, "zg" => 1 / Float64(gravity))
    result = Dict{String,Vector{Float64}}()
    for (field, _, _) in PROFILE_FIELDS, level in LEVELS_HPA
        target = 100.0 * level
        contribution, valid = zeros(n), zeros(n)
        values = sources[field]
        for i in 1:n
            value = MonthlyCore.pressure_temperature(view(values, i, :), sig, Float64(ps[i]), target)
            isnothing(value) && continue
            contribution[i] = scale[field] * value
            valid[i] = 1.0
        end
        name = profile_name(field, level)
        result[name * "_contribution"] = contribution
        result[name * "_valid_fraction"] = valid
    end
    return result
end

end
