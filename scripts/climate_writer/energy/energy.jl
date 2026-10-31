module ClimateEnergy

const SPECS = Dict(
    "open_water_fraction" => (:surface, "1", "open-water fraction of a unique wet ocean cell, zero elsewhere"),
    "sea_ice_surface_fraction" => (:surface, "1", "ice fraction of a unique wet ocean cell, zero elsewhere"),
    "open_water_latent_heat_contribution" => (:surface, "W m-2", "upward open-water latent heat per ocean-cell area, latest cached interface flux"),
    "sea_ice_latent_heat_contribution" => (:surface, "W m-2", "upward sea-ice latent heat per ocean-cell area, latest cached interface flux"),
    "ocean_ice_latent_heat_upward" => (:surface, "W m-2", "sum of upward open-water and sea-ice latent heat per ocean-cell area, latest cached interface flux"),
    "land_surface_fraction" => (:atmosphere, "1", "land fraction of an atmospheric cell"),
    "land_latent_heat_contribution" => (:atmosphere, "W m-2", "upward land latent heat per atmospheric-cell area, reconstructed from retained substep-average land humidity flux"),
    "open_water_vapor_contribution" => (:surface, "kg m-2 s-1", "upward open-water vapour mass per ocean-cell area, latest cached interface flux"),
    "sea_ice_vapor_contribution" => (:surface, "kg m-2 s-1", "upward sea-ice vapour mass per ocean-cell area, latest cached interface flux"),
)

function register_metadata!(registry)
    isempty(intersect(keys(registry), keys(SPECS))) || error("energy fields already registered")
    merge!(registry, SPECS)
    return nothing
end

"""Native ocean-cell energy contributions. Read only physically present tiles.

The caller supplies the unique wet-cell mask, including the tripolar fold
exclusion. Dry and duplicate storage may contain undefined component values.
Concentration must be physical on every included cell. Negative latent heat
(condensation/deposition) is valid and is preserved.
"""
function ocean_samples(open_flux, ice_flux, concentration, active)
    ndims(active) == 2 && !isempty(active) || error("empty or non-matrix ocean mask")
    all(size(x) == size(active) for x in (open_flux, ice_flux, concentration)) ||
        throw(DimensionMismatch("ocean latent heat shapes"))
    all(x -> x == 0 || x == 1, active) || error("ocean mask must be binary")
    water_fraction, ice_fraction = zeros(size(active)), zeros(size(active))
    water_heat, ice_heat = zeros(size(active)), zeros(size(active))
    for i in eachindex(active)
        active[i] == 0 && continue
        a = Float64(concentration[i])
        isfinite(a) && 0 <= a <= 1 || error("invalid active sea-ice fraction at $i")
        water_fraction[i], ice_fraction[i] = 1 - a, a
        if a < 1
            value = Float64(open_flux[i])
            isfinite(value) || error("nonfinite active open-water latent heat at $i")
            water_heat[i] = (1 - a) * value
        end
        if a > 0
            value = Float64(ice_flux[i])
            isfinite(value) || error("nonfinite active sea-ice latent heat at $i")
            ice_heat[i] = a * value
        end
    end
    return Dict("open_water_fraction" => water_fraction,
        "sea_ice_surface_fraction" => ice_fraction,
        "open_water_latent_heat_contribution" => water_heat,
        "sea_ice_latent_heat_contribution" => ice_heat,
        "ocean_ice_latent_heat_upward" => water_heat + ice_heat)
end

"Use the independently stored mass fluxes, without dividing energy by a latent heat."
function vapor_samples(open_flux, ice_flux, concentration, active)
    weighted = ocean_samples(open_flux, ice_flux, concentration, active)
    return Dict("open_water_vapor_contribution" => weighted["open_water_latent_heat_contribution"],
        "sea_ice_vapor_contribution" => weighted["sea_ice_latent_heat_contribution"])
end

"""Recover snow-free land energy from the retained coupling field.

ReadyESM divides its interval-average land energy by the Terrarium constant
before storing surface_humidity_flux. This inverse preserves that averaging
interval, subject to the storage precision. It is not the last auxiliary
land flux. Refuse snow-enabled land, which needs separate phase bookkeeping.
"""
function land_samples(humidity_flux, land_fraction, latent_heat; snow_enabled=false)
    snow_enabled && error("snow-enabled land requires a phase-specific energy observer")
    ndims(land_fraction) == 1 && !isempty(land_fraction) || error("empty or non-vector land fraction")
    size(humidity_flux) == size(land_fraction) || throw(DimensionMismatch("land latent heat shapes"))
    isfinite(latent_heat) && latent_heat > 0 || error("invalid land vaporisation heat")
    all(x -> isfinite(x) && 0 <= x <= 1, land_fraction) || error("invalid land fraction")
    result = zeros(length(land_fraction))
    for i in eachindex(land_fraction)
        f = Float64(land_fraction[i])
        iszero(f) && continue
        q = Float64(humidity_flux[i])
        isfinite(q) || error("nonfinite active land humidity flux at $i")
        result[i] = f * q * Float64(latent_heat)
    end
    return Dict("land_surface_fraction" => Float64.(land_fraction),
        "land_latent_heat_contribution" => result)
end

"""Sum native powers using a supplied physical Earth area.

Atmosphere weights sum to one. Ocean areas are unique wet-cell areas in m².
Report covered-area fraction instead of silently renormalizing away any
land/ocean mask mismatch. These cached/interval-mean snapshots are a climate
diagnostic, not a discrete whole-system energy-conservation residual.
"""
function global_summary(ocean, land, ocean_area, atmosphere_weights, earth_area)
    isfinite(earth_area) && earth_area > 0 || error("invalid Earth area")
    all(x -> isfinite(x) && x >= 0, ocean_area) || error("invalid ocean areas")
    all(x -> isfinite(x) && x > 0, atmosphere_weights) || error("invalid atmosphere weights")
    isapprox(sum(atmosphere_weights), 1; atol=2e-6, rtol=0) || error("atmosphere weights must sum to one")
    size(ocean_area) == size(ocean["ocean_ice_latent_heat_upward"]) || throw(DimensionMismatch("ocean area shape"))
    size(atmosphere_weights) == size(land["land_latent_heat_contribution"]) || throw(DimensionMismatch("land weight shape"))
    ocean_power = sum(Float64.(ocean_area) .* ocean["ocean_ice_latent_heat_upward"])
    land_mean = sum(Float64.(atmosphere_weights) .* land["land_latent_heat_contribution"])
    covered = sum(Float64.(ocean_area) .* (ocean["open_water_fraction"] + ocean["sea_ice_surface_fraction"])) / earth_area +
        sum(Float64.(atmosphere_weights) .* land["land_surface_fraction"])
    return (; ocean_power_w=ocean_power, land_power_w=land_mean * earth_area,
        global_latent_heat_upward=ocean_power / earth_area + land_mean,
        represented_surface_fraction=covered)
end

end
