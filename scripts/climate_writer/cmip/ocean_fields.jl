module OceanCMIPFields

# Ocean fields that ClimateEval reads, as pure functions of host arrays (Claude, 2026-10-29, DESIGN.md).
# Arrays are (i, j, k) with k = 1 at the bottom, as Oceananigans stores them. Depths are positive downward from the
# moving free surface.

# Sums of these bands give ClimateEval's ranges (0-300 m, 0-2000 m, full column) and the common 0-700 m range.
const BANDS = (("0_300m", 0.0, 300.0), ("300_700m", 300.0, 700.0), ("700_2000m", 700.0, 2000.0),
               ("below_2000m", 2000.0, Inf))

const SPECS = let specs = Dict{String,Tuple{Symbol,String,String}}(
        "cmip_zos" => (:surface, "m", "sea-surface height (free-surface displacement); the global mean is not removed"),
        "cmip_mlotst" => (:surface, "m",
            "mixed-layer thickness, sigma-t criterion 0.125 kg m-3 from 10 m (NumericalEarth MixedLayerDepthField)"),
        "cmip_ocean_layer_heat_content" => (:volume, "J m-2",
            "rho_0 c_p theta dz per layer, from the current cell thickness at each step; zero in land cells"),
        "cmip_ocean_wet" => (:volume, "1", "wet cell mask, 1 ocean and 0 land, static (for regridding each level)"))
    for (band, top, bottom) in BANDS
        specs["cmip_ocean_heat_content_$band"] = (:surface, "J m-2",
            "rho_0 c_p theta integrated from $(top) to $(isinf(bottom) ? "the bottom" : "$(bottom) m") below the " *
            "free surface, with partial layers by thickness fraction")
    end
    specs
end

function register_metadata!(registry)
    isempty(intersect(keys(registry), keys(SPECS))) || error("ocean CMIP fields already registered")
    merge!(registry, SPECS)
    return nothing
end

function _check(theta, dz, depth_top, wet)
    size(theta) == size(dz) == size(depth_top) == size(wet) || throw(DimensionMismatch("ocean field shapes"))
    for (name, x) in (("theta", theta), ("dz", dz), ("depth_top", depth_top))
        all(isfinite(x[I]) for I in eachindex(x) if wet[I]) || error("nonfinite $name in a wet cell")
    end
    all(dz[I] > 0 for I in eachindex(dz) if wet[I]) || error("nonpositive wet-cell thickness")
    return nothing
end

"""Heat content per layer, rho_0 c_p theta dz (J m-2), zero in land cells."""
function layer_heat_content(theta, dz, wet, rho_cp)
    _check(theta, dz, zeros(size(theta)), wet)
    isfinite(rho_cp) && rho_cp > 0 || error("invalid rho_0 c_p")
    return [wet[I] ? Float64(rho_cp) * Float64(theta[I]) * Float64(dz[I]) : 0.0 for I in CartesianIndices(theta)]
end

"""
Heat content (J m-2) in each depth band. A layer that crosses a band edge contributes the fraction of its thickness
inside the band, with uniform temperature in the layer. `depth_top` is the depth of each cell's upper face below the
free surface.
"""
function band_heat_content(theta, dz, depth_top, wet, rho_cp)
    _check(theta, dz, depth_top, wet)
    isfinite(rho_cp) && rho_cp > 0 || error("invalid rho_0 c_p")
    nx, ny, nz = size(theta)
    result = Dict(band => zeros(nx, ny) for (band, _, _) in BANDS)
    for k in 1:nz, j in 1:ny, i in 1:nx
        wet[i, j, k] || continue
        top = Float64(depth_top[i, j, k])
        bottom = top + Float64(dz[i, j, k])
        heat_per_metre = Float64(rho_cp) * Float64(theta[i, j, k])
        for (band, upper, lower) in BANDS
            overlap = min(bottom, lower) - max(top, upper)
            overlap > 0 && (result[band][i, j] += heat_per_metre * overlap)
        end
    end
    return Dict("cmip_ocean_heat_content_$band" => values for (band, values) in result)
end

end
