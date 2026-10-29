# Oceananigans 0.110.15 sums the surface TKE flux for tuples of lengths 1–3.
# Longer tuples reach its generic zero-flux fallback. The production ocean
# uses four closures, or five with localized estuary mixing.
# Keep the existing short-tuple methods and recursively handle longer tuples.
# The tuple length and each closure type are known when the kernel compiles.
@inline function Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities.top_tke_flux(
    i, j, grid, clock, fields, parameters, closures::Tuple, buoyancy,
)
    isempty(closures) && return zero(grid)
    flux = Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities.top_tke_flux
    return flux(i, j, grid, clock, fields, parameters, first(closures), buoyancy) +
           flux(i, j, grid, clock, fields, parameters, Base.tail(closures), buoyancy)
end
