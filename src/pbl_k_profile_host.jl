#=
Opt-in boundary-layer K-profile mixing in SpeedyWeather's vertical_diffusion slot (Stage 2 host wiring).

Selected by `atmosphere_vertical_diffusion = :holtslag_boville_k_profile_v1`. Production keeps the zero operator. The
column physics is in pbl_k_profile.jl. Host conventions, agreed with Codex's interface review:
- vertical_diffusion runs first in each column, before convection and the surface fluxes, and SpeedyWeather zeroes the
  top-level surface flux totals before column physics. The surface buoyancy flux therefore comes from the persistent
  prescribed coupled inputs, `prognostic.ocean` and `prognostic.land` sensible heat (W/m2) and humidity (kg/m2/s)
  fluxes, positive upward, conditional means weighted once by the land fraction (a zero-area tile is never read);
- the friction velocity is the stress equivalent of SpeedyWeather's surface momentum flux,
  u*^2 = C_D V0 f |V1|, with the lagged stored drag (or the fixed drag the momentum scheme would use), V0 and the
  surface air density computed locally in the SurfaceCondition form. The shared drag and stress are unchanged. At
  start-up C_D can be zero; then u* = 0 and the calm-safe closure applies (no neutral fallback);
- dry static energy s = cp T + Phi and humidity are mixed from the `_prev` state over 2 dt, and the flux-form tendency
  of the backward-Euler solution is added (the leapfrog update over 2 dt then equals the implicit step). Momentum is
  not mixed in v1. Heights above the surface are hydrostatic from the column's virtual temperature.
=#

const PBL_K_PROFILE_PROVENANCE = "holtslag_boville_1993_k_profile_backward_euler_sigma_flux_scalars_v1"
const PBL_VIRTUAL_FACTOR = 0.61

"""Holtslag-Boville K-profile mixing of dry static energy and humidity (opt-in vertical diffusion)."""
struct KProfileBoundaryLayerMixing{NF} <: SpeedyWeather.AbstractVerticalDiffusion
    critical_richardson::NF
    maximum_height::NF
end

Adapt.@adapt_structure KProfileBoundaryLayerMixing

KProfileBoundaryLayerMixing(spectral_grid::SpeedyWeather.SpectralGrid; critical_richardson = PBL_CRITICAL_RICHARDSON,
                            maximum_height = 5000) =
    KProfileBoundaryLayerMixing{spectral_grid.NF}(spectral_grid.NF(critical_richardson),
                                                  spectral_grid.NF(maximum_height))

_atmosphere_vertical_diffusion_provenance(::KProfileBoundaryLayerMixing) = PBL_K_PROFILE_PROVENANCE

SpeedyWeather.variables(::KProfileBoundaryLayerMixing) = (
    SpeedyWeather.ParameterizationVariable(:pbl_height, SpeedyWeather.Grid2D();
        desc = "K-profile boundary-layer depth", units = "m"),
    SpeedyWeather.ParameterizationVariable(:pbl_friction_velocity, SpeedyWeather.Grid2D();
        desc = "Stress-equivalent friction velocity used by the K-profile", units = "m/s"),
    SpeedyWeather.ParameterizationVariable(:pbl_buoyancy_flux, SpeedyWeather.Grid2D();
        desc = "Surface kinematic virtual heat flux used by the K-profile", units = "K m/s"),
    SpeedyWeather.ParameterizationVariable(:pbl_drag_coefficient, SpeedyWeather.Grid2D();
        desc = "Lagged momentum drag coefficient used by the K-profile", units = "1"),
    SpeedyWeather.ParameterizationVariable(:pbl_lowest_face_diffusivity, SpeedyWeather.Grid2D();
        desc = "K-profile eddy diffusivity at the top face of the lowest layer", units = "m^2/s"),
    SpeedyWeather.ParameterizationVariable(:pbl_height_work, SpeedyWeather.Grid3D();
        desc = "K-profile full-level heights above the surface (work column)", units = "m"),
    SpeedyWeather.ParameterizationVariable(:pbl_solution_work, SpeedyWeather.Grid3D();
        desc = "K-profile implicit solution (work column)", units = "?"),
)

SpeedyWeather.initialize!(::KProfileBoundaryLayerMixing, ::SpeedyWeather.AbstractModel) = nothing

"""Lazy column of virtual potential temperature referenced to the surface pressure, Tv sigma^-kappa."""
struct _PBLThetaV{T, Q, S, A, NF} <: AbstractVector{NF}
    temperature::T
    humidity::Q
    sigma::S
    atmosphere::A
    kappa::NF
    ij::Int
end
Base.size(x::_PBLThetaV) = (length(x.sigma),)
Base.@propagate_inbounds Base.getindex(x::_PBLThetaV, k::Int) =
    SpeedyWeather.virtual_temperature(x.temperature[x.ij, k], x.humidity[x.ij, k], x.atmosphere) *
    x.sigma[k]^(-x.kappa)

"""Lazy column of resolved squared wind speed."""
struct _PBLWind2{U, V, NF} <: AbstractVector{NF}
    u::U
    v::V
    ij::Int
    n::Int
end
Base.size(x::_PBLWind2) = (x.n,)
Base.@propagate_inbounds Base.getindex(x::_PBLWind2{U, V, NF}, k::Int) where {U, V, NF} =
    NF(x.u[x.ij, k]^2 + x.v[x.ij, k]^2)

"""Lazy column of dry static energy cp T + Phi."""
struct _PBLDryStaticEnergy{T, P, NF} <: AbstractVector{NF}
    temperature::T
    geopotential::P
    heat_capacity::NF
    ij::Int
    n::Int
end
Base.size(x::_PBLDryStaticEnergy) = (x.n,)
Base.@propagate_inbounds Base.getindex(x::_PBLDryStaticEnergy, k::Int) =
    x.heat_capacity * x.temperature[x.ij, k] + x.geopotential[x.ij, k]

"""Lazy column of face coefficients a[k] (1/s), k = 1:n-1, from the K-profile at the hydrostatic face heights."""
struct _PBLFaceCoefficients{Z, T, Q, S, H, A, NF} <: AbstractVector{NF}
    heights::Z
    temperature::T
    humidity::Q
    sigma_full::S
    sigma_half::H
    atmosphere::A
    ij::Int
    h::NF
    ustar::NF
    buoyancy_flux::NF
    thetav_reference::NF
    gravity::NF
    gas_constant::NF
end
Base.size(x::_PBLFaceCoefficients) = (length(x.sigma_full) - 1,)
Base.@propagate_inbounds function Base.getindex(x::_PBLFaceCoefficients{Z, T, Q, S, H, A, NF}, k::Int) where {Z, T, Q, S, H, A, NF}
    ij = x.ij
    tv_below = SpeedyWeather.virtual_temperature(x.temperature[ij, k + 1], x.humidity[ij, k + 1], x.atmosphere)
    tv_above = SpeedyWeather.virtual_temperature(x.temperature[ij, k], x.humidity[ij, k], x.atmosphere)
    z_face = x.heights[ij, k + 1] + x.gas_constant * tv_below / x.gravity * log(x.sigma_full[k + 1] / x.sigma_half[k + 1])
    K = pbl_diffusivity(z_face, x.h, x.ustar, x.buoyancy_flux, x.thetav_reference, x.gravity)
    tv_face = (tv_above + tv_below) / 2
    metric = (x.sigma_half[k + 1] * x.gravity / (x.gas_constant * tv_face))^2
    return NF(metric * K / (x.sigma_full[k + 1] - x.sigma_full[k]))
end

"""Momentum drag coefficient exactly as the surface momentum flux scheme selects it, read from the lagged state."""
Base.@propagate_inbounds function _pbl_drag_coefficient(ij, vars, momentum_flux, land_fraction)
    if hasproperty(momentum_flux, :use_boundary_layer_drag) && !momentum_flux.use_boundary_layer_drag
        return land_fraction * momentum_flux.drag_land + (1 - land_fraction) * momentum_flux.drag_ocean
    end
    return vars.parameterizations.boundary_layer_drag[ij]
end

"""Conditional ocean/land flux, weighted once by area, without reading a zero-area tile."""
@inline _pbl_area_weighted(land_value, ocean_value, land_fraction) =
    (land_fraction > 0 ? land_fraction * land_value : zero(land_value)) +
    (land_fraction < 1 ? (1 - land_fraction) * ocean_value : zero(ocean_value))

"""Add the flux-form tendency of the solved column `x_new` (times `scale`) to `tendency[ij, :]`."""
Base.@propagate_inbounds function _pbl_add_flux_tendency!(tendency, ij, x_new, a, dsigma, scale)
    n = length(dsigma)
    flux_above = zero(eltype(x_new))
    for k in 1:n
        flux_below = k < n ? a[k] * (x_new[k + 1] - x_new[k]) : zero(flux_above)
        tendency[ij, k] += scale * (flux_below - flux_above) / dsigma[k]
        flux_above = flux_below
    end
    return nothing
end

Base.@propagate_inbounds function SpeedyWeather.parameterization!(ij, vars, scheme::KProfileBoundaryLayerMixing, model)
    (; geometry, planet, atmosphere, time_stepping) = model
    temperature = vars.grid.temperature_prev
    humidity = vars.grid.humidity_prev
    u = vars.grid.u_prev
    v = vars.grid.v_prev
    NF = eltype(temperature)
    sigma = geometry.σ_levels_full
    sigma_half = geometry.σ_levels_half
    dsigma = geometry.σ_levels_thick
    n = length(sigma)
    gravity = NF(planet.gravity)
    gas_constant = NF(atmosphere.R_dry)
    kappa = NF(atmosphere.κ)
    heat_capacity = NF(atmosphere.heat_capacity)
    surface_pressure = vars.grid.pressure_prev[ij]

    # Hydrostatic full-level heights above the surface.
    heights = vars.parameterizations.pbl_height_work
    tv_lowest = SpeedyWeather.virtual_temperature(temperature[ij, n], humidity[ij, n], atmosphere)
    heights[ij, n] = gas_constant * tv_lowest / gravity * log(1 / sigma[n])
    for k in (n - 1):-1:1
        tv_mean = (SpeedyWeather.virtual_temperature(temperature[ij, k], humidity[ij, k], atmosphere) +
                   SpeedyWeather.virtual_temperature(temperature[ij, k + 1], humidity[ij, k + 1], atmosphere)) / 2
        heights[ij, k] = heights[ij, k + 1] + gas_constant * tv_mean / gravity * log(sigma[k + 1] / sigma[k])
    end
    z_lowest_face = heights[ij, n] + gas_constant * tv_lowest / gravity * log(sigma[n] / sigma_half[n])

    # Surface quantities in the SurfaceCondition form, without calling the mutating scheme early.
    condition = model.surface_condition
    slowdown = NF(condition.wind_slowdown)
    gust = NF(condition.gust_speed)
    u1, v1 = u[ij, n], v[ij, n]
    surface_wind = sqrt((slowdown * u1)^2 + (slowdown * v1)^2 + gust^2)
    density = surface_pressure / (gas_constant * tv_lowest * sigma[n]^(-kappa))
    land_fraction = NF(model.land_sea_mask.mask[ij])
    drag = NF(_pbl_drag_coefficient(ij, vars, model.surface_momentum_flux, land_fraction))
    momentum_slowdown = NF(model.surface_momentum_flux.wind_slowdown)
    ustar = sqrt(max(drag * surface_wind * momentum_slowdown * sqrt(u1^2 + v1^2), zero(NF)))
    sensible = NF(_pbl_area_weighted(vars.prognostic.land.sensible_heat_flux[ij],
                                     vars.prognostic.ocean.sensible_heat_flux[ij], land_fraction))
    evaporation = NF(_pbl_area_weighted(vars.prognostic.land.surface_humidity_flux[ij],
                                        vars.prognostic.ocean.surface_humidity_flux[ij], land_fraction))
    buoyancy_flux = sensible / (density * heat_capacity) +
                    NF(PBL_VIRTUAL_FACTOR) * temperature[ij, n] * evaporation / density

    thetav = _PBLThetaV(temperature, humidity, sigma, atmosphere, kappa, ij)
    wind2 = _PBLWind2{typeof(u), typeof(v), NF}(u, v, ij, n)
    z_full = view(heights.data, ij, :)
    thetav_reference = thetav[n]
    h = pbl_height(z_full, thetav, wind2, ustar, buoyancy_flux, thetav_reference, gravity, z_lowest_face,
                   scheme.maximum_height; critical = scheme.critical_richardson)
    a = _PBLFaceCoefficients(heights, temperature, humidity, sigma, sigma_half, atmosphere, ij, NF(h), ustar,
                             buoyancy_flux, thetav_reference, gravity, gas_constant)

    vars.parameterizations.pbl_height[ij] = h
    vars.parameterizations.pbl_friction_velocity[ij] = ustar
    vars.parameterizations.pbl_buoyancy_flux[ij] = buoyancy_flux
    vars.parameterizations.pbl_drag_coefficient[ij] = drag
    vars.parameterizations.pbl_lowest_face_diffusivity[ij] =
        pbl_diffusivity(z_lowest_face, h, ustar, buoyancy_flux, thetav_reference, gravity)

    dt = 2 * NF(time_stepping.Δt_sec)
    c = view(vars.scratch.grid.a.data, ij, :)
    d = view(vars.scratch.grid.b.data, ij, :)
    x_new = view(vars.parameterizations.pbl_solution_work.data, ij, :)
    dry_static_energy = _PBLDryStaticEnergy(temperature, vars.grid.geopotential, heat_capacity, ij, n)
    pbl_implicit_diffuse!(x_new, dry_static_energy, a, dsigma, dt, c, d)
    _pbl_add_flux_tendency!(vars.tendencies.grid.temperature, ij, x_new, a, dsigma, inv(heat_capacity))
    pbl_implicit_diffuse!(x_new, view(humidity.data, ij, :), a, dsigma, dt, c, d)
    _pbl_add_flux_tendency!(vars.tendencies.grid.humidity, ij, x_new, a, dsigma, one(NF))
    return nothing
end
