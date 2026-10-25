#=
Opt-in boundary-layer K-profile mixing in SpeedyWeather's vertical_diffusion slot (Stage 2 host wiring).

Selected by `atmosphere_vertical_diffusion = :holtslag_boville_k_profile_v1`. Production keeps the zero operator. The
column physics is in pbl_k_profile.jl. Host conventions, agreed with Codex's interface review:
- vertical_diffusion runs first in each column, before convection and the surface fluxes, and SpeedyWeather zeroes the
  top-level surface flux totals before column physics. The surface buoyancy flux therefore comes from the persistent
  prescribed coupled inputs, `prognostic.ocean` and `prognostic.land` sensible heat (W/m2) and humidity (kg/m2/s)
  fluxes, positive upward, conditional means weighted once by the land fraction (a zero-area tile is never read).
  B_s is the surface flux of the virtual potential temperature that the Richardson number uses, Tv sigma^-kappa:
  B_s = (1 + mu q_1) H / (rho cp) + mu T_s E / rho, with T_s = T_1 sigma_n^-kappa and mu = atmosphere.μ_virt_temp;
- the friction velocity is the stress equivalent of SpeedyWeather's surface momentum flux,
  u*^2 = C_D V0 f |V1|, with the lagged stored drag (or the fixed drag the momentum scheme would use), V0 and the
  surface air density computed locally in the SurfaceCondition form. The shared drag and stress are unchanged. At
  start-up C_D can be zero; then u* = 0 and the calm-safe closure applies (no neutral fallback);
- dry static energy s = cp T + Phi and humidity are mixed from the `_prev` state over 2 dt, and the flux-form tendency
  of the backward-Euler solution is added (the leapfrog update over 2 dt then equals the implicit step). Momentum is
  not mixed in v1. Heights above the surface are hydrostatic from the column's virtual temperature;
- the depth h is at least the top face of layer 1 and at most the height of sigma = 0.5 in that column. The bulk
  Richardson number has a 0.1 m/s calm-wind floor (PBL_CALM_WIND_FLOOR), part of the provenance below.
- strength-screen parameters (diagnostic experiments only, defaults keep the scheme bitwise): `diffusivity_scale`
  in (0, 1] multiplies every face coefficient after it is computed (and the reported lowest-face diffusivity), so a
  scale of 1 is an exact multiply. `mix_heat` and `mix_humidity` skip the solve and tendency of that scalar. The
  boundary-layer depth, Ric and the depth cap do not change. A non-default value adds a provenance suffix.
=#

const PBL_K_PROFILE_PROVENANCE =
    "holtslag_boville_1993_k_profile_ric0p3_calm_wind_floor_0p1ms_sigma0p5_cap_backward_euler_sigma_flux_scalars_v1"
const PBL_MAXIMUM_SIGMA = 0.5

"""
Holtslag-Boville K-profile mixing of dry static energy and humidity (opt-in vertical diffusion).

On the GPU, SpeedyWeather passes each column scheme a device copy of the model that holds only the global fields
(geometry, planet, atmosphere, time stepping, land-sea mask); the other parameterizations are separate kernel
arguments. The scheme therefore carries the scalar parameters it needs from the surface condition and the surface
momentum flux. `initialize!` checks them against the model's own schemes, so they cannot drift apart.
"""
struct KProfileBoundaryLayerMixing{NF} <: SpeedyWeather.AbstractVerticalDiffusion
    critical_richardson::NF
    maximum_sigma::NF
    surface_wind_slowdown::NF    # SurfaceCondition.wind_slowdown
    gust_speed::NF               # SurfaceCondition.gust_speed (m/s)
    momentum_wind_slowdown::NF   # SurfaceMomentumFlux.wind_slowdown
    use_boundary_layer_drag::Bool
    drag_ocean::NF
    drag_land::NF
    diffusivity_scale::NF        # strength screen: multiplies the face coefficients (default 1, exact)
    mix_heat::Bool               # strength screen: mix dry static energy (default true)
    mix_humidity::Bool           # strength screen: mix humidity (default true)
    mix_momentum::Bool           # opt-in: mix u and v with the same coefficients (default false)
end

Adapt.@adapt_structure KProfileBoundaryLayerMixing

function KProfileBoundaryLayerMixing(spectral_grid::SpeedyWeather.SpectralGrid;
                                     critical_richardson = PBL_CRITICAL_RICHARDSON, maximum_sigma = PBL_MAXIMUM_SIGMA,
                                     surface_condition = SpeedyWeather.SurfaceCondition(spectral_grid),
                                     surface_momentum_flux = SpeedyWeather.SurfaceMomentumFlux(spectral_grid),
                                     diffusivity_scale = 1, mix_heat = true, mix_humidity = true,
                                     mix_momentum = false)
    0 < diffusivity_scale <= 1 || throw(ArgumentError("diffusivity_scale must be in (0, 1], got $diffusivity_scale"))
    (mix_heat || mix_humidity || mix_momentum) ||
        throw(ArgumentError("the K-profile must mix at least one of heat, humidity and momentum"))
    NF = spectral_grid.NF
    return KProfileBoundaryLayerMixing{NF}(NF(critical_richardson), NF(maximum_sigma),
        NF(surface_condition.wind_slowdown), NF(surface_condition.gust_speed),
        NF(surface_momentum_flux.wind_slowdown), Bool(surface_momentum_flux.use_boundary_layer_drag),
        NF(surface_momentum_flux.drag_ocean), NF(surface_momentum_flux.drag_land),
        NF(diffusivity_scale), Bool(mix_heat), Bool(mix_humidity), Bool(mix_momentum))
end

"""Provenance: the v1 string, plus a suffix for each non-default strength-screen parameter."""
function _atmosphere_vertical_diffusion_provenance(scheme::KProfileBoundaryLayerMixing)
    suffix = ""
    scheme.diffusivity_scale == 1 ||
        (suffix *= "_scale" * replace(string(round(Float64(scheme.diffusivity_scale); digits = 4)), "." => "p"))
    if scheme.mix_momentum && !scheme.mix_heat && !scheme.mix_humidity
        suffix *= "_momentum_only"
    else
        scheme.mix_humidity || (suffix *= "_heat_only")
        scheme.mix_heat || (suffix *= "_humidity_only")
        scheme.mix_momentum && (suffix *= "_with_momentum")
    end
    return PBL_K_PROFILE_PROVENANCE * suffix
end

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
    SpeedyWeather.ParameterizationVariable(:pbl_lowest_face_height, SpeedyWeather.Grid2D();
        desc = "Height of the top face of the lowest layer, the lower bound of the K-profile depth", units = "m"),
    SpeedyWeather.ParameterizationVariable(:pbl_momentum_dissipation_face, SpeedyWeather.Grid2D();
        desc = "Kinetic-energy loss rate of the K-profile momentum update, face term (zero if momentum is off)",
        units = "W/m^2"),
    SpeedyWeather.ParameterizationVariable(:pbl_momentum_dissipation_implicit, SpeedyWeather.Grid2D();
        desc = "Kinetic-energy loss rate of the K-profile momentum update, backward-Euler term (zero if off)",
        units = "W/m^2"),
    SpeedyWeather.ParameterizationVariable(:pbl_height_work, SpeedyWeather.Grid3D();
        desc = "K-profile full-level heights above the surface (work column)", units = "m"),
    SpeedyWeather.ParameterizationVariable(:pbl_solution_work, SpeedyWeather.Grid3D();
        desc = "K-profile implicit solution (work column)", units = "?"),
)

function SpeedyWeather.initialize!(scheme::KProfileBoundaryLayerMixing, model::SpeedyWeather.AbstractModel)
    condition, momentum = model.surface_condition, model.surface_momentum_flux
    condition isa SpeedyWeather.SurfaceCondition && momentum isa SpeedyWeather.SurfaceMomentumFlux ||
        throw(ArgumentError("the K-profile scheme supports SpeedyWeather's SurfaceCondition and SurfaceMomentumFlux"))
    matches = scheme.surface_wind_slowdown == condition.wind_slowdown && scheme.gust_speed == condition.gust_speed &&
        scheme.momentum_wind_slowdown == momentum.wind_slowdown &&
        scheme.use_boundary_layer_drag == momentum.use_boundary_layer_drag &&
        scheme.drag_ocean == momentum.drag_ocean && scheme.drag_land == momentum.drag_land
    matches || throw(ArgumentError("K-profile surface parameters differ from the model's surface condition or " *
                                   "momentum flux; build the scheme from the same objects"))
    return nothing
end

const PBL_K_PROFILE_RESTART_SCHEMA = "readiesm_pbl_k_profile_restart_v1"

"""
Restart record of the K-profile settings, versioned, with exact values (not the rounded provenance suffix). It is
present only when the K-profile is the vertical diffusion, so checkpoints of the production zero operator and of the
bulk Richardson option do not change.
"""
function _pbl_k_profile_restart_state(simulation)
    scheme = simulation.model.vertical_diffusion
    scheme isa KProfileBoundaryLayerMixing || return (;)
    return (; pbl_k_profile_restart_schema = PBL_K_PROFILE_RESTART_SCHEMA,
            pbl_diffusivity_scale = Float64(scheme.diffusivity_scale), pbl_mix_heat = scheme.mix_heat,
            pbl_mix_humidity = scheme.mix_humidity, pbl_mix_momentum = scheme.mix_momentum)
end

"""
Refuse a checkpoint whose K-profile record differs from the reconstructed scheme: K-profile on to off, off to on, an
unknown record version, or any change of the scale or of the heat, humidity or momentum flags. A K-profile
configuration also refuses a checkpoint without the record. Callers run this before any component state is copied.
"""
function _validate_pbl_k_profile_restart(simulation, state)
    scheme = simulation.model.vertical_diffusion
    active = scheme isa KProfileBoundaryLayerMixing
    saved = hasproperty(state, :pbl_k_profile_restart_schema)
    if !active
        saved && throw(ArgumentError(
            "checkpoint was written with the K-profile vertical diffusion, but the configuration does not use it"))
        return nothing
    end
    saved || throw(ArgumentError(
        "the configuration uses the K-profile vertical diffusion, but the checkpoint has no K-profile restart record"))
    state.pbl_k_profile_restart_schema == PBL_K_PROFILE_RESTART_SCHEMA || throw(ArgumentError(
        "unknown K-profile restart record $(state.pbl_k_profile_restart_schema)"))
    expected = _pbl_k_profile_restart_state(simulation)
    for key in (:pbl_diffusivity_scale, :pbl_mix_heat, :pbl_mix_humidity, :pbl_mix_momentum)
        hasproperty(state, key) || throw(ArgumentError("K-profile restart record lacks $key"))
        value = getproperty(state, key)
        typeof(value) == typeof(getproperty(expected, key)) && value == getproperty(expected, key) ||
            throw(ArgumentError("checkpoint and configuration differ in $key " *
                                "(checkpoint $value, configuration $(getproperty(expected, key)))"))
    end
    return nothing
end

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
    scale::NF
end
_PBLFaceCoefficients(heights, temperature, humidity, sigma_full, sigma_half, atmosphere, ij, h::NF, ustar,
                     buoyancy_flux, thetav_reference, gravity, gas_constant) where {NF} =
    _PBLFaceCoefficients(heights, temperature, humidity, sigma_full, sigma_half, atmosphere, ij, h, ustar,
                         buoyancy_flux, thetav_reference, gravity, gas_constant, one(NF))
Base.size(x::_PBLFaceCoefficients) = (length(x.sigma_full) - 1,)
Base.@propagate_inbounds function Base.getindex(x::_PBLFaceCoefficients{Z, T, Q, S, H, A, NF}, k::Int) where {Z, T, Q, S, H, A, NF}
    ij = x.ij
    tv_below = SpeedyWeather.virtual_temperature(x.temperature[ij, k + 1], x.humidity[ij, k + 1], x.atmosphere)
    tv_above = SpeedyWeather.virtual_temperature(x.temperature[ij, k], x.humidity[ij, k], x.atmosphere)
    z_face = x.heights[ij, k + 1] + x.gas_constant * tv_below / x.gravity * log(x.sigma_full[k + 1] / x.sigma_half[k + 1])
    K = pbl_diffusivity(z_face, x.h, x.ustar, x.buoyancy_flux, x.thetav_reference, x.gravity)
    tv_face = (tv_above + tv_below) / 2
    metric = (x.sigma_half[k + 1] * x.gravity / (x.gas_constant * tv_face))^2
    # The scale multiplies the finished coefficient, so a scale of 1 leaves it bitwise unchanged.
    return NF(metric * K / (x.sigma_full[k + 1] - x.sigma_full[k])) * x.scale
end

"""Momentum drag coefficient exactly as the surface momentum flux scheme selects it, read from the lagged state."""
Base.@propagate_inbounds function _pbl_drag_coefficient(ij, vars, scheme::KProfileBoundaryLayerMixing, land_fraction)
    scheme.use_boundary_layer_drag && return vars.parameterizations.boundary_layer_drag[ij]
    return land_fraction * scheme.drag_land + (1 - land_fraction) * scheme.drag_ocean
end

"""Conditional ocean/land flux at `ij`, weighted once by area. A tile with zero area is not read."""
Base.@propagate_inbounds function _pbl_area_weighted(land_field, ocean_field, ij, land_fraction)
    value = zero(land_fraction)
    land_fraction > 0 && (value += land_fraction * land_field[ij])
    land_fraction < 1 && (value += (1 - land_fraction) * ocean_field[ij])
    return value
end

"""Height (m) above the surface of the sigma level `target`, with the same hydrostatic layer-mean virtual temperatures
as the full-level heights in `heights` (log-sigma form between full levels)."""
Base.@propagate_inbounds function _pbl_sigma_height(heights, temperature, humidity, sigma, atmosphere, ij, target,
                                                    gas_constant, gravity)
    n = length(sigma)
    tv(k) = SpeedyWeather.virtual_temperature(temperature[ij, k], humidity[ij, k], atmosphere)
    target >= sigma[n] && return gas_constant * tv(n) / gravity * log(1 / target)
    for k in (n - 1):-1:1
        if target >= sigma[k]
            return heights[ij, k + 1] + gas_constant * (tv(k) + tv(k + 1)) / 2 / gravity * log(sigma[k + 1] / target)
        end
    end
    return heights[ij, 1] + gas_constant * tv(1) / gravity * log(sigma[1] / target)
end

"""
Mix one wind component: backward-Euler solve of `wind` into `x_new`, add the flux-form tendency (m/s2) to
`tendency[ij, :]`, and return the face sum sum_k a[k] (x_new[k+1] - x_new[k])^2 and the implicit sum
sum_k dsigma[k] (x_new[k] - wind[k])^2 used by the kinetic-energy diagnostics.
"""
Base.@propagate_inbounds function _pbl_mix_wind!(tendency, wind, ij, x_new, a, dsigma, dt, c, d)
    pbl_implicit_diffuse!(x_new, wind, a, dsigma, dt, c, d)
    _pbl_add_flux_tendency!(tendency, ij, x_new, a, dsigma, one(eltype(x_new)))
    n = length(dsigma)
    face_sum = zero(eltype(x_new))
    implicit_sum = zero(eltype(x_new))
    for k in 1:n
        k < n && (face_sum += a[k] * (x_new[k + 1] - x_new[k])^2)
        implicit_sum += dsigma[k] * (x_new[k] - wind[k])^2
    end
    return face_sum, implicit_sum
end

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

    # Surface quantities in the SurfaceCondition form, without calling the mutating scheme early. The parameters come
    # from the scheme (see the struct docstring), not from model.surface_condition, which the GPU kernel lacks.
    slowdown = scheme.surface_wind_slowdown
    gust = scheme.gust_speed
    u1, v1 = u[ij, n], v[ij, n]
    surface_wind = sqrt((slowdown * u1)^2 + (slowdown * v1)^2 + gust^2)
    density = surface_pressure / (gas_constant * tv_lowest * sigma[n]^(-kappa))
    land_fraction = NF(model.land_sea_mask.mask[ij])
    drag = NF(_pbl_drag_coefficient(ij, vars, scheme, land_fraction))
    momentum_slowdown = scheme.momentum_wind_slowdown
    ustar = sqrt(max(drag * surface_wind * momentum_slowdown * sqrt(u1^2 + v1^2), zero(NF)))
    sensible = NF(_pbl_area_weighted(vars.prognostic.land.sensible_heat_flux,
                                     vars.prognostic.ocean.sensible_heat_flux, ij, land_fraction))
    evaporation = NF(_pbl_area_weighted(vars.prognostic.land.surface_humidity_flux,
                                        vars.prognostic.ocean.surface_humidity_flux, ij, land_fraction))
    # Surface flux of Tv sigma^-kappa (the Richardson number's variable) at sigma = 1, from the extrapolated T_s.
    mu = NF(atmosphere.μ_virt_temp)
    surface_temperature = temperature[ij, n] * sigma[n]^(-kappa)
    buoyancy_flux = (1 + mu * humidity[ij, n]) * sensible / (density * heat_capacity) +
                    mu * surface_temperature * evaporation / density

    thetav = _PBLThetaV(temperature, humidity, sigma, atmosphere, kappa, ij)
    wind2 = _PBLWind2{typeof(u), typeof(v), NF}(u, v, ij, n)
    z_full = view(heights, ij, :)
    thetav_reference = thetav[n]
    z_maximum = _pbl_sigma_height(heights, temperature, humidity, sigma, atmosphere, ij, scheme.maximum_sigma,
                                  gas_constant, gravity)
    h = pbl_height(z_full, thetav, wind2, ustar, buoyancy_flux, thetav_reference, gravity, z_lowest_face, z_maximum;
                   critical = scheme.critical_richardson)
    # One rounded depth and face height for the coefficients and every stored diagnostic, so that the reported
    # lowest-face K is the applied K (h == z_face gives K == 0 in both).
    h_stored = NF(h)
    z_face_stored = NF(z_lowest_face)
    a = _PBLFaceCoefficients(heights, temperature, humidity, sigma, sigma_half, atmosphere, ij, h_stored, ustar,
                             buoyancy_flux, thetav_reference, gravity, gas_constant, scheme.diffusivity_scale)

    vars.parameterizations.pbl_height[ij] = h_stored
    vars.parameterizations.pbl_lowest_face_height[ij] = z_face_stored
    vars.parameterizations.pbl_friction_velocity[ij] = ustar
    vars.parameterizations.pbl_buoyancy_flux[ij] = buoyancy_flux
    vars.parameterizations.pbl_drag_coefficient[ij] = drag
    vars.parameterizations.pbl_lowest_face_diffusivity[ij] =
        NF(pbl_diffusivity(z_face_stored, h_stored, ustar, buoyancy_flux, thetav_reference, gravity)) *
        scheme.diffusivity_scale

    dt = 2 * NF(time_stepping.Δt_sec)
    # Column views of the arrays themselves: Fields on the CPU, plain device matrices in the GPU kernel.
    c = view(vars.scratch.grid.a, ij, :)
    d = view(vars.scratch.grid.b, ij, :)
    x_new = view(vars.parameterizations.pbl_solution_work, ij, :)
    if scheme.mix_heat
        dry_static_energy = _PBLDryStaticEnergy(temperature, vars.grid.geopotential, heat_capacity, ij, n)
        pbl_implicit_diffuse!(x_new, dry_static_energy, a, dsigma, dt, c, d)
        _pbl_add_flux_tendency!(vars.tendencies.grid.temperature, ij, x_new, a, dsigma, inv(heat_capacity))
    end
    if scheme.mix_humidity
        pbl_implicit_diffuse!(x_new, view(humidity, ij, :), a, dsigma, dt, c, d)
        _pbl_add_flux_tendency!(vars.tendencies.grid.humidity, ij, x_new, a, dsigma, one(NF))
    end
    if scheme.mix_momentum
        # Momentum uses the same coefficients a (HB93 K_m = K, which the scalars already use: Pr = 1), with no flux
        # at the top and at the surface. SpeedyWeather's surface drag stays the only column momentum sink.
        face_sum, implicit_sum = _pbl_mix_wind!(vars.tendencies.grid.u, view(u, ij, :), ij, x_new, a, dsigma, dt, c, d)
        face_v, implicit_v = _pbl_mix_wind!(vars.tendencies.grid.v, view(v, ij, :), ij, x_new, a, dsigma, dt, c, d)
        # Kinetic-energy loss of this frozen-coefficient solve of length dt, per unit area (W/m2):
        # E_old - E_new = dt * D_face + (p_s / 2g) sum dsigma (x_new - x_old)^2 = dt * (D_face + D_implicit).
        mass = surface_pressure / gravity
        vars.parameterizations.pbl_momentum_dissipation_face[ij] = mass * (face_sum + face_v)
        vars.parameterizations.pbl_momentum_dissipation_implicit[ij] = mass * (implicit_sum + implicit_v) / (2 * dt)
    end
    return nothing
end
