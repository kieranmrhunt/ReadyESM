# Tiedtke-Bechtold convection, stage 1a: the pure column function (no host wiring).
#
# Frozen contract: artifacts/development_wave_20261103/tiedtke_bechtold_claude_v1/frozen_v4/FROZEN_CONTRACT_v4.sha256
# (sha256 b07062d259455a207fb45afd1d2a0283387c167f10f8e7a57bb150d36b8c0e76), which pins DESIGN.md 42b7a0d6,
# STAGE1_CHOICES_DRAFT.md 211627e7 (r15) and SPEC_CHAPTER6.md 57d18c70. "Dnn" and "Fn" below are the decision rows
# of STAGE1_CHOICES_DRAFT.md. Stage 1a scope (F1-F7): deep-only trigger from k = n-2, updraught only, liquid only,
# rain only under the water contract instantaneous_rain_no_rup_v1 (no rain loading), PCAPE_bl = 0, strict CFL, no
# downdraughts, no shallow or mid-level convection, no momentum and no surface-flux input.
#
# Levels run from k = 1 (top) to n (surface). Half-level arrays have n+1 entries and index j holds half level j-1/2,
# so layer k lies between p_half[k] (top) and p_half[k+1] (bottom), and the flux "at k+1/2" is stored at j = k+1.
# Everything here is allocation-free, isbits and branch-only (no dynamic dispatch), so the same function can run
# inside a CPU or GPU column kernel and in the offline replay and SCM harness.
#
# Stage 1B (frozen_v5, shallow branch and surface-flux parcel, column science only) is the second method of
# `tiedtke_bechtold_column!`, after the stage 1a method; see the section "Stage 1B, column science" below.

const TIEDTKE_BECHTOLD_CONTRACT_MANIFEST_SHA256 =
    "b07062d259455a207fb45afd1d2a0283387c167f10f8e7a57bb150d36b8c0e76"
const TIEDTKE_BECHTOLD_WATER_CONTRACT = "instantaneous_rain_no_rup_v1"   # D28
const TIEDTKE_BECHTOLD_KE_VERSION = "tb_ke_v1"                           # D20, D21
const TIEDTKE_BECHTOLD_ORGDET_VERSION = "tb_orgdet_v1"                   # D19
const TIEDTKE_BECHTOLD_TRIGGER_VARIANT = "deep_from_n_minus_2"           # D07, T4 contract value
const TIEDTKE_BECHTOLD_PARCEL_ADAPTER = "none"                           # stage 1a, T4 contract value
const TIEDTKE_BECHTOLD_CLOSURE_NAME =
    "pre-Cy40r1 equilibrium CAPE closure (IFS Cy41r2 Eq. 6.25 with PCAPE_bl = 0)"  # D31, F5

"""
Parameters of the stage 1a Tiedtke-Bechtold column (STAGE1_CHOICES_DRAFT.md r15, Sect. 9). Source tags per field:
IFS = IFS Cy41r2 Part IV Chapter 6 as printed, IFS-R = Cy41r2 with a stated reading, RE = ReadyESM choice,
HOST = SpeedyWeather 0.21.1 value, always read from the live planet and atmosphere objects (see
`tiedtke_bechtold_parameters`), never a literal copy. The struct is isbits, so it passes to GPU kernels unchanged.

Rows of Sect. 9 that stage 1a does not use are deliberately absent: c_BF, T_BF and T_ice (D26 not applied), the
(6.19)-(6.20) stage 1b values and floor, and rho_w (used only by the host wrapper, D43).
"""
struct TiedtkeBechtoldParameters{NF <: AbstractFloat}
    # HOST constants (SpeedyWeather 0.21.1 EarthAtmosphere and Earth; rules 5, D02)
    gravity::NF                    # g, m s-2
    heat_capacity::NF              # c_p, J kg-1 K-1
    latent_heat_condensation::NF   # L_v, J kg-1
    R_dry::NF                      # R_d, J kg-1 K-1
    R_vapor::NF                    # R_v, J kg-1 K-1
    mol_ratio::NF                  # eps = R_d/R_v
    mu_virt_temp::NF               # host mu = 0.6078, host hydrostatics only (D04)
    saturation_vapor_pressure::NF  # e_0, Pa (host q_sat)
    temperature_freezing::NF       # T_0, K (host q_sat)
    # Updraught entrainment and detrainment (IFS (6.7), (6.8), (6.12), p. 77-78)
    eps1::NF                       # eps^(1) = 1.75e-3 m-1, IFS
    f_eps::NF                      # f_eps (deep) = 1, IFS
    rh_offset_entrainment::NF      # 1.3 in (6.7), IFS
    rh_offset_detrainment::NF      # 1.6 in (6.8) and (6.12), IFS
    f_scale_exponent::Int          # 3 in (6.7), IFS
    delta1::NF                     # delta^(1) = 0.75e-4 m-1, IFS
    # Updraught kinetic energy (IFS p. 78)
    gamma::NF                      # virtual mass coefficient 0.5, IFS
    f_turbulence::NF               # f = 2, IFS
    C_d::NF                        # drag coefficient 0.506, IFS
    beta::NF                       # 1.875, IFS
    # Deep test parcel (IFS U3, (6.21), p. 80-81)
    test_entrainment_factor::NF    # 0.4, IFS
    test_entrainment_exponent::Int # 3, IFS
    test_condensate_removal::NF    # 0.5 (50% removed at each level), IFS
    delta_T_deep::NF               # 0.2 K, IFS (6.21)
    delta_q_deep::NF               # 1e-4 kg/kg, IFS (6.21)
    w_departure::NF                # 1 m s-1 at the departure level, IFS
    mixed_layer_depth::NF          # 6000 Pa, IFS
    average_depth::NF              # 5000 Pa (50 hPa average), IFS
    window_half_width::NF          # 2500 Pa, RE (D09)
    departure_range::NF            # 35000 Pa, IFS
    deep_depth::NF                 # 20000 Pa, IFS
    # Deep closure (IFS (6.24)-(6.26), p. 81-82)
    mstar_base_factor::NF          # 0.1 in U6, IFS
    alpha_x_coefficient::NF        # 1.66, IFS
    dx_ref::NF                     # 125e3 m, IFS
    tau_min::NF                    # 720 s, IFS
    tau_max::NF                    # 10800 s, IFS
    PCAPE_bl::NF                   # 0 J m-3, IFS (F5, pre-Cy40r1)
    subcloud_exponent::Int         # m = 1 (deep), IFS (6.33)-(6.34)
    grid_spacing::NF               # dx = pi a / N_trunc, RE (D35)
    # Updraught microphysics and rain (IFS (6.36), (6.44), p. 84-85)
    c0::NF                         # 1.4e-3 s-1, IFS
    l_crit::NF                     # 5e-4 kg/kg, IFS
    w_factor::NF                   # 0.75, IFS
    w_cap::NF                      # 10 m s-1, IFS
    l_threshold_water::NF          # 3e-4 kg/kg, IFS
    l_threshold_land::NF           # 5e-4 kg/kg, IFS
    alpha1::NF                     # 5.44e-4 s-1, IFS
    alpha2::NF                     # 5.09e-3, IFS
    alpha3::NF                     # 0.5777, IFS
    C_conv::NF                     # 0.05, IFS
    RH_cr_water::NF                # 0.9, IFS-R
    RH_cr_land::NF                 # 0.7, IFS-R
    virtual_coefficient::NF        # 0.608 in (6.56), IFS
    cfl_factor::NF                 # 1, IFS (F4 strict)
    limiter_timestep_factor::Int   # 2, RE (D41 r15): the host leapfrog applies tendencies over 2 dt; 1 for forward Euler
    # ReadyESM choices
    newton_iterations::Int         # 3, RE (D03)
    rh_min::NF                     # 0, RE (D05 clip)
    rh_max::NF                     # 1, RE (D05 clip)
    land_threshold::NF             # 0.5 land fraction, RE (D30)
    rain_loading::Bool             # false, D28 instantaneous_rain_no_rup_v1
end

"""
    TiedtkeBechtoldParameters(NF; gravity, heat_capacity, latent_heat_condensation, R_dry, R_vapor, mol_ratio,
                              mu_virt_temp, saturation_vapor_pressure, temperature_freezing, radius,
                              truncation = 31, kwargs...)

Stage 1a parameters with the Sect. 9 values. The host constants (g, c_p, L_v, R_d, R_v, eps, the host mu and the
q_sat constants e_0 and T_0) and the planet radius are required keywords with no defaults: they must come from the
live host objects, as the production kernels read them (`model.planet.gravity`, `model.atmosphere.heat_capacity`,
...), so the column ledger closes with the host's own constants (rule 5). Derived host quantities (for example
kappa = R_d/c_p, used only by the test profiles) are always computed from these values. Use
`tiedtke_bechtold_parameters` to build them from a model or from `SpeedyWeather.Earth`/`EarthAtmosphere`.
The other keyword arguments override single Sect. 9 values for sensitivity work. The fixed stage 1a decisions
(F4 cfl_factor = 1, F5 PCAPE_bl = 0, deep m = 1, D28 rain_loading = false) are refused if changed.
"""
function TiedtkeBechtoldParameters(
    ::Type{NF};
    gravity::Real, heat_capacity::Real, latent_heat_condensation::Real, R_dry::Real, R_vapor::Real,
    mol_ratio::Real, mu_virt_temp::Real, saturation_vapor_pressure::Real, temperature_freezing::Real,
    radius::Real,
    truncation::Integer = 31,
    eps1 = 1.75e-3, f_eps = 1, rh_offset_entrainment = 1.3, rh_offset_detrainment = 1.6,
    f_scale_exponent::Integer = 3, delta1 = 0.75e-4,
    gamma = 0.5, f_turbulence = 2, C_d = 0.506, beta = 1.875,
    test_entrainment_factor = 0.4, test_entrainment_exponent::Integer = 3, test_condensate_removal = 0.5,
    delta_T_deep = 0.2, delta_q_deep = 1.0e-4, w_departure = 1,
    mixed_layer_depth = 6000, average_depth = 5000, window_half_width = 2500,
    departure_range = 35000, deep_depth = 20000,
    mstar_base_factor = 0.1, alpha_x_coefficient = 1.66, dx_ref = 125.0e3, tau_min = 720, tau_max = 10800,
    PCAPE_bl = 0, subcloud_exponent::Integer = 1,
    c0 = 1.4e-3, l_crit = 5.0e-4, w_factor = 0.75, w_cap = 10,
    l_threshold_water = 3.0e-4, l_threshold_land = 5.0e-4,
    alpha1 = 5.44e-4, alpha2 = 5.09e-3, alpha3 = 0.5777, C_conv = 0.05,
    RH_cr_water = 0.9, RH_cr_land = 0.7, virtual_coefficient = 0.608, cfl_factor = 1, limiter_timestep_factor::Integer = 2,
    newton_iterations::Integer = 3, rh_min = 0, rh_max = 1, land_threshold = 0.5,
    rain_loading::Bool = false,
) where {NF <: AbstractFloat}
    truncation > 0 || throw(ArgumentError("TB grid spacing needs a positive spectral truncation"))
    PCAPE_bl == 0 || throw(ArgumentError("TB stage 1a fixes PCAPE_bl = 0 (F5)"))
    cfl_factor == 1 || throw(ArgumentError("TB stage 1a uses the strict CFL factor 1 (F4)"))
    limiter_timestep_factor in (1, 2) || throw(ArgumentError(
        "TB D41 limiter_timestep_factor is 2 (host leapfrog) or 1 (forward-Euler driver), contract r15"))
    subcloud_exponent == 1 || throw(ArgumentError("TB stage 1a deep sub-cloud scaling uses m = 1 (D40)"))
    rain_loading && throw(ArgumentError("TB stage 1a contract instantaneous_rain_no_rup_v1 has no rain loading (D28)"))
    f_scale_exponent >= 1 && test_entrainment_exponent >= 1 ||
        throw(ArgumentError("TB entrainment exponents must be positive integers"))
    newton_iterations >= 1 || throw(ArgumentError("TB saturation adjustment needs at least one iteration"))
    grid_spacing = NF(pi) * NF(radius) / NF(truncation)   # D35: dx = pi a / N_trunc
    return TiedtkeBechtoldParameters{NF}(
        NF(gravity), NF(heat_capacity), NF(latent_heat_condensation), NF(R_dry), NF(R_vapor), NF(mol_ratio),
        NF(mu_virt_temp), NF(saturation_vapor_pressure), NF(temperature_freezing),
        NF(eps1), NF(f_eps), NF(rh_offset_entrainment), NF(rh_offset_detrainment), Int(f_scale_exponent),
        NF(delta1),
        NF(gamma), NF(f_turbulence), NF(C_d), NF(beta),
        NF(test_entrainment_factor), Int(test_entrainment_exponent), NF(test_condensate_removal),
        NF(delta_T_deep), NF(delta_q_deep), NF(w_departure),
        NF(mixed_layer_depth), NF(average_depth), NF(window_half_width), NF(departure_range), NF(deep_depth),
        NF(mstar_base_factor), NF(alpha_x_coefficient), NF(dx_ref), NF(tau_min), NF(tau_max),
        NF(PCAPE_bl), Int(subcloud_exponent), grid_spacing,
        NF(c0), NF(l_crit), NF(w_factor), NF(w_cap), NF(l_threshold_water), NF(l_threshold_land),
        NF(alpha1), NF(alpha2), NF(alpha3), NF(C_conv), NF(RH_cr_water), NF(RH_cr_land),
        NF(virtual_coefficient), NF(cfl_factor), Int(limiter_timestep_factor),
        Int(newton_iterations), NF(rh_min), NF(rh_max), NF(land_threshold), rain_loading,
    )
end

"""
    tiedtke_bechtold_host_constants(planet, atmosphere) -> NamedTuple

The host constants of the TB column read from a SpeedyWeather planet and wet atmosphere (for a constructed model:
`model.planet` and `model.atmosphere`), field for field, with no literal copies.
"""
tiedtke_bechtold_host_constants(planet, atmosphere) = (;
    gravity = planet.gravity,
    heat_capacity = atmosphere.heat_capacity,
    latent_heat_condensation = atmosphere.latent_heat_condensation,
    R_dry = atmosphere.R_dry,
    R_vapor = atmosphere.R_vapor,
    mol_ratio = atmosphere.mol_ratio,
    mu_virt_temp = atmosphere.μ_virt_temp,
    saturation_vapor_pressure = atmosphere.saturation_vapor_pressure,
    temperature_freezing = atmosphere.temperature_freezing,
    radius = planet.radius,
)

"""
    tiedtke_bechtold_parameters(NF; planet = SpeedyWeather.Earth(NF),
                                atmosphere = SpeedyWeather.EarthAtmosphere(NF), kwargs...)
    tiedtke_bechtold_parameters(model; kwargs...)

Build `TiedtkeBechtoldParameters` with the host constants of `planet` and `atmosphere`, or of a constructed
SpeedyWeather model (its `planet`, `atmosphere`, number format and spectral truncation). Other keyword arguments
pass through to the constructor (they may also override a host constant explicitly).
"""
function tiedtke_bechtold_parameters(::Type{NF}; planet = SpeedyWeather.Earth(NF),
                                     atmosphere = SpeedyWeather.EarthAtmosphere(NF), kwargs...) where {NF}
    return TiedtkeBechtoldParameters(NF; merge(tiedtke_bechtold_host_constants(planet, atmosphere),
                                               values(kwargs))...)
end

function tiedtke_bechtold_parameters(model::SpeedyWeather.AbstractModel; kwargs...)
    spectral_grid = model.spectral_grid
    return tiedtke_bechtold_parameters(spectral_grid.NF; planet = model.planet, atmosphere = model.atmosphere,
                                       truncation = spectral_grid.trunc, kwargs...)
end

"""
    tiedtke_bechtold_column_work(NF, n)

Host-side allocation of the caller-owned work arrays of one column (half-level arrays of length n+1, layer arrays
of length n). The column call writes every entry before it reads it. After a call the arrays hold the final (closed
and limited) profiles: half-level environment T, q, geopotential and D38-stabilised s, updraught mass flux
(kg m-2 s-1), total-water and liquid-static-energy fluxes (kg m-2 s-1 and W m-2), layer entrainment and turbulent
and organized detrainment (kg m-2 s-1 per layer, D21 top detrainment counted as organized), rain generation G and
rain evaporation e (kg kg-1 s-1) and the detrained-condensate diagnostic (kg kg-1 s-1, re-evaporated in 1a).
"""
function tiedtke_bechtold_column_work(::Type{NF}, n::Integer) where {NF <: AbstractFloat}
    half_level() = zeros(NF, n + 1)
    layer() = zeros(NF, n)
    return (;
        temperature_half = half_level(),
        humidity_half = half_level(),
        geopotential_half = half_level(),
        dry_static_energy_half = half_level(),
        mass_flux_half = half_level(),
        total_water_flux_half = half_level(),
        liquid_static_energy_flux_half = half_level(),
        entrainment = layer(),
        turbulent_detrainment = layer(),
        organized_detrainment = layer(),
        generation = layer(),
        rain_evaporation = layer(),
        detrained_condensate = layer(),
    )
end

# ----------------------------------------------------------------------------------------------------------------
# Thermodynamics (D01-D06)
# ----------------------------------------------------------------------------------------------------------------

# D02 HOST: the SpeedyWeather 0.21.1 Clausius-Clapeyron saturation vapour pressure over liquid, same operation order.
@inline function tb_saturation_vapor_pressure(T, par::TiedtkeBechtoldParameters)
    e0 = par.saturation_vapor_pressure
    Lv = par.latent_heat_condensation
    Rv = par.R_vapor
    T0 = par.temperature_freezing
    return e0 * exp(Lv / Rv * (inv(T0) - inv(T)))
end

# D02 HOST: q_s = eps e_s(T)/p (liquid at all temperatures, D01).
@inline tb_saturation_humidity(T, pressure, par::TiedtkeBechtoldParameters) =
    par.mol_ratio * tb_saturation_vapor_pressure(T, par) / pressure

# D37: dq_s/dT from de_s/dT = L_v e_s/(R_v T^2).
@inline tb_saturation_humidity_slope(T, q_sat, par::TiedtkeBechtoldParameters) =
    q_sat * par.latent_heat_condensation / (par.R_vapor * T^2)

# D04: TB virtual temperature with 0.608 (IFS 6.56).
@inline tb_virtual_temperature(T, q, par::TiedtkeBechtoldParameters) = T * (1 + par.virtual_coefficient * q)

# D04: host mu = 0.6078 stays in the host hydrostatics only (half-level geopotential, D37).
@inline _tb_host_virtual_temperature(T, q, par::TiedtkeBechtoldParameters) = T * (1 + par.mu_virt_temp * q)

# D05: RH as a fraction at the full level, clipped to [0, 1].
@inline tb_relative_humidity(T, q, pressure, par::TiedtkeBechtoldParameters) =
    clamp(q / tb_saturation_humidity(T, pressure, par), par.rh_min, par.rh_max)

# Integer power by repeated multiplication (x^3 = x*x*x), GPU-safe for the fixed exponents of Sect. 9.
@inline function _tb_ipow(x, n::Int)
    y = x
    for _ in 2:n
        y *= x
    end
    return y
end

# D06: f_scale = (q_sat(T_bar_{k+1/2}, p_{k+1/2}) / q_sat(T_bar_base, p_base))^3, saturation specific humidity.
@inline tb_f_scale(T_half_k, p_half_k, T_half_base, p_half_base, par::TiedtkeBechtoldParameters) =
    _tb_ipow(tb_saturation_humidity(T_half_k, p_half_k, par) /
             tb_saturation_humidity(T_half_base, p_half_base, par), par.f_scale_exponent)

"""
    tb_saturation_adjustment(T, q, l, pressure, parameters) -> (T, q, l)

D03: isobaric Newton adjustment at fixed geopotential that conserves s_l = s - L_v l (so c_p T + L_v q) and
q_t = q + l. Two-way: condenses when supersaturated and evaporates carried l when subsaturated, bounded by l.
A fixed number of iterations (parameters.newton_iterations = 3) keeps the GPU kernel free of data-dependent loops.
"""
@inline function tb_saturation_adjustment(T, q, l, pressure, par::TiedtkeBechtoldParameters)
    lv_over_cp = par.latent_heat_condensation / par.heat_capacity
    for _ in 1:par.newton_iterations
        q_sat = tb_saturation_humidity(T, pressure, par)
        slope = tb_saturation_humidity_slope(T, q_sat, par)
        condensed = (q - q_sat) / (1 + lv_over_cp * slope)
        condensed = max(condensed, -l)        # D03: evaporation is bounded by the carried condensate
        T += lv_over_cp * condensed
        q -= condensed
        l += condensed
    end
    return T, q, l
end

# D04, (6.56): B = T_up (1 + 0.608 q_up - l_up - r_up) - T_bar (1 + 0.608 q_bar), q_e := q_bar. The column always
# passes r_up = 0 (D28 instantaneous_rain_no_rup_v1, rain_loading = false). The r_up argument exists for T11.
@inline tb_buoyancy(T_up, q_up, l_up, r_up, T_env, q_env, par::TiedtkeBechtoldParameters) =
    T_up * (1 + par.virtual_coefficient * q_up - l_up - r_up) - tb_virtual_temperature(T_env, q_env, par)

# D32, (6.22) integrand: beta = (T_v,up - T_bar_v)/T_bar_v - l_up (unloaded T_v, loading through -l_up), minus r_up
# (zero in the column, D28; the argument exists for the T11 counterfactual).
@inline function tb_pcape_integrand(T_up, q_up, l_up, r_up, T_env, q_env, par::TiedtkeBechtoldParameters)
    Tv_env = tb_virtual_temperature(T_env, q_env, par)
    return (tb_virtual_temperature(T_up, q_up, par) - Tv_env) / Tv_env - l_up - r_up
end

"""
    tb_kinetic_energy_step(K_lower, c, b, dz)

D20 tb_ke_v1: K_{k-1/2} = [K_{k+1/2}(1 - c dz) + b dz]/(1 + c dz) with the drag damping c >= 0. If c dz >= 1 the
layer-exact solution K_l exp(-2 c dz) + (b/2c)(1 - exp(-2 c dz)) is used instead (the overshooting guard).
"""
@inline function tb_kinetic_energy_step(K_lower, c, b, dz)
    cdz = c * dz
    if cdz < 1
        return (K_lower * (1 - cdz) + b * dz) / (1 + cdz)
    else
        x = 2 * cdz
        return K_lower * exp(-x) + b / (2 * c) * (-expm1(-x))
    end
end

# D30: land if the host land fraction is at least 0.5.
@inline tb_land_flag(land_fraction, par::TiedtkeBechtoldParameters) = land_fraction >= par.land_threshold

# ----------------------------------------------------------------------------------------------------------------
# Trigger: pressure gate and deep test parcel (D07-D16)
# ----------------------------------------------------------------------------------------------------------------

"""
    tb_departure_bounds(p_full, p_half, parameters) -> (first, last)

D07 pressure gate: the admissible set A = {k <= n-2 : p_s - p_k <= 35000 Pa} is scanned k = first, first-1, ...,
last with first = n-2. Levels n and n-1 are never departures. The scan stops at the first k outside A. Level 1 is
never a departure because the first step would need the environment at half level 1/2 (not in A on any real grid).
An empty set gives last = first + 1.
"""
@inline function tb_departure_bounds(p_full, p_half, par::TiedtkeBechtoldParameters)
    n = length(p_full)
    surface_pressure = p_half[n + 1]
    first = n - 2
    last = first + 1
    k = first
    @inbounds while k >= 2 && surface_pressure - p_full[k] <= par.departure_range
        last = k
        k -= 1
    end
    return first, last
end

# D14: deep iff the test cloud (base to top) is thicker than 20000 Pa.
@inline tb_is_deep_cloud(p_base, p_top, par::TiedtkeBechtoldParameters) = p_base - p_top > par.deep_depth

"""
    tb_departure_window(p_full, p_half, k, parameters) -> (mixed_branch, lower, upper)

D09: mixed-layer branch iff p_s - p_k <= 6000 Pa. Window [p_k - 2500, p_k + 2500] Pa, or [p_s - 5000, p_s] Pa if
p_k + 2500 > p_s. (lower, upper) are pressures, lower < upper.
"""
@inline function tb_departure_window(p_full, p_half, k, par::TiedtkeBechtoldParameters)
    n = length(p_full)
    @inbounds surface_pressure = p_half[n + 1]
    @inbounds p_k = p_full[k]
    mixed = surface_pressure - p_k <= par.mixed_layer_depth
    lower = p_k - par.window_half_width
    upper = p_k + par.window_half_width
    if upper > surface_pressure
        lower = surface_pressure - par.average_depth
        upper = surface_pressure
    end
    return mixed, lower, upper
end

# D09: |W intersect [p_{j-1/2}, p_{j+1/2}]|, the pressure overlap of the window with layer j.
@inline _tb_overlap(lower, upper, layer_top, layer_bottom) =
    max(zero(lower), min(upper, layer_bottom) - max(lower, layer_top))

"""
    tb_departure_source(T, q, geopotential, p_full, p_half, k, parameters) -> (mixed_branch, s_up, q_up)
    tb_departure_source(T, q, geopotential, p_full, p_half, k, delta_T, delta_q, parameters)

D08, D09: departure values of the test parcel from level k. Mixed-layer branch: s_tilde + c_p dT and
q_tilde + dq, with s_tilde = sum_j s_bar_j |W intersect layer j| / 5000 Pa (same for q). Otherwise s_bar_k + c_p dT
and q_bar_k + dq. dT = 0.2 K and dq = 1e-4 kg/kg at every stage 1a departure (first form). The second form takes the
perturbation explicitly: stage 1B passes the D46 surface-parcel excess at k = n-1 (D59, D60).
"""
@inline tb_departure_source(T, q, geopotential, p_full, p_half, k, par::TiedtkeBechtoldParameters) =
    tb_departure_source(T, q, geopotential, p_full, p_half, k, par.delta_T_deep, par.delta_q_deep, par)

@inline function tb_departure_source(T, q, geopotential, p_full, p_half, k, delta_T, delta_q,
                                     par::TiedtkeBechtoldParameters)
    cp = par.heat_capacity
    mixed, lower, upper = tb_departure_window(p_full, p_half, k, par)
    if mixed
        s_sum = zero(lower)
        q_sum = zero(lower)
        @inbounds for j in eachindex(p_full)
            weight = _tb_overlap(lower, upper, p_half[j], p_half[j + 1])
            s_sum += (cp * T[j] + geopotential[j]) * weight
            q_sum += q[j] * weight
        end
        s_source = s_sum / par.average_depth
        q_source = q_sum / par.average_depth
    else
        @inbounds s_source = cp * T[k] + geopotential[k]
        @inbounds q_source = q[k]
    end
    return mixed, s_source + cp * delta_T, q_source + delta_q
end

"""
    tb_test_parcel_step(s, q, l, s_env, q_env, eps_ini, dz, geopotential_up, p_up, parameters)

One test-parcel layer step. D10: phi_{k-1/2} = phi_bar_k + (phi_{k+1/2} - phi_bar_k) exp(-eps_ini dz) for s, q and
l (l_bar = 0). D13: saturation is tested before adjustment. D03: adjustment at p_{k-1/2}. D11: the total carried
condensate is halved after the adjustment (l_out = 0.5 (l_in + dl)). The removed water is discarded.
"""
@inline function tb_test_parcel_step(s, q, l, s_env, q_env, eps_ini, dz, geopotential_up, p_up,
                                     par::TiedtkeBechtoldParameters)
    cp = par.heat_capacity
    decay = exp(-eps_ini * dz)
    s_mixed = s_env + (s - s_env) * decay
    q_mixed = q_env + (q - q_env) * decay
    l_mixed = l * decay
    T_mixed = (s_mixed - geopotential_up) / cp
    saturated = q_mixed >= tb_saturation_humidity(T_mixed, p_up, par)
    T_adjusted, q_adjusted, l_adjusted = tb_saturation_adjustment(T_mixed, q_mixed, l_mixed, p_up, par)
    l_out = par.test_condensate_removal * l_adjusted          # D11 (a)
    s_out = cp * T_adjusted + geopotential_up
    return (; s = s_out, q = q_adjusted, l = l_out, T = T_adjusted, saturated,
              l_mixed, q_mixed, l_adjusted)
end

# D07-D15: one deep test-parcel ascent from departure level k_dep. Returns the cloud base (D13), the cloud top
# (D14: highest half level with K > 0) and the parcel state at the base (D15). The first form uses the (6.21)
# perturbation of every stage 1a departure; the second takes it explicitly (stage 1B, D59 at k = n-1).
@inline _tb_test_parcel_ascent(k_dep, T, q, p_full, p_half, geopotential, T_half, q_half, geopotential_half,
                               par::TiedtkeBechtoldParameters) =
    _tb_test_parcel_ascent(k_dep, par.delta_T_deep, par.delta_q_deep, T, q, p_full, p_half, geopotential,
                           T_half, q_half, geopotential_half, par)

@inline function _tb_test_parcel_ascent(k_dep, delta_T, delta_q, T, q, p_full, p_half, geopotential,
                                        T_half, q_half, geopotential_half, par::TiedtkeBechtoldParameters)
    NF = eltype(T)
    n = length(T)
    g = par.gravity
    cp = par.heat_capacity
    ke_buoyancy_factor = g / (par.f_turbulence * (1 + par.gamma))
    drag_factor = 1 + par.beta * par.C_d
    mixed, s, qp = tb_departure_source(T, q, geopotential, p_full, p_half, k_dep, delta_T, delta_q, par)
    l = zero(NF)
    K = par.w_departure^2 / 2                                  # D08, D12: K = 0.5 (1 m/s)^2 at departure
    @inbounds q_sat_surface = tb_saturation_humidity(T[n], p_full[n], par)   # D10: T_bar_surf, p_surf := level n
    @inbounds begin
        T_parcel = (s - geopotential[k_dep]) / cp
        b_lower = tb_buoyancy(T_parcel, qp, l, zero(NF), T[k_dep], q[k_dep], par) /
                  tb_virtual_temperature(T[k_dep], q[k_dep], par)
    end
    has_base = false
    j_base = 0
    j_top = 0
    s_base = zero(NF)
    q_base = zero(NF)
    l_base = zero(NF)
    K_base = zero(NF)
    @inbounds for k in k_dep:-1:2
        j_up = k                                                # half level k-1/2
        # D10: the first step runs from p_{k_dep} to p_{k_dep-1/2}, later steps through whole layers.
        dz = k == k_dep ? (geopotential_half[k] - geopotential[k]) / g :
                          (geopotential_half[k] - geopotential_half[k + 1]) / g
        s_env = cp * T[k] + geopotential[k]
        # D10, U3: eps_ini = 0.4 eps^(1) (q_sat(T_bar)/q_sat(T_bar_surf))^3, numerator at the layer's full level
        eps_ini = par.test_entrainment_factor * par.eps1 *
                  _tb_ipow(tb_saturation_humidity(T[k], p_full[k], par) / q_sat_surface,
                           par.test_entrainment_exponent)
        step = tb_test_parcel_step(s, qp, l, s_env, q[k], eps_ini, dz,
                                   geopotential_half[j_up], p_half[j_up], par)
        # D12, D20: tb_ke_v1 with mu/M := eps_ini and the loaded buoyancy (6.56) of the carried state
        b_upper = tb_buoyancy(step.T, step.q, step.l, zero(NF), T_half[j_up], q_half[j_up], par) /
                  tb_virtual_temperature(T_half[j_up], q_half[j_up], par)
        K = tb_kinetic_energy_step(K, drag_factor * eps_ini, ke_buoyancy_factor * (b_upper + b_lower) / 2, dz)
        s, qp, l = step.s, step.q, step.l
        b_lower = b_upper
        K > 0 || break                                          # w vanishes: the parcel stops
        if !has_base
            if step.saturated                                   # D13: saturated before adjustment, K > 0
                has_base = true
                j_base = j_up
                j_top = j_up
                s_base, q_base, l_base, K_base = s, qp, l, K
            end
        else
            j_top = j_up                                        # D14: highest half level with K > 0
        end
    end
    return (; has_base, j_base, j_top, s_base, q_base, l_base, K_base, mixed)
end

# ----------------------------------------------------------------------------------------------------------------
# Updraught layer (D17-D27)
# ----------------------------------------------------------------------------------------------------------------

"""
    tb_updraught_layer(M_lower, s_lower, q_lower, l_lower, K_lower, B_lower, lfc_passed, rh, f_scale, dz,
                       T_env_lower, q_env_lower, geopotential_lower, T_env_upper, q_env_upper, geopotential_upper,
                       p_upper, land_flag, parameters)

One full-updraught layer k, from half level k+1/2 ("lower") to k-1/2 ("upper"), for the unit mass flux. Order
(D22, Sect. 8): E and D1 (D17, D18), mixing (6.53), adjustment (D03), conversion (D25-D27), B and K (D04, D20), then
the top (D21) or organized detrainment (D19). Mass exchanges E, D1, D2 and the conversion G are layer-integrated
(kg m-2 s-1). `top = true` means K_{k-1/2} <= 0 (D21): the provisional ascent is discarded by the caller.
"""
@inline function tb_updraught_layer(M_lower, s_lower, q_lower, l_lower, K_lower, B_lower, lfc_passed, rh,
                                    f_scale, dz, T_env_lower, q_env_lower, geopotential_lower,
                                    T_env_upper, q_env_upper, geopotential_upper, p_upper, land_flag::Bool,
                                    par::TiedtkeBechtoldParameters)
    NF = typeof(M_lower)
    cp = par.heat_capacity
    # D17: E dp/g = eps^(1) f_eps (1.3 - RH) f_scale M_{k+1/2} dz if B_{k+1/2} > 0, else 0
    entrainment_rate = B_lower > 0 ?
        par.eps1 * par.f_eps * (par.rh_offset_entrainment - rh) * f_scale : zero(NF)
    # D18: D1 dp/g = delta^(1) (1.6 - RH) M_{k+1/2} dz in every cloud layer
    detrainment_rate = par.delta1 * (par.rh_offset_detrainment - rh)
    E = entrainment_rate * M_lower * dz
    D1 = min(detrainment_rate * M_lower * dz, M_lower + E)     # guard: never more than the available mass
    M_mixed = M_lower + E - D1                                 # M' (D23)
    # D22, (6.53): M' phi' = M_{k+1/2} phi_{k+1/2} + E phi_bar_{k+1/2} - D1 phi_{k+1/2}, l_bar = 0
    s_env_lower = cp * T_env_lower + geopotential_lower
    if M_mixed > 0
        s_mixed = ((M_lower - D1) * s_lower + E * s_env_lower) / M_mixed
        q_mixed = ((M_lower - D1) * q_lower + E * q_env_lower) / M_mixed
        l_mixed = (M_lower - D1) * l_lower / M_mixed
    else
        s_mixed, q_mixed, l_mixed = s_lower, q_lower, l_lower
    end
    # Dry adiabatic ascent conserves s, then D03 at p_{k-1/2}
    T_mixed = (s_mixed - geopotential_upper) / cp
    T_adjusted, q_adjusted, l_adjusted = tb_saturation_adjustment(T_mixed, q_mixed, l_mixed, p_upper, par)
    condensation = M_mixed * (q_mixed - q_adjusted)           # D23: c_up dp/g with M' (diagnostic only)
    # D25-D27: conversion with a frozen at the entering mixed l_0 = l_mixed and b dz = l_adjusted - l_mixed
    # (= q_up - q_hat_up, RE-fix of the printed b), w = min(sqrt(2 K_{k+1/2}), 10 m/s); D26 not applied.
    w = min(sqrt(2 * K_lower), par.w_cap)
    threshold = land_flag ? par.l_threshold_land : par.l_threshold_water   # D25, D30
    if l_mixed > threshold && w > 0
        a = par.c0 / (par.w_factor * w) * (1 - exp(-(l_mixed / par.l_crit)^2))   # m-1
    else
        a = zero(NF)
    end
    x = a * dz
    if x > 0
        l_analytic = l_mixed * exp(-x) + (l_adjusted - l_mixed) * (-expm1(-x) / x)
        l_new = clamp(l_analytic, zero(NF), l_adjusted)        # D27: removal cannot exceed the available l
    else
        l_new = l_adjusted
    end
    G = M_mixed * (l_adjusted - l_new)                         # D27: G dp/g = M'(l_0 + b dz - l), >= 0
    s_upper = cp * T_adjusted + geopotential_upper
    # D04, D20 tb_ke_v1: b from the loaded buoyancy (r_up = 0), c = (1 + beta C_d) max(e, d) without D2
    B_upper = tb_buoyancy(T_adjusted, q_adjusted, l_new, zero(NF), T_env_upper, q_env_upper, par)
    b = par.gravity / (par.f_turbulence * (1 + par.gamma)) *
        (B_upper / tb_virtual_temperature(T_env_upper, q_env_upper, par) +
         B_lower / tb_virtual_temperature(T_env_lower, q_env_lower, par)) / 2
    c = (1 + par.beta * par.C_d) * max(entrainment_rate, detrainment_rate)
    K_upper = tb_kinetic_energy_step(K_lower, c, b, dz)
    top = !(K_upper > 0) || !(M_mixed > 0)                     # D21 (and the mass guard)
    # D19 tb_orgdet_v1: active iff the LFC is passed (B > 0 at or below k+1/2), B_{k+1/2} <= 0 and K falls.
    if !top && lfc_passed && !(B_lower > 0) && K_upper < K_lower
        R = (par.rh_offset_detrainment - rh) * sqrt(K_lower / K_upper)
        M_upper = M_mixed / max(R, one(NF))
    else
        M_upper = M_mixed
    end
    D2 = M_mixed - M_upper                                     # 0 <= D2 <= M'
    return (; top, E, D1, D2, M_mixed, M_upper, s = s_upper, q = q_adjusted, l = l_new, T = T_adjusted,
              K = K_upper, B = B_upper, G, condensation, entrainment_rate, detrainment_rate)
end

# D29 (6.44): e = C_conv alpha_1 (RH_cr q_sat - q)[(sqrt(p/p_surf)/alpha_2)(P/C_conv)]^alpha_3, P in kg m-2 s-1
# entering the layer, e in kg kg-1 s-1. Unlimited rate.
@inline function tb_rain_evaporation_rate(P, deficit, pressure, surface_pressure, par::TiedtkeBechtoldParameters)
    x = sqrt(pressure / surface_pressure) / par.alpha2 * (P / par.C_conv)
    return par.C_conv * par.alpha1 * deficit * x^par.alpha3
end

"""
    tb_rain_evaporation_mass(P, q_env, q_sat, dp, dt, land_flag, pressure, surface_pressure, parameters)

D29 limited sub-cloud evaporation of one layer, layer-integrated (kg m-2 s-1): active only where
q_bar < RH_cr q_sat (RH_cr 0.9 water, 0.7 land), limited so that e dp/g <= P and e dt <= RH_cr q_sat - q_bar.
"""
@inline function tb_rain_evaporation_mass(P, q_env, q_sat, dp, dt, land_flag::Bool, pressure,
                                          surface_pressure, par::TiedtkeBechtoldParameters)
    rh_cr = land_flag ? par.RH_cr_land : par.RH_cr_water
    deficit = rh_cr * q_sat - q_env
    (deficit > 0 && P > 0) || return zero(P)
    g = par.gravity
    e = tb_rain_evaporation_rate(P, deficit, pressure, surface_pressure, par)
    return min(e * dp / g, P, deficit * dp / (g * dt))
end

# ----------------------------------------------------------------------------------------------------------------
# The column (Sect. 8 order)
# ----------------------------------------------------------------------------------------------------------------

@inline function _tb_result(::Type{NF}, rain_surface, top_level, base_half_level, convection_type,
                            closure_factor, f_CFL, f_pos, final_factor, PCAPE, tau, base_mass_flux,
                            unit_base_mass_flux, departure_level, departures_scanned, departure_first,
                            departure_last, mixed_layer_departure, test_parcel_base_w,
                            updraught_base_kinetic_energy, deep_test_depth, max_test_cloud_depth,
                            model_top_reached, condensation_column) where {NF}
    return (;
        rain_mass_flux_surface = NF(rain_surface),       # kg m-2 s-1 = rain_mass_flux_half[n+1]
        tb_top_full_level = Int(top_level),              # D21 top layer k; n+1 = none (D43)
        cloud_base_half_level = Int(base_half_level),    # j with p_base = p_half[j]; 0 = none
        convection_type = Int(convection_type),          # 0 none, 1 deep
        closure_factor = NF(closure_factor),             # M_base/M*_base before the limiters (D31)
        f_CFL = NF(f_CFL),                               # D41, min(1, ...)
        f_pos = NF(f_pos),                               # D41, min(1, ...)
        final_factor = NF(final_factor),                 # closure_factor * min(1, f_CFL, f_pos) (D36, D41)
        PCAPE = NF(PCAPE),                               # J m-3, unit-flux plume (D32)
        tau = NF(tau),                                   # s (D35); 0 when no cloud layer
        cloud_base_mass_flux = NF(base_mass_flux),       # kg m-2 s-1, final
        unit_cloud_base_mass_flux = NF(unit_base_mass_flux),  # M*_base (D34)
        departure_level = Int(departure_level),          # winning departure k; 0 = none
        departures_scanned = Int(departures_scanned),
        departure_first = Int(departure_first),          # first scanned k (n-2)
        departure_last = Int(departure_last),            # last scanned k
        mixed_layer_departure = Bool(mixed_layer_departure),
        test_parcel_base_w = NF(test_parcel_base_w),     # sqrt(2 K_test) at the cloud base (D15)
        updraught_base_kinetic_energy = NF(updraught_base_kinetic_energy),  # full-plume K_base (D15)
        deep_test_depth = NF(deep_test_depth),           # p_base - p_top of the winning test parcel
        max_test_cloud_depth = NF(max_test_cloud_depth), # largest test-cloud depth over the scanned departures
        model_top_reached = Bool(model_top_reached),     # D21 flag: the plume reached layer 1
        condensation_column = NF(condensation_column),   # sum of c_up dp/g, final (D23 diagnostic)
    )
end

"""
    tiedtke_bechtold_column!(dT, dq, rain_mass_flux_half, work, T, q, p_full, p_half, geopotential, dt,
                             land_flag, latitude, parameters) -> NamedTuple

Stage 1a Tiedtke-Bechtold column (DESIGN.md "Column interface v1", normative signature r2). Inputs (read only):
`T[k]` (K), `q[k]` (kg/kg), `p_full[k]` (Pa), `p_half[j]` (Pa, j = 1..n+1), `geopotential[k]` (m2 s-2), `dt` (s),
`land_flag` (D30), `latitude` (diagnostics only, unused). Outputs written into the caller's arrays: `dT[k]` (K s-1),
`dq[k]` (kg kg-1 s-1) and `rain_mass_flux_half[j]` (kg m-2 s-1, downward, j = 1..n+1). The caller-owned `work`
arrays (see `tiedtke_bechtold_column_work`) are fully written by the call and hold the final profiles afterwards.
Returns an isbits NamedTuple of scalar diagnostics. No allocation, no exceptions, no dynamic dispatch.
"""
function tiedtke_bechtold_column!(
    dT::AbstractVector{NF}, dq::AbstractVector{NF}, rain_mass_flux_half::AbstractVector{NF}, work::NamedTuple,
    T::AbstractVector{NF}, q::AbstractVector{NF}, p_full::AbstractVector{NF}, p_half::AbstractVector{NF},
    geopotential::AbstractVector{NF}, dt::Real, land_flag::Bool, latitude::Real,
    par::TiedtkeBechtoldParameters{NF},
) where {NF <: AbstractFloat}
    n = length(T)
    g = par.gravity
    cp = par.heat_capacity
    Lv = par.latent_heat_condensation
    timestep = convert(NF, dt)                                  # D42, D43: host step (900 s)
    z = zero(NF)
    (; temperature_half, humidity_half, geopotential_half, dry_static_energy_half, mass_flux_half,
       total_water_flux_half, liquid_static_energy_flux_half, entrainment, turbulent_detrainment,
       organized_detrainment, generation, rain_evaporation, detrained_condensate) = work

    # Every output and work entry is written before any read (scratch can hold anything, including NaN).
    @inbounds for k in 1:n
        dT[k] = z
        dq[k] = z
        entrainment[k] = z
        turbulent_detrainment[k] = z
        organized_detrainment[k] = z
        generation[k] = z
        rain_evaporation[k] = z
        detrained_condensate[k] = z
    end
    @inbounds for j in 1:(n + 1)
        rain_mass_flux_half[j] = z
        mass_flux_half[j] = z
        total_water_flux_half[j] = z
        liquid_static_energy_flux_half[j] = z
    end

    # Step 1 (Sect. 8). D37: half level k+1/2 extrapolated from full level k along h_sat. Half level 1/2 has no
    # full level above it, so it gets a placeholder copy of level 1 and is never read by the physics.
    surface_pressure = p_half[n + 1]
    @inbounds begin
        temperature_half[1] = T[1]
        humidity_half[1] = q[1]
        geopotential_half[1] = geopotential[1]
        for k in 1:n
            q_sat = tb_saturation_humidity(T[k], p_full[k], par)
            slope = tb_saturation_humidity_slope(T[k], q_sat, par)
            dTdp = (par.R_dry * tb_virtual_temperature(T[k], q[k], par) / p_full[k] + Lv * q_sat / p_full[k]) /
                   (cp + Lv * slope)
            dqdp = slope * dTdp - q_sat / p_full[k]
            dp_half = p_half[k + 1] - p_full[k]
            temperature_half[k + 1] = T[k] + dTdp * dp_half
            humidity_half[k + 1] = q[k] + dqdp * dp_half
        end
        # D37: half-level geopotential from the host hydrostatic discretisation (host mu, D04)
        for j in 2:n
            geopotential_half[j] = geopotential[j] +
                par.R_dry * _tb_host_virtual_temperature(T[j], q[j], par) * log(p_full[j] / p_half[j])
        end
        geopotential_half[n + 1] = geopotential[n] -
            par.R_dry * _tb_host_virtual_temperature(T[n], q[n], par) * log(p_half[n + 1] / p_full[n])
        # D38: bottom-up s_bar_{k-1/2} := max(s_bar_{k-1/2}, s_bar_{k+1/2}) on the half levels the fluxes use
        for j in 2:(n + 1)
            dry_static_energy_half[j] = cp * temperature_half[j] + geopotential_half[j]
        end
        dry_static_energy_half[1] = dry_static_energy_half[2]          # placeholder, never read
        for j in n:-1:2
            dry_static_energy_half[j] = max(dry_static_energy_half[j], dry_static_energy_half[j + 1])
        end
    end

    # Step 2. D07-D16: scan the departures bottom up; the first deep test cloud wins.
    first, last = tb_departure_bounds(p_full, p_half, par)
    scanned = 0
    found = false
    k_dep = 0
    j_base = 0
    j_top_test = 0
    mixed_departure = false
    s_base = z
    q_base = z
    l_base = z
    K_base = z
    max_depth = z
    deep_depth_found = z
    k = first
    @inbounds while k >= last && !found
        scanned += 1
        ascent = _tb_test_parcel_ascent(k, T, q, p_full, p_half, geopotential, temperature_half, humidity_half,
                                        geopotential_half, par)
        if ascent.has_base
            depth = p_half[ascent.j_base] - p_half[ascent.j_top]
            max_depth = max(max_depth, depth)
            if tb_is_deep_cloud(p_half[ascent.j_base], p_half[ascent.j_top], par)
                found = true
                k_dep = k
                j_base = ascent.j_base
                j_top_test = ascent.j_top
                mixed_departure = ascent.mixed
                s_base, q_base, l_base, K_base = ascent.s_base, ascent.q_base, ascent.l_base, ascent.K_base
                deep_depth_found = depth
            end
        end
        k -= 1
    end
    scan_last = scanned > 0 ? first - scanned + 1 : first + 1

    # Step 3. D16: no deep cloud gives exact zero tendencies and rain (already written above).
    if !found
        return _tb_result(NF, z, n + 1, 0, 0, z, one(NF), one(NF), z, z, z, z, z, 0, scanned, first, scan_last,
                          false, z, z, z, max_depth, false, z)
    end

    # Step 4. Full updraught from the base half level j_base for the unit mass flux.
    # D34: M*_base = 0.1 dp_base/(g timestep), dp_base of the first cloud layer (layer j_base - 1).
    unit_base_mass_flux = par.mstar_base_factor * (p_half[j_base] - p_half[j_base - 1]) / (g * timestep)
    p_base = p_half[j_base]
    M = unit_base_mass_flux
    # D15: test-parcel w, s, q and l at the cloud base start the full updraught (K_base = 0.5 w_test^2).
    s_up = s_base
    q_up = q_base
    l_up = l_base
    K_up = K_base
    test_parcel_base_w = sqrt(2 * K_base)
    @inbounds begin
        T_up = (s_up - geopotential_half[j_base]) / cp
        B_lower = tb_buoyancy(T_up, q_up, l_up, z, temperature_half[j_base], humidity_half[j_base], par)
        beta_lower = tb_pcape_integrand(T_up, q_up, l_up, z, temperature_half[j_base], humidity_half[j_base], par)
        mass_flux_half[j_base] = M
        # D39 conserved-variable fluxes F_qt = M(q_up + l_up - q_bar), F_sl = M(s_up - L_v l_up - s_bar)
        total_water_flux_half[j_base] = M * (q_up + l_up - humidity_half[j_base])
        liquid_static_energy_flux_half[j_base] = M * (s_up - Lv * l_up - dry_static_energy_half[j_base])
        q_sat_base = tb_saturation_humidity(temperature_half[j_base], p_base, par)   # D06 denominator
    end
    lfc_passed = B_lower > 0
    PCAPE = z
    Den = z
    w_dz_sum = z
    condensation_column = z
    k_top = 0
    model_top = false
    @inbounds for k in (j_base - 1):-1:1
        dp = p_half[k + 1] - p_half[k]
        if k == 1
            # D21 (r14) at the model top: layer 1 is never integrated (no D37 environment at half level 1/2), so a
            # plume reaching half level 3/2 with K > 0 has entered the top layer: layer 1 is the top layer and the
            # model_top ("entered top layer") flag is set.
            organized_detrainment[1] = M
            detrained_condensate[1] = g / dp * M * l_up
            k_top = 1
            model_top = true
            break
        end
        dz = (geopotential_half[k] - geopotential_half[k + 1]) / g
        rh = tb_relative_humidity(T[k], q[k], p_full[k], par)                       # D05
        f_scale = _tb_ipow(tb_saturation_humidity(temperature_half[k + 1], p_half[k + 1], par) / q_sat_base,
                           par.f_scale_exponent)                                     # D06
        layer = tb_updraught_layer(M, s_up, q_up, l_up, K_up, B_lower, lfc_passed, rh, f_scale, dz,
                                   temperature_half[k + 1], humidity_half[k + 1], geopotential_half[k + 1],
                                   temperature_half[k], humidity_half[k], geopotential_half[k], p_half[k],
                                   land_flag, par)
        if layer.top
            # D21: provisional ascent discarded, M_{k-1/2} = 0, G_k = 0, layer k detrains all of M_{k+1/2}
            organized_detrainment[k] = M
            detrained_condensate[k] = g / dp * M * l_up
            k_top = k
            break
        end
        entrainment[k] = layer.E
        turbulent_detrainment[k] = layer.D1
        organized_detrainment[k] = layer.D2
        generation[k] = g * layer.G / dp                                             # kg kg-1 s-1 (D25)
        # Detrained condensate: D1 leaves with the k+1/2 properties, D2 with the adjusted k-1/2 properties.
        detrained_condensate[k] = g / dp * (layer.D1 * l_up + layer.D2 * layer.l)
        condensation_column += layer.condensation
        # Step 5 integrands over the cloud layers between the base and top half levels.
        beta_upper = tb_pcape_integrand(layer.T, layer.q, layer.l, z, temperature_half[k], humidity_half[k], par)
        PCAPE += (beta_lower + beta_upper) / 2 * dp                                  # D32
        Tv_upper = tb_virtual_temperature(temperature_half[k], humidity_half[k], par)
        Tv_lower = tb_virtual_temperature(temperature_half[k + 1], humidity_half[k + 1], par)
        Den += g / tb_virtual_temperature(T[k], q[k], par) * (M + layer.M_upper) / 2 *
               (Tv_upper - Tv_lower + (geopotential_half[k] - geopotential_half[k + 1]) / cp)   # D33
        w_dz_sum += (sqrt(2 * K_up) + sqrt(2 * layer.K)) / 2 * dz                   # D35 w_bar
        M = layer.M_upper
        s_up, q_up, l_up, K_up = layer.s, layer.q, layer.l, layer.K
        B_lower = layer.B
        beta_lower = beta_upper
        lfc_passed = lfc_passed || layer.B > 0
        mass_flux_half[k] = M
        total_water_flux_half[k] = M * (q_up + l_up - humidity_half[k])
        liquid_static_energy_flux_half[k] = M * (s_up - Lv * l_up - dry_static_energy_half[k])
    end
    j_top = k_top + 1                                                                # cloud-top half level

    # D40 (contract r13, a stage-1a adaptation): the unit cloud-base NET conserved-variable fluxes are scaled by
    # Z_j = (p_s - p_j)/(p_s - p_base) (m = 1, linear in pressure) to zero at the surface, F(j) = Z_j F(j_base), and
    # M(j) = Z_j M_base. With the elevated n-2 departure, forming the net flux against the local half-level
    # environment (r12) could reverse its sign near a moist surface and moisten layer n (the C08 failure).
    @inbounds begin
        total_water_base_flux = total_water_flux_half[j_base]
        liquid_static_energy_base_flux = liquid_static_energy_flux_half[j_base]
        for j in (j_base + 1):n
            Z = (surface_pressure - p_half[j]) / (surface_pressure - p_base)
            mass_flux_half[j] = Z * unit_base_mass_flux
            total_water_flux_half[j] = Z * total_water_base_flux
            liquid_static_energy_flux_half[j] = Z * liquid_static_energy_base_flux
        end
        mass_flux_half[n + 1] = z
        total_water_flux_half[n + 1] = z
        liquid_static_energy_flux_half[n + 1] = z
    end

    # Step 5. D35 tau and D31 closure (PCAPE_bl = 0, F5).
    @inbounds H = (geopotential_half[j_top] - geopotential_half[j_base]) / g
    tau = z
    closure_factor = z
    if H > 0 && w_dz_sum > 0
        w_bar = w_dz_sum / H
        alpha_x = 1 + par.alpha_x_coefficient * par.grid_spacing / par.dx_ref
        tau = clamp(H / w_bar * alpha_x, par.tau_min, par.tau_max)
        if PCAPE - par.PCAPE_bl > 0 && Den > 0
            closure_factor = (PCAPE - par.PCAPE_bl) / (tau * Den)
        end
    end

    # Step 6. D41 limiters on the closed unit-flux profiles: one column factor f = min(1, f_CFL, f_pos).
    # Contract r15: both limits use the step the host applies, limiter_timestep_factor * dt (2 dt for the leapfrog).
    limiter_step = par.limiter_timestep_factor * timestep
    f_CFL = one(NF)
    f_pos = one(NF)
    if closure_factor > 0
        @inbounds for j in 2:n
            M_closed = closure_factor * mass_flux_half[j]
            if M_closed > 0
                dp_j = min(p_half[j] - p_half[j - 1], p_half[j + 1] - p_half[j])
                f_CFL = min(f_CFL, par.cfl_factor * dp_j / (g * limiter_step * M_closed))
            end
        end
        @inbounds for k in 1:n
            dp = p_half[k + 1] - p_half[k]
            dq_closed = closure_factor *
                (g / dp * (total_water_flux_half[k + 1] - total_water_flux_half[k]) - generation[k])
            if dq_closed < 0
                f_pos = min(f_pos, max(q[k], z) / (-dq_closed * limiter_step))
            end
        end
    end
    final_factor = closure_factor * min(f_CFL, f_pos)

    if !(final_factor > 0)
        # D31 zero closure (or a zero limiter): exact zero outputs and work profiles. No convection acts, so the
        # type is 0 and the top is n+1 (D43 "none"); the trigger diagnostics (departure, base, PCAPE) are kept.
        @inbounds for j in 1:(n + 1)
            mass_flux_half[j] = z
            total_water_flux_half[j] = z
            liquid_static_energy_flux_half[j] = z
        end
        @inbounds for k in 1:n
            entrainment[k] = z
            turbulent_detrainment[k] = z
            organized_detrainment[k] = z
            generation[k] = z
            detrained_condensate[k] = z
        end
        return _tb_result(NF, z, n + 1, j_base, 0, closure_factor, f_CFL, f_pos, z, PCAPE, tau, z,
                          unit_base_mass_flux, k_dep, scanned, first, scan_last, mixed_departure,
                          test_parcel_base_w, K_base, deep_depth_found, max_depth, model_top, z)
    end

    # Step 7. D36: mass, F_qt, F_sl and G (and the diagnostics) times the final factor.
    @inbounds for j in 1:(n + 1)
        mass_flux_half[j] *= final_factor
        total_water_flux_half[j] *= final_factor
        liquid_static_energy_flux_half[j] *= final_factor
    end
    @inbounds for k in 1:n
        entrainment[k] *= final_factor
        turbulent_detrainment[k] *= final_factor
        organized_detrainment[k] *= final_factor
        generation[k] *= final_factor
        detrained_condensate[k] *= final_factor
    end

    # Steps 7-8. D28 rain flux from the top down (adds G in cloud layers), D29 sub-cloud evaporation (full levels
    # with p_k > p_base), D39 tendencies: dq = (g/dp)(F_qt,k+1/2 - F_qt,k-1/2) - G + e,
    # dT = [(g/dp)(F_sl,k+1/2 - F_sl,k-1/2) + L_v (G - e)]/c_p (geopotential held).
    P = z
    @inbounds rain_mass_flux_half[1] = z
    @inbounds for k in 1:n
        dp = p_half[k + 1] - p_half[k]
        e_k = z
        if p_full[k] > p_base
            q_sat_k = tb_saturation_humidity(T[k], p_full[k], par)
            evaporated = tb_rain_evaporation_mass(P, q[k], q_sat_k, dp, timestep, land_flag, p_full[k],
                                                  surface_pressure, par)
            e_k = evaporated * g / dp
            P -= evaporated                                    # >= 0 exactly: evaporated <= P
        else
            P += generation[k] * dp / g
        end
        rain_evaporation[k] = e_k
        rain_mass_flux_half[k + 1] = P
        dq[k] = g / dp * (total_water_flux_half[k + 1] - total_water_flux_half[k]) - generation[k] + e_k
        dT[k] = (g / dp * (liquid_static_energy_flux_half[k + 1] - liquid_static_energy_flux_half[k]) +
                 Lv * (generation[k] - e_k)) / cp
    end

    return _tb_result(NF, P, k_top, j_base, 1, closure_factor, f_CFL, f_pos, final_factor, PCAPE, tau,
                      final_factor * unit_base_mass_flux, unit_base_mass_flux, k_dep, scanned, first, scan_last,
                      mixed_departure, test_parcel_base_w, K_base, deep_depth_found, max_depth, model_top,
                      final_factor * condensation_column)
end

# ================================================================================================================
# Stage 1B, column science: the Cy41r2 shallow branch and the surface-flux test parcel
# ================================================================================================================
#
# Frozen contract: artifacts/development_wave_20261103/tiedtke_bechtold_claude_v1/frozen_v5/FROZEN_CONTRACT_v5.sha256
# (sha256 c6327c89e97ee1b4676a0d737f46f07cafa87351622b37fcff7db927284eb898), which pins STAGE1B_CHOICES.md b725f19f
# (r4) and DESIGN_STAGE2_v1.md 30a65da1 on top of the stage 1a manifest above. Rows from D45 on and tests T1B, T16-T29
# are those of STAGE1B_CHOICES.md (D45-D48 restate the stage 1a Sect. 7 rows). This is the scientific contract only
# (F10): the persistent lagged-flux state, its restart and reset immunity, and the host wiring are stage 2.
#
# Selection (F10, D81): stage 1B is the method of `tiedtke_bechtold_column!` with three extra arguments (shallow
# parameters, lagged surface flux, closure). The stage 1a method above is unchanged and stays the default. With the
# lagged state unset, stage 1B calls that stage 1a method itself, so the n = 0 fallback is the stage 1a column bit for
# bit. With the lagged state set, stage 1B runs the Sect. 8 order of STAGE1B_CHOICES.md.
#
# Code layout: the stage 1a column body is not refactored, so its bytes stay frozen. Its blocks are repeated below as
# helpers (`_tb_initialize_column!`, `_tb_deep_updraught!`, `_tb_subcloud_fluxes!`, `_tb_deep_closure`,
# `_tb_limiter_factors`, `_tb_zero_profiles!`, `_tb_scale_and_assemble!`) with the same expressions in the same order.
# test/tiedtke_bechtold_stage1b.jl checks that a stage 1B deep cloud from n-2 equals the stage 1a column bit for bit.

const TIEDTKE_BECHTOLD_STAGE1B_CONTRACT_MANIFEST_SHA256 =
    "c6327c89e97ee1b4676a0d737f46f07cafa87351622b37fcff7db927284eb898"
const TIEDTKE_BECHTOLD_STAGE1B_TRIGGER_VARIANT = "shallow_n_deep_n1_v1"   # D82 proposed name (fixed in stage 2)
const TIEDTKE_BECHTOLD_TRIGGER_PATH_STAGE1A_FALLBACK = 1                   # D81, D83
const TIEDTKE_BECHTOLD_TRIGGER_PATH_STAGE1B = 2                            # D83

"""
Parameters that stage 1B adds (STAGE1B_CHOICES.md Sect. 9). Every other value is read from the stage 1a
`TiedtkeBechtoldParameters`: the U3 deep test parcel, w_dep = 1 m/s at n-1 (D59), the 50% test condensate removal
(D52 decides that the shallow parcel uses the D11 rule, so it shares `test_condensate_removal`), the 200 hPa depth
(D55 is the exact complement of D14), m = 1 (D67 = D40) and the D41 limiter step factor (D76). The struct is isbits.
"""
struct TiedtkeBechtoldShallowParameters{NF <: AbstractFloat}
    # (6.19) excesses (IFS p. 80, D46)
    excess_factor::NF              # 1.5, IFS
    delta_T_cap::NF                # 3 K, IFS
    delta_q_cap::NF                # 2e-3 kg/kg, IFS
    # (6.20) convective velocity scale (IFS p. 80, D47)
    wstar_factor::NF               # 1.2, IFS
    wstar_flux_factor::NF          # 1.5, IFS
    wstar_virtual_coefficient::NF  # 0.61 as printed in (6.20) (not the (6.56) 0.608), IFS
    von_karman::NF                 # kappa = 0.4, IFS
    friction_velocity::NF          # u* = 0.1 m/s, IFS; the argument floor u*^3 is RE (D47)
    # U2 shallow test-parcel entrainment eps = b/z + a (IFS p. 80, Ch. 3 (3.37), D51)
    test_entrainment_b::NF         # 0.8 (z in m), IFS
    test_entrainment_a::NF         # 2e-4 m-1, IFS
    # Shallow updraught (IFS (6.7), (6.9), D62, D63)
    f_eps::NF                      # 2, IFS
    rh_offset_detrainment::NF      # 1.6 in (6.9), IFS
end

"""
    TiedtkeBechtoldShallowParameters(NF; kwargs...)
    tiedtke_bechtold_shallow_parameters(NF; kwargs...)

Stage 1B parameters with the Sect. 9 values. Keyword arguments override single values for sensitivity work only.
"""
function TiedtkeBechtoldShallowParameters(
    ::Type{NF};
    excess_factor = 1.5, delta_T_cap = 3, delta_q_cap = 2.0e-3,
    wstar_factor = 1.2, wstar_flux_factor = 1.5, wstar_virtual_coefficient = 0.61, von_karman = 0.4,
    friction_velocity = 0.1, test_entrainment_b = 0.8, test_entrainment_a = 2.0e-4,
    f_eps = 2, rh_offset_detrainment = 1.6,
) where {NF <: AbstractFloat}
    friction_velocity > 0 ||
        throw(ArgumentError("TB D47 needs u* > 0: the w* argument floor is u*^3"))
    (test_entrainment_b >= 0 && test_entrainment_a >= 0) ||
        throw(ArgumentError("TB D51 test entrainment coefficients must be non-negative"))
    return TiedtkeBechtoldShallowParameters{NF}(
        NF(excess_factor), NF(delta_T_cap), NF(delta_q_cap),
        NF(wstar_factor), NF(wstar_flux_factor), NF(wstar_virtual_coefficient), NF(von_karman),
        NF(friction_velocity), NF(test_entrainment_b), NF(test_entrainment_a),
        NF(f_eps), NF(rh_offset_detrainment),
    )
end

tiedtke_bechtold_shallow_parameters(::Type{NF}; kwargs...) where {NF <: AbstractFloat} =
    TiedtkeBechtoldShallowParameters(NF; kwargs...)

"""
    TiedtkeBechtoldLaggedSurfaceFlux(sensible_heat_flux, humidity_flux)   # set
    TiedtkeBechtoldLaggedSurfaceFlux{NF}()                                 # the declared unset state (D81)

The D45-D48 adapter values that the stage 1B column reads, as plain inputs: the previous-step area-weighted host
sensible heat flux H_host (W m-2) and surface humidity flux F_q,host (kg m-2 s-1), both upward positive (the
SpeedyWeather convention), and the `set` flag. The column applies the D45 signs once (J_s = -H_host,
J_q,IFS = -F_q,host). The persistent snapshot, its timestamp and adapter version, and its restart are stage 2 (D48).
"""
struct TiedtkeBechtoldLaggedSurfaceFlux{NF <: AbstractFloat}
    set::Bool
    sensible_heat_flux::NF   # H_host, W m-2, upward positive
    humidity_flux::NF        # F_q,host, kg m-2 s-1, upward positive
end

TiedtkeBechtoldLaggedSurfaceFlux(sensible_heat_flux::NF, humidity_flux::NF) where {NF <: AbstractFloat} =
    TiedtkeBechtoldLaggedSurfaceFlux{NF}(true, sensible_heat_flux, humidity_flux)
TiedtkeBechtoldLaggedSurfaceFlux{NF}() where {NF <: AbstractFloat} =
    TiedtkeBechtoldLaggedSurfaceFlux{NF}(false, zero(NF), zero(NF))

"""
Shallow closure source S_sub of (6.29) (D68). Two separate variants with separate tests (F10, G6): the HOST closure
(D71) and the SCM closure (D72). The host build constructs only `TiedtkeBechtoldHostClosure`.
"""
abstract type AbstractTiedtkeBechtoldShallowClosure end

"""
    TiedtkeBechtoldHostClosure()

D71 HOST closure: S_sub := F_h,surf = H_host + L_v F_q,host from the lagged adapter, W m-2, upward positive.
Radiation and advection are omitted (declared ReadyESM host approximation). It has no field, so it carries no S_LS.
"""
struct TiedtkeBechtoldHostClosure <: AbstractTiedtkeBechtoldShallowClosure end

"""
    TiedtkeBechtoldSCMClosure(temperature_tendency, humidity_tendency)

D72 SCM closure: S_sub := F_h,surf (lagged, as D71) + S_LS. The two profiles are the SCM harness's non-convective
large-scale tendencies on the host full levels, dT_LS[k] = Pi (dtheta_l/dt)_rad + Pi (-omega dtheta/dp) (K s-1) and
dq_LS[k] = (dq/dt)_adv - omega dq/dp (kg kg-1 s-1), which the harness evaluates from the pre-step state before the TB
call (D72 ORDER). The column computes no forcing: it only integrates the given profiles over the layers below the
cloud base it finds, S_LS = (1/g) sum_{k >= j_base} (c_p dT_LS[k] + L_v dq_LS[k]) Delta p_k (W m-2). This type exists
for the SCM harness path only (G6). CPU only (it holds vectors).
"""
struct TiedtkeBechtoldSCMClosure{V <: AbstractVector} <: AbstractTiedtkeBechtoldShallowClosure
    temperature_tendency::V   # K s-1 on full levels 1..n
    humidity_tendency::V      # kg kg-1 s-1 on full levels 1..n
end
# OPEN(D72): D72 defines S_LS as an integral "over the host layers below j_base", and j_base is found inside the column,
# so the SCM input is read here as the two per-layer tendency profiles, not a precomputed scalar. The harness-side
# terms (the Pi conversion, the omega = -rho g wa subsidence derivatives, the conservative and ln-p remaps and the
# interp_time step-start rule, T22 r4 items (2)-(3) and the time rule) belong to barry_scm_harness_v1 and are not
# column code. They are not implemented or tested here.

"""
    tb_scm_large_scale_source(closure::TiedtkeBechtoldSCMClosure, j_base, p_half, parameters) -> S_LS (W m-2)

D72: S_LS = (1/g) sum over the layers k = j_base..n below the cloud-base half level of (c_p dT_LS + L_v dq_LS) Delta p_k.
"""
@inline function tb_scm_large_scale_source(closure::TiedtkeBechtoldSCMClosure, j_base, p_half,
                                           par::TiedtkeBechtoldParameters)
    n = length(p_half) - 1
    cp = par.heat_capacity
    Lv = par.latent_heat_condensation
    g = par.gravity
    S_LS = zero(eltype(p_half))
    @inbounds for k in j_base:n
        S_LS += (cp * closure.temperature_tendency[k] + Lv * closure.humidity_tendency[k]) *
                (p_half[k + 1] - p_half[k]) / g
    end
    return S_LS
end

# D71, D72: (S_sub, S_LS). The HOST method has no S_LS term at all; the SCM method adds the integrated forcing.
@inline _tb_shallow_closure_source(::TiedtkeBechtoldHostClosure, F_h, j_base, p_half,
                                   par::TiedtkeBechtoldParameters) = (F_h, zero(F_h))
@inline function _tb_shallow_closure_source(closure::TiedtkeBechtoldSCMClosure, F_h, j_base, p_half,
                                            par::TiedtkeBechtoldParameters)
    S_LS = tb_scm_large_scale_source(closure, j_base, p_half, par)
    return F_h + S_LS, S_LS
end

# ----------------------------------------------------------------------------------------------------------------
# Surface parcel (D45-D53)
# ----------------------------------------------------------------------------------------------------------------

"""
    tb_surface_parcel_excess(lagged, T_n, q_n, p_n, z_n, shallow, parameters) -> NamedTuple

D45-D47 at the lowest full level n, with the lagged upward host fluxes H (W m-2) and F_q (kg m-2 s-1):
- D45 signs: J_s = -H and J_q,IFS = J_q/L = -F_q, so -[J_s/c_p + 0.61 T_bar J_q/L] = H/c_p + 0.61 T_bar F_q is the
  upward surface buoyancy flux (K kg m-2 s-1). L_v never enters (J_q is a mass flux here).
- D47 (6.20): w* = 1.2 max(u*^3 + 1.5 (g z kappa/(rho_bar T_bar))(H/c_p + 0.61 T_bar F_q), u*^3)^(1/3), so
  w* >= 1.2 u* (argument floor RE).
- D46 (6.19), OpenIFS-lineage sign: Delta T = min(3, 1.5 H/(rho_bar c_p w*)), Delta q = min(2e-3, 1.5 F_q/(rho_bar w*)),
  caps only, no floors (a downward flux gives a negative excess).
z_n, rho_bar and T_bar are those of level n (D47). `nonupward` marks a zero or downward buoyancy flux (D53 counter).
"""
@inline function tb_surface_parcel_excess(lagged::TiedtkeBechtoldLaggedSurfaceFlux, T_n, q_n, p_n, z_n,
                                          shp::TiedtkeBechtoldShallowParameters, par::TiedtkeBechtoldParameters)
    g = par.gravity
    cp = par.heat_capacity
    H = lagged.sensible_heat_flux
    F_q = lagged.humidity_flux
    # OPEN(D47): the contract fixes rho_bar "of the lowest full level n" but not its virtual-temperature form. Read
    # here as the host state density p_n/(R_d T_n (1 + mu q_n)) with the host mu of the D04 hydrostatics (the same
    # relation that gives z_n). The (6.56) 0.608 form differs by about 1e-6 relative. Reviewer to confirm.
    density = p_n / (par.R_dry * _tb_host_virtual_temperature(T_n, q_n, par))
    buoyancy_flux = H / cp + shp.wstar_virtual_coefficient * T_n * F_q
    ustar_cubed = _tb_ipow(shp.friction_velocity, 3)
    argument = ustar_cubed + shp.wstar_flux_factor * g * z_n * shp.von_karman / (density * T_n) * buoyancy_flux
    wstar = shp.wstar_factor * cbrt(max(argument, ustar_cubed))
    delta_T = min(shp.delta_T_cap, shp.excess_factor * H / (density * cp * wstar))
    delta_q = min(shp.delta_q_cap, shp.excess_factor * F_q / (density * wstar))
    return (; wstar, delta_T, delta_q, density, buoyancy_flux, argument,
              floored = !(argument > ustar_cubed), nonupward = !(buoyancy_flux > 0))
end

# D51 (U2): eps = 0.8/z + 2e-4 (m-1), z = height of the full level of the layer being crossed above the surface (m).
@inline tb_shallow_test_entrainment(z, shp::TiedtkeBechtoldShallowParameters) =
    shp.test_entrainment_b / z + shp.test_entrainment_a

"""
    tb_surface_parcel_departure(T, q, geopotential, excess, parameters) -> (s, q, K)

D49, D50: the shallow test parcel leaves full level n with s_bar_n + c_p Delta T, q_bar_n + Delta q, l = 0 and
K_dep = 0.5 w*^2. No U4 layer average (Cy41r2 states it only for deep departures).
"""
@inline function tb_surface_parcel_departure(T, q, geopotential, excess, par::TiedtkeBechtoldParameters)
    n = length(T)
    cp = par.heat_capacity
    @inbounds s = cp * T[n] + geopotential[n] + cp * excess.delta_T
    @inbounds q_departure = q[n] + excess.delta_q
    K = excess.wstar^2 / 2
    return s, q_departure, K
end

# D58: A_1B = {k <= n-1 : p_s - p_k <= 35000 Pa}, scanned k = n-1, n-2, ..., stopping at the first k outside A_1B.
# Level 1 is never a departure (as D07). An empty set gives last = first + 1.
@inline function tb_departure_bounds_stage1b(p_full, p_half, par::TiedtkeBechtoldParameters)
    n = length(p_full)
    surface_pressure = p_half[n + 1]
    first = n - 1
    last = first + 1
    k = first
    @inbounds while k >= 2 && surface_pressure - p_full[k] <= par.departure_range
        last = k
        k -= 1
    end
    return first, last
end

# D59: the n-1 departure reuses the D46 values of level n (no floor); every higher departure keeps (6.21) (D08).
@inline tb_stage1b_departure_perturbation(k, n, excess, par::TiedtkeBechtoldParameters) =
    k == n - 1 ? (excess.delta_T, excess.delta_q) : (par.delta_T_deep, par.delta_q_deep)

# D55: shallow iff the test cloud is at most 20000 Pa deep, the exact complement of the D14 deep test.
@inline tb_is_shallow_cloud(p_base, p_top, par::TiedtkeBechtoldParameters) = !tb_is_deep_cloud(p_base, p_top, par)

# D49-D55: the shallow test-parcel ascent from full level n. D51 U2 entrainment with z_k of the layer being crossed
# (z_n for the first half step from p_n to p_{n-1/2}), D10 exponential step, D03, D52 = D11 removal (the stage 1a
# `tb_test_parcel_step`), D54 = D12 K (tb_ke_v1 coefficients, loaded buoyancy, mu/M := eps_k) from K_dep = 0.5 w*^2,
# D13 base (saturated before adjustment, K > 0), top = highest half level with K > 0 (the D14 measure).
@inline function _tb_shallow_test_parcel_ascent(excess, T, q, p_full, p_half, geopotential, T_half, q_half,
                                                geopotential_half, shp::TiedtkeBechtoldShallowParameters,
                                                par::TiedtkeBechtoldParameters)
    NF = eltype(T)
    n = length(T)
    g = par.gravity
    cp = par.heat_capacity
    ke_buoyancy_factor = g / (par.f_turbulence * (1 + par.gamma))
    drag_factor = 1 + par.beta * par.C_d
    s, qp, K = tb_surface_parcel_departure(T, q, geopotential, excess, par)
    s_departure, q_departure, K_departure = s, qp, K
    l = zero(NF)
    @inbounds begin
        geopotential_surface = geopotential_half[n + 1]            # D37 host hydrostatic surface value
        T_parcel = (s - geopotential[n]) / cp
        b_lower = tb_buoyancy(T_parcel, qp, l, zero(NF), T[n], q[n], par) / tb_virtual_temperature(T[n], q[n], par)
    end
    has_base = false
    j_base = 0
    j_top = 0
    s_base = zero(NF)
    q_base = zero(NF)
    l_base = zero(NF)
    K_base = zero(NF)
    @inbounds for k in n:-1:2
        j_up = k                                                    # half level k-1/2
        dz = k == n ? (geopotential_half[k] - geopotential[k]) / g :
                      (geopotential_half[k] - geopotential_half[k + 1]) / g
        s_env = cp * T[k] + geopotential[k]
        eps_ini = tb_shallow_test_entrainment((geopotential[k] - geopotential_surface) / g, shp)   # D51
        step = tb_test_parcel_step(s, qp, l, s_env, q[k], eps_ini, dz, geopotential_half[j_up], p_half[j_up], par)
        b_upper = tb_buoyancy(step.T, step.q, step.l, zero(NF), T_half[j_up], q_half[j_up], par) /
                  tb_virtual_temperature(T_half[j_up], q_half[j_up], par)
        K = tb_kinetic_energy_step(K, drag_factor * eps_ini, ke_buoyancy_factor * (b_upper + b_lower) / 2, dz)
        s, qp, l = step.s, step.q, step.l
        b_lower = b_upper
        K > 0 || break                                              # w vanishes: the parcel stops
        if !has_base
            if step.saturated                                       # D13 / D54
                has_base = true
                j_base = j_up
                j_top = j_up
                s_base, q_base, l_base, K_base = s, qp, l, K
            end
        else
            j_top = j_up                                            # D55: highest half level with K > 0
        end
    end
    return (; has_base, j_base, j_top, s_base, q_base, l_base, K_base, s_departure, q_departure, K_departure)
end

# ----------------------------------------------------------------------------------------------------------------
# Shallow updraught layer (D61-D66)
# ----------------------------------------------------------------------------------------------------------------

"""
    tb_shallow_updraught_layer(M_lower, s_lower, q_lower, l_lower, K_lower, B_lower, lfc_passed, rh, f_scale, dz,
                               T_env_lower, q_env_lower, geopotential_lower, T_env_upper, q_env_upper,
                               geopotential_upper, p_upper, land_flag, shallow, parameters)

One shallow-updraught layer k (unit mass flux), as `tb_updraught_layer` with three stage 1B changes:
- D62: E dp/g = eps^(1) f_eps (1.3 - RH) f_scale M_{k+1/2} dz with f_eps = 2, if B_{k+1/2} > 0, else 0.
- D63 (replaces D18): D1 dp/g = (1.6 - RH) E dp/g (zero where E is zero).
- D65 (Sect. 8 step 6, before D22, D03, D25-D27, D04 and D20): E and D1 only form M' = M_{k+1/2} + E - D1. If
  M' <= 0 the layer is the top layer at once (`mass_guard = top = true`): nothing is mixed, adjusted, converted or
  integrated, G = 0, and the caller detrains all of M_{k+1/2} with the k+1/2 plume properties.
Otherwise D22 mixing, D03, D25-D27 conversion, D04 B and D20 K with c = (1 + beta C_d) max(e, d) of the fractional
rates of D62-D63, then the D21 top or D19 tb_orgdet_v1 (D64, unchanged), as for deep.
"""
@inline function tb_shallow_updraught_layer(M_lower, s_lower, q_lower, l_lower, K_lower, B_lower, lfc_passed, rh,
                                            f_scale, dz, T_env_lower, q_env_lower, geopotential_lower,
                                            T_env_upper, q_env_upper, geopotential_upper, p_upper,
                                            land_flag::Bool, shp::TiedtkeBechtoldShallowParameters,
                                            par::TiedtkeBechtoldParameters)
    NF = typeof(M_lower)
    cp = par.heat_capacity
    z = zero(NF)
    # D62: positively buoyant convection only (B at k+1/2), f_eps = 2
    entrainment_rate = B_lower > 0 ?
        par.eps1 * shp.f_eps * (par.rh_offset_entrainment - rh) * f_scale : z
    # D63: d = (1.6 - RH) e as fractional rates, D1 = (1.6 - RH) E as layer exchanges
    detrainment_rate = (shp.rh_offset_detrainment - rh) * entrainment_rate
    E = entrainment_rate * M_lower * dz
    D1 = (shp.rh_offset_detrainment - rh) * E
    M_mixed = M_lower + E - D1                                  # M' (D23, D65)
    if !(M_mixed > 0)
        # D65: not cloud. No D22, D03, D25-D27, D04 or D20 step in this layer.
        return (; mass_guard = true, top = true, E, D1, D2 = z, M_mixed, M_upper = z, s = s_lower, q = q_lower,
                  l = l_lower, T = (s_lower - geopotential_lower) / cp, K = K_lower, B = B_lower, G = z,
                  condensation = z, entrainment_rate, detrainment_rate)
    end
    # D22, (6.53): M' phi' = M_{k+1/2} phi_{k+1/2} + E phi_bar_{k+1/2} - D1 phi_{k+1/2}, l_bar = 0
    s_env_lower = cp * T_env_lower + geopotential_lower
    s_mixed = ((M_lower - D1) * s_lower + E * s_env_lower) / M_mixed
    q_mixed = ((M_lower - D1) * q_lower + E * q_env_lower) / M_mixed
    l_mixed = (M_lower - D1) * l_lower / M_mixed
    # Dry adiabatic ascent conserves s, then D03 at p_{k-1/2}
    T_mixed = (s_mixed - geopotential_upper) / cp
    T_adjusted, q_adjusted, l_adjusted = tb_saturation_adjustment(T_mixed, q_mixed, l_mixed, p_upper, par)
    condensation = M_mixed * (q_mixed - q_adjusted)           # D23 diagnostic
    # D75: D25-D27 unchanged for shallow
    w = min(sqrt(2 * K_lower), par.w_cap)
    threshold = land_flag ? par.l_threshold_land : par.l_threshold_water
    if l_mixed > threshold && w > 0
        a = par.c0 / (par.w_factor * w) * (1 - exp(-(l_mixed / par.l_crit)^2))
    else
        a = z
    end
    x = a * dz
    if x > 0
        l_analytic = l_mixed * exp(-x) + (l_adjusted - l_mixed) * (-expm1(-x) / x)
        l_new = clamp(l_analytic, z, l_adjusted)
    else
        l_new = l_adjusted
    end
    G = M_mixed * (l_adjusted - l_new)
    s_upper = cp * T_adjusted + geopotential_upper
    # D66 = D20 tb_ke_v1 with the shallow fractional rates
    B_upper = tb_buoyancy(T_adjusted, q_adjusted, l_new, z, T_env_upper, q_env_upper, par)
    b = par.gravity / (par.f_turbulence * (1 + par.gamma)) *
        (B_upper / tb_virtual_temperature(T_env_upper, q_env_upper, par) +
         B_lower / tb_virtual_temperature(T_env_lower, q_env_lower, par)) / 2
    c = (1 + par.beta * par.C_d) * max(entrainment_rate, detrainment_rate)
    K_upper = tb_kinetic_energy_step(K_lower, c, b, dz)
    top = !(K_upper > 0)                                       # D21
    # D64 = D19 tb_orgdet_v1 unchanged
    if !top && lfc_passed && !(B_lower > 0) && K_upper < K_lower
        R = (par.rh_offset_detrainment - rh) * sqrt(K_lower / K_upper)
        M_upper = M_mixed / max(R, one(NF))
    else
        M_upper = M_mixed
    end
    D2 = M_mixed - M_upper
    return (; mass_guard = false, top, E, D1, D2, M_mixed, M_upper, s = s_upper, q = q_adjusted, l = l_new,
              T = T_adjusted, K = K_upper, B = B_upper, G, condensation, entrainment_rate, detrainment_rate)
end

# ----------------------------------------------------------------------------------------------------------------
# Column blocks shared by the stage 1B deep and shallow paths (repeats of the stage 1a column body, same order)
# ----------------------------------------------------------------------------------------------------------------

# Stage 1a: every output and work entry written before any read, then Sect. 8 step 1 (D37, D38).
@inline function _tb_initialize_column!(dT, dq, rain_mass_flux_half, work, T, q, p_full, p_half, geopotential,
                                        par::TiedtkeBechtoldParameters)
    NF = eltype(T)
    n = length(T)
    cp = par.heat_capacity
    Lv = par.latent_heat_condensation
    z = zero(NF)
    (; temperature_half, humidity_half, geopotential_half, dry_static_energy_half, mass_flux_half,
       total_water_flux_half, liquid_static_energy_flux_half, entrainment, turbulent_detrainment,
       organized_detrainment, generation, rain_evaporation, detrained_condensate) = work
    @inbounds for k in 1:n
        dT[k] = z
        dq[k] = z
        entrainment[k] = z
        turbulent_detrainment[k] = z
        organized_detrainment[k] = z
        generation[k] = z
        rain_evaporation[k] = z
        detrained_condensate[k] = z
    end
    @inbounds for j in 1:(n + 1)
        rain_mass_flux_half[j] = z
        mass_flux_half[j] = z
        total_water_flux_half[j] = z
        liquid_static_energy_flux_half[j] = z
    end
    @inbounds begin
        temperature_half[1] = T[1]
        humidity_half[1] = q[1]
        geopotential_half[1] = geopotential[1]
        for k in 1:n
            q_sat = tb_saturation_humidity(T[k], p_full[k], par)
            slope = tb_saturation_humidity_slope(T[k], q_sat, par)
            dTdp = (par.R_dry * tb_virtual_temperature(T[k], q[k], par) / p_full[k] + Lv * q_sat / p_full[k]) /
                   (cp + Lv * slope)
            dqdp = slope * dTdp - q_sat / p_full[k]
            dp_half = p_half[k + 1] - p_full[k]
            temperature_half[k + 1] = T[k] + dTdp * dp_half
            humidity_half[k + 1] = q[k] + dqdp * dp_half
        end
        for j in 2:n
            geopotential_half[j] = geopotential[j] +
                par.R_dry * _tb_host_virtual_temperature(T[j], q[j], par) * log(p_full[j] / p_half[j])
        end
        geopotential_half[n + 1] = geopotential[n] -
            par.R_dry * _tb_host_virtual_temperature(T[n], q[n], par) * log(p_half[n + 1] / p_full[n])
        for j in 2:(n + 1)
            dry_static_energy_half[j] = cp * temperature_half[j] + geopotential_half[j]
        end
        dry_static_energy_half[1] = dry_static_energy_half[2]
        for j in n:-1:2
            dry_static_energy_half[j] = max(dry_static_energy_half[j], dry_static_energy_half[j + 1])
        end
    end
    return nothing
end

# Stage 1a step 4 (D15, D17-D27, D32, D33, D35 integrands): the deep unit-flux updraught from j_base.
@inline function _tb_deep_updraught!(work, j_base, s_base, q_base, l_base, K_base, T, q, p_full, p_half, timestep,
                                     land_flag::Bool, par::TiedtkeBechtoldParameters)
    NF = eltype(T)
    g = par.gravity
    cp = par.heat_capacity
    Lv = par.latent_heat_condensation
    z = zero(NF)
    (; temperature_half, humidity_half, geopotential_half, dry_static_energy_half, mass_flux_half,
       total_water_flux_half, liquid_static_energy_flux_half, entrainment, turbulent_detrainment,
       organized_detrainment, generation, detrained_condensate) = work
    unit_base_mass_flux = par.mstar_base_factor * (p_half[j_base] - p_half[j_base - 1]) / (g * timestep)
    p_base = p_half[j_base]
    M = unit_base_mass_flux
    s_up = s_base
    q_up = q_base
    l_up = l_base
    K_up = K_base
    test_parcel_base_w = sqrt(2 * K_base)
    @inbounds begin
        T_up = (s_up - geopotential_half[j_base]) / cp
        B_lower = tb_buoyancy(T_up, q_up, l_up, z, temperature_half[j_base], humidity_half[j_base], par)
        beta_lower = tb_pcape_integrand(T_up, q_up, l_up, z, temperature_half[j_base], humidity_half[j_base], par)
        mass_flux_half[j_base] = M
        total_water_flux_half[j_base] = M * (q_up + l_up - humidity_half[j_base])
        liquid_static_energy_flux_half[j_base] = M * (s_up - Lv * l_up - dry_static_energy_half[j_base])
        q_sat_base = tb_saturation_humidity(temperature_half[j_base], p_base, par)
    end
    lfc_passed = B_lower > 0
    PCAPE = z
    Den = z
    w_dz_sum = z
    condensation_column = z
    k_top = 0
    model_top = false
    @inbounds for k in (j_base - 1):-1:1
        dp = p_half[k + 1] - p_half[k]
        if k == 1
            organized_detrainment[1] = M
            detrained_condensate[1] = g / dp * M * l_up
            k_top = 1
            model_top = true
            break
        end
        dz = (geopotential_half[k] - geopotential_half[k + 1]) / g
        rh = tb_relative_humidity(T[k], q[k], p_full[k], par)
        f_scale = _tb_ipow(tb_saturation_humidity(temperature_half[k + 1], p_half[k + 1], par) / q_sat_base,
                           par.f_scale_exponent)
        layer = tb_updraught_layer(M, s_up, q_up, l_up, K_up, B_lower, lfc_passed, rh, f_scale, dz,
                                   temperature_half[k + 1], humidity_half[k + 1], geopotential_half[k + 1],
                                   temperature_half[k], humidity_half[k], geopotential_half[k], p_half[k],
                                   land_flag, par)
        if layer.top
            organized_detrainment[k] = M
            detrained_condensate[k] = g / dp * M * l_up
            k_top = k
            break
        end
        entrainment[k] = layer.E
        turbulent_detrainment[k] = layer.D1
        organized_detrainment[k] = layer.D2
        generation[k] = g * layer.G / dp
        detrained_condensate[k] = g / dp * (layer.D1 * l_up + layer.D2 * layer.l)
        condensation_column += layer.condensation
        beta_upper = tb_pcape_integrand(layer.T, layer.q, layer.l, z, temperature_half[k], humidity_half[k], par)
        PCAPE += (beta_lower + beta_upper) / 2 * dp
        Tv_upper = tb_virtual_temperature(temperature_half[k], humidity_half[k], par)
        Tv_lower = tb_virtual_temperature(temperature_half[k + 1], humidity_half[k + 1], par)
        Den += g / tb_virtual_temperature(T[k], q[k], par) * (M + layer.M_upper) / 2 *
               (Tv_upper - Tv_lower + (geopotential_half[k] - geopotential_half[k + 1]) / cp)
        w_dz_sum += (sqrt(2 * K_up) + sqrt(2 * layer.K)) / 2 * dz
        M = layer.M_upper
        s_up, q_up, l_up, K_up = layer.s, layer.q, layer.l, layer.K
        B_lower = layer.B
        beta_lower = beta_upper
        lfc_passed = lfc_passed || layer.B > 0
        mass_flux_half[k] = M
        total_water_flux_half[k] = M * (q_up + l_up - humidity_half[k])
        liquid_static_energy_flux_half[k] = M * (s_up - Lv * l_up - dry_static_energy_half[k])
    end
    return (; unit_base_mass_flux, PCAPE, Den, w_dz_sum, condensation_column, k_top, model_top, test_parcel_base_w)
end

# Stage 1B Sect. 8 step 6 (D61-D66): the shallow unit-flux updraught from j_base, M*_base of D34 and the D15 state.
# No PCAPE, Den or tau (D68: D31-D35 are deep only).
@inline function _tb_shallow_updraught!(work, j_base, s_base, q_base, l_base, K_base, T, q, p_full, p_half,
                                        timestep, land_flag::Bool, shp::TiedtkeBechtoldShallowParameters,
                                        par::TiedtkeBechtoldParameters)
    NF = eltype(T)
    g = par.gravity
    cp = par.heat_capacity
    Lv = par.latent_heat_condensation
    z = zero(NF)
    (; temperature_half, humidity_half, geopotential_half, dry_static_energy_half, mass_flux_half,
       total_water_flux_half, liquid_static_energy_flux_half, entrainment, turbulent_detrainment,
       organized_detrainment, generation, detrained_condensate) = work
    # D61, D34: unit flux M*_base = 0.1 dp_base/(g timestep), dp_base of the first cloud layer
    unit_base_mass_flux = par.mstar_base_factor * (p_half[j_base] - p_half[j_base - 1]) / (g * timestep)
    M = unit_base_mass_flux
    s_up = s_base                                               # D15: test-parcel w, s, q and l at the base
    q_up = q_base
    l_up = l_base
    K_up = K_base
    test_parcel_base_w = sqrt(2 * K_base)
    @inbounds begin
        T_up = (s_up - geopotential_half[j_base]) / cp
        B_lower = tb_buoyancy(T_up, q_up, l_up, z, temperature_half[j_base], humidity_half[j_base], par)
        mass_flux_half[j_base] = M
        total_water_flux_half[j_base] = M * (q_up + l_up - humidity_half[j_base])
        liquid_static_energy_flux_half[j_base] = M * (s_up - Lv * l_up - dry_static_energy_half[j_base])
        q_sat_base = tb_saturation_humidity(temperature_half[j_base], p_half[j_base], par)   # D06 (shallow base)
    end
    lfc_passed = B_lower > 0
    condensation_column = z
    k_top = 0
    model_top = false
    mass_guard_top = false
    @inbounds for k in (j_base - 1):-1:1
        dp = p_half[k + 1] - p_half[k]
        if k == 1
            # D21 r14 model top, as for deep
            organized_detrainment[1] = M
            detrained_condensate[1] = g / dp * M * l_up
            k_top = 1
            model_top = true
            break
        end
        dz = (geopotential_half[k] - geopotential_half[k + 1]) / g
        rh = tb_relative_humidity(T[k], q[k], p_full[k], par)                       # D05
        f_scale = _tb_ipow(tb_saturation_humidity(temperature_half[k + 1], p_half[k + 1], par) / q_sat_base,
                           par.f_scale_exponent)                                     # D06
        layer = tb_shallow_updraught_layer(M, s_up, q_up, l_up, K_up, B_lower, lfc_passed, rh, f_scale, dz,
                                           temperature_half[k + 1], humidity_half[k + 1], geopotential_half[k + 1],
                                           temperature_half[k], humidity_half[k], geopotential_half[k], p_half[k],
                                           land_flag, shp, par)
        if layer.top
            # D65 (M' <= 0) or D21 (K <= 0): M_{k-1/2} = 0, G_k = 0, no E or D1 is applied, and layer k detrains all
            # of M_{k+1/2} with the k+1/2 plume properties (counted as organized, as the D21 top).
            organized_detrainment[k] = M
            detrained_condensate[k] = g / dp * M * l_up
            k_top = k
            mass_guard_top = layer.mass_guard
            break
        end
        entrainment[k] = layer.E
        turbulent_detrainment[k] = layer.D1
        organized_detrainment[k] = layer.D2
        generation[k] = g * layer.G / dp
        detrained_condensate[k] = g / dp * (layer.D1 * l_up + layer.D2 * layer.l)
        condensation_column += layer.condensation
        M = layer.M_upper
        s_up, q_up, l_up, K_up = layer.s, layer.q, layer.l, layer.K
        B_lower = layer.B
        lfc_passed = lfc_passed || layer.B > 0
        mass_flux_half[k] = M
        total_water_flux_half[k] = M * (q_up + l_up - humidity_half[k])
        liquid_static_energy_flux_half[k] = M * (s_up - Lv * l_up - dry_static_energy_half[k])
    end
    return (; unit_base_mass_flux, condensation_column, k_top, model_top, mass_guard_top, test_parcel_base_w)
end

# D40 r13 (deep) and D67 (shallow, unchanged): F(j) = Z_j F(j_base), M(j) = Z_j M*_base, zero at the surface.
@inline function _tb_subcloud_fluxes!(work, j_base, unit_base_mass_flux, p_half)
    (; mass_flux_half, total_water_flux_half, liquid_static_energy_flux_half) = work
    n = length(p_half) - 1
    z = zero(unit_base_mass_flux)
    @inbounds begin
        surface_pressure = p_half[n + 1]
        p_base = p_half[j_base]
        total_water_base_flux = total_water_flux_half[j_base]
        liquid_static_energy_base_flux = liquid_static_energy_flux_half[j_base]
        for j in (j_base + 1):n
            Z = (surface_pressure - p_half[j]) / (surface_pressure - p_base)
            mass_flux_half[j] = Z * unit_base_mass_flux
            total_water_flux_half[j] = Z * total_water_base_flux
            liquid_static_energy_flux_half[j] = Z * liquid_static_energy_base_flux
        end
        mass_flux_half[n + 1] = z
        total_water_flux_half[n + 1] = z
        liquid_static_energy_flux_half[n + 1] = z
    end
    return nothing
end

# Stage 1a step 5: D35 tau and the D31 closure factor (PCAPE_bl = 0, F5).
@inline function _tb_deep_closure(work, j_base, j_top, PCAPE, Den, w_dz_sum, par::TiedtkeBechtoldParameters)
    geopotential_half = work.geopotential_half
    z = zero(PCAPE)
    @inbounds H = (geopotential_half[j_top] - geopotential_half[j_base]) / par.gravity
    tau = z
    closure_factor = z
    if H > 0 && w_dz_sum > 0
        w_bar = w_dz_sum / H
        alpha_x = 1 + par.alpha_x_coefficient * par.grid_spacing / par.dx_ref
        tau = clamp(H / w_bar * alpha_x, par.tau_min, par.tau_max)
        if PCAPE - par.PCAPE_bl > 0 && Den > 0
            closure_factor = (PCAPE - par.PCAPE_bl) / (tau * Den)
        end
    end
    return tau, closure_factor
end

# Stage 1a step 6, D41 r15 (and D76): one column factor from f_CFL and f_pos over limiter_timestep_factor * dt.
@inline function _tb_limiter_factors(work, closure_factor, q, p_half, timestep, par::TiedtkeBechtoldParameters)
    (; mass_flux_half, total_water_flux_half, generation) = work
    NF = typeof(closure_factor)
    n = length(q)
    g = par.gravity
    z = zero(NF)
    limiter_step = par.limiter_timestep_factor * timestep
    f_CFL = one(NF)
    f_pos = one(NF)
    if closure_factor > 0
        @inbounds for j in 2:n
            M_closed = closure_factor * mass_flux_half[j]
            if M_closed > 0
                dp_j = min(p_half[j] - p_half[j - 1], p_half[j + 1] - p_half[j])
                f_CFL = min(f_CFL, par.cfl_factor * dp_j / (g * limiter_step * M_closed))
            end
        end
        @inbounds for k in 1:n
            dp = p_half[k + 1] - p_half[k]
            dq_closed = closure_factor *
                (g / dp * (total_water_flux_half[k + 1] - total_water_flux_half[k]) - generation[k])
            if dq_closed < 0
                f_pos = min(f_pos, max(q[k], z) / (-dq_closed * limiter_step))
            end
        end
    end
    return f_CFL, f_pos
end

# Stage 1a zero closure: exact zero mass, flux and exchange profiles (the half-level environment is kept).
@inline function _tb_zero_profiles!(work, n)
    (; mass_flux_half, total_water_flux_half, liquid_static_energy_flux_half, entrainment, turbulent_detrainment,
       organized_detrainment, generation, detrained_condensate) = work
    z = zero(eltype(mass_flux_half))
    @inbounds for j in 1:(n + 1)
        mass_flux_half[j] = z
        total_water_flux_half[j] = z
        liquid_static_energy_flux_half[j] = z
    end
    @inbounds for k in 1:n
        entrainment[k] = z
        turbulent_detrainment[k] = z
        organized_detrainment[k] = z
        generation[k] = z
        detrained_condensate[k] = z
    end
    return nothing
end

# Stage 1a steps 7-8: D36 scaling by the final factor, D28 rain flux top down, D29 sub-cloud evaporation (p_k > p_base),
# D39 tendencies. Returns the surface rain mass flux (kg m-2 s-1).
@inline function _tb_scale_and_assemble!(dT, dq, rain_mass_flux_half, work, final_factor, T, q, p_full, p_half,
                                         p_base, timestep, land_flag::Bool, par::TiedtkeBechtoldParameters)
    (; mass_flux_half, total_water_flux_half, liquid_static_energy_flux_half, entrainment, turbulent_detrainment,
       organized_detrainment, generation, rain_evaporation, detrained_condensate) = work
    NF = eltype(T)
    n = length(T)
    g = par.gravity
    cp = par.heat_capacity
    Lv = par.latent_heat_condensation
    z = zero(NF)
    surface_pressure = p_half[n + 1]
    @inbounds for j in 1:(n + 1)
        mass_flux_half[j] *= final_factor
        total_water_flux_half[j] *= final_factor
        liquid_static_energy_flux_half[j] *= final_factor
    end
    @inbounds for k in 1:n
        entrainment[k] *= final_factor
        turbulent_detrainment[k] *= final_factor
        organized_detrainment[k] *= final_factor
        generation[k] *= final_factor
        detrained_condensate[k] *= final_factor
    end
    P = z
    @inbounds rain_mass_flux_half[1] = z
    @inbounds for k in 1:n
        dp = p_half[k + 1] - p_half[k]
        e_k = z
        if p_full[k] > p_base
            q_sat_k = tb_saturation_humidity(T[k], p_full[k], par)
            evaporated = tb_rain_evaporation_mass(P, q[k], q_sat_k, dp, timestep, land_flag, p_full[k],
                                                  surface_pressure, par)
            e_k = evaporated * g / dp
            P -= evaporated
        else
            P += generation[k] * dp / g
        end
        rain_evaporation[k] = e_k
        rain_mass_flux_half[k + 1] = P
        dq[k] = g / dp * (total_water_flux_half[k + 1] - total_water_flux_half[k]) - generation[k] + e_k
        dT[k] = (g / dp * (liquid_static_energy_flux_half[k + 1] - liquid_static_energy_flux_half[k]) +
                 Lv * (generation[k] - e_k)) / cp
    end
    return P
end

# D83 outputs: the stage 1a diagnostics plus the stage 1B fields. Counters are per-call increments (0 or 1); the host
# sums them in stage 2. For a deep column every shallow field is zero (D57), except `shallow_candidate_replaced`.
@inline function _tb_stage1b_result(::Type{NF}, base::NamedTuple, trigger_path, lagged_flux_unset, nonupward,
                                    nonpositive_jump, nonpositive_source, candidate_replaced, test_base_half_level,
                                    test_top_half_level, test_depth, test_base_w, wstar, delta_T, delta_q, source,
                                    large_scale_source, jump, closure_mass_flux, mass_guard_top) where {NF}
    return merge(base, (;
        trigger_path = Int(trigger_path),                          # 1 stage1a_fallback (D81), 2 stage1b
        lagged_flux_unset = Int(lagged_flux_unset),                # D81: tb_lagged_flux_unset_steps increment
        shallow_nonupward_buoyancy_flux = Int(nonupward),          # D53 counter increment
        shallow_nonpositive_jump = Int(nonpositive_jump),          # D73 counter increment (Delta h_base <= 0)
        shallow_nonpositive_source = Int(nonpositive_source),      # D73 counter increment (S_sub <= 0)
        shallow_candidate_replaced = Bool(candidate_replaced),     # D57: a shallow candidate lost to a deep cloud
        shallow_test_base_half_level = Int(test_base_half_level),  # D54 surface-parcel cloud base j; 0 = none
        shallow_test_top_half_level = Int(test_top_half_level),    # D55 highest half level with K > 0; 0 = none
        shallow_test_depth = NF(test_depth),                       # D55 p_base - p_top of the surface parcel (Pa)
        shallow_test_base_w = NF(test_base_w),                     # sqrt(2 K) of the surface parcel at its base
        surface_parcel_wstar = NF(wstar),                          # D47 w* (m s-1); 0 when unset
        surface_parcel_delta_T = NF(delta_T),                      # D46 Delta T (K), also the n-1 value (D59)
        surface_parcel_delta_q = NF(delta_q),                      # D46 Delta q (kg/kg), also the n-1 value (D59)
        shallow_source = NF(source),                               # S_sub (W m-2, upward positive), D68
        shallow_large_scale_source = NF(large_scale_source),       # S_LS (W m-2), D72; 0 for the HOST closure
        shallow_moist_static_energy_jump = NF(jump),               # Delta h_base (J kg-1), D69
        shallow_closure_mass_flux = NF(closure_mass_flux),         # M_base = S_sub/Delta h_base before D41 (kg m-2 s-1)
        shallow_mass_guard_top = Bool(mass_guard_top),             # D65 set the shallow cloud top
    ))
end

"""
    tiedtke_bechtold_column!(dT, dq, rain_mass_flux_half, work, T, q, p_full, p_half, geopotential, dt,
                             land_flag, latitude, parameters, shallow, lagged, closure) -> NamedTuple

Stage 1B column (STAGE1B_CHOICES.md r4, Sect. 8). The inputs and outputs of the stage 1a method, plus
`shallow::TiedtkeBechtoldShallowParameters`, `lagged::TiedtkeBechtoldLaggedSurfaceFlux` (the D45-D48 adapter values)
and `closure` (`TiedtkeBechtoldHostClosure()` for the host, D71, or `TiedtkeBechtoldSCMClosure` for the SCM
harness only, D72).

- `lagged.set == false` (D81): the stage 1a method runs unchanged (deep from n-2, no shallow test, no n-1 departure),
  and `lagged_flux_unset = 1`, `trigger_path = 1` (stage1a_fallback).
- Otherwise (`trigger_path = 2`): the shallow test from level n (D49-D56), the deep scan from n-1 (D58-D60), the
  first deep cloud replaces the shallow candidate (D57). Deep columns then follow stage 1a steps 4-8. Shallow columns
  follow D61-D76: unit-flux plume with the D65 guard, D67 sub-cloud scaling, M_base = S_sub/Delta h_base with the D73
  guard, D41 limiter, rain D75. Neither: exact zeros (D16, D56, D77).

`convection_type` is 0 none, 1 deep, 2 shallow. For shallow the column only reports the type and `tb_top_full_level`;
it never touches a host `cloud_top` (D79 (a), host wiring in stage 2). Same contract as stage 1a: no allocation, no
exceptions, isbits result (for the HOST closure), every work entry written before it is read.
"""
function tiedtke_bechtold_column!(
    dT::AbstractVector{NF}, dq::AbstractVector{NF}, rain_mass_flux_half::AbstractVector{NF}, work::NamedTuple,
    T::AbstractVector{NF}, q::AbstractVector{NF}, p_full::AbstractVector{NF}, p_half::AbstractVector{NF},
    geopotential::AbstractVector{NF}, dt::Real, land_flag::Bool, latitude::Real,
    par::TiedtkeBechtoldParameters{NF}, shp::TiedtkeBechtoldShallowParameters{NF},
    lagged::TiedtkeBechtoldLaggedSurfaceFlux{NF}, closure::AbstractTiedtkeBechtoldShallowClosure,
) where {NF <: AbstractFloat}
    z = zero(NF)
    n = length(T)

    # Sect. 8 step 2, D81: an unset lagged state runs the stage 1a column itself (exact stage 1a path).
    if !lagged.set
        stage1a = tiedtke_bechtold_column!(dT, dq, rain_mass_flux_half, work, T, q, p_full, p_half, geopotential,
                                           dt, land_flag, latitude, par)
        return _tb_stage1b_result(NF, stage1a, TIEDTKE_BECHTOLD_TRIGGER_PATH_STAGE1A_FALLBACK, true, false, false,
                                  false, false, 0, 0, z, z, z, z, z, z, z, z, z, false)
    end

    g = par.gravity
    Lv = par.latent_heat_condensation
    timestep = convert(NF, dt)
    # Step 1 (D37, D38), as stage 1a.
    _tb_initialize_column!(dT, dq, rain_mass_flux_half, work, T, q, p_full, p_half, geopotential, par)
    (; temperature_half, humidity_half, geopotential_half, dry_static_energy_half) = work

    # Step 3. Shallow test from level n (D49-D56). D47, D51: z_n is the height of level n above the surface.
    @inbounds z_n = (geopotential[n] - geopotential_half[n + 1]) / g
    @inbounds excess = tb_surface_parcel_excess(lagged, T[n], q[n], p_full[n], z_n, shp, par)
    shallow_test = _tb_shallow_test_parcel_ascent(excess, T, q, p_full, p_half, geopotential, temperature_half,
                                                  humidity_half, geopotential_half, shp, par)
    shallow_test_depth = z
    shallow_candidate = false
    if shallow_test.has_base
        @inbounds p_test_base = p_half[shallow_test.j_base]
        @inbounds p_test_top = p_half[shallow_test.j_top]
        shallow_test_depth = p_test_base - p_test_top
        shallow_candidate = tb_is_shallow_cloud(p_test_base, p_test_top, par)          # D55, D56
    end
    shallow_test_base_w = shallow_test.has_base ? sqrt(2 * shallow_test.K_base) : z

    # Step 4. Deep scan k = n-1, n-2, ... (D58-D60), the first deep cloud wins.
    first, last = tb_departure_bounds_stage1b(p_full, p_half, par)
    scanned = 0
    found = false
    k_dep = 0
    j_base = 0
    mixed_departure = false
    s_base = z
    q_base = z
    l_base = z
    K_base = z
    max_depth = z
    deep_depth_found = z
    k = first
    @inbounds while k >= last && !found
        scanned += 1
        delta_T, delta_q = tb_stage1b_departure_perturbation(k, n, excess, par)         # D59, D60
        ascent = _tb_test_parcel_ascent(k, delta_T, delta_q, T, q, p_full, p_half, geopotential, temperature_half,
                                        humidity_half, geopotential_half, par)
        if ascent.has_base
            depth = p_half[ascent.j_base] - p_half[ascent.j_top]
            max_depth = max(max_depth, depth)
            if tb_is_deep_cloud(p_half[ascent.j_base], p_half[ascent.j_top], par)
                found = true
                k_dep = k
                j_base = ascent.j_base
                mixed_departure = ascent.mixed
                s_base, q_base, l_base, K_base = ascent.s_base, ascent.q_base, ascent.l_base, ascent.K_base
                deep_depth_found = depth
            end
        end
        k -= 1
    end
    scan_last = scanned > 0 ? first - scanned + 1 : first + 1
    stage1b = TIEDTKE_BECHTOLD_TRIGGER_PATH_STAGE1B
    nonupward = excess.nonupward

    # Step 5, D57. Deep: the shallow candidate is replaced completely; stage 1a steps 4-8 follow. D65 for deep layers
    # (D18): the stage 1a layer already makes a layer with M' <= 0 the top layer with the D65 outcome (M_{k-1/2} = 0,
    # G_k = 0, all of M_{k+1/2} detrained with the k+1/2 properties), and it can bind only if 1.6 delta^(1) dz >= 1.
    if found
        deep = _tb_deep_updraught!(work, j_base, s_base, q_base, l_base, K_base, T, q, p_full, p_half, timestep,
                                   land_flag, par)
        j_top = deep.k_top + 1
        _tb_subcloud_fluxes!(work, j_base, deep.unit_base_mass_flux, p_half)
        tau, closure_factor = _tb_deep_closure(work, j_base, j_top, deep.PCAPE, deep.Den, deep.w_dz_sum, par)
        f_CFL, f_pos = _tb_limiter_factors(work, closure_factor, q, p_half, timestep, par)
        final_factor = closure_factor * min(f_CFL, f_pos)
        if !(final_factor > 0)
            _tb_zero_profiles!(work, n)
            base = _tb_result(NF, z, n + 1, j_base, 0, closure_factor, f_CFL, f_pos, z, deep.PCAPE, tau, z,
                              deep.unit_base_mass_flux, k_dep, scanned, first, scan_last, mixed_departure,
                              deep.test_parcel_base_w, K_base, deep_depth_found, max_depth, deep.model_top, z)
        else
            @inbounds p_base = p_half[j_base]
            P = _tb_scale_and_assemble!(dT, dq, rain_mass_flux_half, work, final_factor, T, q, p_full, p_half,
                                        p_base, timestep, land_flag, par)
            base = _tb_result(NF, P, deep.k_top, j_base, 1, closure_factor, f_CFL, f_pos, final_factor, deep.PCAPE,
                              tau, final_factor * deep.unit_base_mass_flux, deep.unit_base_mass_flux, k_dep, scanned,
                              first, scan_last, mixed_departure, deep.test_parcel_base_w, K_base, deep_depth_found,
                              max_depth, deep.model_top, final_factor * deep.condensation_column)
        end
        return _tb_stage1b_result(NF, base, stage1b, false, nonupward, false, false, shallow_candidate, 0, 0, z, z,
                                  excess.wstar, excess.delta_T, excess.delta_q, z, z, z, z, false)
    end

    # No deep and no shallow cloud (D16, D56, D77): exact zeros, already written. The surface-parcel test
    # diagnostics are kept, so a test cloud deeper than 200 hPa (D56) can be seen.
    if !shallow_candidate
        base = _tb_result(NF, z, n + 1, 0, 0, z, one(NF), one(NF), z, z, z, z, z, 0, scanned, first, scan_last,
                          false, z, z, z, max_depth, false, z)
        return _tb_stage1b_result(NF, base, stage1b, false, nonupward, false, false, false, shallow_test.j_base,
                                  shallow_test.j_top, shallow_test_depth, shallow_test_base_w, excess.wstar,
                                  excess.delta_T, excess.delta_q, z, z, z, z, false)
    end

    # Step 6. Shallow full updraught (D61-D66), then D67 sub-cloud scaling.
    j_base = shallow_test.j_base
    plume = _tb_shallow_updraught!(work, j_base, shallow_test.s_base, shallow_test.q_base, shallow_test.l_base,
                                   shallow_test.K_base, T, q, p_full, p_half, timestep, land_flag, shp, par)
    _tb_subcloud_fluxes!(work, j_base, plume.unit_base_mass_flux, p_half)

    # Step 7. Shallow closure. D69: Delta h_base = (c_p T_up + L_v q_up + Phi) - (s_bar + L_v q_bar) at j_base, with
    # c_p T_up + Phi = s_up and the D39 half-level s_bar (after D37, D38) and q_bar, so that F_sl + L_v F_qt at the base
    # is M Delta h_base. D71 / D72: S_sub. D45: L_v applied once. D73 guard. D68: M_base = S_sub / Delta h_base.
    @inbounds jump = (shallow_test.s_base - dry_static_energy_half[j_base]) +
                     Lv * (shallow_test.q_base - humidity_half[j_base])
    surface_source = lagged.sensible_heat_flux + Lv * lagged.humidity_flux
    source, large_scale_source = _tb_shallow_closure_source(closure, surface_source, j_base, p_half, par)
    nonpositive_jump = !(jump > 0)
    nonpositive_source = !(source > 0)
    if nonpositive_jump || nonpositive_source
        # D73: M_base = 0, non-convective column (exact zeros, D16), type none, matching counter(s) incremented.
        _tb_zero_profiles!(work, n)
        base = _tb_result(NF, z, n + 1, j_base, 0, z, one(NF), one(NF), z, z, z, z, plume.unit_base_mass_flux, n,
                          scanned, first, scan_last, false, plume.test_parcel_base_w, shallow_test.K_base, z,
                          max_depth, plume.model_top, z)
        return _tb_stage1b_result(NF, base, stage1b, false, nonupward, nonpositive_jump, nonpositive_source, false,
                                  shallow_test.j_base, shallow_test.j_top, shallow_test_depth, shallow_test_base_w,
                                  excess.wstar, excess.delta_T, excess.delta_q, source, large_scale_source, jump, z,
                                  plume.mass_guard_top)
    end
    closure_mass_flux = source / jump
    closure_factor = closure_mass_flux / plume.unit_base_mass_flux                       # M_base / M*_base (D61)

    # Step 8. D41 / D76 limiter: final factor = (M_base/M*_base) f (D74).
    f_CFL, f_pos = _tb_limiter_factors(work, closure_factor, q, p_half, timestep, par)
    final_factor = closure_factor * min(f_CFL, f_pos)
    if !(final_factor > 0)
        _tb_zero_profiles!(work, n)
        base = _tb_result(NF, z, n + 1, j_base, 0, closure_factor, f_CFL, f_pos, z, z, z, z,
                          plume.unit_base_mass_flux, n, scanned, first, scan_last, false, plume.test_parcel_base_w,
                          shallow_test.K_base, z, max_depth, plume.model_top, z)
        return _tb_stage1b_result(NF, base, stage1b, false, nonupward, false, false, false, shallow_test.j_base,
                                  shallow_test.j_top, shallow_test_depth, shallow_test_base_w, excess.wstar,
                                  excess.delta_T, excess.delta_q, source, large_scale_source, jump,
                                  closure_mass_flux, plume.mass_guard_top)
    end

    # Steps 9-10. D36 scaling, D28 rain, D29 sub-cloud evaporation, D39 tendencies (D75: unchanged for shallow).
    @inbounds p_base = p_half[j_base]
    P = _tb_scale_and_assemble!(dT, dq, rain_mass_flux_half, work, final_factor, T, q, p_full, p_half, p_base,
                                timestep, land_flag, par)
    base = _tb_result(NF, P, plume.k_top, j_base, 2, closure_factor, f_CFL, f_pos, final_factor, z, z,
                      final_factor * plume.unit_base_mass_flux, plume.unit_base_mass_flux, n, scanned, first,
                      scan_last, false, plume.test_parcel_base_w, shallow_test.K_base, z, max_depth, plume.model_top,
                      final_factor * plume.condensation_column)
    return _tb_stage1b_result(NF, base, stage1b, false, nonupward, false, false, false, shallow_test.j_base,
                              shallow_test.j_top, shallow_test_depth, shallow_test_base_w, excess.wstar,
                              excess.delta_T, excess.delta_q, source, large_scale_source, jump, closure_mass_flux,
                              plume.mass_guard_top)
end
