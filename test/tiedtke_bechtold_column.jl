# Tests of the stage 1a Tiedtke-Bechtold pure column (src/tiedtke_bechtold_column.jl) against the frozen contract
# FROZEN_CONTRACT_v1.sha256 b2953c07 (STAGE1_CHOICES_DRAFT.md r12, Sect. 10). Pure-column tests only: T1, T2, T3,
# T5, T6, T6b, T7, T8, T9, T11 and T13. T4 (default-off and restart identity), T10 (stage 1b), T12 and T14 (host
# wrapper) belong to later parts. Every test runs in Float64 and Float32, and every value is built in that format.
module TiedtkeBechtoldColumnTests
using ReadyESM, Test
const R = ReadyESM
const SW = R.SpeedyWeather

# Declared tolerances. Column identities (T2, T11, T13 totals): |residual| <= 64 n eps(NF) times the gross exchange
# (Float64 3.8e-13, Float32 2.1e-4 at n = 27). Per-layer and per-level identities: 16 eps(NF) relative.
closure_tolerance(::Type{NF}, n) where {NF} = 64 * n * eps(NF)
local_tolerance(::Type{NF}) where {NF} = 16 * eps(NF)

# --------------------------------------------------------------------------------------------------------------
# Synthetic profiles, built analytically in NF. sigma = p/p_s.
# - sigma >= sigma_m: well-mixed layer, dry adiabat T = T_s sigma^kappa, constant q = q_m;
# - above it: constant lapse rate lapse1 (K/m) in height, T = T_m (sigma/sigma_m)^(R_d lapse1/g);
# - above an optional inversion at sigma_i (0 = none): a jump of `inversion` K, then lapse rate lapse2;
# - T floored at T_trop (isothermal stratosphere); q = RH(sigma) q_sat(T, p) with RH piecewise linear in sigma
#   between `nodes` (host SpeedyWeather q_sat). Geopotential from the host hydrostatic discretisation.
# --------------------------------------------------------------------------------------------------------------

# Tropical ocean: 300 K, 6.5 K/km to a 195 K tropopause near 100 hPa, RH 85% below 900 hPa falling to 60% at
# 600 hPa and 30% at 200 hPa. Conditionally unstable (large CAPE). It triggers deep convection from k = n-2: the
# 0.2 K / 0.1 g/kg parcel at 925 hPa (RH 85%) saturates within one layer while still positively buoyant (base
# half level 23.5, test-parcel base speed 1.71 m/s), and the moist ascent stays warmer than the 6.5 K/km
# environment up to about 300 hPa, so the test cloud is about 710 hPa deep (> 200 hPa). Checked with an offline
# line-by-line mirror in Float64 and Float32 (decision margins: K >= 0.8 m2 s-2, saturation 7.5%, B >= 0.26 K).
const TROPICAL = (T_s = 300.0, sigma_m = 1.0, q_m = 0.0, lapse1 = 6.5e-3, sigma_i = 0.0, inversion = 0.0,
    lapse2 = 6.5e-3, T_trop = 195.0, nodes = ((0.9, 0.85), (0.6, 0.6), (0.2, 0.3)))
# Mid-latitude summer land: 298 K well-mixed boundary layer to 850 hPa (q 10.5 g/kg, RH 96% at its top), then
# 6.5 K/km to 215 K, RH 70% falling to 50% at 500 hPa. Deep (test cloud about 630 hPa, base speed 3.5 m/s).
const MIDLATITUDE = (T_s = 298.0, sigma_m = 0.85, q_m = 0.0105, lapse1 = 6.5e-3, sigma_i = 0.0, inversion = 0.0,
    lapse2 = 6.5e-3, T_trop = 215.0, nodes = ((0.85, 0.7), (0.5, 0.5), (0.2, 0.3)))
# High plateau (p_s = 600 hPa, 4 km): 285 K, 6.5 K/km, RH 90%. Level n-2 is 55.6 hPa above the surface, so the
# departure uses the 60 hPa mixed-layer branch (D09). Deep (test cloud about 330 hPa).
const PLATEAU = (T_s = 285.0, sigma_m = 1.0, q_m = 0.0, lapse1 = 6.5e-3, sigma_i = 0.0, inversion = 0.0,
    lapse2 = 6.5e-3, T_trop = 200.0, nodes = ((0.9, 0.9), (0.5, 0.6), (0.2, 0.3)))
# Absolutely stable: isothermal 260 K, RH 50%. No parcel reaches saturation.
const STABLE = (T_s = 260.0, sigma_m = 1.0, q_m = 0.0, lapse1 = 0.0, sigma_i = 0.0, inversion = 0.0,
    lapse2 = 0.0, T_trop = 100.0, nodes = ((1.0, 0.5),))
# Trade-wind family: 300 K well-mixed layer to sigma 0.89 (12.8 g/kg), a conditionally unstable cloud layer
# (5.5 K/km, RH 90%) capped by a 6 K inversion at sigma_i with dry air above (RH 25% then 20%). The test cloud
# runs from the mixed-layer top to the inversion: 7.5, 11.3 and 18.8 kPa deep for sigma_i = 0.75, 0.70, 0.65
# (no deep cloud) and 26.3 kPa for sigma_i = 0.60 (deep).
trade(sigma_i) = (T_s = 300.0, sigma_m = 0.89, q_m = 0.0128, lapse1 = 5.5e-3, sigma_i = sigma_i, inversion = 6.0,
    lapse2 = 5.0e-3, T_trop = 200.0,
    nodes = ((0.89, 0.9), (sigma_i + 0.01, 0.9), (sigma_i - 0.01, 0.25), (0.2, 0.2)))

function relative_humidity_profile(sigma::NF, nodes) where {NF}
    s0, r0 = NF(nodes[1][1]), NF(nodes[1][2])
    sigma >= s0 && return r0
    for i in 2:length(nodes)
        s1, r1 = NF(nodes[i][1]), NF(nodes[i][2])
        sigma >= s1 && return r0 + (r1 - r0) * (sigma - s0) / (s1 - s0)
        s0, r0 = s1, r1
    end
    return r0
end

function analytic_profile(::Type{NF}, p_full, p_s, spec, atmosphere, gravity) where {NF}
    kappa = atmosphere.R_dry / atmosphere.heat_capacity
    exponent1 = atmosphere.R_dry * NF(spec.lapse1) / gravity
    exponent2 = atmosphere.R_dry * NF(spec.lapse2) / gravity
    T_s, sigma_m, sigma_i, T_trop = NF(spec.T_s), NF(spec.sigma_m), NF(spec.sigma_i), NF(spec.T_trop)
    T_m = T_s * sigma_m^kappa
    T_i = sigma_i > 0 ? max(T_m * (sigma_i / sigma_m)^exponent1, T_trop) + NF(spec.inversion) : zero(NF)
    T = zeros(NF, length(p_full))
    q = zeros(NF, length(p_full))
    for k in eachindex(p_full)
        sigma = p_full[k] / NF(p_s)
        if sigma >= sigma_m
            T[k] = T_s * sigma^kappa
            q[k] = NF(spec.q_m)
        else
            T[k] = sigma >= sigma_i ? max(T_m * (sigma / sigma_m)^exponent1, T_trop) :
                                      max(T_i * (sigma / sigma_i)^exponent2, T_trop)
            q[k] = relative_humidity_profile(sigma, spec.nodes) *
                   SW.saturation_humidity(T[k], p_full[k], atmosphere)
        end
    end
    return T, q
end

# Host hydrostatics (SpeedyWeather 0.21.1 geopotential_compute!, written out with pressures).
function host_geopotential(T::Vector{NF}, q, p_full, p_half, surface_geopotential, atmosphere) where {NF}
    n = length(T)
    R_dry = atmosphere.R_dry
    mu = atmosphere.μ_virt_temp
    geopotential = zeros(NF, n)
    geopotential[n] = surface_geopotential + R_dry * T[n] * (1 + mu * q[n]) * log(p_half[n + 1] / p_full[n])
    for k in (n - 1):-1:1
        geopotential[k] = geopotential[k + 1] +
            R_dry * T[k + 1] * (1 + mu * q[k + 1]) * log(p_full[k + 1] / p_half[k + 1]) +
            R_dry * T[k] * (1 + mu * q[k]) * log(p_half[k + 1] / p_full[k])
    end
    return geopotential
end

# Grids: production L27 (equal sigma steps of 1/27, SpeedyWeather default_sigma_coordinates), a deliberately
# stretched 20-level grid (sigma_half = 1 - (1 - j/20)^2, fine near the surface) and a synthetic low-top column
# (ten 50 hPa layers from 500 to 1000 hPa) that the tropical plume overshoots.
equal_sigma_half(::Type{NF}, n) where {NF} = [NF(j) / NF(n) for j in 0:n]
stretched_sigma_half(::Type{NF}, n) where {NF} = [1 - (1 - NF(j) / NF(n))^2 for j in 0:n]
function sigma_grid(sigma_half::Vector{NF}, p_s) where {NF}
    p_half = sigma_half .* NF(p_s)
    p_full = [(sigma_half[k] + sigma_half[k + 1]) / 2 * NF(p_s) for k in 1:(length(sigma_half) - 1)]
    return p_full, p_half
end
function low_top_grid(::Type{NF}) where {NF}
    p_half = [NF(50000) + NF(5000) * NF(j) for j in 0:10]
    p_full = [(p_half[k] + p_half[k + 1]) / 2 for k in 1:10]
    return p_full, p_half
end
function make_grid(::Type{NF}, name::Symbol, p_s) where {NF}
    name === :L27 && return sigma_grid(equal_sigma_half(NF, 27), p_s)
    name === :stretched && return sigma_grid(stretched_sigma_half(NF, 20), p_s)
    return low_top_grid(NF)
end

const DEEP_CASES = (
    tropical = (grid = :L27, p_s = 101325.0, profile = TROPICAL, z_s = 0.0, land = false, dt = 900.0),
    midlatitude = (grid = :L27, p_s = 100000.0, profile = MIDLATITUDE, z_s = 0.0, land = true, dt = 900.0),
    plateau = (grid = :L27, p_s = 60000.0, profile = PLATEAU, z_s = 4000.0, land = true, dt = 900.0),
    stretched = (grid = :stretched, p_s = 100000.0, profile = TROPICAL, z_s = 0.0, land = false, dt = 900.0),
    low_top = (grid = :low_top, p_s = 100000.0, profile = TROPICAL, z_s = 0.0, land = false, dt = 900.0),
)

function build_column(::Type{NF}, case) where {NF}
    atmosphere = SW.EarthAtmosphere(NF)
    gravity = SW.Earth(NF).gravity
    p_full, p_half = make_grid(NF, case.grid, case.p_s)
    T, q = analytic_profile(NF, p_full, case.p_s, case.profile, atmosphere, gravity)
    geopotential = host_geopotential(T, q, p_full, p_half, gravity * NF(case.z_s), atmosphere)
    return (; T, q, p_full, p_half, geopotential, dt = NF(case.dt), land = case.land)
end

# Replace the humidity profile and rebuild the hydrostatic geopotential over the given surface geopotential.
function with_humidity(column, q, surface_geopotential)
    atmosphere = SW.EarthAtmosphere(eltype(q))
    geopotential = host_geopotential(column.T, q, column.p_full, column.p_half, surface_geopotential, atmosphere)
    return merge(column, (; q, geopotential))
end

function run_column(column, parameters; dt = column.dt, latitude = zero(eltype(column.T)),
                    prefill = zero(eltype(column.T)))
    NF = eltype(column.T)
    n = length(column.T)
    dT = fill(prefill, n)
    dq = fill(prefill, n)
    P = fill(prefill, n + 1)
    work = R.tiedtke_bechtold_column_work(NF, n)
    foreach(a -> fill!(a, prefill), values(work))
    res = R.tiedtke_bechtold_column!(dT, dq, P, work, column.T, column.q, column.p_full, column.p_half,
                                     column.geopotential, NF(dt), column.land, latitude, parameters)
    return (; dT, dq, P, work, res)
end

# T2: both named identities at the exact discrete dp[k] = p_half[k+1] - p_half[k] used by the column.
function ledgers(out, column, par)
    NF = eltype(out.dT)
    g, cp, Lv = par.gravity, par.heat_capacity, par.latent_heat_condensation
    rain = out.res.rain_mass_flux_surface
    water = rain                    # water_mass_closure_v1: sum dq dp/g + rain_mass_flux_surface
    water_gross = rain
    transform = zero(NF)            # latent_sensible_transform_v1: sum (c_p dT + L_v dq) dp/g
    transform_gross = zero(NF)
    for k in eachindex(out.dT)
        mass = (column.p_half[k + 1] - column.p_half[k]) / g
        water += out.dq[k] * mass
        water_gross += abs(out.dq[k]) * mass
        transform += (cp * out.dT[k] + Lv * out.dq[k]) * mass
        transform_gross += (cp * abs(out.dT[k]) + Lv * abs(out.dq[k])) * mass
    end
    return (; water, water_gross, transform, transform_gross)
end

function measure_allocations(dT, dq, P, work, T, q, p_full, p_half, geopotential, dt, land, latitude, par)
    R.tiedtke_bechtold_column!(dT, dq, P, work, T, q, p_full, p_half, geopotential, dt, land, latitude, par)
    return @allocated R.tiedtke_bechtold_column!(dT, dq, P, work, T, q, p_full, p_half, geopotential, dt, land,
                                                 latitude, par)
end

const PRECISIONS = (Float64, Float32)

@testset "Tiedtke-Bechtold stage 1a column" begin

@testset "Parameters, host constants and versions (Sect. 9, D02, D30) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    atmosphere = SW.EarthAtmosphere(NF)
    @test isbits(par)
    @test par.gravity == SW.Earth(NF).gravity == NF(9.81)
    @test par.heat_capacity == NF(1004)
    @test par.latent_heat_condensation == NF(2.501e6)
    @test isapprox(par.R_dry, NF(287.05); rtol = NF(1e-4))
    @test isapprox(par.R_vapor, NF(461.52); rtol = NF(1e-4))
    @test isapprox(par.mol_ratio, NF(0.62197); rtol = NF(1e-4))
    @test par.saturation_vapor_pressure == NF(610.78)
    @test par.temperature_freezing == NF(273.15)
    # Host constants come from the host objects, field for field: the parameters hold no literal copies.
    planet = SW.Earth(NF)
    @test (par.gravity, par.heat_capacity, par.latent_heat_condensation, par.R_dry, par.R_vapor, par.mol_ratio,
           par.mu_virt_temp, par.saturation_vapor_pressure, par.temperature_freezing) ==
          (planet.gravity, atmosphere.heat_capacity, atmosphere.latent_heat_condensation, atmosphere.R_dry,
           atmosphere.R_vapor, atmosphere.mol_ratio, atmosphere.μ_virt_temp, atmosphere.saturation_vapor_pressure,
           atmosphere.temperature_freezing)
    # The explicit constructor with the same constants gives the same parameters; a missing constant is an error;
    # a changed host constant propagates (no hidden default).
    host = R.tiedtke_bechtold_host_constants(planet, atmosphere)
    @test R.TiedtkeBechtoldParameters(NF; host...) === par
    @test_throws UndefKeywordError R.TiedtkeBechtoldParameters(NF; Base.structdiff(host, (; gravity = 0))...)
    changed = R.tiedtke_bechtold_parameters(NF; atmosphere = SW.EarthAtmosphere(NF; heat_capacity = 1005))
    @test changed.heat_capacity == NF(1005)
    # From a constructed model: its planet, atmosphere, number format and spectral truncation (T31).
    model = SW.PrimitiveWetModel(SW.SpectralGrid(; NF, trunc = 31, nlayers = 8))
    @test R.tiedtke_bechtold_parameters(model) === par
    for (name, value) in ((:eps1, 1.75e-3), (:f_eps, 1), (:rh_offset_entrainment, 1.3),
            (:rh_offset_detrainment, 1.6), (:delta1, 0.75e-4), (:gamma, 0.5), (:f_turbulence, 2), (:C_d, 0.506),
            (:beta, 1.875), (:test_entrainment_factor, 0.4), (:test_condensate_removal, 0.5),
            (:delta_T_deep, 0.2), (:delta_q_deep, 1e-4), (:w_departure, 1), (:mixed_layer_depth, 6000),
            (:average_depth, 5000), (:window_half_width, 2500), (:departure_range, 35000), (:deep_depth, 20000),
            (:mstar_base_factor, 0.1), (:alpha_x_coefficient, 1.66), (:dx_ref, 125e3), (:tau_min, 720),
            (:tau_max, 10800), (:PCAPE_bl, 0), (:c0, 1.4e-3), (:l_crit, 5e-4), (:w_factor, 0.75), (:w_cap, 10),
            (:l_threshold_water, 3e-4), (:l_threshold_land, 5e-4), (:alpha1, 5.44e-4), (:alpha2, 5.09e-3),
            (:alpha3, 0.5777), (:C_conv, 0.05), (:RH_cr_water, 0.9), (:RH_cr_land, 0.7),
            (:virtual_coefficient, 0.608), (:cfl_factor, 1), (:rh_min, 0), (:rh_max, 1), (:land_threshold, 0.5))
        @test getfield(par, name) == NF(value)
    end
    @test (par.f_scale_exponent, par.test_entrainment_exponent, par.subcloud_exponent, par.newton_iterations) ==
          (3, 3, 1, 3)
    @test par.rain_loading == false
    @test isapprox(par.grid_spacing, NF(645.6e3); rtol = NF(1e-3))              # D35: pi a / 31
    @test isapprox(1 + par.alpha_x_coefficient * par.grid_spacing / par.dx_ref, NF(9.57); rtol = NF(1e-3))
    # D02: the host saturation humidity, bit for bit.
    for T in NF.((220, 250, 273.15, 300, 310)), p in NF.((20000, 50000, 85000, 101325))
        @test R.tb_saturation_humidity(T, p, par) == SW.saturation_humidity(T, p, atmosphere)
    end
    # Fixed stage 1a decisions cannot be changed through the constructor.
    @test_throws ArgumentError R.tiedtke_bechtold_parameters(NF; rain_loading = true)
    @test_throws ArgumentError R.tiedtke_bechtold_parameters(NF; PCAPE_bl = 1)
    @test_throws ArgumentError R.tiedtke_bechtold_parameters(NF; cfl_factor = 3)
    @test_throws ArgumentError R.tiedtke_bechtold_parameters(NF; subcloud_exponent = 2)
    # D30 land flag, D05 RH clip.
    @test R.tb_land_flag(NF(0.5), par)
    @test !R.tb_land_flag(NF(0.49), par)
    @test R.tb_relative_humidity(NF(290), NF(0.05), NF(85000), par) == 1
    @test R.tb_relative_humidity(NF(290), zero(NF), NF(85000), par) == 0
    @test R.TIEDTKE_BECHTOLD_WATER_CONTRACT == "instantaneous_rain_no_rup_v1"
    @test R.TIEDTKE_BECHTOLD_KE_VERSION == "tb_ke_v1"
    @test R.TIEDTKE_BECHTOLD_ORGDET_VERSION == "tb_orgdet_v1"
    @test R.TIEDTKE_BECHTOLD_CONTRACT_MANIFEST_SHA256 ==
          "b2953c077772f6027eba2c1aa693b07509fa5e6ff922ed107260fe87f068fb95"
end

@testset "Column call: allocation-free, isbits, scratch-safe, latitude unused $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    for case in (DEEP_CASES.tropical, merge(DEEP_CASES.tropical, (; profile = STABLE)))
        column = build_column(NF, case)
        n = length(column.T)
        dT, dq, P = zeros(NF, n), zeros(NF, n), zeros(NF, n + 1)
        work = R.tiedtke_bechtold_column_work(NF, n)
        @test measure_allocations(dT, dq, P, work, column.T, column.q, column.p_full, column.p_half,
                                  column.geopotential, column.dt, column.land, zero(NF), par) == 0
        clean = run_column(column, par)
        @test isbitstype(typeof(clean.res))
        # Every work array is written before it is read: NaN-filled scratch gives bit-identical results.
        dirty = run_column(column, par; prefill = NF(NaN))
        @test isequal(clean.res, dirty.res)
        @test isequal((clean.dT, clean.dq, clean.P), (dirty.dT, dirty.dq, dirty.P))
        @test all(isequal(getfield(clean.work, f), getfield(dirty.work, f)) for f in keys(clean.work))
        # latitude is a diagnostics-only input
        shifted = run_column(column, par; latitude = NF(45))
        @test isequal(clean.res, shifted.res)
        @test isequal((clean.dT, clean.dq, clean.P), (shifted.dT, shifted.dq, shifted.P))
    end
end

@testset "T1 pressure gate (D07, D09, D14) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    # Measured L27 facts: candidate counts 7 / 9 / 14 (k = 19-25, 17-25, 12-25); level n-1 is 55.6 / 47.2 /
    # 33.3 hPa and level n-2 is 92.6 / 78.7 / 55.6 hPa above the surface at p_s = 1000 / 850 / 600 hPa.
    for (p_s, count, range, above_n1, above_n2) in ((100000, 7, 19:25, 5555.6, 9259.3),
                                                    (85000, 9, 17:25, 4722.2, 7870.4),
                                                    (60000, 14, 12:25, 3333.3, 5555.6))
        p_full, p_half = make_grid(NF, :L27, p_s)
        first, last = R.tb_departure_bounds(p_full, p_half, par)
        @test first - last + 1 == count
        @test last:first == range
        @test isapprox(NF(p_s) - p_full[26], NF(above_n1); rtol = NF(1e-4))
        @test isapprox(NF(p_s) - p_full[25], NF(above_n2); rtol = NF(1e-4))
    end
    # The scanned set equals A = {k <= n-2 : p_s - p_k <= 35000 Pa} on L27 and on the stretched grid, the scan
    # starts at n-2 (never n or n-1), the 350 hPa bound follows p_s, and the column scans exactly A when no
    # departure is deep (stable column).
    for grid in (:L27, :stretched), p_s in (100000.0, 85000.0, 60000.0)
        case = (grid = grid, p_s = p_s, profile = STABLE, z_s = (100000 - p_s) / 12, land = false, dt = 900.0)
        column = build_column(NF, case)
        n = length(column.T)
        A = [k for k in 1:(n - 2) if NF(p_s) - column.p_full[k] <= NF(35000)]
        first, last = R.tb_departure_bounds(column.p_full, column.p_half, par)
        @test collect(first:-1:last) == reverse(A)
        @test first == n - 2
        @test NF(p_s) - column.p_full[last] <= NF(35000)
        @test last == 2 || NF(p_s) - column.p_full[last - 1] > NF(35000)
        out = run_column(column, par)
        @test out.res.departures_scanned == length(A)
        @test out.res.departure_first == n - 2
        @test out.res.departure_last == minimum(A)
        @test out.res.convection_type == 0
        # The 60 hPa mixed-layer branch follows pressure, not a level count.
        for k in A
            mixed, lower, upper = R.tb_departure_window(column.p_full, column.p_half, k, par)
            @test mixed == (NF(p_s) - column.p_full[k] <= NF(6000))
            @test upper - lower == NF(5000)
            @test upper <= NF(p_s)
        end
    end
    # Window weights follow pressure. L27 at p_s = 600 hPa: k = 25 is 55.6 hPa above the surface (mixed branch),
    # window [p_25 - 2500, p_25 + 2500] overlaps layers 24, 25, 26 by 5/18, 4/9, 5/18 of 5000 Pa. An indicator
    # humidity in layer j returns its weight (T = 0 and geopotential = 0 make s_bar = 0).
    p_full, p_half = make_grid(NF, :L27, 60000)
    n = length(p_full)
    expected = Dict(24 => NF(5) / 18, 25 => NF(4) / 9, 26 => NF(5) / 18)
    weights = zeros(NF, n)
    for j in 1:n
        indicator = [NF(i == j) for i in 1:n]
        mixed, s_up, q_up = R.tb_departure_source(zeros(NF, n), indicator, zeros(NF, n), p_full, p_half, 25, par)
        @test mixed
        weights[j] = q_up - par.delta_q_deep
    end
    for j in 1:n
        @test isapprox(weights[j], get(expected, j, zero(NF)); atol = NF(1e-5))
    end
    @test isapprox(sum(weights), one(NF); rtol = NF(1e-5))
    mixed, s_up, q_up = R.tb_departure_source(fill(NF(290), n), fill(NF(0.01), n), zeros(NF, n), p_full,
                                              p_half, 24, par)
    @test !mixed                                                 # k = 24 is 77.8 hPa above the surface
    @test s_up == par.heat_capacity * NF(290) + par.heat_capacity * par.delta_T_deep
    @test q_up == NF(0.01) + par.delta_q_deep
    # Stretched grid, p_s = 1000 hPa: k = 18 is 16.25 hPa above the surface, so p_k + 2500 > p_s and the window
    # moves to [p_s - 5000, p_s]; k = 16-18 are inside 60 hPa, k = 15 (76.25 hPa) is not.
    p_full, p_half = make_grid(NF, :stretched, 100000)
    mixed, lower, upper = R.tb_departure_window(p_full, p_half, 18, par)
    @test mixed
    @test (lower, upper) == (NF(95000), NF(100000))
    @test [R.tb_departure_window(p_full, p_half, k, par)[1] for k in 15:18] == [false, true, true, true]
    total = zero(NF)
    for j in eachindex(p_full)
        total += R._tb_overlap(lower, upper, p_half[j], p_half[j + 1]) / par.average_depth
    end
    @test isapprox(total, one(NF); rtol = NF(1e-5))
    # Deep flag at 20000 Pa (D14): exactly 200 hPa is not deep, anything thicker is.
    @test !R.tb_is_deep_cloud(NF(90000), NF(70000), par)
    @test R.tb_is_deep_cloud(NF(90000), NF(69999), par)
    # ... and on a synthetic profile family: the trade inversion is raised until the test cloud is deeper than
    # 200 hPa. Depths 7.5, 11.3, 18.8 kPa give no convection, 26.3 kPa gives deep convection.
    for (sigma_i, deep) in ((0.75, false), (0.70, false), (0.65, false), (0.60, true))
        case = (grid = :L27, p_s = 101325.0, profile = trade(sigma_i), z_s = 0.0, land = false, dt = 900.0)
        out = run_column(build_column(NF, case), par)
        @test (out.res.convection_type == 1) == deep
        @test (out.res.max_test_cloud_depth > NF(20000)) == deep
        @test out.res.max_test_cloud_depth > 0                  # an LCL is found in every member
        @test out.res.departure_level <= 25                     # never n or n-1
    end
    # Deep columns launch from n-2 or above; the plateau departure uses the mixed-layer branch.
    for name in keys(DEEP_CASES)
        column = build_column(NF, DEEP_CASES[name])
        out = run_column(column, par)
        n = length(column.T)
        @test out.res.convection_type == 1
        @test 1 <= out.res.departure_level <= n - 2
        @test out.res.deep_test_depth > NF(20000)
        @test out.res.mixed_layer_departure == (name === :plateau || name === :stretched)
    end
end

@testset "T2 ledger: water_mass_closure_v1 and latent_sensible_transform_v1 (D39, D40) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    g, cp, Lv = par.gravity, par.heat_capacity, par.latent_heat_condensation
    for name in keys(DEEP_CASES)
        column = build_column(NF, DEEP_CASES[name])
        out = run_column(column, par)
        n = length(column.T)
        @test out.res.convection_type == 1
        @test out.res.rain_mass_flux_surface > 0
        l = ledgers(out, column, par)
        @test abs(l.water) <= closure_tolerance(NF, n) * l.water_gross
        @test abs(l.transform) <= closure_tolerance(NF, n) * l.transform_gross
        # The telescoping is visible: the conserved-variable fluxes vanish exactly at the model top (j = 1) and at
        # the surface (j = n+1), and each layer is the flux-form budget of its two half levels.
        F_qt = out.work.total_water_flux_half
        F_sl = out.work.liquid_static_energy_flux_half
        @test F_qt[1] == 0
        @test F_qt[n + 1] == 0
        @test F_sl[1] == 0
        @test F_sl[n + 1] == 0
        @test out.work.mass_flux_half[1] == 0
        @test out.work.mass_flux_half[n + 1] == 0
        G, e = out.work.generation, out.work.rain_evaporation
        for k in 1:n
            mass = (column.p_half[k + 1] - column.p_half[k]) / g
            source = (G[k] - e[k]) * mass
            water_layer = (F_qt[k + 1] - F_qt[k]) - source
            water_scale = abs(F_qt[k + 1] - F_qt[k]) + (G[k] + e[k]) * mass
            @test abs(out.dq[k] * mass - water_layer) <= local_tolerance(NF) * water_scale
            energy_layer = (F_sl[k + 1] - F_sl[k]) + Lv * source
            energy_scale = abs(F_sl[k + 1] - F_sl[k]) + Lv * (G[k] + e[k]) * mass
            @test abs(cp * out.dT[k] * mass - energy_layer) <= local_tolerance(NF) * energy_scale
        end
    end
end

@testset "T3 stable and no-deep columns give exact zeros (D16) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    cases = Any[(grid = grid, p_s = p_s, profile = STABLE, z_s = (100000 - p_s) / 12, land = false, dt = 900.0)
                for grid in (:L27, :stretched) for p_s in (100000.0, 85000.0, 60000.0)]
    for sigma_i in (0.75, 0.70, 0.65)        # shallow trade clouds: an LCL, but thinner than 200 hPa
        push!(cases, (grid = :L27, p_s = 101325.0, profile = trade(sigma_i), z_s = 0.0, land = false, dt = 900.0))
    end
    for case in cases, land in (false, true)
        column = merge(build_column(NF, case), (; land))
        out = run_column(column, par)
        n = length(column.T)
        @test out.res.convection_type == 0
        @test out.res.tb_top_full_level == n + 1
        @test out.res.cloud_base_half_level == 0
        @test out.res.rain_mass_flux_surface == 0
        @test all(iszero, out.dT)
        @test all(iszero, out.dq)
        @test all(iszero, out.P)
        @test all(iszero, out.work.mass_flux_half)
        @test all(iszero, out.work.generation)
        @test out.res.max_test_cloud_depth <= NF(20000)
    end
end

@testset "T5 limiters: forced CFL and positivity (D41) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    g = par.gravity
    cfl_numbers(out, column, dt) = [out.work.mass_flux_half[j] * g * NF(dt) /
        min(column.p_half[j] - column.p_half[j - 1], column.p_half[j + 1] - column.p_half[j])
        for j in 2:length(column.T)]
    column = build_column(NF, DEEP_CASES.tropical)
    n = length(column.T)
    # Not binding at the host step: strict CFL M dt g / dp < 1 everywhere.
    out = run_column(column, par)
    @test out.res.f_CFL == 1
    @test maximum(cfl_numbers(out, column, 900)) < 1
    # Forced CFL: a 6 h step (the closed mass flux does not depend on dt, the limit does).
    out = run_column(column, par; dt = 21600)
    numbers = cfl_numbers(out, column, 21600)
    @test out.res.f_CFL < 1
    @test out.res.f_pos == 1
    @test isapprox(out.res.final_factor, out.res.closure_factor * out.res.f_CFL; rtol = local_tolerance(NF))
    @test maximum(numbers) <= 1 + 64 * eps(NF)
    @test isapprox(maximum(numbers), one(NF); rtol = 64 * eps(NF))
    l = ledgers(out, column, par)
    @test abs(l.water) <= closure_tolerance(NF, n) * l.water_gross
    @test abs(l.transform) <= closure_tolerance(NF, n) * l.transform_gross
    # Forced positivity: the same tropical column with its two lowest layers made very dry (2 g/kg and 2e-4 g/kg,
    # a synthetic limiter case). The sub-cloud flux (D40 r13: the scaled cloud-base net flux) removes moisture from
    # layer n faster than it holds it.
    q = copy(column.q)
    q[n - 1] = NF(2e-3)
    q[n] = NF(2e-7)
    dry = with_humidity(column, q, zero(NF))                       # tropical surface geopotential is 0
    out = run_column(dry, par)
    @test out.res.convection_type == 1
    @test out.res.f_pos < 1
    @test out.res.f_pos < out.res.f_CFL
    @test isapprox(out.res.final_factor, out.res.closure_factor * out.res.f_pos; rtol = local_tolerance(NF))
    F_qt, G = out.work.total_water_flux_half, out.work.generation
    dt = dry.dt
    margins = zeros(NF, n)
    scales = zeros(NF, n)
    for k in 1:n
        dq_before_evaporation = g / (dry.p_half[k + 1] - dry.p_half[k]) * (F_qt[k + 1] - F_qt[k]) - G[k]
        margins[k] = q[k] + dq_before_evaporation * dt
        scales[k] = q[k] + abs(dq_before_evaporation) * dt
        @test margins[k] >= -64 * eps(NF) * scales[k]              # q + f dq dt >= 0 before evaporation
        @test q[k] + out.dq[k] * dt >= -64 * eps(NF) * scales[k]   # evaporation only adds vapour
    end
    @test minimum(abs.(margins) ./ scales) <= 64 * eps(NF)         # the binding level reaches zero
    l = ledgers(out, dry, par)
    @test abs(l.water) <= closure_tolerance(NF, n) * l.water_gross
    @test abs(l.transform) <= closure_tolerance(NF, n) * l.transform_gross
end

@testset "T6 tb_ke_v1 and the cloud top (D15, D20, D21) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    K, dz = NF(2), NF(500)
    # Drag sign: drag only decelerates (a printed "+" drag would accelerate an entraining plume).
    @test R.tb_kinetic_energy_step(K, zero(NF), zero(NF), dz) == K
    @test R.tb_kinetic_energy_step(K, NF(1e-3), zero(NF), dz) < K
    @test isapprox(R.tb_kinetic_energy_step(K, zero(NF), NF(0.01), dz), K + NF(0.01) * dz; rtol = 4 * eps(NF))
    @test R.tb_kinetic_energy_step(K, NF(1e-3), NF(0.01), dz) < K + NF(0.01) * dz
    # Branch c dz < 1 (c dz = 0.5): semi-implicit form.
    c, b = NF(1e-3), NF(0.01)
    @test isapprox(R.tb_kinetic_energy_step(K, c, b, dz), (K * (1 - c * dz) + b * dz) / (1 + c * dz);
                   rtol = 4 * eps(NF))
    # Branch c dz >= 1 (c dz = 2): layer-exact solution, which stays positive where the semi-implicit form
    # would overshoot to a spurious negative K.
    c = NF(4e-3)
    exact = K * exp(-2 * c * dz) + b / (2 * c) * (1 - exp(-2 * c * dz))
    @test isapprox(R.tb_kinetic_energy_step(K, c, b, dz), exact; rtol = 16 * eps(NF))
    @test (K * (1 - c * dz)) / (1 + c * dz) < 0
    @test R.tb_kinetic_energy_step(K, c, zero(NF), dz) > 0
    # Top handling at mid-column: the D21 top layer has M = 0 above it, no conversion, and detrains all of
    # M_{k+1/2}; nothing happens above it.
    for name in (:tropical, :midlatitude, :plateau, :stretched)
        column = build_column(NF, DEEP_CASES[name])
        out = run_column(column, par)
        n = length(column.T)
        k_top = out.res.tb_top_full_level
        @test out.res.convection_type == 1
        @test !out.res.model_top_reached
        @test 1 < k_top < out.res.cloud_base_half_level - 1
        if out.res.convection_type == 1 && 1 < k_top <= n
            M = out.work.mass_flux_half
            @test M[k_top] == 0
            @test M[k_top + 1] > 0
            @test out.work.generation[k_top] == 0
            @test out.work.entrainment[k_top] == 0
            @test out.work.turbulent_detrainment[k_top] == 0
            @test isapprox(out.work.organized_detrainment[k_top], M[k_top + 1]; rtol = local_tolerance(NF))
            @test all(iszero, out.dT[1:(k_top - 1)])
            @test all(iszero, out.dq[1:(k_top - 1)])
            @test all(iszero, M[1:k_top])
            @test all(iszero, out.P[1:(k_top + 1)])
            @test out.dq[k_top] > 0                                   # the detrained top layer is moistened
        end
    end
    # Top handling at the model top: the plume overshoots the 500 hPa lid of the low-top column, layer 1 becomes
    # the top layer and the flag is set.
    column = build_column(NF, DEEP_CASES.low_top)
    out = run_column(column, par)
    @test out.res.convection_type == 1
    @test out.res.model_top_reached
    @test out.res.tb_top_full_level == 1
    @test out.work.mass_flux_half[1] == 0
    @test out.work.mass_flux_half[2] > 0
    @test isapprox(out.work.organized_detrainment[1], out.work.mass_flux_half[2]; rtol = local_tolerance(NF))
    @test out.work.generation[1] == 0
    # D15: the full plume starts from the test parcel's cloud-base speed (1.71 m/s here), not 1 m/s.
    for name in keys(DEEP_CASES)
        out = run_column(build_column(NF, DEEP_CASES[name]), par)
        @test isapprox(out.res.updraught_base_kinetic_energy, out.res.test_parcel_base_w^2 / 2;
                       rtol = local_tolerance(NF))
    end
    out = run_column(build_column(NF, DEEP_CASES.tropical), par)
    @test abs(out.res.test_parcel_base_w - 1) > NF(0.1)
end

@testset "T6b test-parcel condensate removal (D11) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    cp = par.heat_capacity
    rtol = sqrt(eps(NF))
    # Layer 1: a parcel 5% supersaturated at 285 K and 850 hPa, no carried condensate, no environmental contrast.
    p1, geopotential1 = NF(85000), NF(14000)
    s0 = cp * NF(285) + geopotential1
    q0 = NF(1.05) * R.tb_saturation_humidity(NF(285), p1, par)
    step1 = R.tb_test_parcel_step(s0, q0, zero(NF), s0, q0, NF(1e-4), NF(300), geopotential1, p1, par)
    condensed1 = step1.q_mixed - step1.q
    @test step1.saturated
    @test condensed1 > 0
    @test isapprox(step1.l, (step1.l_mixed + condensed1) / 2; rtol)
    # Layer 2: 450 m higher; the carried l (nonzero) is mixed and new condensate forms. l_out = 0.5 (l_in + dl).
    p2, geopotential2 = NF(80000), geopotential1 + NF(9.81) * NF(450)
    step2 = R.tb_test_parcel_step(step1.s, step1.q, step1.l, step1.s, step1.q, NF(1e-4), NF(450), geopotential2,
                                  p2, par)
    condensed2 = step2.q_mixed - step2.q
    @test step1.l > 0
    @test step2.l_mixed > 0
    @test condensed2 > 0
    @test isapprox(step2.l, (step2.l_mixed + condensed2) / 2; rtol)
    @test step2.l == par.test_condensate_removal * step2.l_adjusted
end

@testset "T7 tb_orgdet_v1 (D19) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    tol = local_tolerance(NF)
    # Direct sweep of one layer over buoyancy sign, LFC state, K, RH, f_scale and the upper environment.
    cp = par.heat_capacity
    geopotential_lower = NF(2.6e4)
    geopotential_upper = geopotential_lower + NF(9.81) * NF(400)
    organized = 0
    for K in NF.((0.5, 5, 20)), B in NF.((-0.5, 0.5)), lfc in (false, true), rh in NF.((0.3, 0.9)),
        f_scale in NF.((1, 0.1)), dT_upper in NF.((-2, 0, 2))
        layer = R.tb_updraught_layer(NF(0.04), cp * NF(282) + geopotential_lower, NF(7.5e-3), NF(1e-3), K, B, lfc,
                                     rh, f_scale, NF(400), NF(280), NF(6e-3), geopotential_lower, NF(277) + dT_upper,
                                     NF(5e-3), geopotential_upper, NF(70000), false, par)
        @test (layer.E > 0) == (B > 0)                       # entrainment for positively buoyant air only
        layer.top && continue
        organized += layer.D2 > 0
        @test iszero(layer.E) || iszero(layer.D2)            # E and D2 never active in one layer
        @test lfc || iszero(layer.D2)                         # no D2 below the LFC
        @test layer.D2 >= 0
        @test layer.D2 <= layer.M_mixed
        @test layer.M_upper >= 0
    end
    @test organized > 0
    # The same properties on the deep columns, read from the column's own profiles.
    organized = 0
    for name in keys(DEEP_CASES)
        column = build_column(NF, DEEP_CASES[name])
        out = run_column(column, par)
        @test out.res.convection_type == 1
        w = out.work
        @test all(>=(0), w.mass_flux_half)
        k_top, j_base = out.res.tb_top_full_level, out.res.cloud_base_half_level
        out.res.convection_type == 1 || continue
        for k in (k_top + 1):(j_base - 1)
            E, D1, D2 = w.entrainment[k], w.turbulent_detrainment[k], w.organized_detrainment[k]
            M_lower, M_upper = w.mass_flux_half[k + 1], w.mass_flux_half[k]
            M_mixed = M_lower + E - D1
            scale = M_lower + E + D1 + D2
            organized += D2 > 0
            @test iszero(E) || iszero(D2)
            @test D2 >= 0
            @test D2 <= M_mixed + tol * scale
            @test abs(M_upper - (M_mixed - D2)) <= tol * scale
            # LFC passed: B > 0 at some half level at or below k+1/2 means E > 0 in some layer k' >= k.
            @test iszero(D2) || any(w.entrainment[kk] > 0 for kk in k:(j_base - 1))
        end
    end
    @test organized > 0
end

@testset "T8 f_scale reading and saturation adjustment (D03, D06) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    atmosphere = SW.EarthAtmosphere(NF)
    cp, Lv = par.heat_capacity, par.latent_heat_condensation
    # D06 reading (a): saturation specific humidity at both points (with pressure). At 500 hPa over a 950 hPa
    # base the vapour-pressure reading (b) is smaller by (950/500)^3 = 6.9, the "about x7" of D06.
    T, p, T_base, p_base = NF(265), NF(50000), NF(293), NF(95000)
    f_scale = R.tb_f_scale(T, p, T_base, p_base, par)
    @test isapprox(f_scale, (SW.saturation_humidity(T, p, atmosphere) /
                             SW.saturation_humidity(T_base, p_base, atmosphere))^3; rtol = local_tolerance(NF))
    vapour_reading = (R.tb_saturation_vapor_pressure(T, par) / R.tb_saturation_vapor_pressure(T_base, par))^3
    @test isapprox(f_scale / vapour_reading, (p_base / p)^3; rtol = sqrt(eps(NF)))
    @test f_scale / vapour_reading > 5
    # D03: supersaturated parcels (saturated at T0, then cooled by dTs) are adjusted to saturation within the
    # recorded residual; q + l and c_p T + L_v q are conserved to rounding.
    for (T0, p0) in ((300, 100000), (290, 85000), (250, 40000)), cooling in (2, 5, 8)
        T_cold = NF(T0) - NF(cooling)
        q0 = SW.saturation_humidity(NF(T0), NF(p0), atmosphere)
        T1, q1, l1 = R.tb_saturation_adjustment(T_cold, q0, zero(NF), NF(p0), par)
        q_sat = SW.saturation_humidity(T1, NF(p0), atmosphere)
        residual = abs(q1 - q_sat) / q_sat
        println("TB D03 residual ", NF, " T0=", T0, " p=", p0, " cooling=", cooling, " K: ", residual)
        tolerance = cooling <= 5 && NF === Float64 ? NF(1e-7) : NF(1e-5)
        @test l1 > 0
        @test residual <= tolerance
        @test isapprox(q1 + l1, q0; rtol = 8 * eps(NF))
        @test isapprox(cp * T1 + Lv * q1, cp * T_cold + Lv * q0; rtol = 8 * eps(NF))
    end
    # Two-way: subsaturated parcels evaporate carried l; enough l re-saturates, too little l evaporates fully.
    T0, p0 = NF(291), NF(85000)
    q0 = SW.saturation_humidity(NF(290), p0, atmosphere)
    T1, q1, l1 = R.tb_saturation_adjustment(T0, q0, NF(1e-3), p0, par)
    @test 0 < l1 < NF(1e-3)
    @test abs(q1 - SW.saturation_humidity(T1, p0, atmosphere)) <= NF(1e-5) * q1
    T1, q1, l1 = R.tb_saturation_adjustment(T0, q0, NF(1e-4), p0, par)
    @test l1 == 0
    @test q1 < SW.saturation_humidity(T1, p0, atmosphere)
    @test isapprox(q1 + l1, q0 + NF(1e-4); rtol = 8 * eps(NF))
    # Subsaturated without condensate: unchanged, exactly.
    @test R.tb_saturation_adjustment(T0, q0, zero(NF), p0, par) == (T0, q0, zero(NF))
end

@testset "T9 rain evaporation limits and units (D29) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    g = par.gravity
    # Units (P in kg m-2 s-1 as in (6.42)): P = 1e-3 and a deficit of 3e-3 evaporate about 18% of P per 100 hPa.
    rate = R.tb_rain_evaporation_rate(NF(1e-3), NF(3e-3), NF(1e5), NF(1e5), par)
    @test isapprox(rate * NF(10000) / g / NF(1e-3), NF(0.18338); rtol = NF(1e-3))
    @test isapprox(R.tb_rain_evaporation_rate(NF(2e-3), NF(3e-3), NF(1e5), NF(1e5), par) / rate,
                   NF(2)^par.alpha3; rtol = sqrt(eps(NF)))
    dp, dt = NF(10000), NF(900)
    # Limited by the available rain: e dp/g <= P.
    @test R.tb_rain_evaporation_mass(NF(1e-7), zero(NF), NF(0.02), dp, dt, false, NF(1e5), NF(1e5), par) ==
          NF(1e-7)
    # Limited by the saturation deficit: e dt <= RH_cr q_sat - q.
    deficit = par.RH_cr_water * NF(0.02) - NF(0.017)
    evaporated = R.tb_rain_evaporation_mass(one(NF), NF(0.017), NF(0.02), dp, dt, false, NF(1e5), NF(1e5), par)
    @test isapprox(evaporated, deficit * dp / (g * dt); rtol = 4 * eps(NF))
    @test R.tb_rain_evaporation_rate(one(NF), deficit, NF(1e5), NF(1e5), par) * dp / g > evaporated
    # Activation: only where q < RH_cr q_sat (0.9 over water, 0.7 over land); none without rain.
    @test R.tb_rain_evaporation_mass(NF(1e-4), NF(0.95) * NF(0.02), NF(0.02), dp, dt, false, NF(1e5), NF(1e5),
                                     par) == 0
    @test R.tb_rain_evaporation_mass(NF(1e-4), NF(0.8) * NF(0.02), NF(0.02), dp, dt, true, NF(1e5), NF(1e5),
                                     par) == 0
    @test R.tb_rain_evaporation_mass(NF(1e-4), NF(0.8) * NF(0.02), NF(0.02), dp, dt, false, NF(1e5), NF(1e5),
                                     par) > 0
    @test R.tb_rain_evaporation_mass(zero(NF), zero(NF), NF(0.02), dp, dt, false, NF(1e5), NF(1e5), par) == 0
    # In the columns: evaporation only below the cloud base, never more than the rain or the deficit allows.
    evaporating = 0
    for name in keys(DEEP_CASES)
        column = build_column(NF, DEEP_CASES[name])
        out = run_column(column, par)
        @test out.res.convection_type == 1
        out.res.convection_type == 1 || continue
        p_base = column.p_half[out.res.cloud_base_half_level]
        rh_cr = column.land ? par.RH_cr_land : par.RH_cr_water
        for k in eachindex(column.T)
            e = out.work.rain_evaporation[k]
            evaporating += e > 0
            @test e >= 0
            @test column.p_full[k] > p_base || e == 0
            mass = (column.p_half[k + 1] - column.p_half[k]) / g
            @test e * mass <= out.P[k] * (1 + local_tolerance(NF))
            deficit_k = rh_cr * R.tb_saturation_humidity(column.T[k], column.p_full[k], par) - column.q[k]
            @test e * column.dt <= max(deficit_k, zero(NF)) * (1 + local_tolerance(NF))
        end
    end
    @test evaporating > 0
end

@testset "T11 instantaneous_rain_no_rup_v1 (D28) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    g = par.gravity
    # A fixed column with generated rain closes the atmospheric water loss to the surface convective rain.
    column = build_column(NF, DEEP_CASES.tropical)
    out = run_column(column, par)
    n = length(column.T)
    loss = zero(NF)
    gross = zero(NF)
    for k in 1:n
        loss -= out.dq[k] * (column.p_half[k + 1] - column.p_half[k]) / g
        gross += abs(out.dq[k]) * (column.p_half[k + 1] - column.p_half[k]) / g
    end
    @test out.res.rain_mass_flux_surface > 0
    @test abs(loss - out.res.rain_mass_flux_surface) <= closure_tolerance(NF, n) * gross
    # r_up is zero by construction (no rain state, no loading) and the column returns only convective rain.
    @test R.TIEDTKE_BECHTOLD_WATER_CONTRACT == "instantaneous_rain_no_rup_v1"
    @test !par.rain_loading
    @test_throws ArgumentError R.tiedtke_bechtold_parameters(NF; rain_loading = true)
    @test !any(occursin("large", String(name)) || occursin("snow", String(name)) for name in keys(out.res))
    # Scalar counterfactual for the next increment: a positive r_up lowers the virtual-temperature buoyancy term
    # (6.22) by exactly r_up and the loaded buoyancy (6.56) by T_up r_up.
    T_up, q_up, l_up, T_env, q_env, r_up = NF(290), NF(0.012), NF(1e-3), NF(289), NF(0.011), NF(5e-3)
    term0 = R.tb_pcape_integrand(T_up, q_up, l_up, zero(NF), T_env, q_env, par)
    term1 = R.tb_pcape_integrand(T_up, q_up, l_up, r_up, T_env, q_env, par)
    @test isapprox(term0 - term1, r_up; rtol = sqrt(eps(NF)))
    B0 = R.tb_buoyancy(T_up, q_up, l_up, zero(NF), T_env, q_env, par)
    B1 = R.tb_buoyancy(T_up, q_up, l_up, r_up, T_env, q_env, par)
    @test isapprox(B0 - B1, T_up * r_up; rtol = sqrt(eps(NF)))
end

@testset "T13 local half-level rain flux (D28, D29, D43) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    g = par.gravity
    for name in keys(DEEP_CASES)
        column = build_column(NF, DEEP_CASES[name])
        out = run_column(column, par)
        n = length(column.T)
        P, G, e = out.P, out.work.generation, out.work.rain_evaporation
        @test out.res.convection_type == 1
        @test length(P) == n + 1
        @test P[1] == 0                                         # top interface
        @test P[n + 1] == out.res.rain_mass_flux_surface        # surface interface
        @test all(>=(0), P)
        j_base = out.res.cloud_base_half_level
        atol = local_tolerance(NF) * maximum(P)
        total = zero(NF)
        gross = zero(NF)
        for k in 1:n
            mass = (column.p_half[k + 1] - column.p_half[k]) / g
            @test abs((P[k + 1] - P[k]) - (G[k] - e[k]) * mass) <= atol
            @test iszero(G[k]) || iszero(e[k])
            @test k < j_base || G[k] == 0                       # no generation in sub-cloud layers
            @test k >= j_base || e[k] == 0                      # no evaporation in cloud layers
            total += (G[k] - e[k]) * mass
            gross += (G[k] + e[k]) * mass
        end
        @test abs(total - P[n + 1]) <= closure_tolerance(NF, n) * gross   # telescoped surface identity
    end
end

@testset "D37 half-level environment and D38 stability $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    atmosphere = SW.EarthAtmosphere(NF)
    column = build_column(NF, DEEP_CASES.tropical)
    out = run_column(column, par)
    n = length(column.T)
    cp, Lv, R_dry, R_v = par.heat_capacity, par.latent_heat_condensation, atmosphere.R_dry, atmosphere.R_vapor
    for k in (10, 20, n)
        T, q, p = column.T[k], column.q[k], column.p_full[k]
        q_sat = SW.saturation_humidity(T, p, atmosphere)
        slope = q_sat * Lv / (R_v * T^2)
        dTdp = (R_dry * T * (1 + NF(0.608) * q) / p + Lv * q_sat / p) / (cp + Lv * slope)
        dqdp = slope * dTdp - q_sat / p
        @test isapprox(out.work.temperature_half[k + 1], T + dTdp * (column.p_half[k + 1] - p);
                       rtol = local_tolerance(NF))
        @test isapprox(out.work.humidity_half[k + 1], q + dqdp * (column.p_half[k + 1] - p);
                       rtol = local_tolerance(NF))
    end
    @test all(out.work.dry_static_energy_half[j] >= out.work.dry_static_energy_half[j + 1] for j in 2:n)
end

@testset "T15 sub-cloud flux closure, D40 contract r13 $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    for name in keys(DEEP_CASES)
        column = build_column(NF, DEEP_CASES[name])
        out = run_column(column, par)
        n = length(column.T)
        @test out.res.convection_type == 1
        j_base = out.res.cloud_base_half_level
        @test 2 <= j_base <= n
        p_s = column.p_half[n + 1]
        p_base = column.p_half[j_base]
        F_qt = out.work.total_water_flux_half
        F_sl = out.work.liquid_static_energy_flux_half
        M = out.work.mass_flux_half
        # 1. Every sub-cloud interface carries the scaled cloud-base net flux, F(j) = Z(j) F(j_base), and the same
        # scaling of the mass flux.
        for j in (j_base + 1):n
            Z = (p_s - column.p_half[j]) / (p_s - p_base)
            @test abs(F_qt[j] - Z * F_qt[j_base]) <= local_tolerance(NF) * abs(F_qt[j_base])
            @test abs(F_sl[j] - Z * F_sl[j_base]) <= local_tolerance(NF) * abs(F_sl[j_base])
            @test abs(M[j] - Z * M[j_base]) <= local_tolerance(NF) * M[j_base]
        end
        # 2. Exactly zero at the surface (the telescoping column budgets are T2).
        @test F_qt[n + 1] == 0
        @test F_sl[n + 1] == 0
        @test M[n + 1] == 0
    end
    # 3. The moist-surface C08 regression: the tropical case (RH 85 % below 900 hPa) has an upward cloud-base
    # total-water flux, and the lowest layer dries. Under the r12 reading of D40 this layer was moistened.
    column = build_column(NF, DEEP_CASES.tropical)
    out = run_column(column, par)
    n = length(column.T)
    @test out.work.total_water_flux_half[out.res.cloud_base_half_level] > 0
    @test out.dq[n] < 0
end

end # top-level testset
end # module

println("TIEDTKE_BECHTOLD_COLUMN_TESTS_PASS")
