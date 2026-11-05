# Tests of the stage 1B Tiedtke-Bechtold column science (src/tiedtke_bechtold_column.jl, the stage 1B method of
# `tiedtke_bechtold_column!`) against the frozen contract frozen_v5/FROZEN_CONTRACT_v5.sha256 c6327c89
# (STAGE1B_CHOICES.md r4 b725f19f, Sect. 10). Column-level tests only: T1B, T16-T25, T28, T29 and the column-level part
# of T23. Stage 2 (PART B, host realisation) owns:
#   OPEN(D81, D82, T23): the two-step sentinel, the fresh-process pickup and the refusal of an absent, mismatched or
#     stage 1a checkpoint need the persistent lagged state, which is not column code.
#   OPEN(D79, T26, T27): the host `cloud_top` rule for shallow columns and the Betts-Miller / TB regression pair need the
#     stage 1B host wrapper. The present `tiedtke_bechtold_host_cloud_top` refuses type 2, so stage 2 must extend it.
#     Here the column only reports type 2 and `tb_top_full_level` (tested below).
#   OPEN(D72, T22): S_LS terms (2)-(3) as the harness forms them (Pi conversion, omega = -rho g wa subsidence derivatives,
#     remaps) and the interp_time step-start rule are SCM harness code (barry_scm_harness_v1). Here only the column's
#     sub-cloud integration of given tendency profiles is tested.
# Every test runs in Float64 and Float32, and every value is built in that format. The profile choices below were
# checked with an offline Float64 line-by-line mirror of the stage 1B column (decision margins quoted per case).
module TiedtkeBechtoldStage1BTests
using ReadyESM, Test
const R = ReadyESM
const SW = R.SpeedyWeather

# Declared tolerances, as in test/tiedtke_bechtold_column.jl. Column identities: 64 n eps(NF) times the gross
# exchange. Per-layer and per-level identities: 16 eps(NF) relative.
closure_tolerance(::Type{NF}, n) where {NF} = 64 * n * eps(NF)
local_tolerance(::Type{NF}) where {NF} = 16 * eps(NF)

# --------------------------------------------------------------------------------------------------------------
# Synthetic profiles, copied from test/tiedtke_bechtold_column.jl at 5fe01c3 (so this file does not depend on that
# test module). sigma = p/p_s. Well-mixed dry-adiabatic layer below sigma_m with constant q_m; constant lapse rate
# above; optional inversion at sigma_i; T floored at T_trop; q = RH(sigma) q_sat. Host hydrostatic geopotential.
# --------------------------------------------------------------------------------------------------------------
const TROPICAL = (T_s = 300.0, sigma_m = 1.0, q_m = 0.0, lapse1 = 6.5e-3, sigma_i = 0.0, inversion = 0.0,
    lapse2 = 6.5e-3, T_trop = 195.0, nodes = ((0.9, 0.85), (0.6, 0.6), (0.2, 0.3)))
const MIDLATITUDE = (T_s = 298.0, sigma_m = 0.85, q_m = 0.0105, lapse1 = 6.5e-3, sigma_i = 0.0, inversion = 0.0,
    lapse2 = 6.5e-3, T_trop = 215.0, nodes = ((0.85, 0.7), (0.5, 0.5), (0.2, 0.3)))
const PLATEAU = (T_s = 285.0, sigma_m = 1.0, q_m = 0.0, lapse1 = 6.5e-3, sigma_i = 0.0, inversion = 0.0,
    lapse2 = 6.5e-3, T_trop = 200.0, nodes = ((0.9, 0.9), (0.5, 0.6), (0.2, 0.3)))
const STABLE = (T_s = 260.0, sigma_m = 1.0, q_m = 0.0, lapse1 = 0.0, sigma_i = 0.0, inversion = 0.0,
    lapse2 = 0.0, T_trop = 100.0, nodes = ((1.0, 0.5),))
# Trade-wind family (the stage 1a T1 family): 300 K mixed layer to sigma 0.89 (12.8 g/kg, RH 62% at level n), a
# 5.5 K/km RH 90% cloud layer and a 6 K inversion at sigma_i with dry air above.
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

# Grids: production L27 (equal sigma 1/27), the stage 1a stretched 20-level grid, and an exact grid of 25 layers of
# 4000 Pa (p_s = 100000 Pa) whose half-level pressures are exact in Float32 and Float64, so a test-cloud depth of
# exactly 20000 Pa (five layers) is representable.
equal_sigma_half(::Type{NF}, n) where {NF} = [NF(j) / NF(n) for j in 0:n]
stretched_sigma_half(::Type{NF}, n) where {NF} = [1 - (1 - NF(j) / NF(n))^2 for j in 0:n]
function sigma_grid(sigma_half::Vector{NF}, p_s) where {NF}
    p_half = sigma_half .* NF(p_s)
    p_full = [(sigma_half[k] + sigma_half[k + 1]) / 2 * NF(p_s) for k in 1:(length(sigma_half) - 1)]
    return p_full, p_half
end
function exact_grid(::Type{NF}) where {NF}
    p_half = [NF(4000) * NF(j) for j in 0:25]
    p_full = [(p_half[k] + p_half[k + 1]) / 2 for k in 1:25]
    return p_full, p_half
end
function make_grid(::Type{NF}, name::Symbol, p_s) where {NF}
    name === :L27 && return sigma_grid(equal_sigma_half(NF, 27), p_s)
    name === :stretched && return sigma_grid(stretched_sigma_half(NF, 20), p_s)
    return exact_grid(NF)
end

case(profile; grid = :L27, p_s = 101325.0, z_s = 0.0, land = false, dt = 900.0) = (; grid, p_s, profile, z_s, land, dt)

function build_column(::Type{NF}, c) where {NF}
    atmosphere = SW.EarthAtmosphere(NF)
    gravity = SW.Earth(NF).gravity
    p_full, p_half = make_grid(NF, c.grid, c.p_s)
    T, q = analytic_profile(NF, p_full, c.p_s, c.profile, atmosphere, gravity)
    geopotential = host_geopotential(T, q, p_full, p_half, gravity * NF(c.z_s), atmosphere)
    return (; T, q, p_full, p_half, geopotential, dt = NF(c.dt), land = c.land, z_s = NF(c.z_s))
end

# Change one full level k (humidity set to `q`, temperature raised by `dT`) and rebuild the hydrostatic geopotential.
function with_level(column, k; q = nothing, dT = 0)
    NF = eltype(column.T)
    T = copy(column.T)
    humidity = copy(column.q)
    T[k] += NF(dT)
    q === nothing || (humidity[k] = NF(q))
    atmosphere = SW.EarthAtmosphere(NF)
    gravity = SW.Earth(NF).gravity
    geopotential = host_geopotential(T, humidity, column.p_full, column.p_half, gravity * column.z_s, atmosphere)
    return merge(column, (; T, q = humidity, geopotential))
end

# Lagged surface flux of a warm evaporating trade-wind ocean: H = 15 W m-2 and F_q = 6e-5 kg m-2 s-1 (latent
# 150 W m-2), both upward. Unset: the declared n = 0 state.
lagged(::Type{NF}, H = 15.0, F_q = 6.0e-5) where {NF} = R.TiedtkeBechtoldLaggedSurfaceFlux(NF(H), NF(F_q))
unset(::Type{NF}) where {NF} = R.TiedtkeBechtoldLaggedSurfaceFlux{NF}()
host_closure() = R.TiedtkeBechtoldHostClosure()

function run_stage1a(column, par; dt = column.dt, prefill = zero(eltype(column.T)))
    NF = eltype(column.T)
    n = length(column.T)
    dT, dq, P = fill(prefill, n), fill(prefill, n), fill(prefill, n + 1)
    work = R.tiedtke_bechtold_column_work(NF, n)
    foreach(a -> fill!(a, prefill), values(work))
    res = R.tiedtke_bechtold_column!(dT, dq, P, work, column.T, column.q, column.p_full, column.p_half,
                                     column.geopotential, NF(dt), column.land, zero(NF), par)
    return (; dT, dq, P, work, res)
end

function run_stage1b(column, par, flux; shallow = R.tiedtke_bechtold_shallow_parameters(eltype(column.T)),
                     closure = host_closure(), dt = column.dt, prefill = zero(eltype(column.T)))
    NF = eltype(column.T)
    n = length(column.T)
    dT, dq, P = fill(prefill, n), fill(prefill, n), fill(prefill, n + 1)
    work = R.tiedtke_bechtold_column_work(NF, n)
    foreach(a -> fill!(a, prefill), values(work))
    res = R.tiedtke_bechtold_column!(dT, dq, P, work, column.T, column.q, column.p_full, column.p_half,
                                     column.geopotential, NF(dt), column.land, zero(NF), par, shallow, flux, closure)
    return (; dT, dq, P, work, res)
end

function measure_allocations(dT, dq, P, work, T, q, p_full, p_half, geopotential, dt, land, latitude, par, shallow,
                             flux, closure)
    R.tiedtke_bechtold_column!(dT, dq, P, work, T, q, p_full, p_half, geopotential, dt, land, latitude, par, shallow,
                               flux, closure)
    return @allocated R.tiedtke_bechtold_column!(dT, dq, P, work, T, q, p_full, p_half, geopotential, dt, land,
                                                 latitude, par, shallow, flux, closure)
end

# T2 identities at the exact discrete dp[k] used by the column.
function ledgers(out, column, par)
    NF = eltype(out.dT)
    g, cp, Lv = par.gravity, par.heat_capacity, par.latent_heat_condensation
    rain = out.res.rain_mass_flux_surface
    water, water_gross = rain, rain
    transform, transform_gross = zero(NF), zero(NF)
    for k in eachindex(out.dT)
        mass = (column.p_half[k + 1] - column.p_half[k]) / g
        water += out.dq[k] * mass
        water_gross += abs(out.dq[k]) * mass
        transform += (cp * out.dT[k] + Lv * out.dq[k]) * mass
        transform_gross += (cp * abs(out.dT[k]) + Lv * abs(out.dq[k])) * mass
    end
    return (; water, water_gross, transform, transform_gross)
end

# The surface-parcel excess exactly as the column forms it (z_n from the column's own D37 surface geopotential).
function column_excess(column, work, flux, shallow, par)
    n = length(column.T)
    z_n = (column.geopotential[n] - work.geopotential_half[n + 1]) / par.gravity
    return R.tb_surface_parcel_excess(flux, column.T[n], column.q[n], column.p_full[n], z_n, shallow, par)
end

# The column's own shallow test ascent (D49-D55) on the half-level work arrays of a previous call.
column_shallow_ascent(column, work, excess, shallow, par) =
    R._tb_shallow_test_parcel_ascent(excess, column.T, column.q, column.p_full, column.p_half, column.geopotential,
                                     work.temperature_half, work.humidity_half, work.geopotential_half, shallow, par)

"""
Test-side re-implementation of the D49-D55 ascent from the public layer primitives, with the full-level heights
`heights[k]` (m above the surface) supplied by the caller. It returns the result and the per-step trace.
"""
function traced_shallow_ascent(column, work, excess, shallow, par, heights)
    NF = eltype(column.T)
    n = length(column.T)
    g, cp = par.gravity, par.heat_capacity
    T, q, phi, p_half = column.T, column.q, column.geopotential, column.p_half
    T_half, q_half, phi_half = work.temperature_half, work.humidity_half, work.geopotential_half
    ke_buoyancy_factor = g / (par.f_turbulence * (1 + par.gamma))
    drag_factor = 1 + par.beta * par.C_d
    s, qp, K = R.tb_surface_parcel_departure(T, q, phi, excess, par)
    l = zero(NF)
    b_lower = R.tb_buoyancy((s - phi[n]) / cp, qp, l, zero(NF), T[n], q[n], par) /
              R.tb_virtual_temperature(T[n], q[n], par)
    steps = []
    has_base, j_base, j_top = false, 0, 0
    s_base = q_base = l_base = K_base = zero(NF)
    for k in n:-1:2
        dz = k == n ? (phi_half[k] - phi[k]) / g : (phi_half[k] - phi_half[k + 1]) / g
        eps_k = R.tb_shallow_test_entrainment(heights[k], shallow)
        step = R.tb_test_parcel_step(s, qp, l, cp * T[k] + phi[k], q[k], eps_k, dz, phi_half[k], p_half[k], par)
        b_upper = R.tb_buoyancy(step.T, step.q, step.l, zero(NF), T_half[k], q_half[k], par) /
                  R.tb_virtual_temperature(T_half[k], q_half[k], par)
        K = R.tb_kinetic_energy_step(K, drag_factor * eps_k, ke_buoyancy_factor * (b_upper + b_lower) / 2, dz)
        push!(steps, (; j = k, eps = eps_k, z = heights[k], step, K))
        s, qp, l = step.s, step.q, step.l
        b_lower = b_upper
        K > 0 || break
        if !has_base
            if step.saturated
                has_base, j_base, j_top = true, k, k
                s_base, q_base, l_base, K_base = s, qp, l, K
            end
        else
            j_top = k
        end
    end
    return (; has_base, j_base, j_top, s_base, q_base, l_base, K_base, steps)
end

const STAGE1B_FIELDS = (:trigger_path, :lagged_flux_unset, :shallow_nonupward_buoyancy_flux,
    :shallow_nonpositive_jump, :shallow_nonpositive_source, :shallow_candidate_replaced,
    :shallow_test_base_half_level, :shallow_test_top_half_level, :shallow_test_depth, :shallow_test_base_w,
    :surface_parcel_wstar, :surface_parcel_delta_T, :surface_parcel_delta_q, :shallow_source,
    :shallow_large_scale_source, :shallow_moist_static_energy_jump, :shallow_closure_mass_flux,
    :shallow_mass_guard_top)
# Fields that describe the shallow candidate or closure (all zero in a deep or a fallback column).
const SHALLOW_ONLY_FIELDS = (:shallow_nonpositive_jump, :shallow_nonpositive_source, :shallow_test_base_half_level,
    :shallow_test_top_half_level, :shallow_test_depth, :shallow_test_base_w, :shallow_source,
    :shallow_large_scale_source, :shallow_moist_static_energy_jump, :shallow_closure_mass_flux,
    :shallow_mass_guard_top)

const PRECISIONS = (Float64, Float32)

@testset "Tiedtke-Bechtold stage 1B column" begin

@testset "Stage 1B parameters, inputs, closures and identity (Sect. 9, D45, D71, D72, D82) $NF" for NF in PRECISIONS
    shallow = R.tiedtke_bechtold_shallow_parameters(NF)
    @test isbits(shallow)
    @test R.TiedtkeBechtoldShallowParameters(NF) === shallow
    for (name, value) in ((:excess_factor, 1.5), (:delta_T_cap, 3), (:delta_q_cap, 2e-3), (:wstar_factor, 1.2),
            (:wstar_flux_factor, 1.5), (:wstar_virtual_coefficient, 0.61), (:von_karman, 0.4),
            (:friction_velocity, 0.1), (:test_entrainment_b, 0.8), (:test_entrainment_a, 2e-4), (:f_eps, 2),
            (:rh_offset_detrainment, 1.6))
        @test getfield(shallow, name) == NF(value)
    end
    @test_throws ArgumentError R.tiedtke_bechtold_shallow_parameters(NF; friction_velocity = 0)
    @test_throws ArgumentError R.tiedtke_bechtold_shallow_parameters(NF; test_entrainment_a = -1)
    # Values shared with stage 1a (D52 = D11 removal, D59 w_dep, D55 = D14 depth, D67 = D40 m).
    par = R.tiedtke_bechtold_parameters(NF)
    @test par.test_condensate_removal == NF(0.5)
    @test par.w_departure == 1
    @test par.deep_depth == NF(20000)
    @test par.subcloud_exponent == 1
    # Lagged adapter values: upward-positive H (W m-2) and F_q (kg m-2 s-1), set flag; the declared unset state.
    flux = lagged(NF)
    @test isbits(flux) && flux.set
    @test (flux.sensible_heat_flux, flux.humidity_flux) == (NF(15), NF(6e-5))
    @test !unset(NF).set && unset(NF).sensible_heat_flux == 0 && unset(NF).humidity_flux == 0
    # G6: the HOST closure has no field (no S_LS input), the SCM closure holds the two tendency profiles, and the host
    # wrapper source has no path to the SCM closure.
    @test fieldcount(R.TiedtkeBechtoldHostClosure) == 0
    @test isbits(R.TiedtkeBechtoldHostClosure())
    @test fieldnames(R.TiedtkeBechtoldSCMClosure) == (:temperature_tendency, :humidity_tendency)
    host_source = read(joinpath(R.PROJECT_ROOT, "src", "tiedtke_bechtold_host.jl"), String)
    @test !occursin("SCMClosure", host_source)
    @test !occursin("large_scale_source", host_source)
    F = NF(123.5)
    @test R._tb_shallow_closure_source(R.TiedtkeBechtoldHostClosure(), F, 24, zeros(NF, 28), par) === (F, zero(NF))
    # Identity strings (D82 names are proposals, fixed with Codex in stage 2).
    @test R.TIEDTKE_BECHTOLD_STAGE1B_CONTRACT_MANIFEST_SHA256 ==
          "c6327c89e97ee1b4676a0d737f46f07cafa87351622b37fcff7db927284eb898"
    @test R.TIEDTKE_BECHTOLD_STAGE1B_TRIGGER_VARIANT == "shallow_n_deep_n1_v1"
    @test (R.TIEDTKE_BECHTOLD_TRIGGER_PATH_STAGE1A_FALLBACK, R.TIEDTKE_BECHTOLD_TRIGGER_PATH_STAGE1B) == (1, 2)
    # The stage 1a identity is untouched.
    @test R.TIEDTKE_BECHTOLD_TRIGGER_VARIANT == "deep_from_n_minus_2"
    @test R.TIEDTKE_BECHTOLD_PARCEL_ADAPTER == "none"
end

@testset "Stage 1B call: allocation-free, isbits, scratch-safe, inferred, D83 outputs $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    shallow = R.tiedtke_bechtold_shallow_parameters(NF)
    for (profile, flux) in ((trade(0.70), lagged(NF)), (TROPICAL, lagged(NF)), (STABLE, lagged(NF)),
                            (trade(0.70), unset(NF)))
        column = build_column(NF, case(profile))
        n = length(column.T)
        dT, dq, P = zeros(NF, n), zeros(NF, n), zeros(NF, n + 1)
        work = R.tiedtke_bechtold_column_work(NF, n)
        @test measure_allocations(dT, dq, P, work, column.T, column.q, column.p_full, column.p_half,
                                  column.geopotential, column.dt, column.land, zero(NF), par, shallow, flux,
                                  host_closure()) == 0
        clean = run_stage1b(column, par, flux)
        @test isbitstype(typeof(clean.res))
        @test all(name -> haskey(clean.res, name), STAGE1B_FIELDS)
        @test all(name -> haskey(clean.res, name), keys(run_stage1a(column, par).res))
        dirty = run_stage1b(column, par, flux; prefill = NF(NaN))
        @test isequal(clean.res, dirty.res)
        @test isequal((clean.dT, clean.dq, clean.P), (dirty.dT, dirty.dq, dirty.P))
        @test all(isequal(getfield(clean.work, f), getfield(dirty.work, f)) for f in keys(clean.work))
        @test (@inferred R.tiedtke_bechtold_column!(dT, dq, P, work, column.T, column.q, column.p_full, column.p_half,
                   column.geopotential, column.dt, column.land, zero(NF), par, shallow, flux, host_closure())) isa
              NamedTuple
    end
    # D83 type flag in {0 none, 1 deep, 2 shallow} and the trigger path.
    @test run_stage1b(build_column(NF, case(trade(0.70))), par, lagged(NF)).res.convection_type == 2
    @test run_stage1b(build_column(NF, case(TROPICAL)), par, lagged(NF)).res.convection_type == 1
    @test run_stage1b(build_column(NF, case(STABLE)), par, lagged(NF)).res.convection_type == 0
    @test run_stage1b(build_column(NF, case(STABLE)), par, lagged(NF)).res.trigger_path == 2
    @test run_stage1b(build_column(NF, case(STABLE)), par, unset(NF)).res.trigger_path == 1
end

@testset "T1B pressure gate, stage 1B (D58-D60) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    shallow = R.tiedtke_bechtold_shallow_parameters(NF)
    for grid in (:L27, :stretched), p_s in (100000.0, 85000.0, 60000.0)
        column = build_column(NF, case(STABLE; grid, p_s, z_s = (100000 - p_s) / 12))
        n = length(column.T)
        A_1B = [k for k in 2:(n - 1) if NF(p_s) - column.p_full[k] <= NF(35000)]
        A = [k for k in 2:(n - 2) if NF(p_s) - column.p_full[k] <= NF(35000)]
        first, last = R.tb_departure_bounds_stage1b(column.p_full, column.p_half, par)
        @test first == n - 1                                            # n-1 is scanned first
        @test collect(first:-1:last) == reverse(A_1B)
        @test last == 2 || NF(p_s) - column.p_full[last - 1] > NF(35000)
        # Lagged state set: the stable column scans exactly A_1B and finds nothing.
        out = run_stage1b(column, par, lagged(NF))
        @test out.res.departure_first == n - 1
        @test out.res.departures_scanned == length(A_1B)
        @test out.res.departure_last == minimum(A_1B)
        @test out.res.convection_type == 0
        # Lagged state unset: exactly the stage 1a set A.
        fallback = run_stage1b(column, par, unset(NF))
        @test fallback.res.departure_first == n - 2
        @test fallback.res.departures_scanned == length(A)
        # The 60 hPa branch at n-1 follows pressure (at L27, n-1 is 55.6 / 47.2 / 33.3 hPa above the surface).
        @test R.tb_departure_window(column.p_full, column.p_half, n - 1, par)[1] ==
              (NF(p_s) - column.p_full[n - 1] <= NF(6000))
        grid === :L27 && @test R.tb_departure_window(column.p_full, column.p_half, n - 1, par)[1]
    end
    # D59: the n-1 perturbation is the D46 value of level n, with no floor (here a downward flux makes it negative);
    # every higher departure keeps (6.21). D60: the n-1 source is the D09 window mean plus that perturbation.
    column = build_column(NF, case(trade(0.70)))
    n = length(column.T)
    work = run_stage1b(column, par, lagged(NF)).work
    for flux in (lagged(NF), lagged(NF, -20.0, -1.0e-5))
        excess = column_excess(column, work, flux, shallow, par)
        @test R.tb_stage1b_departure_perturbation(n - 1, n, excess, par) === (excess.delta_T, excess.delta_q)
        @test R.tb_stage1b_departure_perturbation(n - 2, n, excess, par) === (par.delta_T_deep, par.delta_q_deep)
        mixed, s_n1, q_n1 = R.tb_departure_source(column.T, column.q, column.geopotential, column.p_full,
                                                  column.p_half, n - 1, excess.delta_T, excess.delta_q, par)
        mixed_a, s_a, q_a = R.tb_departure_source(column.T, column.q, column.geopotential, column.p_full,
                                                  column.p_half, n - 1, par)
        @test mixed && mixed_a
        cp = par.heat_capacity
        @test abs((s_n1 - s_a) - cp * (excess.delta_T - par.delta_T_deep)) <= 8 * eps(NF) * abs(s_a)
        @test abs((q_n1 - q_a) - (excess.delta_q - par.delta_q_deep)) <= 8 * eps(NF) * abs(q_a)
        if flux.sensible_heat_flux < 0
            @test excess.delta_T < 0 && excess.delta_q < 0                # no floor (D46, D59)
        end
    end
    # With the lagged state set, the deep tropical and mid-latitude clouds now leave from n-1 (w_dep = 1 m/s, the
    # mixed-layer branch at L27); the stage 1a fallback still leaves from n-2.
    for (profile, p_s) in ((TROPICAL, 101325.0), (MIDLATITUDE, 100000.0))
        column = build_column(NF, case(profile; p_s))
        n = length(column.T)
        out = run_stage1b(column, par, lagged(NF))
        @test out.res.convection_type == 1
        @test out.res.departure_level == n - 1
        @test out.res.mixed_layer_departure
        @test run_stage1b(column, par, unset(NF)).res.departure_level == n - 2
    end
end

@testset "T16 shallow test parcel (D49-D51) and D45-D48 sign cases $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    shallow = R.tiedtke_bechtold_shallow_parameters(NF)
    g, cp = par.gravity, par.heat_capacity
    rtol = sqrt(eps(NF))
    column = build_column(NF, case(trade(0.70)))
    n = length(column.T)
    flux = lagged(NF)
    out = run_stage1b(column, par, flux)
    @test out.res.convection_type == 2
    # Independent hydrostatic heights above the surface: the test's own geopotential minus g z_s.
    heights = [(column.geopotential[k] - g * column.z_s) / g for k in 1:n]
    @test isapprox(out.work.geopotential_half[n + 1], g * column.z_s; atol = 64 * eps(NF) * column.geopotential[1])
    # D46, D47 written out with the independent z_n, rho_bar and T_bar of level n (D45 signs, L_v not used).
    T_n, q_n, p_n = column.T[n], column.q[n], column.p_full[n]
    H, F_q = flux.sensible_heat_flux, flux.humidity_flux
    rho = p_n / (par.R_dry * T_n * (1 + par.mu_virt_temp * q_n))
    ustar3 = NF(0.1)^3
    argument = ustar3 + NF(1.5) * g * heights[n] * NF(0.4) / (rho * T_n) * (H / cp + NF(0.61) * T_n * F_q)
    wstar = NF(1.2) * cbrt(max(argument, ustar3))
    @test isapprox(out.res.surface_parcel_wstar, wstar; rtol)
    @test isapprox(out.res.surface_parcel_delta_T, min(NF(3), NF(1.5) * H / (rho * cp * wstar)); rtol)
    @test isapprox(out.res.surface_parcel_delta_q, min(NF(2e-3), NF(1.5) * F_q / (rho * wstar)); rtol)
    # D49, D50: departure state at full level n.
    excess = column_excess(column, out.work, flux, shallow, par)
    s_dep, q_dep, K_dep = R.tb_surface_parcel_departure(column.T, column.q, column.geopotential, excess, par)
    @test isapprox(s_dep, cp * T_n + column.geopotential[n] + cp * excess.delta_T; rtol = 4 * eps(NF))
    @test isapprox(q_dep, q_n + excess.delta_q; rtol = 4 * eps(NF))
    @test isapprox(K_dep, excess.wstar^2 / 2; rtol = 4 * eps(NF))
    # D51: eps = 0.8/z + 2e-4 m-1.
    for z in NF.((50, 165, 1000, 5000))
        @test isapprox(R.tb_shallow_test_entrainment(z, shallow), NF(0.8) / z + NF(2e-4); rtol = 4 * eps(NF))
    end
    # Every layer crossed uses eps_k of its own full-level height: the column's ascent equals the test-side ascent
    # driven by the independent heights (to rounding), and the first half step uses z_n.
    ascent = column_shallow_ascent(column, out.work, excess, shallow, par)
    traced = traced_shallow_ascent(column, out.work, excess, shallow, par, heights)
    @test ascent.has_base && traced.has_base
    @test (ascent.j_base, ascent.j_top) == (traced.j_base, traced.j_top) == (24, 21)
    @test (out.res.shallow_test_base_half_level, out.res.shallow_test_top_half_level) == (ascent.j_base, ascent.j_top)
    @test isapprox(ascent.s_base, traced.s_base; rtol)
    @test isapprox(ascent.q_base, traced.q_base; rtol)
    @test isapprox(ascent.K_base, traced.K_base; rtol)
    @test (ascent.s_departure, ascent.q_departure, ascent.K_departure) == (s_dep, q_dep, K_dep)
    @test traced.steps[1].j == n && traced.steps[1].z == heights[n]
    @test all(isapprox(st.eps, NF(0.8) / heights[st.j] + NF(2e-4); rtol) for st in traced.steps)
    @test length(traced.steps) >= n - ascent.j_top                         # every layer up to the top
    # D48 sign cases: a warm evaporating ocean gives a warmer, moister parcel and a larger w*; no flux gives no excess
    # and the u* floor; the sign-reversed case gives a colder, drier parcel (no floor) at the u* floor.
    excess_of(H, F_q) = column_excess(column, out.work, lagged(NF, H, F_q), shallow, par)
    warm, none, reversed, stronger = (excess_of(15.0, 6e-5), excess_of(0.0, 0.0), excess_of(-15.0, -6e-5),
                                      excess_of(30.0, 1.2e-4))
    floor_w = NF(1.2) * cbrt(NF(0.1)^3)
    @test warm.delta_T > 0 && warm.delta_q > 0 && warm.wstar > floor_w
    @test stronger.wstar > warm.wstar
    @test none.delta_T == 0 && none.delta_q == 0 && isapprox(none.wstar, floor_w; rtol = 4 * eps(NF))
    @test reversed.delta_T < 0 && reversed.delta_q < 0 && isapprox(reversed.wstar, floor_w; rtol = 4 * eps(NF))
    @test !warm.nonupward && none.nonupward && reversed.nonupward
    # Caps (D46): very large fluxes (w* about 11 m/s, uncapped excess about 12 K and 0.12 kg/kg) give 3 K and 2e-3.
    capped = excess_of(1.0e5, 1.0)
    @test capped.delta_T == NF(3) && capped.delta_q == NF(2e-3)
end

@testset "T17 classification (D54-D56) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    shallow = R.tiedtke_bechtold_shallow_parameters(NF)
    # D55 rule: shallow iff p_base - p_top <= 20000 Pa (exact complement of D14).
    @test R.tb_is_shallow_cloud(NF(90000), NF(70001), par)
    @test R.tb_is_shallow_cloud(NF(90000), NF(70000), par)
    @test !R.tb_is_shallow_cloud(NF(90000), NF(69999), par)
    for (p_base, p_top) in ((NF(90000), NF(70000)), (NF(85000), NF(60000)), (NF(95000), NF(94000)))
        @test R.tb_is_shallow_cloud(p_base, p_top, par) == !R.tb_is_deep_cloud(p_base, p_top, par)
    end
    # L27 trade family with the warm-ocean fluxes: the surface parcel saturates at half level 24 and stops below the
    # inversion: 7.5, 11.3 and 15.0 kPa deep (mirror), all shallow, and no deep cloud from n-1 (7.5, 11.3, 18.8 kPa).
    for (sigma_i, j_top) in ((0.75, 22), (0.70, 21), (0.65, 20))
        column = build_column(NF, case(trade(sigma_i)))
        out = run_stage1b(column, par, lagged(NF))
        @test out.res.convection_type == 2
        @test (out.res.shallow_test_base_half_level, out.res.shallow_test_top_half_level) == (24, j_top)
        @test out.res.shallow_test_depth == column.p_half[24] - column.p_half[j_top]
        @test 0 < out.res.shallow_test_depth <= NF(20000)
        @test out.res.cloud_base_half_level == 24
        @test out.res.max_test_cloud_depth <= NF(20000)
    end
    # Just below, at and just above 20000 Pa on the exact 4000 Pa grid. Sensitivity setting: U2 b = 0.1 (instead of
    # 0.8) lets the surface parcel reach the inversion on this grid (with b = 0.8 it stops at its base). Mirror depths
    # 16000 (inversion sigma 0.6375), 20000 (0.60) and 24000 Pa (0.5625), each in the middle of a band of inversion
    # heights with the same depth.
    weak = R.tiedtke_bechtold_shallow_parameters(NF; test_entrainment_b = 0.1)
    for (sigma_i, depth, classified_shallow) in ((0.6375, 16000, true), (0.60, 20000, true), (0.5625, 24000, false))
        column = build_column(NF, case(trade(sigma_i); grid = :exact4000, p_s = 100000.0))
        out = run_stage1b(column, par, lagged(NF); shallow = weak)
        excess = column_excess(column, out.work, lagged(NF), weak, par)
        ascent = column_shallow_ascent(column, out.work, excess, weak, par)
        @test ascent.has_base
        @test column.p_half[ascent.j_base] - column.p_half[ascent.j_top] === NF(depth)
        @test R.tb_is_shallow_cloud(column.p_half[ascent.j_base], column.p_half[ascent.j_top], par) ==
              classified_shallow
        # In the column: a shallow cloud (type 2), or a shallow candidate replaced by a deeper n-1 cloud (D57), or
        # (above 200 hPa) no shallow candidate at all.
        candidate = out.res.convection_type == 2 || out.res.shallow_candidate_replaced
        @test candidate == classified_shallow
        depth == 16000 && @test out.res.convection_type == 2
    end
    # K <= 0 at the LCL gives none: level 23 warmed by 2 K. The surface parcel reaches saturation at half level 24
    # (its LCL, as in the unchanged profile) with K < 0 there (mirror K = -1.9 m2 s-2), so there is no base; no deep
    # cloud either.
    reference = build_column(NF, case(trade(0.70)))
    @test run_stage1b(reference, par, lagged(NF)).res.shallow_test_base_half_level == 24
    warmed = with_level(reference, 23; dT = 2)
    out = run_stage1b(warmed, par, lagged(NF))
    @test out.res.shallow_test_base_half_level == 0
    @test out.res.convection_type == 0
    @test all(iszero, out.dT) && all(iszero, out.dq) && all(iszero, out.P)
    n = length(warmed.T)
    heights = [(warmed.geopotential[k] - out.work.geopotential_half[n + 1]) / par.gravity for k in 1:n]
    traced = traced_shallow_ascent(warmed, out.work, column_excess(warmed, out.work, lagged(NF), shallow, par),
                                   shallow, par, heights)
    last_step = traced.steps[end]
    @test !traced.has_base
    @test last_step.j == 24 && last_step.step.saturated && !(last_step.K > 0)
    # D56: a surface-parcel cloud deeper than 200 hPa with no deep cloud from n-1 gives none. trade(0.70) with a moist
    # lowest level (q_n = 19.5 g/kg, RH 95%): the surface parcel saturates at n-1/2 and reaches the inversion,
    # 22.5 kPa deep, while the drier n-1 parcel stays 18.8 kPa deep (mirror; the same for q_n 19 to 20 g/kg).
    moist = with_level(reference, n; q = 0.0195)
    out = run_stage1b(moist, par, lagged(NF))
    @test out.res.shallow_test_base_half_level == n
    @test out.res.shallow_test_depth > NF(20000)
    @test out.res.max_test_cloud_depth <= NF(20000)
    @test out.res.convection_type == 0
    @test !out.res.shallow_candidate_replaced
    @test all(iszero, out.dT) && all(iszero, out.dq) && all(iszero, out.P) && all(iszero, out.work.mass_flux_half)
end

@testset "T18 order: one cloud type, deep replaces the shallow candidate (D57) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    # trade(0.60): the surface parcel gives an 18.8 kPa shallow candidate, the n-1 parcel a 26.3 kPa deep cloud.
    column = build_column(NF, case(trade(0.60)))
    n = length(column.T)
    out = run_stage1b(column, par, lagged(NF))
    @test out.res.convection_type == 1
    @test out.res.shallow_candidate_replaced
    @test out.res.departure_level == n - 1
    @test out.res.deep_test_depth > NF(20000)
    @test out.res.PCAPE > 0 && out.res.tau > 0                          # the deep closure, not the shallow one
    for name in SHALLOW_ONLY_FIELDS
        @test iszero(getproperty(out.res, name))
    end
    # The deep cloud's own base and top (D57: base, top, w, s, q, l and type all from the deep parcel).
    @test out.res.cloud_base_half_level == 24
    @test out.res.tb_top_full_level < out.res.cloud_base_half_level - 1
    # Exactly one type per column over every column of this file's families.
    columns = Any[build_column(NF, case(trade(s))) for s in (0.75, 0.70, 0.65, 0.60)]
    append!(columns, [build_column(NF, case(TROPICAL)), build_column(NF, case(MIDLATITUDE; p_s = 100000.0)),
                      build_column(NF, case(STABLE)), build_column(NF, case(PLATEAU; p_s = 60000.0, z_s = 4000.0,
                                                                                land = true))])
    for c in columns
        r = run_stage1b(c, par, lagged(NF)).res
        @test r.convection_type in (0, 1, 2)
        if r.convection_type == 1
            @test all(name -> iszero(getproperty(r, name)), SHALLOW_ONLY_FIELDS)
            @test r.departure_level <= n - 1
        elseif r.convection_type == 2
            @test r.departure_level == length(c.T) && !r.shallow_candidate_replaced
            @test r.PCAPE == 0 && r.tau == 0 && r.deep_test_depth == 0
        else
            @test r.cloud_base_mass_flux == 0 && r.rain_mass_flux_surface == 0
        end
    end
end

@testset "T19 shallow updraught: E, D1, D2 and the D65 guard (D62-D65) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    shallow = R.tiedtke_bechtold_shallow_parameters(NF)
    cp = par.heat_capacity
    tol = local_tolerance(NF)
    geopotential_lower = NF(2.6e4)
    geopotential_upper = geopotential_lower + NF(9.81) * NF(400)
    layer_args(B, K, lfc, rh, f_scale, dT_upper) = (NF(0.04), cp * NF(282) + geopotential_lower, NF(7.5e-3),
        NF(1e-3), K, B, lfc, rh, f_scale, NF(400), NF(280), NF(6e-3), geopotential_lower, NF(277) + dT_upper,
        NF(5e-3), geopotential_upper, NF(70000), false)
    # D62: shallow E is twice the deep D17 E for the same inputs; D63: D1 = (1.6 - RH) E, zero where E is zero.
    organized = 0
    for K in NF.((0.5, 5, 20)), B in NF.((-0.5, 0.5)), lfc in (false, true), rh in NF.((0.6, 0.9)),
        f_scale in NF.((1, 0.1)), dT_upper in NF.((-2, 0, 2))
        args = layer_args(B, K, lfc, rh, f_scale, dT_upper)
        deep = R.tb_updraught_layer(args..., par)
        sh = R.tb_shallow_updraught_layer(args..., shallow, par)
        @test isapprox(sh.E, 2 * deep.E; rtol = 4 * eps(NF))
        @test (sh.E > 0) == (B > 0)
        @test isapprox(sh.D1, (NF(1.6) - rh) * sh.E; rtol = 4 * eps(NF))
        @test iszero(sh.E) == iszero(sh.D1)
        @test !sh.mass_guard                                   # RH >= 0.6: M' >= M, no guard
        sh.top && continue
        organized += sh.D2 > 0
        # D64 = D19: E and D2 never in one layer, no D2 below the LFC, 0 <= D2 <= M', M >= 0.
        @test iszero(sh.E) || iszero(sh.D2)
        @test lfc || iszero(sh.D2)
        @test 0 <= sh.D2 <= sh.M_mixed
        @test sh.M_upper >= 0
    end
    @test organized > 0
    # D65 on one layer: RH = 0, f_scale = 1, dz = 500 m, B > 0: E/M = 2 x 1.75e-3 x 1.3 x 500 = 2.275 and
    # M' = M (1 - 0.6 x 2.275) < 0. The layer is the top at once: nothing is mixed, adjusted, converted or integrated.
    M, s, qv, l, K, B = NF(0.04), cp * NF(282) + geopotential_lower, NF(7.5e-3), NF(1e-3), NF(5), NF(0.5)
    guarded = R.tb_shallow_updraught_layer(M, s, qv, l, K, B, true, zero(NF), one(NF), NF(500), NF(280), NF(6e-3),
                                           geopotential_lower, NF(277), NF(5e-3), geopotential_upper, NF(70000),
                                           false, shallow, par)
    @test guarded.mass_guard && guarded.top
    @test guarded.M_mixed <= 0 && guarded.M_upper == 0 && guarded.D2 == 0
    @test guarded.G == 0 && guarded.condensation == 0
    @test (guarded.s, guarded.q, guarded.l, guarded.K) === (s, qv, l, K)
    # In a column: trade(0.65). D1 = (1.6 - RH) E in every cloud layer, M >= 0 everywhere.
    column = build_column(NF, case(trade(0.65)))
    out = run_stage1b(column, par, lagged(NF))
    @test out.res.convection_type == 2
    w = out.work
    @test all(>=(0), w.mass_flux_half)
    k_top, j_base = out.res.tb_top_full_level, out.res.cloud_base_half_level
    entraining = 0
    for k in (k_top + 1):(j_base - 1)
        rh = R.tb_relative_humidity(column.T[k], column.q[k], column.p_full[k], par)
        E, D1, D2 = w.entrainment[k], w.turbulent_detrainment[k], w.organized_detrainment[k]
        entraining += E > 0
        @test abs(D1 - (NF(1.6) - rh) * E) <= tol * (D1 + E)
        @test iszero(E) == iszero(D1)
        @test iszero(E) || iszero(D2)
        scale = w.mass_flux_half[k + 1] + E + D1 + D2
        @test abs(w.mass_flux_half[k] - (w.mass_flux_half[k + 1] + E - D1 - D2)) <= tol * scale
    end
    @test entraining > 0
    # Finding 11 (recorded): the shallow cloud-base buoyancy is negative here, so the first cloud layer neither
    # entrains nor detrains.
    @test w.entrainment[j_base - 1] == 0 && w.turbulent_detrainment[j_base - 1] == 0
    # D65 in a column: trade(0.70) with a dry level 23 (q = 0) above the base half level 24. The dry layer makes the
    # base buoyant (B about +1.8 K), so layer 23 entrains with RH = 0 and f_scale = 1 (E/M = 1.70, mirror) and
    # M' = M - 0.6 E is about -2% of M: layer 23 becomes the top layer before any mixing.
    dry = with_level(build_column(NF, case(trade(0.70))), 23; q = 0)
    out = run_stage1b(dry, par, lagged(NF))
    n = length(dry.T)
    w = out.work
    g = par.gravity
    @test out.res.convection_type == 2
    @test out.res.shallow_mass_guard_top
    @test out.res.cloud_base_half_level == 24
    @test out.res.tb_top_full_level == 23
    @test w.mass_flux_half[23] == 0 && w.mass_flux_half[24] > 0
    @test all(>=(0), w.mass_flux_half)
    @test w.entrainment[23] == 0 && w.turbulent_detrainment[23] == 0       # E and D1 only formed M'
    @test isapprox(w.organized_detrainment[23], w.mass_flux_half[24]; rtol = tol)
    @test w.generation[23] == 0                                            # no condensate converted, no rain
    @test out.P[24] == out.P[23]
    @test w.total_water_flux_half[23] == 0 && w.liquid_static_energy_flux_half[23] == 0
    dp23 = dry.p_half[24] - dry.p_half[23]
    # The only flux change in layer 23 is the detrainment of M_{23+1/2}: dq_23 = (g/dp) F_qt(24), no G.
    @test isapprox(out.dq[23], g / dp23 * w.total_water_flux_half[24]; rtol = tol)
    @test all(iszero, out.dq[1:22]) && all(iszero, out.dT[1:22])
end

@testset "T28 shallow test-parcel condensate removal (D52) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    shallow = R.tiedtke_bechtold_shallow_parameters(NF)
    column = build_column(NF, case(trade(0.70)))
    n = length(column.T)
    out = run_stage1b(column, par, lagged(NF))
    excess = column_excess(column, out.work, lagged(NF), shallow, par)
    heights = [(column.geopotential[k] - out.work.geopotential_half[n + 1]) / par.gravity for k in 1:n]
    # With the column's own heights the test-side ascent is the column's ascent, operation for operation.
    ascent = column_shallow_ascent(column, out.work, excess, shallow, par)
    traced = traced_shallow_ascent(column, out.work, excess, shallow, par, heights)
    @test isequal((ascent.has_base, ascent.j_base, ascent.j_top, ascent.s_base, ascent.q_base, ascent.l_base,
                   ascent.K_base),
                  (traced.has_base, traced.j_base, traced.j_top, traced.s_base, traced.q_base, traced.l_base,
                   traced.K_base))
    # At every level the parcel carries exactly half of its adjusted condensate (carried l plus the new condensate),
    # as the deep D11 test.
    condensing = 0
    for st in traced.steps
        @test st.step.l == par.test_condensate_removal * st.step.l_adjusted
        if st.step.q_mixed > st.step.q
            condensing += 1
            @test isapprox(st.step.l, (st.step.l_mixed + (st.step.q_mixed - st.step.q)) / 2; rtol = sqrt(eps(NF)))
        end
    end
    @test condensing >= 2
    # A variant without removal (all condensate carried) loads the parcel more: its test cloud stops at its base
    # (depth 0) instead of 11.3 kPa (mirror), so the rule acts on the class depth.
    keep_all = R.tiedtke_bechtold_parameters(NF; test_condensate_removal = 1)
    ascent_keep = column_shallow_ascent(column, out.work, excess, shallow, keep_all)
    depth(a) = column.p_half[a.j_base] - column.p_half[a.j_top]
    @test ascent_keep.has_base
    @test depth(ascent_keep) != depth(ascent)
    @test depth(ascent) == out.res.shallow_test_depth
end

@testset "T29 non-upward lagged buoyancy flux (D53) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    shallow = R.tiedtke_bechtold_shallow_parameters(NF)
    column = build_column(NF, case(trade(0.75)))
    n = length(column.T)
    floor_w = NF(1.2) * cbrt(NF(0.1)^3)
    # H < 0, F_q = 0: downward buoyancy flux. The shallow test still runs (it finds its base at half level 24), w* is
    # the u* floor, the excess follows D46 with no floor, and the counter is one. S_sub = H < 0 then stops it (D73).
    out = run_stage1b(column, par, lagged(NF, -5.0, 0.0))
    rho = column.p_full[n] / (par.R_dry * column.T[n] * (1 + par.mu_virt_temp * column.q[n]))
    @test out.res.shallow_nonupward_buoyancy_flux == 1
    @test out.res.shallow_test_base_half_level == 24
    @test isapprox(out.res.surface_parcel_wstar, floor_w; rtol = 8 * eps(NF))
    @test isapprox(out.res.surface_parcel_delta_T, NF(1.5) * NF(-5) / (rho * par.heat_capacity * floor_w);
                   rtol = sqrt(eps(NF)))
    @test out.res.surface_parcel_delta_T < 0
    @test out.res.surface_parcel_delta_q == 0
    @test out.res.shallow_nonpositive_source == 1 && out.res.convection_type == 0
    # A zero buoyancy flux also counts as non-upward; an upward flux leaves the counter unchanged.
    @test run_stage1b(column, par, lagged(NF, 0.0, 0.0)).res.shallow_nonupward_buoyancy_flux == 1
    @test run_stage1b(column, par, lagged(NF)).res.shallow_nonupward_buoyancy_flux == 0
    @test run_stage1b(column, par, unset(NF)).res.shallow_nonupward_buoyancy_flux == 0     # no test when unset
end

@testset "T20 HOST closure (D68, D69, D71, D74) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    g, cp, Lv = par.gravity, par.heat_capacity, par.latent_heat_condensation
    flux = lagged(NF)
    F_h = flux.sensible_heat_flux + Lv * flux.humidity_flux
    # Base at n-1/2: trade(0.775) with q_n = 19.5 g/kg (mirror: base half level 27, 15.0 kPa, no limiter; the same
    # for q_n 19 to 20 g/kg). And a base higher up: trade(0.70) (base half level 24).
    base_n = with_level(build_column(NF, case(trade(0.775))), 27; q = 0.0195)
    for (column, base_at_n) in ((base_n, true), (build_column(NF, case(trade(0.70))), false))
        n = length(column.T)
        out = run_stage1b(column, par, flux)
        r = out.res
        @test r.convection_type == 2
        @test (r.cloud_base_half_level == n) == base_at_n
        @test r.f_CFL == 1 && r.f_pos == 1
        @test r.final_factor == r.closure_factor
        # D71: S_sub = H + L_v F_q (L_v once); D68: M_base Delta h_base = S_sub before limiters; D74: final = M_base.
        @test r.shallow_source == F_h
        @test r.shallow_large_scale_source == 0
        @test r.shallow_moist_static_energy_jump > 0
        @test isapprox(r.shallow_closure_mass_flux * r.shallow_moist_static_energy_jump, F_h; rtol = 4 * eps(NF))
        @test isapprox(r.cloud_base_mass_flux, r.shallow_closure_mass_flux; rtol = 4 * eps(NF))
        @test r.cloud_base_mass_flux > 0
        # D69 with D39: the cloud-base net flux F_sl + L_v F_qt equals M_base Delta h_base.
        j = r.cloud_base_half_level
        F_sl, F_qt = out.work.liquid_static_energy_flux_half[j], out.work.total_water_flux_half[j]
        @test abs(F_sl + Lv * F_qt - r.shallow_closure_mass_flux * r.shallow_moist_static_energy_jump) <=
              64 * eps(NF) * (abs(F_sl) + Lv * abs(F_qt))
        if base_at_n
            # My algebra (D39, D67): with the base at n-1/2, layer n gets the whole sub-cloud divergence, so
            # (c_p dT + L_v dq) dp_n/g = -(H + L_v F_q).
            mass = (column.p_half[n + 1] - column.p_half[n]) / g
            layer_n = (cp * out.dT[n] + Lv * out.dq[n]) * mass
            @test abs(layer_n + F_h) <= 64 * eps(NF) * (cp * abs(out.dT[n]) + Lv * abs(out.dq[n])) * mass
            @test out.dq[n] < 0                                            # layer n is ventilated
        end
    end
end

@testset "T21 guards: non-positive jump or source (D73) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    function assert_zero_column(out)
        @test out.res.convection_type == 0
        @test out.res.tb_top_full_level == length(out.dT) + 1
        @test out.res.rain_mass_flux_surface == 0 && out.res.cloud_base_mass_flux == 0
        @test out.res.shallow_closure_mass_flux == 0
        @test all(iszero, out.dT) && all(iszero, out.dq) && all(iszero, out.P)
        for name in (:mass_flux_half, :total_water_flux_half, :liquid_static_energy_flux_half, :entrainment,
                     :turbulent_detrainment, :organized_detrainment, :generation, :rain_evaporation,
                     :detrained_condensate)
            @test all(iszero, getproperty(out.work, name))
        end
    end
    reference = build_column(NF, case(trade(0.75)))
    # Delta h_base <= 0 with positive F_h: level 23 saturated (RH 100%) raises the base half-level environment h
    # above the parcel's (mirror Delta h_base = -522 J/kg) while the base stays at half level 24.
    T23, p23 = reference.T[23], reference.p_full[23]
    saturated = with_level(reference, 23; q = R.tb_saturation_humidity(T23, p23, par))
    out = run_stage1b(saturated, par, lagged(NF))
    @test out.res.shallow_test_base_half_level == 24
    @test out.res.shallow_source > 0
    @test out.res.shallow_moist_static_energy_jump <= 0
    @test (out.res.shallow_nonpositive_jump, out.res.shallow_nonpositive_source) == (1, 0)
    assert_zero_column(out)
    # F_h <= 0: H = -5 W m-2, F_q = 1e-6 kg m-2 s-1 (F_h = -2.5 W m-2), and F_h = 0.
    for (H, F_q) in ((-5.0, 1.0e-6), (0.0, 0.0))
        out = run_stage1b(reference, par, lagged(NF, H, F_q))
        @test out.res.shallow_test_base_half_level == 24
        @test out.res.shallow_moist_static_energy_jump > 0
        @test (out.res.shallow_nonpositive_jump, out.res.shallow_nonpositive_source) == (0, 1)
        assert_zero_column(out)
    end
    # A positive source and jump: no count.
    out = run_stage1b(reference, par, lagged(NF))
    @test (out.res.shallow_nonpositive_jump, out.res.shallow_nonpositive_source) == (0, 0)
    @test out.res.convection_type == 2
end

@testset "T22 SCM closure separation and S_LS integration (D72, G6) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    g, cp, Lv = par.gravity, par.heat_capacity, par.latent_heat_condensation
    scm(column, dT, dq) = R.TiedtkeBechtoldSCMClosure(fill(NF(dT), length(column.T)), fill(NF(dq), length(column.T)))
    columns = (build_column(NF, case(trade(0.65))),
               with_level(build_column(NF, case(trade(0.775))), 27; q = 0.0195),
               build_column(NF, case(TROPICAL)), build_column(NF, case(STABLE)))
    # S_LS = 0: the SCM closure equals the HOST closure bit for bit.
    for column in columns, flux in (lagged(NF), lagged(NF, -5.0, 1.0e-6), unset(NF))
        host = run_stage1b(column, par, flux)
        zero_scm = run_stage1b(column, par, flux; closure = scm(column, 0, 0))
        @test isequal(host.res, zero_scm.res)
        @test isequal((host.dT, host.dq, host.P), (zero_scm.dT, zero_scm.dq, zero_scm.P))
        @test all(isequal(getfield(host.work, f), getfield(zero_scm.work, f)) for f in keys(host.work))
    end
    column = columns[1]
    n = length(column.T)
    host = run_stage1b(column, par, lagged(NF))
    j_base = host.res.cloud_base_half_level
    p_s, p_base = column.p_half[n + 1], column.p_half[j_base]
    F_h = host.res.shallow_source
    # Term (1): a uniform humidity tendency X gives S_LS = L_v X (p_s - p_base)/g.
    X = NF(1e-8)                                                         # 0.86 g/kg/day
    out = run_stage1b(column, par, lagged(NF); closure = scm(column, 0, X))
    @test out.res.cloud_base_half_level == j_base
    @test isapprox(out.res.shallow_large_scale_source, Lv * X * (p_s - p_base) / g; rtol = closure_tolerance(NF, n))
    @test isapprox(out.res.shallow_source, F_h + out.res.shallow_large_scale_source; rtol = 4 * eps(NF))
    @test isapprox(out.res.shallow_closure_mass_flux * out.res.shallow_moist_static_energy_jump,
                   out.res.shallow_source; rtol = 4 * eps(NF))
    @test out.res.shallow_closure_mass_flux > host.res.shallow_closure_mass_flux
    # Term (2): a uniform radiative theta tendency Y converted with Pi = (p/1e5)^kappa gives c_p (1/g) int Pi Y dp.
    # The column integrates layer by layer, so the closed form holds to the midpoint-rule error (about 1e-5 here).
    Y = NF(-2 / 86400)                                                   # -2 K/day
    kappa = par.R_dry / cp
    radiative = [(column.p_full[k] / NF(1e5))^kappa * Y for k in 1:n]
    out = run_stage1b(column, par, lagged(NF); closure = R.TiedtkeBechtoldSCMClosure(radiative, zeros(NF, n)))
    closed_form = cp * Y / (g * (1 + kappa) * NF(1e5)^kappa) * (p_s^(1 + kappa) - p_base^(1 + kappa))
    @test isapprox(out.res.shallow_large_scale_source, closed_form; rtol = NF(1e-3))
    @test out.res.shallow_large_scale_source < 0
    # Only the layers below the cloud base count: a tendency confined above the base gives S_LS = 0.
    above = zeros(NF, n)
    above[1:(j_base - 1)] .= NF(1e-6)
    out = run_stage1b(column, par, lagged(NF); closure = R.TiedtkeBechtoldSCMClosure(zeros(NF, n), above))
    @test out.res.shallow_large_scale_source == 0
    @test isequal(out.res.shallow_closure_mass_flux, host.res.shallow_closure_mass_flux)
    # A strong large-scale sink makes S_sub <= 0 and the D73 guard applies in the SCM path too.
    out = run_stage1b(column, par, lagged(NF); closure = scm(column, 0, -1e-5))
    @test out.res.shallow_source <= 0 && out.res.shallow_nonpositive_source == 1 && out.res.convection_type == 0
end

@testset "T23 n = 0 fallback equals stage 1a bit for bit (D81), and stage 1B deep from n-2 $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    cases = (case(TROPICAL), case(MIDLATITUDE; p_s = 100000.0, land = true),
             case(PLATEAU; p_s = 60000.0, z_s = 4000.0, land = true), case(STABLE), case(trade(0.65)),
             case(TROPICAL; grid = :stretched, p_s = 100000.0))
    for c in cases, prefill in (zero(NF), NF(NaN))
        column = build_column(NF, c)
        stage1a = run_stage1a(column, par; prefill)
        fallback = run_stage1b(column, par, unset(NF); prefill)
        @test isequal((stage1a.dT, stage1a.dq, stage1a.P), (fallback.dT, fallback.dq, fallback.P))
        @test all(isequal(getfield(stage1a.work, f), getfield(fallback.work, f)) for f in keys(stage1a.work))
        @test all(name -> isequal(getproperty(stage1a.res, name), getproperty(fallback.res, name)),
                  keys(stage1a.res))
        @test fallback.res.lagged_flux_unset == 1
        @test fallback.res.trigger_path == R.TIEDTKE_BECHTOLD_TRIGGER_PATH_STAGE1A_FALLBACK
        @test fallback.res.convection_type in (0, 1)                       # never shallow while unset
        @test all(name -> iszero(getproperty(fallback.res, name)), SHALLOW_ONLY_FIELDS)
        @test fallback.res.surface_parcel_wstar == 0
    end
    @test run_stage1b(build_column(NF, case(trade(0.65))), par, lagged(NF)).res.lagged_flux_unset == 0
    # The stage 1B deep path repeats the stage 1a blocks: with a strongly downward lagged flux (H = -1000 W m-2,
    # Delta T = -10.9 K, w* at its floor) neither the surface parcel nor the n-1 parcel reaches a base (mirror K < -3.5
    # after one step), so the deep cloud leaves from n-2 exactly as in stage 1a and every output matches bit for bit.
    for c in (case(TROPICAL), case(MIDLATITUDE; p_s = 100000.0, land = true))
        column = build_column(NF, c)
        n = length(column.T)
        stage1a = run_stage1a(column, par)
        stage1b = run_stage1b(column, par, lagged(NF, -1000.0, 0.0))
        @test stage1b.res.convection_type == stage1a.res.convection_type == 1
        @test stage1b.res.departure_level == stage1a.res.departure_level == n - 2
        @test stage1b.res.shallow_test_base_half_level == 0
        @test stage1b.res.departures_scanned == stage1a.res.departures_scanned + 1
        @test isequal((stage1a.dT, stage1a.dq, stage1a.P), (stage1b.dT, stage1b.dq, stage1b.P))
        @test all(isequal(getfield(stage1a.work, f), getfield(stage1b.work, f)) for f in keys(stage1a.work))
        for name in keys(stage1a.res)
            name in (:departures_scanned, :departure_first) && continue
            @test isequal(getproperty(stage1a.res, name), getproperty(stage1b.res, name))
        end
    end
end

@testset "T24 limiter on a shallow column (D76, D41 r15) $NF" for NF in PRECISIONS
    g = R.tiedtke_bechtold_parameters(NF).gravity
    column = build_column(NF, case(trade(0.75)))
    n = length(column.T)
    strong = lagged(NF, 50.0, 2.0e-4)                                    # 550 W m-2: M_base about 0.22 kg m-2 s-1
    cfl_numbers(out, factor, dt) = [out.work.mass_flux_half[j] * g * factor * NF(dt) /
        min(column.p_half[j] - column.p_half[j - 1], column.p_half[j + 1] - column.p_half[j]) for j in 2:n]
    results = Dict{Int, Any}()
    for factor in (2, 1)
        par = R.tiedtke_bechtold_parameters(NF; limiter_timestep_factor = factor)
        # At the host step the strong case already binds with the leapfrog factor 2 (mirror f_CFL 0.855) and not
        # with factor 1 (1.71, capped at 1).
        host_step = run_stage1b(column, par, strong)
        @test host_step.res.convection_type == 2
        @test (host_step.res.f_CFL < 1) == (factor == 2)
        # Forced: a 1 h step (the closed mass flux does not depend on dt, the limit does).
        out = run_stage1b(column, par, strong; dt = 3600)
        results[factor] = out
        @test out.res.convection_type == 2
        @test out.res.f_CFL < 1
        @test isapprox(out.res.final_factor, out.res.closure_factor * min(out.res.f_CFL, out.res.f_pos);
                       rtol = local_tolerance(NF))
        numbers = cfl_numbers(out, factor, 3600)
        @test maximum(numbers) <= 1 + 64 * eps(NF)
        out.res.f_CFL < out.res.f_pos && @test isapprox(maximum(numbers), one(NF); rtol = 64 * eps(NF))
        l = ledgers(out, column, par)
        @test abs(l.water) <= closure_tolerance(NF, n) * l.water_gross
        @test abs(l.transform) <= closure_tolerance(NF, n) * l.transform_gross
    end
    # One column factor over limiter_timestep_factor * dt: halving the applied step doubles the factor.
    @test isapprox(results[1].res.final_factor, 2 * results[2].res.final_factor; rtol = 64 * eps(NF))
    @test results[1].res.closure_factor == results[2].res.closure_factor
end

@testset "T25 ledgers on shallow columns with and without rain (D75, D67: T2, T11, T13, T15) $NF" for NF in PRECISIONS
    par = R.tiedtke_bechtold_parameters(NF)
    no_rain = R.tiedtke_bechtold_parameters(NF; l_threshold_water = 1, l_threshold_land = 1)   # conversion off
    g = par.gravity
    columns = (build_column(NF, case(trade(0.75))), build_column(NF, case(trade(0.70))),
               build_column(NF, case(trade(0.65))), with_level(build_column(NF, case(trade(0.775))), 27; q = 0.0195),
               build_column(NF, case(trade(0.70); land = true)))
    raining = 0
    for column in columns, (p, rain_allowed) in ((par, true), (no_rain, false))
        out = run_stage1b(column, p, lagged(NF))
        n = length(column.T)
        @test out.res.convection_type == 2
        G, e, P = out.work.generation, out.work.rain_evaporation, out.P
        if rain_allowed
            raining += sum(G) > 0
        else
            @test all(iszero, G) && all(iszero, P) && out.res.rain_mass_flux_surface == 0
        end
        # T2: both identities.
        l = ledgers(out, column, p)
        @test abs(l.water) <= closure_tolerance(NF, n) * l.water_gross
        @test abs(l.transform) <= closure_tolerance(NF, n) * l.transform_gross
        # T11: the column water loss equals the surface rain.
        loss = -sum(out.dq[k] * (column.p_half[k + 1] - column.p_half[k]) / g for k in 1:n)
        gross = sum(abs(out.dq[k]) * (column.p_half[k + 1] - column.p_half[k]) / g for k in 1:n)
        @test abs(loss - out.res.rain_mass_flux_surface) <= closure_tolerance(NF, n) * gross
        # T13: the local rain-flux profile.
        @test P[1] == 0 && P[n + 1] == out.res.rain_mass_flux_surface && all(>=(0), P)
        j_base = out.res.cloud_base_half_level
        atol = local_tolerance(NF) * max(maximum(P), floatmin(NF))
        for k in 1:n
            mass = (column.p_half[k + 1] - column.p_half[k]) / g
            @test abs((P[k + 1] - P[k]) - (G[k] - e[k]) * mass) <= atol
            @test iszero(G[k]) || iszero(e[k])
            @test k < j_base || G[k] == 0
            @test k >= j_base || e[k] == 0
        end
        # T15 / D67: sub-cloud interfaces carry Z_j times the cloud-base net flux, zero at the surface.
        F_qt, F_sl, M = (out.work.total_water_flux_half, out.work.liquid_static_energy_flux_half,
                         out.work.mass_flux_half)
        p_s, p_base = column.p_half[n + 1], column.p_half[j_base]
        for j in (j_base + 1):n
            Z = (p_s - column.p_half[j]) / (p_s - p_base)
            @test abs(F_qt[j] - Z * F_qt[j_base]) <= local_tolerance(NF) * abs(F_qt[j_base])
            @test abs(F_sl[j] - Z * F_sl[j_base]) <= local_tolerance(NF) * abs(F_sl[j_base])
            @test abs(M[j] - Z * M[j_base]) <= local_tolerance(NF) * M[j_base]
        end
        @test F_qt[n + 1] == 0 && F_sl[n + 1] == 0 && M[n + 1] == 0
        @test F_qt[1] == 0 && F_sl[1] == 0 && M[1] == 0
        # The shallow cloud ventilates the sub-cloud layer: upward cloud-base total-water flux. Without rain layer n
        # dries. With rain, D29 sub-cloud evaporation can outweigh the scaled flux divergence in layer n (mirror: the
        # sea columns trade(0.70) and trade(0.65) moisten layer n), so no sign is asserted there.
        @test F_qt[j_base] > 0
        rain_allowed || @test out.dq[n] < 0
        # D79 at column level: the type flag and a full-level top above the base (host cloud_top is stage 2).
        @test 1 <= out.res.tb_top_full_level < j_base
    end
    @test raining >= 3                                    # mirror: generated rain in every default trade column
end

end # top-level testset
end # module

println("TIEDTKE_BECHTOLD_STAGE1B_TESTS_PASS")
