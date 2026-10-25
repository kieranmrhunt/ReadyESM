module PBLKProfileMomentumTests
# Opt-in K-profile momentum mixing (atmosphere_pbl_mix_momentum): configuration, refusals and provenance, the column
# algebra (momentum conservation, the backward-Euler kinetic-energy identity, a hand-solved column), the coupled
# column behaviour with 27 layers, the CPU form of the GPU interface, and a coupled save, restore and continue.
using ReadyESM, Test, LinearAlgebra, Random
const R = ReadyESM
const S = R.SpeedyWeather
const O = R.Oceananigans

fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
with(config; kwargs...) = R.ExperimentConfig(; merge(fields(config), (; kwargs...))...)
bits(x::Float32) = reinterpret(UInt32, x)
bits(x::Float64) = reinterpret(UInt64, x)
bits(x::Complex) = (bits(real(x)), bits(imag(x)))   # spectral prognostics are ComplexF32
same_bits(a, b) = size(a) == size(b) && eltype(a) == eltype(b) && all(bits(x) == bits(y) for (x, y) in zip(a, b))

function compact_l27(output; steps = 4, kwargs...)
    base = R.load_config(joinpath(R.PROJECT_ROOT, "config/production.yml"); check_input_files = false)
    era5 = joinpath(R.PROJECT_ROOT, "data/era5_initial_state")
    return with(base; device = :cpu, truncation = 7, nlayers = 27, ocean_grid = :idealized,
        ocean_nlongitude = 32, ocean_nlatitude = 16, ocean_nlayers = 4, ocean_initial_conditions = :analytic,
        ecco_initial_conditions_directory = "", duration_days = steps * 3600 / 86_400, analysis_start_day = 0.0,
        output_dir = output,
        era5_pressure_levels_path = joinpath(era5, "era5_pressure_levels_stratosphere_19930101T0000Z.nc"),
        era5_single_levels_path = joinpath(era5, "era5_single_levels_19930101T0000Z.nc"),
        era5_land_state_path = joinpath(era5, "era5_land_state_19930101T0000Z.nc"),
        atmosphere_vertical_diffusion = :holtslag_boville_k_profile_v1, kwargs...)
end

"""Call `scheme` on every column with zeroed T, q, u and v tendencies; return copies of the four tendencies."""
function column_tendencies(vars, model, scheme)
    t = vars.tendencies.grid
    for x in (t.temperature, t.humidity, t.u, t.v)
        fill!(x.data, 0)
    end
    for ij in 1:size(vars.grid.temperature_prev.data, 1)
        S.parameterization!(ij, vars, scheme, model)
    end
    return map(x -> copy(x.data), (t.temperature, t.humidity, t.u, t.v))
end

"""Hand-built backward-Euler matrix for the flux form with zero boundary fluxes."""
function implicit_matrix(a, dsigma, dt)
    n = length(dsigma)
    M = Matrix{Float64}(I, n, n)
    for k in 1:(n - 1)
        M[k, k] += dt * a[k] / dsigma[k]
        M[k, k + 1] -= dt * a[k] / dsigma[k]
        M[k + 1, k + 1] += dt * a[k] / dsigma[k + 1]
        M[k + 1, k] -= dt * a[k] / dsigma[k + 1]
    end
    return M
end

prognostic(vars) = map(n -> copy(Array(S.get_step(getproperty(vars.prognostic, n), 1).data)),
                       (:temperature, :humidity, :vorticity, :divergence, :pressure))

@testset "PBL K-profile momentum mixing" begin
    production = R.load_config(joinpath(R.PROJECT_ROOT, "config/production.yml"); check_input_files = false)
    grid = S.SpectralGrid(; NF = Float32, trunc = 7, nlayers = 8)
    on = with(production; atmosphere_vertical_diffusion = :holtslag_boville_k_profile_v1)

    @testset "bitwise comparison helper" begin
        a = ComplexF32[1 + 2im, -3.5f0 + 0.25f0im, 0 + 0im]
        @test same_bits(a, copy(a))
        b = copy(a); b[2] = complex(real(b[2]), nextfloat(imag(b[2])))
        @test !same_bits(a, b)                                            # changed imaginary part
        @test !same_bits(ComplexF32[complex(0f0, 0f0)], ComplexF32[complex(-0f0, 0f0)])   # signed zero, real
        @test !same_bits(ComplexF32[complex(0f0, 0f0)], ComplexF32[complex(0f0, -0f0)])   # signed zero, imaginary
        @test same_bits(ComplexF32[complex(NaN32, 1f0)], ComplexF32[complex(NaN32, 1f0)]) # same NaN bits
        @test !same_bits(a, a[1:2])                                       # different length
        @test !same_bits(reshape([a; a], 3, 2), reshape([a; a], 2, 3))    # different shape
        @test !same_bits(a, ComplexF64.(a))                               # different element type
        @test !same_bits(Float32[1, 2], Float64[1, 2])
        @test same_bits(Float32[1, -0f0], Float32[1, -0f0]) && !same_bits(Float32[0], Float32[-0f0])
    end

    @testset "configuration, provenance and refusals" begin
        @test production.atmosphere_pbl_mix_momentum == false
        default = R._atmosphere_vertical_diffusion(on, grid)
        @test default.mix_momentum == false
        @test R._atmosphere_vertical_diffusion_provenance(default) == R.PBL_K_PROFILE_PROVENANCE
        base = R.PBL_K_PROFILE_PROVENANCE
        for (kwargs, suffix) in (
                ((; atmosphere_pbl_mix_momentum = true), "_with_momentum"),
                ((; atmosphere_pbl_mix_momentum = true, atmosphere_pbl_mix_heat = false,
                    atmosphere_pbl_mix_humidity = false), "_momentum_only"),
                ((; atmosphere_pbl_mix_momentum = true, atmosphere_pbl_mix_humidity = false), "_heat_only_with_momentum"),
                ((; atmosphere_pbl_mix_humidity = false), "_heat_only"),
                ((; atmosphere_pbl_mix_heat = false), "_humidity_only"))
            config = with(on; kwargs...)
            @test R.validate(config; check_input_files = false) === config
            scheme = R._atmosphere_vertical_diffusion(config, grid)
            @test R._atmosphere_vertical_diffusion_provenance(scheme) == base * suffix
            @test scheme.mix_momentum == get(kwargs, :atmosphere_pbl_mix_momentum, false)
            @test isbitstype(typeof(scheme))
        end
        @test_throws ArgumentError R.validate(with(on; atmosphere_pbl_mix_heat = false, atmosphere_pbl_mix_humidity = false);
                                              check_input_files = false)
        @test_throws ArgumentError R.KProfileBoundaryLayerMixing(grid; mix_heat = false, mix_humidity = false)
        @test R.KProfileBoundaryLayerMixing(grid; mix_heat = false, mix_humidity = false, mix_momentum = true).mix_momentum
        # Momentum mixing needs the K-profile: the zero-operator control and the bulk Richardson option refuse it.
        @test_throws ArgumentError R.validate(with(production; atmosphere_pbl_mix_momentum = true); check_input_files = false)
        @test_throws ArgumentError R.validate(with(production; atmosphere_vertical_diffusion = :metric_bulk_richardson,
                                                   atmosphere_pbl_mix_momentum = true); check_input_files = false)
        @test R.validate(production; check_input_files = false) === production
        mktempdir() do dir
            lines = readlines(joinpath(R.PROJECT_ROOT, "config/production.yml"))
            index = findall(l -> startswith(l, "atmosphere_vertical_diffusion:"), lines)
            @test length(index) == 1
            lines[only(index)] = "atmosphere_vertical_diffusion: holtslag_boville_k_profile_v1"
            push!(lines, "atmosphere_pbl_mix_momentum: true", "atmosphere_pbl_mix_heat: false",
                  "atmosphere_pbl_mix_humidity: false")
            path = joinpath(dir, "momentum.yml"); write(path, join(lines, "\n") * "\n")
            parsed = R.load_config(path; check_input_files = false)
            @test parsed.atmosphere_pbl_mix_momentum && !parsed.atmosphere_pbl_mix_heat && !parsed.atmosphere_pbl_mix_humidity
        end
    end

    @testset "column algebra" begin
        rng = Random.Xoshiro(20261025)
        for NF in (Float32, Float64), trial in 1:20
            n = trial == 1 ? 3 : 27
            dsigma = NF.(rand(rng, n) .+ 0.2); dsigma ./= sum(dsigma)
            a = NF.(rand(rng, n - 1) .* (trial % 4 == 0 ? 0 : 1e-3))   # every fourth trial: a = 0
            dt = NF(1800)
            wind = NF.(10 .* randn(rng, n))
            c, d, x_new = zeros(NF, n), zeros(NF, n), zeros(NF, n)
            tendency = zeros(NF, 1, n)
            face, implicit = R._pbl_mix_wind!(tendency, wind, 1, x_new, a, dsigma, dt, c, d)
            tol = NF == Float32 ? 1e-5 : 1e-12
            scale = sum(abs.(tendency[1, :]) .* dsigma) + eps(NF)
            # Column momentum: the sigma-weighted sum of the tendency vanishes (zero boundary fluxes).
            @test abs(sum(tendency[1, :] .* dsigma)) <= tol * scale
            # Flux form equals (x_new - x_old) / dt, and the solve matches the hand-built matrix.
            exact = implicit_matrix(Float64.(a), Float64.(dsigma), Float64(dt)) \ Float64.(wind)
            @test maximum(abs, Float64.(x_new) .- exact) <= (NF == Float32 ? 1e-4 : 1e-10) * (maximum(abs, wind) + 1)
            @test maximum(abs, Float64.(tendency[1, :]) .- (exact .- Float64.(wind)) ./ Float64(dt)) <=
                  (NF == Float32 ? 1e-6 : 1e-12) * (maximum(abs, wind) + 1)
            # Kinetic energy: E_old - E_new = dt * face + implicit / 2 (per unit p_s/g), and both terms are >= 0.
            loss = sum(Float64.(dsigma) .* (Float64.(wind) .^ 2 .- exact .^ 2)) / 2
            @test face >= 0 && implicit >= 0
            @test abs(loss - (Float64(dt) * face + implicit / 2)) <= (NF == Float32 ? 1e-3 : 1e-9) * (abs(loss) + 1e-6)
            @test loss >= -tol
            # Bottom face, exactly: the lowest-layer tendency is a[n-1] (x_new[n-1] - x_new[n]) / dsigma[n], so its
            # sign is that of the solved gradient across the lowest face (down the solved gradient).
            bottom = a[n - 1] * (x_new[n - 1] - x_new[n]) / dsigma[n]
            @test tendency[1, n] == bottom
            @test sign(tendency[1, n]) == sign(x_new[n - 1] - x_new[n]) || a[n - 1] == 0
            if all(iszero, a)
                @test all(iszero, tendency) && x_new == wind && face == 0 && implicit == 0
            end
            # A uniform wind is unchanged.
            fill!(tendency, 0)
            uniform = fill(NF(7), n)
            R._pbl_mix_wind!(tendency, uniform, 1, x_new, a, dsigma, dt, c, d)
            @test maximum(abs, tendency) <= tol
        end
    end

    @testset "lowest-face rounding contract" begin
        z = 331f0
        h_raw = Float64(z) + 1e-7
        @test Float32(h_raw) == z
        @test R.pbl_diffusivity(Float64(z), h_raw, 0.3, 0.01, 300.0, 9.81) > 0        # unrounded: tiny but nonzero
        @test R.pbl_diffusivity(z, Float32(h_raw), 0.3f0, 0.01f0, 300f0, 9.81f0) == 0  # stored values: exactly zero
    end

    @testset "compact coupled model, 27 layers" begin
        mktempdir() do dir
            config = compact_l27(joinpath(dir, "out"))
            R.validate(config)
            simulation = R.build_dynamic_esm(config; steps = 4)
            model = simulation.model.atmosphere.model
            vars = simulation.model.atmosphere.variables
            O.run!(simulation)
            @test simulation.model.clock.iteration == 4
            sg, p = model.spectral_grid, vars.parameterizations
            default = model.vertical_diffusion
            @test default.mix_momentum == false
            dT0, dq0, du0, dv0 = column_tendencies(vars, model, default)
            @test all(iszero, du0) && all(iszero, dv0)
            @test all(iszero, Array(p.pbl_momentum_dissipation_face.data))
            @test all(iszero, Array(p.pbl_momentum_dissipation_implicit.data))
            with_momentum = R.KProfileBoundaryLayerMixing(sg; mix_momentum = true)
            @test S.initialize!(with_momentum, model) === nothing
            dT1, dq1, du1, dv1 = column_tendencies(vars, model, with_momentum)
            # The scalar tendencies are bitwise those of the scalar-only scheme.
            @test dT1 == dT0 && dq1 == dq0
            momentum_only = R.KProfileBoundaryLayerMixing(sg; mix_heat = false, mix_humidity = false, mix_momentum = true)
            dT2, dq2, du2, dv2 = column_tendencies(vars, model, momentum_only)
            @test all(iszero, dT2) && all(iszero, dq2) && du2 == du1 && dv2 == dv1
            npoints = size(du1, 1)
            active = [ij for ij in 1:npoints if any(!=(0), du1[ij, :])]
            println("PBL_MOMENTUM_ACTIVE_COLUMNS ", length(active), " of ", npoints)
            @test length(active) > npoints ÷ 10
            dsigma = Float64.(model.geometry.σ_levels_thick)
            Dface = Float64.(Array(p.pbl_momentum_dissipation_face.data))
            Dimpl = Float64.(Array(p.pbl_momentum_dissipation_implicit.data))
            @test all(isfinite, Dface) && all(>=(0), Dface) && all(isfinite, Dimpl) && all(>=(0), Dimpl)
            u, v = Float64.(Array(vars.grid.u_prev.data)), Float64.(Array(vars.grid.v_prev.data))
            ps = Float64.(Array(vars.grid.pressure_prev.data))
            gravity = Float64(model.planet.gravity)
            dt = 2 * Float64(model.time_stepping.Δt_sec)
            worst_momentum = worst_energy = 0.0
            for ij in active
                for (x, dx) in ((u, du1), (v, dv1))
                    total = sum(Float64.(dx[ij, :]) .* dsigma)
                    magnitude = sum(abs.(Float64.(dx[ij, :])) .* dsigma)
                    if magnitude == 0
                        @test total == 0      # a component with no tendency in this column
                    else
                        worst_momentum = max(worst_momentum, abs(total) / magnitude)
                    end
                end
                # Kinetic-energy loss of the isolated update, from the flux-form tendencies, against D_face + D_implicit.
                mass = ps[ij] / gravity
                unew, vnew = u[ij, :] .+ dt .* Float64.(du1[ij, :]), v[ij, :] .+ dt .* Float64.(dv1[ij, :])
                loss = mass / 2 * sum(dsigma .* (u[ij, :] .^ 2 .+ v[ij, :] .^ 2 .- unew .^ 2 .- vnew .^ 2)) / dt
                expected = Dface[ij] + Dimpl[ij]
                worst_energy = max(worst_energy, abs(loss - expected) / max(abs(expected), 1e-6))
            end
            println("PBL_MOMENTUM_COLUMN_CHECKS worst_momentum=$worst_momentum worst_energy=$worst_energy ",
                    "point_mean_D_face=", sum(Dface) / npoints, " point_mean_D_implicit=", sum(Dimpl) / npoints,
                    " (unweighted means over reduced-grid points)")
            @test worst_momentum < 1e-4
            @test worst_energy < 1e-2
            n = model.geometry.nlayers
            # Lowest-face diagnostics (v3): one rounded depth and face height for the coefficients and the reports.
            hs = Array(p.pbl_height.data); zf = Array(p.pbl_lowest_face_height.data)
            Kf = Array(p.pbl_lowest_face_diffusivity.data)
            @test all(isfinite, zf) && all(>(0), zf) && all(hs .>= zf)
            @test all(Kf[ij] == 0 for ij in eachindex(hs) if hs[ij] == zf[ij])     # implication only, not iff
            ustar_s, buoy_s = Array(p.pbl_friction_velocity.data), Array(p.pbl_buoyancy_flux.data)
            g32 = Float32(model.planet.gravity)
            for ij in eachindex(hs)
                thetav = R._PBLThetaV(vars.grid.temperature_prev, vars.grid.humidity_prev, model.geometry.σ_levels_full,
                                      model.atmosphere, Float32(model.atmosphere.κ), ij)
                @test Kf[ij] == Float32(R.pbl_diffusivity(zf[ij], hs[ij], ustar_s[ij], buoy_s[ij], thetav[n], g32)) *
                               with_momentum.diffusivity_scale
            end
            println("PBL_MOMENTUM_LOWEST_FACE active_fraction=", count(>(0), Kf) / length(Kf),
                    " clamped_fraction=", count(hs .== zf) / length(hs))
            # Fixture diagnostic only (not an invariant of the implicit solve, which the column algebra tests exactly):
            # the fraction of columns whose lowest-layer tendency has the sign of the old bottom-face difference.
            slower = [ij for ij in active if abs(u[ij, n]) < abs(u[ij, n - 1]) && sign(u[ij, n]) == sign(u[ij, n - 1])]
            agree = isempty(slower) ? NaN :
                count(sign(du1[ij, n]) == sign(u[ij, n - 1] - u[ij, n]) for ij in slower) / length(slower)
            println("PBL_MOMENTUM_OLD_GRADIENT_SIGN_FRACTION ", agree, " of ", length(slower))
            # The CPU form of the GPU interface: device model fields and plain arrays give bitwise the same result.
            device_model = (; geometry = model.geometry, planet = model.planet, atmosphere = model.atmosphere,
                            time_stepping = model.time_stepping, land_sea_mask = model.land_sea_mask)
            plain(x) = copy(Array(x.data))
            g = vars.grid
            device_vars = (
                grid = (; temperature_prev = plain(g.temperature_prev), humidity_prev = plain(g.humidity_prev),
                        u_prev = plain(g.u_prev), v_prev = plain(g.v_prev), pressure_prev = plain(g.pressure_prev),
                        geopotential = plain(g.geopotential)),
                parameterizations = (; (name => plain(getproperty(p, name)) for name in
                    (:pbl_height, :pbl_friction_velocity, :pbl_buoyancy_flux, :pbl_drag_coefficient,
                     :pbl_lowest_face_diffusivity, :pbl_lowest_face_height, :pbl_momentum_dissipation_face, :pbl_momentum_dissipation_implicit,
                     :pbl_height_work, :pbl_solution_work, :boundary_layer_drag))...),
                prognostic = (; land = (; sensible_heat_flux = plain(vars.prognostic.land.sensible_heat_flux),
                                        surface_humidity_flux = plain(vars.prognostic.land.surface_humidity_flux)),
                              ocean = (; sensible_heat_flux = plain(vars.prognostic.ocean.sensible_heat_flux),
                                         surface_humidity_flux = plain(vars.prognostic.ocean.surface_humidity_flux))),
                scratch = (; grid = (; a = plain(vars.scratch.grid.a), b = plain(vars.scratch.grid.b))),
                tendencies = (; grid = (; temperature = zero(dT1), humidity = zero(dq1), u = zero(du1), v = zero(dv1))))
            for (scheme, expected) in ((with_momentum, (dT1, dq1, du1, dv1)), (momentum_only, (dT2, dq2, du2, dv2)))
                t = device_vars.tendencies.grid
                foreach(x -> fill!(x, 0), (t.temperature, t.humidity, t.u, t.v))
                for ij in 1:npoints
                    S.parameterization!(ij, device_vars, scheme, device_model)
                end
                @test (t.temperature, t.humidity, t.u, t.v) == expected
                @test device_vars.parameterizations.pbl_momentum_dissipation_face == Array(p.pbl_momentum_dissipation_face.data)
            end
        end
    end

    @testset "K-profile restart record" begin
        sim(vd) = (; model = (; vertical_diffusion = vd))
        zero_operator = R._atmosphere_vertical_diffusion(production, grid)
        @test !(zero_operator isa R.KProfileBoundaryLayerMixing)
        @test R._pbl_k_profile_restart_state(sim(zero_operator)) == (;)     # production checkpoints do not change
        @test R._validate_pbl_k_profile_restart(sim(zero_operator), (;)) === nothing
        scalar = R.KProfileBoundaryLayerMixing(grid)
        momentum = R.KProfileBoundaryLayerMixing(grid; mix_momentum = true)
        only = R.KProfileBoundaryLayerMixing(grid; mix_heat = false, mix_humidity = false, mix_momentum = true)
        scaled = R.KProfileBoundaryLayerMixing(grid; diffusivity_scale = 0.3f0)
        record(x) = R._pbl_k_profile_restart_state(sim(x))
        @test record(momentum) == (; pbl_k_profile_restart_schema = R.PBL_K_PROFILE_RESTART_SCHEMA,
            pbl_diffusivity_scale = 1.0, pbl_mix_heat = true, pbl_mix_humidity = true, pbl_mix_momentum = true)
        @test record(scaled).pbl_diffusivity_scale === Float64(0.3f0)      # exact value, not the rounded suffix
        for x in (scalar, momentum, only, scaled)
            @test R._validate_pbl_k_profile_restart(sim(x), record(x)) === nothing    # same mode
        end
        for (target, source) in ((scalar, momentum), (momentum, scalar), (momentum, only), (scalar, scaled))
            @test_throws ArgumentError R._validate_pbl_k_profile_restart(sim(target), record(source))
        end
        @test_throws ArgumentError R._validate_pbl_k_profile_restart(sim(zero_operator), record(scalar))  # on -> off
        @test_throws ArgumentError R._validate_pbl_k_profile_restart(sim(momentum), (;))                 # off -> on
        @test_throws ArgumentError R._validate_pbl_k_profile_restart(sim(momentum),
            merge(record(momentum), (; pbl_k_profile_restart_schema = "readiesm_pbl_k_profile_restart_v0")))
        @test_throws ArgumentError R._validate_pbl_k_profile_restart(sim(momentum),
            Base.structdiff(record(momentum), NamedTuple{(:pbl_mix_momentum,)}))
        @test_throws ArgumentError R._validate_pbl_k_profile_restart(sim(scaled),
            merge(record(scaled), (; pbl_diffusivity_scale = 0.3)))         # 0.3 differs from Float64(0.3f0)
    end

    @testset "coupled run, and save, restore and continue, momentum only" begin
        # The save happens inside the run (callback at iteration 2). The continuation follows the initialized-restart
        # stepping path of scripts/validate_coupled_restart.jl, because native run! would initialize a second time.
        mktempdir() do dir
            settings = (; atmosphere_pbl_mix_momentum = true, atmosphere_pbl_mix_heat = false,
                        atmosphere_pbl_mix_humidity = false, atmosphere_vertical_advection = :centered_full_temperature)
            reference = R.build_dynamic_esm(compact_l27(joinpath(dir, "reference"); settings...); steps = 4)
            @test reference.model.atmosphere.model.vertical_diffusion.mix_momentum
            path = joinpath(dir, "restart_step2.jld2")
            function save_at_two(sim)
                sim.model.clock.iteration == 2 && R.save_dynamic_restart_state(sim, path)
                return nothing
            end
            O.add_callback!(reference, save_at_two, O.IterationInterval(1); name = :momentum_restart_save)
            O.run!(reference)
            @test reference.model.clock.iteration == 4 && isfile(path)
            final = prognostic(reference.model.atmosphere.variables)
            @test all(x -> all(isfinite, x), final)    # a restart identity on a nonfinite state is no evidence
            continued = R.build_dynamic_esm(compact_l27(joinpath(dir, "continued"); settings...); steps = 4)
            R.restore_dynamic_restart_state!(continued, path)
            @test continued.model.clock.iteration == 2
            continued.stop_iteration = 4
            continued.running = true
            while continued.running
                O.time_step!(continued)
            end
            @test continued.model.clock.iteration == 4
            @test all(same_bits(x, y) for (x, y) in zip(prognostic(continued.model.atmosphere.variables), final))
            # On to off: a scalar-only K-profile configuration refuses the momentum checkpoint before any copy.
            scalar_settings = merge(settings, (; atmosphere_pbl_mix_momentum = false, atmosphere_pbl_mix_heat = true,
                                               atmosphere_pbl_mix_humidity = true))
            scalar_run = R.build_dynamic_esm(compact_l27(joinpath(dir, "scalar"); scalar_settings...); steps = 4)
            scalar_path = joinpath(dir, "scalar_step1.jld2")
            function save_at_one(sim)
                sim.model.clock.iteration == 1 && R.save_dynamic_restart_state(sim, scalar_path)
                return nothing
            end
            O.add_callback!(scalar_run, save_at_one, O.IterationInterval(1); name = :scalar_restart_save)
            scalar_run.stop_iteration = 1
            O.run!(scalar_run)
            @test isfile(scalar_path)
            before = prognostic(scalar_run.model.atmosphere.variables)
            @test_throws ArgumentError R.restore_dynamic_restart_state!(scalar_run, path)
            @test scalar_run.model.clock.iteration == 1
            @test all(same_bits(x, y) for (x, y) in zip(prognostic(scalar_run.model.atmosphere.variables), before))
            # Off to on: the momentum configuration refuses the scalar-only checkpoint, and keeps its state.
            fresh = R.build_dynamic_esm(compact_l27(joinpath(dir, "fresh"); settings...); steps = 4)
            initial = prognostic(fresh.model.atmosphere.variables)
            @test_throws ArgumentError R.restore_dynamic_restart_state!(fresh, scalar_path)
            @test fresh.model.clock.iteration == 0
            @test all(same_bits(x, y) for (x, y) in zip(prognostic(fresh.model.atmosphere.variables), initial))
        end
    end
end
end
println("PBL_K_PROFILE_MOMENTUM_UNIT_PASS")
