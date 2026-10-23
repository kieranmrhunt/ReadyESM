module PBLKProfileTests
using ReadyESM, Test, Random
const R = ReadyESM

const G = 9.81
const RD = 287.0
const THETA = 300.0

# A 27-layer sigma column like production's (top to bottom), with layer thickness and half levels.
sigma_half = collect(range(0.0, 1.0; length = 28))
sigma_full = [(sigma_half[k] + sigma_half[k + 1]) / 2 for k in 1:27]
dsigma = diff(sigma_half)

@testset "PBL K-profile column functions" begin
    rng = MersenneTwister(20261023)
    n = length(sigma_full)
    c, d, out = zeros(n), zeros(n), zeros(n)

    @testset "conservation, monotonicity and positivity" begin
        worst = 0.0
        for trial in 1:200
            x = 3.0e5 .+ 2.0e4 .* rand(rng, n)
            q = 0.02 .* rand(rng, n) .^ 3
            a = 10.0 .^ (4 .* rand(rng, n - 1) .- 6)       # 1e-6 to 1e-2 1/s
            dt = 1800.0
            R.pbl_implicit_diffuse!(out, x, a, dsigma, dt, c, d)
            change = sum(dsigma .* (out .- x))
            worst = max(worst, abs(change) / sum(dsigma .* abs.(x)))
            @test minimum(x) - 1e-9 <= minimum(out) && maximum(out) <= maximum(x) + 1e-9
            qn = R.pbl_implicit_diffuse!(zeros(n), q, a, dsigma, dt, c, d)
            @test all(qn .>= 0)
            @test abs(sum(dsigma .* (qn .- q))) <= 1e-14 * sum(dsigma .* q) + 1e-20
        end
        println("PBL_CONSERVATION_WORST_RELATIVE ", worst)
        @test worst < 1e-14
        # Float32, as in the model.
        x32 = Float32.(3.0e5 .+ 2.0e4 .* rand(rng, n))
        a32 = Float32.(10.0 .^ (4 .* rand(rng, n - 1) .- 6))
        o32 = R.pbl_implicit_diffuse!(zeros(Float32, n), x32, a32, Float32.(dsigma), 1800.0f0, zeros(Float32, n),
                                      zeros(Float32, n))
        @test abs(sum(Float64.(dsigma) .* (Float64.(o32) .- Float64.(x32)))) / sum(Float64.(dsigma) .* Float64.(x32)) < 1e-6
    end

    @testset "analytic decay of discrete Neumann modes" begin
        m_layers = 20
        ds = fill(1 / m_layers, m_layers)
        aval = 3.0e-4
        a = fill(aval, m_layers - 1)
        dt = 1800.0
        cc, dd, oo = zeros(m_layers), zeros(m_layers), zeros(m_layers)
        for mode in 1:5
            x = [cos(pi * mode * (k - 0.5) / m_layers) for k in 1:m_layers]
            lambda = aval / ds[1] * 2 * (1 - cos(pi * mode / m_layers))
            R.pbl_implicit_diffuse!(oo, x, a, ds, dt, cc, dd)
            @test maximum(abs.(oo .- x ./ (1 + dt * lambda))) < 1e-12
        end
        # The constant mode is untouched.
        @test maximum(abs.(R.pbl_implicit_diffuse!(oo, ones(m_layers), a, ds, dt, cc, dd) .- 1)) < 1e-14
    end

    @testset "stress: coefficient a = 1e3 s^-1 stays finite and mixes toward the column mean" begin
        x = 3.0e5 .+ 2.0e4 .* rand(rng, n)
        a = fill(1.0e3, n - 1)
        out2 = R.pbl_implicit_diffuse!(zeros(n), x, a, dsigma, 1800.0, c, d)
        mean = sum(dsigma .* x) / sum(dsigma)
        @test all(isfinite, out2)
        @test maximum(abs.(out2 .- mean)) < 1e-3 * (maximum(x) - minimum(x))
    end

    @testset "calm columns stay finite in all three stability branches (Codex reproducer)" begin
        z3, th3 = [1500.0, 500.0, 100.0], [310.0, 303.0, 300.0]
        for flux in (0.0, -0.02, 0.1)
            h3 = R.pbl_height(z3, th3, zeros(3), 0.0, flux, 300.0, G, 100.0, 3000.0)
            @test isfinite(h3) && 100.0 <= h3 <= 3000.0
            @test isfinite(R.pbl_velocity_scale(50.0, h3, 0.0, flux, 300.0, G))
            @test isfinite(R.pbl_diffusivity(50.0, h3, 0.0, flux, 300.0, G))
        end
        @test R.pbl_velocity_scale(50.0, 500.0, 0.0, 0.0, 300.0, G) == 0
        @test R.pbl_velocity_scale(50.0, 500.0, 0.0, -0.02, 300.0, G) == 0
        @test R.pbl_velocity_scale(50.0, 500.0, 0.0, 0.1, 300.0, G) ≈ cbrt(0.6 * (G / 300 * 0.1 * 500))
        # Stable calm: the layer stays at its floor; unstable calm: the thermal excess deepens it.
        @test 100.0 <= R.pbl_height(z3, th3, zeros(3), 0.0, -0.02, 300.0, G, 100.0, 3000.0) < 101.0
        @test R.pbl_height(z3, th3, zeros(3), 0.0, 0.1, 300.0, G, 100.0, 3000.0) > 100.0
    end

    @testset "partially intersected top cell (coefficients from the closure)" begin
        # Codex's case: centres at 1500/500/100 m, faces at 1000/300 m, h = 400 m. Only the 300 m face is active,
        # so the cell centred at 500 m changes through its lower face while the top cell is unchanged.
        z_half4 = [3000.0, 1000.0, 300.0, 0.0]
        a4 = zeros(2)
        R.pbl_interface_coefficients!(a4, z_half4, [0.75, 0.85, 0.94], [0.7, 0.8, 0.9, 1.0], fill(300.0, 4), 400.0,
                                      0.2, 0.01, 300.0, G, RD)
        @test a4[1] == 0 && a4[2] > 0
        x4 = [320.0, 310.0, 300.0]
        o4 = R.pbl_implicit_diffuse!(zeros(3), x4, a4, [0.1, 0.1, 0.1], 1800.0, zeros(3), zeros(3))
        @test o4[1] == x4[1] && o4[2] != x4[2] && o4[3] != x4[3]
    end

    @testset "flux-form tendency: conservative and equal to the state difference" begin
        x = 3.0e5 .+ 2.0e4 .* rand(rng, n)
        a = 10.0 .^ (4 .* rand(rng, n - 1) .- 6)
        tend, xn = zeros(n), zeros(n)
        R.pbl_implicit_tendency!(tend, x, a, dsigma, 1800.0, c, d, xn)
        @test maximum(abs.(tend .- (xn .- x) ./ 1800.0)) < 1e-9 * maximum(abs.(tend))
        @test abs(sum(dsigma .* tend)) < 1e-13 * sum(dsigma .* abs.(tend))
        # Float32, 1000 steps, physical K = 1000 m2/s through the sigma metric on the lowest ten faces.
        z_half = vcat([-8000 * log(max(s, 1e-6)) for s in sigma_half[1:end-1]], 0.0)
        a_phys = zeros(n - 1)
        for k in (n - 10):(n - 1)
            a_phys[k] = (sigma_half[k + 1] * G / (RD * THETA))^2 * 1000.0 / (sigma_full[k + 1] - sigma_full[k])
        end
        drift(tendency_form) = begin
            x32 = Float32.(3.0e5 .+ 2.0e4 .* rand(MersenneTwister(7), n))
            initial = sum(Float64.(dsigma) .* Float64.(x32))
            a32, ds32 = Float32.(a_phys), Float32.(dsigma)
            t32, n32, c32, d32 = zeros(Float32, n), zeros(Float32, n), zeros(Float32, n), zeros(Float32, n)
            for step in 1:1000
                if tendency_form
                    R.pbl_implicit_tendency!(t32, x32, a32, ds32, 1800.0f0, c32, d32, n32)
                else
                    R.pbl_implicit_diffuse!(n32, x32, a32, ds32, 1800.0f0, c32, d32)
                    t32 .= (n32 .- x32) ./ 1800.0f0
                end
                x32 .+= 1800.0f0 .* t32
            end
            (sum(Float64.(dsigma) .* Float64.(x32)) - initial) / initial
        end
        state_drift, flux_drift = drift(false), drift(true)
        println("PBL_FLOAT32_1000STEP_DRIFT state_difference=$state_drift flux_form=$flux_drift")
        @test abs(flux_drift) <= abs(state_drift)
        @test abs(flux_drift) < 1e-5
    end

    @testset "no tendency above the boundary layer" begin
        x = 3.0e5 .+ 2.0e4 .* rand(rng, n)
        a = zeros(n - 1)
        a[end - 5:end] .= 1.0e-3                  # mixing only across the lowest six interfaces
        out3 = R.pbl_implicit_diffuse!(zeros(n), x, a, dsigma, 1800.0, c, d)
        @test out3[1:(n - 7)] == x[1:(n - 7)]
        @test out3[end] != x[end]
    end

    @testset "velocity scale and diffusivity" begin
        ustar = 0.3
        @test R.pbl_obukhov_length(ustar, 0.0, THETA, G) == Inf
        @test R.pbl_velocity_scale(100.0, 1000.0, ustar, 0.0, THETA, G) == ustar
        stable_flux = -0.02
        L = R.pbl_obukhov_length(ustar, stable_flux, THETA, G)
        @test L ≈ ustar^3 * THETA / (0.4 * G * 0.02)
        @test R.pbl_velocity_scale(50.0, 200.0, ustar, stable_flux, THETA, G) ≈ ustar / (1 + 5 * 50 / L)
        unstable_flux = 0.1
        h = 1200.0
        wstar = cbrt(G / THETA * unstable_flux * h)
        @test R.pbl_convective_velocity(h, unstable_flux, THETA, G) ≈ wstar
        @test R.pbl_velocity_scale(300.0, h, ustar, unstable_flux, THETA, G) ≈ cbrt(ustar^3 + 0.6 * wstar^3)
        @test R.pbl_diffusivity(300.0, h, ustar, unstable_flux, THETA, G) ≈
              0.4 * cbrt(ustar^3 + 0.6 * wstar^3) * 300 * (1 - 300 / h)^2
        @test R.pbl_diffusivity(h, h, ustar, unstable_flux, THETA, G) == 0
        @test R.pbl_diffusivity(1500.0, h, ustar, unstable_flux, THETA, G) == 0
        @test R.pbl_diffusivity(0.0, h, ustar, unstable_flux, THETA, G) == 0
    end

    @testset "boundary-layer depth" begin
        # Heights of full levels for a 27-layer column with an 8 km scale height, and a stable layer above 1 km.
        z_full = [-8000 * log(s) for s in sigma_full]
        thetav = [THETA + (z > 1000 ? 0.008 * (z - 1000) : 0.0) for z in z_full]
        wind2 = fill(25.0, n)
        ustar = 0.3
        h = R.pbl_height(z_full, thetav, wind2, ustar, 0.0, THETA, G, z_full[end], 5000.0)
        # Independent linear interpolation of the same Ri_b profile.
        ri = [G * z_full[k] * (thetav[k] - thetav[end]) /
              (THETA * (wind2[k] + 100 * ustar^2 + R.PBL_CALM_WIND_FLOOR^2)) for k in 1:n]
        k = findlast(r -> r >= 0.3, ri)
        expected = z_full[k + 1] + (0.3 - ri[k + 1]) / (ri[k] - ri[k + 1]) * (z_full[k] - z_full[k + 1])
        @test h ≈ expected
        # Linear interpolation of Ri_b between the bracketing full levels: h lies between them (with a sharp
        # inversion and coarse levels it can sit below the inversion base).
        @test z_full[k + 1] <= h <= z_full[k]
        # An upward buoyancy flux deepens the layer (thermal excess); the clamps hold.
        @test R.pbl_height(z_full, thetav, wind2, ustar, 0.15, THETA, G, z_full[end], 5000.0) > h
        @test R.pbl_height(z_full, thetav, wind2, ustar, 0.0, THETA, G, 2500.0, 5000.0) == 2500.0
        @test R.pbl_height(z_full, fill(THETA, n), wind2, ustar, 0.0, THETA, G, z_full[end], 3000.0) == 3000.0
        # Interface coefficients: zero at and above h, positive below, and dimensionally the hydrostatic metric.
        z_half = vcat([-8000 * log(max(s, 1e-6)) for s in sigma_half[1:end-1]], 0.0)
        tv_half = fill(THETA, n + 1)
        a = zeros(n - 1)
        R.pbl_interface_coefficients!(a, z_half, sigma_full, sigma_half, tv_half, h, ustar, 0.0, THETA, G, RD)
        for k in 1:(n - 1)
            z = z_half[k + 1]
            if z >= h
                @test a[k] == 0
            else
                K = R.pbl_diffusivity(z, h, ustar, 0.0, THETA, G)
                @test a[k] ≈ (sigma_half[k + 1] * G / (RD * THETA))^2 * K / (sigma_full[k + 1] - sigma_full[k])
                @test a[k] > 0
            end
        end
    end
end
end
println("PBL_K_PROFILE_UNIT_PASS")
