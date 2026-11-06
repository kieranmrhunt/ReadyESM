module RRTMGPBandOpticsHookTests
using ReadyESM, Test
const R = ReadyESM
const G = R.RRTMGP
const S = R.SpeedyWeather
const OP = G.Optics

# Hand-written two-stream combination and f = g^2 delta scaling, independent of the RRTMGP source.
combine(τ1, ω1, g1, τ2, ω2, g2) = (τ1 + τ2, (τ1 * ω1 + τ2 * ω2) / (τ1 + τ2),
                                   (τ1 * ω1 * g1 + τ2 * ω2 * g2) / (τ1 * ω1 + τ2 * ω2))
delta(τ, ω, g) = (τ * (1 - ω * g^2), ω * (1 - g^2) / (1 - ω * g^2), g / (1 + g))

"""Band optics for one column with `nbnd` bands and `nlay` layers, zero except where `set!` writes."""
function band_optics(nbnd_sw, nbnd_lw, nlay; T = Float64)
    z(n) = zeros(T, n, nlay, 1)
    return (tau_sw = z(nbnd_sw), ssa_sw = z(nbnd_sw), g_sw = z(nbnd_sw),
            tau_lw = z(nbnd_lw), ssa_lw = z(nbnd_lw), g_lw = z(nbnd_lw))
end

@testset "RRTMGP band-optics hook (vendored RRTMGP 0.22.2 ReadyESM patch)" begin
    @testset "analytic one-column values" begin
        bo = band_optics(14, 16, 3)
        bo.tau_lw[6, 2, 1], bo.ssa_lw[6, 2, 1], bo.g_lw[6, 2, 1] = 0.3, 0.4, 0.5
        bo.tau_sw[11, 1, 1], bo.ssa_sw[11, 1, 1], bo.g_sw[11, 1, 1] = 0.2, 0.9, 0.7
        # LW absorption only: τ += τ_b (1 - ssa_b) in the chosen band and layer, nothing elsewhere.
        τ = [1.0, 2.0, 3.0]
        OP.add_band_optics_lw_1scalar!(τ, bo, 6, 1)
        @test τ == [1.0, 2.0 + 0.3 * (1 - 0.4), 3.0]
        τ = [1.0, 2.0, 3.0]
        OP.add_band_optics_lw_1scalar!(τ, bo, 5, 1)
        @test τ == [1.0, 2.0, 3.0]
        OP.add_band_optics_lw_1scalar!(τ, nothing, 6, 1)
        @test τ == [1.0, 2.0, 3.0]
        # LW two-stream: the τ-weighted combination, no delta scaling.
        τ, ω, g = [1.0, 2.0, 3.0], [0.1, 0.2, 0.3], [0.0, 0.1, 0.2]
        OP.add_band_optics_lw_2stream!(τ, ω, g, bo, 6, 1)
        @test all(isapprox.((τ[2], ω[2], g[2]), combine(2.0, 0.2, 0.1, 0.3, 0.4, 0.5); rtol = 1e-14))
        @test (τ[[1, 3]], ω[[1, 3]], g[[1, 3]]) == ([1.0, 3.0], [0.1, 0.3], [0.0, 0.2])
        # SW two-stream: delta-scaled before combining; the 550 nm diagnostics get the unscaled values.
        τ, ω, g = [1.0, 2.0, 3.0], [0.5, 0.6, 0.7], [0.1, 0.2, 0.3]
        aod_ext, aod_sca = fill(0.05), fill(0.04)
        OP.add_band_optics_sw_2stream!(τ, ω, g, aod_ext, aod_sca, bo, 11, 1, 11)
        @test all(isapprox.((τ[1], ω[1], g[1]), combine(1.0, 0.5, 0.1, delta(0.2, 0.9, 0.7)...); rtol = 1e-14))
        @test aod_ext[] ≈ 0.05 + 0.2 && aod_sca[] ≈ 0.04 + 0.2 * 0.9
        aod_ext, aod_sca = fill(0.05), fill(0.04)
        OP.add_band_optics_sw_2stream!([1.0, 2.0, 3.0], [0.5, 0.6, 0.7], [0.1, 0.2, 0.3], aod_ext, aod_sca, bo, 11, 1, 4)
        @test aod_ext[] == 0.05 && aod_sca[] == 0.04          # not the 550 nm band: diagnostics untouched
        τ, ω, g = [1.0, 2.0, 3.0], [0.5, 0.6, 0.7], [0.1, 0.2, 0.3]
        OP.add_band_optics_sw_2stream!(τ, ω, g, nothing, nothing, nothing, 11, 1, 11)
        @test (τ, ω, g) == ([1.0, 2.0, 3.0], [0.5, 0.6, 0.7], [0.1, 0.2, 0.3])
        # Zero prescribed optics leave every value bit for bit unchanged.
        z = band_optics(14, 16, 3)
        τ, ω, g = [1.0, 2.0, 3.0], [0.5, 0.6, 0.7], [0.1, 0.2, 0.3]
        OP.add_band_optics_sw_2stream!(τ, ω, g, fill(0.0), fill(0.0), z, 11, 1, 11)
        OP.add_band_optics_lw_2stream!(τ, ω, g, z, 6, 1)
        OP.add_band_optics_lw_1scalar!(τ, z, 6, 1)
        @test (τ, ω, g) == ([1.0, 2.0, 3.0], [0.5, 0.6, 0.7], [0.1, 0.2, 0.3])
    end

    @testset "solver: zero prescribed optics are bitwise neutral, nonzero ones act" begin
        grid = S.SpectralGrid(; NF = Float32, trunc = 7, nlayers = 8)
        base = R.RRTMGPRadiation(grid).solver
        plain = R._enable_sulfate_aerosol(base, 0.0f0)
        hooked = R._enable_sulfate_aerosol(base, 0.0f0; band_optics = true)
        @test isnothing(plain.as.aerosol_state.band_optics)
        bo = hooked.as.aerosol_state.band_optics
        @test size(bo.tau_sw, 1) == hooked.lookups.nbnd_sw == 14 && size(bo.tau_lw, 1) == hooked.lookups.nbnd_lw == 16
        fluxes(s) = [Array(f(s)) for f in (G.sw_flux_up, G.sw_flux_dn, G.lw_flux_up, G.lw_flux_dn, G.net_flux)]
        R._reset_solver_cloud_sampler!(plain, 0); G.update_fluxes!(plain)
        R._reset_solver_cloud_sampler!(hooked, 0); G.update_fluxes!(hooked)
        @test fluxes(plain) == fluxes(hooked)                                      # bitwise
        @test Array(G.aod_sw_extinction(plain)) == Array(G.aod_sw_extinction(hooked))
        # A pure scatterer of optical depth 0.1 in the 550 nm band, spread over the layers: the AOD diagnostic rises
        # by exactly that amount and the SW fluxes change; the LW fluxes do not.
        iband = hooked.lookups.lookup_sw_aero.iband_550nm
        nlay = size(bo.tau_sw, 2)
        bo.tau_sw[iband, :, :] .= 0.1f0 / nlay; bo.ssa_sw[iband, :, :] .= 1; bo.g_sw[iband, :, :] .= 0.7f0
        before = Array(G.aod_sw_extinction(hooked))
        G.update_fluxes!(hooked)
        after = Array(G.aod_sw_extinction(hooked))
        @test maximum(abs.(after .- before .- 0.1f0)) < 1.0f-5
        f0, f1 = fluxes(plain), fluxes(hooked)
        @test f0[1] != f1[1] && f0[2] != f1[2] && f0[3] == f1[3] && f0[4] == f1[4]
        # An LW absorber changes the LW fluxes.
        bo.tau_lw[6, :, :] .= 0.05f0; bo.ssa_lw[6, :, :] .= 0
        G.update_fluxes!(hooked)
        f2 = fluxes(hooked)
        @test f2[3] != f1[3] && all(isfinite, f2[3]) && all(isfinite, f2[2])
    end
end
end
println("RRTMGP_BAND_OPTICS_HOOK_PASS")
