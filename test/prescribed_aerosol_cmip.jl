module PrescribedAerosolCMIPTests
using ReadyESM, Test, Dates
const R = ReadyESM
const G = R.RRTMGP
const S = R.SpeedyWeather

# The prepared MACv2 background is development-only data (not redistributed), so the test reads the pinned files from
# the paths the CPU packet sets; a missing path is a failure, not a skip.
const BACKGROUND = get(ENV, "READYESM_PA2_BACKGROUND", "")
const PLUMES = get(ENV, "READYESM_PA2_PLUMES", "")

fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
with_forcing(config; kwargs...) = R.ExperimentConfig(; merge(fields(config),
    (; forcing = R.ForcingConfig(; merge(fields(config.forcing), (; kwargs...))...)))...)
slab(; kwargs...) = R.ExperimentConfig(; name = "pa2", device = :cpu, truncation = 7, nlayers = 8,
    duration_days = 0.5, mixed_layer_depth_m = 30.0, start_date = DateTime(1993, 1, 1),
    forcing = R.ForcingConfig(; radiation = :rrtmgp_clear_sky, kwargs...))

@testset "Prescribed CMIP aerosol (PA2: MACv2 background + MACv2-SP)" begin
    @test isfile(BACKGROUND) && isfile(PLUMES)

    @testset "helpers" begin
        # Monthly midpoints at (m - 0.5)/12, cyclic.
        @test R._month_weights(0.125) == (2, 3, 0.0)            # 12 x 0.125 - 0.5 = 1, exactly
        m0, m1, w = R._month_weights(0.0)
        @test (m0, m1) == (12, 1) && w ≈ 0.5
        # Periodic bilinear weights on a -179.5..179.5 grid.
        glon, glat = collect(-179.5:1.0:179.5), collect(-89.5:1.0:89.5)
        w1 = R._bilinear_weights(glon, glat, 180.0, 0.0)
        @test w1[1] == 360 && w1[2] == 1 && w1[5] ≈ 0.5
        w2 = R._bilinear_weights(glon, glat, 10.25, 89.9)
        @test w2[6] ≈ 1.0 && w2[3] == 179
        # Optics combination and the ETH tropopause-layer rule.
        @test all(R._combine_optics(((0.2, 0.9, 0.7), (0.1, 0.5, 0.4))) .≈
                  (0.3, (0.18 + 0.05) / 0.3, (0.2 * 0.9 * 0.7 + 0.1 * 0.5 * 0.4) / (0.18 + 0.05)))
        @test R._combine_optics(((0.0, 1.0, 0.0),)) == (0.0, 1.0, 0.0)
        @test all(R._tropopause_mean((0.2, 0.9, 0.6), (0.04, 0.99, 0.7)) .≈ (0.12, 0.945, 0.65))
        @test R._tropopause_mean((0.2, 0.9, 0.6), (0.0, 1.0, 0.0)) == (0.1, 0.9, 0.6)
        # WMO lapse-rate tropopause on a synthetic column: 6.5 K/km to 11 km, isothermal above.
        z = reshape(collect(500.0:1000.0:19500.0), :, 1)
        t = reshape([zz <= 11000 ? 288 - 6.5e-3 * zz : 288 - 6.5e-3 * 11000 for zz in vec(z)], :, 1)
        p = reshape([101325 * exp(-zz / 7500) for zz in vec(z)], :, 1)
        tp = R._tropopause_layers(t, p, z)
        # 10500 -> 11500 m still cools at 3.25 K/km; the first layer with <= 2 K/km held over 2 km is at 11500 m.
        @test z[tp[1], 1] ≈ 11500.0
    end

    @testset "config validation" begin
        on = slab(; aerosol_scheme = :cmip_prescribed, aerosol_cmip_background_path = BACKGROUND,
                  aerosol_cmip_plumes_path = PLUMES)
        @test R.validate(on; check_input_files = false) === on
        @test R._validate_cmip_aerosol_options(on, true, true) === nothing
        @test_throws ArgumentError R.validate(with_forcing(on; aerosol_optical_depth_550nm = 0.1); check_input_files = false)
        @test_throws ArgumentError R.validate(with_forcing(on; radiation = :speedy_simplified); check_input_files = false)
        @test_throws ArgumentError R.validate(with_forcing(on; aerosol_cmip_plumes_path = ""); check_input_files = false)
        @test_throws ArgumentError R.validate(with_forcing(on; aerosol_cmip_year = 1700); check_input_files = false)
        @test_throws ArgumentError R.validate(slab(; aerosol_cmip_year = 1993); check_input_files = false)
        @test_throws ArgumentError R.validate(slab(; aerosol_scheme = :macv2); check_input_files = false)
        production = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)
        @test production.forcing.aerosol_scheme == :none
    end

    @testset "model: off path has no scheme; on path writes validated optics" begin
        off = R.build_rrtmgp(slab())
        S.run!(off; steps = 2)
        @test isnothing(off.model.longwave_radiation.cmip_aerosol)
        on = slab(; aerosol_scheme = :cmip_prescribed, aerosol_cmip_background_path = BACKGROUND,
                  aerosol_cmip_plumes_path = PLUMES)
        simulation = R.build_rrtmgp(on)
        S.run!(simulation; steps = 3)
        radiation = simulation.model.longwave_radiation
        a = radiation.cmip_aerosol
        @test a isa R.CMIPAerosol && a.writes >= 1
        solver = radiation.solver
        @test all(iszero, solver.as.aerosol_state.aero_mass)
        bo = solver.as.aerosol_state.band_optics
        @test R._validate_band_optics(solver) === nothing
        @test maximum(bo.tau_sw) > 0 && maximum(bo.tau_lw) > 0
        aod = Float64.(Array(G.aod_sw_extinction(solver)))
        @test all(isfinite, aod) && all(>=(0), aod)
        weights = Float64.(Array(R._global_point_weights(simulation.model.spectral_grid)))
        gm = sum(weights .* aod) / sum(weights)
        println("PA2_MODEL_GLOBAL_MEAN_AOD_550_BAND ", gm)
        @test 0.03 < gm < 0.5
        nlay = size(bo.tau_sw, 2)
        @test all(1 .<= a.tropopause_layer .<= nlay)
        @test all(isfinite, Array(S.get_step(simulation.variables.prognostic.temperature, 1).data))
        provenance = R._forcing_provenance(radiation, on)
        @test provenance.aerosol_forcing_active && provenance.cmip_aerosol == R.CMIP_AEROSOL_LABEL_TROPOSPHERIC
        @test provenance.cmip_aerosol_validated_writes == a.writes
        # Public-path admission: invalid source optics are refused before anything is copied into the solver.
        before = Array(bo.tau_sw)
        old = a.background.tau_fine_sw[1, 1, 1, 1]
        a.background.tau_fine_sw[:, :, :, 11] .= NaN32
        @test_throws ArgumentError R._apply_cmip_aerosol!(a, solver, simulation.variables, simulation.model)
        @test Array(bo.tau_sw) == before
        a.background.tau_fine_sw[:, :, :, 11] .= R.MACv2Background(BACKGROUND).tau_fine_sw[:, :, :, 11]
        writes = a.writes
        @test R._apply_cmip_aerosol!(a, solver, simulation.variables, simulation.model) === nothing
        @test a.writes == writes + 1
        # A nonzero MERRA mass alongside the prescribed optics is refused through the same path.
        mass = G.aerosol_column_mass_density(solver, "sulfate")
        mass[1, 1] = 1.0f-6
        @test_throws ArgumentError R._apply_cmip_aerosol!(a, solver, simulation.variables, simulation.model)
        mass[1, 1] = 0
        # Restart identity: a checkpoint with the scheme does not restore into radiation without it, and back.
        state = R._rrtmgp_radiation_state(radiation)
        @test state.cmip_aerosol == (background_sha256 = a.background.sha256, plumes_sha256 = a.plumes.sha256,
                                     fixed_year = 0)
        plain = off.model.longwave_radiation
        @test_throws ArgumentError R._restore_rrtmgp_radiation_state!(plain, state)
        @test_throws ArgumentError R._restore_rrtmgp_radiation_state!(radiation, R._rrtmgp_radiation_state(plain))
    end
end
end
println("PRESCRIBED_AEROSOL_CMIP_PASS")
