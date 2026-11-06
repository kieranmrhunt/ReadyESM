module PrescribedAerosolStratTests
using ReadyESM, Test, Dates
const R = ReadyESM
const S = R.SpeedyWeather
const NC = R.NCDatasets

# Component C of PA2: the prepared CMIP7 stratospheric optics (UOEXETER-CMIP-2-2-1). The prepared files are
# development-only data (not redistributed), so the test reads the pinned files from the paths the CPU packet sets; a
# missing path is a failure, not a skip. Synthetic files written here give the exact answers.
const BACKGROUND = get(ENV, "READYESM_PA2_BACKGROUND", "")
const PLUMES = get(ENV, "READYESM_PA2_PLUMES", "")
const STRAT = get(ENV, "READYESM_PA2_STRAT", "")

fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
with_forcing(config; kwargs...) = R.ExperimentConfig(; merge(fields(config),
    (; forcing = R.ForcingConfig(; merge(fields(config.forcing), (; kwargs...))...)))...)
slab(; kwargs...) = R.ExperimentConfig(; name = "pa2c", device = :cpu, truncation = 7, nlayers = 8,
    duration_days = 0.5, mixed_layer_depth_m = 30.0, start_date = DateTime(1993, 1, 1),
    forcing = R.ForcingConfig(; radiation = :rrtmgp_clear_sky, kwargs...))

days1850(d) = Dates.value(d - DateTime(1850, 1, 1)) / 86_400_000
const MONTHS = [DateTime(1992, 12, 1), DateTime(1993, 1, 1), DateTime(1993, 2, 1), DateTime(1993, 3, 1)]
const E0 = 1.0e-6

"""
Synthetic prepared file: 14 SW and 16 LW bands, 70 layers of 500 m from 5 km, latitudes -60, 0, 60, months December
1992 to February 1993. Extinction E0 x (latitude index) x (1, 2, 4 by month) between 15 and 25 km, zero elsewhere;
ssa 0.9 and g 0.6. `order` permutes the moment variables' dimensions in the file.
"""
function write_strat(path; order = (1, 2, 3, 4), schema = R.CMIP_AEROSOL_STRAT_SCHEMA, negative = false,
                     g_case = :none)
    nsw, nlw, nh, nlat, nt = 14, 16, 70, 3, 3
    bottom = collect(5000.0:500.0:39500.0)
    NC.NCDataset(path, "c") do ds
        ds.attrib["schema"] = schema
        for (n, l) in (("time", nt), ("lat", nlat), ("height", nh), ("sw_band", nsw), ("lw_band", nlw), ("bnds", 2))
            NC.defDim(ds, n, l)
        end
        tb = NC.defVar(ds, "time_bnds", Float64, ("bnds", "time"))
        tb.var[:, :] = [days1850(MONTHS[i + k - 1]) for k in 1:2, i in 1:nt]
        tb.attrib["units"] = "days since 1850-01-01"
        tb.attrib["calendar"] = "proleptic_gregorian"
        NC.defVar(ds, "lat", [-60.0, 0.0, 60.0], ("lat",))
        NC.defVar(ds, "height_bnds", permutedims(hcat(bottom, bottom .+ 500)), ("bnds", "height"))
        for (spec, nb) in (("sw", nsw), ("lw", nlw))
            E = zeros(nb, nh, nlat, nt)
            for t in 1:nt, j in 1:nlat, h in 1:nh, b in 1:nb
                15000 <= bottom[h] < 25000 && (E[b, h, j, t] = E0 * j * (1, 2, 4)[t])
            end
            negative && (E[1, 25, 1, 1] = -E0)
            names = ("$(spec)_band", "height", "lat", "time")
            Sm, Gm = 0.9 .* E, 0.6 .* 0.9 .* E
            # Emitted moments the reader must refuse: negative G, G = S (g = 1), G > 0 where S = 0.
            g_case == :negative && (Gm[1, 25, 2, 2] = -1.0e-9)
            g_case == :equal && (Gm[1, 25, 2, 2] = Sm[1, 25, 2, 2])
            g_case == :orphan && (Gm[1, 1, 2, 2] = 1.0e-9)
            for (m, A) in (("E", E), ("S", Sm), ("G", Gm))
                NC.defVar(ds, "$(m)_$spec", permutedims(A, order), names[collect(order)])
            end
        end
    end
    return path
end

@testset "Prescribed CMIP aerosol component C (stratospheric, UOEXETER-CMIP-2-2-1)" begin
    @test isfile(BACKGROUND) && isfile(PLUMES) && isfile(STRAT)
    dir = mktempdir()
    syn = write_strat(joinpath(dir, "strat.nc"))

    @testset "reader: axes by dimension name; refusals" begin
        st = R.StratosphericAerosol(syn)
        perm = R.StratosphericAerosol(write_strat(joinpath(dir, "strat_perm.nc"); order = (2, 4, 1, 3)))
        @test size(st.E_sw) == (14, 70, 3, 3) && size(st.E_lw) == (16, 70, 3, 3)
        @test st.E_sw == perm.E_sw && st.G_lw == perm.G_lw
        @test st.month_mid[2] == DateTime(1993, 1, 16, 12)
        @test st.height_bottom_m[1] == 5000 && st.height_top_m[end] == 40000
        @test_throws ArgumentError R.StratosphericAerosol(write_strat(joinpath(dir, "bad_schema.nc"); schema = "x"))
        @test_throws ArgumentError R.StratosphericAerosol(write_strat(joinpath(dir, "negative.nc"); negative = true))
        for g_case in (:negative, :equal, :orphan)
            @test_throws ArgumentError R.StratosphericAerosol(write_strat(joinpath(dir, "g_$g_case.nc"); g_case))
        end
        @test_throws ArgumentError R.StratosphericAerosol(joinpath(dir, "missing.nc"))
    end

    @testset "time and latitude weights" begin
        st = R.StratosphericAerosol(syn)
        @test R._strat_time_weights(st, DateTime(1993, 1, 16, 12)) == (2, 3, 0.0)
        @test R._strat_time_weights(st, st.month_mid[end]) == (2, 3, 1.0)
        i0, i1, w = R._strat_time_weights(st, DateTime(1993, 1, 1))
        @test (i0, i1) == (1, 2) && 0 < w < 1
        @test_throws ArgumentError R._strat_time_weights(st, DateTime(1992, 12, 1))
        @test_throws ArgumentError R._strat_time_weights(st, DateTime(1993, 3, 1))
        @test R._strat_lat_weights(st.lat, 30.0) == (2, 3, 0.5)
        @test R._strat_lat_weights(st.lat, 80.0) == (2, 3, 1.0)
        @test R._strat_lat_weights(st.lat, -90.0) == (1, 2, 0.0)
    end

    @testset "column moments are conservative in height" begin
        st = R.StratosphericAerosol(syn)
        tw = R._strat_time_weights(st, DateTime(1993, 1, 16, 12))          # January: month factor 2
        lw = R._strat_lat_weights(st.lat, 30.0)                            # latitude factor 2.5
        zbc = [0.0, 14750.0, 15250.0, 20000.0, 40000.0, 50000.0]
        M = zeros(3, 14, 5)
        R._strat_column_moments!(M, st.E_sw, st.S_sw, st.G_sw, st, lw, tw, zbc)
        k = E0 * 2.5 * 2
        @test M[1, 1, 1] == 0 && M[1, 1, 5] == 0
        @test M[1, 1, 2] ≈ k * 250 && M[1, 1, 3] ≈ k * 4750 && M[1, 1, 4] ≈ k * 5000
        @test sum(M[1, 7, :]) ≈ k * 10000
        @test all(M[2, :, 2:4] .≈ 0.9 .* M[1, :, 2:4]) && all(M[3, :, 2:4] .≈ 0.6 .* M[2, :, 2:4])
    end

    @testset "optics: tropospheric below the tropopause, stratospheric above, ETH mean in its layer" begin
        on = slab(; aerosol_scheme = :cmip_prescribed, aerosol_cmip_background_path = BACKGROUND,
                  aerosol_cmip_plumes_path = PLUMES)
        simulation = R.build_rrtmgp(on)
        S.run!(simulation; steps = 1)
        vars, model = simulation.variables, simulation.model
        solver = model.longwave_radiation.solver
        nlay = solver.grid_params.nlay
        zc, zb = R._layer_heights(Float64.(Array(vars.grid.geopotential.data)),
                                  Float64.(vec(Array(model.orography.orography.data))), model.planet.gravity, nlay)
        oro = Float64.(vec(Array(model.orography.orography.data)))
        t_lay = Float64.(Array(R.RRTMGP.layer_temperature(solver)))
        p_lay = Float64.(Array(R.RRTMGP.layer_pressure(solver)))
        size(t_lay, 1) == nlay || (t_lay = vcat(t_lay, repeat(t_lay[end:end, :], nlay - size(t_lay, 1))))
        size(p_lay, 1) == nlay || (p_lay = vcat(p_lay, repeat(p_lay[end:end, :], nlay - size(p_lay, 1))))
        time = DateTime(1993, 1, 16, 12)
        trop = R.CMIPAerosol(BACKGROUND, PLUMES, model.spectral_grid)
        full = R.CMIPAerosol(BACKGROUND, PLUMES, model.spectral_grid; stratospheric_path = syn)
        nsw, nlw = solver.lookups.nbnd_sw, solver.lookups.nbnd_lw
        ot = R._cmip_aerosol_optics(trop, time, zc, zb, oro, t_lay, p_lay, nsw, nlw)
        of = R._cmip_aerosol_optics(full, time, zc, zb, oro, t_lay, p_lay, nsw, nlw)
        @test trop.tropopause_layer == full.tropopause_layer
        st = full.stratospheric
        tw = R._strat_time_weights(st, time)
        M = zeros(3, nsw, nlay)
        below = above = at = true
        strat_column = 0.0
        for c in axes(zc, 2)
            tp = full.tropopause_layer[c]
            R._strat_column_moments!(M, st.E_sw, st.S_sw, st.G_sw, st, full.strat_lat_weights[c], tw, view(zb, :, c))
            for b in 1:nsw
                for r in 1:(tp - 1)
                    below &= of.tau_sw[b, r, c] == ot.tau_sw[b, r, c] && of.ssa_sw[b, r, c] == ot.ssa_sw[b, r, c]
                end
                for r in (tp + 1):nlay
                    above &= ot.tau_sw[b, r, c] == 0 && isapprox(of.tau_sw[b, r, c], M[1, b, r]; rtol = 1e-12) &&
                             (M[1, b, r] == 0 || isapprox(of.ssa_sw[b, r, c], 0.9; rtol = 1e-12))
                end
                # tau at the tropopause layer: (trop + strat) / 2 with the tropospheric-only run holding trop / 2.
                at &= isapprox(of.tau_sw[b, tp, c], ot.tau_sw[b, tp, c] + M[1, b, tp] / 2; rtol = 1e-12, atol = 1e-15)
            end
            strat_column += sum(M[1, 11, :])
        end
        @test below && above && at
        @test strat_column > 0
        # Without the file the optics are those of A + B alone, unchanged from before component C.
        again = R._cmip_aerosol_optics(trop, time, zc, zb, oro, t_lay, p_lay, nsw, nlw)
        @test again.tau_sw == ot.tau_sw && again.tau_lw == ot.tau_lw
        @test_throws ArgumentError R._cmip_aerosol_optics(full, DateTime(1993, 6, 1), zc, zb, oro, t_lay, p_lay,
                                                          nsw, nlw)
        @test_throws ArgumentError R.CMIPAerosol(BACKGROUND, PLUMES, model.spectral_grid; fixed_year = 1993,
                                                 stratospheric_path = syn)
    end

    @testset "config validation" begin
        on = slab(; aerosol_scheme = :cmip_prescribed, aerosol_cmip_background_path = BACKGROUND,
                  aerosol_cmip_plumes_path = PLUMES, aerosol_cmip_stratospheric_path = STRAT)
        @test R.validate(on; check_input_files = false) === on
        @test R._validate_cmip_aerosol_options(on, true, true) === nothing
        @test_throws ArgumentError R.validate(with_forcing(on; aerosol_cmip_year = 1993); check_input_files = false)
        @test_throws ArgumentError R.validate(slab(; aerosol_cmip_stratospheric_path = STRAT);
                                              check_input_files = false)
        late = R.ExperimentConfig(; merge(fields(on), (; start_date = DateTime(1995, 1, 1)))...)
        @test R.validate(late; check_input_files = false) === late
        @test_throws ArgumentError R._validate_cmip_aerosol_options(late, true, true)
        production = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)
        @test production.forcing.aerosol_cmip_stratospheric_path == ""
    end

    @testset "real prepared file: Pinatubo-era column AOD and a model run" begin
        st = R.StratosphericAerosol(STRAT)
        @test size(st.E_sw, 1) == 14 && size(st.E_lw, 1) == 16
        @test st.month_mid[1] < DateTime(1993, 1, 1) && st.month_mid[end] > DateTime(1994, 1, 1)
        # Zonal-mean column optical depth in the 550 nm band (SW band 11), January 1993, area-weighted.
        tw = R._strat_time_weights(st, DateTime(1993, 1, 16, 12))
        thick = st.height_top_m .- st.height_bottom_m
        col = [sum(((1 - tw[3]) * st.E_sw[11, h, j, tw[1]] + tw[3] * st.E_sw[11, h, j, tw[2]]) * thick[h]
                   for h in eachindex(thick)) for j in eachindex(st.lat)]
        gm = sum(col .* cosd.(st.lat)) / sum(cosd.(st.lat))
        println("PA2C_STRAT_GLOBAL_MEAN_AOD_550_BAND_JAN1993 ", gm)
        @test 0.02 < gm < 0.3
        on = slab(; aerosol_scheme = :cmip_prescribed, aerosol_cmip_background_path = BACKGROUND,
                  aerosol_cmip_plumes_path = PLUMES, aerosol_cmip_stratospheric_path = STRAT)
        # The public slab builder starts at the configured date (job 65718246 found it started on 2000-01-01, which the
        # stratospheric file correctly refused).
        simulation = R.build_rrtmgp(on)
        @test simulation.variables.prognostic.clock.time == DateTime(1993, 1, 1)
        S.run!(simulation; steps = 3)
        @test simulation.variables.prognostic.clock.time > DateTime(1993, 1, 1)
        radiation = simulation.model.longwave_radiation
        a = radiation.cmip_aerosol
        @test a.stratospheric isa R.StratosphericAerosol && a.writes >= 1
        @test R._validate_band_optics(radiation.solver) === nothing
        @test all(isfinite, Array(S.get_step(simulation.variables.prognostic.temperature, 1).data))
        provenance = R._forcing_provenance(radiation, on)
        @test provenance.cmip_aerosol == R.CMIP_AEROSOL_LABEL_FULL
        @test provenance.cmip_aerosol_stratospheric_sha256 == a.stratospheric.sha256
        state = R._rrtmgp_radiation_state(radiation)
        @test state.cmip_aerosol.stratospheric_sha256 == a.stratospheric.sha256
        # A configuration outside the file's months is refused at the first radiation application (whether that is in
        # the build or the first step), with no optics written: the clock is at the configured date from the start.
        late = R.ExperimentConfig(; merge(fields(on), (; start_date = DateTime(1995, 1, 1)))...)
        late_refused = try
            sim = R.build_rrtmgp(late)
            @test sim.variables.prognostic.clock.time == DateTime(1995, 1, 1)
            try
                S.run!(sim; steps = 1)
                false
            catch err
                err isa ArgumentError && occursin("outside the stratospheric aerosol months", err.msg) &&
                    sim.model.longwave_radiation.cmip_aerosol.writes == 0
            end
        catch err
            err isa ArgumentError && occursin("outside the stratospheric aerosol months", err.msg)
        end
        @test late_refused
        # A checkpoint with the stratospheric file does not restore into radiation without it.
        trop = R.build_rrtmgp(with_forcing(on; aerosol_cmip_stratospheric_path = ""))
        @test_throws ArgumentError R._restore_rrtmgp_radiation_state!(trop.model.longwave_radiation, state)
    end
end
end
println("PRESCRIBED_AEROSOL_STRAT_PASS")
