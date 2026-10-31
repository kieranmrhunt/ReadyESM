module PrescribedOzoneTests
using ReadyESM, Test, Dates
const R = ReadyESM
const G = R.RRTMGP
const S = R.SpeedyWeather

"""Decimal year of the first instant of (year, month) in a 365-day calendar, as prepare_ozone.py writes noleap files."""
noleap(year, month) = year + (cumsum([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])[month]) / 365

"""
Write a schema-conforming monthly ozone file. `f(p, lat, t)` gives the mixing ratio at pressure p (Pa), latitude and
midpoint decimal year t. Months run from `first` for `nmonths`, in the 365-day calendar.
"""
function write_ozone(path, f; first = (2000, 1), nmonths = 36, pressure = [1.0e5, 5.0e4, 2.0e4, 5.0e3, 1.0e3, 1.0e2, 1.0],
                     latitude = collect(-90.0:30.0:90.0), schema = R.PRESCRIBED_OZONE_SCHEMA, mutate = nothing)
    starts = Tuple{Int, Int}[]
    y, m = first
    for _ in 1:nmonths
        push!(starts, (y, m)); (y, m) = m == 12 ? (y + 1, 1) : (y, m + 1)
    end
    lower = [noleap(s...) for s in starts]
    upper = [m == 12 ? y + 1.0 : noleap(y, m + 1) for (y, m) in starts]
    time = (lower .+ upper) ./ 2
    values = [f(p, φ, t) for φ in latitude, p in pressure, t in time]  # Julia order (lat, plev, time)
    mutate == :negative && (values[1, 1, 1] = -1e-9)
    mutate == :time && (time[3] = time[2])
    R.NCDatasets.NCDataset(path, "c") do ds
        ds.attrib["schema"] = schema
        ds.attrib["family"] = "test"
        ds.attrib["scenario"] = "synthetic"
        R.NCDatasets.defDim(ds, "time", nmonths); R.NCDatasets.defDim(ds, "plev", length(pressure))
        R.NCDatasets.defDim(ds, "lat", length(latitude)); R.NCDatasets.defDim(ds, "bnds", 2)
        R.NCDatasets.defVar(ds, "time", time, ("time",))
        R.NCDatasets.defVar(ds, "time_bnds", permutedims(hcat(lower, upper)), ("bnds", "time"))
        R.NCDatasets.defVar(ds, "year", Int32[s[1] for s in starts], ("time",))
        R.NCDatasets.defVar(ds, "month", Int32[s[2] for s in starts], ("time",))
        R.NCDatasets.defVar(ds, "plev", pressure, ("plev",); attrib = Dict("units" => "Pa"))
        R.NCDatasets.defVar(ds, "lat", latitude, ("lat",))
        R.NCDatasets.defVar(ds, "vmro3", Float32.(values), ("lat", "plev", "time"); attrib = Dict("units" => "mole mole-1"))
    end
    return path, time, lower, upper
end

# Linear in ln p, latitude and time, so every interpolation in the model is exact (up to Float32 storage).
linear_field(p, φ, t) = 5e-6 + 2e-7 * log(p) + 1e-8 * φ + 1e-7 * (t - 2000)
monthly_field(p, φ, t) = 1e-6 * (1 + round(Int, 12 * (t - floor(t)) + 0.5))  # depends on the month only

fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
with_forcing(config; kwargs...) = R.ExperimentConfig(; merge(fields(config),
    (; forcing = R.ForcingConfig(; merge(fields(config.forcing), (; kwargs...))...)))...)

slab_config(path; start = DateTime(2001, 3, 1), kwargs...) = R.ExperimentConfig(; name = "ozone", device = :cpu,
    truncation = 7, nlayers = 8, duration_days = 0.5, mixed_layer_depth_m = 30.0, start_date = start,
    forcing = R.ForcingConfig(; radiation = :rrtmgp_clear_sky, ozone_path = path, kwargs...))

@testset "Prescribed zonal-mean ozone" begin
    dir = mktempdir()
    path, time, lower, upper = write_ozone(joinpath(dir, "linear.nc"), linear_field)
    grid = S.SpectralGrid(; NF = Float32, trunc = 7, nlayers = 8)
    _, model_latitude = R.RG.get_londlatds(grid.grid)

    @testset "reader" begin
        data = R.PrescribedOzoneData(path)
        @test data.time ≈ time && data.lower == lower[1] && data.upper == upper[end]
        @test size(data.values) == (7, 7, 36) && data.family == "test"
        @test data.values[3, 2, 5] ≈ linear_field(2.0e4, -60.0, time[5]) rtol = 1e-6
        @test length(data.sha256) == 64
        @test_throws ArgumentError R.PrescribedOzoneData(write_ozone(joinpath(dir, "s.nc"), linear_field; schema = "x")[1])
        @test_throws ArgumentError R.PrescribedOzoneData(write_ozone(joinpath(dir, "n.nc"), linear_field; mutate = :negative)[1])
        @test_throws ArgumentError R.PrescribedOzoneData(write_ozone(joinpath(dir, "t.nc"), linear_field; mutate = :time)[1])
        @test_throws ArgumentError R.PrescribedOzoneData(joinpath(dir, "missing.nc"))
    end

    @testset "reader: labels, bounds and latitude count (Codex audit 61852576)" begin
        # A 360-day-calendar year as decimal years (the schema can express it), with optional defects.
        function audit_fixture(path, variant)
            lower = [2001 + (m - 1) / 12 for m in 1:12]; upper = [2001 + m / 12 for m in 1:12]
            time = (lower .+ upper) ./ 2
            years = fill(2001, 12); months = collect(1:12)
            variant == :shifted_months && (months = circshift(months, 6))
            variant == :shifted_one_month && (months = circshift(months, 1))
            variant == :shifted_years && (years .= 2000)
            variant == :repeated && (months[5] = 4)
            variant == :infinite_bounds && (lower[1] = -Inf; upper[end] = Inf)
            latitude = variant == :single_latitude ? [0.0] : [-90.0, 90.0]
            R.NCDatasets.NCDataset(path, "c") do ds
                ds.attrib["schema"] = R.PRESCRIBED_OZONE_SCHEMA
                for (name, n) in (("time", 12), ("plev", 2), ("lat", length(latitude)), ("bnds", 2))
                    R.NCDatasets.defDim(ds, name, n)
                end
                R.NCDatasets.defVar(ds, "time", time, ("time",))
                R.NCDatasets.defVar(ds, "time_bnds", permutedims(hcat(lower, upper)), ("bnds", "time"))
                R.NCDatasets.defVar(ds, "year", Int32.(years), ("time",))
                R.NCDatasets.defVar(ds, "month", Int32.(months), ("time",))
                R.NCDatasets.defVar(ds, "plev", [1e5, 1e3], ("plev",); attrib = Dict("units" => "Pa"))
                R.NCDatasets.defVar(ds, "lat", latitude, ("lat",))
                values = [Float64(m) * 1e-6 for φ in latitude, p in 1:2, m in 1:12]
                R.NCDatasets.defVar(ds, "vmro3", values, ("lat", "plev", "time"); attrib = Dict("units" => "mole mole-1"))
            end
            return path
        end
        valid = R.PrescribedOzoneData(audit_fixture(joinpath(dir, "a_valid.nc"), :valid))
        @test R._ozone_time_weights(valid, 2001, DateTime(2001, 1, 1)) == (12, 1, 0.5)  # 360-day phases accepted
        for variant in (:shifted_months, :shifted_one_month, :shifted_years, :repeated, :infinite_bounds,
                        :single_latitude)
            @test_throws ArgumentError R.PrescribedOzoneData(audit_fixture(joinpath(dir, "a_$variant.nc"), variant))
        end
    end

    @testset "time, latitude and ln p interpolation against an analytic field" begin
        ozone = R.PrescribedOzone(path, grid)
        for t in (DateTime(2001, 3, 10, 7), DateTime(2000, 7, 1), DateTime(2002, 12, 1))
            R._update_prescribed_ozone!(ozone, t)
            x = R._decimal_year(t)
            expect = [linear_field(p, Float64(φ), x) for p in ozone.data.pressure, φ in model_latitude]
            @test maximum(abs.(ozone.host_values .- expect) ./ expect) < 1e-5
            @test Array(ozone.columns.values) == ozone.host_values
            for (ij, p) in ((1, 7.3e4), (50, 3.1e3), (100, 12.0), (288, 9.9e4))
                @test R._layer_ozone(ozone.columns, ij, Float32(p), 1f5) ≈
                      linear_field(p, Float64(model_latitude[ij]), x) rtol = 1e-5
            end
            # Outside the source pressure range the end levels are held.
            @test R._layer_ozone(ozone.columns, 7, 1.02f5, 1.02f5) == ozone.host_values[1, 7]
            @test R._layer_ozone(ozone.columns, 7, 0.5f0, 1f5) == ozone.host_values[end, 7]
        end
    end

    @testset "coverage: endpoint hold inside the bounds, refusal outside" begin
        data = R.PrescribedOzoneData(path)
        @test R._ozone_time_weights(data, 0, DateTime(2000, 1, 1)) == (1, 1, 0.0)           # first instant
        @test R._ozone_time_weights(data, 0, DateTime(2002, 12, 31, 23)) == (36, 36, 0.0)   # last instant inside
        @test_throws ArgumentError R._ozone_time_weights(data, 0, DateTime(1999, 12, 31, 12))
        @test_throws ArgumentError R._ozone_time_weights(data, 0, DateTime(2003, 1, 2))
        i, j, w = R._ozone_time_weights(data, 0, DateTime(2001, 6, 20))
        @test j == i + 1 && 0 <= w <= 1 && data.time[i] <= R._decimal_year(DateTime(2001, 6, 20)) <= data.time[j]
    end

    @testset "fixed year: cyclic December-to-January interpolation" begin
        mpath = write_ozone(joinpath(dir, "monthly.nc"), monthly_field)[1]
        data = R.PrescribedOzoneData(mpath)
        ozone = R.PrescribedOzone(data, grid; fixed_year = 2001)
        idx = R._ozone_year_indices(data, 2001)
        phase = data.time[idx] .- 2001
        jan, dec = data.values[1, 1, idx[1]], data.values[1, 1, idx[12]]
        @test jan != dec
        for t in (DateTime(2001, 1, 1), DateTime(2003, 1, 1), DateTime(1990, 1, 1))  # any year: only the phase matters
            R._update_prescribed_ozone!(ozone, t)
            f = R._decimal_year(t) - year(t)
            w = (f - (phase[12] - 1)) / (phase[1] - (phase[12] - 1))
            @test ozone.host_values[1, 1] ≈ (1 - w) * dec + w * jan rtol = 1e-6
        end
        R._update_prescribed_ozone!(ozone, DateTime(2001, 12, 31, 18))  # after the December midpoint: towards January
        @test dec < ozone.host_values[1, 1] < jan || jan < ozone.host_values[1, 1] < dec
        R._update_prescribed_ozone!(ozone, DateTime(2001, 7, 1))
        i, j, _ = R._ozone_time_weights(data, 2001, DateTime(2001, 7, 1))
        @test (data.months[i], data.months[j]) == (6, 7)
        @test_throws ArgumentError R.PrescribedOzone(data, grid; fixed_year = 2005)
    end

    @testset "validation" begin
        base = slab_config(path)
        @test R.validate(base; check_input_files = false) === base
        @test_throws ArgumentError R.validate(with_forcing(base; ozone_path = "", ozone_fixed_year = 2001);
                                              check_input_files = false)
        @test_throws ArgumentError R.validate(with_forcing(base; radiation = :speedy_simplified); check_input_files = false)
        # With the file checked, the run must stay inside the coverage, or the fixed year must be complete.
        @test R._validate_prescribed_ozone_options(base, true, true) === nothing
        # A half-day run from 18:00 on 31 December 2002 ends after the coverage end (1 January 2003).
        @test R._validate_prescribed_ozone_options(slab_config(path; start = DateTime(2002, 12, 31, 11)), true, true) ===
              nothing
        @test_throws ArgumentError R._validate_prescribed_ozone_options(
            slab_config(path; start = DateTime(2002, 12, 31, 18)), true, true)
        @test_throws ArgumentError R._validate_prescribed_ozone_options(slab_config(path; start = DateTime(1999, 12, 1)),
                                                                        true, true)
        @test R._validate_prescribed_ozone_options(slab_config(path; start = DateTime(1950, 1, 1), ozone_fixed_year = 2001),
                                                   true, true) === nothing
        @test_throws ArgumentError R._validate_prescribed_ozone_options(slab_config(path; ozone_fixed_year = 2007),
                                                                        true, true)
        mktempdir() do ydir
            raw = R.YAML.load_file(joinpath(R.PROJECT_ROOT, "config", "production.yml"))
            raw["forcing"]["ozone_path"] = path
            raw["forcing"]["ozone_fixed_year"] = 2001
            ypath = joinpath(ydir, "c.yml")
            R.YAML.write_file(ypath, raw)
            parsed = R.load_config(ypath; check_input_files = false)
            @test parsed.forcing.ozone_path == path && parsed.forcing.ozone_fixed_year == 2001
        end
        production = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)
        @test production.forcing.ozone_path == "" && production.forcing.ozone_fixed_year == 0
        fixed = with_forcing(production; ozone_path = path, ozone_fixed_year = 2001)
        @test R._forcing_time_dependence(fixed) == R._forcing_time_dependence(production) * "_and_zonal_mean_ozone_fixed_at_year_2001" ||
              R._forcing_time_dependence(fixed) == "zonal_mean_ozone_fixed_at_year_2001"
    end

    @testset "model: off path is the analytic profile, on path follows the file" begin
        off = R.ExperimentConfig(; name = "ozone_off", device = :cpu, truncation = 7, nlayers = 8, duration_days = 0.5,
            mixed_layer_depth_m = 30.0, start_date = DateTime(2001, 3, 1),
            forcing = R.ForcingConfig(; radiation = :rrtmgp_clear_sky))
        simulation = R.build_rrtmgp(off)
        S.run!(simulation; steps = 3)
        radiation = simulation.model.longwave_radiation
        @test isnothing(radiation.ozone)
        p_lay = Array(G.layer_pressure(radiation.solver)); o3 = Array(G.volume_mixing_ratio(radiation.solver, "o3"))
        analytic = map(p -> typeof(p)(3.0e-8) + typeof(p)(7.5e-6) * exp(-(log(p / typeof(p)(1200)))^2 / typeof(p)(2.88)),
                       p_lay)
        @test o3 == analytic  # bitwise
        @test R._forcing_time_dependence(off) == "constant_for_entire_integration"

        on = slab_config(path)
        simulation = R.build_rrtmgp(on)
        S.run!(simulation; steps = 3)
        radiation = simulation.model.longwave_radiation
        ozone = radiation.ozone
        _, column_latitude = R.RG.get_londlatds(simulation.model.spectral_grid.grid)  # the model's own grid
        @test ozone isa R.PrescribedOzone
        p_lay = Array(G.layer_pressure(radiation.solver)); o3 = Array(G.volume_mixing_ratio(radiation.solver, "o3"))
        nlay = size(o3, 1)
        @test all(o3[r, ij] == R._layer_ozone(ozone.columns, ij, p_lay[r, ij], p_lay[r, ij])
                  for r in 1:nlay, ij in axes(o3, 2))
        now = simulation.variables.prognostic.clock.time
        dt = Millisecond(simulation.model.time_stepping.Δt_millisec)
        # The slab builder's clock starts on 2000-01-01 whatever start_date says, which is before the first monthly
        # midpoint: the first month is then held, so the expected field uses the clamped time.
        x = [clamp(R._decimal_year(t), ozone.data.time[1], ozone.data.time[end]) for t in (now, now - dt)]
        # Inside the source pressure range [1 Pa, 1e5 Pa] the interpolation of the linear field is exact.
        ok = [maximum(abs(o3[r, ij] - linear_field(Float64(p_lay[r, ij]), Float64(column_latitude[ij]), xi)) /
                      linear_field(Float64(p_lay[r, ij]), Float64(column_latitude[ij]), xi)
                      for r in 1:nlay, ij in axes(o3, 2) if 1.0 <= p_lay[r, ij] <= 1.0e5) for xi in x]
        println("OZONE_MODEL_LAYER_MAX_REL_ERROR ", minimum(ok))
        @test minimum(ok) < 1e-5
        @test all(isfinite, Array(S.get_step(simulation.variables.prognostic.temperature, 1).data))
        provenance = R._forcing_provenance(radiation, on)
        @test provenance.ozone_sha256 == ozone.data.sha256 && provenance.ozone_fixed_year == 0
        @test isfinite(provenance.ozone_final_model_column_du_global_mean)
        @test R._forcing_time_dependence(on) == R.PRESCRIBED_OZONE_TIME_DEPENDENCE
    end

    @testset "column ozone: constant mixing ratio and the same pressure interval" begin
        cpath = write_ozone(joinpath(dir, "const.nc"), (p, φ, t) -> 1.0e-6)[1]
        simulation = R.build_rrtmgp(slab_config(cpath))
        S.run!(simulation; steps = 2)
        solver = simulation.model.longwave_radiation.solver
        o3 = Array(G.volume_mixing_ratio(solver, "o3")); p_lev = Float64.(Array(G.level_pressure(solver)))
        @test all(==(1.0f-6), o3)
        column = R._model_column_ozone_du(o3, p_lev, 9.80665, 0.0289644)
        analytic = Float64(1.0f-6) .* (p_lev[1, :] .- p_lev[end, :]) ./ 9.80665 ./ 0.0289644 .* 6.02214076e23 ./ 2.6867e20
        @test column ≈ analytic rtol = 1e-12
        # Linear-in-ln p profile: model-layer column against the exact integral over the same pressure interval.
        simulation = R.build_rrtmgp(slab_config(path))
        S.run!(simulation; steps = 2)
        radiation = simulation.model.longwave_radiation
        o3 = Float64.(Array(G.volume_mixing_ratio(radiation.solver, "o3")))
        p_lev = Float64.(Array(G.level_pressure(radiation.solver)))
        layer = R._model_column_ozone_du(o3, p_lev, 9.80665, 0.0289644)
        x = clamp(R._decimal_year(simulation.variables.prognostic.clock.time), radiation.ozone.data.time[1],
                  radiation.ozone.data.time[end])
        _, column_latitude = R.RG.get_londlatds(simulation.model.spectral_grid.grid)
        # The implemented profile: a + b ln p inside [1 Pa, 1e5 Pa], held at the end levels outside.
        # ∫ (a + b ln p) dp = a p + b (p ln p - p), with a including the latitude and time terms.
        exact = map(axes(o3, 2)) do ij
            a = 5e-6 + 1e-8 * Float64(column_latitude[ij]) + 1e-7 * (x - 2000); b = 2e-7
            F(p) = a * p + b * (p * log(p) - p)
            bottom, top = p_lev[1, ij], max(p_lev[end, ij], 1.0)
            inside = F(min(bottom, 1.0e5)) - F(top)
            held = max(bottom - 1.0e5, 0.0) * (a + b * log(1.0e5))
            (inside + held) / 9.80665 / 0.0289644 * 6.02214076e23 / 2.6867e20
        end
        quadrature = maximum(abs.(layer .- exact) ./ exact)
        println("OZONE_COLUMN_QUADRATURE_MAX_REL_ERROR_L8 ", quadrature)
        @test quadrature < 0.05  # midpoint quadrature over 8 layers and the time step; recorded, not tuned
    end

    @testset "restart record" begin
        small = S.SpectralGrid(; NF = Float32, trunc = 1, nlayers = 8)
        mpath = joinpath(dir, "monthly.nc")
        rad = R.RRTMGPRadiation(small; prescribed_ozone = R.PrescribedOzone(path, small))
        other = R.RRTMGPRadiation(small; prescribed_ozone = R.PrescribedOzone(mpath, small))
        fixed = R.RRTMGPRadiation(small; prescribed_ozone = R.PrescribedOzone(path, small; fixed_year = 2001))
        base = R.RRTMGPRadiation(small)
        saved = R._rrtmgp_radiation_state(rad)
        @test saved.prescribed_ozone == (sha256 = rad.ozone.data.sha256, fixed_year = 0)
        @test R._restore_rrtmgp_radiation_state!(rad, saved) === rad
        legacy = (; (k => v for (k, v) in pairs(saved) if k != :prescribed_ozone)...)
        for (radiation, state) in ((rad, legacy), (other, saved), (fixed, saved), (base, saved),
                                   (rad, R._rrtmgp_radiation_state(base)))
            before = radiation.call_counter
            @test_throws ArgumentError R._restore_rrtmgp_radiation_state!(radiation, merge(state, (; call_counter = 99)))
            @test radiation.call_counter == before
        end
        @test !hasproperty(R._rrtmgp_radiation_state(base), :prescribed_ozone)
        @test R._restore_rrtmgp_radiation_state!(base, R._rrtmgp_radiation_state(base)) === base
    end
end
end
println("PRESCRIBED_OZONE_UNIT_PASS")
