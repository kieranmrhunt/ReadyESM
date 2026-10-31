module WellMixedGasScheduleTests
using ReadyESM, Test, Dates
const R = ReadyESM
const G = R.RRTMGP

"""Write a schedule file; `values` maps each gas to a vector over `years` (mol mol-1)."""
function write_schedule(path, years, values; schema = R.WELL_MIXED_GAS_SCHEMA, units = "mol mol-1", scenario = "test")
    R.NCDatasets.NCDataset(path, "c") do ds
        ds.attrib["schema"] = schema
        ds.attrib["scenario"] = scenario
        R.NCDatasets.defDim(ds, "year", length(years))
        R.NCDatasets.defVar(ds, "year", Int32.(years), ("year",))
        for name in R.WELL_MIXED_GAS_NAMES
            v = R.NCDatasets.defVar(ds, String(name), Float64.(values[name]), ("year",))
            v.attrib["units"] = units
            v.attrib["rrtmgp_gas"] = R.WELL_MIXED_GAS_RRTMGP_NAMES[name]
        end
    end
    return path
end

linear(years, a, b) = [a + b * (y - first(years)) for y in years]
fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
with_forcing(config; kwargs...) = R.ExperimentConfig(; merge(fields(config),
    (; forcing = R.ForcingConfig(; merge(fields(config.forcing), (; kwargs...))...)))...)
with(config; kwargs...) = R.ExperimentConfig(; merge(fields(config), (; kwargs...))...)
# v0.1.16 and later production turns on the seasonal solar distance, which the time-dependence label appends.
solar_suffix(config) = config.forcing.solar_distance_correction ? "_and_seasonal_solar_distance_from_model_clock" : ""

@testset "Time-varying well-mixed greenhouse gases" begin
    dir = mktempdir()
    years = 1990:2010
    values = (co2 = linear(years, 350e-6, 2e-6), ch4 = linear(years, 1700e-9, 5e-9), n2o = linear(years, 300e-9, 1e-9),
              cfc12 = linear(years, 500e-12, 1e-12), cfc11eq = linear(years, 700e-12, 2e-12))
    path = write_schedule(joinpath(dir, "schedule.nc"), years, values)
    production = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)

    @testset "off: production gases and objects unchanged" begin
        @test production.forcing.well_mixed_gas_path == "" && production.forcing.well_mixed_gas_fixed_year == 0
        @test isnothing(R._well_mixed_gas_schedule(production))
        grid = R.SpeedyWeather.SpectralGrid(; NF = Float32, trunc = 7, nlayers = 8)
        base = R.RRTMGPRadiation(grid)
        @test isnothing(base.gas_schedule) && isnothing(base.cfc12_vmr) && isnothing(base.cfc11_vmr)
        @test G.volume_mixing_ratio(base.solver, "cfc11") == 0 && G.volume_mixing_ratio(base.solver, "cfc12") == 0
        @test R._forcing_time_dependence(production) == (production.forcing.solar_distance_correction ?
            "seasonal_solar_distance_from_model_clock" : "constant_for_entire_integration")
        @test isempty(R._well_mixed_gas_provenance(base, production, true))
        base_state = R._rrtmgp_radiation_state(base)
        @test keys(base_state) == (:call_counter, :ch4_vmr, :n2o_vmr, :solar_distance_correction,
            :heating_rate, :applied_heating_rate, :lw_up, :lw_down, :lw_net, :sw_up, :sw_down,
            :sw_net, :sw_direct_down, :net, :clear, :sulfate_mass, :cloud_sampler,
            :prognostic_cloud_condensate)
        @test !hasproperty(base_state, :well_mixed_gas_schedule)
        @test !hasproperty(base_state, :well_mixed_gas_values)
    end

    @testset "reader" begin
        s = R.WellMixedGasSchedule(path)
        @test s.years == collect(years) && s.scenario == "test" && length(s.sha256) == 64
        @test s.values[1, 1] == 350e-6 && s.values[end, 5] == values.cfc11eq[end]
        bad_schema = write_schedule(joinpath(dir, "bad_schema.nc"), years, values; schema = "other")
        bad_units = write_schedule(joinpath(dir, "bad_units.nc"), years, values; units = "ppm")
        negative = write_schedule(joinpath(dir, "negative.nc"), years, merge(values, (; cfc12 = -values.cfc12)))
        zero_co2 = write_schedule(joinpath(dir, "zero_co2.nc"), years, merge(values, (; co2 = 0 .* values.co2)))
        gap = write_schedule(joinpath(dir, "gap.nc"), [1990, 1992], map(v -> v[1:2], values))
        for bad in (bad_schema, bad_units, negative, zero_co2, gap, joinpath(dir, "missing.nc"))
            @test_throws ArgumentError R.WellMixedGasSchedule(bad)
        end
        @test_throws ArgumentError R.WellMixedGasSchedule(path; fixed_year = 1980)
    end

    @testset "interpolation" begin
        s = R.WellMixedGasSchedule(path)
        # 1993 has 365 days: 2 July 12:00 is decimal year 1993.5 exactly, the 1993 annual mean.
        at_mid = R.well_mixed_gas_values(s, DateTime(1993, 7, 2, 12))
        @test all(at_mid[g] == values[g][1993 - 1989] for g in R.WELL_MIXED_GAS_NAMES)
        # 1 January 1994 is halfway between the 1993 and 1994 annual means.
        new_year = R.well_mixed_gas_values(s, DateTime(1994))
        @test all(isapprox(new_year[g], (values[g][4] + values[g][5]) / 2; rtol = 1e-14) for g in R.WELL_MIXED_GAS_NAMES)
        # Linear in time between mid-years.
        t = DateTime(1994, 3, 15)
        w = R._decimal_year(t) - 1993.5
        @test isapprox(R.well_mixed_gas_values(s, t).co2, (1 - w) * values.co2[4] + w * values.co2[5]; rtol = 1e-14)
        @test R._decimal_year(DateTime(2000)) == 2000.0
        @test R.well_mixed_gas_values(s, DateTime(2010, 7, 2, 12)).co2 == values.co2[end]
        @test_throws ArgumentError R.well_mixed_gas_values(s, DateTime(1990, 6, 1))
        @test_throws ArgumentError R.well_mixed_gas_values(s, DateTime(2010, 12, 1))
        fixed = R.WellMixedGasSchedule(path; fixed_year = 1995)
        @test R.well_mixed_gas_values(fixed, DateTime(1700)).co2 == values.co2[6]
        @test R.well_mixed_gas_values(fixed, DateTime(2500)).cfc11eq == values.cfc11eq[6]
    end

    @testset "validation" begin
        on = with_forcing(production; well_mixed_gas_path = path)
        # Production initializes from January 1993 ERA5/ECCO, so keep that start month.
        on = with(on; start_date = DateTime(1993), duration_days = 365.0, analysis_start_day = 0.0)
        @test R.validate(on; check_input_files = false) === on
        @test_throws ArgumentError R.validate(with_forcing(on; ch4_ppb = 1800.0); check_input_files = false)
        @test_throws ArgumentError R.validate(with_forcing(on; radiation = :speedy_simplified); check_input_files = false)
        @test_throws ArgumentError R.validate(with_forcing(production; well_mixed_gas_fixed_year = 1850);
                                              check_input_files = false)
        @test_throws ArgumentError R._validate_well_mixed_gas_options(with(on; carbon = R.CarbonConfig(; mode = :prescribed)),
                                                                      true, false)
        # The run must stay inside the interpolable span when the file is checked.
        @test R._validate_well_mixed_gas_options(on, true, true) === nothing
        @test_throws ArgumentError R._validate_well_mixed_gas_options(with(on; duration_days = 7000.0), true, true)
        @test_throws ArgumentError R._validate_well_mixed_gas_options(with(on; start_date = DateTime(1990, 3, 1)), true, true)
        fixed = with_forcing(on; well_mixed_gas_fixed_year = 1995)
        @test R._validate_well_mixed_gas_options(with(fixed; duration_days = 40000.0), true, true) === nothing
        @test R._forcing_time_dependence(on) == R.WELL_MIXED_GAS_TIME_DEPENDENCE * solar_suffix(on)
        @test R._forcing_time_dependence(fixed) == "well_mixed_greenhouse_gases_fixed_at_year_1995" * solar_suffix(fixed)
        mktempdir() do ydir
            raw = R.YAML.load_file(joinpath(R.PROJECT_ROOT, "config", "production.yml"))
            raw["forcing"]["well_mixed_gas_path"] = path
            raw["forcing"]["well_mixed_gas_fixed_year"] = 1995
            ypath = joinpath(ydir, "c.yml")
            R.YAML.write_file(ypath, raw)
            parsed = R.load_config(ypath; check_input_files = false)
            @test parsed.forcing.well_mixed_gas_path == path && parsed.forcing.well_mixed_gas_fixed_year == 1995
        end
    end

    @testset "radiation: halocarbons and CO2 reach RRTMGP" begin
        grid = R.SpeedyWeather.SpectralGrid(; NF = Float32, trunc = 7, nlayers = 8)
        s = R.WellMixedGasSchedule(path; fixed_year = 2000)
        rad = R.RRTMGPRadiation(grid; well_mixed_gas_schedule = s, well_mixed_gas_start = DateTime(2000))
        @test rad.co2_ppm == Float32(values.co2[11] * 1e6)
        for (gas, field) in (("ch4", :ch4_vmr), ("n2o", :n2o_vmr), ("cfc12", :cfc12_vmr), ("cfc11", :cfc11_vmr))
            @test G.volume_mixing_ratio(rad.solver, gas) == getfield(rad, field)
        end
        @test G.volume_mixing_ratio(rad.solver, "co2") == Float32(values.co2[11])
        @test_throws ArgumentError R.RRTMGPRadiation(grid; well_mixed_gas_schedule = s)
        @test_throws ArgumentError R.RRTMGPRadiation(grid; well_mixed_gas_schedule = s,
                                                     well_mixed_gas_start = DateTime(2000), ch4_ppb = 1800.0)
        # Same standard profile: removing the halocarbons raises outgoing longwave (they absorb), by a few tenths W/m2.
        toa_up(radiation) = sum(Array(G.lw_flux_up(radiation.solver))[end, :]) / size(G.lw_flux_up(radiation.solver), 2)
        with_halocarbons = toa_up(rad)
        G.set_volume_mixing_ratio!(rad.solver, "cfc11", 0.0f0)
        G.set_volume_mixing_ratio!(rad.solver, "cfc12", 0.0f0)
        G.update_fluxes!(rad.solver)
        without = toa_up(rad)
        println("WELL_MIXED_GAS_HALOCARBON_OLR_EFFECT_Wm2 ", without - with_halocarbons)
        @test 0.05 < without - with_halocarbons < 1.0
        # Doubling CO2 in the schedule lowers outgoing longwave by more than 1 W/m2 on the same profile.
        doubled = write_schedule(joinpath(dir, "doubled.nc"), years, merge(values, (; co2 = 2 .* values.co2)))
        rad2 = R.RRTMGPRadiation(grid; well_mixed_gas_schedule = R.WellMixedGasSchedule(doubled; fixed_year = 2000),
                                 well_mixed_gas_start = DateTime(2000))
        rad1 = R.RRTMGPRadiation(grid; well_mixed_gas_schedule = s, well_mixed_gas_start = DateTime(2000))
        println("WELL_MIXED_GAS_CO2_DOUBLING_OLR_EFFECT_Wm2 ", toa_up(rad1) - toa_up(rad2))
        @test toa_up(rad1) - toa_up(rad2) > 1.0
    end

    @testset "model: gases follow the clock through the radiation solve" begin
        config = R.ExperimentConfig(; name = "wmghg", device = :cpu, truncation = 7, nlayers = 8, duration_days = 0.5,
            mixed_layer_depth_m = 30.0, start_date = DateTime(2000, 6, 30),
            forcing = R.ForcingConfig(; radiation = :rrtmgp_clear_sky, well_mixed_gas_path = path))
        @test R.validate(config; check_input_files = false) === config
        simulation = R.build_rrtmgp(config)
        R.SpeedyWeather.run!(simulation; steps = 4)
        radiation = simulation.model.longwave_radiation
        @test radiation.gas_schedule isa R.WellMixedGasSchedule
        now = simulation.variables.prognostic.clock.time
        dt = Millisecond(simulation.model.time_stepping.Δt_millisec)
        candidates = [R.well_mixed_gas_values(radiation.gas_schedule, t) for t in (now, now - dt)]
        @test any(radiation.co2_ppm == Float32(c.co2 * 1e6) && radiation.cfc11_vmr == Float32(c.cfc11eq) for c in candidates)
        @test G.volume_mixing_ratio(radiation.solver, "cfc11") == radiation.cfc11_vmr
        @test G.volume_mixing_ratio(radiation.solver, "n2o") == radiation.n2o_vmr
        @test radiation.co2_ppm != Float32(R.well_mixed_gas_values(radiation.gas_schedule, config.start_date).co2 * 1e6)
        provenance = R._forcing_provenance(radiation, config)
        @test provenance.well_mixed_gas_sha256 == radiation.gas_schedule.sha256
        @test provenance.well_mixed_gas_final_co2_ppm == Float64(radiation.co2_ppm)
        @test provenance.requested_co2_ppm_ignored
        @test all(isfinite, Array(R.SpeedyWeather.get_step(simulation.variables.prognostic.temperature, 1).data))
    end

    @testset "restart record" begin
        grid = R.SpeedyWeather.SpectralGrid(; NF = Float32, trunc = 1, nlayers = 8)
        s = R.WellMixedGasSchedule(path)
        other_path = write_schedule(joinpath(dir, "other.nc"), years, merge(values, (; ch4 = 1.01 .* values.ch4)))
        rad = R.RRTMGPRadiation(grid; well_mixed_gas_schedule = s, well_mixed_gas_start = DateTime(2000))
        other = R.RRTMGPRadiation(grid; well_mixed_gas_schedule = R.WellMixedGasSchedule(other_path),
                                  well_mixed_gas_start = DateTime(2000))
        fixed = R.RRTMGPRadiation(grid; well_mixed_gas_schedule = R.WellMixedGasSchedule(path; fixed_year = 2000),
                                  well_mixed_gas_start = DateTime(2000))
        base = R.RRTMGPRadiation(grid)
        saved = R._rrtmgp_radiation_state(rad)
        @test saved.well_mixed_gas_schedule == (sha256 = s.sha256, fixed_year = 0)
        # Gas values at the checkpoint differ from the reconstructed start values; the identity check still passes.
        later = merge(saved, (; ch4_vmr = rad.ch4_vmr * 1.1f0, n2o_vmr = rad.n2o_vmr * 1.1f0))
        @test R._restore_rrtmgp_radiation_state!(rad, later) === rad
        # The schedule-active record also carries the current gas values, and a restore reinstates them exactly.
        @test keys(saved.well_mixed_gas_values) == (:co2_ppm, :ch4_vmr, :n2o_vmr, :cfc12_vmr, :cfc11_vmr)
        @test !hasproperty(R._rrtmgp_radiation_state(base), :well_mixed_gas_values)
        shifted = map(v -> v * 1.05f0, saved.well_mixed_gas_values)
        @test R._restore_rrtmgp_radiation_state!(rad, merge(saved, (; well_mixed_gas_values = shifted))) === rad
        @test (rad.co2_ppm, rad.ch4_vmr, rad.n2o_vmr, rad.cfc12_vmr, rad.cfc11_vmr) == Tuple(shifted)
        @test R._rrtmgp_radiation_state(rad).well_mixed_gas_values == shifted
        # A schedule-active checkpoint without the values is refused before any state changes.
        no_values = (; (k => v for (k, v) in pairs(saved) if k != :well_mixed_gas_values)...)
        before_counter = rad.call_counter
        @test_throws ArgumentError R._restore_rrtmgp_radiation_state!(rad, merge(no_values, (; call_counter = 77)))
        @test rad.call_counter == before_counter && rad.co2_ppm == shifted.co2_ppm
        legacy = (; (k => v for (k, v) in pairs(saved) if k != :well_mixed_gas_schedule)...)
        for (radiation, state) in ((rad, legacy), (other, saved), (fixed, saved), (base, saved),
                                   (rad, R._rrtmgp_radiation_state(base)))
            before = radiation.call_counter
            @test_throws ArgumentError R._restore_rrtmgp_radiation_state!(radiation, merge(state, (; call_counter = 99)))
            @test radiation.call_counter == before
        end
        # A constant-gas checkpoint without the record still restores into a constant-gas configuration.
        @test !hasproperty(R._rrtmgp_radiation_state(base), :well_mixed_gas_schedule)
        @test R._restore_rrtmgp_radiation_state!(base, R._rrtmgp_radiation_state(base)) === base
    end

    prepared = "/gws/ssde/j25a/ncas_climate/vol1/users/kieran/readiesm/artifacts/development_wave_20261023/ssp_forcing_v1/inputs/wmghg_cmip6_ssp245_1850-2100.nc"
    if isfile(prepared)
        @testset "prepared CMIP6 SSP2-4.5 file" begin
            s = R.WellMixedGasSchedule(prepared)
            @test s.years == collect(1850:2100) && s.scenario == "ssp245"
            @test isapprox(R.well_mixed_gas_values(s, DateTime(2100, 7, 2, 12)).co2 * 1e6, 602.78; atol = 0.01)
            @test isapprox(R.well_mixed_gas_values(s, DateTime(1850, 7, 2, 12)).co2 * 1e6, 284.32; atol = 0.01)
            @test R.well_mixed_gas_values(s, DateTime(1850, 7, 2, 12)).cfc12 == 0
        end
    end
end
end
println("WELL_MIXED_GAS_SCHEDULE_UNIT_PASS")
