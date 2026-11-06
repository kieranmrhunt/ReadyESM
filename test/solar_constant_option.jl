module SolarConstantOptionTests
# Opt-in total solar irradiance (design artifacts/development_wave_20261101/tsi_option_claude/DESIGN.md v2): settings,
# provenance, the shared planet in every PrimitiveWetModel, restart identity, and a short coupled run with the default
# and with 1360.747 W/m2. The exact on/off check is on the irradiance RRTMGP receives (toa_flux = S0 f(t), typed in
# the model's number format, 2 eps). The global incoming ratio and the global incoming against S0 f(t)/4 are printed:
# at T7 the grid integral of the cosine of the zenith angle is too coarse for the 0.05 W/m2 bound of the design, which
# is checked at T31 in a run, not here.
using ReadyESM, Test, Dates
const R = ReadyESM
const S = R.SpeedyWeather
const O = R.Oceananigans

fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
with(config; kwargs...) = R.ExperimentConfig(; merge(fields(config), (; kwargs...))...)
with_forcing(config; kwargs...) =
    with(config; forcing = R.ForcingConfig(; merge(fields(config.forcing), (; kwargs...))...))
include(joinpath(@__DIR__, "analytic_coupled_fixture.jl"))

const CMIP6 = 1360.747
const LABEL = "CMIP6 piControl (Matthes et al. 2017)"

"Text of the call starting at `start` up to its matching closing parenthesis."
function call_text(source, start)
    depth = 0
    for i in start:lastindex(source)
        c = source[i]
        c == '(' && (depth += 1)
        if c == ')'
            depth -= 1
            depth == 0 && return source[start:i]
        end
    end
    error("unbalanced call")
end

@testset "Solar constant option: settings, provenance, planet" begin
    default = R.ExperimentConfig()
    @test default.forcing.solar_constant_w_m2 == 1365.0 == R.DEFAULT_SOLAR_CONSTANT_W_M2
    @test default.forcing.solar_constant_source == ""
    @test R._solar_constant_source(default) == "SpeedyWeather Earth default"
    production = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)
    @test production.forcing.solar_constant_w_m2 == 1365.0
    on = with_forcing(production; solar_constant_w_m2 = CMIP6, solar_constant_source = LABEL,
        solar_constant_entry = "picontrol")
    @test R.validate(on; check_input_files = false) === on
    @test R._solar_constant_source(on) == LABEL
    @test_throws ArgumentError R.validate(with_forcing(production; solar_constant_w_m2 = CMIP6); check_input_files = false)
    for bad in (0.0, -1.0, NaN, Inf)
        @test_throws ArgumentError R.validate(with_forcing(on; solar_constant_w_m2 = bad); check_input_files = false)
    end
    # The default value carries no source label of its own (Barry, port read note 2).
    @test_throws ArgumentError R.validate(with_forcing(production; solar_constant_source = "stated default");
        check_input_files = false)
    @test R.validate(with_forcing(production; solar_constant_source = "  "); check_input_files = false) isa
          R.ExperimentConfig
    # A non-default value must name a pinned entry, agree with its value and carry its label (design v2 amendment 2).
    invalid(; kwargs...) = @test_throws ArgumentError R.validate(with_forcing(on; kwargs...); check_input_files = false)
    invalid(solar_constant_entry = "")                                  # a free label alone is not enough
    invalid(solar_constant_entry = "piControl")                         # unknown entry
    invalid(solar_constant_w_m2 = 1361.079)                             # the 1993 value under the piControl entry
    invalid(solar_constant_w_m2 = CMIP6 + 2e-3)                         # outside the table's three-decimal precision
    invalid(solar_constant_source = "CMIP6 piControl")                  # not the entry's label
    @test_throws ArgumentError R.validate(with_forcing(production; solar_constant_entry = "picontrol");
        check_input_files = false)                                      # an entry with the default value
    historical = R.SOLAR_CONSTANT_ENTRIES["historical_1993"]
    hist = with_forcing(production; solar_constant_w_m2 = 1361.079, solar_constant_source = historical.label,
        solar_constant_entry = "historical_1993")
    @test R.validate(hist; check_input_files = false) === hist
    @test R._solar_constant_provenance(hist).solar_constant_year == "1993"
    # The constants in src/solar_forcing.jl are the pinned table's values.
    table_path = joinpath(R.PROJECT_ROOT, R.SOLAR_CONSTANT_TABLE)
    @test bytes2hex(R.SHA.sha256(read(table_path))) == R.SOLAR_CONSTANT_TABLE_SHA256
    table = read(table_path, String)
    @test sort(collect(keys(R.SOLAR_CONSTANT_ENTRIES))) == ["historical_1993", "picontrol"]
    @test occursin("tsi_w_m2 = $(repr(R.SOLAR_CONSTANT_ENTRIES["picontrol"].value_w_m2))", table)
    @test occursin("annual_mean_tsi_w_m2 = $(repr(historical.value_w_m2))", table)
    @test occursin("reduction = \"$(historical.reduction)\"", table)
    for (name, pinned) in R.SOLAR_CONSTANT_ENTRIES
        @test occursin("label = \"$(pinned.label)\"", table)
        @test occursin("path = \"$(pinned.archive_path)\"", table)
        @test occursin("sha256 = \"$(pinned.archive_sha256)\"", table)
        @test abs(round(pinned.value_w_m2; digits = 3) - pinned.value_w_m2) <= R.SOLAR_CONSTANT_ENTRY_TOLERANCE_W_M2
    end
    provenance = R._solar_constant_provenance(on)
    @test provenance.solar_constant_entry == "picontrol" && provenance.solar_constant_year == "" &&
          provenance.solar_constant_table_sha256 == R.SOLAR_CONSTANT_TABLE_SHA256 &&
          provenance.solar_constant_archive_sha256 == R.SOLAR_CONSTANT_ENTRIES["picontrol"].archive_sha256
    @test R._solar_constant_entry_provenance(production) == (;)
    # The output check accepts matching provenance, and refuses a default with an entry or a changed SHA-256.
    meta(config) = (; solar_constant_w_m2 = config.forcing.solar_constant_w_m2,
        solar_constant_source = R._solar_constant_source(config), R._solar_constant_entry_provenance(config)...)
    @test R._validate_solar_constant_provenance(meta(production)) === nothing
    @test R._validate_solar_constant_provenance(meta(on)) === nothing
    @test R._validate_solar_constant_provenance(meta(hist)) === nothing
    @test_throws ErrorException R._validate_solar_constant_provenance(
        merge(meta(production), (; solar_constant_entry = "picontrol")))
    @test_throws ErrorException R._validate_solar_constant_provenance(
        merge(meta(production), (; solar_constant_source = "stated default")))
    @test_throws ErrorException R._validate_solar_constant_provenance(
        merge(meta(on), (; solar_constant_archive_sha256 = "0"^64)))
    @test_throws ErrorException R._validate_solar_constant_provenance(
        (; solar_constant_w_m2 = CMIP6, solar_constant_source = LABEL))
    # Restart identity: the setting is part of the configuration hash, so a restart across a change is refused.
    @test R._dynamic_config_identity_hash(on) != R._dynamic_config_identity_hash(production)

    mktempdir() do directory
        raw = R.YAML.load_file(joinpath(R.PROJECT_ROOT, "config", "production.yml"))
        raw["forcing"]["solar_constant_w_m2"] = CMIP6
        raw["forcing"]["solar_constant_source"] = LABEL
        raw["forcing"]["solar_constant_entry"] = "picontrol"
        path = joinpath(directory, "tsi.yml")
        R.YAML.write_file(path, raw)
        loaded = R.load_config(path; check_input_files = false)
        @test loaded.forcing.solar_constant_w_m2 == CMIP6 && loaded.forcing.solar_constant_source == LABEL &&
              loaded.forcing.solar_constant_entry == "picontrol"
    end

    for NF in (Float32, Float64)
        grid = S.SpectralGrid(; NF, trunc = 7, nlayers = 8)
        reference = S.Earth(grid)
        planet = R._readyesm_planet(production, grid)
        @test all(isequal(getfield(planet, n), getfield(reference, n)) for n in fieldnames(typeof(reference)))
        @test R._readyesm_planet(on, grid).solar_constant === NF(CMIP6)
    end

    # Every PrimitiveWetModel built in src must take the configured planet.
    for file in readdir(joinpath(R.PROJECT_ROOT, "src"); join = true)
        endswith(file, ".jl") || continue
        source = read(file, String)
        for hit in findall("PrimitiveWetModel(", source)
            @test occursin("planet = _readyesm_planet(config, spectral_grid)", call_text(source, last(hit)))
        end
    end
    @test count("PrimitiveWetModel(", read(joinpath(R.PROJECT_ROOT, "src", "baseline.jl"), String)) == 4
end

function run_case(config, steps)
    simulation = R.build_dynamic_esm(config; steps)
    model = simulation.model.atmosphere.model
    solver = model.longwave_radiation.solver
    weights = Float64.(Array(R._global_point_weights(model.spectral_grid)))
    record = Dict{Int, NamedTuple}()
    function sample(sim)
        i = sim.model.clock.iteration
        i > 0 || return nothing
        toa = Array(R.RRTMGP.toa_flux(solver))
        incoming = Float64.(Array(R.RRTMGP.sw_flux_dn(solver))[end, :])
        record[i] = (; toa, global_incoming = sum(weights .* incoming) / sum(weights),
            time = sim.model.atmosphere.variables.prognostic.clock.time)
        return nothing
    end
    O.add_callback!(simulation, sample, O.IterationInterval(1); name = :tsi_sample)
    O.run!(simulation)
    # Collect, save and validate the production diagnostics with the setting (the M1 lesson), and read its provenance.
    diagnostics = R.collect_dynamic_diagnostics(simulation, config)
    R.save_dynamic_diagnostics(diagnostics, config; render_figure = false)
    R.validate_dynamic_diagnostics(diagnostics)
    return model, record, diagnostics.metadata
end

@testset "Solar constant option: coupled analytic run, default and 1360.747" begin
    production = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)
    base = with(analytic_coupled_test_config(production; output_dir = mktempdir());
        atmosphere_timestep_at_t31_minutes = 3.75, duration_days = 8 * 900 / 86_400, start_date = DateTime(1993, 1, 1))
    startswith(String(base.forcing.radiation), "rrtmgp") || error("production radiation is not RRTMGP")
    base.forcing.solar_distance_correction || error("production solar distance correction is off")
    on = with_forcing(base; solar_constant_w_m2 = CMIP6, solar_constant_source = LABEL, solar_constant_entry = "picontrol")
    model_off, off, meta_off = run_case(base, 8)
    model_on, on_record, meta_on = run_case(with(on; output_dir = mktempdir()), 8)
    @test meta_off.solar_constant_w_m2 == 1365.0 && meta_off.solar_constant_source == "SpeedyWeather Earth default"
    @test meta_on.solar_constant_w_m2 == CMIP6 && meta_on.solar_constant_source == LABEL
    # The saved provenance carries the pinned entry only when the setting is used.
    @test !hasproperty(meta_off, :solar_constant_entry)
    @test meta_on.solar_constant_entry == "picontrol" &&
          meta_on.solar_constant_archive_sha256 == R.SOLAR_CONSTANT_ENTRIES["picontrol"].archive_sha256 &&
          meta_on.solar_constant_table_sha256 == R.SOLAR_CONSTANT_TABLE_SHA256
    NF = typeof(model_on.planet.solar_constant)
    @test model_off.planet.solar_constant === NF(1365)
    @test model_on.planet.solar_constant === NF(CMIP6)
    @test sort(collect(keys(off))) == sort(collect(keys(on_record))) == collect(1:8)
    expected = NF(CMIP6) / NF(1365)
    # Steps after a radiation solve carry a set irradiance; require at least two of them.
    active = [i for i in 1:8 if all(>(0), off[i].toa) && all(>(0), on_record[i].toa)]
    println("TSI_ACTIVE_STEPS ", active)
    @test length(active) >= 2
    for i in active
        a, b = on_record[i], off[i]
        @test a.time == b.time
        @test all(isfinite, a.toa) && all(isfinite, b.toa)
        @test all(abs.(a.toa ./ b.toa .- expected) .<= 2 * eps(NF))       # the typed irradiance RRTMGP receives
        f = R._solar_distance_factor(a.time)
        println("TSI_INSTANT step=$i global_incoming_ratio=$(a.global_incoming / b.global_incoming) ",
            "expected=$(Float64(expected)) on_minus_S0f4=$(a.global_incoming - CMIP6 * f / 4) ",
            "off_minus_S0f4=$(b.global_incoming - 1365 * f / 4)")
        @test isapprox(a.global_incoming / b.global_incoming, Float64(expected); rtol = 1e-5)
    end
end
println("SOLAR_CONSTANT_OPTION_TESTS_PASS")
end
