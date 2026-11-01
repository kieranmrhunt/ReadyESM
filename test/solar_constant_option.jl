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
    on = with_forcing(production; solar_constant_w_m2 = CMIP6, solar_constant_source = LABEL)
    @test R.validate(on; check_input_files = false) === on
    @test R._solar_constant_source(on) == LABEL
    @test_throws ArgumentError R.validate(with_forcing(production; solar_constant_w_m2 = CMIP6); check_input_files = false)
    for bad in (0.0, -1.0, NaN, Inf)
        @test_throws ArgumentError R.validate(with_forcing(on; solar_constant_w_m2 = bad); check_input_files = false)
    end
    @test R.validate(with_forcing(production; solar_constant_source = "stated default"); check_input_files = false) isa
          R.ExperimentConfig
    # Restart identity: the setting is part of the configuration hash, so a restart across a change is refused.
    @test R._dynamic_config_identity_hash(on) != R._dynamic_config_identity_hash(production)

    mktempdir() do directory
        raw = R.YAML.load_file(joinpath(R.PROJECT_ROOT, "config", "production.yml"))
        raw["forcing"]["solar_constant_w_m2"] = CMIP6
        raw["forcing"]["solar_constant_source"] = LABEL
        path = joinpath(directory, "tsi.yml")
        R.YAML.write_file(path, raw)
        loaded = R.load_config(path; check_input_files = false)
        @test loaded.forcing.solar_constant_w_m2 == CMIP6 && loaded.forcing.solar_constant_source == LABEL
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
    on = with_forcing(base; solar_constant_w_m2 = CMIP6, solar_constant_source = LABEL)
    model_off, off, meta_off = run_case(base, 8)
    model_on, on_record, meta_on = run_case(with(on; output_dir = mktempdir()), 8)
    @test meta_off.solar_constant_w_m2 == 1365.0 && meta_off.solar_constant_source == "SpeedyWeather Earth default"
    @test meta_on.solar_constant_w_m2 == CMIP6 && meta_on.solar_constant_source == LABEL
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
