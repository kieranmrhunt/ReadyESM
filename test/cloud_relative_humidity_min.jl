@testset "Cloud relative-humidity threshold option" begin
    R = ReadyESM
    default = R.ExperimentConfig()
    @test default.forcing.cloud_relative_humidity_min == 0.3

    fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
    with_forcing(c; kw...) = R.ExperimentConfig(;
        merge(fields(c), (; forcing = R.ForcingConfig(; merge(fields(c.forcing), (; kw...))...)))...)

    # The released production file still selects the SPEEDY value.
    production = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)
    @test production.forcing.cloud_relative_humidity_min == 0.3

    # A YAML value reaches the config.
    mktempdir() do dir
        text = read(joinpath(R.PROJECT_ROOT, "config", "production.yml"), String)
        path = joinpath(dir, "rh.yml")
        write(path, replace(text, "  cloud_scheme: diagnostic\n" =>
            "  cloud_scheme: diagnostic\n  cloud_relative_humidity_min: 0.7\n"))
        loaded = R.load_config(path; check_input_files = false)
        @test loaded.forcing.cloud_relative_humidity_min == 0.7
    end

    # Validation: inside [0, 1), and only with an all-sky cloud scheme.
    diag = with_forcing(production; cloud_relative_humidity_min = 0.7)
    @test R.validate(diag; check_input_files = false) === diag
    @test Float32(1.0 - 1e-9) == 1.0f0
    for bad in (-0.1, 1.0, 1.5, 1.0 - 1e-9)
        @test_throws ArgumentError R.validate(with_forcing(production; cloud_relative_humidity_min = bad);
            check_input_files = false)
    end
    clear = with_forcing(production; cloud_scheme = :none, cloud_relative_humidity_min = 0.7)
    @test_throws ArgumentError R.validate(clear; check_input_files = false)

    # The radiation adapter builds the cloud law with the requested threshold, and the default is unchanged.
    grid = R.SpeedyWeather.SpectralGrid(; NF = Float32, trunc = 7, nlayers = 8)
    tuned = R.RRTMGPRadiation(grid; cloud_scheme = :diagnostic, cloud_relative_humidity_min = 0.7)
    @test tuned.clouds.relative_humidity_threshold_min == 0.7f0
    released = R.RRTMGPRadiation(grid; cloud_scheme = :diagnostic)
    @test released.clouds.relative_humidity_threshold_min == 0.3f0
    @test released.clouds.relative_humidity_threshold_min ==
          R.SpeedyWeather.DiagnosticClouds(grid).relative_humidity_threshold_min
    @test_throws ArgumentError R.RRTMGPRadiation(grid; cloud_scheme = :diagnostic, cloud_relative_humidity_min = 1.0)
    @test_throws ArgumentError R.RRTMGPRadiation(grid; cloud_scheme = :diagnostic,
        cloud_relative_humidity_min = 1.0 - 1e-9)

    # Provenance records the requested and the applied threshold.
    provenance = R._forcing_provenance(tuned, diag)
    @test provenance.cloud_relative_humidity_min == 0.7
    @test provenance.applied_cloud_relative_humidity_min == Float64(0.7f0)
    @test R._forcing_provenance(released, production).applied_cloud_relative_humidity_min == Float64(0.3f0)
    # A clear-sky radiation object (NoClouds) has no cloud law: the applied value is NaN, never an error,
    # even when the config names an all-sky scheme (as in rrtmgp_trace_gases.jl).
    clear_radiation = R.RRTMGPRadiation(grid)
    @test clear_radiation.clouds isa R.SpeedyWeather.NoClouds
    @test isnan(R._forcing_provenance(clear_radiation, production).applied_cloud_relative_humidity_min)
end
println("CLOUD_RELATIVE_HUMIDITY_MIN_UNIT_PASS")
