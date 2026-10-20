const R = ReadyESM

function _valid_tracer_only_config(; carbon = R.CarbonConfig(
    mode = :tracer_only,
), forcing = R.ForcingConfig(
    radiation = :rrtmgp_clear_sky,
), kwargs...)
    settings = merge((
        ; device = :cpu,
        ocean_model = :oceananigans,
        forcing,
        carbon,
    ), (; kwargs...))
    return R.ExperimentConfig(; settings...)
end

@testset "BS10 Stage 0a configuration" begin
    default = R.ExperimentConfig()
    @test default.carbon == R.CarbonConfig()
    @test isempty(keys(R._carbon_ocean_kwargs(default)))

    active = _valid_tracer_only_config()
    @test R.validate(active; check_input_files = false) === active
    @test R._carbon_ocean_kwargs(active) == (
        tracers = (:T, :S, :DIC, :TA),
    )

    @test_throws ArgumentError R.validate(R.ExperimentConfig(
        carbon = R.CarbonConfig(mode = :interactive),
    ); check_input_files = false)
    @test_throws ArgumentError R.validate(_valid_tracer_only_config(
        device = :gpu,
    ); check_input_files = false)
    @test_throws ArgumentError R.validate(_valid_tracer_only_config(
        ocean_model = :slab,
    ); check_input_files = false)
    @test_throws ArgumentError R.validate(_valid_tracer_only_config(
        forcing = R.ForcingConfig(radiation = :speedy_simplified),
    ); check_input_files = false)
    @test_throws ArgumentError R.validate(_valid_tracer_only_config(
        forcing = R.ForcingConfig(
            radiation = :rrtmgp_clear_sky,
            co2_ppm = 560,
        ),
    ); check_input_files = false)
    @test_throws ArgumentError R.validate(_valid_tracer_only_config(
        carbon = R.CarbonConfig(
            mode = :tracer_only,
            alkalinity_at_salinity_35_umol_kg = 0,
        ),
    ); check_input_files = false)
    @test_throws ArgumentError R.validate(_valid_tracer_only_config(
        carbon = R.CarbonConfig(
            mode = :tracer_only,
            chemistry_salinity_floor = 0,
        ),
    ); check_input_files = false)
    @test_throws ArgumentError R.validate(_valid_tracer_only_config(
        ocean_river_mouth_mixing = :localized_estuary_v3,
    ); check_input_files = false)
    @test_throws ArgumentError R.validate(R.ExperimentConfig(
        carbon = R.CarbonConfig(ocean_equilibrium_co2_ppm = 300),
    ); check_input_files = false)

    mktempdir() do directory
        path = joinpath(directory, "unknown_carbon_key.yml")
        write(path, "carbon:\n  mode: off\n  unknown_key: 1\n")
        @test_throws ArgumentError R.load_config(
            path;
            check_input_files = false,
        )
    end
end

@testset "BS10 Stage 0a carbonate chemistry" begin
    constants = R._carbon_carbonate_constants(298.15, 35.0)
    @test isapprox(log(constants.K0), -3.5617; atol = 5e-5, rtol = 0)
    @test isapprox(log10(constants.K1), -5.8472; atol = 5e-5, rtol = 0)
    @test isapprox(log10(constants.K2), -8.9660; atol = 5e-5, rtol = 0)
    @test isapprox(log(constants.KB), -19.7964; atol = 5e-5, rtol = 0)
    @test isapprox(log(constants.Kw), -30.434; atol = 5e-4, rtol = 0)
    @test constants.total_boron == 4.157e-4

    dic = R._carbon_equilibrium_dic_umol_kg(
        280.0,
        2300.0,
        298.15,
        35.0,
    )
    @test isapprox(
        R._carbon_pco2_uatm(dic, 2300.0, 298.15, 35.0),
        280.0;
        rtol = 1e-12,
    )

    original_dic = 2050.0
    pco2 = R._carbon_pco2_uatm(
        original_dic,
        2300.0,
        298.15,
        35.0,
    )
    recovered_dic = R._carbon_equilibrium_dic_umol_kg(
        pco2,
        2300.0,
        298.15,
        35.0,
    )
    @test isapprox(recovered_dic, original_dic; rtol = 1e-9)
end

@testset "BS10 Stage 0a controller state" begin
    config = _valid_tracer_only_config()
    @test R._carbon_provenance(config).schema ==
          "readiesm_carbon_stage0a_v2"
    ppm = config.carbon.atmospheric_co2_ppm
    controller = R.CarbonController(
        :tracer_only,
        ppm * R._CARBON_PGC_PER_PPM,
        ppm,
        38_000.0,
        0.0,
        0.0,
        0.0,
        0.0,
        NaN,
        -1,
        2400.0,
        R._carbon_provenance(config),
    )
    state = R.Oceananigans.prognostic_state(controller)
    controller.atmospheric_carbon_pgc = 1.0
    controller.atmospheric_co2_ppm = 1.0
    R.Oceananigans.restore_prognostic_state!(controller, state)
    @test controller.atmospheric_carbon_pgc ==
          controller.atmospheric_co2_ppm * R._CARBON_PGC_PER_PPM

    active_model = Ref(:active_carbon_model)
    R._register_carbon_controller!(active_model, controller)
    @test isnothing(R._validate_carbon_restart_compatibility(
        active_model,
        state,
    ))
    @test_throws ErrorException R._validate_carbon_restart_compatibility(
        active_model,
        nothing,
    )
    @test_throws ErrorException R._validate_carbon_restart_compatibility(
        Ref(:off_carbon_model),
        state,
    )

    changed = _valid_tracer_only_config(carbon = R.CarbonConfig(
        mode = :tracer_only,
        ocean_equilibrium_co2_ppm = 300.0,
    ))
    changed_controller = R.CarbonController(
        :tracer_only,
        ppm * R._CARBON_PGC_PER_PPM,
        ppm,
        38_000.0,
        0.0,
        0.0,
        0.0,
        0.0,
        NaN,
        -1,
        2400.0,
        R._carbon_provenance(changed),
    )
    @test_throws ErrorException R.Oceananigans.restore_prognostic_state!(
        changed_controller,
        state,
    )

    legacy_provenance = merge(
        controller.provenance,
        (; schema = "readiesm_carbon_stage0a_v1"),
    )
    legacy_v1_state = (
        mode = "tracer_only",
        provenance = legacy_provenance,
    )
    legacy_error = try
        R.Oceananigans.restore_prognostic_state!(controller, legacy_v1_state)
        nothing
    catch error
        error
    end
    @test legacy_error isa ErrorException
    @test sprint(showerror, legacy_error) ==
          "checkpoint carbon provenance differs from the reconstructed controller"
end
