const R = ReadyESM
const O = R.Oceananigans

function _valid_prescribed_carbon_config(; device = :cpu)
    return R.ExperimentConfig(
        device = device,
        ocean_model = :oceananigans,
        ocean_dynamics = true,
        forcing = R.ForcingConfig(radiation = :rrtmgp_clear_sky),
        carbon = R.CarbonConfig(mode = :prescribed),
    )
end

function _stage0c_test_controller(config)
    ppm = config.carbon.atmospheric_co2_ppm
    return R.CarbonController(
        :prescribed,
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
end

@testset "BS10 Stage 0c configuration and storage" begin
    cpu = _valid_prescribed_carbon_config()
    gpu = _valid_prescribed_carbon_config(device = :gpu)
    @test R.validate(cpu; check_input_files = false) === cpu
    @test R.validate(gpu; check_input_files = false) === gpu
    @test R._carbon_ocean_kwargs(cpu) == (
        tracers = (:T, :S, :DIC, :TA),
    )
    @test R._carbon_storage_offset(:prescribed, Val(:DIC), Float32) ==
          2000f0
    @test R._carbon_storage_offset(:prescribed, Val(:TA), Float32) ==
          2300f0
    @test iszero(R._carbon_storage_offset(
        :tracer_only,
        Val(:DIC),
        Float32,
    ))
    @test cpu.carbon.mode == :prescribed
end

@testset "BS10 Stage 0c absolute freshwater carry" begin
    grid = O.LatitudeLongitudeGrid(
        O.CPU(),
        Float64;
        size = (4, 4, 2),
        longitude = (0, 360),
        latitude = (-60, 60),
        z = (-100, 0),
        halo = (1, 1, 1),
    )
    carrying_flux = O.Field{
        O.Center,
        O.Center,
        Nothing,
    }(grid)
    stored_dic = O.CenterField(grid)
    stored_ta = O.CenterField(grid)
    O.set!(carrying_flux, 2e-5)
    O.set!(stored_dic, 125.0)
    O.set!(stored_ta, -25.0)
    dic_exchange = R._CarbonFreshwaterExchange{
        :DIC,
        typeof(carrying_flux),
        Float64,
    }(carrying_flux, 2000.0)
    ta_exchange = R._CarbonFreshwaterExchange{
        :TA,
        typeof(carrying_flux),
        Float64,
    }(carrying_flux, 2300.0)
    fields = (; DIC = stored_dic, TA = stored_ta)
    @test dic_exchange(2, 2, grid, nothing, fields) == 2125.0 * 2e-5
    @test ta_exchange(2, 2, grid, nothing, fields) == 2275.0 * 2e-5
end

@testset "BS10 carbon boundary kernel representation" begin
    host_exchange = R._CarbonFreshwaterExchange{:DIC, Tuple{Float32}, Float32}(
        (2f-5,), 2000f0,
    )
    @test !isbitstype(typeof(host_exchange))
    kernel_exchange = R.Adapt.adapt(nothing, host_exchange)
    @test kernel_exchange isa R._CarbonFreshwaterKernelExchange{:DIC}
    @test isbitstype(typeof(kernel_exchange))
    @test kernel_exchange.carrying_flux == host_exchange.carrying_flux
    @test kernel_exchange.storage_offset == host_exchange.storage_offset
    @test R.Adapt.adapt(nothing, kernel_exchange) === kernel_exchange
    host_exchange.carrying_flux = (3f-5,)
    @test kernel_exchange.carrying_flux == (2f-5,)
    @test R.Adapt.adapt(nothing, host_exchange).carrying_flux == (3f-5,)
    array_exchange = R._CarbonFreshwaterExchange{:TA, Vector{Float32}, Float32}(
        [2f-5], 2300f0,
    )
    @test R.Adapt.adapt(nothing, array_exchange) isa R._CarbonFreshwaterExchange{:TA}
end

@testset "BS10 Stage 0c prescribed ledger order" begin
    config = _valid_prescribed_carbon_config()
    controller = _stage0c_test_controller(config)
    initial_invariant = controller.atmospheric_carbon_pgc +
        controller.ocean_carbon_pgc

    R._advance_carbon_ledger!(controller, 0, 0.0, 1e-8)
    @test controller.ledger_iteration == 0
    @test controller.previous_air_sea_uptake_rate_pgc_s == 1e-8
    R._advance_carbon_ledger!(controller, 0, 0.0, 9e-8)
    @test controller.previous_air_sea_uptake_rate_pgc_s == 1e-8

    R._advance_carbon_ledger!(controller, 1, 2400.0, 2e-8)
    expected_uptake = 1e-8 * 2400.0
    @test controller.cumulative_air_sea_uptake_pgc == expected_uptake
    @test controller.cumulative_implied_emissions_pgc == expected_uptake
    @test controller.ocean_carbon_pgc == 38_000.0 + expected_uptake
    @test controller.previous_air_sea_uptake_rate_pgc_s == 2e-8
    @test isapprox(
        controller.atmospheric_carbon_pgc + controller.ocean_carbon_pgc -
        controller.cumulative_implied_emissions_pgc,
        initial_invariant;
        rtol = 64eps(Float64),
        atol = 0,
    )
    @test_throws ErrorException R._advance_carbon_ledger!(
        controller,
        3,
        7200.0,
        3e-8,
    )

    state = O.prognostic_state(controller)
    restored = _stage0c_test_controller(config)
    O.restore_prognostic_state!(restored, state)
    @test O.prognostic_state(restored) == state
end

@testset "BS10 Stage 0c Float64 transport offset identity" begin
    offset = R._CARBON_DIC_STORAGE_OFFSET_UMOL_KG
    grid = O.RectilinearGrid(
        O.CPU(),
        Float64;
        size = 16,
        x = (0, 1),
        topology = (O.Periodic, O.Flat, O.Flat),
        halo = 4,
    )
    function transport_model()
        return O.NonhydrostaticModel(
            grid;
            tracers = :c,
            advection = O.WENO(order = 5),
            buoyancy = nothing,
            coriolis = nothing,
        )
    end

    uniform = transport_model()
    absolute = transport_model()
    anomaly = transport_model()
    profile(x) = 2050.0 + 10sin(2pi * x)
    anomaly_profile(x) = profile(x) - offset
    O.set!(uniform, u = 1.0, c = 2050.0)
    O.set!(absolute, u = 1.0, c = profile)
    O.set!(anomaly, u = 1.0, c = anomaly_profile)
    O.time_step!(uniform, 1e-4)
    O.time_step!(absolute, 1e-4)
    O.time_step!(anomaly, 1e-4)

    uniform_values = Array(O.interior(uniform.tracers.c))
    absolute_values = Array(O.interior(absolute.tracers.c))
    reconstructed_values = Array(O.interior(anomaly.tracers.c)) .+ offset
    @test all(==(2050.0), uniform_values)
    @test isapprox(
        reconstructed_values,
        absolute_values;
        rtol = 2e-14,
        atol = 2e-12,
    )
end
