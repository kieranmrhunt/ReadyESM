const R = ReadyESM

function _reference_constants_stage0b(temperature, salinity)
    temperature_100 = temperature / 100
    root_salinity = sqrt(salinity)
    log_k0 = -60.2409 + 93.4517 / temperature_100 +
             23.3585 * log(temperature_100) + salinity * (
                 0.023517 - 0.023656 * temperature_100 +
                 0.0047036 * temperature_100^2
             )
    pk1 = 3633.86 / temperature - 61.2172 + 9.6777 * log(temperature) -
          0.011555 * salinity + 0.0001152 * salinity^2
    pk2 = 471.78 / temperature + 25.9290 - 3.16967 * log(temperature) -
          0.01781 * salinity + 0.0001122 * salinity^2
    log_kb = (
        -8966.90 - 2890.53 * root_salinity - 77.942 * salinity +
        1.728 * salinity * root_salinity - 0.0996 * salinity^2
    ) / temperature + 148.0248 + 137.1942 * root_salinity +
        1.62142 * salinity - (
            24.4344 + 25.085 * root_salinity + 0.2474 * salinity
        ) * log(temperature) + 0.053105 * root_salinity * temperature
    log_kw = 148.9652 - 13847.26 / temperature - 23.6521 * log(temperature) +
             (-5.977 + 118.67 / temperature + 1.0495 * log(temperature)) *
             root_salinity - 0.01615 * salinity
    return (
        K0 = exp(log_k0),
        K1 = 10.0^(-pk1),
        K2 = 10.0^(-pk2),
        KB = exp(log_kb),
        Kw = exp(log_kw),
        total_boron = 4.157e-4 * salinity / 35,
    )
end

function _reference_pco2_stage0b(dic_umol_kg, alkalinity_umol_kg, temperature, salinity)
    constants = _reference_constants_stage0b(temperature, salinity)
    dic = dic_umol_kg * 1e-6
    alkalinity = alkalinity_umol_kg * 1e-6
    residual(hydrogen) = dic * (
        hydrogen * constants.K1 + 2 * constants.K1 * constants.K2
    ) / (
        hydrogen^2 + hydrogen * constants.K1 + constants.K1 * constants.K2
    ) + constants.total_boron * constants.KB / (constants.KB + hydrogen) +
        constants.Kw / hydrogen - hydrogen - alkalinity
    low = 1e-10
    high = 1e-6
    for _ in 1:200
        midpoint = sqrt(low * high)
        if residual(midpoint) > 0
            low = midpoint
        else
            high = midpoint
        end
        high / low < 1 + 1e-14 && break
    end
    hydrogen = sqrt(low * high)
    dissolved_co2 = dic * hydrogen^2 / (
        hydrogen^2 + hydrogen * constants.K1 + constants.K1 * constants.K2
    )
    return dissolved_co2 / constants.K0 * 1e6
end

function _reference_hydrogen_stage0b(
    dic_umol_kg,
    alkalinity_umol_kg,
    temperature,
    salinity,
)
    constants = _reference_constants_stage0b(temperature, salinity)
    dic = dic_umol_kg * 1e-6
    alkalinity = alkalinity_umol_kg * 1e-6
    residual(hydrogen) = dic * (
        hydrogen * constants.K1 + 2 * constants.K1 * constants.K2
    ) / (
        hydrogen^2 + hydrogen * constants.K1 + constants.K1 * constants.K2
    ) + constants.total_boron * constants.KB / (constants.KB + hydrogen) +
        constants.Kw / hydrogen - hydrogen - alkalinity
    low = 1e-12
    high = 1e-4
    for _ in 1:200
        midpoint = sqrt(low * high)
        if residual(midpoint) > 0
            low = midpoint
        else
            high = midpoint
        end
    end
    return sqrt(low * high)
end

function _reference_water_vapour_pressure_stage0b(temperature, salinity)
    return exp(
        24.4543 - 67.4509 * (100 / temperature) -
        4.8489 * log(temperature / 100) - 0.000544 * salinity,
    )
end

function _reference_flux_stage0b(
    dic,
    alkalinity,
    temperature,
    salinity,
    wind,
    ice,
    atmosphere,
    surface_pressure,
    density,
)
    pco2 = _reference_pco2_stage0b(dic, alkalinity, temperature, salinity)
    constants = _reference_constants_stage0b(temperature, salinity)
    temperature_celsius = temperature - 273.15
    schmidt = 2116.8 - 136.25 * temperature_celsius +
              4.7353 * temperature_celsius^2 -
              0.092307 * temperature_celsius^3 +
              0.0007555 * temperature_celsius^4
    transfer_velocity = 0.251 * wind^2 * sqrt(660 / schmidt) * 0.01 / 3600
    dry_air_pressure = surface_pressure / 101325 -
                       _reference_water_vapour_pressure_stage0b(temperature, salinity)
    air_pco2 = atmosphere * dry_air_pressure
    return transfer_velocity * density * constants.K0 * (1 - ice) *
           (air_pco2 - pco2) * 1e-6
end

@testset "BS10 Stage 0b device-safe component" begin
    state_grid = Tuple{Float64, Float64, Float64, Float64}[]
    for temperature in (271.15, 283.15, 298.15, 313.15),
        salinity in (5.0, 20.0, 35.0, 42.0),
        alkalinity in (1900.0, 2300.0, 2700.0),
        dic_fraction in (0.78, 0.88)
        push!(state_grid, (
            dic_fraction * alkalinity,
            alkalinity,
            temperature,
            salinity,
        ))
    end

    maximum_float64_error = 0.0
    maximum_float32_error = 0.0
    for (dic, alkalinity, temperature, salinity) in state_grid
        reference = R._carbon_pco2_uatm(dic, alkalinity, temperature, salinity)
        pco2_64, valid_64 = R._carbon_device_pco2_uatm(
            dic,
            alkalinity,
            temperature,
            salinity,
        )
        pco2_32, valid_32 = R._carbon_device_pco2_uatm(
            Float32(dic),
            Float32(alkalinity),
            Float32(temperature),
            Float32(salinity),
        )
        @test valid_64
        @test valid_32
        error_64 = abs(pco2_64 - reference) / reference
        error_32 = abs(Float64(pco2_32) - reference) / reference
        maximum_float64_error = max(maximum_float64_error, error_64)
        maximum_float32_error = max(maximum_float32_error, error_32)
    end
    @test maximum_float64_error <= 2e-12
    @test maximum_float32_error <= 8e-5

    for invalid in (
        (NaN32, 2300f0, 298.15f0, 35f0),
        (2050f0, 2300f0, 270f0, 35f0),
        (2050f0, 2300f0, 298.15f0, 4f0),
        (5000f0, 100f0, 298.15f0, 35f0),
    )
        pco2, valid = R._carbon_device_pco2_uatm(invalid...)
        @test !valid
        @test iszero(pco2)
    end

    one_cell = (
        dic = 2050.0,
        alkalinity = 2300.0,
        temperature = 291.15,
        salinity = 35.0,
        wind = 7.0,
        ice = 0.25,
        atmosphere = 420.0,
        surface_pressure = 101325.0,
        density = 1025.0,
    )
    reference_flux = _reference_flux_stage0b(
        one_cell.dic,
        one_cell.alkalinity,
        one_cell.temperature,
        one_cell.salinity,
        one_cell.wind,
        one_cell.ice,
        one_cell.atmosphere,
        one_cell.surface_pressure,
        one_cell.density,
    )
    flux64 = R._carbon_air_sea_flux_mol_m2_s(
        one_cell.dic,
        one_cell.alkalinity,
        one_cell.temperature,
        one_cell.salinity,
        one_cell.wind,
        one_cell.ice,
        one_cell.atmosphere,
        one_cell.surface_pressure,
        one_cell.density,
    )
    flux32 = R._carbon_air_sea_flux_mol_m2_s(
        Float32(one_cell.dic),
        Float32(one_cell.alkalinity),
        Float32(one_cell.temperature),
        Float32(one_cell.salinity),
        Float32(one_cell.wind),
        Float32(one_cell.ice),
        Float32(one_cell.atmosphere),
        Float32(one_cell.surface_pressure),
        Float32(one_cell.density),
    )
    @test flux64.valid
    @test flux32.valid
    @test isapprox(flux64.flux, reference_flux; rtol = 2e-12, atol = 0)
    @test isapprox(Float64(flux32.flux), reference_flux; rtol = 8e-5, atol = 0)
    reference_air_pco2 = one_cell.atmosphere * (
        one_cell.surface_pressure / 101325 -
        _reference_water_vapour_pressure_stage0b(
            one_cell.temperature,
            one_cell.salinity,
        )
    )
    @test isapprox(flux64.air_pco2, reference_air_pco2; rtol = 2e-15, atol = 0)
    @test isapprox(Float64(flux32.air_pco2), reference_air_pco2; rtol = 8e-5, atol = 0)

    boundary_states = [
        (temperature = prevfloat(271.15f0), salinity = 35f0, valid = false),
        (temperature = 271.15f0, salinity = 35f0, valid = true),
        (temperature = nextfloat(271.15f0), salinity = 35f0, valid = true),
        (temperature = prevfloat(313.15f0), salinity = 35f0, valid = true),
        (temperature = 313.15f0, salinity = 35f0, valid = true),
        (temperature = nextfloat(313.15f0), salinity = 35f0, valid = false),
        (temperature = 298.15f0, salinity = prevfloat(5f0), valid = false),
        (temperature = 298.15f0, salinity = 5f0, valid = true),
        (temperature = 298.15f0, salinity = nextfloat(5f0), valid = true),
        (temperature = 298.15f0, salinity = prevfloat(42f0), valid = true),
        (temperature = 298.15f0, salinity = 42f0, valid = true),
        (temperature = 298.15f0, salinity = nextfloat(42f0), valid = false),
    ]
    for state in boundary_states
        pco2, pco2_valid = R._carbon_device_pco2_uatm(
            2050f0,
            2300f0,
            state.temperature,
            state.salinity,
        )
        flux = R._carbon_air_sea_flux_mol_m2_s(
            2050f0,
            2300f0,
            state.temperature,
            state.salinity,
            7f0,
            0.25f0,
            420f0,
            101325f0,
            1025f0,
        )
        @test pco2_valid == state.valid
        @test flux.valid == pco2_valid
        if state.valid
            @test isfinite(pco2)
            @test isfinite(flux.flux)
            @test isfinite(flux.pco2)
            @test isfinite(flux.air_pco2)
        else
            @test iszero(pco2)
            @test iszero(flux.flux)
        end
    end

    ph_edge_states = [
        (alkalinity = 4932.7354f0, hydrogen_edge = 1e-10, valid = false),
        (alkalinity = 1205.4783f0, hydrogen_edge = 1e-6, valid = false),
    ]
    for state in ph_edge_states
        reference_hydrogen = _reference_hydrogen_stage0b(
            Float64(2050f0),
            Float64(state.alkalinity),
            Float64(298.15f0),
            Float64(35f0),
        )
        @test abs(reference_hydrogen / state.hydrogen_edge - 1) <= 1e-6
        _, pco2_valid = R._carbon_device_pco2_uatm(
            2050f0,
            state.alkalinity,
            298.15f0,
            35f0,
        )
        flux = R._carbon_air_sea_flux_mol_m2_s(
            2050f0,
            state.alkalinity,
            298.15f0,
            35f0,
            7f0,
            0.25f0,
            420f0,
            101325f0,
            1025f0,
        )
        @test pco2_valid == state.valid
        @test flux.valid == pco2_valid
    end

    ph_guard_states = [
        (alkalinity = 4932.714f0, hydrogen_edge = 1f-10,
         guard_multiple = 0.5f0, valid = false),
        (alkalinity = 4932.6396f0, hydrogen_edge = 1f-10,
         guard_multiple = 2f0, valid = true),
        (alkalinity = 1205.4883f0, hydrogen_edge = 1f-6,
         guard_multiple = 0.5f0, valid = false),
        (alkalinity = 1205.5063f0, hydrogen_edge = 1f-6,
         guard_multiple = 2f0, valid = true),
    ]
    for state in ph_guard_states
        constants = R._carbon_device_constants(298.15f0, 35f0)
        alkalinity_mol_kg = state.alkalinity * 1f-6
        endpoint_residual = R._carbon_device_alkalinity_residual(
            state.hydrogen_edge,
            2050f0 * 1f-6,
            alkalinity_mol_kg,
            constants,
        )
        relative_residual = abs(endpoint_residual) / alkalinity_mol_kg
        @test isapprox(
            relative_residual,
            state.guard_multiple *
                Float32(R._CARBON_STAGE0B_BRACKET_RESIDUAL_GUARD);
            rtol = 0.05f0,
            atol = 0f0,
        )
        pco2, pco2_valid = R._carbon_device_pco2_uatm(
            2050f0,
            state.alkalinity,
            298.15f0,
            35f0,
        )
        flux = R._carbon_air_sea_flux_mol_m2_s(
            2050f0,
            state.alkalinity,
            298.15f0,
            35f0,
            7f0,
            0.25f0,
            420f0,
            101325f0,
            1025f0,
        )
        @test pco2_valid == state.valid
        @test flux.valid == pco2_valid
        if state.valid
            reference_pco2 = _reference_pco2_stage0b(
                Float64(2050f0),
                Float64(state.alkalinity),
                Float64(298.15f0),
                Float64(35f0),
            )
            @test abs(Float64(pco2) - reference_pco2) / reference_pco2 <= 8e-5
        end
    end

    zero_wind = R._carbon_air_sea_flux_mol_m2_s(
        one_cell.dic,
        one_cell.alkalinity,
        one_cell.temperature,
        one_cell.salinity,
        0.0,
        one_cell.ice,
        one_cell.atmosphere,
        one_cell.surface_pressure,
        one_cell.density,
    )
    full_ice = R._carbon_air_sea_flux_mol_m2_s(
        one_cell.dic,
        one_cell.alkalinity,
        one_cell.temperature,
        one_cell.salinity,
        one_cell.wind,
        1.0,
        one_cell.atmosphere,
        one_cell.surface_pressure,
        one_cell.density,
    )
    dry_air_pressure = one_cell.surface_pressure / 101325 -
                       _reference_water_vapour_pressure_stage0b(
        one_cell.temperature,
        one_cell.salinity,
    )
    equilibrium = R._carbon_air_sea_flux_mol_m2_s(
        one_cell.dic,
        one_cell.alkalinity,
        one_cell.temperature,
        one_cell.salinity,
        one_cell.wind,
        one_cell.ice,
        flux64.pco2 / dry_air_pressure,
        one_cell.surface_pressure,
        one_cell.density,
    )
    @test iszero(zero_wind.flux)
    @test iszero(full_ice.flux)
    @test isapprox(equilibrium.air_pco2, equilibrium.pco2; rtol = 8eps(Float64), atol = 0)
    @test abs(equilibrium.flux) <= 8eps(Float64) * abs(flux64.flux)
    @test flux64.flux * (flux64.air_pco2 - flux64.pco2) > 0

    invalid_dry_air = R._carbon_air_sea_flux_mol_m2_s(
        one_cell.dic,
        one_cell.alkalinity,
        one_cell.temperature,
        one_cell.salinity,
        one_cell.wind,
        one_cell.ice,
        one_cell.atmosphere,
        1000.0,
        one_cell.density,
    )
    @test !invalid_dry_air.valid
    @test iszero(invalid_dry_air.flux)

    cell = R._carbon_apply_surface_cell_flux(
        one_cell.dic,
        -1.0,
        3600.0,
        1e10,
        1.0,
        one_cell.density,
    )
    @test cell.valid
    @test cell.limited
    @test isapprox(cell.dic_umol_kg, 0.1 * one_cell.dic; rtol = 32eps(Float64))

    values = [1e-12, 1e12, -1e12, 3.0, -2.0]
    @test R._carbon_fixed_order_sum(values) == 1.000000000001

    for mode in (:emissions_driven, :prescribed)
        initial = R._Stage0bCarbonLedger(mode, 892.08, 38_000.0)
        state = initial
        for step in 1:48
            uptake = isodd(step) ? 0.03 : -0.01
            emissions = mode === :emissions_driven ? 0.02 : 0.0
            state = R._carbon_ledger_step(state, uptake, emissions)
        end
        residual = R._carbon_ledger_invariant(state) -
                   R._carbon_ledger_invariant(initial)
        tolerance = 64eps(Float64) * (
            state.atmospheric_carbon_pgc + state.ocean_carbon_pgc
        )
        @test abs(residual) <= tolerance
        @test state.cumulative_air_sea_uptake_pgc > 0
        if mode === :prescribed
            @test state.atmospheric_carbon_pgc == initial.atmospheric_carbon_pgc
            @test state.cumulative_implied_emissions_pgc ==
                  state.cumulative_air_sea_uptake_pgc
        else
            @test isapprox(
                state.cumulative_emissions_pgc,
                48 * 0.02;
                rtol = 8eps(Float64),
                atol = 0,
            )
            @test iszero(state.cumulative_implied_emissions_pgc)
        end

        checkpoint = (
            mode = state.mode,
            atmospheric_carbon_pgc = state.atmospheric_carbon_pgc,
            ocean_carbon_pgc = state.ocean_carbon_pgc,
            cumulative_emissions_pgc = state.cumulative_emissions_pgc,
            cumulative_air_sea_uptake_pgc = state.cumulative_air_sea_uptake_pgc,
            cumulative_implied_emissions_pgc = state.cumulative_implied_emissions_pgc,
        )
        restored = R._Stage0bCarbonLedger(
            checkpoint.mode,
            checkpoint.atmospheric_carbon_pgc,
            checkpoint.ocean_carbon_pgc,
            checkpoint.cumulative_emissions_pgc,
            checkpoint.cumulative_air_sea_uptake_pgc,
            checkpoint.cumulative_implied_emissions_pgc,
        )
        continued = R._carbon_ledger_step(restored, 0.015, 0.0)
        direct = R._carbon_ledger_step(state, 0.015, 0.0)
        @test isequal(continued, direct)
    end
end
