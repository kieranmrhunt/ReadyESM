module MonthlyUpperRadiativeCapTests
using ReadyESM, Test, Dates, SHA
const R = ReadyESM
const G = R.RRTMGP
const S = R.SpeedyWeather
const TABLE = joinpath(R.PROJECT_ROOT, "data", "upper_radiative_cap", "era5_1993_sampled_monthly_upper_cap.csv")

@testset "Sampled monthly upper radiative cap" begin
    @test bytes2hex(open(SHA.sha256, TABLE)) == R.UPPER_CAP_MONTHLY1993_SHA256
    data = R._monthly_upper_cap_data(TABLE)
    @test size(data) == (5, 12, 181)
    @test all(isfinite, data)

    for (year, month) in ((1993, 1), (1993, 12), (2024, 2), (2025, 2))
        @test R._upper_cap_month_bracket(R._upper_cap_month_centre(year, month)) ==
            (month, month == 12 ? 1 : month + 1, 0.0)
    end
    @test R._upper_cap_month_bracket(DateTime(2024, 1, 1))[1:2] == (12, 1)
    @test R._upper_cap_month_bracket(DateTime(2024, 12, 31))[1:2] == (12, 1)
    @test R._upper_cap_month_centre(2024, 2) == DateTime(2024, 2, 15, 12)
    @test R._upper_cap_month_centre(2025, 2) == DateTime(2025, 2, 15)
    for t in (DateTime(1993, 1, 1), DateTime(2024, 2, 29, 12), DateTime(2024, 12, 31, 23))
        left, right, f = R._upper_cap_month_bracket(t)
        @test 1 <= left <= 12 && 1 <= right <= 12 && 0 <= f <= 1
    end

    grid = S.SpectralGrid(; NF = Float32, trunc = 7, nlayers = 8)
    radiation = R.RRTMGPRadiation(grid; upper_radiative_cap_path = TABLE)
    cap = radiation.upper_cap
    @test cap isa R.MonthlyUpperRadiativeCap
    @test size(cap.temperatures) == (5, 12, grid.npoints)
    @test cap.sha256 == R.UPPER_CAP_MONTHLY1993_SHA256
    @test R._upper_cap_identity(cap).method == R.UPPER_CAP_MONTHLY_METHOD
    @test R._validate_upper_cap_restart(radiation, (; upper_radiative_cap = R._upper_cap_identity(cap))) === nothing
    @test_throws ArgumentError R._validate_upper_cap_restart(radiation, (; upper_radiative_cap = nothing))

    solver = radiation.solver
    latitude = Array(G.latitude(solver))
    col = argmin(abs.(latitude))
    surface_pressure = fill(1.0f5, grid.npoints)
    nmodel = grid.nlayers
    for (time, left, right, f) in ((R._upper_cap_month_centre(1993, 1), 1, 2, 0.0),
                                    (DateTime(1993, 2, 1), 1, 2, R._upper_cap_month_bracket(DateTime(1993, 2, 1))[3]),
                                    (R._upper_cap_month_centre(1993, 12), 12, 1, 0.0))
        @test R._fill_upper_radiative_cap!(cap, solver, surface_pressure, nothing, nmodel, time) === nothing
        temperature = Array(G.layer_temperature(solver))
        boundary = Array(G.level_temperature(solver))
        expected(field) = (1 - f) * cap.temperatures[field, left, col] + f * cap.temperatures[field, right, col]
        @test temperature[nmodel + 1, col] ≈ expected(1) atol = 5e-5
        @test temperature[nmodel + 2, col] ≈ expected(2) atol = 5e-5
        @test boundary[nmodel + 1, col] ≈ expected(3) atol = 5e-5
        @test boundary[nmodel + 2, col] ≈ expected(4) atol = 5e-5
        @test boundary[nmodel + 3, col] ≈ expected(5) atol = 5e-5
    end

    allsky = R.RRTMGPRadiation(grid; cloud_scheme = :diagnostic,
                              upper_radiative_cap_path = TABLE)
    fraction = G.cloud_fraction(allsky.solver)
    fill!(fraction, 0)
    fraction[1, col] = 0.25f0
    model = (; geometry = (; npoints = grid.npoints, nlayers = nmodel),
              longwave_radiation = allsky)
    profile = R._atmosphere_cloud_layer_diagnostics(model)
    @test size(profile.fraction) == (grid.npoints, nmodel)
    @test profile.fraction[col, end] == 0.25
    @test profile.column_fraction[col] == 0.25
    fraction[nmodel + 1, col] = 0.1f0
    @test_throws ArgumentError R._atmosphere_cloud_layer_diagnostics(model)

    uncapped = R.RRTMGPRadiation(grid; cloud_scheme = :diagnostic)
    uncapped_profile = R._atmosphere_cloud_layer_diagnostics((;
        geometry = model.geometry, longwave_radiation = uncapped))
    @test size(uncapped_profile.fraction) == (grid.npoints, nmodel)
    @test R._upper_cap_identity(uncapped.upper_cap) === nothing
end
end
