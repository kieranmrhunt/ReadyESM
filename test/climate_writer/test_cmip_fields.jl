# Unit tests for CMIPFields (Claude, 2026-10-29). Pure Julia: Dates and Test only, no model.
# usage: julia --startup-file=no test_cmip_fields.jl
using Test
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CMIP_DIR = joinpath(PROJECT_ROOT, "scripts", "climate_writer", "cmip")
include(joinpath(CMIP_DIR, "cmip_fields.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "climate_writer", "profiles", "profiles.jl"))
using .CMIPFields: max_random_cover, grid_mean_path, radiative_surface_temperature, surface_samples,
                   surface_stress, selected_drag,
                   profile_samples, profile_name, LEVELS_HPA, SPECS

@testset "CMIP fields" begin
    @testset "maximum-random overlap" begin
        @test max_random_cover(Float64[]) == 0
        @test max_random_cover([0.0, 0.0, 0.0]) == 0
        @test max_random_cover([0.0, 0.4, 0.0]) ≈ 0.4
        @test max_random_cover([0.3, 0.6, 0.2]) ≈ 0.6          # one block: maximum overlap
        @test max_random_cover([0.3, 0.0, 0.6]) ≈ 1 - 0.7 * 0.4  # two blocks: random overlap
        @test max_random_cover([0.6, 0.3, 0.0, 0.5]) ≈ 1 - 0.4 * 0.5
        @test max_random_cover([0.2, 1.0, 0.1]) == 1
        @test max_random_cover(reverse([0.3, 0.0, 0.6, 0.6, 0.1])) ≈ max_random_cover([0.3, 0.0, 0.6, 0.6, 0.1])
    end
    @testset "paths and surface temperature" begin
        @test grid_mean_path([0.5, 0.0, 1.0], [40.0, 99.0, 10.0]) ≈ 0.030
        @test radiative_surface_temperature(0.0, 280.0, 250.0) ≈ 280.0
        @test radiative_surface_temperature(1.0, 280.0, 250.0) ≈ 250.0
        @test radiative_surface_temperature(0.5, 280.0, 250.0)^4 ≈ (280.0^4 + 250.0^4) / 2
        cf = [0.0 0.5; 0.4 0.5; 0.4 0.0]                       # 3 layers, 2 columns
        lw = [0.0 30.0; 20.0 30.0; 20.0 0.0]
        iw = [0.0 10.0; 5.0 0.0; 5.0 0.0]
        s = surface_samples([270.0, 290.0], [0.0, 1.0], [271.0, 300.0], [260.0, 288.0], cf, lw, iw,
                            [0.1, -0.2], [0.0, 0.3])
        @test s["cmip_tas"] == [270.0, 290.0]
        @test s["cmip_ts"] ≈ [271.0, 288.0]
        @test s["cmip_clt"] ≈ [40.0, 50.0]
        @test s["cmip_clivi"] ≈ [0.004, 0.005]
        @test s["cmip_clwvi"] ≈ [0.020, 0.035]
        @test s["cmip_tauu"] == [0.1, -0.2] && s["cmip_tauv"] == [0.0, 0.3]
        @test_throws ErrorException surface_samples([270.0, NaN], [0.0, 1.0], [271.0, 300.0], [260.0, 288.0], cf, lw, iw,
                                                    [0.1, -0.2], [0.0, 0.3])
        @test_throws DimensionMismatch surface_samples([270.0], [0.0, 1.0], [271.0, 300.0], [260.0, 288.0], cf, lw, iw,
                                                       [0.1, -0.2], [0.0, 0.3])
    end
    @testset "surface stress" begin
        # Matches SpeedyWeather's flux_u_upward = -rho * drag * V0 * f * u with the opposite (downward) sign.
        rho, drag, V0, f, u, v = [1.2, 1.0], [1.5e-3, 2.0e-3], [8.0, 5.0], 0.95, [6.0, -3.0], [-2.0, 4.0]
        tu, tv = surface_stress(rho, drag, V0, f, u, v)
        @test tu ≈ -(-rho .* drag .* V0 .* f .* u)
        @test tv ≈ rho .* drag .* V0 .* f .* v
        @test tu[1] > 0 && tu[2] < 0                       # eastward wind gives eastward stress on the surface
        @test selected_drag(true, [1e-3, 2e-3], [0.0, 1.0], 2.4e-3, 0.8e-3) == [1e-3, 2e-3]
        @test selected_drag(false, [1e-3, 2e-3], [0.0, 0.25], 2.4e-3, 0.8e-3) ≈ [0.8e-3, 0.25 * 2.4e-3 + 0.75 * 0.8e-3]
        @test_throws ErrorException surface_stress(rho, drag, V0, f, [NaN, 1.0], v)
    end
    @testset "pressure levels" begin
        nz = 27
        sigma = [(k - 0.5) / nz for k in 1:nz]
        dsigma = fill(1 / nz, nz)
        n = 3
        ps = [101325.0, 95000.0, 60000.0]
        T = [220.0 + 70 * sigma[k]^0.3 + 3i for i in 1:n, k in 1:nz]
        u = [10.0 + 20 * (1 - sigma[k]) - i for i in 1:n, k in 1:nz]
        v = [-2.0 + 5 * sigma[k] * i for i in 1:n, k in 1:nz]
        q = [1e-5 + 0.012 * sigma[k]^3 / i for i in 1:n, k in 1:nz]
        R, g = 287.05, 9.81
        T0 = 250.0
        phis = [0.0, 500.0 * g, 4000.0 * g]
        # Isothermal hydrostatic geopotential: exactly linear in log pressure.
        phi = [phis[i] + R * T0 * log(1 / sigma[k]) for i in 1:n, k in 1:nz]
        result = profile_samples(T, u, v, q, phi, ps, sigma, g)
        old = ClimateProfiles.samples(T, q, u, ps, sigma, dsigma, g)
        for (level, field, oldname) in ((200, "ta", "temperature_200hpa"), (500, "ta", "temperature_500hpa"),
                                        (850, "ta", "temperature_850hpa"), (200, "ua", "zonal_wind_200hpa"),
                                        (850, "ua", "zonal_wind_850hpa"))
            name = profile_name(field, level)
            @test result[name * "_contribution"] == old[oldname * "_contribution"]
            @test result[name * "_valid_fraction"] == old[oldname * "_valid_fraction"]
        end
        for level in LEVELS_HPA
            name = profile_name("zg", level)
            p = 100.0 * level
            for i in 1:n
                inside = sigma[1] * ps[i] <= p <= sigma[end] * ps[i]
                @test result[name * "_valid_fraction"][i] == (inside ? 1.0 : 0.0)
                if inside
                    @test isapprox(result[name * "_contribution"][i], phis[i] / g + R * T0 / g * log(ps[i] / p); rtol=1e-12)
                end
            end
        end
        @test result[profile_name("ta", 1000) * "_valid_fraction"] == [0.0, 0.0, 0.0]   # below the lowest full level
        @test result[profile_name("hus", 20) * "_valid_fraction"] == [1.0, 1.0, 1.0]
        @test_throws ErrorException profile_samples(fill(NaN, n, nz), u, v, q, phi, ps, sigma, g)
    end
    @testset "metadata" begin
        @test length(SPECS) == 7 + 5 * length(LEVELS_HPA) * 2
        registry = Dict{String,Tuple{Symbol,String,String}}()
        CMIPFields.register_metadata!(registry)
        @test_throws ErrorException CMIPFields.register_metadata!(registry)
        @test isempty(intersect(keys(SPECS), keys(ClimateProfiles.SPECS)))
    end
end
println("CMIP_FIELDS_TESTS_PASS")
