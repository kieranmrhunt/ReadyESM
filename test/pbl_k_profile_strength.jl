module PBLKProfileStrengthTests
# Strength-screen parameters of the K-profile (diffusivity_scale, mix_heat, mix_humidity): configuration, refusals,
# provenance, the exact default path, and the column behaviour on a compact coupled model with 27 layers (at 8 layers
# only one column in 288 mixes, which is too weak a test).
using ReadyESM, Test
const R = ReadyESM
const S = R.SpeedyWeather

fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
with(config; kwargs...) = R.ExperimentConfig(; merge(fields(config), (; kwargs...))...)

function compact_l27(output)
    base = R.load_config(joinpath(R.PROJECT_ROOT, "config/production.yml"); check_input_files = false)
    era5 = joinpath(R.PROJECT_ROOT, "data/era5_initial_state")
    return with(base; device = :cpu, truncation = 7, nlayers = 27, ocean_grid = :idealized,
        ocean_nlongitude = 32, ocean_nlatitude = 16, ocean_nlayers = 4, ocean_initial_conditions = :analytic,
        ecco_initial_conditions_directory = "", duration_days = 4 * 3600 / 86_400, analysis_start_day = 0.0,
        output_dir = output,
        era5_pressure_levels_path = joinpath(era5, "era5_pressure_levels_stratosphere_19930101T0000Z.nc"),
        era5_single_levels_path = joinpath(era5, "era5_single_levels_19930101T0000Z.nc"),
        era5_land_state_path = joinpath(era5, "era5_land_state_19930101T0000Z.nc"),
        atmosphere_vertical_diffusion = :holtslag_boville_k_profile_v1)
end

"""Call `scheme` on every column with zeroed temperature/humidity tendencies; return the tendencies."""
function column_tendencies(vars, model, scheme)
    fill!(vars.tendencies.grid.temperature.data, 0)
    fill!(vars.tendencies.grid.humidity.data, 0)
    for ij in 1:size(vars.grid.temperature_prev.data, 1)
        S.parameterization!(ij, vars, scheme, model)
    end
    return copy(vars.tendencies.grid.temperature.data), copy(vars.tendencies.grid.humidity.data)
end

"""Largest column-relative heat and water imbalance of a set of tendencies."""
function imbalance(dT, dq, model)
    dsigma = Float64.(model.geometry.σ_levels_thick)
    worst_heat = worst_water = 0.0
    for ij in axes(dT, 1)
        h, hs = sum(Float64.(dT[ij, :]) .* dsigma), sum(abs.(Float64.(dT[ij, :])) .* dsigma)
        w, ws = sum(Float64.(dq[ij, :]) .* dsigma), sum(abs.(Float64.(dq[ij, :])) .* dsigma)
        hs > 0 && (worst_heat = max(worst_heat, abs(h) / hs))
        ws > 0 && (worst_water = max(worst_water, abs(w) / ws))
    end
    return worst_heat, worst_water
end

face_coefficients(vars, model, ij, scale...) = begin
    p, n = vars.parameterizations, model.geometry.nlayers
    thetav = R._PBLThetaV(vars.grid.temperature_prev, vars.grid.humidity_prev, model.geometry.σ_levels_full,
                          model.atmosphere, Float32(model.atmosphere.κ), ij)
    a = R._PBLFaceCoefficients(p.pbl_height_work, vars.grid.temperature_prev, vars.grid.humidity_prev,
        model.geometry.σ_levels_full, model.geometry.σ_levels_half, model.atmosphere, ij, p.pbl_height[ij],
        p.pbl_friction_velocity[ij], p.pbl_buoyancy_flux[ij], thetav[n], Float32(model.planet.gravity),
        Float32(model.atmosphere.R_dry), scale...)
    [a[k] for k in 1:(n - 1)]
end

@testset "PBL K-profile strength-screen parameters" begin
    production = R.load_config(joinpath(R.PROJECT_ROOT, "config/production.yml"); check_input_files = false)
    grid = S.SpectralGrid(; NF = Float32, trunc = 7, nlayers = 8)
    on = with(production; atmosphere_vertical_diffusion = :holtslag_boville_k_profile_v1)

    @testset "configuration, provenance and refusals" begin
        @test production.atmosphere_pbl_diffusivity_scale == 1.0
        @test production.atmosphere_pbl_mix_heat && production.atmosphere_pbl_mix_humidity
        default = R._atmosphere_vertical_diffusion(on, grid)
        @test default.diffusivity_scale === 1.0f0 && default.mix_heat && default.mix_humidity
        @test R._atmosphere_vertical_diffusion_provenance(default) == R.PBL_K_PROFILE_PROVENANCE
        @test isbitstype(typeof(default))
        base = R.PBL_K_PROFILE_PROVENANCE
        for (kwargs, suffix) in (
                ((; atmosphere_pbl_diffusivity_scale = 0.3), "_scale0p3"),
                ((; atmosphere_pbl_diffusivity_scale = 0.1), "_scale0p1"),
                ((; atmosphere_pbl_mix_humidity = false), "_heat_only"),
                ((; atmosphere_pbl_mix_heat = false), "_humidity_only"),
                ((; atmosphere_pbl_diffusivity_scale = 0.3, atmosphere_pbl_mix_heat = false), "_scale0p3_humidity_only"))
            config = with(on; kwargs...)
            @test R.validate(config; check_input_files = false) === config
            scheme = R._atmosphere_vertical_diffusion(config, grid)
            @test R._atmosphere_vertical_diffusion_provenance(scheme) == base * suffix
            @test isbitstype(typeof(scheme))
            @test scheme.critical_richardson == default.critical_richardson && scheme.maximum_sigma == default.maximum_sigma
        end
        for bad in (0.0, -0.1, 1.5, NaN, Inf)
            @test_throws ArgumentError R.validate(with(on; atmosphere_pbl_diffusivity_scale = bad); check_input_files = false)
            @test_throws ArgumentError R.KProfileBoundaryLayerMixing(grid; diffusivity_scale = bad)
        end
        @test_throws ArgumentError R.validate(with(on; atmosphere_pbl_mix_heat = false, atmosphere_pbl_mix_humidity = false);
                                              check_input_files = false)
        @test_throws ArgumentError R.KProfileBoundaryLayerMixing(grid; mix_heat = false, mix_humidity = false)
        # Non-default values need the K-profile: the zero-operator control refuses each of them.
        for kwargs in ((; atmosphere_pbl_diffusivity_scale = 0.3), (; atmosphere_pbl_mix_heat = false),
                       (; atmosphere_pbl_mix_humidity = false))
            @test_throws ArgumentError R.validate(with(production; kwargs...); check_input_files = false)
        end
        @test R.validate(production; check_input_files = false) === production
        # YAML keys are parsed.
        mktempdir() do dir
            lines = readlines(joinpath(R.PROJECT_ROOT, "config/production.yml"))
            index = findall(l -> startswith(l, "atmosphere_vertical_diffusion:"), lines)
            @test length(index) == 1
            lines[only(index)] = "atmosphere_vertical_diffusion: holtslag_boville_k_profile_v1"
            push!(lines, "atmosphere_pbl_diffusivity_scale: 0.3", "atmosphere_pbl_mix_heat: false")
            path = joinpath(dir, "screen.yml"); write(path, join(lines, "\n") * "\n")
            parsed = R.load_config(path; check_input_files = false)
            @test parsed.atmosphere_pbl_diffusivity_scale == 0.3
            @test !parsed.atmosphere_pbl_mix_heat && parsed.atmosphere_pbl_mix_humidity
        end
    end

    @testset "compact coupled model, 27 layers" begin
        mktempdir() do dir
            config = compact_l27(joinpath(dir, "out"))
            R.validate(config)
            simulation = R.build_dynamic_esm(config; steps = 4)
            model = simulation.model.atmosphere.model
            vars = simulation.model.atmosphere.variables
            R.Oceananigans.run!(simulation)
            @test simulation.model.clock.iteration == 4
            sg = model.spectral_grid
            default = model.vertical_diffusion
            p = vars.parameterizations
            dT1, dq1 = column_tendencies(vars, model, default)
            K1, h1 = copy(Array(p.pbl_lowest_face_diffusivity.data)), copy(Array(p.pbl_height.data))
            npoints = size(dT1, 1)
            active = [ij for ij in 1:npoints if any(!=(0), dq1[ij, :])]
            println("PBL_STRENGTH_ACTIVE_COLUMNS ", length(active), " of ", npoints)
            @test length(active) > npoints ÷ 10
            heat1, water1 = imbalance(dT1, dq1, model)
            println("PBL_STRENGTH_DEFAULT_IMBALANCE heat=$heat1 water=$water1")
            @test heat1 < 1e-4 && water1 < 1e-4
            a1 = Dict(ij => face_coefficients(vars, model, ij) for ij in active[1:min(end, 20)])
            # An explicitly constructed default scheme gives bitwise the same tendencies.
            explicit = R.KProfileBoundaryLayerMixing(sg; diffusivity_scale = 1, mix_heat = true, mix_humidity = true)
            @test S.initialize!(explicit, model) === nothing
            @test column_tendencies(vars, model, explicit) == (dT1, dq1)
            for s in (0.1f0, 0.3f0)
                scheme = R.KProfileBoundaryLayerMixing(sg; diffusivity_scale = s)
                @test S.initialize!(scheme, model) === nothing
                dT, dq = column_tendencies(vars, model, scheme)
                @test all(isfinite, dT) && all(isfinite, dq)
                heat, water = imbalance(dT, dq, model)
                println("PBL_STRENGTH_SCALE $s imbalance heat=$heat water=$water abs_dq_ratio=",
                        sum(abs, Float64.(dq)) / sum(abs, Float64.(dq1)))
                @test heat < 1e-4 && water < 1e-4
                # Depth, u* and B do not change. The lowest-face diffusivity and every face coefficient scale exactly.
                @test Array(p.pbl_height.data) == h1
                @test Array(p.pbl_lowest_face_diffusivity.data) == K1 .* s
                for (ij, a) in a1
                    @test face_coefficients(vars, model, ij, s) == a .* s
                end
                # Weaker mixing: smaller total tendencies, and different from the default somewhere.
                @test sum(abs, Float64.(dq)) < sum(abs, Float64.(dq1)) && sum(abs, Float64.(dT)) < sum(abs, Float64.(dT1))
                @test dq != dq1
            end
            heat_only = R.KProfileBoundaryLayerMixing(sg; mix_humidity = false)
            dT, dq = column_tendencies(vars, model, heat_only)
            @test dT == dT1 && all(==(0), dq)
            humidity_only = R.KProfileBoundaryLayerMixing(sg; mix_heat = false)
            dT, dq = column_tendencies(vars, model, humidity_only)
            @test dq == dq1 && all(==(0), dT)
            # GPU interface on the CPU (only the device model's fields, plain arrays): bitwise equal for every variant.
            device_model = (; geometry = model.geometry, planet = model.planet, atmosphere = model.atmosphere,
                            time_stepping = model.time_stepping, land_sea_mask = model.land_sea_mask)
            plain(x) = copy(Array(x.data))
            g = vars.grid
            device_vars = (
                grid = (; temperature_prev = plain(g.temperature_prev), humidity_prev = plain(g.humidity_prev),
                        u_prev = plain(g.u_prev), v_prev = plain(g.v_prev), pressure_prev = plain(g.pressure_prev),
                        geopotential = plain(g.geopotential)),
                parameterizations = (; (name => plain(getproperty(p, name)) for name in
                    (:pbl_height, :pbl_friction_velocity, :pbl_buoyancy_flux, :pbl_drag_coefficient,
                     :pbl_lowest_face_diffusivity, :pbl_lowest_face_height, :pbl_height_work, :pbl_solution_work, :boundary_layer_drag))...),
                prognostic = (; land = (; sensible_heat_flux = plain(vars.prognostic.land.sensible_heat_flux),
                                        surface_humidity_flux = plain(vars.prognostic.land.surface_humidity_flux)),
                              ocean = (; sensible_heat_flux = plain(vars.prognostic.ocean.sensible_heat_flux),
                                         surface_humidity_flux = plain(vars.prognostic.ocean.surface_humidity_flux))),
                scratch = (; grid = (; a = plain(vars.scratch.grid.a), b = plain(vars.scratch.grid.b))),
                tendencies = (; grid = (; temperature = zero(dT1), humidity = zero(dq1))))
            for scheme in (default, R.KProfileBoundaryLayerMixing(sg; diffusivity_scale = 0.3f0), heat_only, humidity_only)
                expected = column_tendencies(vars, model, scheme)
                fill!(device_vars.tendencies.grid.temperature, 0); fill!(device_vars.tendencies.grid.humidity, 0)
                for ij in 1:npoints
                    S.parameterization!(ij, device_vars, scheme, device_model)
                end
                @test (device_vars.tendencies.grid.temperature, device_vars.tendencies.grid.humidity) == expected
                @test device_vars.parameterizations.pbl_lowest_face_diffusivity == Array(p.pbl_lowest_face_diffusivity.data)
            end
        end
    end
end
end
println("PBL_K_PROFILE_STRENGTH_UNIT_PASS")
