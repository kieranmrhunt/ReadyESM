module FullTemperatureVerticalAdvectionTests
# Port of upstream SpeedyWeather PR #1286 (full_temperature_vertical_advection.jl): configuration, dispatch, the
# upstream reference-profile testset adapted to SpeedyWeather 0.21.1, an independent check that the fix adds exactly
# the vertical advection of the reference profile, and the coupled ReadyESM path.
using ReadyESM, Test
const R = ReadyESM
const S = R.SpeedyWeather

fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
with(config; kwargs...) = R.ExperimentConfig(; merge(fields(config), (; kwargs...))...)

"""Grid temperature tendency of the dynamical core alone (as in the upstream test)."""
function grid_temperature_tendency(vars, model; method = nothing)
    S.reset_tendencies!(vars)
    if isnothing(method)
        S.dynamics_tendencies!(vars, 2, model)
    else
        invoke(S.dynamics_tendencies!, method, vars, 2, model)
    end
    return copy(vars.tendencies.grid.temperature.data)
end

"""The explicit kernel's vertical advection of a horizontally uniform profile Tk (centred, order 2, clamped stencil)."""
function reference_profile_advection(w, Tk, dsigma)
    npoints, n = size(w)
    out = zeros(Float64, npoints, n)
    for k in 1:n, ij in 1:npoints
        km, kp = max(1, k - 1), k
        wm, wp = Float64(w[ij, km]), Float64(w[ij, kp])
        s = (clamp(k - 1, 1, n), k, clamp(k + 1, 1, n))
        face_p = (Float64(Tk[s[2]]) + Float64(Tk[s[3]])) / 2
        face_m = (Float64(Tk[s[1]]) + Float64(Tk[s[2]])) / 2
        out[ij, k] = -(wp * face_p - wm * face_m - Float64(Tk[k]) * (wp - wm)) / Float64(dsigma[k])
    end
    return out
end

@testset "full-temperature vertical advection (upstream PR #1286 port)" begin
    production = R.load_config(joinpath(R.PROJECT_ROOT, "config/production.yml"); check_input_files = false)
    grid = S.SpectralGrid(; NF = Float32, trunc = 7, nlayers = 8)

    @testset "configuration and dispatch" begin
        @test production.atmosphere_vertical_advection == :centered_full_temperature
        @test R.ExperimentConfig().atmosphere_vertical_advection == :speedyweather_0_21_1_centered
        # Explicit legacy settings and configurations without the key retain the old dynamics.
        legacy = with(production; atmosphere_vertical_advection = :speedyweather_0_21_1_centered)
        @test R._atmosphere_vertical_advection_scheme(legacy, grid) === S.CenteredVerticalAdvection(grid)
        on = with(production; atmosphere_vertical_advection = :centered_full_temperature)
        @test R.validate(on; check_input_files = false) === on
        scheme = R._atmosphere_vertical_advection_scheme(on, grid)
        @test scheme isa R.FullTemperatureCenteredVerticalAdvection{Float32, 1}
        @test scheme isa S.DispersiveVerticalAdvection && isbitstype(typeof(scheme))
        @test R._atmosphere_vertical_advection_provenance(scheme) == R.FULL_TEMPERATURE_VERTICAL_ADVECTION_PROVENANCE
        @test R._atmosphere_vertical_advection_provenance(S.CenteredVerticalAdvection(grid)) ==
              R.DEFAULT_VERTICAL_ADVECTION_PROVENANCE
        @test_throws ArgumentError R.validate(with(production; atmosphere_vertical_advection = :upwind);
                                              check_input_files = false)
        mktempdir() do dir
            lines = filter(line -> !startswith(line, "atmosphere_vertical_advection:"),
                           readlines(joinpath(R.PROJECT_ROOT, "config/production.yml")))
            path = joinpath(dir, "fix.yml"); write(path, join(lines, "\n") * "\n")
            @test R.load_config(path; check_input_files = false).atmosphere_vertical_advection == :speedyweather_0_21_1_centered
            for option in (:centered_full_temperature, :speedyweather_0_21_1_centered)
                write(path, join([lines; "atmosphere_vertical_advection: $option"], "\n") * "\n")
                @test R.load_config(path; check_input_files = false).atmosphere_vertical_advection == option
            end
        end
        # Face reconstruction is exactly the centred one.
        xi = reshape(Float32[1, 4, 9, 16, 25], 1, 5)
        for B in (1, 2), stencil in ((1, 2, 3, 4, 5)[1:2B],)
            @test S.reconstructed_at_face(xi, 1, stencil, 1.0f0, R.FullTemperatureCenteredVerticalAdvection{Float32, B}()) ===
                  S.reconstructed_at_face(xi, 1, stencil, 1.0f0, S.CenteredVerticalAdvection{Float32, B}())
        end
        # Dispatch: default models keep SpeedyWeather's method, the fix gets the copy, both options get the combined copy.
        m(va, sp) = S.PrimitiveWetModel(grid; vertical_advection = va, stochastic_physics = sp, dynamics_only = true)
        V = S.Variables
        default_model = S.PrimitiveWetModel(grid; dynamics_only = true)
        fixed_model = m(R.FullTemperatureCenteredVerticalAdvection(grid), default_model.stochastic_physics)
        moist_model = m(S.CenteredVerticalAdvection(grid), R.MoistMassCoupling())
        both_model = m(R.FullTemperatureCenteredVerticalAdvection(grid), R.MoistMassCoupling())
        method(model) = which(S.dynamics_tendencies!, Tuple{V, Int, typeof(model)})
        file(model) = basename(String(method(model).file))
        @test method(default_model).module === S && file(default_model) == "tendencies.jl"
        @test fixed_model isa R.FullTemperatureAdvectionModel && !(default_model isa R.FullTemperatureAdvectionModel)
        @test file(fixed_model) == "full_temperature_vertical_advection.jl"
        @test file(moist_model) == "moist_mass_coupling.jl"
        @test both_model isa R.MoistMassFullTemperatureModel && both_model isa R.MoistMassCoupledModel
        @test file(both_model) == "full_temperature_vertical_advection.jl" && method(both_model) !== method(fixed_model)
        @test method(both_model).line > method(fixed_model).line   # the combined copy, defined second
    end

    @testset "restart record" begin
        sim(scheme) = (; model = (; vertical_advection = scheme))
        fixed, default = sim(R.FullTemperatureCenteredVerticalAdvection(grid)), sim(S.CenteredVerticalAdvection(grid))
        @test R._full_temperature_restart_state(default) == (;)            # default checkpoints are unchanged
        @test R._full_temperature_restart_state(fixed) == (; full_temperature_vertical_advection = true)
        @test R._validate_full_temperature_restart(default, (; variables = 1)) === nothing      # legacy checkpoint
        @test R._validate_full_temperature_restart(fixed, R._full_temperature_restart_state(fixed)) === nothing
        @test_throws ArgumentError R._validate_full_temperature_restart(fixed, (; variables = 1))
        @test_throws ArgumentError R._validate_full_temperature_restart(default, (; full_temperature_vertical_advection = true))
        @test_throws ArgumentError R._validate_full_temperature_restart(fixed, (; full_temperature_vertical_advection = 1))
    end

    @testset "upstream reference-profile test, adapted to 0.21.1" begin
        sg = S.SpectralGrid(; NF = Float32, trunc = 15, nlayers = 8)
        for (label, scheme, fixed) in (("centred (0.21.1 default)", S.CenteredVerticalAdvection(sg), false),
                                       ("upwind", S.UpwindVerticalAdvection(sg), true),
                                       ("centred full temperature", R.FullTemperatureCenteredVerticalAdvection(sg), true))
            model = S.PrimitiveWetModel(sg; vertical_advection = scheme, dynamics_only = true)
            model.feedback.verbose = false
            simulation = S.initialize!(model)
            S.run!(simulation; period = S.Hour(6))      # develop a nonzero vertical velocity
            vars = simulation.variables
            tendency = grid_temperature_tendency(vars, model)
            change = 10 .* model.geometry.σ_levels_full    # change the stratification of the reference profile
            model.implicit.temp_profile .+= change
            changed = grid_temperature_tendency(vars, model)
            divergence = vars.grid.divergence.data
            passes = isapprox(changed, tendency .- change' .* divergence)
            println("PBL_PR1286_UPSTREAM_TEST $label passes=$passes max_dev=",
                    maximum(abs.(changed .- (tendency .- change' .* divergence))))
            @test passes == fixed     # the 0.21.1 centred scheme fails, as upstream reports
        end
    end

    @testset "the fix adds exactly the vertical advection of the reference profile" begin
        sg = S.SpectralGrid(; NF = Float32, trunc = 15, nlayers = 8)
        model = S.PrimitiveWetModel(sg; vertical_advection = R.FullTemperatureCenteredVerticalAdvection(sg),
                                    dynamics_only = true)
        model.feedback.verbose = false
        simulation = S.initialize!(model)
        S.run!(simulation; period = S.Hour(6))
        vars = simulation.variables
        fixed = grid_temperature_tendency(vars, model)
        u_fixed = copy(vars.tendencies.grid.u.data); q_fixed = copy(vars.tendencies.grid.humidity.data)
        original = grid_temperature_tendency(vars, model; method = Tuple{S.Variables, Integer, S.PrimitiveEquation})
        @test vars.tendencies.grid.u.data == u_fixed && vars.tendencies.grid.humidity.data == q_fixed
        expected = reference_profile_advection(vars.dynamics.w.data, model.implicit.temp_profile,
                                               model.geometry.σ_levels_thick)
        difference = Float64.(fixed) .- Float64.(original)
        scale = maximum(abs, expected)
        worst = maximum(abs, difference .- expected) / scale
        println("PBL_PR1286_REFERENCE_ADVECTION max_relative_deviation=$worst scale=$scale")
        @test scale > 0 && worst < 1e-3
    end

    @testset "coupled ReadyESM model with the fix" begin
        mktempdir() do dir
            base = production
            era5 = joinpath(R.PROJECT_ROOT, "data/era5_initial_state")
            config = with(base; device = :cpu, truncation = 7, nlayers = 8, ocean_grid = :idealized,
                ocean_nlongitude = 32, ocean_nlatitude = 16, ocean_nlayers = 4, ocean_initial_conditions = :analytic,
                ecco_initial_conditions_directory = "", duration_days = 4 * 3600 / 86_400, analysis_start_day = 0.0,
                output_dir = joinpath(dir, "out"),
                era5_pressure_levels_path = joinpath(era5, "era5_pressure_levels_stratosphere_19930101T0000Z.nc"),
                era5_single_levels_path = joinpath(era5, "era5_single_levels_19930101T0000Z.nc"),
                era5_land_state_path = joinpath(era5, "era5_land_state_19930101T0000Z.nc"),
                atmosphere_vertical_advection = :centered_full_temperature)
            R.validate(config)
            simulation = R.build_dynamic_esm(config; steps = 4)
            model = simulation.model.atmosphere.model
            @test model.vertical_advection isa R.FullTemperatureCenteredVerticalAdvection
            @test model isa R.FullTemperatureAdvectionModel
            R.Oceananigans.run!(simulation)
            @test simulation.model.clock.iteration == 4
            vars = simulation.model.atmosphere.variables
            @test all(isfinite, Array(S.get_step(vars.prognostic.temperature, 1).data))
            diagnostics = R.collect_dynamic_diagnostics(simulation, config)
            @test diagnostics.metadata.atmosphere_vertical_advection == R.FULL_TEMPERATURE_VERTICAL_ADVECTION_PROVENANCE
        end
    end
end
end
println("PBL_PR1286_PORT_UNIT_PASS")
