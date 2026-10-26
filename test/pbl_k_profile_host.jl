module PBLKProfileHostTests
using ReadyESM, Test
const R = ReadyESM
const S = R.SpeedyWeather

fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
with(config; kwargs...) = R.ExperimentConfig(; merge(fields(config), (; kwargs...))...)
include("pbl_test_config.jl")

"""A tile field that must never be read (Codex: a zero-area tile is not read, not only weighted by zero)."""
struct AbsentTile <: AbstractVector{Float32} end
Base.size(::AbsentTile) = (4,)
Base.getindex(::AbsentTile, ::Int) = error("absent tile read")

"""Height cap of column `ij`: the height of sigma = scheme.maximum_sigma."""
column_cap(vars, model, ij) = R._pbl_sigma_height(vars.parameterizations.pbl_height_work, vars.grid.temperature_prev,
    vars.grid.humidity_prev, model.geometry.σ_levels_full, model.atmosphere, ij, model.vertical_diffusion.maximum_sigma,
    Float32(model.atmosphere.R_dry), Float32(model.planet.gravity))

"""Compact analytic CPU fixture: T7/L27 with a 32x16x4 ocean."""
function compact_coupled(output)
    base = analytic_pbl_test_config()
    return with(base; device = :cpu, truncation = 7, nlayers = 27, ocean_grid = :idealized,
        ocean_nlongitude = 32, ocean_nlatitude = 16, ocean_nlayers = 4, ocean_initial_conditions = :analytic,
        ecco_initial_conditions_directory = "", duration_days = 4 * 3600 / 86_400, analysis_start_day = 0.0,
        output_dir = output,
        atmosphere_vertical_diffusion = :holtslag_boville_k_profile_v1,
        # The analytic L8 launch has no active mixing face after four steps.
        # L27 resolves the boundary layer without weakening the K > 0 check.
        # Keep this legacy-advection fixture separate from the full-temperature
        # strength and momentum tests, which also use L27.
        atmosphere_vertical_advection = :speedyweather_0_21_1_centered)
end

"""Call the scheme on every column with zeroed temperature/humidity tendencies; return the tendencies."""
function column_tendencies(vars, model)
    fill!(vars.tendencies.grid.temperature.data, 0)
    fill!(vars.tendencies.grid.humidity.data, 0)
    for ij in 1:size(vars.grid.temperature_prev.data, 1)
        S.parameterization!(ij, vars, model.vertical_diffusion, model)
    end
    return copy(vars.tendencies.grid.temperature.data), copy(vars.tendencies.grid.humidity.data)
end

function check_columns(vars, model; label)
    dT, dq = column_tendencies(vars, model)
    @test all(isfinite, dT) && all(isfinite, dq)
    dsigma = Float64.(model.geometry.σ_levels_thick)
    cp = Float64(model.atmosphere.heat_capacity)
    npoints, n = size(dT)
    heat_error = water_error = 0.0
    inactive_ok = true
    scheme = model.vertical_diffusion
    for ij in 1:npoints
        heat = sum(cp .* Float64.(dT[ij, :]) .* dsigma)
        heat_scale = sum(cp .* abs.(Float64.(dT[ij, :])) .* dsigma)
        water = sum(Float64.(dq[ij, :]) .* dsigma)
        water_scale = sum(abs.(Float64.(dq[ij, :])) .* dsigma)
        heat_scale > 0 && (heat_error = max(heat_error, abs(heat) / heat_scale))
        water_scale > 0 && (water_error = max(water_error, abs(water) / water_scale))
        # Rebuild this column's face coefficients; cells with both faces inactive must have no tendency.
        p = vars.parameterizations
        a = R._PBLFaceCoefficients(p.pbl_height_work, vars.grid.temperature_prev, vars.grid.humidity_prev,
            model.geometry.σ_levels_full, model.geometry.σ_levels_half, model.atmosphere, ij, p.pbl_height[ij],
            p.pbl_friction_velocity[ij], p.pbl_buoyancy_flux[ij],
            R._PBLThetaV(vars.grid.temperature_prev, vars.grid.humidity_prev, model.geometry.σ_levels_full,
                         model.atmosphere, Float32(model.atmosphere.κ), ij)[n],
            Float32(model.planet.gravity), Float32(model.atmosphere.R_dry))
        for k in 1:n
            above = k > 1 ? a[k - 1] : 0.0f0
            below = k < n ? a[k] : 0.0f0
            if above == 0 && below == 0
                inactive_ok &= dT[ij, k] == 0 && dq[ij, k] == 0
            end
        end
    end
    println("PBL_HOST_COLUMN_CHECK $label max_relative_heat=$heat_error max_relative_water=$water_error")
    @test heat_error < 1e-4
    @test water_error < 1e-4
    @test inactive_ok
    h = Array(vars.parameterizations.pbl_height.data)
    caps = [column_cap(vars, model, ij) for ij in 1:npoints]
    @test all(isfinite, h) && all(h .<= caps .* (1 + 1f-6)) && all(3000 .< caps .< 8000)
    return dT, dq
end

@testset "PBL K-profile host wiring" begin
    production = legacy_pbl_test_config()
    grid = S.SpectralGrid(; NF = Float32, trunc = 7, nlayers = 8)

    @testset "off and validation" begin
        @test production.atmosphere_vertical_diffusion == :speedyweather_0_21_1_zero_operator_control
        @test R._atmosphere_vertical_diffusion(production, grid) isa S.BulkRichardsonDiffusion
        on = with(production; atmosphere_vertical_diffusion = :holtslag_boville_k_profile_v1)
        @test R.validate(on; check_input_files = false) === on
        scheme = R._atmosphere_vertical_diffusion(on, grid)
        @test scheme isa R.KProfileBoundaryLayerMixing && scheme.critical_richardson == 0.3f0
        @test scheme.maximum_sigma == 0.5f0
        # Surface parameters carried by the scheme (the GPU kernel's model has no surface_condition).
        condition, momentum = S.SurfaceCondition(grid), S.SurfaceMomentumFlux(grid)
        @test scheme.surface_wind_slowdown == condition.wind_slowdown && scheme.gust_speed == condition.gust_speed
        @test scheme.momentum_wind_slowdown == momentum.wind_slowdown && scheme.drag_land == momentum.drag_land &&
              scheme.drag_ocean == momentum.drag_ocean && scheme.use_boundary_layer_drag == momentum.use_boundary_layer_drag
        @test isbitstype(typeof(scheme))
        @test occursin("calm_wind_floor_0p1ms", R.PBL_K_PROFILE_PROVENANCE)
        @test R._atmosphere_vertical_diffusion_provenance(scheme) == R.PBL_K_PROFILE_PROVENANCE
        @test_throws ArgumentError R.validate(with(on; land_model = :speedy_bucket); check_input_files = false)
        @test_throws ArgumentError R.validate(with(on; ocean_model = :slab); check_input_files = false)
        @test_throws ArgumentError R.validate(with(on; atmosphere_moist_mass_coupling = true); check_input_files = false)
    end

    @testset "absent tiles and the sigma cap" begin
        @test R._pbl_area_weighted(AbsentTile(), fill(2.0f0, 4), 3, 0.0f0) == 2.0f0
        @test R._pbl_area_weighted(fill(5.0f0, 4), AbsentTile(), 3, 1.0f0) == 5.0f0
        @test R._pbl_area_weighted(fill(5.0f0, 4), fill(2.0f0, 4), 3, 0.25f0) ≈ 0.25f0 * 5 + 0.75f0 * 2
        @test_throws ErrorException R._pbl_area_weighted(AbsentTile(), fill(2.0f0, 4), 3, 0.5f0)
        # Isothermal dry column: z(sigma) = R T / g log(1 / sigma) exactly, for any target between or below levels.
        atmosphere = S.EarthAtmosphere(grid)
        sigma = Float64[0.05, 0.2, 0.45, 0.7, 0.9, 0.99]
        Rd, g, T0 = 287.0, 9.81, 260.0
        T = fill(T0, 1, 6); q = zeros(1, 6)
        heights = reshape(Rd * T0 / g .* log.(1 ./ sigma), 1, 6)
        for target in (0.5, 0.995, 0.3, 0.02)
            @test R._pbl_sigma_height(heights, T, q, sigma, atmosphere, 1, target, Rd, g) ≈ Rd * T0 / g * log(1 / target)
        end
    end

    @testset "compact coupled model" begin
        mktempdir() do dir
            config = compact_coupled(joinpath(dir, "out"))
            R.validate(config)
            simulation = R.build_dynamic_esm(config; steps = 4)
            model = simulation.model.atmosphere.model
            vars = simulation.model.atmosphere.variables
            @test model.vertical_diffusion isa R.KProfileBoundaryLayerMixing
            # Start-up: the stored drag can be zero; the calm-safe closure must still be finite and conservative.
            check_columns(vars, model; label = "startup")
            R.Oceananigans.run!(simulation)
            @test simulation.model.clock.iteration == 4
            @test all(isfinite, Array(S.get_step(vars.prognostic.temperature, 1).data))
            p = vars.parameterizations
            ustar = Array(p.pbl_friction_velocity.data)
            drag = Array(p.pbl_drag_coefficient.data)
            println("PBL_HOST_AFTER_STEPS ustar_mean=$(sum(ustar) / length(ustar)) drag_mean=$(sum(drag) / length(drag)) " *
                    "h_mean=$(sum(Array(p.pbl_height.data)) / length(ustar)) " *
                    "B_range=$(extrema(Array(p.pbl_buoyancy_flux.data)))")
            @test all(isfinite, ustar) && all(>=(0), ustar) && sum(ustar) > 0
            @test all(isfinite, Array(p.pbl_buoyancy_flux.data))
            @test any(>(0), Array(p.pbl_lowest_face_diffusivity.data))
            dT, dq = check_columns(vars, model; label = "after_4_steps")
            # GPU interface on the CPU: a model holding only the fields of SpeedyWeather's device model, and plain
            # arrays instead of Fields. The tendencies must be bitwise those of the full-model, Field-based call.
            device_model = (; geometry = model.geometry, planet = model.planet, atmosphere = model.atmosphere,
                            time_stepping = model.time_stepping, land_sea_mask = model.land_sea_mask)
            plain(x) = copy(Array(x.data))
            g, p = vars.grid, vars.parameterizations
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
                tendencies = (; grid = (; temperature = zero(dT), humidity = zero(dq))))
            for ij in axes(dT, 1)
                S.parameterization!(ij, device_vars, model.vertical_diffusion, device_model)
            end
            @test device_vars.tendencies.grid.temperature == dT
            @test device_vars.tendencies.grid.humidity == dq
            @test device_vars.parameterizations.pbl_height == Array(vars.parameterizations.pbl_height.data)
            # initialize! refuses a scheme whose surface parameters differ from the model's, for every copied value.
            sg = model.spectral_grid
            for (condition, momentum) in (
                    (S.SurfaceCondition(sg; gust_speed = 3), S.SurfaceMomentumFlux(sg)),
                    (S.SurfaceCondition(sg; wind_slowdown = 0.9), S.SurfaceMomentumFlux(sg)),
                    (S.SurfaceCondition(sg), S.SurfaceMomentumFlux(sg; wind_slowdown = 0.9)),
                    (S.SurfaceCondition(sg), S.SurfaceMomentumFlux(sg; use_boundary_layer_drag = false)),
                    (S.SurfaceCondition(sg), S.SurfaceMomentumFlux(sg; drag_ocean = 1.5e-3)),
                    (S.SurfaceCondition(sg), S.SurfaceMomentumFlux(sg; drag_land = 2.0e-3)))
                other = R.KProfileBoundaryLayerMixing(sg; surface_condition = condition, surface_momentum_flux = momentum)
                @test_throws ArgumentError S.initialize!(other, model)
            end
            @test S.initialize!(model.vertical_diffusion, model) === nothing
            # Fixed-drag mode (use_boundary_layer_drag = false) through both interfaces: bitwise equal, and the drag
            # is the area-weighted fixed land and ocean coefficient.
            fixed = R.KProfileBoundaryLayerMixing(sg; surface_momentum_flux = S.SurfaceMomentumFlux(sg; use_boundary_layer_drag = false))
            fill!(vars.tendencies.grid.temperature.data, 0); fill!(vars.tendencies.grid.humidity.data, 0)
            for ij in axes(dT, 1)
                S.parameterization!(ij, vars, fixed, model)
            end
            fixed_dT = copy(vars.tendencies.grid.temperature.data); fixed_dq = copy(vars.tendencies.grid.humidity.data)
            fill!(device_vars.tendencies.grid.temperature, 0); fill!(device_vars.tendencies.grid.humidity, 0)
            for ij in axes(dT, 1)
                S.parameterization!(ij, device_vars, fixed, device_model)
            end
            @test device_vars.tendencies.grid.temperature == fixed_dT && device_vars.tendencies.grid.humidity == fixed_dq
            f = Float32.(Array(model.land_sea_mask.mask.data))
            @test device_vars.parameterizations.pbl_drag_coefficient ≈ f .* fixed.drag_land .+ (1 .- f) .* fixed.drag_ocean
            @test Array(vars.parameterizations.pbl_drag_coefficient.data) == device_vars.parameterizations.pbl_drag_coefficient
            # B_s is the surface flux of Tv sigma^-kappa: finite difference of T_s (1 + mu q) along (H/(rho cp), E/rho).
            atmosphere = model.atmosphere
            mu, kappa, Rd = Float64(atmosphere.μ_virt_temp), Float64(atmosphere.κ), Float64(atmosphere.R_dry)
            cp, n = Float64(atmosphere.heat_capacity), model.geometry.nlayers
            sigma_n = Float64(model.geometry.σ_levels_full[n])
            worst = 0.0
            for ij in axes(dT, 1)
                T1 = Float64(vars.grid.temperature_prev[ij, n]); q1 = Float64(vars.grid.humidity_prev[ij, n])
                Ts = T1 * sigma_n^(-kappa)
                rho = Float64(vars.grid.pressure_prev[ij]) / (Rd * T1 * (1 + mu * q1) * sigma_n^(-kappa))
                f = Float64(model.land_sea_mask.mask[ij])
                H = (f > 0 ? f * vars.prognostic.land.sensible_heat_flux[ij] : 0.0) +
                    (f < 1 ? (1 - f) * vars.prognostic.ocean.sensible_heat_flux[ij] : 0.0)
                E = (f > 0 ? f * vars.prognostic.land.surface_humidity_flux[ij] : 0.0) +
                    (f < 1 ? (1 - f) * vars.prognostic.ocean.surface_humidity_flux[ij] : 0.0)
                eps = 1e-3
                thetav(T, q) = T * (1 + mu * q)
                fd = (thetav(Ts + eps * H / (rho * cp), q1 + eps * E / rho) - thetav(Ts, q1)) / eps
                B = Float64(vars.parameterizations.pbl_buoyancy_flux[ij])
                worst = max(worst, abs(B - fd) / (abs(fd) + 1e-4))
            end
            println("PBL_HOST_BUOYANCY_CONVENTION max_relative=$worst")
            @test worst < 1e-3
            # Mixing must actually act (at L8 usually only the lowest face is inside the boundary layer, which
            # makes the two-layer exchange conserve to the bit).
            println("PBL_HOST_ACTIVE_COLUMNS ", count(ij -> any(!=(0), dq[ij, :]), axes(dq, 1)), " of ", size(dq, 1))
            @test count(!=(0), dq) > 0 && count(!=(0), dT) > 0
        end
    end
end
end
println("PBL_K_PROFILE_HOST_UNIT_PASS")
