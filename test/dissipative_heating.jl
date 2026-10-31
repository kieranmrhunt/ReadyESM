module DissipativeHeatingTests
# Dissipative heating (atmosphere_dissipative_heating): configuration, the legacy components when off, the
# kinetic-energy accounting of each wrapped operator on a compact coupled model, diagnose_v1 bitwise equal to off, and
# the restart record.
using ReadyESM, Test
const R = ReadyESM
const S = R.SpeedyWeather
const O = R.Oceananigans

fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
with(config; kwargs...) = R.ExperimentConfig(; merge(fields(config), (; kwargs...))...)
include("pbl_test_config.jl")
bits(x::Float32) = reinterpret(UInt32, x)
bits(x::Complex) = (bits(real(x)), bits(imag(x)))
same_bits(a, b) = size(a) == size(b) && all(bits(x) == bits(y) for (x, y) in zip(a, b))

function compact_l27(output, mode; steps = 4)
    return with(analytic_pbl_test_config(); device = :cpu, truncation = 7, nlayers = 27, ocean_grid = :idealized,
        ocean_nlongitude = 32, ocean_nlatitude = 16, ocean_nlayers = 4, ocean_initial_conditions = :analytic,
        ecco_initial_conditions_directory = "", duration_days = steps * 3600 / 86_400, analysis_start_day = 0.0,
        output_dir = output, atmosphere_vertical_diffusion = :holtslag_boville_k_profile_v1,
        atmosphere_pbl_mix_momentum = true, atmosphere_dissipative_heating = mode)
end

function run_compact(dir, mode)
    config = compact_l27(joinpath(dir, String(mode)), mode)
    R.validate(config)
    simulation = R.build_dynamic_esm(config; steps = 4)
    O.run!(simulation)
    @test simulation.model.clock.iteration == 4
    return config, simulation
end

prognostic(vars) = map(n -> copy(Array(S.get_step(getproperty(vars.prognostic, n), 1).data)),
                       (:temperature, :humidity, :vorticity, :divergence, :pressure))

zero_tendencies!(vars) = foreach(x -> fill!(x.data, 0),
    (vars.tendencies.grid.temperature, vars.tendencies.grid.humidity, vars.tendencies.grid.u, vars.tendencies.grid.v))

relative(a, b; floor = 1e-6) = maximum(abs.(a .- b) ./ max.(abs.(b), floor))

@testset "dissipative heating" begin
    production = legacy_pbl_test_config()
    grid = S.SpectralGrid(; NF = Float32, trunc = 7, nlayers = 8)

    @testset "configuration" begin
        released = R.load_config(joinpath(R.PROJECT_ROOT, "config/production.yml"); check_input_files = false)
        @test released.atmosphere_dissipative_heating == :kinetic_energy_v1
        @test R.ExperimentConfig().atmosphere_dissipative_heating == :off
        mktempdir() do dir
            path = joinpath(dir, "config.yml")
            production_text = read(joinpath(R.PROJECT_ROOT, "config/production.yml"), String)
            write(path, replace(production_text,
                "atmosphere_dissipative_heating: kinetic_energy_v1" =>
                "atmosphere_dissipative_heating: diagnose_v1"))
            @test R.load_config(path; check_input_files = false).atmosphere_dissipative_heating == :diagnose_v1
        end
        for mode in R.DISSIPATIVE_HEATING_MODES
            config = with(production; atmosphere_dissipative_heating = mode)
            @test R.validate(config; check_input_files = false) === config
        end
        @test_throws ArgumentError R.validate(with(production; atmosphere_dissipative_heating = :on);
                                              check_input_files = false)
        @test_throws ArgumentError R.KProfileBoundaryLayerMixing(grid; dissipative_heating = true)
        # The speed-limit wrapper needs ReadyESM's coupled drag, which has the surface branch.
        @test_throws MethodError R.DissipativeSpeedLimitDrag(S.SpeedLimitDrag(grid), true)
    end

    @testset "K-profile heating allocation on sheared columns" begin
        n = 12
        dsigma = fill(1.0 / n, n)
        a = [0.02 + 0.01 * sin(k) for k in 1:(n - 1)]
        dt = 1800.0
        own_cools = false
        for wind in ([10 * exp(-((k - 4) / 2)^2) for k in 1:n], [2.0 * k for k in 1:n],
                     [k == 9 ? 18.0 : 3.0 for k in 1:n])
            x_new, c, d = zeros(n), zeros(n), zeros(n)
            R._pbl_mix_wind!(zeros(1, n), wind, 1, x_new, a, dsigma, dt, c, d)
            heat = zeros(1, n)
            R._pbl_heat_from_wind!(heat, wind, 1, x_new, a, dsigma, dt, 1.0)
            loss = sum(dsigma .* (wind .^ 2 .- x_new .^ 2)) / (2 * dt)
            @test loss > 0
            @test all(>=(0), heat)
            @test sum(heat[1, :] .* dsigma) ≈ loss rtol = 1e-10
            own_cools |= any(<(0), wind .^ 2 .- x_new .^ 2)   # a layer's own KE change would cool it
        end
        @test own_cools
    end

    @testset "off builds the legacy components" begin
        @test R._atmosphere_horizontal_diffusion(production, grid) isa S.HyperDiffusion
        @test R._atmosphere_speed_limit_drag(production, grid) isa R.CoupledSurfaceSpeedLimitDrag
        @test R._dissipative_surface_momentum_flux(production, grid) === (;)
        @test R._dissipative_heating_provenance(production) === (;)
        for (mode, heat) in ((:diagnose_v1, false), (:kinetic_energy_v1, true))
            on = with(production; atmosphere_dissipative_heating = mode)
            diffusion = R._atmosphere_horizontal_diffusion(on, grid)
            @test diffusion isa R.DissipativeHyperDiffusion && diffusion.heat == heat
            @test R._plain_horizontal_diffusion(diffusion) isa S.HyperDiffusion
            drag = R._atmosphere_speed_limit_drag(on, grid)
            @test drag isa R.DissipativeSpeedLimitDrag && drag.heat == heat
            flux = R._dissipative_surface_momentum_flux(on, grid).surface_momentum_flux
            @test flux isa R.DissipativeSurfaceMomentumFlux && flux.heat == heat
            @test R._plain_surface_momentum_flux(flux) isa S.SurfaceMomentumFlux
            @test R._dissipative_heating_provenance(on).atmosphere_dissipative_heating == String(mode)
        end
    end

    mktempdir() do dir
        _, off = run_compact(dir, :off)
        off_atmosphere = off.model.atmosphere
        off_model = off_atmosphere.model
        @test off_model.horizontal_diffusion isa S.HyperDiffusion
        @test off_model.drag isa R.CoupledSurfaceSpeedLimitDrag
        @test off_model.surface_momentum_flux isa S.SurfaceMomentumFlux
        @test !off_model.vertical_diffusion.dissipative_heating
        @test R._dissipative_heating_restart_state(off_atmosphere) === (;)
        @test !haskey(off_model.callbacks, :dissipative_heating_budget)
        @test !hasproperty(R._speedy_callback_state(off_model.callbacks), :dissipative_heating_budget)
        @test !hasproperty(O.prognostic_state(off_atmosphere), :atmosphere_dissipative_heating)
        @test_throws ArgumentError R._validate_dissipative_heating_restart(off_atmosphere,
            (; atmosphere_dissipative_heating = "kinetic_energy_v1"))
        off_state = prognostic(off_atmosphere.variables)

        for (mode, heats) in ((:diagnose_v1, false), (:kinetic_energy_v1, true))
            @testset "$mode on a compact coupled model" begin
                config, simulation = run_compact(dir, mode)
                atmosphere = simulation.model.atmosphere
                model, vars = atmosphere.model, atmosphere.variables
                p, tend, sg = vars.parameterizations, vars.tendencies.grid, model.spectral_grid
                n = model.geometry.nlayers
                npoints = size(tend.u.data, 1)
                dsigma = Float64.(Array(model.geometry.σ_levels_thick))
                gravity = Float64(model.planet.gravity)
                cp = Float64(model.atmosphere.heat_capacity)
                half_step = Float64(model.time_stepping.Δt_sec)
                ps = Float64.(Array(vars.grid.pressure_prev.data))
                mass = ps ./ gravity

                state = prognostic(vars)
                if mode == :diagnose_v1
                    # Recording the losses changes no prognostic bit.
                    @test all(same_bits(a, b) for (a, b) in zip(state, off_state))
                else
                    @test !same_bits(state[1], off_state[1])
                end

                @testset "restart record" begin
                    record = O.prognostic_state(atmosphere)
                    @test record.atmosphere_dissipative_heating == String(mode)
                    @test R._validate_dissipative_heating_restart(atmosphere, record) === nothing
                    other = mode == :diagnose_v1 ? "kinetic_energy_v1" : "diagnose_v1"
                    @test_throws ArgumentError R._validate_dissipative_heating_restart(atmosphere,
                        merge(record, (; atmosphere_dissipative_heating = other)))
                    @test_throws ArgumentError R._validate_dissipative_heating_restart(atmosphere, (;))
                    @test R._check_dissipative_heating_host(model, config) === nothing
                    @test_throws ErrorException R._check_dissipative_heating_host(off_model, config)
                end

                @testset "native-step budget" begin
                    budget = model.callbacks[:dissipative_heating_budget]
                    @test budget isa R.DissipativeHeatingBudgetCallback && budget.steps > 0
                    means = R.dissipative_heating_budget_means(budget)
                    @test means.seconds ≈ budget.steps * half_step
                    @test means.surface_drag > 0 && means.pbl_face > 0 && means.hyperdiffusion > 0
                    @test means.speed_limit_drag >= 0 && isfinite(means.filter)
                    @test length(means.hyperdiffusion_layers) == n && all(>(0), means.hyperdiffusion_layers)
                    println("DISSIPATIVE_HEATING_BUDGET mode=$mode ", means)
                    saved = R._speedy_callback_state(model.callbacks).dissipative_heating_budget
                    @test saved.steps == budget.steps && saved.hyperdiffusion_layers == budget.hyperdiffusion_layers
                    fresh = R.DissipativeHeatingBudgetCallback(sg)
                    R._restore_dissipative_heating_budget!(fresh, saved)
                    @test R.dissipative_heating_budget_means(fresh) == means
                    without = Base.structdiff(R._speedy_callback_state(model.callbacks),
                                              NamedTuple{(:dissipative_heating_budget,)})
                    @test_throws ArgumentError R._restore_speedy_callback_state!(model.callbacks, without)
                end

                @testset "surface wind stress" begin
                    flux = model.surface_momentum_flux
                    @test flux isa R.DissipativeSurfaceMomentumFlux && flux.heat == heats
                    zero_tendencies!(vars)
                    foreach(ij -> S.parameterization!(ij, vars, flux, model), 1:npoints)
                    du, dv, dT = copy(tend.u.data), copy(tend.v.data), copy(tend.temperature.data)
                    D = Float64.(Array(p.surface_drag_dissipation.data))
                    zero_tendencies!(vars)
                    foreach(ij -> S.surface_wind_stress!(ij, vars, flux.flux, model), 1:npoints)
                    # The wrapper leaves the momentum tendencies bitwise unchanged.
                    @test same_bits(du, tend.u.data) && same_bits(dv, tend.v.data)
                    Fu, Fv = Float64.(du[:, n]), Float64.(dv[:, n])
                    up, vp = Float64.(Array(vars.grid.u_prev.data))[:, n], Float64.(Array(vars.grid.v_prev.data))[:, n]
                    expected = mass .* dsigma[n] .* .-((up .+ half_step .* Fu) .* Fu .+ (vp .+ half_step .* Fv) .* Fv)
                    @test all(>=(0), D) && count(>(0), D) > npoints ÷ 2
                    @test relative(D, expected) < 1e-4
                    if heats
                        @test all(iszero, dT[:, 1:(n - 1)])
                        @test relative(cp .* Float64.(dT[:, n]) .* mass .* dsigma[n], D) < 1e-4
                    else
                        @test all(iszero, dT)
                    end
                    println("DISSIPATIVE_HEATING_SURFACE_DRAG mode=$mode point_mean_W_m2=", sum(D) / npoints)
                end

                @testset "K-profile momentum mixing" begin
                    pbl = model.vertical_diffusion
                    @test pbl isa R.KProfileBoundaryLayerMixing && pbl.mix_momentum && pbl.dissipative_heating == heats
                    # Momentum only, so that the heated scheme's temperature tendency is the returned heat alone.
                    settings = (; mix_heat = false, mix_humidity = false, mix_momentum = true)
                    plain = R.KProfileBoundaryLayerMixing(sg; settings...)
                    heated = R.KProfileBoundaryLayerMixing(sg; settings..., dissipative_heating = true)
                    @test S.initialize!(plain, model) === nothing && S.initialize!(heated, model) === nothing
                    results = map((plain, heated)) do scheme
                        zero_tendencies!(vars)
                        foreach(ij -> S.parameterization!(ij, vars, scheme, model), 1:npoints)
                        map(x -> copy(x.data), (tend.temperature, tend.humidity, tend.u, tend.v))
                    end
                    (dT0, dq0, du0, dv0), (dT1, dq1, du1, dv1) = results
                    @test same_bits(dq0, dq1) && same_bits(du0, du1) && same_bits(dv0, dv1)
                    @test all(iszero, dT0)
                    # Dissipation heats every layer: no layer cools when it gains momentum from a neighbour.
                    @test all(>=(0), dT1)
                    @test count(>(0), dT1) > npoints ÷ 10   # mixing reaches only the lowest layers
                    Dsum = Float64.(Array(p.pbl_momentum_dissipation_face.data)) .+
                           Float64.(Array(p.pbl_momentum_dissipation_implicit.data))
                    column = [mass[ij] * cp * sum(Float64.(dT1[ij, :]) .* dsigma) for ij in 1:npoints]
                    active = findall(>(1e-3), Dsum)
                    @test length(active) > npoints ÷ 10
                    # Float32 columns: the same 1e-2 bound as the K-profile KE identity test, and 1e-3 on the total.
                    @test relative(column[active], Dsum[active]) < 1e-2
                    @test abs(sum(column) - sum(Dsum)) < 1e-3 * sum(Dsum)
                    println("DISSIPATIVE_HEATING_K_PROFILE mode=$mode point_mean_W_m2=", sum(Dsum) / npoints,
                            " worst_relative=", relative(column[active], Dsum[active]))
                end

                @testset "hyperdiffusion" begin
                    diffusion = model.horizontal_diffusion
                    @test diffusion isa R.DissipativeHyperDiffusion && diffusion.heat == heats
                    spectral = vars.tendencies
                    foreach(x -> fill!(x.data, 0),
                            (spectral.vorticity, spectral.divergence, spectral.temperature, spectral.humidity))
                    vorticity = copy(S.get_step(vars.prognostic.vorticity, 1).data)
                    divergence = copy(S.get_step(vars.prognostic.divergence, 1).data)
                    S.horizontal_diffusion!(vars, diffusion, model)
                    loss = Float64.(Array(diffusion.layer_loss))
                    filter_change = Float64.(Array(diffusion.filter_change))
                    dt = 2 * model.time_stepping.Δt
                    Tv, Td = copy(spectral.vorticity.data), copy(spectral.divergence.data)
                    # Independent checks on the grid. The transform is linear, so the exact KE change of an increment
                    # δ from a is the sphere mean of (wind of a + δ/2) . (wind of δ).
                    weights = Float64.(Array(R._global_point_weights(sg)))
                    function grid_wind(vor, div)
                        S.get_step(vars.prognostic.vorticity, 2).data .= vor
                        S.get_step(vars.prognostic.divergence, 2).data .= div
                        S.transform!(vars, 2, model; initialize = true)
                        return Float64.(Array(vars.grid.u.data)), Float64.(Array(vars.grid.v.data))
                    end
                    function ke_change(vor, div, dvor, ddiv)
                        um, vm = grid_wind(vor .+ dvor ./ 2, div .+ ddiv ./ 2)
                        du, dv = grid_wind(dvor, ddiv)
                        return vec(sum(weights .* (um .* du .+ vm .* dv); dims = 1))
                    end
                    # The filter diagnostic against SpeedyWeather's own leapfrog on copies of the two time levels.
                    levels = map((:vorticity, :divergence)) do name
                        variable = getproperty(vars.prognostic, name)
                        old, now = copy(S.get_step(variable, 1)), copy(S.get_step(variable, 2))
                        tendency = getproperty(spectral, name)
                        unfiltered = old.data .+ dt .* tendency.data
                        before = copy(now.data)
                        S.leapfrog!(old, now, tendency, dt, 2, model.time_stepping)   # old <- filtered now, now <- new
                        (; now = before, unfiltered, filtered_now = copy(old.data), filtered_new = copy(now.data))
                    end
                    vor, div = levels
                    filter_expected = (ke_change(vor.now, div.now, vor.filtered_now .- vor.now,
                                                 div.filtered_now .- div.now) .+
                                       ke_change(vor.unfiltered, div.unfiltered, vor.filtered_new .- vor.unfiltered,
                                                 div.filtered_new .- div.unfiltered)) ./ (2 * half_step)
                    @test maximum(abs, filter_expected) > 0
                    @test maximum(abs.(filter_change .- filter_expected)) < 0.02 * maximum(abs, filter_expected)
                    println("DISSIPATIVE_HEATING_FILTER mode=$mode layer_W_kg=", filter_change)
                    heating = copy(spectral.temperature.data[1, :])
                    expected = -ke_change(vorticity, divergence, dt .* Tv, dt .* Td) ./ (2 * half_step)
                    @test all(>(0), loss)
                    # The grid winds are U / cos(latitude), so Gaussian quadrature of their squares is not exact: at
                    # T7 it differs from the exact spectral energy by about 0.3 % (0.1 % at T31).
                    @test relative(loss, expected; floor = 1e-9) < 1e-2
                    if heats
                        # One l = 0 coefficient per layer: a uniform heating of loss / cp (in the scaled time).
                        uniform = zeros(eltype(spectral.temperature.data), sg.spectrum, n)
                        uniform.data[1, :] .= heating
                        field = Float64.(Array(S.transform(uniform, model.spectral_transform).data))
                        target = loss ./ cp .* (half_step / Float64(model.time_stepping.Δt))
                        @test maximum(abs.(field .- target') ./ abs.(target')) < 1e-4
                    else
                        @test all(iszero, heating)
                    end
                    column = sum(loss .* dsigma) * sum(weights .* ps) / gravity
                    println("DISSIPATIVE_HEATING_HYPERDIFFUSION mode=$mode global_W_m2=$column worst_relative=",
                            relative(loss, expected; floor = 1e-9))
                end

                @testset "speed-limit drag" begin
                    drag = model.drag
                    @test drag isa R.DissipativeSpeedLimitDrag && drag.heat == heats
                    # Winds 1 m/s above each limit: 61 m/s aloft (all-level limit 60) and 41 m/s in the lowest layer
                    # (surface limit 40), with no v, so the step changes the wind only slightly.
                    fill!(vars.grid.u.data, 61); fill!(vars.grid.v.data, 0)
                    fill!(vars.grid.u_prev.data, 60.5); fill!(vars.grid.v_prev.data, 0)
                    vars.grid.u.data[:, n] .= 41
                    vars.grid.u_prev.data[:, n] .= 40.5
                    zero_tendencies!(vars)
                    S.drag!(vars, drag, 2, model)
                    du, dv, dT = copy(tend.u.data), copy(tend.v.data), copy(tend.temperature.data)
                    D = Float64.(Array(p.speed_limit_drag_dissipation.data))
                    zero_tendencies!(vars)
                    S.drag!(vars, drag.drag, 2, model)
                    @test same_bits(du, tend.u.data) && same_bits(dv, tend.v.data)
                    scale = Float64(vars.prognostic.scale[])
                    Fu, Fv = Float64.(du) ./ scale, Float64.(dv) ./ scale
                    up = Float64.(Array(vars.grid.u_prev.data))
                    # The hyperdiffusion check's transform!(...; initialize = true) copied the current state into the
                    # _prev fields, including surface pressure, so take the layer mass from the state used here.
                    mass_now = Float64.(Array(vars.grid.pressure_prev.data)) ./ gravity
                    @test all(iszero, Fv)
                    @test all(<(0), Fu[:, 1:(n - 1)]) && all(<(0), Fu[:, n])   # every layer is above its limit
                    loss = .-((up .+ half_step .* Fu) .* Fu)
                    expected = mass_now .* vec(sum(loss .* dsigma'; dims = 2))
                    @test all(>(0), D)
                    @test relative(D, expected) < 1e-4
                    if heats
                        @test relative(cp .* Float64.(dT) ./ scale, loss) < 1e-4
                    else
                        @test all(iszero, dT)
                    end
                end
            end
        end
    end
end

end
