module AtmosphereSpectralTruncationTests
using ReadyESM, Test, Dates
const R = ReadyESM
const SW = R.SpeedyWeather
const ST = SW.SpeedyTransforms

with(config; kwargs...) = R.ExperimentConfig(; merge(
    (; (n => getfield(config, n) for n in fieldnames(typeof(config)))...), (; kwargs...))...)

lm_indices(x) = (Array(x.spectrum.l_indices), Array(x.spectrum.m_indices))
"""Power in the l = lmax + 1 row and in l = lmax of a spectral field level (any number of layers)."""
function top_rows(x)
    l, m = lm_indices(x); d = reshape(Array(x.data), length(l), :); w = ifelse.(m .== 1, 1, 2)
    pw(L) = sum(abs2.(d[l .== L, :]) .* w[l .== L])
    return pw(maximum(l)), pw(maximum(l) - 1)
end
spectra(vars) = Dict(name => [copy(Array(SW.get_step(getproperty(vars.prognostic, name), lf).data)) for lf in 1:2]
                     for name in (:vorticity, :divergence, :temperature, :humidity, :pressure))
step_of(var, lf) = SW.get_step(var, lf)

@testset "atmosphere initial spectral truncation" begin
    production = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"))
    @test production.atmosphere_initial_spectral_truncation                     # selected in production
    @test !R.ExperimentConfig().atmosphere_initial_spectral_truncation           # generic compatibility default

    # ERA5 launch exactly as ReadyESM builds it (T31/L27, CPU), without the coupled components
    config = with(production; atmosphere_initial_conditions = :era5_instantaneous,
        atmosphere_initial_spectral_truncation = false)                         # retain the unfixed reproducer
    sg = SW.SpectralGrid(; NF = Float32, trunc = config.truncation, Grid = SW.RingGrids.FullGaussianGrid,
                         nlayers = config.nlayers)
    model = SW.PrimitiveWetModel(sg; initial_conditions = R._atmosphere_initial_conditions(config, sg))
    sim = SW.initialize!(model; time = config.start_date)
    vars = sim.variables

    @testset "the ERA5 launch leaves l = lmax + 1 content in ln ps, T and q, and it stays frozen (the bug)" begin
        @test top_rows(step_of(vars.prognostic.pressure, 1))[1] > 0
        @test top_rows(step_of(vars.prognostic.temperature, 1))[1] > 0
        @test top_rows(step_of(vars.prognostic.humidity, 1))[1] > 0
        @test top_rows(step_of(vars.prognostic.vorticity, 1))[1] == 0          # curl! output is truncated
        # a separate launch, integrated: the l = lmax + 1 row of ln ps does not change (only tendencies are
        # truncated), while the l = lmax row is damped
        sim2 = SW.initialize!(SW.PrimitiveWetModel(sg; initial_conditions = R._atmosphere_initial_conditions(config, sg),
                                  time_stepping = SW.Leapfrog(sg; Δt_at_T31 = SW.Minute(15)),
                                  horizontal_diffusion = R._atmosphere_horizontal_diffusion(config, sg));
                              time = config.start_date)
        p0 = top_rows(step_of(sim2.variables.prognostic.pressure, 1))
        l, _ = lm_indices(step_of(sim2.variables.prognostic.pressure, 1)); top = l .== maximum(l)
        c0 = Array(step_of(sim2.variables.prognostic.pressure, 1).data)[top]
        SW.run!(sim2; steps = 96)
        c1 = Array(step_of(sim2.variables.prognostic.pressure, 2).data)[top]
        p1 = top_rows(step_of(sim2.variables.prognostic.pressure, 2))
        @test all(isfinite, Array(step_of(sim2.variables.prognostic.pressure, 2).data))
        @test c1 == c0                                                          # frozen, bitwise
        @test p1[2] < p0[2] / 10                                                # l = lmax damped within a day
    end

    before = spectra(vars)
    R.truncate_prognostic_lmax_plus_1!(vars)
    after = spectra(vars)
    @testset "truncation zeroes exactly the l = lmax + 1 row, bitwise elsewhere, idempotent" begin
        l, _ = lm_indices(step_of(vars.prognostic.pressure, 1)); top = l .== maximum(l)
        for (name, levels) in after, lf in 1:2
            A = reshape(levels[lf], length(l), :); Bm = reshape(before[name][lf], length(l), :)
            @test all(iszero, A[top, :])
            @test A[.!top, :] == Bm[.!top, :]
        end
        R.truncate_prognostic_lmax_plus_1!(vars)
        @test spectra(vars) == after
    end

    current_names = (:vorticity, :divergence, :temperature, :pressure, :humidity, :u, :v)
    previous_names = (:u_prev, :v_prev, :temperature_prev, :pressure_prev, :humidity_prev)
    grab(names) = Dict(n => copy(Array(getproperty(vars.grid, n).data)) for n in names)
    @testset "restored-state cleaning: transform-only identity, row removal, previous grids preserved" begin
        # step as the coupled model does (initialize!, then timestep! without finalize), so vorticity and divergence
        # stay radius-scaled; the current grid is then the transform of level 2
        SW.initialize!(sim; steps = 3)
        for _ in 1:3; SW.timestep!(sim); end
        @test vars.prognostic.scale[] == model.planet.radius
        fake = (; model = (; atmosphere = sim))
        # identity: on an already-truncated state the recompute reproduces the current grid bitwise and, with the
        # default preserve_previous, leaves every *_prev grid bitwise unchanged
        current0 = grab(current_names); previous0 = grab(previous_names)
        @test all(n -> previous0[n] != current0[n === :pressure_prev ? :pressure : Symbol(replace(String(n), "_prev" => ""))],
                  (:u_prev, :temperature_prev, :humidity_prev))                 # prev really is a different (older) state
        R.truncate_restored_atmosphere_lmax_plus_1!(fake)
        @test grab(current_names) == current0
        @test grab(previous_names) == previous0
        # contaminate the top row, then clean: row zeroed, current grid recomputed from it, *_prev untouched
        T2 = step_of(vars.prognostic.temperature, 2)
        l, m = lm_indices(T2); i = findfirst((l .== maximum(l)) .& (m .== 1))
        T2.data[i, 1] += 1.0f0                                                  # l = lmax + 1, m = 0, top layer
        @test top_rows(T2)[1] > 0
        R.truncate_restored_atmosphere_lmax_plus_1!(fake)
        @test top_rows(step_of(vars.prognostic.temperature, 2))[1] == 0
        @test grab(previous_names) == previous0
        @test grab(current_names) == current0                                   # the added row is gone again
        regrid = ST.transform(step_of(vars.prognostic.temperature, 2), model.spectral_transform)
        @test isapprox(Array(regrid.data), Array(vars.grid.temperature.data); rtol = 1e-6)
        # preserve_previous = false keeps the initialisation semantics (*_prev = current)
        R.truncate_restored_atmosphere_lmax_plus_1!(fake; preserve_previous = false)
        @test Array(vars.grid.temperature.data) == Array(vars.grid.temperature_prev.data)
        @test Array(vars.grid.humidity.data) == Array(vars.grid.humidity_prev.data)
        # an unscaled state (after run!'s finalize) would give u/radius: refused
        SW.unscale!(vars)
        @test_throws ArgumentError R.truncate_restored_atmosphere_lmax_plus_1!(fake)
    end

    @testset "builder: the flag zeroes only the l = lmax + 1 row of the launch state" begin
        small = R.ExperimentConfig(device = :cpu, truncation = 7, nlayers = 8,
            forcing = R.ForcingConfig(radiation = :rrtmgp_all_sky))
        a = R.build_rrtmgp(small)                                               # default: unchanged launch
        b = R.build_rrtmgp(with(small; atmosphere_initial_spectral_truncation = true))
        lb, _ = lm_indices(step_of(b.variables.prognostic.pressure, 1)); top = lb .== maximum(lb)
        sa = spectra(a.variables); sb = spectra(b.variables)
        for name in keys(sb), lf in 1:2
            A = reshape(sa[name][lf], length(lb), :); Bm = reshape(sb[name][lf], length(lb), :)
            @test all(iszero, Bm[top, :])
            @test A[.!top, :] == Bm[.!top, :]
        end
        # SpeedyWeather's own analytic launch also leaves l = lmax + 1 temperature content (only ln ps is truncated
        # there), so the default path differs from the truncated one exactly in that row
        @test any(!iszero, reshape(sa[:temperature][1], length(lb), :)[top, :])
    end

    # Contract on the ordering the launch flag relies on (Codex review 2026-10-14 17:38): the builder truncates the
    # spectra, and SpeedyWeather's simulation initialize! (scale, then transform! with initialize = true) rebuilds the
    # first grid state from them. After it, for both flag values: the current T, q, ln ps, u, v grids equal fresh
    # transforms of the selected leapfrog level's spectra, and every *_prev grid equals the initialised current grid.
    @testset "post-initialisation grids are the transforms of the (truncated) spectra; *_prev = current" begin
        small = R.ExperimentConfig(device = :cpu, truncation = 7, nlayers = 8, atmosphere_timestep_at_t31_minutes = 3.75,
            atmosphere_convection_timescale_seconds = 7200, forcing = R.ForcingConfig(radiation = :rrtmgp_all_sky))
        grids = Dict{Bool, Dict{Symbol, Array}}()
        for flag in (false, true)
            sim = R.build_rrtmgp(with(small; atmosphere_initial_spectral_truncation = flag))
            model = sim.model; v = sim.variables; S = model.spectral_transform
            @test model.time_stepping.Δt_sec == 900
            SW.initialize!(sim; steps = 4)
            @test v.prognostic.scale[] == model.planet.radius
            lf = model.time_stepping.first_step_euler ? 1 : 2
            P = v.prognostic; G = v.grid
            for name in (:temperature, :humidity, :pressure)
                spec = step_of(getproperty(P, name), lf)
                flag && @test top_rows(spec)[1] == 0                            # the first grid state comes from clean spectra
                expected = ST.transform(spec, S)
                name === :humidity && SW.hole_filling!(expected, model.hole_filling, model)   # as transform! does
                @test isapprox(Array(getproperty(G, name).data), Array(expected.data); rtol = 1e-6, atol = 1e-12)
            end
            vor = step_of(P.vorticity, lf); div = step_of(P.divergence, lf)
            U = similar(vor); V = similar(vor)
            ST.UV_from_vordiv!(U, V, vor, div, S)
            @test isapprox(Array(G.u.data), Array(ST.transform(U, S; unscale_coslat = true).data); rtol = 1e-6, atol = 1e-6)
            @test isapprox(Array(G.v.data), Array(ST.transform(V, S; unscale_coslat = true).data); rtol = 1e-6, atol = 1e-6)
            @test Array(G.u_prev.data) == Array(G.u.data) && Array(G.v_prev.data) == Array(G.v.data)
            @test Array(G.temperature_prev.data) == Array(G.temperature.data)
            @test Array(G.humidity_prev.data) == Array(G.humidity.data)
            @test Array(G.pressure_prev.data) == exp.(Array(G.pressure.data))       # stored in Pa
            grids[flag] = Dict(n => copy(Array(getproperty(G, n).data)) for n in (:temperature, :humidity, :pressure))
            for _ in 1:4
                SW.timestep!(sim)
                @test all(all(isfinite, Array(getproperty(G, n).data)) for n in (:temperature, :humidity, :pressure, :u, :v))
            end
        end
        # the flag changes the first grid state (the imprint is not only in unused spectra)
        @test grids[true][:temperature] != grids[false][:temperature]
    end
end
end
