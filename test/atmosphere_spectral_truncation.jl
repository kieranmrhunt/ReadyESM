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
    @test production.atmosphere_initial_spectral_truncation == false            # opt-in

    # ERA5 launch exactly as ReadyESM builds it (T31/L27, CPU), without the coupled components
    config = with(production; atmosphere_initial_conditions = :era5_instantaneous)
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

    @testset "restored-state cleaning recomputes a consistent grid" begin
        SW.run!(sim; steps = 3)                                               # evolve, then contaminate the top row
        T2 = step_of(vars.prognostic.temperature, 2)
        l, m = lm_indices(T2); i = findfirst((l .== maximum(l)) .& (m .== 1))
        T2.data[i, 1] += 1.0f0                                                  # l = lmax + 1, m = 0, top layer
        @test top_rows(T2)[1] > 0
        fake = (; model = (; atmosphere = sim))
        R.truncate_restored_atmosphere_lmax_plus_1!(fake)
        @test top_rows(step_of(vars.prognostic.temperature, 2))[1] == 0
        @test Array(vars.grid.temperature.data) == Array(vars.grid.temperature_prev.data)
        regrid = ST.transform(step_of(vars.prognostic.temperature, 2), model.spectral_transform)
        @test isapprox(Array(regrid.data), Array(vars.grid.temperature.data); rtol = 1e-6)
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
end
end
