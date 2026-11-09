module LargeScaleMeltCoolingTests
using ReadyESM, Test, Dates
const R = ReadyESM
const S = R.SpeedyWeather
const O = R.Oceananigans

fields(value) = (; (k => getfield(value, k) for k in fieldnames(typeof(value)))...)
with(config; kwargs...) = R.ExperimentConfig(; merge(fields(config), (; kwargs...))...)
include(joinpath(@__DIR__, "analytic_coupled_fixture.jl"))

# Column heat input of the large-scale scheme (W/m²): the sum over layers of (cₚ dT/dt + Lᵥ dq/dt) times the layer mass.
function column_heat(vars, model, ij)
    atmosphere, Δσ = model.atmosphere, model.geometry.σ_levels_thick
    mass = Float64(vars.grid.pressure_prev[ij]) / Float64(model.planet.gravity)
    return sum((Float64(atmosphere.heat_capacity) * Float64(vars.tendencies.grid.temperature[ij, k]) +
                Float64(atmosphere.latent_heat_condensation) * Float64(vars.tendencies.grid.humidity[ij, k])) *
               Float64(Δσ[k]) * mass for k in eachindex(Δσ))
end

latent_fusion_flux(model, rate) =
    Float64(model.atmosphere.latent_heat_fusion) * Float64(model.atmosphere.water_density) * Float64(rate)

function reset_column!(vars)
    fill!(vars.tendencies.grid.temperature.data, 0)
    fill!(vars.tendencies.grid.humidity.data, 0)
    p = vars.parameterizations
    for name in (:rain_large_scale, :snow_large_scale, :rain_rate_large_scale, :snow_rate_large_scale,
                 :rain_rate, :snow_rate, :large_scale_snowmelt_rate)
        fill!(getproperty(p, name).data, 0)
    end
end

# Upper layers cold and supersaturated (large-scale snow forms), lower layers at `warm` K and moist but unsaturated.
function set_column!(vars, model, ij, warm)
    σ = model.geometry.σ_levels_full
    ps = 1.0f5
    vars.grid.pressure_prev[ij] = ps
    for k in eachindex(σ)
        cold = σ[k] < 0.5f0
        T = cold ? 250.0f0 : Float32(warm)
        vars.grid.temperature_prev[ij, k] = T
        qsat = S.saturation_humidity(T, σ[k] * ps, model.atmosphere)
        vars.grid.humidity_prev[ij, k] = (cold ? 1.2f0 : 0.8f0) * qsat
    end
end

function run_column!(vars, model, scheme, ij)
    reset_column!(vars)
    S.parameterization!(ij, vars, scheme, model)
    p = vars.parameterizations
    return (heat = column_heat(vars, model, ij), snow = p.snow_rate_large_scale[ij],
            rain = p.rain_rate_large_scale[ij], melt = p.large_scale_snowmelt_rate[ij],
            dT = copy(vars.tendencies.grid.temperature.data[ij, :]),
            dq = copy(vars.tendencies.grid.humidity.data[ij, :]))
end

@testset "Large-scale snowmelt cooling (opt-in)" begin
    default = R.ExperimentConfig()
    @test default.atmosphere_large_scale_melt_cooling == false
    base = fields(default)
    on = R.ExperimentConfig(; merge(base, (; atmosphere_large_scale_precipitation = :net_column_conservative,
        atmosphere_large_scale_melt_cooling = true))...)
    @test R.validate(on; check_input_files = false) === on
    bad = R.ExperimentConfig(; merge(base, (; atmosphere_large_scale_precipitation = :upstream_observed,
        atmosphere_large_scale_melt_cooling = true))...)
    @test_throws ArgumentError R.validate(bad; check_input_files = false)
    off = R.ExperimentConfig(; merge(base, (; atmosphere_large_scale_precipitation = :net_column_conservative))...)
    grid = S.SpectralGrid(; NF = Float32, trunc = 7, nlayers = 8)
    @test R._atmosphere_large_scale_condensation(on, grid) isa
        R.NetColumnImplicitCondensation{<:R.MeltCoolingImplicitCondensation{<:S.ImplicitCondensation}}
    # Off: the same object type and parameters as the release.
    off_scheme = R._atmosphere_large_scale_condensation(off, grid)
    @test off_scheme isa R.NetColumnImplicitCondensation{<:S.ImplicitCondensation}
    @test typeof(off_scheme) == typeof(R.NetColumnImplicitCondensation(grid))
    @test fields(off_scheme.condensation) == fields(S.ImplicitCondensation(grid))
    @test fields(R._atmosphere_large_scale_condensation(on, grid).condensation.condensation) ==
        fields(S.ImplicitCondensation(grid))

    # The YAML key is read.
    mktempdir() do directory
        source = joinpath(R.PROJECT_ROOT, "config", "production.yml")
        path = joinpath(directory, "melt.yml")
        @test R.load_config(source; check_input_files = false).atmosphere_large_scale_melt_cooling
        original = read(source, String)
        @test occursin("atmosphere_large_scale_melt_cooling: true", original)
        write(path, replace(original, "atmosphere_large_scale_melt_cooling: true" =>
            "atmosphere_large_scale_melt_cooling: false"))
        @test !R.load_config(path; check_input_files = false).atmosphere_large_scale_melt_cooling
    end

    # Columns on a real CPU model with the option on, so the snowmelt output exists.
    config = R.ExperimentConfig(; name = "melt_cooling_column", device = :cpu, truncation = 7, nlayers = 8,
        duration_days = 0.1, mixed_layer_depth_m = 30., atmosphere_large_scale_precipitation = :net_column_conservative,
        atmosphere_large_scale_melt_cooling = true)
    simulation, _ = R.run_baseline!(config; steps = 2)
    vars, model = simulation.variables, simulation.model
    wrapper = model.large_scale_condensation
    @test wrapper isa R.NetColumnImplicitCondensation{<:R.MeltCoolingImplicitCondensation}
    fixed = wrapper.condensation
    upstream = fixed.condensation
    ij = 1

    for (label, warm) in (("full melt", 285.0), ("partial melt", 278.02))
        set_column!(vars, model, ij, warm)
        a = run_column!(vars, model, fixed, ij)
        set_column!(vars, model, ij, warm)
        b = run_column!(vars, model, upstream, ij)
        scale = latent_fusion_flux(model, a.snow + a.melt)
        @test a.melt > 0
        @test scale > 1  # W/m²: a substantial snow signal
        # Fixed: the heat input equals Lᵢ ρ_w times the snow leaving the column (gross, before the wrapper).
        @test isapprox(a.heat, latent_fusion_flux(model, a.snow); atol = 1e-4 * scale)
        # Upstream: the excess equals Lᵢ ρ_w times the melted snow (negative control: the defect is present).
        @test b.snow == a.snow && b.melt == 0  # upstream does not write the melt output; reset left it at zero
        @test isapprox(b.heat - latent_fusion_flux(model, b.snow), latent_fusion_flux(model, a.melt);
            atol = 1e-3 * scale)
        @test b.heat - latent_fusion_flux(model, b.snow) > 0.5 * latent_fusion_flux(model, a.melt)
        println("MELT_COOLING_COLUMN case=$(replace(label, ' ' => '_')) snow_Wm2=$(latent_fusion_flux(model, a.snow)) ",
            "melt_Wm2=$(latent_fusion_flux(model, a.melt)) fixed_residual=$(a.heat - latent_fusion_flux(model, a.snow)) ",
            "upstream_excess=$(b.heat - latent_fusion_flux(model, b.snow))")
    end

    # No melting (all layers cold): bit-identical tendencies and rates.
    set_column!(vars, model, ij, 250.0)
    a = run_column!(vars, model, fixed, ij)
    set_column!(vars, model, ij, 250.0)
    b = run_column!(vars, model, upstream, ij)
    @test a.melt == 0 && a.snow > 0
    @test a.dT == b.dT && a.dq == b.dq && a.snow == b.snow && a.rain == b.rain

    # Through the production wrapper: it runs, records the melt, and caps the surface snow at the net-column total.
    set_column!(vars, model, ij, 278.02)
    reset_column!(vars)
    S.parameterization!(ij, vars, wrapper, model)
    p = vars.parameterizations
    @test p.large_scale_snowmelt_rate[ij] > 0
    @test 0 <= p.snow_rate_large_scale[ij]
    @test all(isfinite, vars.tendencies.grid.temperature.data[ij, :])
end

@testset "Large-scale snowmelt cooling: coupled CPU build and run, off and on" begin
    released = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)
    @test released.atmosphere_large_scale_precipitation == :net_column_conservative
    states = Dict{Bool, Any}()
    for melt in (false, true)
        base = analytic_coupled_test_config(released; output_dir = mktempdir())
        config = with(base; atmosphere_timestep_at_t31_minutes = 3.75, duration_days = 0.25,
            start_date = DateTime(1993, 1, 1), atmosphere_large_scale_melt_cooling = melt)
        R.validate(config; check_input_files = false)
        simulation = R.build_dynamic_esm(config)
        lsc = simulation.model.atmosphere.model.large_scale_condensation
        @test lsc isa R.NetColumnImplicitCondensation
        @test (lsc.condensation isa R.MeltCoolingImplicitCondensation) == melt
        O.run!(simulation)
        vars = simulation.model.atmosphere.variables
        T, q = Array(vars.grid.temperature_prev.data), Array(vars.grid.humidity_prev.data)
        @test all(isfinite, T) && all(isfinite, q)
        if melt
            rate = Array(vars.parameterizations.large_scale_snowmelt_rate.data)
            @test all(x -> isfinite(x) && x >= 0, rate)
            println("MELT_COOLING_COUPLED max_snowmelt_rate_m_s=$(maximum(rate)) columns_melting=$(count(>(0), rate))")
        else
            @test !hasproperty(vars.parameterizations, :large_scale_snowmelt_rate)
        end
        states[melt] = (T, q)
    end
    println("MELT_COOLING_COUPLED atmosphere_identical_off_on=$(states[false] == states[true])")
end
# Reached only if both top-level test sets passed (a failing set throws when it ends).
println("LARGE_SCALE_MELT_COOLING_TESTS_PASS")
end
