module MoistMassCouplingTests
using ReadyESM, Test, Random
const R = ReadyESM
const S = R.SpeedyWeather
include(joinpath(@__DIR__, "moist_mass_column_reference.jl"))
using .MoistMassColumn

fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
with(config; kwargs...) = R.ExperimentConfig(; merge(fields(config), (; kwargs...))...)
with_forcing(config; kwargs...) =
    with(config; forcing = R.ForcingConfig(; merge(fields(config.forcing), (; kwargs...))...))

function mock(rng, npts, K; E_sign = 1)
    Δσ = rand(rng, K) .+ 0.2
    Δσ ./= sum(Δσ)
    g, a = 9.81f0, 6.371f6
    ps = Float32.(9.5e4 .+ 1e4 .* rand(rng, npts))
    q = Float32.(1e-3 .+ 1.5e-2 .* rand(rng, npts, K))
    T = Float32.(200 .+ 100 .* rand(rng, npts, K))
    u = Float32.(20 .* randn(rng, npts, K))
    v = Float32.(20 .* randn(rng, npts, K))
    E = Float32.(E_sign .* 1e-4 .* rand(rng, npts))
    lsc = Float32.(1e-7 .* randn(rng, npts, K))                       # condensation and re-evaporation, kg/kg/s
    conv = Float32.(1e-7 .* randn(rng, npts, K))                      # Betts-Miller: net source plus redistribution
    tq = lsc .+ conv
    for ij in 1:npts
        tq[ij, K] += g * E[ij] / (ps[ij] * Float32(Δσ[K]))           # surface_flux_to_tendency
    end
    zeros2 = () -> zeros(Float32, npts)
    vars = (; parameterizations = (; large_scale_condensation_humidity_tendency = lsc, surface_humidity_flux = E,
                mmc_column_source = zeros2(), mmc_large_scale_source = zeros2(), mmc_surface_source = zeros2(),
                mmc_remainder_source = zeros2(), mmc_source_enthalpy = zeros2()),
            grid = (; pressure_prev = ps, humidity_prev = q, temperature_prev = T, u_prev = u, v_prev = v),
            tendencies = (; grid = (; humidity = copy(tq), temperature = Float32.(1e-5 .* randn(rng, npts, K)),
                u = Float32.(1e-4 .* randn(rng, npts, K)), v = Float32.(1e-4 .* randn(rng, npts, K)), pressure = zeros2())))
    model = (; geometry = (; σ_levels_thick = Float32.(Δσ)), planet = (; gravity = g, radius = a),
             atmosphere = (; heat_capacity = 1004f0, latent_heat_condensation = 2.501f6))
    return vars, model, (; Δσ = Float32.(Δσ), conv)
end
snapshot(vars) = deepcopy(vars.tendencies.grid)

"Text of one top-level method starting at `header` and ending at the first line that is exactly `end`."
function method_text(source, header)
    start = findfirst(header, source)
    isnothing(start) && error("missing method header: $header")
    count(header, source) == 1 || error("method header not unique: $header")
    stop = findnext("\nend\n", source, first(start))
    return source[first(start):last(stop)-1]
end

copy_text(source, name) = begin
    a = findfirst("# BEGIN MOIST MASS COPY: $name\n", source)
    b = findfirst("# END MOIST MASS COPY: $name", source)
    body = source[last(a)+1:first(b)-1]
    body[findfirst("function ", body)[1]:end]
end

@testset "Moist-mass coupling option" begin
    production_path = joinpath(R.PROJECT_ROOT, "config", "production.yml")
    production = load_config(production_path; check_input_files = false)
    @test R.ExperimentConfig().atmosphere_moist_mass_coupling == false
    @test production.atmosphere_moist_mass_coupling == false
    raw = R.YAML.load_file(production_path)
    raw["atmosphere_moist_mass_coupling"] = true
    path = tempname() * ".yml"
    R.YAML.write_file(path, raw)
    coupled = load_config(path; check_input_files = false)
    @test coupled.atmosphere_moist_mass_coupling == true
    @test R.validate(coupled; check_input_files = false) === coupled
    @test R.validate(with(coupled; atmosphere_convection = :betts_miller_sigma_rh_v1); check_input_files = false) isa
          R.ExperimentConfig
    for bad in (with(coupled; device = :cpu, atmosphere_convection = :finite_multiplume_c08_v1),
                with(coupled; device = :cpu, atmosphere_convection = :finite_multiplume_c07_v1),
                with(coupled; atmosphere_large_scale_precipitation = :upstream_observed),
                with(coupled; atmosphere_vertical_diffusion = :metric_bulk_richardson),
                with_forcing(coupled; radiation = :speedy_simplified),
                with_forcing(coupled; cloud_scheme = :prognostic_condensate))
        @test_throws ArgumentError R.validate(bad; check_input_files = false)
    end
    @test R._atmosphere_stochastic_physics(production) === nothing
    @test R._atmosphere_stochastic_physics(coupled) isa R.MoistMassCoupling
    @test R._moist_mass_coupling_provenance(production) == (;)
    provenance = R._moist_mass_coupling_provenance(coupled)
    @test provenance.requested_atmosphere_moist_mass_coupling === true
    @test provenance.atmosphere_moist_mass_coupling_scheme == "betts_miller_sigma_spread_remainder_v1"
    @test provenance.atmosphere_pressure_zero_mode_resets == "removed_surface_and_implicit_v1"
end

@testset "Moist-mass coupling column scheme against the column reference" begin
    rng = Xoshiro(20261021)
    for E_sign in (1, -1), trial in 1:3
        npts, K = 7, 27
        vars, model, extra = mock(rng, npts, K; E_sign)
        before = snapshot(vars)
        for ij in 1:npts
            S.parameterization!(ij, vars, R.MoistMassCoupling(), model)
        end
        after = vars.tendencies.grid
        g = Float64(model.planet.gravity); a = Float64(model.planet.radius)
        Δσ = Float64.(extra.Δσ)
        for ij in 1:npts
            ps = Float64(vars.grid.pressure_prev[ij])
            m = Δσ .* (ps / g)
            E = Float64(vars.parameterizations.surface_humidity_flux[ij])
            L = m .* Float64.(vars.parameterizations.large_scale_condensation_humidity_tendency[ij, :])
            Cref = sum(m .* Float64.(extra.conv[ij, :]))
            Σref = sum(m .* Float64.(before.humidity[ij, :]))
            p = vars.parameterizations
            @test isapprox(p.mmc_column_source[ij], Σref; rtol = 1e-5, atol = 1e-12)
            @test isapprox(p.mmc_remainder_source[ij], Cref; rtol = 1e-3, atol = 1e-9)
            @test isapprox(p.mmc_large_scale_source[ij] + p.mmc_surface_source[ij] + p.mmc_remainder_source[ij],
                           p.mmc_column_source[ij]; rtol = 1e-5, atol = 1e-12)
            Sk = L .+ Δσ .* Cref .+ [k == length(Δσ) ? E : 0.0 for k in eachindex(Δσ)]
            Href = sum(Sk .* (1004.0 .* Float64.(vars.grid.temperature_prev[ij, :]) .+ Float64(2.501f6)))
            @test isapprox(p.mmc_source_enthalpy[ij], Href; rtol = 1e-3, atol = 1e-3 * sum(abs.(Sk)) * 2.8e6)
            Σ = Float64(p.mmc_column_source[ij])
            nonE = L .+ Δσ .* (Σ - sum(L) - E)
            S_add = max.(nonE, 0.0); S_remove = max.(-nonE, 0.0)
            E > 0 ? (S_add[end] += E) : (S_remove[end] -= E)
            for (component, ucol, du) in (("u", vars.grid.u_prev, after.u .- before.u), ("v", vars.grid.v_prev, after.v .- before.v))
                s = ColumnState(ps, Float64.(vars.grid.humidity_prev[ij, :]), Float64.(vars.grid.temperature_prev[ij, :]),
                                Float64.(ucol[ij, :]))
                u_add = copy(s.u)
                E > 0 && (u_add[end] = ((S_add[end] - E) * s.u[end]) / S_add[end])          # evaporated vapour at rest
                ref = mass_coupling_correction(s, S_add, S_remove, Δσ, g; T_add = copy(s.T), u_add)
                scale(x) = 1e-5 * maximum(abs, x) + 1e-12
                @test maximum(abs, du[ij, :] .- ref.du) <= scale(ref.du) + 4 * maximum(abs, before.u) * eps(Float32)
                if component == "u"
                    dq = after.humidity[ij, :] .- before.humidity[ij, :]
                    dT = after.temperature[ij, :] .- before.temperature[ij, :]
                    @test maximum(abs, dq .- ref.dq) <= scale(ref.dq) + 4 * maximum(abs, before.humidity) * eps(Float32)
                    @test maximum(abs, dT .- ref.dT) <= scale(ref.dT) + 4 * maximum(abs, before.temperature) * eps(Float32)
                    @test isapprox(after.pressure[ij], -a * ref.dlnps; rtol = 1e-5)
                    # The Float32 layer thicknesses do not sum to exactly 1, so sum(S) - Σ = (sum(Δσ) - 1) C exactly.
                    Sref = sum(S_add .- S_remove)
                    @test isapprox(ref.dlnps, g * Sref / ps; rtol = 1e-12)
                    @test abs((Sref - Σ) - (sum(Δσ) - 1) * (Σ - sum(L) - E)) <= 1e-12 * (abs(Σ) + sum(abs.(nonE)) + abs(E))
                end
            end
        end
    end
    # The remainder is spread in sigma, so moving the Betts-Miller net source between layers (same column total)
    # changes only the dilution, never the T, u or v tendencies.
    vars1, model, extra = mock(Xoshiro(7), 3, 27)
    vars2 = deepcopy(vars1)
    for ij in 1:3
        vars2.tendencies.grid.humidity[ij, 3] += 1f-7 / extra.Δσ[3]
        vars2.tendencies.grid.humidity[ij, 20] -= 1f-7 / extra.Δσ[20]
    end
    b1, b2 = snapshot(vars1), snapshot(vars2)
    for ij in 1:3
        S.parameterization!(ij, vars1, R.MoistMassCoupling(), model)
        S.parameterization!(ij, vars2, R.MoistMassCoupling(), model)
    end
    for name in (:temperature, :u, :v)
        d1 = getproperty(vars1.tendencies.grid, name) .- getproperty(b1, name)
        d2 = getproperty(vars2.tendencies.grid, name) .- getproperty(b2, name)
        @test maximum(abs, d1 .- d2) <= 1e-3 * maximum(abs, d1) + 1e-12
    end
    # Zero sources give zero corrections, bitwise.
    vars, model, _ = mock(Xoshiro(13), 2, 27)
    vars.tendencies.grid.humidity .= 0
    vars.parameterizations.large_scale_condensation_humidity_tendency .= 0
    vars.parameterizations.surface_humidity_flux .= 0
    before = snapshot(vars)
    for ij in 1:2
        S.parameterization!(ij, vars, R.MoistMassCoupling(), model)
    end
    for name in (:humidity, :temperature, :u, :v, :pressure)
        @test getproperty(vars.tendencies.grid, name) == getproperty(before, name)
    end
end

@testset "Moist-mass coupling copies of the upstream pressure methods" begin
    sw = dirname(pathof(S))
    tendencies = read(joinpath(sw, "dynamics", "tendencies.jl"), String)
    implicit = read(joinpath(sw, "time_stepping", "implicit.jl"), String)
    ours = read(joinpath(R.PROJECT_ROOT, "src", "moist_mass_coupling.jl"), String)
    strip_qualifier(text) = replace(text, "SpeedyWeather." => "")
    surface_reset = "    pres_tend.data[1:1] .= 0                # for mass conservation\n"
    implicit_reset = "    pres_tend.data[1:1] .= 0    # mass conservation\n\n"
    # Source-scaling path and the two reset removals are separate checks.
    up_surface = method_text(tendencies, "function surface_pressure_tendency!(\n")
    @test count(surface_reset, up_surface) == 1
    expected_surface = replace(up_surface, surface_reset => "",
        "function surface_pressure_tendency!(" => "function _moist_mass_surface_pressure_tendency!(")
    @test strip_qualifier(copy_text(ours, "surface_pressure_tendency!")) == expected_surface * "\n"
    up_dynamics = method_text(tendencies, "function dynamics_tendencies!(\n        vars::Variables,\n        " *
        "lf::Integer,                # leapfrog index for tendencies\n        model::PrimitiveEquation,")
    @test count("    surface_pressure_tendency!(vars, spectral_transform)\n", up_dynamics) == 1
    expected_dynamics = replace(up_dynamics,
        "        model::PrimitiveEquation,\n" => "        model::MoistMassCoupledModel,\n",
        "    surface_pressure_tendency!(vars, spectral_transform)\n" =>
            "    _moist_mass_surface_pressure_tendency!(vars, spectral_transform)\n")
    @test strip_qualifier(copy_text(ours, "dynamics_tendencies!")) == expected_dynamics * "\n"
    up_implicit = method_text(implicit, "function implicit_correction!(\n        vars::Variables,\n        " *
        "implicit::ImplicitPrimitiveEquation,\n        model::PrimitiveEquation,")
    @test count(implicit_reset, up_implicit) == 1
    expected_implicit = replace(up_implicit,
        "        model::PrimitiveEquation,\n" => "        model::MoistMassCoupledModel,\n",
        implicit_reset => "")
    @test strip_qualifier(copy_text(ours, "implicit_correction!")) == expected_implicit * "\n"
    for name in ("surface_pressure_tendency!", "dynamics_tendencies!", "implicit_correction!")
        @test !occursin("pres_tend.data[1:1]", copy_text(ours, name))
    end
end

@testset "Moist-mass coupling model selection and l = 0 pressure tendency" begin
    base = R.ExperimentConfig(; truncation = 7, nlayers = 8, duration_days = 1,
        atmosphere_large_scale_precipitation = :net_column_conservative,
        forcing = R.ForcingConfig(radiation = :rrtmgp_clear_sky))
    results = Dict{Bool, Any}()
    for on in (false, true)
        simulation = R.build_rrtmgp(with(base; atmosphere_moist_mass_coupling = on))
        model, vars = simulation.model, simulation.variables
        @test (model.stochastic_physics isa R.MoistMassCoupling) == on
        @test (model isa R.MoistMassCoupledModel) == on
        selected(f, implicit_type) = which(f, Tuple{typeof(vars), implicit_type, typeof(model)}).module
        @test selected(S.dynamics_tendencies!, Int) == (on ? R : S)
        @test selected(S.implicit_correction!, typeof(model.implicit)) == (on ? R : S)
        S.run!(simulation; steps = 3)
        pres_tend = vars.tendencies.pressure.data
        @test all(isfinite, Array(vars.grid.temperature.data)) && all(isfinite, Array(vars.grid.pressure.data))
        results[on] = (; l0 = pres_tend[1],
            column = on ? sum(Float64.(Array(vars.parameterizations.mmc_column_source.data))) : 0.0)
    end
    @test results[false].l0 == 0                         # both upstream resets act on the uncoupled model
    @test results[true].l0 != 0 && results[true].column != 0
    # Explicit time stepping (α = 0) leaves implicit_correction! early in both methods.
    sim = R.build_rrtmgp(with(base; atmosphere_moist_mass_coupling = true))
    sim.model.implicit.α = 0
    S.run!(sim; steps = 2)
    @test all(isfinite, Array(sim.variables.grid.pressure.data))
end

@testset "Moist-mass coupling restart record" begin
    coupled = (; model = (; stochastic_physics = R.MoistMassCoupling()))
    plain = (; model = (; stochastic_physics = nothing))
    @test R._moist_mass_restart_state(coupled) == (; moist_mass_coupling = true)
    @test R._moist_mass_restart_state(plain) == (;)
    @test isnothing(R._validate_moist_mass_restart(coupled, (; moist_mass_coupling = true)))
    @test isnothing(R._validate_moist_mass_restart(plain, (; radiation = 1)))           # legacy: uncoupled
    @test_throws ArgumentError R._validate_moist_mass_restart(coupled, (; radiation = 1))
    @test_throws ArgumentError R._validate_moist_mass_restart(plain, (; moist_mass_coupling = true))
    @test_throws ArgumentError R._validate_moist_mass_restart(coupled, (; moist_mass_coupling = 1))
end
println("MOIST_MASS_COUPLING_UNIT_PASS")
end
