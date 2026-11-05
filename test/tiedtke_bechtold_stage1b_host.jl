module TiedtkeBechtoldStage1BHostTests
# Tiedtke-Bechtold stage 1B, PART B: host realisation (src/tiedtke_bechtold_stage1b_host.jl and its config, factory,
# metadata and restart hooks) against frozen_v6 (D48, D79, D81, D82, D83, G1-G6) and DESIGN_PARTB_v0.md P1-P9.
# Tests: the configuration name, the shared stage 1a refusals and the factory; declarations; T26 host cloud top; the
# D82 restart record and every cross-stage refusal (standalone); the P2 surface-flux refusal; and a short coupled
# analytic CPU run with stage 1B selected: the D81 n = 0 rule (step 1 equals stage 1a bit for bit), the D48 adapter
# read (a two-step sentinel on the four prescribed inputs), the column water ledger and diagnostics, restart pickup in
# a fresh build (bit for bit to step 4) and the real-checkpoint refusals. A fresh-PROCESS pickup and GPU are not
# covered here (later gates).
using ReadyESM, Test, Dates
const R = ReadyESM
const S = R.SpeedyWeather
const O = R.Oceananigans

fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
with(config; kwargs...) = R.ExperimentConfig(; merge(fields(config), (; kwargs...))...)
with_forcing(config; kwargs...) = with(config; forcing = R.ForcingConfig(; merge(fields(config.forcing), (; kwargs...))...))
include(joinpath(@__DIR__, "analytic_coupled_fixture.jl"))

const TB1A = :tiedtke_bechtold_eqcape_stage1
const TB1B = :tiedtke_bechtold_eqcape_stage1b
const MANIFEST_1B = "8d0f295c947a3466877bb25f26a1710d857c210f9c8e944804c48630d676108c"
const PRECISIONS = (Float64, Float32)
closure_tolerance(::Type{NF}, n) where {NF} = 64 * n * eps(NF)

bits(x::Float32) = reinterpret(UInt32, x)
bits(x::Float64) = reinterpret(UInt64, x)
bits(x::Complex) = (bits(real(x)), bits(imag(x)))
same_bits(a, b) = size(a) == size(b) && eltype(a) == eltype(b) && all(bits(x) == bits(y) for (x, y) in zip(a, b))
prognostic(vars) = map(n -> copy(Array(S.get_step(getproperty(vars.prognostic, n), 1).data)),
                       (:temperature, :humidity, :vorticity, :divergence, :pressure))
flag(atmosphere) = Array(atmosphere.variables.prognostic.convection.tiedtke_bechtold_lagged_flux_set.data)
const INPUTS = ((:ocean, :sensible_heat_flux), (:land, :sensible_heat_flux),
                (:ocean, :surface_humidity_flux), (:land, :surface_humidity_flux))
input_field(atmosphere, (space, name)) = getproperty(getproperty(atmosphere.variables.prognostic, space), name)

@testset "Tiedtke-Bechtold stage 1B host integration" begin

production = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)
tb_production = with(production; atmosphere_convection = TB1B, atmosphere_saturation_law = :liquid,
                     atmosphere_convection_timescale_seconds = 14400)

@testset "Configuration: name, shared stage 1a refusals, factory, provenance" begin
    @test R.validate(tb_production; check_input_files = false) === tb_production
    @test R._is_tiedtke_bechtold_convection(TB1A) && R._is_tiedtke_bechtold_convection(TB1B)
    @test !R._is_tiedtke_bechtold_convection(:betts_miller_sigma_rh_v1)
    @test R.TIEDTKE_BECHTOLD_STAGE1B_CONVECTION_SYMBOL == TB1B
    # The stage 1a refusals apply unchanged, with the stage 1B name in the message.
    @test_throws r"atmosphere_convection=tiedtke_bechtold_eqcape_stage1b requires atmosphere_saturation_law=liquid" R.validate(
        with(tb_production; atmosphere_saturation_law = :mixed_phase); check_input_files = false)
    @test_throws ArgumentError R.validate(with(tb_production; atmosphere_moist_mass_coupling = true);
                                          check_input_files = false)
    @test_throws r"prognostic cloud condensate" R.validate(
        with_forcing(tb_production; cloud_scheme = :prognostic_condensate, cloud_condensate_retention_fraction = 0.5);
        check_input_files = false)
    for name in R._TIEDTKE_BECHTOLD_INERT_CONVECTION_SETTINGS
        value = getfield(R.ExperimentConfig(), name)
        changed = value isa Integer ? value + 1 : value * 0.5 + 0.1
        @test_throws ArgumentError R.validate(with(tb_production; NamedTuple{(name,)}((changed,))...);
                                              check_input_files = false)
    end
    # Factory: the stage 1B wrapper with the stage 1a scheme and the Sect. 9 shallow parameters.
    grid = S.SpectralGrid(; NF = Float32, trunc = 7, nlayers = 8)
    convection = R._atmosphere_convection(tb_production, grid)
    @test convection isa R.TiedtkeBechtoldStage1BHostConvection
    @test convection.convection.parameters === R.tiedtke_bechtold_parameters(Float32; truncation = 7)
    @test convection.shallow === R.tiedtke_bechtold_shallow_parameters(Float32)
    @test R._atmosphere_convection(with(tb_production; atmosphere_convection = TB1A), grid) isa
          R.TiedtkeBechtoldHostConvection                                          # stage 1a unchanged
    # Restart identity binds the symbol; provenance carries the stage 1B values.
    @test R._dynamic_config_identity_hash(tb_production) !=
          R._dynamic_config_identity_hash(with(tb_production; atmosphere_convection = TB1A))
    meta = R._tiedtke_bechtold_metadata(tb_production)
    @test meta.atmosphere_convection_stage == "stage1b"
    @test meta.atmosphere_convection_trigger_variant == "shallow_n_deep_n1_v1"
    @test meta.atmosphere_convection_parcel_adapter == "coupled_prescribed_inputs_v1"
    @test meta.atmosphere_convection_contract_manifest_sha256 == MANIFEST_1B
    @test R._tiedtke_bechtold_metadata(with(tb_production; atmosphere_convection = TB1A)).atmosphere_convection_stage ==
          "stage1a"
    @test R._tiedtke_bechtold_metadata(production) == (;)
end

@testset "Declarations: stage 1a fields plus the D81 flag, D83 diagnostics and counters, isbits $NF" for NF in PRECISIONS
    grid = S.SpectralGrid(; NF, trunc = 7, nlayers = 8)
    convection = R.TiedtkeBechtoldStage1BHostConvection(grid)
    @test isbits(convection)
    @test R.Adapt.adapt(Array, convection) === convection
    declared = S.variables(convection)
    names = [v.name for v in declared]
    @test allunique(names)
    stage1a = [v.name for v in S.variables(R.TiedtkeBechtoldHostConvection(convection.convection))]
    @test issubset(stage1a, names)
    set_flag = only(filter(v -> v.name == :tiedtke_bechtold_lagged_flux_set, declared))
    @test set_flag isa S.PrognosticVariable && set_flag.dims isa S.Grid2D && set_flag.namespace == :convection
    for (name, _, _) in (R._TIEDTKE_BECHTOLD_STAGE1B_DIAGNOSTICS..., R._TIEDTKE_BECHTOLD_STAGE1B_COUNTERS...)
        v = only(filter(v -> v.name == name, declared))
        @test v isa S.ParameterizationVariable && v.dims isa S.Grid2D
    end
    @test length(names) == length(stage1a) + 1 + length(R._TIEDTKE_BECHTOLD_STAGE1B_DIAGNOSTICS) +
                           length(R._TIEDTKE_BECHTOLD_STAGE1B_COUNTERS)
    @test hasmethod(S.parameterization!, Tuple{Int, Any, typeof(convection), Any})
    @test R._uses_tiedtke_bechtold_convection(convection)
end

@testset "T26 host cloud top: shallow leaves it unchanged, deep as D43 (D79 (a)) $NF" for NF in PRECISIONS
    n = 27
    top = R.tiedtke_bechtold_stage1b_host_cloud_top
    @test top(NF(n + 1), n + 1, 0, n) === NF(n + 1)
    @test top(NF(2), n + 1, 0, n) === NF(2)
    for preset in (n + 1, 5, 2)
        @test top(NF(preset), 3, 2, n) === NF(preset)                  # shallow: unchanged, whatever its top
        @test top(NF(preset), 3, 1, n) === NF(min(preset, 3))           # deep: the D43 min rule
    end
    for (tb_top, kind) in ((0, 1), (n + 2, 1), (n + 1, 1), (0, 2), (n + 1, 2), (n + 2, 2), (3, 0), (3, 3), (3, -1))
        @test_throws ArgumentError top(NF(n + 1), tb_top, kind, n)
    end
    @test_throws ArgumentError R.tiedtke_bechtold_host_cloud_top(NF(n + 1), 3, 2, n)   # stage 1a still refuses 2
end

@testset "D82 restart record and cross-stage refusals (standalone)" begin
    grid = S.SpectralGrid(; NF = Float32, trunc = 7, nlayers = 8)
    sim(convection) = (; model = (; convection))
    tb1b = R.TiedtkeBechtoldStage1BHostConvection(grid)
    tb1a = R.TiedtkeBechtoldHostConvection(grid)
    bm = R.NetColumnBettsMillerConvection(S.BettsMillerConvection(grid))
    contract = R.tiedtke_bechtold_restart_contract(tb1b)
    @test contract == (; schema = "readiesm_convection_restart_contract_v1", scheme = "tiedtke_bechtold_eqcape_stage1",
        stage_version = "stage1b", water_contract = "instantaneous_rain_no_rup_v1",
        trigger_variant = "shallow_n_deep_n1_v1", parcel_adapter = "coupled_prescribed_inputs_v1",
        contract_manifest_sha256 = MANIFEST_1B)
    @test R.tiedtke_bechtold_restart_contract(tb1a) == R.tiedtke_bechtold_restart_contract()   # stage 1a unchanged
    record_1b = (; convection_restart_contract = contract)
    record_1a = (; convection_restart_contract = R.tiedtke_bechtold_restart_contract())
    @test R._convection_restart_state(sim(tb1b)) == record_1b
    @test R._convection_restart_state(sim(tb1a)) == record_1a
    @test R._convection_restart_state(sim(bm)) == (;)
    @test R._validate_convection_restart(sim(tb1b), record_1b) === nothing
    @test R._validate_convection_restart(sim(tb1b), (; convection_restart_contract = (; reverse(pairs(contract))...))) ===
          nothing                                                          # JLD2 key order does not matter
    @test_throws r"stage_version" R._validate_convection_restart(sim(tb1b), record_1a)     # 1a -> 1B
    @test_throws r"stage_version" R._validate_convection_restart(sim(tb1a), record_1b)     # 1B -> 1a
    @test_throws ArgumentError R._validate_convection_restart(sim(tb1b), (;))              # BM -> 1B
    @test_throws ArgumentError R._validate_convection_restart(sim(bm), record_1b)          # 1B -> BM
    for key in keys(contract)
        @test_throws ArgumentError R._validate_convection_restart(sim(tb1b),
            (; convection_restart_contract = merge(contract, NamedTuple{(key,)}(("changed",)))))
    end
    @test_throws ArgumentError R.with_dynamic_restart_identity_override(
        () -> R._validate_convection_restart(sim(tb1b), record_1a))
end

@testset "P2: stage 1B refuses non-prescribed surface fluxes" begin
    grid = S.SpectralGrid(; NF = Float32, trunc = 7, nlayers = 8)
    prescribed = (; surface_heat_flux = S.SurfaceHeatFlux(grid; ocean = S.PrescribedOceanHeatFlux(grid),
                                                           land = S.PrescribedLandHeatFlux()),
                    surface_humidity_flux = S.SurfaceHumidityFlux(grid; ocean = S.PrescribedOceanHumidityFlux(grid),
                                                                   land = S.PrescribedLandHumidityFlux()))
    @test R._check_tiedtke_bechtold_stage1b_surface_fluxes(prescribed) === nothing
    computed_land = merge(prescribed, (; surface_heat_flux = S.SurfaceHeatFlux(grid;
        ocean = S.PrescribedOceanHeatFlux(grid), land = S.SurfaceLandHeatFlux(grid))))
    @test_throws r"land component must be PrescribedLandHeatFlux" R._check_tiedtke_bechtold_stage1b_surface_fluxes(
        computed_land)
    computed_ocean = merge(prescribed, (; surface_humidity_flux = S.SurfaceHumidityFlux(grid;
        ocean = S.SurfaceOceanHumidityFlux(grid), land = S.PrescribedLandHumidityFlux())))
    @test_throws r"ocean component must be PrescribedOceanHumidityFlux" R._check_tiedtke_bechtold_stage1b_surface_fluxes(
        computed_ocean)
end

@testset "Coupled analytic CPU run with stage 1B: n = 0 rule, D48 read, ledger, pickup, refusals" begin
    mktempdir() do dir
        output(name) = mkpath(joinpath(dir, name))
        base = with(analytic_coupled_test_config(production; output_dir = output("base"));
            atmosphere_timestep_at_t31_minutes = 3.75, duration_days = 4 * 900 / 86_400,
            start_date = DateTime(1993, 1, 1))
        tb_config = with(base; atmosphere_convection = TB1B, atmosphere_saturation_law = :liquid,
                         atmosphere_convection_timescale_seconds = 14400, output_dir = output("reference"))
        @test R.validate(tb_config; check_input_files = false) === tb_config
        reference = R.build_dynamic_esm(tb_config)
        atmosphere = reference.model.atmosphere
        convection = atmosphere.model.convection
        @test convection isa R.TiedtkeBechtoldStage1BHostConvection
        @test convection.convection.parameters === R.tiedtke_bechtold_parameters(atmosphere.model)
        @test all(iszero, flag(atmosphere))                                     # D81: unset at a fresh start
        model = atmosphere.model
        Lv = convection.convection.parameters.latent_heat_condensation
        land_fraction = Array(model.land_sea_mask.mask.data)
        # Sentinel (T23 (ii)): after step 3 every input is replaced by a distinct uniform value, so step 4 must read
        # exactly these values through the D48 adapter.
        sentinel = (21.25f0, 7.5f0, 6.25f-5, 1.875f-5)
        apply_sentinel!(atm) = foreach(((name, value),) -> fill!(input_field(atm, name).data, value),
                                       zip(INPUTS, sentinel))
        expected_lagged(ij) = (R._pbl_area_weighted(input_field(atmosphere, INPUTS[2]), input_field(atmosphere, INPUTS[1]),
                                                    ij, Float32(land_fraction[ij])),
                               R._pbl_area_weighted(input_field(atmosphere, INPUTS[4]), input_field(atmosphere, INPUTS[3]),
                                                    ij, Float32(land_fraction[ij])))
        path = joinpath(dir, "tb1b_step2.jld2")
        activity = NamedTuple[]
        function observe(sim)
            k = sim.model.clock.iteration
            p = sim.model.atmosphere.variables.parameterizations
            type = Array(p.tiedtke_bechtold_convection_type.data)
            push!(activity, (; iteration = k, flag = copy(flag(sim.model.atmosphere)),
                trigger_path = copy(Array(p.tiedtke_bechtold_trigger_path.data)),
                unset_steps = copy(Array(p.tiedtke_bechtold_lagged_flux_unset_steps.data)),
                lagged_H = copy(Array(p.tiedtke_bechtold_lagged_sensible_heat_flux.data)),
                lagged_E = copy(Array(p.tiedtke_bechtold_lagged_humidity_flux.data)),
                source = copy(Array(p.tiedtke_bechtold_shallow_source.data)), type = copy(type),
                shallow_columns = count(==(2), type), deep_columns = count(==(1), type),
                shallow_tests = count(>(0), Array(p.tiedtke_bechtold_shallow_test_base_half_level.data))))
            if k == 2
                R.save_dynamic_restart_state(sim, path)
                push!(activity, (; iteration = -2,
                    inputs = map(name -> copy(Array(input_field(sim.model.atmosphere, name).data)), INPUTS)))
            end
            if k == 3
                apply_sentinel!(sim.model.atmosphere)
                push!(activity, (; iteration = -3, expected = [expected_lagged(ij) for ij in eachindex(land_fraction)]))
            end
            return nothing
        end
        O.add_callback!(reference, observe, O.IterationInterval(1); name = :tiedtke_bechtold_stage1b_observer)
        O.run!(reference)
        @test reference.model.clock.iteration == 4 && isfile(path)
        final = prognostic(atmosphere.variables)
        @test all(x -> all(isfinite, x), final)
        stepped = Dict(a.iteration => a for a in activity if a.iteration >= 1)
        @test Set(keys(stepped)) == Set(1:4)
        # D81: step 1 ran the stage 1a fallback everywhere and set every flag; later steps run stage 1B.
        @test all(==(1), stepped[1].trigger_path) && all(==(1), stepped[1].unset_steps)
        @test all(==(1), stepped[1].flag)
        @test all(iszero, stepped[1].lagged_H) && all(iszero, stepped[1].lagged_E)
        for k in 2:4
            @test all(==(2), stepped[k].trigger_path) && all(==(1), stepped[k].unset_steps)
        end
        # T23 (ii): step 4 read exactly the sentinel inputs written after step 3 (area weighting as the K-profile).
        expected = only(a for a in activity if a.iteration == -3).expected
        @test all(bits(stepped[4].lagged_H[ij]) == bits(expected[ij][1]) for ij in eachindex(expected))
        @test all(bits(stepped[4].lagged_E[ij]) == bits(expected[ij][2]) for ij in eachindex(expected))
        shallow = findall(==(2), stepped[4].type)
        @test all(stepped[4].source[ij] == expected[ij][1] + Lv * expected[ij][2] for ij in shallow)   # D71
        # Diagnostics and the column water ledger of the last step (as the stage 1a host test).
        diagnostics = R.collect_dynamic_diagnostics(reference, tb_config)
        R.save_dynamic_diagnostics(diagnostics, tb_config; render_figure = false)
        R.validate_dynamic_diagnostics(diagnostics)
        @test diagnostics.metadata.atmosphere_convection == "tiedtke_bechtold_eqcape_stage1b"
        @test diagnostics.metadata.atmosphere_convection_stage == "stage1b"
        p = atmosphere.variables.parameterizations
        @test all(iszero, Array(p.convective_precipitation_mass_correction.data))
        dq = Float32.(Array(p.convective_humidity_tendency.data))
        ps = Array(p.convective_precipitation_surface_pressure.data)
        dsigma = Array(model.geometry.σ_levels_thick)
        gravity = model.planet.gravity
        flux = Array(p.tiedtke_bechtold_rain_mass_flux_surface.data)
        n = length(dsigma)
        worst = maximum(abs(flux[ij] + ps[ij] / gravity * sum(dq[ij, :] .* dsigma)) /
                        max(flux[ij] + ps[ij] / gravity * sum(abs.(dq[ij, :]) .* dsigma), floatmin(Float32))
                        for ij in eachindex(flux))
        @test worst <= closure_tolerance(Float32, n)
        println("TIEDTKE_BECHTOLD_STAGE1B_COUPLED_ACTIVITY ",
                join(("step$(k): shallow=$(stepped[k].shallow_columns) deep=$(stepped[k].deep_columns) " *
                      "shallow_tests=$(stepped[k].shallow_tests)" for k in 1:4), " | "),
                " columns=$(length(flux)) capped_total=$(sum(Array(p.tiedtke_bechtold_shallow_cfl_capped_steps.data)))")
        # T23 (i): step 1 equals a stage 1a run bit for bit (the n = 0 fallback in the host).
        stage1a = R.build_dynamic_esm(with(tb_config; atmosphere_convection = TB1A, output_dir = output("stage1a")))
        stage1a.stop_iteration = 1
        O.run!(stage1a)
        fresh_1b = R.build_dynamic_esm(with(tb_config; output_dir = output("step1")))
        fresh_1b.stop_iteration = 1
        O.run!(fresh_1b)
        @test all(same_bits(x, y) for (x, y) in zip(prognostic(fresh_1b.model.atmosphere.variables),
                                                       prognostic(stage1a.model.atmosphere.variables)))
        # T23 (iii) in a fresh build: restore the step-2 checkpoint, apply the same sentinel after step 3, and
        # continue to step 4 bit for bit. The flag and the counters come back from the checkpoint.
        continued = R.build_dynamic_esm(with(tb_config; output_dir = output("continued")))
        R.restore_dynamic_restart_state!(continued, path)
        @test continued.model.clock.iteration == 2
        # Paul's T23 (iii) requirement: right after restore, before any step, the four persistent inputs are the saved
        # step-2 values bit for bit (so no restore-time update changes them and a later step compensates), and so is
        # the flag.
        saved_inputs = only(a for a in activity if a.iteration == -2).inputs
        for (name, saved) in zip(INPUTS, saved_inputs)
            @test same_bits(Array(input_field(continued.model.atmosphere, name).data), saved)
        end
        @test same_bits(flag(continued.model.atmosphere), stepped[2].flag)
        @test all(==(1), flag(continued.model.atmosphere))
        @test all(==(1), Array(continued.model.atmosphere.variables.parameterizations.tiedtke_bechtold_lagged_flux_unset_steps.data))
        continued.stop_iteration = 4
        continued.running = true
        while continued.running
            O.time_step!(continued)
            continued.model.clock.iteration == 3 && apply_sentinel!(continued.model.atmosphere)
        end
        @test continued.model.clock.iteration == 4
        @test all(same_bits(x, y) for (x, y) in zip(prognostic(continued.model.atmosphere.variables), final))
        # T23 (iv) with real checkpoints: stage 1B into stage 1a is refused, also with the identity check waived.
        @test_throws ErrorException R.restore_dynamic_restart_state!(stage1a, path)
        @test_throws r"stage_version" R.restore_dynamic_restart_state!(stage1a, path; allow_identity_mismatch = true)
        @test stage1a.model.clock.iteration == 1
    end
end

end
println("TIEDTKE_BECHTOLD_STAGE1B_HOST_TESTS_PASS")
end
