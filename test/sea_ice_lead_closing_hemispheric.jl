module SeaIceLeadClosingHemisphericTests
# Opt-in southern lead-closing thickness (sea_ice_lead_closing_thickness_south_m, vendored ClimaSeaIce
# HemisphericLeadClosingEvolution; design artifacts/development_wave_20261104/lead_closing_hemispheric_claude_v1/DESIGN.md):
# settings and refusals, the restart identity (unchanged while unset), the per-cell rule on a two-hemisphere grid
# against the global and released rules bit for bit, and a short coupled analytic run whose live sea-ice model carries
# the rule and whose diagnostics validate.
using ReadyESM, Test, Dates
const R = ReadyESM
const CSI = R.ClimaSeaIce.SeaIceThermodynamics
const O = R.Oceananigans

fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
with(config; kwargs...) = R.ExperimentConfig(; merge(fields(config), (; kwargs...))...)
include(joinpath(@__DIR__, "analytic_coupled_fixture.jl"))

# Independently construct the v0.1.21 typed identity with the optional south
# key omitted. An unset key must not change the identity of an old checkpoint.
function hash_without_key(config)
    io = IOBuffer()
    print(io, "readiesm-experiment-config-v2")
    for name in fieldnames(typeof(config))
        name in R._DYNAMIC_RESTART_RUN_CONTROL_FIELDS && continue
        name == :sea_ice_lead_closing_thickness_south_m && continue
        R._dynamic_restart_write_field!(io, name, getfield(config, name), :experiment)
    end
    return bytes2hex(R.SHA.sha256(take!(io)))
end

@testset "Hemispheric lead closing: settings and identity" begin
    @test isnothing(R.ExperimentConfig().sea_ice_lead_closing_thickness_south_m)
    production = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)
    @test isnothing(production.sea_ice_lead_closing_thickness_south_m)
    global_on = with(production; sea_ice_lead_closing_thickness_m = 0.5)
    north = with(production; sea_ice_lead_closing_thickness_m = 0.5, sea_ice_lead_closing_thickness_south_m = 0.0)
    south = with(production; sea_ice_lead_closing_thickness_m = 0.0, sea_ice_lead_closing_thickness_south_m = 0.5)
    split = with(production; sea_ice_lead_closing_thickness_m = 0.5, sea_ice_lead_closing_thickness_south_m = 0.3)
    for config in (north, south, split)
        @test R.validate(config; check_input_files = false) === config
    end
    for bad in (-0.5, 0.01, 2.5, NaN, Inf)
        @test_throws ArgumentError R.validate(with(north; sea_ice_lead_closing_thickness_south_m = bad);
            check_input_files = false)
    end
    # Equal values are the global rule (or the release), so the key must stay unset.
    for h in (0.0, 0.5)
        @test_throws ArgumentError R.validate(with(production; sea_ice_lead_closing_thickness_m = h,
            sea_ice_lead_closing_thickness_south_m = h); check_input_files = false)
    end
    # Production now selects global lead closing at 0.5 m (LC1), so it equals global_on. The released rule is 0.0.
    released = with(production; sea_ice_lead_closing_thickness_m = 0.0)
    @test production.sea_ice_lead_closing_thickness_m == 0.5
    @test R._dynamic_config_identity_hash(production) == R._dynamic_config_identity_hash(global_on)
    # Unset, the identity is the pre-key identity, so existing checkpoints keep their identity. Set, it is bound.
    for config in (released, production, global_on)
        @test R._dynamic_config_identity_hash(config) == hash_without_key(config)
    end
    # The released v0.1.21 serializer and production config have these exact
    # portable identities; output paths and checkout directories are excluded.
    @test R._dynamic_config_identity_hash(released) ==
        "42f78154105e1810160d5a97886ae3e788adf488e9df40463ebf6b42759da899"
    @test R._dynamic_config_identity_hash(production) ==
        "fe4f0d8c32835a40530a1294a3465adc733a48127ab57369b9d5abb017a5d3b4"
    hashes = R._dynamic_config_identity_hash.((released, global_on, north, south, split))
    @test allunique(hashes)

    @test isnothing(R._sea_ice_concentration_evolution(released, Float32))
    @test R._sea_ice_concentration_evolution(production, Float32) === CSI.LeadClosingEvolution(0.5f0)
    @test R._sea_ice_concentration_evolution(global_on, Float32) === CSI.LeadClosingEvolution(0.5f0)
    @test R._sea_ice_concentration_evolution(north, Float32) === CSI.HemisphericLeadClosingEvolution(0.5f0, 0.0f0)
    @test R._sea_ice_concentration_evolution(south, Float64) === CSI.HemisphericLeadClosingEvolution(0.0, 0.5)
    @test R._sea_ice_lead_closing_provenance(global_on) == (;
        sea_ice_concentration_evolution = "proportional_with_lead_closing_thickness_v1",
        sea_ice_lead_closing_thickness_m = 0.5, sea_ice_consolidation_thickness_m = 0.05)
    provenance = R._sea_ice_lead_closing_provenance(north)
    @test provenance.sea_ice_concentration_evolution == "proportional_with_hemispheric_lead_closing_thickness_v1"
    @test provenance.sea_ice_lead_closing_thickness_m == 0.5 && provenance.sea_ice_lead_closing_thickness_south_m == 0.0
end

@testset "Hemispheric lead closing: per-cell rule" begin
    prop = CSI.ProportionalEvolution()
    for FT in (Float32, Float64)
        # Rows 1-2 are south of the equator and rows 3-4 north of it (centres 67.5S, 22.5S, 22.5N, 67.5N). The kernel
        # tests construct HemisphericLeadClosingEvolution directly, so (0.5, 0.5) is allowed here although the
        # configuration refuses equal values.
        grid = O.Grids.LatitudeLongitudeGrid(O.CPU(), FT; size = (4, 4, 1), halo = (1, 1, 1), longitude = (0, 360),
                                             latitude = (-90, 90), z = (-1, 0))
        C = O.Grids.Center()
        @test [sign(O.Grids.φnode(1, j, 1, grid, C, C, C)) for j in 1:4] == [-1, -1, 1, 1]
        S, N = 1, 4
        hc, dt = FT(0.05), FT(900)
        lead = CSI.LeadClosingEvolution(FT(0.5))
        # The released and global rules are the same object in every cell.
        for i in 1:4, j in 1:4
            @test CSI.local_concentration_evolution(prop, i, j, grid, hc) === prop
            @test CSI.local_concentration_evolution(lead, i, j, grid, hc) === lead
        end
        north_only = CSI.HemisphericLeadClosingEvolution(FT(0.5), zero(FT))
        south_only = CSI.HemisphericLeadClosingEvolution(zero(FT), FT(0.5))
        both = CSI.HemisphericLeadClosingEvolution(FT(0.5), FT(0.5))
        thickness(rule, i, j) = CSI.freezing_concentration_thickness(CSI.local_concentration_evolution(rule, i, j, grid, hc), hc)
        for i in 1:4, j in 1:4
            north_row = j > 2
            @test thickness(north_only, i, j) === (north_row ? FT(0.5) : hc)
            @test thickness(south_only, i, j) === (north_row ? hc : FT(0.5))
        end
        for A in FT.((0.0, 0.2, 0.55, 0.9, 1.0)), h in FT.((0.1, 0.6, 2.0)), dV in FT.((-3e-7, -1e-8, 0, 1e-8, 3e-7))
            conc(rule) = CSI.concentration_thermodynamic_step(rule, dV, A, h, hc, dt)
            volume(rule) = CSI.ice_volume_update((; concentration_evolution = rule), dV, h, A, hc, dt)
            local_rule(rule, j) = CSI.local_concentration_evolution(rule, 1, j, grid, hc)
            # North-only: the global rule in the north and the released rule in the south, bit for bit.
            @test conc(local_rule(north_only, N)) === conc(lead)
            @test conc(local_rule(north_only, S)) === conc(prop)
            @test volume(local_rule(north_only, N)) === volume(lead)
            @test volume(local_rule(north_only, S)) === volume(prop)
            # South-only is the mirror. Equal values are the global rule everywhere.
            @test conc(local_rule(south_only, S)) === conc(lead)
            @test conc(local_rule(south_only, N)) === conc(prop)
            @test conc(local_rule(both, S)) === conc(local_rule(both, N)) === conc(lead)
        end
    end
end

@testset "Hemispheric lead closing: hemisphere on the production tripolar grid" begin
    # The horizontal production grid (one layer), as built in _build_ocean. Each column runs south to north across one
    # sign change, the fold rows are all north, and the rule picks the northern value exactly where latitude >= 0.
    production = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)
    @test production.ocean_grid == :tripolar_1degree
    Nx, Ny = production.ocean_nlongitude, production.ocean_nlatitude
    grid = O.TripolarGrid(O.CPU(), Float32; size = (Nx, Ny, 1), z = (-1.0f0, 0.0f0), halo = (5, 5, 1))
    C = O.Grids.Center()
    φ = [O.Grids.φnode(i, j, 1, grid, C, C, C) for i in 1:Nx, j in 1:Ny]
    @test all(isfinite, φ)
    @test all(i -> issorted(sign.(φ[i, :])), 1:Nx)
    @test all(>(0), φ[:, (Ny - 5):Ny])
    @test any(<(0), φ) && any(>(0), φ)
    rule = CSI.HemisphericLeadClosingEvolution(0.5f0, 0.0f0)
    hc = 0.05f0
    picked = [CSI.freezing_concentration_thickness(CSI.local_concentration_evolution(rule, i, j, grid, hc), hc)
              for i in 1:Nx, j in 1:Ny]
    @test picked == ifelse.(φ .>= 0, 0.5f0, hc)
    println("HEMISPHERIC_TRIPOLAR_SIGN nx=$Nx ny=$Ny north_cells=$(count(>=(0), φ)) south_cells=$(count(<(0), φ)) " *
            "min_north_lat=$(minimum(φ[φ .>= 0])) max_south_lat=$(maximum(φ[φ .< 0]))")
end

@testset "Hemispheric lead closing: coupled analytic run carries the rule" begin
    production = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)
    base = with(analytic_coupled_test_config(production; output_dir = mktempdir());
        atmosphere_timestep_at_t31_minutes = 3.75, duration_days = 4 * 900 / 86_400, start_date = DateTime(1993, 1, 1))
    north = with(base; sea_ice_lead_closing_thickness_m = 0.5, sea_ice_lead_closing_thickness_south_m = 0.0,
                 output_dir = mktempdir())
    simulation = R.build_dynamic_esm(north)
    evolution = simulation.model.sea_ice.model.ice_thermodynamics.concentration_evolution
    @test evolution isa CSI.HemisphericLeadClosingEvolution && evolution.north == 0.5 && evolution.south == 0
    O.run!(simulation)
    @test simulation.model.clock.iteration == 4
    A = Array(O.interior(simulation.model.sea_ice.model.ice_concentration))
    @test all(isfinite, A) && all(x -> 0 <= x <= 1, A)
    diagnostics = R.collect_dynamic_diagnostics(simulation, north)
    R.save_dynamic_diagnostics(diagnostics, north; render_figure = false)
    R.validate_dynamic_diagnostics(diagnostics)
    @test diagnostics.metadata.sea_ice_concentration_evolution ==
          "proportional_with_hemispheric_lead_closing_thickness_v1"
    @test diagnostics.metadata.sea_ice_lead_closing_thickness_m == 0.5
    @test diagnostics.metadata.sea_ice_lead_closing_thickness_south_m == 0.0
end
println("SEA_ICE_LEAD_CLOSING_HEMISPHERIC_TESTS_PASS")
end
