module SeaIceLeadClosingTests
# Opt-in sea-ice lead-closing thickness (src/sea_ice_lead_closing.jl, vendored ClimaSeaIce LeadClosingEvolution; design
# artifacts/development_wave_20261103/sea_ice_lead_closing_claude/DESIGN.md): settings, the kernel rule against the
# released ProportionalEvolution, the freeze/melt ratchet replayed through the vendored volume update, restart identity,
# and a short coupled analytic run whose live sea-ice model carries the rule and whose diagnostics validate.
using ReadyESM, Test, Dates
const R = ReadyESM
const CSI = R.ClimaSeaIce.SeaIceThermodynamics
const O = R.Oceananigans

fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
with(config; kwargs...) = R.ExperimentConfig(; merge(fields(config), (; kwargs...))...)
include(joinpath(@__DIR__, "analytic_coupled_fixture.jl"))

# 31 days at 15-minute steps of a cell-net volume tendency (m per s per ice area) with a day/night cycle, through the
# vendored ice_volume_update. Illustrative forcing: daily mean -1 cm/day, night +0.5 cm/day for 8 h.
function replay(rule, ::Type{FT}; A = FT(0.5), h = FT(0.6), hc = FT(0.05), dt = FT(900)) where {FT}
    therm = (; concentration_evolution = rule)
    night, mean, hours = FT(0.5), FT(-1.0), 8
    day = (mean * 24 - night * hours) / (24 - hours)
    for d in 1:31, k in 0:95
        rate = (k * 900 / 3600 < hours ? night : day) / FT(100) / FT(86400)
        h, A = CSI.ice_volume_update(therm, rate * A, h, A, hc, dt)
        A > 0 || return (A, h)
    end
    return (A, h)
end

@testset "Sea-ice lead closing: settings and kernel rule" begin
    default = R.ExperimentConfig()
    @test default.sea_ice_lead_closing_thickness_m == 0.0
    production = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)
    @test production.sea_ice_lead_closing_thickness_m == 0.0
    on = with(production; sea_ice_lead_closing_thickness_m = 0.5)
    @test R.validate(on; check_input_files = false) === on
    for bad in (-0.5, 0.01, 2.5, NaN, Inf)
        @test_throws ArgumentError R.validate(with(production; sea_ice_lead_closing_thickness_m = bad);
            check_input_files = false)
    end
    @test R._dynamic_config_identity_hash(on) != R._dynamic_config_identity_hash(production)
    @test isnothing(R._sea_ice_concentration_evolution(production, Float32))
    @test R._sea_ice_concentration_evolution(on, Float32) === CSI.LeadClosingEvolution(0.5f0)
    @test R._sea_ice_lead_closing_provenance(production) == (;)
    @test R._sea_ice_lead_closing_provenance(on).sea_ice_lead_closing_thickness_m == 0.5

    prop = CSI.ProportionalEvolution()
    for FT in (Float32, Float64)
        hc = FT(0.05)
        @test CSI.freezing_concentration_thickness(prop, hc) === hc
        @test CSI.freezing_concentration_thickness(CSI.LeadClosingEvolution(FT(0.5)), hc) === FT(0.5)
        same = CSI.LeadClosingEvolution(hc)
        for A in FT.((0.0, 0.2, 0.55, 0.9, 1.0)), h in FT.((0.1, 0.6, 2.0)), dV in FT.((-3e-7, -1e-8, 0, 1e-8, 3e-7))
            dt = FT(900)
            # At h_0 = h_c the new rule is the released rule, bit for bit, in both branches.
            @test CSI.concentration_thermodynamic_step(same, dV, A, h, hc, dt) ===
                  CSI.concentration_thermodynamic_step(prop, dV, A, h, hc, dt)
            lead = CSI.concentration_thermodynamic_step(CSI.LeadClosingEvolution(FT(0.5)), dV, A, h, hc, dt)
            if dV < 0      # melting is unchanged
                @test lead === CSI.concentration_thermodynamic_step(prop, dV, A, h, hc, dt)
            elseif dV > 0  # freezing closes leads at (1 - A)/h_0
                @test lead ≈ max(zero(FT), A + dt * (1 - A) / FT(0.5) * dV)
            end
        end
    end
    # The ratchet: under net melt the released rule raises A, while h_0 = 0.5 m lets it fall.
    for FT in (Float32, Float64)
        A_release, h_release = replay(prop, FT)
        A_lead, _ = replay(CSI.LeadClosingEvolution(FT(0.5)), FT)
        println("LEAD_CLOSING_RATCHET FT=$FT release_A=$A_release release_h=$h_release lead_0.5m_A=$A_lead")
        @test A_release > FT(0.5) > A_lead
    end
end

@testset "Sea-ice lead closing: coupled analytic run carries the rule" begin
    production = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)
    base = with(analytic_coupled_test_config(production; output_dir = mktempdir());
        atmosphere_timestep_at_t31_minutes = 3.75, duration_days = 4 * 900 / 86_400, start_date = DateTime(1993, 1, 1))
    off = R.build_dynamic_esm(base)
    @test off.model.sea_ice.model.ice_thermodynamics.concentration_evolution isa CSI.ProportionalEvolution
    on = with(base; sea_ice_lead_closing_thickness_m = 0.5, output_dir = mktempdir())
    simulation = R.build_dynamic_esm(on)
    evolution = simulation.model.sea_ice.model.ice_thermodynamics.concentration_evolution
    @test evolution isa CSI.LeadClosingEvolution && evolution.lead_closing_thickness == 0.5
    O.run!(simulation)
    @test simulation.model.clock.iteration == 4
    A = Array(O.interior(simulation.model.sea_ice.model.ice_concentration))
    @test all(isfinite, A) && all(x -> 0 <= x <= 1, A)
    diagnostics = R.collect_dynamic_diagnostics(simulation, on)
    R.save_dynamic_diagnostics(diagnostics, on; render_figure = false)
    R.validate_dynamic_diagnostics(diagnostics)
    @test diagnostics.metadata.sea_ice_lead_closing_thickness_m == 0.5
    @test diagnostics.metadata.sea_ice_concentration_evolution == "proportional_with_lead_closing_thickness_v1"
end
println("SEA_ICE_LEAD_CLOSING_TESTS_PASS")
end
