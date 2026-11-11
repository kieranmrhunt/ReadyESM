module SeaIceLeadClosingSnowFixedPointTests
# Snow-covered freezing fixed point of the opt-in sea-ice lead-closing thickness h_0 (src/sea_ice_lead_closing.jl,
# vendored ClimaSeaIce LeadClosingEvolution, design artifacts/development_wave_20261103/sea_ice_lead_closing_claude/
# DESIGN.md). test/sea_ice_lead_closing.jl covers the scalar rule, the no-snow ratchet and a short coupled run. This
# file runs the real layered snow and ice kernel (_layered_thermodynamic_time_step! in
# vendor/ClimaSeaIce/src/SeaIceThermodynamics/thermodynamic_time_step.jl) for one call on four CPU cells with
# prescribed inputs. A Float64 scalar oracle written here, which calls no ClimaSeaIce helper, checks:
#   1. a snow-covered freezing fixture with snow melt (Q_s > 0) at h_0 = 0.5 and 0.2 m, and at the release h_c;
#   2. the closed form A_next = (A + K α)/(1 - K β) with K = Δt (1 - A)/h_0, and its dV/dt, against the returned
#      concentration, ice and snow thickness, snow-ice interface temperature, and ice, snow and snowfall mass fluxes;
#   3. the residual dV/dt - (α + β A_next) of the returned state at rounding, and the column energy closure;
#   4. the snow rebase hs A / A_next with melt and snowfall, flooding inactive, and the snow volume accounting;
#   5. h_0 = h_c (0.05 m) against the released ProportionalEvolution, bit for bit in every output.
# The oracle needs no value read back from the kernel. The snow surface temperature is prescribed (as in ReadyESM), so
# Q_s follows from the inputs. The kernel's stored snow-ice interface temperature is read back only to cross-check.
#
# Latent heat. The closed form uses ρ_i ℒ₀ (the reference latent heat). The final dV/dt comes from
# ice_melt_freeze_tendency (slab_thermodynamics_tendencies.jl):
#   dV/dt = (Q_u + Q_s A_tmp - Q_ii)/(ρ_i ℒ(T_si)) + (Q_ii - Q_b)/(ρ_i ℒ(T_b)),
#   ℒ(T) = ℒ₀ + (ρ_ℓ c_ℓ/ρ - c)(T - T₀),
# where Q_ii is the ice-interior conductive flux at T_si and A_tmp is the closed-form concentration. Q_ii cancels and the
# identity in 3 is exact only when ℒ does not depend on T. Fixture "constant_latent_heat" sets ρ_ℓ c_ℓ = ρ c exactly
# (1000 x 1890 = 900 x 2100, exact in Float32 and Float64), so the identity must hold at rounding. Fixture
# "production_phase" keeps ReadyESM's phase transitions (ρ_ℓ 999.8, c_ℓ 4186, c 2100). There the oracle reproduces the
# kernel's dV/dt with ℒ(T) and Q_ii exactly, and the residual must equal the oracle's residual. That residual is about
# 1.4% of dV/dt here. It is a property of the released closed form, not of h_0.
#
# Fixture. Production values from src/dynamic_ocean.jl (_build_sea_ice) and NumericalEarth's default snow: k_i 2 and
# k_s 0.31 W m-1 K-1, ρ_i 900 and ρ_s 330 kg m-3, h_c 0.05 m, ice salinity 4, Δt 900 s. Ocean salinity 34 gives
# T_b = -1.836 C. The snow surface is held at 0 C (melting snow). Heat fluxes are per cell area and positive upward.
# Q_u = -90 or -45 W m-2 at A = 0.6 or 0.3 puts 150 W m-2 per ice area into the snow surface. Downward conduction is only
# 3.2 or 1.3 W m-2, so the snow melts at Q_s = 147 or 149 W m-2 per ice area, far below its energy capacity (1.2e4 and
# 3.1e4 W m-2 for one step). Q_b = -150 or -120 W m-2 draws heat from the ice base (as frazil heat does), so
# α = (Q_u - Q_b)/(ρ_i ℒ₀) > 0. With α > 0 and β = Q_s/(ρ_i ℒ₀) > 0, both closed-form candidates give dV/dt > 0. The
# kernel therefore takes the freeze branch, and A rises. No clip, ridging cap or flooding acts (the test asserts this).
# Cells, in Oceananigans order (1,1) (2,1) (1,2) (2,2): cells 2 and 4 repeat cells 1 and 3 with snowfall.
# Production runs this kernel once per RK3 substep with the substep Δt. One call at 900 s tests the same arithmetic.
using ReadyESM, Test
const R = ReadyESM
const O = R.Oceananigans
const C = R.ClimaSeaIce
const CSI = C.SeaIceThermodynamics

const DT_SECONDS = 900
const HC = 0.05
const ICE_DENSITY = 900
const SNOW_DENSITY = 330
const ICE_CONDUCTIVITY = 2
const SNOW_CONDUCTIVITY = 0.31
const ICE_SALINITY = 4
const OCEAN_SALINITY = 34
const SNOW_SURFACE_TEMPERATURE = 0
const LIQUIDUS_SLOPE = 0.054 # ClimaSeaIce LinearLiquidus default, with a freshwater melting point of 0 C
const CELLS = (
    A = (0.6, 0.6, 0.3, 0.3),
    h = (0.5, 0.5, 1.2, 1.2),
    hs = (0.1, 0.1, 0.25, 0.25),
    Qu = (-90, -90, -45, -45),     # W m-2 per cell area, positive upward
    Qb = (-150, -150, -120, -120), # W m-2 per cell area, positive upward
    Ps = (0, 5e-5, 0, 5e-5),       # snowfall, kg m-2 s-1
)
const PRODUCTION_PHASE = (density = 900, heat_capacity = 2100, liquid_density = 999.8, liquid_heat_capacity = 4186,
    reference_latent_heat = 334e3, reference_temperature = 0)
const CONSTANT_LATENT_PHASE = merge(PRODUCTION_PHASE, (liquid_density = 1000, liquid_heat_capacity = 1890))
const FIXTURES = (("constant_latent_heat", CONSTANT_LATENT_PHASE, true), ("production_phase", PRODUCTION_PHASE, false))
const COMBOS = ((Float64, Float64), (Float32, Float32), (Float32, Float64)) # (grid number format, Δt type)

grid_for(FT) = O.RectilinearGrid(O.CPU(), FT; size = (2, 2, 1), extent = (2, 2, 1),
    topology = (O.Periodic, O.Periodic, O.Bounded))

# The same construction as ReadyESM's _sea_ice_simulation_with_momentum_forcing (src/dynamic_ocean.jl), the constructor
# that the option uses, with dynamics and advection off. That constructor takes the bottom salinity from the ocean, so
# this one takes it as an argument, with a choice of phase transitions. The second testset compares the two.
function build(FT, evolution, phase; salinity = OCEAN_SALINITY)
    grid = grid_for(FT)
    top_surface_temperature = O.Field{O.Center, O.Center, Nothing}(grid)
    ice_thermodynamics = C.sea_ice_slab_thermodynamics(grid;
        internal_heat_flux = CSI.ConductiveFlux(FT; conductivity = FT(ICE_CONDUCTIVITY)),
        top_heat_boundary_condition = CSI.PrescribedTemperature(top_surface_temperature.data),
        bottom_heat_boundary_condition = CSI.IceWaterThermalEquilibrium(FT(salinity)),
        concentration_evolution = evolution)
    top_heat_flux = O.Field{O.Center, O.Center, Nothing}(grid)
    return C.SeaIceModel(grid;
        clock = O.Clock(grid),
        ice_salinity = FT(ICE_SALINITY),
        advection = nothing,
        ice_consolidation_thickness = FT(HC),
        sea_ice_density = FT(ICE_DENSITY),
        snow_density = FT(SNOW_DENSITY),
        phase_transitions = C.PhaseTransitions(FT; map(FT, phase)...),
        ice_thermodynamics,
        snow_thermodynamics = R.NumericalEarth.SeaIces.default_snow_thermodynamics(grid),
        snowfall = similar(top_heat_flux),
        dynamics = nothing,
        bottom_heat_flux = similar(top_heat_flux),
        top_heat_flux,
        forcing = NamedTuple(),
        boundary_conditions = NamedTuple())
end

function check_build(model, FT, evolution, phase)
    pt = model.phase_transitions
    @test model.ice_thermodynamics.concentration_evolution === evolution
    @test model.ice_thermodynamics.heat_boundary_conditions.top isa CSI.PrescribedTemperature
    @test model.snow_thermodynamics.heat_boundary_conditions.top isa CSI.PrescribedTemperature
    @test model.ice_thermodynamics.internal_heat_flux.conductivity === FT(ICE_CONDUCTIVITY)
    @test model.snow_thermodynamics.internal_heat_flux.conductivity === FT(SNOW_CONDUCTIVITY)
    for key in keys(phase)
        @test getproperty(pt, key) === FT(getproperty(phase, key))
    end
    @test pt.liquidus.slope === FT(LIQUIDUS_SLOPE)
    @test iszero(pt.liquidus.freshwater_melting_temperature)
    return nothing
end

cells!(field, values) = (O.interior(field) .= reshape(collect(values), 2, 2, 1); field)
host(field) = vec(Array(O.interior(field)))

function set_fixture!(model, FT; Tu = SNOW_SURFACE_TEMPERATURE)
    cells!(model.ice_concentration, FT.(CELLS.A))
    cells!(model.ice_thickness, FT.(CELLS.h))
    cells!(model.snow_thickness, FT.(CELLS.hs))
    cells!(model.external_heat_fluxes.top, FT.(CELLS.Qu))
    cells!(model.external_heat_fluxes.bottom, FT.(CELLS.Qb))
    cells!(model.snowfall, FT.(CELLS.Ps))
    cells!(model.snow_thermodynamics.top_surface_temperature, ntuple(_ -> FT(Tu), 4))
    return model
end

# One call of the real layered kernel, then every output it writes.
function step!(model, Δt)
    CSI.thermodynamic_time_step!(model, model.ice_thermodynamics, model.snow_thermodynamics, Δt)
    fluxes = model.mass_fluxes
    return (A = host(model.ice_concentration), h = host(model.ice_thickness), hs = host(model.snow_thickness),
        ice = host(fluxes.thermodynamics.ice), snow = host(fluxes.thermodynamics.snow),
        intercepted = host(fluxes.intercepted_snowfall), Tsi = host(model.ice_thermodynamics.top_surface_temperature))
end

# Float64 scalar oracle for cell k. Inputs are the model's own values (FT), converted exactly to Float64.
function oracle(FT, phase, hl, k, Δt; Tu = SNOW_SURFACE_TEMPERATURE, salinity = OCEAN_SALINITY)
    x(v) = Float64(FT(v))
    A, hi, hs = x(CELLS.A[k]), x(CELLS.h[k]), x(CELLS.hs[k])
    Qu, Qb, Ps = x(CELLS.Qu[k]), x(CELLS.Qb[k]), x(CELLS.Ps[k])
    ρi, ρs, ki, ks, hc = x(ICE_DENSITY), x(SNOW_DENSITY), x(ICE_CONDUCTIVITY), x(SNOW_CONDUCTIVITY), x(HC)
    ℒ₀, T₀ = x(phase.reference_latent_heat), x(phase.reference_temperature)
    ρ, c, ρℓ, cℓ = x(phase.density), x(phase.heat_capacity), x(phase.liquid_density), x(phase.liquid_heat_capacity)
    Δt, hl, Tus = Float64(Δt), Float64(hl), x(Tu)
    Tb = -x(LIQUIDUS_SLOPE) * x(salinity)
    # Conduction through snow and ice in series (consolidated ice, hi ≥ hc), with the snow surface at Tus.
    Ri, Rs = hi / ki, hs / ks
    Tsi = Tb + (Tus - Tb) * Ri / (Rs + Ri)
    Qis = (Tb - Tus) / (Rs + Ri)    # per ice area, positive upward
    melt = max(0.0, Qis - Qu / A)   # energy left to melt snow, per ice area
    capacity = ρs * ℒ₀ * hs / Δt
    Qs = min(melt, capacity)
    Gs_melt, Gs_fall = Qs / (ρs * ℒ₀), Ps / ρs
    # Closed form: dV/dt = α + β A_next and A_next = A + K dV/dt, for each branch.
    E₀ = ρi * ℒ₀
    α, β = (Qu - Qb) / E₀, Qs / E₀
    Kf, Km = Δt * (1 - A) / hl, Δt * A / (2hi)
    A_freeze = (A + Kf * α) / (1 - Kf * β)
    A_melt = (A + Km * α) / (1 - Km * β)
    # Final dV/dt as ice_melt_freeze_tendency forms it, with A_tmp = A_freeze (see the header).
    Qii = -ki * (Tsi - Tb) / hi
    ℒ(T) = ℒ₀ + (ρℓ * cℓ / ρ - c) * (T - T₀)
    dVdt = (Qu + Qs * A_freeze - Qii) / (ρi * ℒ(Tsi)) + (Qii - Qb) / (ρi * ℒ(Tb))
    # Volume, then the freeze branch of the concentration rule at h_l, then the snow rebase, melt and snowfall.
    V0 = hi * A
    V1 = V0 + Δt * dVdt
    A1 = A + Δt * (1 - A) / hl * dVdt
    h1 = V1 / A1
    hs1 = hs * A / A1 + Δt * (Gs_fall - Gs_melt)
    freeboard = h1 * (1 - ρi / ρℓ) - hs1 * ρs / ρℓ
    return (; A, hi, hs, Qu, Qb, Ps, ρi, ρs, ℒ₀, hc, Tb, Tus, Tsi, Qis, Qii, Qs, capacity, Gs_melt, Gs_fall, α, β,
        Kf, Km, A_freeze, A_melt, dVdt, V0, V1, A1, h1, hs1, freeboard, Lmin = min(ℒ(Tsi), ℒ(Tb), ℒ₀),
        ice = ρi * dVdt, snow = ρs * (hs1 * A1 - hs * A) / Δt - Ps * A1, intercepted = Ps * A1,
        residual = dVdt - (α + β * A1))
end

# Rounding bounds: 16 eps of the model's number format times the size of the terms that enter each quantity. The volume
# update subtracts V0 from V1, so a volume rate carries an error of order eps V/Δt, and the concentration step
# multiplies that by K. In Float64 these bounds are many orders of magnitude below the effects tested (h_0, the
# fixed-point term and the branch). In Float32 they are below the h_0 effect.
function tolerances(FT, o, Δt)
    u = 16 * eps(FT)
    Δt = Float64(Δt)
    rate = (o.V0 + o.V1) / Δt + (abs(o.Qu) + abs(o.Qb) + 2 * abs(o.Qii) + 2 * o.Qs) / (o.ρi * o.Lmin)
    A = 1 + o.Kf * rate
    hs = o.hs + o.hs * A / o.A1 + Δt * (o.Gs_fall + o.Gs_melt)
    snow_volume = o.hs * o.A + o.hs1 * o.A1 + Δt * (o.Gs_fall + o.Gs_melt)
    return (; rate = u * rate, A = u * A, h = u * (o.h1 + (Δt * rate + o.h1 * A) / o.A1), hs = u * hs,
        ice = u * o.ρi * rate, snow = u * (o.ρs * (hs + o.hs * A) / Δt + o.Ps * (1 + A)),
        intercepted = u * o.Ps * (1 + A), Tsi = u * 2 * (abs(o.Tb) + abs(o.Tus)), residual = u * (rate + o.β * A),
        snow_volume = u * snow_volume, energy = u * o.ℒ₀ * (o.ρi * rate + o.ρs * snow_volume / Δt + o.Ps))
end

function check_cell(FT, out, o, Δt, k; exact_identity, report)
    tol = tolerances(FT, o, Δt)
    Δt = Float64(Δt)
    m(key) = Float64(getproperty(out, key)[k])
    A1, h1, hs1, ice, snow, intercepted, Tsi = m(:A), m(:h), m(:hs), m(:ice), m(:snow), m(:intercepted), m(:Tsi)
    # 1. The fixture: consolidated ice, snow melt below capacity, an unambiguous freeze branch, no clip, cap or flood.
    @test o.hi >= o.hc
    @test o.Qs > 0
    @test o.Qs < o.capacity
    @test o.α > 0
    @test o.β > 0
    @test o.α + o.β * o.A_melt > 0   # the melt-branch candidate gives dV/dt > 0, so the kernel rejects it
    @test o.α + o.β * o.A_freeze > 0 # the freeze-branch candidate is self-consistent
    @test o.dVdt > 0
    @test o.A < o.A1 < 1
    @test o.hs1 > 0
    @test o.freeboard > 0
    # 2. Returned state and mass fluxes against the oracle.
    @test abs(Tsi - o.Tsi) <= tol.Tsi
    @test abs(A1 - o.A1) <= tol.A
    @test abs(h1 - o.h1) <= tol.h
    @test abs(hs1 - o.hs1) <= tol.hs
    @test abs(ice - o.ice) <= tol.ice
    @test abs(snow - o.snow) <= tol.snow
    @test abs(intercepted - o.intercepted) <= tol.intercepted
    # 3. Fixed-point residual of the returned state, with dV/dt from the returned ice mass flux.
    dVdt = ice / o.ρi
    residual = dVdt - (o.α + o.β * A1)
    if exact_identity
        @test abs(residual) <= tol.residual
        # Column energy: latent heat of the ice gain less the snow melt equals the net heat loss Q_u - Q_b.
        @test abs(o.ℒ₀ * (ice + snow) - (o.Qu - o.Qb)) <= tol.energy
    else
        @test abs(residual - o.residual) <= tol.residual
    end
    # 4. Snow volume accounting: the rebase keeps hs A, so the change is melt and snowfall over the new concentration.
    @test abs((hs1 * A1 - o.hs * o.A) - Δt * (o.Gs_fall - o.Gs_melt) * A1) <= tol.snow_volume
    @test abs(snow + o.ρs * o.Gs_melt * A1) <= tol.snow
    @test abs(intercepted - o.Ps * A1) <= tol.intercepted
    # The fixed-point term and the branch choice are far above the bounds (resolved in Float64 only).
    naive = o.A + o.Kf * (o.α + o.β * o.A)
    if FT == Float64
        @test abs(A1 - naive) > 1000 * tol.A
        @test abs(dVdt - (o.α + o.β * o.A_melt)) > 1000 * tol.rate
    end
    println("SNOW_FIXED_POINT ", report, " A0=", o.A, " A1=", A1, " dA=", A1 - o.A, " dVdt=", dVdt, " residual=",
        residual, " residual_over_dVdt=", residual / dVdt, " fixed_point_shift=", A1 - naive)
    return nothing
end

@testset "Sea-ice lead closing: snow-covered freezing fixed point" begin
    @testset "Kernel against the oracle" begin
        for (name, phase, exact_identity) in FIXTURES, (FT, DT) in COMBOS
            @testset "$name FT=$FT Δt::$DT" begin
                Δt = DT(DT_SECONDS)
                if exact_identity
                    # ρ_ℓ c_ℓ/ρ - c is exactly zero in the model's number format, so ℒ(T) = ℒ₀.
                    @test FT(phase.liquid_density) * FT(phase.liquid_heat_capacity) / FT(phase.density) -
                          FT(phase.heat_capacity) === zero(FT)
                end
                rules = (("release", CSI.ProportionalEvolution(), FT(HC)),
                    ("h0=0.2", CSI.LeadClosingEvolution(FT(0.2)), FT(0.2)),
                    ("h0=0.5", CSI.LeadClosingEvolution(FT(0.5)), FT(0.5)),
                    ("h0=hc", CSI.LeadClosingEvolution(FT(HC)), FT(HC)))
                outputs = Dict{String, Any}()
                for (label, evolution, hl) in rules
                    model = build(FT, evolution, phase)
                    check_build(model, FT, evolution, phase)
                    outputs[label] = step!(set_fixture!(model, FT), Δt)
                end
                # 5. h_0 = h_c gives the released snow-path result, bit for bit, in every output.
                for key in keys(outputs["release"])
                    @test all(outputs["release"][key] .=== outputs["h0=hc"][key])
                end
                # Snowfall does not change the ice: cells 2 and 4 repeat the ice of cells 1 and 3.
                for label in ("release", "h0=0.2", "h0=0.5"), key in (:A, :h, :ice, :Tsi), (a, b) in ((1, 2), (3, 4))
                    @test outputs[label][key][a] === outputs[label][key][b]
                end
                # 1-4 for h_0 = 0.5 and 0.2 m, and for the release (h_l = h_c).
                for (label, evolution, hl) in rules[1:3], k in 1:4
                    @testset "$label cell $k" begin
                        o = oracle(FT, phase, hl, k, Δt)
                        check_cell(FT, outputs[label], o, Δt, k; exact_identity,
                            report = "fixture=$name FT=$FT DT=$DT rule=$label cell=$k")
                    end
                end
                # A larger h_0 closes leads more slowly: the concentration rise falls as h_0 grows.
                for k in 1:4
                    A0 = FT(CELLS.A[k])
                    @test outputs["release"].A[k] - A0 > outputs["h0=0.2"].A[k] - A0
                    @test outputs["h0=0.2"].A[k] - A0 > outputs["h0=0.5"].A[k] - A0
                    @test outputs["h0=0.5"].A[k] > A0
                end
            end
        end
    end

    @testset "Fixture construction matches ReadyESM's sea-ice constructor" begin
        for (FT, DT) in COMBOS
            @testset "FT=$FT Δt::$DT" begin
                evolution = CSI.LeadClosingEvolution(FT(0.5))
                simulation = R._sea_ice_simulation_with_momentum_forcing(grid_for(FT), nothing;
                    dynamics = nothing, advection = nothing, forcing = NamedTuple(),
                    conductivity = FT(ICE_CONDUCTIVITY),
                    internal_heat_flux = CSI.ConductiveFlux(FT; conductivity = FT(ICE_CONDUCTIVITY)),
                    ice_salinity = FT(ICE_SALINITY), ice_heat_capacity = FT(PRODUCTION_PHASE.heat_capacity),
                    ice_consolidation_thickness = FT(HC), sea_ice_density = FT(ICE_DENSITY),
                    snow_density = FT(SNOW_DENSITY), concentration_evolution = evolution)
                reference = simulation.model
                check_build(reference, FT, evolution, PRODUCTION_PHASE)
                # Without an ocean that constructor uses zero bottom salinity. A cold snow surface keeps conduction on.
                expected = step!(set_fixture!(reference, FT; Tu = -0.5), DT(DT_SECONDS))
                mirror = build(FT, evolution, PRODUCTION_PHASE; salinity = 0)
                check_build(mirror, FT, evolution, PRODUCTION_PHASE)
                actual = step!(set_fixture!(mirror, FT; Tu = -0.5), DT(DT_SECONDS))
                for key in keys(expected)
                    @test all(expected[key] .=== actual[key])
                end
            end
        end
    end
end
println("SEA_ICE_LEAD_SNOW_FIXED_POINT_TESTS_PASS")
end
