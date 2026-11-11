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

# --------------------------------------------------------------------------------------------------------------
# Fixed-column wrapper fixture (GPU-interface form), copied from test/tiedtke_bechtold_host.jl at e9a4133 and the
# stage 1B column test at afb714d (so this file does not depend on those test modules), with the stage 1B additions:
# the D81 flag in `prognostic.convection`, the four prescribed flux inputs in `prognostic.ocean` and `.land`, and a
# trade-wind column that convects shallow.
# --------------------------------------------------------------------------------------------------------------
const TROPICAL = (T_s = 300.0, sigma_m = 1.0, q_m = 0.0, lapse1 = 6.5e-3, sigma_i = 0.0, inversion = 0.0,
    lapse2 = 6.5e-3, T_trop = 195.0, nodes = ((0.9, 0.85), (0.6, 0.6), (0.2, 0.3)))
const MIDLATITUDE = (T_s = 298.0, sigma_m = 0.85, q_m = 0.0105, lapse1 = 6.5e-3, sigma_i = 0.0, inversion = 0.0,
    lapse2 = 6.5e-3, T_trop = 215.0, nodes = ((0.85, 0.7), (0.5, 0.5), (0.2, 0.3)))
const STABLE = (T_s = 260.0, sigma_m = 1.0, q_m = 0.0, lapse1 = 0.0, sigma_i = 0.0, inversion = 0.0,
    lapse2 = 0.0, T_trop = 100.0, nodes = ((1.0, 0.5),))
trade(sigma_i) = (T_s = 300.0, sigma_m = 0.89, q_m = 0.0128, lapse1 = 5.5e-3, sigma_i = sigma_i, inversion = 6.0,
    lapse2 = 5.0e-3, T_trop = 200.0,
    nodes = ((0.89, 0.9), (sigma_i + 0.01, 0.9), (sigma_i - 0.01, 0.25), (0.2, 0.2)))

function relative_humidity_profile(sigma::NF, nodes) where {NF}
    s0, r0 = NF(nodes[1][1]), NF(nodes[1][2])
    sigma >= s0 && return r0
    for i in 2:length(nodes)
        s1, r1 = NF(nodes[i][1]), NF(nodes[i][2])
        sigma >= s1 && return r0 + (r1 - r0) * (sigma - s0) / (s1 - s0)
        s0, r0 = s1, r1
    end
    return r0
end

function analytic_profile(::Type{NF}, p_full, p_s, spec, atmosphere, gravity) where {NF}
    kappa = atmosphere.R_dry / atmosphere.heat_capacity
    exponent1 = atmosphere.R_dry * NF(spec.lapse1) / gravity
    exponent2 = atmosphere.R_dry * NF(spec.lapse2) / gravity
    T_s, sigma_m, sigma_i, T_trop = NF(spec.T_s), NF(spec.sigma_m), NF(spec.sigma_i), NF(spec.T_trop)
    T_m = T_s * sigma_m^kappa
    T_i = sigma_i > 0 ? max(T_m * (sigma_i / sigma_m)^exponent1, T_trop) + NF(spec.inversion) : zero(NF)
    T = zeros(NF, length(p_full))
    q = zeros(NF, length(p_full))
    for k in eachindex(p_full)
        sigma = p_full[k] / NF(p_s)
        if sigma >= sigma_m
            T[k] = T_s * sigma^kappa
            q[k] = NF(spec.q_m)
        else
            T[k] = sigma >= sigma_i ? max(T_m * (sigma / sigma_m)^exponent1, T_trop) :
                                      max(T_i * (sigma / sigma_i)^exponent2, T_trop)
            q[k] = relative_humidity_profile(sigma, spec.nodes) * S.saturation_humidity(T[k], p_full[k], atmosphere)
        end
    end
    return T, q
end

function host_geopotential(T::Vector{NF}, q, p_full, p_half, surface_geopotential, atmosphere) where {NF}
    n = length(T)
    R_dry = atmosphere.R_dry
    mu = atmosphere.μ_virt_temp
    geopotential = zeros(NF, n)
    geopotential[n] = surface_geopotential + R_dry * T[n] * (1 + mu * q[n]) * log(p_half[n + 1] / p_full[n])
    for k in (n - 1):-1:1
        geopotential[k] = geopotential[k + 1] +
            R_dry * T[k + 1] * (1 + mu * q[k + 1]) * log(p_full[k + 1] / p_half[k + 1]) +
            R_dry * T[k] * (1 + mu * q[k]) * log(p_half[k + 1] / p_full[k])
    end
    return geopotential
end

# Columns: ocean and land inputs differ, so the D48 weighting is exercised in the partial-land column.
const HOST_COLUMNS_1B = (
    tropical = (profile = TROPICAL, p_s = 101325.0, land = 0.0),
    midlatitude = (profile = MIDLATITUDE, p_s = 100000.0, land = 1.0),
    trade = (profile = trade(0.75), p_s = 101325.0, land = 0.0),
    coast = (profile = trade(0.75), p_s = 101325.0, land = 0.4),
    stable = (profile = STABLE, p_s = 100000.0, land = 0.3),
)
const OCEAN_INPUT = (H = 15.0, E = 6.0e-5)     # W m-2, kg m-2 s-1, upward
const LAND_INPUT = (H = 40.0, E = 3.0e-5)

function host_fixture_1b(::Type{NF}, convection; prefill = zero(NF), set = 1) where {NF}
    n = 27
    sigma_half = [NF(j) / NF(n) for j in 0:n]
    sigma_full = [NF(0.5) * (sigma_half[k + 1] + sigma_half[k]) for k in 1:n]
    sigma_thick = [sigma_half[k + 1] - sigma_half[k] for k in 1:n]
    atmosphere = S.EarthAtmosphere(NF)
    gravity = S.Earth(NF).gravity
    npoints = length(HOST_COLUMNS_1B)
    T, q, phi = zeros(NF, npoints, n), zeros(NF, npoints, n), zeros(NF, npoints, n)
    ps, land, latitude = zeros(NF, npoints), zeros(NF, npoints), zeros(NF, npoints)
    for (ij, column) in enumerate(HOST_COLUMNS_1B)
        ps[ij] = NF(column.p_s)
        p_full = sigma_full .* ps[ij]
        p_half = sigma_half .* ps[ij]
        Tc, qc = analytic_profile(NF, p_full, column.p_s, column.profile, atmosphere, gravity)
        T[ij, :] = Tc
        q[ij, :] = qc
        phi[ij, :] = host_geopotential(Tc, qc, p_full, p_half, zero(NF), atmosphere)
        land[ij] = NF(column.land)
    end
    allocate(v, value) = v.dims isa S.Grid3D ? fill(value, npoints, n) :
                         v.dims isa S.Grid2D ? fill(value, npoints) : error("unexpected dims $(v.dims)")
    declared = S.variables(convection)
    scratch = (; (v.name => allocate(v, prefill) for v in declared if v isa S.ScratchVariable)...)
    parameterizations = (; (v.name => allocate(v, zero(NF)) for v in declared if v isa S.ParameterizationVariable)...)
    fill!(parameterizations.cloud_top, NF(n + 1))
    fill!(parameterizations.convective_temperature_tendency, prefill)
    fill!(parameterizations.convective_humidity_tendency, prefill)
    prognostic = (;
        convection = (; tiedtke_bechtold_lagged_flux_set = fill(NF(set), npoints)),
        ocean = (; sensible_heat_flux = fill(NF(OCEAN_INPUT.H), npoints),
                   surface_humidity_flux = fill(NF(OCEAN_INPUT.E), npoints)),
        land = (; sensible_heat_flux = fill(NF(LAND_INPUT.H), npoints),
                  surface_humidity_flux = fill(NF(LAND_INPUT.E), npoints)),
    )
    vars = (;
        grid = (; temperature_prev = T, humidity_prev = q, geopotential = phi, pressure_prev = ps),
        tendencies = (; grid = (; temperature = zeros(NF, npoints, n), humidity = zeros(NF, npoints, n))),
        parameterizations,
        prognostic,
        scratch = (; grid = scratch),
    )
    model = (;
        geometry = (; σ_levels_full = sigma_full, σ_levels_half = sigma_half, σ_levels_thick = sigma_thick,
                    latds = latitude),
        time_stepping = (; Δt_sec = NF(900)),
        land_sea_mask = (; mask = land),
    )
    return (; vars, model, n)
end

"""The stage 1B column on fixture column `ij` with plain vectors (the reference), and the lagged value it used."""
function direct_column_1b(vars, model, ij, convection)
    parameters = convection.convection.parameters
    NF = typeof(parameters.gravity)
    n = length(model.geometry.σ_levels_full)
    ps = vars.grid.pressure_prev[ij]
    f = NF(model.land_sea_mask.mask[ij])
    # The D48 weighting written out independently of `_pbl_area_weighted` (land first, then ocean, zero-area tiles
    # skipped), so the reference does not reuse the code under test.
    weighted(land_value, ocean_value) = (f > 0 ? f * NF(land_value) : zero(NF)) + (f < 1 ? (1 - f) * NF(ocean_value) : zero(NF))
    lagged = vars.prognostic.convection.tiedtke_bechtold_lagged_flux_set[ij] > 0 ?
        R.TiedtkeBechtoldLaggedSurfaceFlux(weighted(LAND_INPUT.H, OCEAN_INPUT.H), weighted(LAND_INPUT.E, OCEAN_INPUT.E)) :
        R.TiedtkeBechtoldLaggedSurfaceFlux{NF}()
    dT, dq, P = zeros(NF, n), zeros(NF, n), zeros(NF, n + 1)
    work = R.tiedtke_bechtold_column_work(NF, n)
    result = R.tiedtke_bechtold_column!(dT, dq, P, work,
        NF[vars.grid.temperature_prev[ij, k] for k in 1:n], NF[vars.grid.humidity_prev[ij, k] for k in 1:n],
        NF[model.geometry.σ_levels_full[k] * ps for k in 1:n], NF[model.geometry.σ_levels_half[j] * ps for j in 1:(n + 1)],
        NF[vars.grid.geopotential[ij, k] for k in 1:n], model.time_stepping.Δt_sec,
        R.tb_land_flag(model.land_sea_mask.mask[ij], parameters), model.geometry.latds[ij], parameters,
        convection.shallow, lagged, R.TiedtkeBechtoldHostClosure())
    return (; dT, dq, P, result, lagged)
end
row(field, ij) = [field[ij, k] for k in 1:size(field, 2)]
function wrapper_allocations(vars, convection, model, ij)
    S.parameterization!(ij, vars, convection, model)
    return @allocated S.parameterization!(ij, vars, convection, model)
end

@testset "Stage 1B wrapper equals the plain column; D81 set/unset; NaN scratch; no allocation; inference $NF" for NF in PRECISIONS
    convection = R.TiedtkeBechtoldStage1BHostConvection(R.TiedtkeBechtoldConvection(R.tiedtke_bechtold_parameters(NF)),
                                                        R.tiedtke_bechtold_shallow_parameters(NF))
    types_seen = Set{Int}()
    for set in (0, 1)
        clean = host_fixture_1b(NF, convection; set)
        poisoned = host_fixture_1b(NF, convection; set, prefill = NF(NaN))
        for ij in 1:length(HOST_COLUMNS_1B)
            reference = direct_column_1b(clean.vars, clean.model, ij, convection)
            r = reference.result
            set == 1 && push!(types_seen, r.convection_type)
            for fixture in (clean, poisoned)
                S.parameterization!(ij, fixture.vars, convection, fixture.model)
                p = fixture.vars.parameterizations
                @test same_bits(row(p.convective_temperature_tendency, ij), reference.dT)
                @test same_bits(row(p.convective_humidity_tendency, ij), reference.dq)
                @test same_bits(row(fixture.vars.tendencies.grid.temperature, ij), reference.dT)   # tendencies start at 0
                @test same_bits(row(fixture.vars.tendencies.grid.humidity, ij), reference.dq)
                @test bits(p.tiedtke_bechtold_rain_mass_flux_surface[ij]) == bits(r.rain_mass_flux_surface)
                @test p.tiedtke_bechtold_convection_type[ij] == r.convection_type
                @test p.tiedtke_bechtold_trigger_path[ij] == r.trigger_path == (set == 1 ? 2 : 1)
                @test bits(p.tiedtke_bechtold_lagged_sensible_heat_flux[ij]) == bits(reference.lagged.sensible_heat_flux)
                @test bits(p.tiedtke_bechtold_lagged_humidity_flux[ij]) == bits(reference.lagged.humidity_flux)
                @test bits(p.tiedtke_bechtold_shallow_source[ij]) == bits(r.shallow_source)
                @test bits(p.tiedtke_bechtold_shallow_closure_mass_flux[ij]) == bits(r.shallow_closure_mass_flux)
                @test bits(p.tiedtke_bechtold_cloud_base_mass_flux[ij]) == bits(r.cloud_base_mass_flux)
                @test p.tiedtke_bechtold_lagged_flux_unset_steps[ij] == r.lagged_flux_unset          # one call
                @test p.tiedtke_bechtold_shallow_cfl_capped_steps[ij] == r.shallow_cfl_capped
                @test p.cloud_top[ij] == (r.convection_type == 1 ? r.tb_top_full_level : clean.n + 1)  # D79 (a)
                @test fixture.vars.prognostic.convection.tiedtke_bechtold_lagged_flux_set[ij] == 1  # D81: set after
                # A second call on a set column runs stage 1B and adds the counters once more.
                fixture === clean && set == 0 && begin
                    S.parameterization!(ij, fixture.vars, convection, fixture.model)
                    @test p.tiedtke_bechtold_trigger_path[ij] == 2
                    @test p.tiedtke_bechtold_lagged_flux_unset_steps[ij] == 1
                end
            end
        end
    end
    @test 2 in types_seen && 1 in types_seen          # the fixture exercises shallow and deep columns
    # The CPU form of the GPU interface: no allocation and fully inferred.
    clean = host_fixture_1b(NF, convection)
    for ij in 1:length(HOST_COLUMNS_1B)
        @test wrapper_allocations(clean.vars, convection, clean.model, ij) == 0
    end
    @test (@inferred R._tiedtke_bechtold_lagged_flux(1, clean.vars, clean.model, NF)) isa
          R.TiedtkeBechtoldLaggedSurfaceFlux{NF}
    @test (@inferred S.parameterization!(1, clean.vars, convection, clean.model)) === nothing
    @test isbits(convection)
end

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
    @test R._validate_convection_restart(sim(tb1b), (; convection_restart_contract = NamedTuple{reverse(keys(contract))}(reverse(values(contract))))) ===
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
        # D81. At a fresh start the first coupled iteration holds two atmosphere steps (SpeedyWeather's first Euler half
        # step, then the leapfrog step; the land coupling period is split the same way, terrarium_gpu_coupling.jl
        # l. 1178-1190). So the first TB call is unset (stage 1a fallback) and sets the flag, and the second TB call in
        # iteration 1 already runs stage 1B. 65411434 showed exactly this; v1 of this test wrongly expected the fallback
        # in the iteration-1 record. Before any step (the iteration-0 callback, if it runs) every flag and counter is 0.
        initial = filter(a -> a.iteration == 0, activity)
        for a in initial
            @test all(iszero, a.flag) && all(iszero, a.unset_steps)
        end
        @test all(==(1), stepped[1].unset_steps)                  # exactly one unset call per column, ever
        @test all(==(1), stepped[1].flag)
        @test all(==(2), stepped[1].trigger_path)                 # the record is the second (stage 1B) call
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
