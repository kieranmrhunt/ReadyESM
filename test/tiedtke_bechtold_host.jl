module TiedtkeBechtoldHostTests
# Tiedtke-Bechtold stage 1a, part 2: host integration (src/tiedtke_bechtold_host.jl and its config, factory and restart
# hooks) against the frozen contract tiedtke_bechtold_claude_v1/frozen_v4/FROZEN_CONTRACT_v4.sha256 b07062d2
# (DESIGN.md "Integration gate after Paul's code audit" items 1-6; STAGE1_CHOICES_DRAFT.md r15 D43, T4, T12, T14).
# Tests: configuration name and every refusal (with the host's prognostic-cloud selector); the declared work arrays
# and TB-named diagnostics; the wrapper on fixed columns in Float64 and Float32 (T12 unit boundary, accumulation, the
# host water budget in kg m-2, the never-written correction, bit identity with the plain column call, NaN-filled
# scratch, no allocation, inference); T14 host cloud top; a real constructed model (live constants, declared arrays,
# a real column); T4 restart (contract, refusals with and without the identity override, TB to TB acceptance and
# identity); and a short coupled analytic CPU run with TB selected whose diagnostics validate.
using ReadyESM, Test, Dates
const R = ReadyESM
const S = R.SpeedyWeather
const O = R.Oceananigans

fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
with(config; kwargs...) = R.ExperimentConfig(; merge(fields(config), (; kwargs...))...)
with_forcing(config; kwargs...) = with(config; forcing = R.ForcingConfig(; merge(fields(config.forcing), (; kwargs...))...))
include(joinpath(@__DIR__, "analytic_coupled_fixture.jl"))

const TB = :tiedtke_bechtold_eqcape_stage1
const PRECISIONS = (Float64, Float32)
const MANIFEST = "b07062d259455a207fb45afd1d2a0283387c167f10f8e7a57bb150d36b8c0e76"
# Declared tolerance of column identities (as in the part-1 test): 64 n eps(NF) times the gross exchange.
closure_tolerance(::Type{NF}, n) where {NF} = 64 * n * eps(NF)

bits(x::Float32) = reinterpret(UInt32, x)
bits(x::Float64) = reinterpret(UInt64, x)
bits(x::Complex) = (bits(real(x)), bits(imag(x)))
same_bits(a, b) = size(a) == size(b) && eltype(a) == eltype(b) && all(bits(x) == bits(y) for (x, y) in zip(a, b))
prognostic(vars) = map(n -> copy(Array(S.get_step(getproperty(vars.prognostic, n), 1).data)),
                       (:temperature, :humidity, :vorticity, :divergence, :pressure))

# --------------------------------------------------------------------------------------------------------------
# Analytic columns, copied from test/tiedtke_bechtold_column.jl at 42712b4 (so this file does not depend on that test
# module). On the equal-sigma L27 grid below they are the part-1 tropical (deep, water) and mid-latitude (deep, land)
# cases bit for bit, and a stable column (no deep cloud).
# --------------------------------------------------------------------------------------------------------------
const TROPICAL = (T_s = 300.0, sigma_m = 1.0, q_m = 0.0, lapse1 = 6.5e-3, sigma_i = 0.0, inversion = 0.0,
    lapse2 = 6.5e-3, T_trop = 195.0, nodes = ((0.9, 0.85), (0.6, 0.6), (0.2, 0.3)))
const MIDLATITUDE = (T_s = 298.0, sigma_m = 0.85, q_m = 0.0105, lapse1 = 6.5e-3, sigma_i = 0.0, inversion = 0.0,
    lapse2 = 6.5e-3, T_trop = 215.0, nodes = ((0.85, 0.7), (0.5, 0.5), (0.2, 0.3)))
const STABLE = (T_s = 260.0, sigma_m = 1.0, q_m = 0.0, lapse1 = 0.0, sigma_i = 0.0, inversion = 0.0,
    lapse2 = 0.0, T_trop = 100.0, nodes = ((1.0, 0.5),))

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
            q[k] = relative_humidity_profile(sigma, spec.nodes) *
                   S.saturation_humidity(T[k], p_full[k], atmosphere)
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

const HOST_COLUMNS = (
    tropical = (profile = TROPICAL, p_s = 101325.0, z_s = 0.0, land = 0.0, latitude = 5.0, deep = true),
    midlatitude = (profile = MIDLATITUDE, p_s = 100000.0, z_s = 0.0, land = 1.0, latitude = 45.0, deep = true),
    stable = (profile = STABLE, p_s = 100000.0, z_s = 0.0, land = 0.3, latitude = -60.0, deep = false),
)

host_wrapper(::Type{NF}) where {NF} =
    R.TiedtkeBechtoldHostConvection(R.TiedtkeBechtoldConvection(R.tiedtke_bechtold_parameters(NF)))

"""
Host fixture with the layout of SpeedyWeather 0.21.1 Variables as the GPU kernel sees them (npoints x nlayers
matrices and npoints vectors in place of the Fields): equal-sigma L27 geometry, the columns above, every declared
scratch and parameterization variable allocated from the wrapper's own `variables` declaration (`prefill` for the
scratch and the process tendencies, which the column must overwrite), host cloud top n+1 as after
SpeedyWeather's reset, and a 900 s step.
"""
function host_fixture(::Type{NF}, convection; prefill = zero(NF)) where {NF}
    n = 27
    sigma_half = [NF(j) / NF(n) for j in 0:n]
    sigma_full = [NF(0.5) * (sigma_half[k + 1] + sigma_half[k]) for k in 1:n]   # SpeedyWeather Geometry formula
    sigma_thick = [sigma_half[k + 1] - sigma_half[k] for k in 1:n]
    atmosphere = S.EarthAtmosphere(NF)
    gravity = S.Earth(NF).gravity
    npoints = length(HOST_COLUMNS)
    T, q, phi = zeros(NF, npoints, n), zeros(NF, npoints, n), zeros(NF, npoints, n)
    ps, land, latitude = zeros(NF, npoints), zeros(NF, npoints), zeros(NF, npoints)
    for (ij, column) in enumerate(HOST_COLUMNS)
        ps[ij] = NF(column.p_s)
        p_full = sigma_full .* ps[ij]
        p_half = sigma_half .* ps[ij]
        Tc, qc = analytic_profile(NF, p_full, column.p_s, column.profile, atmosphere, gravity)
        T[ij, :] = Tc
        q[ij, :] = qc
        phi[ij, :] = host_geopotential(Tc, qc, p_full, p_half, gravity * NF(column.z_s), atmosphere)
        land[ij] = NF(column.land)
        latitude[ij] = NF(column.latitude)
    end
    allocate(v, value) = v.dims isa S.Grid3D ? fill(value, npoints, n) :
                         v.dims isa S.Grid2D ? fill(value, npoints) : error("unexpected dims $(v.dims)")
    declared = S.variables(convection)
    scratch = (; (v.name => allocate(v, prefill) for v in declared if v isa S.ScratchVariable)...)
    parameterizations = (; (v.name => allocate(v, zero(NF)) for v in declared if v isa S.ParameterizationVariable)...)
    fill!(parameterizations.cloud_top, NF(n + 1))
    fill!(parameterizations.convective_temperature_tendency, prefill)
    fill!(parameterizations.convective_humidity_tendency, prefill)
    vars = (;
        grid = (; temperature_prev = T, humidity_prev = q, geopotential = phi, pressure_prev = ps),
        tendencies = (; grid = (; temperature = zeros(NF, npoints, n), humidity = zeros(NF, npoints, n))),
        parameterizations,
        scratch = (; grid = scratch),
    )
    model = (;
        geometry = (; σ_levels_full = sigma_full, σ_levels_half = sigma_half, σ_levels_thick = sigma_thick,
                    latds = latitude),
        time_stepping = (; Δt_sec = NF(900)),
        land_sea_mask = (; mask = land),
    )
    return (; vars, model, sigma_full, sigma_half, sigma_thick, n)
end

"""The pure column on host column `ij` with plain vectors and the part-1 work arrays (the reference)."""
function direct_column(vars, model, ij, parameters; sigma_full = model.geometry.σ_levels_full,
                       sigma_half = model.geometry.σ_levels_half)
    NF = eltype(parameters.gravity)
    n = length(sigma_full)
    ps = vars.grid.pressure_prev[ij]
    dT, dq, P = zeros(NF, n), zeros(NF, n), zeros(NF, n + 1)
    work = R.tiedtke_bechtold_column_work(NF, n)
    result = R.tiedtke_bechtold_column!(dT, dq, P, work,
        NF[vars.grid.temperature_prev[ij, k] for k in 1:n], NF[vars.grid.humidity_prev[ij, k] for k in 1:n],
        NF[sigma_full[k] * ps for k in 1:n], NF[sigma_half[j] * ps for j in 1:(n + 1)],
        NF[vars.grid.geopotential[ij, k] for k in 1:n], model.time_stepping.Δt_sec,
        R.tb_land_flag(model.land_sea_mask.mask[ij], parameters), model.geometry.latds[ij], parameters)
    return (; dT, dq, P, work, result)
end

"""A half-level work array of column ij read back from the scratch pair (top half level, then the lower ones)."""
half_level(scratch, name, ij) = vcat(getproperty(scratch, R.tiedtke_bechtold_scratch_top_name(name))[ij],
    [getproperty(scratch, R.tiedtke_bechtold_scratch_name(name))[ij, k]
     for k in 1:size(getproperty(scratch, R.tiedtke_bechtold_scratch_name(name)), 2)])
row(field, ij) = [field[ij, k] for k in 1:size(field, 2)]

function wrapper_allocations(vars, convection, model, ij)
    S.parameterization!(ij, vars, convection, model)
    return @allocated S.parameterization!(ij, vars, convection, model)
end

# Host fields the wrapper writes, per column.
const HOST_FIELDS_2D = (:rain_convection, :rain_rate_convection, :rain_rate, :cloud_top,
    :convective_precipitation_surface_pressure, :convective_precipitation_mass_correction)
const HOST_FIELDS_3D = (:convective_humidity_tendency, :convective_temperature_tendency)

@testset "Tiedtke-Bechtold stage 1a host integration" begin

production = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)
# Production uses mixed-phase saturation and a 7200 s Betts-Miller time scale, which TB refuses (see below).
tb_production = with(production; atmosphere_convection = TB, atmosphere_saturation_law = :liquid,
                     atmosphere_convection_timescale_seconds = 14400)

@testset "Configuration: name, refusals, prognostic-cloud selector, factory" begin
    @test R.validate(tb_production; check_input_files = false) === tb_production
    @test R.ExperimentConfig().atmosphere_convection == :betts_miller_constant_rh     # TB is opt-in
    @test production.atmosphere_convection == :betts_miller_sigma_rh_v1
    # The name parses from YAML.
    mktempdir() do dir
        lines = readlines(joinpath(R.PROJECT_ROOT, "config", "production.yml"))
        replacements = ("atmosphere_convection:" => "atmosphere_convection: tiedtke_bechtold_eqcape_stage1",
                        "atmosphere_saturation_law:" => "atmosphere_saturation_law: liquid",
                        "atmosphere_convection_timescale_seconds:" => "atmosphere_convection_timescale_seconds: 14400")
        for (prefix, line) in replacements
            index = findall(l -> startswith(l, prefix), lines)
            @test length(index) == 1
            lines[only(index)] = line
        end
        path = joinpath(dir, "tb.yml")
        write(path, join(lines, "\n") * "\n")
        @test R.load_config(path; check_input_files = false).atmosphere_convection == TB
    end
    # Every other unknown value is still rejected.
    for name in (:tiedtke_bechtold, :tiedtke_bechtold_v1, :tiedtke_bechtold_eqcape_stage1b, :tiedtke_bechtold_stage1,
                 :dilute_cape_adjustment_v1, :Tiedtke_Bechtold_eqcape_stage1)
        @test_throws ArgumentError R.validate(with(tb_production; atmosphere_convection = name);
                                              check_input_files = false)
    end
    # The host prognostic-cloud selector (release 4caab68 has one member) is used by the builder and the refusal.
    @test R._PROGNOSTIC_CLOUD_SCHEMES == (:prognostic_condensate,)
    @test !R._uses_prognostic_cloud_condensate(:diagnostic) && !R._uses_prognostic_cloud_condensate(:none)
    for scheme in R._PROGNOSTIC_CLOUD_SCHEMES
        @test R._uses_prognostic_cloud_condensate(scheme)
        cloudy = (; cloud_scheme = scheme, cloud_condensate_retention_fraction = 0.5)
        bm = with_forcing(production; cloudy...)
        @test R.validate(bm; check_input_files = false) === bm     # the same cloud settings are valid with BM
        @test_throws ArgumentError R.validate(with_forcing(tb_production; cloudy...); check_input_files = false)
        @test_throws r"prognostic cloud condensate" R.validate(with_forcing(tb_production; cloudy...);
                                                                check_input_files = false)
    end
    for scheme in (:diagnostic, :none)
        config = with_forcing(tb_production; cloud_scheme = scheme)
        @test R.validate(config; check_input_files = false) === config
    end
    # Liquid saturation only (D01, D02), no moist-mass coupling, and no Betts-Miller or C08/C07 setting.
    @test_throws ArgumentError R.validate(with(tb_production; atmosphere_saturation_law = :mixed_phase);
                                          check_input_files = false)
    @test_throws ArgumentError R.validate(with(tb_production; atmosphere_moist_mass_coupling = true);
                                          check_input_files = false)
    changed = (atmosphere_convection_timescale_seconds = 7200, atmosphere_convection_upper_relative_humidity = 0.4,
        atmosphere_convection_lower_relative_humidity = 0.8, atmosphere_convection_transition_top_sigma = 0.25,
        atmosphere_convection_transition_bottom_sigma = 0.75, atmosphere_c08_entrainment_scale = 0.5,
        atmosphere_c08_precipitation_conversion_m_inv = 0.002, atmosphere_c08_condensate_threshold_kg_kg = 0.0004,
        atmosphere_c07_momentum_strength = 0.25, atmosphere_c07_maximum_subcycles = 8)
    @test Set(keys(changed)) == Set(R._TIEDTKE_BECHTOLD_INERT_CONVECTION_SETTINGS)
    for name in keys(changed)
        @test getfield(R.ExperimentConfig(), name) != changed[name]
        @test_throws ArgumentError R.validate(with(tb_production; NamedTuple{(name,)}((changed[name],))...);
                                              check_input_files = false)
    end
    # Factory: TB gets its own wrapper; the Betts-Miller paths are unchanged.
    grid = S.SpectralGrid(; NF = Float32, trunc = 7, nlayers = 8)
    convection = R._atmosphere_convection(tb_production, grid)
    @test convection isa R.TiedtkeBechtoldHostConvection
    @test convection.convection isa R.TiedtkeBechtoldConvection{Float32}
    @test !occursin("BettsMiller", string(typeof(convection)))
    @test convection.convection.parameters === R.tiedtke_bechtold_parameters(Float32; truncation = 7)
    sigma_bm = R._atmosphere_convection(production, grid)
    @test sigma_bm isa R.NetColumnBettsMillerConvection{<:R.SigmaProfileBettsMillerConvection}
    constant_bm = R._atmosphere_convection(with(production; atmosphere_convection = :betts_miller_constant_rh), grid)
    @test constant_bm isa R.NetColumnBettsMillerConvection{<:S.BettsMillerConvection}
    bm_names = Dict(v.name => v.desc for v in S.variables(constant_bm))
    @test bm_names[:convective_humidity_tendency] == "Incremental Betts-Miller specific-humidity tendency"
    @test !any(name -> startswith(String(name), "tiedtke_bechtold"), keys(bm_names))
    # Restart identity: the convection symbol is bound by the config hash (not a run-control field).
    @test :atmosphere_convection ∉ R._DYNAMIC_RESTART_RUN_CONTROL_FIELDS
    @test R._dynamic_config_identity_hash(tb_production) !=
          R._dynamic_config_identity_hash(with(tb_production; atmosphere_convection = :betts_miller_constant_rh))
    # Provenance: none for Betts-Miller, TB-named values for TB, no Betts-Miller label.
    @test R._tiedtke_bechtold_metadata(production) == (;)
    meta = R._tiedtke_bechtold_metadata(tb_production)
    @test meta.atmosphere_convection_contract_manifest_sha256 == MANIFEST
    @test meta.atmosphere_convection_water_contract == "instantaneous_rain_no_rup_v1"
    @test !any(value -> occursin("Betts", string(value)), values(meta))
end

@testset "Declarations: scratch work arrays, host fields, TB names, isbits $NF" for NF in PRECISIONS
    convection = host_wrapper(NF)
    @test isbits(convection)
    @test isbits(R.Adapt.adapt(Array, convection)) && R.Adapt.adapt(Array, convection) === convection
    declared = S.variables(convection)
    names = [v.name for v in declared]
    @test allunique(names)
    find(name) = only(filter(v -> v.name == name, declared))
    scratch = filter(v -> v isa S.ScratchVariable, declared)
    @test length(scratch) == length(R._TIEDTKE_BECHTOLD_LAYER_WORK) + 2 * length(R._TIEDTKE_BECHTOLD_HALF_LEVEL_WORK)
    @test all(v -> v.namespace == :grid, scratch)
    for (name, _, _) in R._TIEDTKE_BECHTOLD_LAYER_WORK
        v = find(R.tiedtke_bechtold_scratch_name(name))
        @test v isa S.ScratchVariable && v.dims isa S.Grid3D
    end
    for (name, _, _) in R._TIEDTKE_BECHTOLD_HALF_LEVEL_WORK
        lower, top = find(R.tiedtke_bechtold_scratch_name(name)), find(R.tiedtke_bechtold_scratch_top_name(name))
        @test lower isa S.ScratchVariable && lower.dims isa S.Grid3D
        @test top isa S.ScratchVariable && top.dims isa S.Grid2D
    end
    # The scratch covers exactly the part-1 column work arrays plus the rain-flux profile.
    work_names = Set([first(x) for x in (R._TIEDTKE_BECHTOLD_LAYER_WORK..., R._TIEDTKE_BECHTOLD_HALF_LEVEL_WORK...)])
    @test work_names == union(Set(keys(R.tiedtke_bechtold_column_work(NF, 4))), Set([:rain_mass_flux_half]))
    for name in (HOST_FIELDS_2D..., (first(x) for x in R._TIEDTKE_BECHTOLD_DIAGNOSTICS)...)
        @test find(name) isa S.ParameterizationVariable && find(name).dims isa S.Grid2D
    end
    for name in HOST_FIELDS_3D
        @test find(name) isa S.ParameterizationVariable && find(name).dims isa S.Grid3D
    end
    @test !any(v -> occursin("Betts", v.desc), declared)          # TB-named, no Betts-Miller labels
    @test :snow_rate_convection ∉ names && :snow_convection ∉ names   # stage 1a is rain only
    # The scheme alone is not a model component: only the host wrapper has a column method.
    @test !hasmethod(S.parameterization!, Tuple{Int, Any, typeof(convection.convection), Any})
    @test hasmethod(S.parameterization!, Tuple{Int, Any, typeof(convection), Any})
    @test R._uses_tiedtke_bechtold_convection(convection)
    @test R._uses_tiedtke_bechtold_convection(R.PrognosticCondensateConvection(convection, nothing))
    @test !R._uses_tiedtke_bechtold_convection(convection.convection)
end

@testset "Wrapper on fixed columns: T12 unit boundary, water budget, never-written correction $NF" for NF in PRECISIONS
    convection = host_wrapper(NF)
    parameters = convection.convection.parameters
    fixture = host_fixture(NF, convection)
    (; vars, model, n) = fixture
    p = vars.parameterizations
    rho = NF(R._FRESHWATER_DENSITY_KG_M3)
    @test rho == 1000
    dt = model.time_stepping.Δt_sec
    g = parameters.gravity
    tol = closure_tolerance(NF, n)
    large_scale_rain = NF(3e-6)
    fill!(p.convective_precipitation_mass_correction, NF(7))   # sentinel: TB must never write it
    fill!(p.rain_rate, large_scale_rain)                       # large-scale condensation runs before convection
    for (ij, name) in enumerate(keys(HOST_COLUMNS))
        deep = HOST_COLUMNS[name].deep
        reference = direct_column(vars, model, ij, parameters)
        S.parameterization!(ij, vars, convection, model)
        flux = p.tiedtke_bechtold_rain_mass_flux_surface[ij]
        rate = p.rain_rate_convection[ij]
        @test (p.tiedtke_bechtold_convection_type[ij] == 1) == deep
        @test bits(flux) == bits(reference.result.rain_mass_flux_surface)
        @test (flux > 0) == deep
        @test rate == flux / rho                                   # D43: the single unit boundary
        @test abs(rate * rho - flux) <= 4 * eps(NF) * flux         # T12: flux = host rate x _FRESHWATER_DENSITY_KG_M3
        @test deep ? rate < flux : rate == 0                       # the density factor is applied
        @test p.rain_convection[ij] == dt * rate                   # T12: accumulation from zero
        @test p.rain_rate[ij] == large_scale_rain + rate           # host total rain gains exactly the convective rate
        @test p.convective_precipitation_surface_pressure[ij] == vars.grid.pressure_prev[ij]
        @test p.convective_precipitation_mass_correction[ij] == 7  # never written
        # The tendencies started at zero, so they equal the column's dT and dq, which are also the process diagnostics.
        @test same_bits(row(vars.tendencies.grid.temperature, ij), zero(NF) .+ reference.dT)
        @test same_bits(row(vars.tendencies.grid.humidity, ij), zero(NF) .+ reference.dq)
        @test same_bits(row(p.convective_temperature_tendency, ij), reference.dT)
        @test same_bits(row(p.convective_humidity_tendency, ij), reference.dq)
        # T12: the host water budget over the step closes in kg m-2 (host sigma thickness, host fields only).
        ps = vars.grid.pressure_prev[ij]
        dq_host = row(p.convective_humidity_tendency, ij)
        vapour_loss = -ps / g * sum(dq_host .* fixture.sigma_thick) * dt
        surface_rain = rho * p.rain_convection[ij]
        budget_gross = ps / g * sum(abs.(dq_host) .* fixture.sigma_thick) * dt + surface_rain
        @test abs(vapour_loss - surface_rain) <= tol * budget_gross
        # DESIGN item 3: the scheme's own surface rain equals the net column water loss (column dp), no repair.
        dp = [fixture.sigma_half[k + 1] * ps - fixture.sigma_half[k] * ps for k in 1:n]
        ledger = flux + sum(dq_host .* dp) / g
        @test abs(ledger) <= tol * (flux + sum(abs.(dq_host) .* dp) / g)
        # From a nonzero accumulator the increment is dt times the rate, bit for bit.
        before = NF(0.25)
        p.rain_convection[ij] = before
        S.parameterization!(ij, vars, convection, model)
        @test p.rain_convection[ij] == before + dt * rate
        @test p.rain_rate_convection[ij] == rate
        @test p.convective_precipitation_mass_correction[ij] == 7
    end
end

@testset "Wrapper equals the plain column; NaN-filled scratch; no allocation; inference $NF" for NF in PRECISIONS
    convection = host_wrapper(NF)
    parameters = convection.convection.parameters
    clean = host_fixture(NF, convection)
    poisoned = host_fixture(NF, convection; prefill = NF(NaN))
    for ij in 1:length(HOST_COLUMNS)
        reference = direct_column(clean.vars, clean.model, ij, parameters)
        for fixture in (clean, poisoned)
            S.parameterization!(ij, fixture.vars, convection, fixture.model)
            scratch = fixture.vars.scratch.grid
            p = fixture.vars.parameterizations
            @test same_bits(half_level(scratch, :rain_mass_flux_half, ij), reference.P)
            for (name, _, _) in R._TIEDTKE_BECHTOLD_HALF_LEVEL_WORK
                name == :rain_mass_flux_half && continue
                @test same_bits(half_level(scratch, name, ij), getproperty(reference.work, name))
            end
            for (name, _, _) in R._TIEDTKE_BECHTOLD_LAYER_WORK
                @test same_bits(row(getproperty(scratch, R.tiedtke_bechtold_scratch_name(name)), ij),
                                getproperty(reference.work, name))
            end
            r = reference.result
            @test p.tiedtke_bechtold_top_full_level[ij] == r.tb_top_full_level
            @test p.tiedtke_bechtold_base_half_level[ij] == r.cloud_base_half_level
            @test p.tiedtke_bechtold_departure_level[ij] == r.departure_level
            @test p.tiedtke_bechtold_convection_type[ij] == r.convection_type
            @test bits(p.tiedtke_bechtold_cloud_base_mass_flux[ij]) == bits(r.cloud_base_mass_flux)
            @test bits(p.tiedtke_bechtold_closure_factor[ij]) == bits(r.closure_factor)
            @test bits(p.tiedtke_bechtold_cfl_factor[ij]) == bits(r.f_CFL)
            @test bits(p.tiedtke_bechtold_positivity_factor[ij]) == bits(r.f_pos)
            @test bits(p.tiedtke_bechtold_final_factor[ij]) == bits(r.final_factor)
            @test bits(p.tiedtke_bechtold_pcape[ij]) == bits(r.PCAPE)
            @test bits(p.tiedtke_bechtold_timescale[ij]) == bits(r.tau)
            @test bits(p.tiedtke_bechtold_max_test_cloud_depth[ij]) == bits(r.max_test_cloud_depth)
            @test p.tiedtke_bechtold_model_top_reached[ij] == r.model_top_reached
            @test bits(p.tiedtke_bechtold_condensation[ij]) == bits(r.condensation_column)
        end
        # NaN-filled scratch and process arrays give bit-identical tendencies and host fields (write before read).
        for name in (:temperature, :humidity)
            @test same_bits(row(getproperty(clean.vars.tendencies.grid, name), ij),
                            row(getproperty(poisoned.vars.tendencies.grid, name), ij))
        end
        for name in HOST_FIELDS_2D
            @test bits(getproperty(clean.vars.parameterizations, name)[ij]) ==
                  bits(getproperty(poisoned.vars.parameterizations, name)[ij])
        end
        for name in HOST_FIELDS_3D
            @test same_bits(row(getproperty(clean.vars.parameterizations, name), ij),
                            row(getproperty(poisoned.vars.parameterizations, name), ij))
        end
    end
    # The CPU form of the GPU interface: no allocation and fully inferred.
    for ij in 1:length(HOST_COLUMNS)
        @test wrapper_allocations(clean.vars, convection, clean.model, ij) == 0
    end
    dT = view(clean.vars.parameterizations.convective_temperature_tendency, 1, :)
    dq = view(clean.vars.parameterizations.convective_humidity_tendency, 1, :)
    @test (@inferred R._tiedtke_bechtold_column_call!(dT, dq, 1, clean.vars, convection.convection, clean.model)) isa
          NamedTuple
    @test (@inferred S.parameterization!(1, clean.vars, convection, clean.model)) === nothing
end

@testset "T14 host cloud top $NF" for NF in PRECISIONS
    n = 27
    top = R.tiedtke_bechtold_host_cloud_top
    @test top(NF(n + 1), n + 1, 0, n) === NF(n + 1)     # no deep cloud leaves n+1
    @test top(NF(2), n + 1, 0, n) === NF(2)             # and leaves any existing top unchanged
    @test top(NF(n + 1), 3, 1, n) === NF(3)             # a TB top of 3 sets 3
    @test top(NF(2), 3, 1, n) === NF(2)                 # an existing higher top of 2 stays 2
    @test top(NF(5), 3, 1, n) === NF(3)
    for (tb_top, kind) in ((0, 1), (n + 2, 1), (0, 0), (n + 2, 0), (n + 1, 1), (3, 0), (3, 2), (3, -1))
        @test_throws ArgumentError top(NF(n + 1), tb_top, kind, n)   # 0 and n+2 (and inconsistent pairs) refused
    end
    # Through the wrapper: the deep tropical column and the stable column.
    convection = host_wrapper(NF)
    fixture = host_fixture(NF, convection)
    (; vars, model) = fixture
    p = vars.parameterizations
    S.parameterization!(1, vars, convection, model)
    tb_top = Int(p.tiedtke_bechtold_top_full_level[1])
    @test 1 <= tb_top <= n
    @test tb_top < p.tiedtke_bechtold_base_half_level[1]   # a top above the base, not the base
    @test p.cloud_top[1] == tb_top
    for preset in (2, tb_top + 1, n + 1)
        p.cloud_top[1] = NF(preset)
        S.parameterization!(1, vars, convection, model)
        @test p.cloud_top[1] == min(preset, tb_top)
    end
    stable = findfirst(==(:stable), keys(HOST_COLUMNS))
    for preset in (n + 1, 2)
        p.cloud_top[stable] = NF(preset)
        S.parameterization!(stable, vars, convection, model)
        @test p.tiedtke_bechtold_convection_type[stable] == 0
        @test p.tiedtke_bechtold_top_full_level[stable] == n + 1
        @test p.cloud_top[stable] == preset
    end
end

@testset "Constructed model: live constants, declared arrays, a real column" begin
    config = R.ExperimentConfig(; truncation = 7, nlayers = 27, duration_days = 1, atmosphere_convection = TB,
        forcing = R.ForcingConfig(radiation = :rrtmgp_clear_sky))
    @test R.validate(config; check_input_files = false) === config
    simulation = R.build_rrtmgp(config)
    model, vars = simulation.model, simulation.variables
    convection = model.convection
    @test convection isa R.TiedtkeBechtoldHostConvection
    @test convection.convection.parameters === R.tiedtke_bechtold_parameters(model)
    @test S.initialize!(convection, model) === nothing
    changed = R.tiedtke_bechtold_parameters(Float32; atmosphere = S.EarthAtmosphere(Float32; heat_capacity = 1005),
                                            truncation = 7)
    @test_throws ArgumentError S.initialize!(R.TiedtkeBechtoldHostConvection(R.TiedtkeBechtoldConvection(changed)), model)
    npoints, n = model.spectral_grid.npoints, model.geometry.nlayers
    @test n == 27
    scratch = vars.scratch.grid
    for (name, _, _) in R._TIEDTKE_BECHTOLD_LAYER_WORK
        @test size(getproperty(scratch, R.tiedtke_bechtold_scratch_name(name))) == (npoints, n)
    end
    for (name, _, _) in R._TIEDTKE_BECHTOLD_HALF_LEVEL_WORK
        @test size(getproperty(scratch, R.tiedtke_bechtold_scratch_name(name))) == (npoints, n)
        @test size(getproperty(scratch, R.tiedtke_bechtold_scratch_top_name(name))) == (npoints,)
    end
    for (name, _, _) in R._TIEDTKE_BECHTOLD_DIAGNOSTICS
        @test haskey(vars.parameterizations, name)
    end
    @test !haskey(vars.parameterizations, :snow_rate_convection)
    # A real column: the tropical profile on the model's own sigma levels, written into column 1.
    ij = 1
    sigma_full, sigma_half = model.geometry.σ_levels_full, model.geometry.σ_levels_half
    ps = Float32(101325)
    p_full, p_half = sigma_full .* ps, sigma_half .* ps
    T, q = analytic_profile(Float32, p_full, 101325.0, TROPICAL, model.atmosphere, model.planet.gravity)
    phi = host_geopotential(T, q, p_full, p_half, 0.0f0, model.atmosphere)
    g, p = vars.grid, vars.parameterizations
    g.temperature_prev.data[ij, :] .= T
    g.humidity_prev.data[ij, :] .= q
    g.geopotential.data[ij, :] .= phi
    g.pressure_prev.data[ij] = ps
    vars.tendencies.grid.temperature.data[ij, :] .= 0
    vars.tendencies.grid.humidity.data[ij, :] .= 0
    p.cloud_top.data[ij] = n + 1
    p.rain_rate.data[ij] = 0
    p.rain_convection.data[ij] = 0
    p.convective_precipitation_mass_correction.data[ij] = 7
    reference = direct_column(vars, model, ij, convection.convection.parameters)
    S.parameterization!(ij, vars, convection, model)
    @test reference.result.convection_type == 1 && p.tiedtke_bechtold_convection_type[ij] == 1
    @test same_bits(row(vars.tendencies.grid.temperature, ij), 0.0f0 .+ reference.dT)
    @test same_bits(row(vars.tendencies.grid.humidity, ij), 0.0f0 .+ reference.dq)
    @test same_bits(half_level(scratch, :rain_mass_flux_half, ij), reference.P)
    flux = p.tiedtke_bechtold_rain_mass_flux_surface[ij]
    @test bits(flux) == bits(reference.result.rain_mass_flux_surface) && flux > 0
    @test p.rain_rate_convection[ij] == flux / Float32(R._FRESHWATER_DENSITY_KG_M3)
    @test p.rain_rate[ij] == p.rain_rate_convection[ij]
    @test p.rain_convection[ij] == model.time_stepping.Δt_sec * p.rain_rate_convection[ij]
    @test p.cloud_top[ij] == reference.result.tb_top_full_level
    @test p.convective_precipitation_mass_correction[ij] == 7
    println("TIEDTKE_BECHTOLD_REAL_COLUMN top=", reference.result.tb_top_full_level, " base=",
            reference.result.cloud_base_half_level, " rain_mm_per_day=", 86_400 * flux, " field_view_allocations=",
            wrapper_allocations(vars, convection, model, ij))
end

@testset "T4 restart contract: record, legacy rule, refusals not waived by the override (standalone)" begin
    contract = R.tiedtke_bechtold_restart_contract()
    @test contract == (; schema = "readiesm_convection_restart_contract_v1", scheme = "tiedtke_bechtold_eqcape_stage1",
        stage_version = "stage1a", water_contract = "instantaneous_rain_no_rup_v1",
        trigger_variant = "deep_from_n_minus_2", parcel_adapter = "none", contract_manifest_sha256 = MANIFEST)
    grid = S.SpectralGrid(; NF = Float32, trunc = 7, nlayers = 8)
    sim(convection) = (; model = (; convection))
    tb = R.TiedtkeBechtoldHostConvection(grid)
    bm = R.NetColumnBettsMillerConvection(S.BettsMillerConvection(grid))
    record = (; convection_restart_contract = contract)
    @test R._convection_restart_state(sim(bm)) == (;)            # Betts-Miller checkpoints do not change
    @test R._convection_restart_state(sim(tb)) == record
    @test R._validate_convection_restart(sim(tb), record) === nothing
    @test R._validate_convection_restart(sim(bm), (;)) === nothing          # legacy rule: no record with Betts-Miller
    @test_throws ArgumentError R._validate_convection_restart(sim(bm), record)   # TB -> BM
    @test_throws ArgumentError R._validate_convection_restart(sim(tb), (;))      # BM (or legacy) -> TB
    for key in keys(contract)
        @test_throws ArgumentError R._validate_convection_restart(sim(tb),
            (; convection_restart_contract = merge(contract, NamedTuple{(key,)}(("changed",)))))
        @test_throws ArgumentError R._validate_convection_restart(sim(tb),
            (; convection_restart_contract = Base.structdiff(contract, NamedTuple{(key,)})))
    end
    @test_throws ArgumentError R._validate_convection_restart(sim(tb),
        (; convection_restart_contract = merge(contract, (; extra = "x"))))
    @test_throws ArgumentError R._validate_convection_restart(sim(tb),
        (; convection_restart_contract = merge(contract, (; scheme = TB))))      # a Symbol is not the string record
    reordered = NamedTuple{reverse(keys(contract))}(reverse(values(contract)))
    @test R._validate_convection_restart(sim(tb), (; convection_restart_contract = reordered)) === nothing
    @test_throws ArgumentError R.with_dynamic_restart_identity_override(() -> R._validate_convection_restart(sim(tb), (;)))
    @test_throws ArgumentError R.with_dynamic_restart_identity_override(() -> R._validate_convection_restart(sim(bm), record))
    # Real standalone atmospheres: save, reload through the checkpoint serializer, restore.
    tb_config = R.ExperimentConfig(; truncation = 7, nlayers = 8, duration_days = 1, atmosphere_convection = TB,
        forcing = R.ForcingConfig(radiation = :rrtmgp_clear_sky))
    bm_config = with(tb_config; atmosphere_convection = :betts_miller_constant_rh)
    tb_sim, bm_sim = R.build_rrtmgp(tb_config), R.build_rrtmgp(bm_config)
    tb_state, bm_state = O.prognostic_state(tb_sim), O.prognostic_state(bm_sim)
    @test tb_state.convection_restart_contract == contract
    @test !hasproperty(bm_state, :convection_restart_contract)
    @test Set(keys(tb_state)) == union(Set(keys(bm_state)), Set([:convection_restart_contract]))
    @test haskey(tb_state.variables.scratch.grid, :tiedtke_bechtold_mass_flux_half)   # TB work arrays are saved
    @test !haskey(bm_state.variables.scratch.grid, :tiedtke_bechtold_mass_flux_half)
    marker = tb_sim.variables.scratch.grid.tiedtke_bechtold_mass_flux_half
    fill!(marker.data, 0.125f0)
    mktempdir() do dir
        tb_path = R.save_dynamic_restart_state(tb_sim, joinpath(dir, "tb.jld2"))
        bm_path = R.save_dynamic_restart_state(bm_sim, joinpath(dir, "bm.jld2"))
        tb_loaded = O.OutputWriters.load_checkpoint_state(tb_path; base_path = "simulation")
        bm_loaded = O.OutputWriters.load_checkpoint_state(bm_path; base_path = "simulation")
        @test R._validate_convection_restart(tb_sim, tb_loaded) === nothing      # the reloaded record validates
        @test R._validate_convection_restart(bm_sim, bm_loaded) === nothing
        receiver_tb, receiver_bm = R.build_rrtmgp(tb_config), R.build_rrtmgp(bm_config)
        for (receiver, source, pattern) in ((receiver_bm, tb_loaded, r"carries a convection restart contract"),
                                            (receiver_tb, bm_loaded, r"no convection restart contract"))
            before = prognostic(receiver.variables)
            @test_throws ArgumentError O.restore_prognostic_state!(receiver, source)
            @test_throws ArgumentError R.with_dynamic_restart_identity_override(
                () -> O.restore_prognostic_state!(receiver, source))
            @test_throws pattern O.restore_prognostic_state!(receiver, source)
            @test all(same_bits(x, y) for (x, y) in zip(prognostic(receiver.variables), before))
        end
        @test all(iszero, receiver_tb.variables.scratch.grid.tiedtke_bechtold_mass_flux_half.data)
        @test O.restore_prognostic_state!(receiver_tb, tb_loaded) === receiver_tb          # TB -> TB accepted
        @test all(==(0.125f0), receiver_tb.variables.scratch.grid.tiedtke_bechtold_mass_flux_half.data)
        @test all(same_bits(x, y) for (x, y) in zip(prognostic(receiver_tb.variables), prognostic(tb_sim.variables)))
        @test O.restore_prognostic_state!(receiver_bm, bm_loaded) === receiver_bm          # BM -> BM unchanged
    end
end

@testset "Coupled analytic CPU run with TB: diagnostics, rain, restart identity and refusals" begin
    mktempdir() do dir
        base = with(analytic_coupled_test_config(production; output_dir = mkpath(joinpath(dir, "base")));
            atmosphere_timestep_at_t31_minutes = 3.75, duration_days = 4 * 900 / 86_400,
            start_date = DateTime(1993, 1, 1))
        output(name) = mkpath(joinpath(dir, name))
        tb_config = with(base; atmosphere_convection = TB, atmosphere_saturation_law = :liquid,
                         atmosphere_convection_timescale_seconds = 14400, output_dir = output("reference"))
        @test R.validate(tb_config; check_input_files = false) === tb_config
        reference = R.build_dynamic_esm(tb_config)
        atmosphere = reference.model.atmosphere
        convection = atmosphere.model.convection
        @test convection isa R.TiedtkeBechtoldHostConvection
        @test convection.convection.parameters === R.tiedtke_bechtold_parameters(atmosphere.model)
        path = joinpath(dir, "tb_step2.jld2")
        activity = NamedTuple[]
        function observe(sim)
            p = sim.model.atmosphere.variables.parameterizations
            push!(activity, (; iteration = sim.model.clock.iteration,
                deep_columns = count(==(1), Array(p.tiedtke_bechtold_convection_type.data)),
                deep_test_clouds = count(>(0), Array(p.tiedtke_bechtold_departure_level.data)),
                max_test_cloud_depth = Float64(maximum(Array(p.tiedtke_bechtold_max_test_cloud_depth.data))),
                max_rain_rate = Float64(maximum(Array(p.rain_rate_convection.data)))))
            sim.model.clock.iteration == 2 && R.save_dynamic_restart_state(sim, path)
            return nothing
        end
        O.add_callback!(reference, observe, O.IterationInterval(1); name = :tiedtke_bechtold_observer)
        O.run!(reference)
        @test reference.model.clock.iteration == 4 && isfile(path)
        final = prognostic(atmosphere.variables)
        @test all(x -> all(isfinite, x), final)
        diagnostics = R.collect_dynamic_diagnostics(reference, tb_config)
        R.save_dynamic_diagnostics(diagnostics, tb_config; render_figure = false)
        R.validate_dynamic_diagnostics(diagnostics)       # includes the convective rain against the layer tendency
        meta = diagnostics.metadata
        @test meta.atmosphere_convection == "tiedtke_bechtold_eqcape_stage1"
        @test meta.atmosphere_convective_precipitation_scheme ==
              "Tiedtke-Bechtold column surface rain mass flux over water density"
        @test meta.atmosphere_convection_contract_manifest_sha256 == MANIFEST
        @test meta.atmosphere_convection_reference_humidity_profile == "not_applicable"
        @test isnan(meta.atmosphere_convection_upper_relative_humidity)
        @test !occursin("Betts", meta.atmosphere_precipitation_mass_correction_diagnostics)
        @test all(iszero, diagnostics.atmosphere_convective_precipitation_mass_correction_flux)
        p = atmosphere.variables.parameterizations
        rho = Float32(R._FRESHWATER_DENSITY_KG_M3)
        @test all(iszero, Array(p.convective_precipitation_mass_correction.data))     # never written
        @test same_bits(Array(p.rain_rate_convection.data), Array(p.tiedtke_bechtold_rain_mass_flux_surface.data) ./ rho)
        # Column water ledger of the last step at every point, with the surface pressure the column used.
        dq = Float32.(Array(p.convective_humidity_tendency.data))
        ps = Array(p.convective_precipitation_surface_pressure.data)
        dsigma = Array(atmosphere.model.geometry.σ_levels_thick)
        gravity = atmosphere.model.planet.gravity
        flux = Array(p.tiedtke_bechtold_rain_mass_flux_surface.data)
        n = length(dsigma)
        worst = maximum(abs(flux[ij] + ps[ij] / gravity * sum(dq[ij, :] .* dsigma)) /
                        max(flux[ij] + ps[ij] / gravity * sum(abs.(dq[ij, :]) .* dsigma), floatmin(Float32))
                        for ij in eachindex(flux))
        @test worst <= closure_tolerance(Float32, n)
        stepped = filter(a -> a.iteration >= 1, activity)
        @test Set(a.iteration for a in stepped) == Set(1:4)
        deep_column_steps = sum(a.deep_columns for a in stepped)
        deep_test_cloud_steps = sum(a.deep_test_clouds for a in stepped)
        accumulated = Array(p.rain_convection.data)
        deep_depth = Float64(convection.convection.parameters.deep_depth)
        summary = "deep_column_steps=$deep_column_steps deep_test_cloud_steps=$deep_test_cloud_steps " *
            "columns_with_rain=$(count(>(0), accumulated)) of $(length(accumulated)) " *
            "max_rain_mm_per_day=$(86_400 * 1000 * maximum(a.max_rain_rate for a in stepped)) " *
            "max_test_cloud_depth_pa=$(maximum(a.max_test_cloud_depth for a in stepped)) deep_depth_pa=$deep_depth"
        if any(>(0), accumulated)
            println("TIEDTKE_BECHTOLD_COUPLED_RAIN branch=positive ", summary)
            @test deep_column_steps > 0
            @test maximum(a.max_rain_rate for a in stepped) > 0
        elseif deep_column_steps == 0 && deep_test_cloud_steps == 0
            # Stated explanation (HOST_NOTES.md): the analytic Jablonowski-Williamson state (about 5 K/km, RH 70%) is
            # stable for the 0.2 K / 0.1 g/kg deep test parcel from the admissible L8 departures, so no test cloud is
            # deeper than the 20000 Pa deep threshold (D14) in any step and TB returns exact zeros (D16).
            println("TIEDTKE_BECHTOLD_COUPLED_RAIN branch=no_deep_test_cloud ", summary)
            @test all(a.max_test_cloud_depth <= deep_depth for a in stepped)
            @test all(a.max_rain_rate == 0 for a in stepped)
            @test all(iszero, Array(p.convective_humidity_tendency.data))
            @test all(iszero, Array(p.convective_temperature_tendency.data))
        else
            # Not foreseen by the contract (deep test clouds with a zero closure, or deep convection without surface
            # rain): reported and failed, not explained away.
            println("TIEDTKE_BECHTOLD_COUPLED_RAIN branch=unforeseen ", summary)
            @test deep_column_steps == 0 && deep_test_cloud_steps == 0
        end
        # T4 restart identity of a TB run: restore the step-2 checkpoint into a fresh build and continue to step 4.
        continued = R.build_dynamic_esm(with(tb_config; output_dir = output("continued")))
        R.restore_dynamic_restart_state!(continued, path)
        @test continued.model.clock.iteration == 2
        continued.stop_iteration = 4
        continued.running = true
        while continued.running
            O.time_step!(continued)
        end
        @test continued.model.clock.iteration == 4
        @test all(same_bits(x, y) for (x, y) in zip(prognostic(continued.model.atmosphere.variables), final))
        # TB -> BM and BM -> TB: the identity check refuses first; with it waived, the convection contract refuses.
        bm_run = R.build_dynamic_esm(with(base; output_dir = output("bm")))
        @test bm_run.model.atmosphere.model.convection isa R.NetColumnBettsMillerConvection
        bm_path = joinpath(dir, "bm_step1.jld2")
        function save_bm(sim)
            sim.model.clock.iteration == 1 && R.save_dynamic_restart_state(sim, bm_path)
            return nothing
        end
        O.add_callback!(bm_run, save_bm, O.IterationInterval(1); name = :bm_restart_save)
        bm_run.stop_iteration = 1
        O.run!(bm_run)
        @test isfile(bm_path)
        before = prognostic(bm_run.model.atmosphere.variables)
        @test_throws ErrorException R.restore_dynamic_restart_state!(bm_run, path)
        @test_throws ArgumentError R.restore_dynamic_restart_state!(bm_run, path; allow_identity_mismatch = true)
        @test_throws r"carries a convection restart contract" R.restore_dynamic_restart_state!(bm_run, path;
                                                                                             allow_identity_mismatch = true)
        @test bm_run.model.clock.iteration == 1
        @test all(same_bits(x, y) for (x, y) in zip(prognostic(bm_run.model.atmosphere.variables), before))
        R.restore_dynamic_restart_state!(bm_run, bm_path)        # the Betts-Miller restart path still accepts its own
        @test bm_run.model.clock.iteration == 1
        fresh = R.build_dynamic_esm(with(tb_config; output_dir = output("fresh")))
        initial = prognostic(fresh.model.atmosphere.variables)
        @test_throws ErrorException R.restore_dynamic_restart_state!(fresh, bm_path)
        @test_throws ArgumentError R.restore_dynamic_restart_state!(fresh, bm_path; allow_identity_mismatch = true)
        @test_throws r"no convection restart contract" R.restore_dynamic_restart_state!(fresh, bm_path;
                                                                                      allow_identity_mismatch = true)
        @test fresh.model.clock.iteration == 0
        @test all(same_bits(x, y) for (x, y) in zip(prognostic(fresh.model.atmosphere.variables), initial))
    end
end

end
println("TIEDTKE_BECHTOLD_HOST_TESTS_PASS")
end
