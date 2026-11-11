module SeaIceBranchObserverTests
# Read-only sea-ice branch observer v1 (src/sea_ice_branch_observer.jl; bundle
# artifacts/development_wave_20261103/sea_ice_branch_observer_claude_v1/README.md), on the 4-step analytic coupled CPU
# case of test/sea_ice_lead_closing.jl:
# 1. settings: default off, validated, a run-control field outside the restart identity, inert hooks;
# 2. passivity: no observer, observer off and observer on end bit-identical in every prognostic field, and the
#    fifteen v019 prognostic hashes agree at every step. The same holds on the h0 = 0.5 m lead-closing arm;
# 3. one sample per coupled step (final Runge-Kutta stage only), finite states with 0 <= A <= 1 and h, hs >= 0, the
#    boundary records equal the live post-step state, the hemisphere sums equal an independent host recomputation
#    from analytic cell areas, and the ice and two-term snow mass ledgers close at rounding;
# 4. with SEA_ICE_BRANCH_OBSERVER_PARENT_REFERENCE set (the launcher sets it): observer off and on, 24 steps, equal
#    the sealed off-path hashes of the parent source 95459f5 (job 64431839), which has no observer code at all.
using ReadyESM, Test, Dates, SHA
const R = ReadyESM
const O = R.Oceananigans
const NC = R.NCDatasets

fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
with(config; kwargs...) = R.ExperimentConfig(; merge(fields(config), (; kwargs...))...)
include(joinpath(@__DIR__, "analytic_coupled_fixture.jl"))

const PRODUCTION = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)
analytic(steps; kwargs...) = with(analytic_coupled_test_config(PRODUCTION; output_dir = mktempdir());
    atmosphere_timestep_at_t31_minutes = 3.75, duration_days = steps * 900 / 86_400,
    start_date = DateTime(1993, 1, 1), kwargs...)

host(field) = Array(Base.has_offset_axes(field.data) ? parent(field.data) : field.data)
digest(values) = bytes2hex(sha256(reinterpret(UInt8, vec(values))))
interior_digest(field) = digest(Array(O.interior(field)))
function surface(field)
    values = Array(O.interior(field))
    return reshape(values, size(values, 1), size(values, 2))
end

# The off-path receipt of sea_ice_lead_closing_claude/offpath_lead.jl: the fifteen standard fields, every Terrarium
# land field and the grid temperature and humidity.
function offpath_hashes(simulation)
    p = simulation.model.atmosphere.variables.prognostic
    hashes = Dict{String, String}()
    for k in (:temperature, :humidity, :vorticity, :divergence, :pressure)
        hashes[String(k)] = digest(host(getproperty(p, k)))
    end
    ocean, ice = simulation.model.ocean.model, simulation.model.sea_ice.model
    for (k, f) in (("ocean_T", ocean.tracers.T), ("ocean_S", ocean.tracers.S), ("ocean_u", ocean.velocities.u),
                   ("ocean_v", ocean.velocities.v), ("ocean_eta", ocean.free_surface.displacement),
                   ("ice_h", ice.ice_thickness), ("ice_A", ice.ice_concentration), ("snow_h", ice.snow_thickness),
                   ("ice_u", ice.velocities.u), ("ice_v", ice.velocities.v))
        hashes[k] = interior_digest(f)
    end
    land = simulation.model.atmosphere.variables.prognostic.land.terrarium
    is_field(name) = try getproperty(land, name) isa O.Fields.AbstractField catch; false end
    for name in propertynames(land)
        is_field(name) && (hashes["land_" * String(name)] = interior_digest(getproperty(land, name)))
    end
    grid = simulation.model.atmosphere.variables.grid
    hashes["grid_temperature"] = digest(host(grid.temperature_prev))
    hashes["grid_humidity"] = digest(host(grid.humidity_prev))
    return hashes
end

# The off-path receipt plus every ocean and sea-ice prognostic field, the sea-ice mass fluxes, the ice top
# temperature and the clock.
function state_hashes(simulation)
    hashes = offpath_hashes(simulation)
    ocean, ice = simulation.model.ocean.model, simulation.model.sea_ice.model
    for (prefix, model) in (("ocean_prognostic_", ocean), ("ice_prognostic_", ice))
        prognostic = O.prognostic_fields(model)
        for name in keys(prognostic)
            hashes[prefix * String(name)] = interior_digest(getproperty(prognostic, name))
        end
    end
    hashes["ice_mass_flux_ice"] = interior_digest(ice.mass_fluxes.thermodynamics.ice)
    hashes["ice_mass_flux_snow"] = interior_digest(ice.mass_fluxes.thermodynamics.snow)
    hashes["ice_intercepted_snowfall"] = interior_digest(ice.mass_fluxes.intercepted_snowfall)
    hashes["ice_top_surface_temperature"] = interior_digest(ice.ice_thermodynamics.top_surface_temperature)
    hashes["clock"] = string(simulation.model.clock.iteration, ":", simulation.model.clock.time)
    return hashes
end

# Every arm gets the same passivity recorder (a read-only step-boundary callback).
function run_arm(config; install = nothing)
    simulation = R.build_dynamic_esm(config)
    isnothing(install) || install(simulation)
    registered = !isnothing(R.sea_ice_branch_observer(simulation))
    recorder = R.install_sea_ice_branch_passivity_hashes!(simulation)
    ice = simulation.model.sea_ice.model
    initial = (A = surface(ice.ice_concentration), h = surface(ice.ice_thickness), hs = surface(ice.snow_thickness))
    O.run!(simulation)
    return (simulation = simulation, config = config, recorder = recorder, initial = initial,
            registered = registered, final = state_hashes(simulation))
end

# Independent geometry of the analytic grid (32 x 16, longitude 0-360, latitude -80 to 80, Oceananigans radius).
const NX, NY = 32, 16
const RADIUS = 6.371e6
latitude_center(j) = -85.0 + 10.0 * j
analytic_area(j) = RADIUS^2 * (2π / NX) * (sind(-80.0 + 10.0 * j) - sind(-90.0 + 10.0 * j))
longitude_center(i) = 5.625 + 11.25 * (i - 1)
const AREA = [analytic_area(j) for i in 1:NX, j in 1:NY]
const SOUTH = [latitude_center(j) <= -50 for i in 1:NX, j in 1:NY]
const NORTH = .!SOUTH
const SECTOR = [SOUTH[i, j] && 120 <= longitude_center(i) < 180 for i in 1:NX, j in 1:NY]
region_sum(mask, values) = sum(AREA[I] * Float64(values[I]) for I in eachindex(values) if mask[I]; init = 0.0)

function read_record(path)
    return NC.NCDataset(joinpath(path, "sea_ice_branch_observer_v1.nc")) do ds
        values = Dict{String, Any}(name => Array(ds[name]) for name in keys(ds))
        attributes = Dict{String, Any}(name => ds.attrib[name] for name in keys(ds.attrib))
        (values = values, attributes = attributes, regions = String.(split(attributes["region_names"], ",")))
    end
end

function check_saved_files(path)
    lines = collect(eachline(joinpath(path, "SHA256SUMS")))
    @test length(lines) == 3
    for line in lines
        sha, name = split(line, "  ")
        @test bytes2hex(open(sha256, joinpath(path, name))) == sha
    end
    receipt = Dict{String, String}()
    for line in eachline(joinpath(path, "receipt.txt"))
        key, value = split(line, " = "; limit = 2)
        receipt[key] = value
    end
    @test receipt["schema"] == "sea_ice_branch_observer_v1"
    @test receipt["field_schema_version"] == "1"
    @test receipt["observer_mutates_live_state"] == "0"
    @test receipt["observer_code_sha256"] ==
          bytes2hex(open(sha256, joinpath(R.PROJECT_ROOT, "src", "sea_ice_branch_observer.jl")))
    @test receipt["observer_hook_code_sha256"] ==
          bytes2hex(open(sha256, joinpath(R.PROJECT_ROOT, "src", "conservative_sea_ice_advection.jl")))
    @test occursin("final Runge-Kutta stage", receipt["cadence"])
    @test occursin("frazil", receipt["not_in_v1"])
    @test receipt["result_pass"] == "true"
    return receipt
end

@testset "Sea-ice branch observer: settings, restart identity and inert hooks" begin
    default = R.ExperimentConfig()
    @test default.sea_ice_branch_observer === false
    @test PRODUCTION.sea_ice_branch_observer === false
    @test :sea_ice_branch_observer in R._DYNAMIC_RESTART_RUN_CONTROL_FIELDS
    on = with(PRODUCTION; sea_ice_branch_observer = true)
    @test R.validate(on; check_input_files = false) === on
    @test R._dynamic_config_identity_hash(on) == R._dynamic_config_identity_hash(PRODUCTION)
    bad = with(PRODUCTION; sea_ice_branch_observer = true, sea_ice_dynamics = false)
    message = try
        R.validate(bad; check_input_files = false)
        "accepted"
    catch err
        err isa ArgumentError ? err.msg : "unexpected $(typeof(err))"
    end
    @test occursin("sea_ice_branch_observer", message)
    @test isempty(R._SEA_ICE_BRANCH_OBSERVERS)
    @test R._observe_sea_ice_branch_before_thermodynamics!(nothing, 900.0) === nothing
    @test R._observe_sea_ice_branch_after_thermodynamics!(nothing, 900.0) === nothing
    @test R.SEA_ICE_BRANCH_OBSERVER_SCHEMA == "sea_ice_branch_observer_v1"
    @test length(R.SEA_ICE_BRANCH_OBSERVER_NOT_IN_V1) == 6
end

const NONE = run_arm(analytic(4))
const OFF = run_arm(analytic(4; sea_ice_branch_observer = false))
const ON = run_arm(analytic(4; sea_ice_branch_observer = true))
const ON_SAVED = R.save_sea_ice_branch_observer(ON.simulation, ON.config)

@testset "Sea-ice branch observer: passivity off and on" begin
    @test !NONE.registered && !OFF.registered && ON.registered
    @test isempty(R._SEA_ICE_BRANCH_OBSERVERS)
    @test NONE.simulation.model.clock.iteration == 4
    @test all(k -> haskey(NONE.final, k), ("ice_A", "ice_h", "snow_h", "ice_mass_flux_ice", "ice_mass_flux_snow",
                                          "ice_intercepted_snowfall", "ice_top_surface_temperature", "grid_humidity"))
    @test count(startswith("land_"), keys(NONE.final)) > 0
    @test count(startswith("ice_prognostic_"), keys(NONE.final)) >= 3
    @test count(startswith("ocean_prognostic_"), keys(NONE.final)) >= 4
    # Run-to-run determinism in one process (no observer, and observer explicitly off).
    @test OFF.final == NONE.final
    @test OFF.recorder.hashes == NONE.recorder.hashes
    # Passivity when enabled: every prognostic field at the end, and the fifteen v019 hashes after every step.
    @test ON.final == NONE.final
    @test ON.recorder.iterations == [0, 1, 2, 3, 4]
    @test ON.recorder.hashes == NONE.recorder.hashes
    # The observer adds no callback, so checkpoints have the same structure.
    @test Set(keys(ON.simulation.callbacks)) == Set(keys(NONE.simulation.callbacks))
    @test R._dynamic_config_identity_hash(ON.config) == R._dynamic_config_identity_hash(NONE.config)
end

@testset "Sea-ice branch observer: cadence, states, regions and ledgers (config path)" begin
    s = ON_SAVED
    observer = s.observer
    @test s.pass
    @test s.stages == 4
    @test s.boundaries == 5
    @test observer.stage_iteration == [1, 2, 3, 4]
    @test observer.boundary_iteration == [0, 1, 2, 3, 4]
    @test all(==(900.0), observer.stage_dt)
    @test s.invalid_cells == 0
    @test s.ice_ledger_failed_cells == 0
    @test s.snow_ledger_failed_cells == 0
    @test s.max_ice_ledger_ratio <= 1
    @test s.max_snow_ledger_ratio <= 1
    # The analytic ice has A = 0.8f0, which is above 0.80 in the contract's Float64 comparison, so the contract
    # loose-pack mask is empty on this case. The supplied-mask arm below exercises the per-cell path.
    @test Float64(0.8f0) > 0.80
    @test s.fixed_mask_cells == 0
    @test s.fixed_mask_source == "iteration_0_state"
    # Each boundary record equals the live post-step state seen by the passivity recorder.
    for (k, n) in enumerate(observer.boundary_iteration)
        live = ON.recorder.hashes[findfirst(==(n), ON.recorder.iterations)]
        @test observer.boundary_hashes[k] == [live["ice_A"], live["ice_h"], live["snow_h"]]
    end
    record = read_record(s.path)
    v = record.values
    r(name) = findfirst(==(name), record.regions)
    @test record.regions == ["north_of_50S", "south_of_50S", "fixed_loose_pack",
                             "sector_120E_180E_south_of_50S", "current_loose_pack"]
    @test record.attributes["schema"] == "sea_ice_branch_observer_v1"
    @test record.attributes["validation_pass"] == 1
    @test record.attributes["map_every_n_steps"] == 96
    # Hemisphere sums at the first and last boundary against an independent host recomputation.
    ice = ON.simulation.model.sea_ice.model
    final = (A = surface(ice.ice_concentration), h = surface(ice.ice_thickness), hs = surface(ice.snow_thickness))
    for (state, k) in ((ON.initial, 1), (final, 5))
        A, h, hs = Float64.(state.A), Float64.(state.h), Float64.(state.hs)
        for (mask, name) in ((SOUTH, "south_of_50S"), (NORTH, "north_of_50S"))
            @test isapprox(v["boundary_ice_area"][r(name), k], region_sum(mask, A); rtol = 1e-5)
            @test isapprox(v["boundary_ice_volume"][r(name), k], region_sum(mask, A .* h); rtol = 1e-5)
            @test isapprox(v["boundary_snow_volume"][r(name), k], region_sum(mask, A .* hs); rtol = 1e-5, atol = 1e-3)
            @test isapprox(v["boundary_cell_area"][r(name), k], region_sum(mask, ones(NX, NY)); rtol = 1e-5)
        end
    end
    @test v["boundary_ice_area"][r("south_of_50S"), 1] > 0
    @test v["boundary_ice_area"][r("north_of_50S"), 1] > 0
    @test count(SECTOR) == 15
    @test v["stage_cells"][r("sector_120E_180E_south_of_50S"), 1] == count(SECTOR)
    @test isapprox(v["boundary_cell_area"][r("sector_120E_180E_south_of_50S"), 1],
                   region_sum(SECTOR, ones(NX, NY)); rtol = 1e-5)
    # Ledgers in every region and step.
    @test all(abs.(v["stage_ice_mass_ledger_residual"]) .<= v["stage_ice_mass_ledger_tolerance"])
    @test all(abs.(v["stage_snow_mass_ledger_residual"]) .<= v["stage_snow_mass_ledger_tolerance"])
    hemispheres = [r("north_of_50S"), r("south_of_50S")]
    @test sum(abs, v["stage_net_volume_tendency"][hemispheres, :]) > 0
    # Stage maps at the first step only (one-day map cadence) and boundary maps at the first and last boundary.
    @test v["stage_map_iteration"] == [1]
    @test v["boundary_map_iteration"] == [0, 4]
    @test v["boundary_map_ice_concentration"][:, :, 1] == ON.initial.A
    check_saved_files(s.path)
end

function south_ice_mask(simulation)
    A = surface(simulation.model.sea_ice.model.ice_concentration)
    return BitMatrix([A[i, j] > 0 && latitude_center(j) <= -50 for i in 1:NX, j in 1:NY])
end
const CELLS = run_arm(analytic(4); install = simulation -> R.install_sea_ice_branch_observer!(simulation;
    fixed_loose_pack_mask = south_ice_mask(simulation), fixed_mask_source = "test_south_ice_cells",
    map_every_n_steps = 1))
const CELLS_SAVED = R.save_sea_ice_branch_observer(CELLS.simulation, joinpath(CELLS.config.output_dir, "observer"))

@testset "Sea-ice branch observer: per-cell rows, maps, independent sums and ledgers (supplied mask)" begin
    s = CELLS_SAVED
    observer = s.observer
    @test CELLS.registered
    @test isempty(R._SEA_ICE_BRANCH_OBSERVERS)
    @test CELLS.final == NONE.final
    @test CELLS.recorder.hashes == NONE.recorder.hashes
    @test s.pass
    @test s.fixed_mask_cells == 64
    @test s.fixed_mask_source == "test_south_ice_cells"
    # The contract's representative-cell rule on this mask (all cells wet, rows j = 1 and 2).
    @test [(role, Tuple(I)) for (role, I) in observer.selected_cells] ==
          [("interior_first", (1, 1)), ("interior_last", (32, 2)), ("edge_first", (1, 2))]
    record = read_record(s.path)
    v = record.values
    r(name) = findfirst(==(name), record.regions)
    @test v["stage_map_iteration"] == [1, 2, 3, 4]
    @test v["boundary_map_iteration"] == [0, 1, 2, 3, 4]
    @test v["boundary_map_ice_concentration"][:, :, 1] == CELLS.initial.A
    @test v["boundary_map_ice_thickness"][:, :, 1] == CELLS.initial.h
    fixed = SOUTH .& (CELLS.initial.A .> 0)
    @test count(fixed) == 64
    ci, cj = v["cell_i"], v["cell_j"]
    @test length(ci) == 64
    @test all(j -> j in (1, 2), cj)
    total_change = 0.0
    for k in 1:4
        A0, h0, z0 = (Float64.(v["map_pre_" * n][:, :, k]) for n in ("ice_concentration", "ice_thickness", "snow_thickness"))
        A1, h1, z1 = (Float64.(v["map_post_" * n][:, :, k]) for n in ("ice_concentration", "ice_thickness", "snow_thickness"))
        Fi = Float64.(v["map_ice_thermodynamic_mass_flux"][:, :, k])
        Fs = Float64.(v["map_snow_thermodynamic_mass_flux"][:, :, k])
        P = Float64.(v["map_intercepted_snowfall"][:, :, k])
        dt = v["stage_dt_seconds"][k]
        @test dt == 900.0
        @test all(isfinite, A0) && all(isfinite, A1) && all(isfinite, Fi) && all(isfinite, Fs) && all(isfinite, P)
        @test all(x -> 0 <= x <= 1, A0) && all(x -> 0 <= x <= 1, A1)
        @test all(>=(0), h0) && all(>=(0), h1) && all(>=(0), z0) && all(>=(0), z1)
        # Each stage hash is the hash of the saved map.
        @test observer.stage_hashes[k][1] == digest(v["map_pre_ice_concentration"][:, :, k])
        @test observer.stage_hashes[k][7] == digest(v["map_ice_thermodynamic_mass_flux"][:, :, k])
        # Per-cell rows equal the maps at the mask cells.
        @test [v["map_pre_ice_concentration"][ci[c], cj[c], k] for c in eachindex(ci)] ==
              v["cell_pre_ice_concentration"][:, k]
        @test [v["map_post_ice_thickness"][ci[c], cj[c], k] for c in eachindex(ci)] ==
              v["cell_post_ice_thickness"][:, k]
        @test [v["map_ice_thermodynamic_mass_flux"][ci[c], cj[c], k] for c in eachindex(ci)] ==
              v["cell_ice_thermodynamic_mass_flux"][:, k]
        for (mask, name) in ((SOUTH, "south_of_50S"), (NORTH, "north_of_50S"))
            # Hemisphere sums against the independent analytic-area recomputation.
            @test isapprox(v["stage_pre_ice_area"][r(name), k], region_sum(mask, A0); rtol = 1e-5)
            @test isapprox(v["stage_pre_ice_volume"][r(name), k], region_sum(mask, A0 .* h0); rtol = 1e-5)
            @test isapprox(v["stage_pre_snow_volume"][r(name), k], region_sum(mask, A0 .* z0); rtol = 1e-5, atol = 1e-3)
            @test isapprox(v["stage_post_ice_area"][r(name), k], region_sum(mask, A1); rtol = 1e-5)
            @test isapprox(v["stage_post_ice_volume"][r(name), k], region_sum(mask, A1 .* h1); rtol = 1e-5)
            @test isapprox(v["stage_post_snow_volume"][r(name), k], region_sum(mask, A1 .* z1); rtol = 1e-5, atol = 1e-3)
            positive = region_sum(mask, max.(A1 .* h1 .- A0 .* h0, 0.0)) / dt
            @test isapprox(v["stage_positive_volume_tendency"][r(name), k], positive; rtol = 1e-5, atol = 1e-9)
            positive_area = region_sum(mask, max.(A1 .- A0, 0.0)) / dt
            @test isapprox(v["stage_positive_concentration_tendency"][r(name), k], positive_area; rtol = 1e-5,
                           atol = 1e-9)
            # Ice ledger: sum of F_ice * area * dt equals the change in sum of rho_i * A * h * area, at rounding.
            flux_mass = region_sum(mask, Fi .* dt)
            state_mass = region_sum(mask, 900.0 .* (A1 .* h1 .- A0 .* h0))
            scale = region_sum(mask, 900.0 .* (abs.(A1 .* h1) .+ abs.(A0 .* h0)))
            @test abs(flux_mass - state_mass) <= 8 * eps(Float32) * scale
            # Two-term snow ledger: (F_snow + P_intercepted) * dt equals rho_s * d(A * hs).
            snow_flux_mass = region_sum(mask, (Fs .+ P) .* dt)
            snow_state_mass = region_sum(mask, 330.0 .* (A1 .* z1 .- A0 .* z0))
            snow_scale = region_sum(mask, 330.0 .* (abs.(A1 .* z1) .+ abs.(A0 .* z0)) .+ abs.(P) .* dt)
            @test abs(snow_flux_mass - snow_state_mass) <= 8 * eps(Float32) * snow_scale
            total_change += region_sum(mask, abs.(A1 .* h1 .- A0 .* h0))
        end
        # The fixed-mask region is the supplied mask, fixed for the run.
        @test isapprox(v["stage_pre_ice_volume"][r("fixed_loose_pack"), k], region_sum(fixed, A0 .* h0); rtol = 1e-5)
        @test isapprox(v["stage_post_ice_area"][r("fixed_loose_pack"), k], region_sum(fixed, A1); rtol = 1e-5)
    end
    # The thermodynamics changed the ice, so the ledger is not trivially zero.
    @test total_change > 0
    println("SEA_ICE_BRANCH_OBSERVER_LEDGER_NOTE thermodynamic_abs_volume_change_m3=", total_change)
    check_saved_files(s.path)
end

# Production now selects lead closing at h0 = 0.5 m (LC1), so the base arms above already run it. The contrast arm
# is therefore the released rule, h0 = 0.0 m. (The LEAD_* names are kept; they now mean "the other arm".)
const LEAD_NONE = run_arm(analytic(4; sea_ice_lead_closing_thickness_m = 0.0))
const LEAD_ON = run_arm(analytic(4; sea_ice_lead_closing_thickness_m = 0.0, sea_ice_branch_observer = true))
const LEAD_SAVED = R.save_sea_ice_branch_observer(LEAD_ON.simulation, LEAD_ON.config)

@testset "Sea-ice branch observer: the same observer on the released h0 = 0.0 m arm" begin
    @test PRODUCTION.sea_ice_lead_closing_thickness_m == 0.5
    @test ON.simulation.model.sea_ice.model.ice_thermodynamics.concentration_evolution isa
          R.ClimaSeaIce.SeaIceThermodynamics.LeadClosingEvolution
    evolution = LEAD_ON.simulation.model.sea_ice.model.ice_thermodynamics.concentration_evolution
    @test !(evolution isa R.ClimaSeaIce.SeaIceThermodynamics.LeadClosingEvolution)
    @test !LEAD_NONE.registered && LEAD_ON.registered
    @test LEAD_ON.final == LEAD_NONE.final
    @test LEAD_ON.recorder.hashes == LEAD_NONE.recorder.hashes
    @test LEAD_SAVED.pass
    @test LEAD_SAVED.stages == 4
    # One observer: the same code, hook, settings and schema on both arms. Only the arm identity differs.
    for key in ("observer_code_sha256", "observer_hook_code_sha256", "observer_settings_sha256", "schema",
                "source_commit", "source_tree", "cadence")
        @test LEAD_SAVED.observer.provenance[key] == ON_SAVED.observer.provenance[key]
    end
    @test LEAD_SAVED.observer.provenance["sea_ice_lead_closing_thickness_m"] == 0.0
    @test ON_SAVED.observer.provenance["sea_ice_lead_closing_thickness_m"] == 0.5
    @test LEAD_SAVED.observer.provenance["config_identity_sha256"] !=
          ON_SAVED.observer.provenance["config_identity_sha256"]
    check_saved_files(LEAD_SAVED.path)
    changed = sort!([k for k in keys(NONE.final) if NONE.final[k] != LEAD_NONE.final[k]])
    println("SEA_ICE_BRANCH_OBSERVER_LEAD_ARM_NOTE fields_changed_by_h0_0.0m_vs_production_in_4_steps=", length(changed))
end

# Parse the sealed off-path TOML without the TOML package (not a ReadyESM dependency).
function parse_offpath_reference(path)
    hashes = Dict{String, String}()
    top = Dict{String, String}()
    section = ""
    for raw in eachline(path)
        line = strip(raw)
        isempty(line) && continue
        if startswith(line, "[")
            section = String(line)
            continue
        end
        quoted = match(r"^([A-Za-z0-9_]+) = \"([^\"]*)\"$", line)
        if quoted !== nothing
            target = section == "[hashes]" ? hashes : top
            target[String(quoted.captures[1])] = String(quoted.captures[2])
            continue
        end
        number = match(r"^([A-Za-z0-9_]+) = ([0-9]+)$", line)
        number === nothing || (top[String(number.captures[1])] = String(number.captures[2]))
    end
    return (hashes = hashes, top = top)
end

const PARENT_REFERENCE = get(ENV, "SEA_ICE_BRANCH_OBSERVER_PARENT_REFERENCE", "")
if isempty(PARENT_REFERENCE)
    println("SEA_ICE_BRANCH_OBSERVER_PARENT_CHECK_SKIPPED")
else
    bytes2hex(open(sha256, PARENT_REFERENCE)) == get(ENV, "SEA_ICE_BRANCH_OBSERVER_PARENT_REFERENCE_SHA256", "") ||
        error("the parent reference does not match its declared SHA-256")
    reference = parse_offpath_reference(PARENT_REFERENCE)
    # The same build and run as offpath_lead.jl mode "off": no extra callback.
    parent_off = R.build_dynamic_esm(analytic(24))
    O.run!(parent_off)
    parent_off_hashes = offpath_hashes(parent_off)
    parent_on_config = analytic(24; sea_ice_branch_observer = true)
    parent_on = R.build_dynamic_esm(parent_on_config)
    O.run!(parent_on)
    parent_on_hashes = offpath_hashes(parent_on)
    parent_on_saved = R.save_sea_ice_branch_observer(parent_on, parent_on_config)
    for (label, hashes) in (("off", parent_off_hashes), ("on", parent_on_hashes))
        differ = sort!([k for k in keys(reference.hashes) if get(hashes, k, "") != reference.hashes[k]])
        println("SEA_ICE_BRANCH_OBSERVER_PARENT_COMPARE arm=$label fields=$(length(reference.hashes)) differ=$differ")
    end
    @testset "Sea-ice branch observer: observer off and on equal the parent source over 24 steps" begin
        @test reference.top["mode"] == "off"
        @test reference.top["source_commit"] == "95459f537da2fafb070e3a8a1514df4996284c5b"
        @test reference.top["iterations"] == "24"
        @test length(reference.hashes) == 50
        @test parent_off.model.clock.iteration == 24
        @test parent_on.model.clock.iteration == 24
        @test Set(keys(parent_off_hashes)) == Set(keys(reference.hashes))
        @test parent_off_hashes == reference.hashes
        @test parent_on_hashes == reference.hashes
        @test parent_on_saved.pass
        @test parent_on_saved.stages == 24
        @test parent_on_saved.boundaries == 25
    end
    println("SEA_ICE_BRANCH_OBSERVER_PARENT_EQUIVALENCE_PASS")
end

# v1.1 ledger floor (LC1 pair v3, job 64575899): subnormal snow and ice states whose fluxes are flushed to zero pass;
# normal cells that close at rounding pass; a normal cell broken by 1e-4 of its snow mass still fails.
@testset "Sea-ice branch observer: v1.1 ledger floor for subnormal states, $FT" for FT in (Float32, Float64)
    ρi, ρs, Δτ = 900.0, 330.0, 900.0
    tiny = floatmin(FT)
    # Cells: 1 subnormal snow melting to zero with a flushed snow flux (the pair v3 case), 2 subnormal ice lost with a
    # flushed ice flux, 3 normal and closed, 4 normal with a snow flux broken by 1e-4 of its snow mass.
    A0 = FT[0.4770869, 0.4, 0.5, 0.5]
    A1 = FT[0.4778271, 0.4, 0.52, 0.52]
    h0 = FT[1.0, tiny / 64, 1.2, 1.2]
    h1 = FT[1.0, 0, 1.25, 1.25]
    hs0 = FT[tiny / 64, 0, 0.2, 0.2]
    hs1 = FT[0, 0, 0.19, 0.19]
    exact_ice(i) = ρi * (Float64(A1[i]) * Float64(h1[i]) - Float64(A0[i]) * Float64(h0[i])) / Δτ
    exact_snow(i) = ρs * (Float64(A1[i]) * Float64(hs1[i]) - Float64(A0[i]) * Float64(hs0[i])) / Δτ
    broken = 1e-4 * ρs * (Float64(A1[4]) * Float64(hs1[4]) + Float64(A0[4]) * Float64(hs0[4])) / Δτ
    Fi = FT[exact_ice(1), 0, exact_ice(3), exact_ice(4)]
    Fs = FT[0, 0, exact_snow(3), exact_snow(4) + broken]
    column(v) = reshape(v, 4, 1)
    s = (; A0 = column(A0), h0 = column(h0), hs0 = column(hs0), A1 = column(A1), h1 = column(h1), hs1 = column(hs1),
         Fi = column(Fi), Fs = column(Fs), P = column(zeros(FT, 4)))
    observer = (; weight = ones(4, 1), ice_density = fill(ρi, 4, 1), snow_density = fill(ρs, 4, 1),
                ledger_epsilon = Float64(eps(FT)), float_type = FT)
    # The reduction returns its values in the order of _SEA_ICE_BRANCH_STAGE_QUANTITIES.
    v = R._sea_ice_branch_stage_values(observer, trues(4, 1), s, Δτ)
    q(name) = findfirst(x -> first(x) === name, R._SEA_ICE_BRANCH_STAGE_QUANTITIES)
    @test length(v) == length(R._SEA_ICE_BRANCH_STAGE_QUANTITIES)
    @test v[q(:invalid_cells)] == 0
    @test v[q(:ice_ledger_failed_cells)] == 0
    @test v[q(:snow_ledger_failed_cells)] == 1
    # Without the floor, cells 1 and 2 would fail: their residual exceeds the relative bound alone.
    relative = 8 * eps(FT)
    r_snow_1 = -ρs * (Float64(A1[1]) * Float64(hs1[1]) - Float64(A0[1]) * Float64(hs0[1]))
    @test abs(r_snow_1) > relative * ρs * Float64(A0[1]) * Float64(hs0[1])
    r_ice_2 = -ρi * (Float64(A1[2]) * Float64(h1[2]) - Float64(A0[2]) * Float64(h0[2]))
    @test abs(r_ice_2) > relative * ρi * Float64(A0[2]) * Float64(h0[2])
    # The floor is negligible: at most 8 floatmin(FT) (2 rho + 2 dt) per cell.
    @test 8 * Float64(tiny) * (2ρs + 2Δτ) < 1e-30
end
println("SEA_ICE_BRANCH_OBSERVER_TESTS_PASS")
end
