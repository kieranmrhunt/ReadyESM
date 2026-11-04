# Sea-ice branch observer v1 (schema sea_ice_branch_observer_v1). An opt-in, read-only record of the final-stage sea-ice
# thermodynamic update, for the paired January lead-closing screen. Measurement contract:
# artifacts/development_wave_20261103/v019_ice_branch_trace_contract_codex_v1/README.md (manifest 897467e3).
# Bundle and requirement map: artifacts/development_wave_20261103/sea_ice_branch_observer_claude_v1/README.md.
#
# What it records
# ---------------
# Inside `rk_substep!` (conservative_sea_ice_advection.jl), in the final Runge-Kutta stage only, immediately before
# and immediately after `thermodynamic_time_step!` and before `_consolidate_conditional_ice_state!`:
#   - before: the post-advection state A, h and hs, and the step-start state from the timestepper cache Ψ⁻, which is
#     the state at the previous step boundary;
#   - after: A, h and hs, and the thermodynamic ice and snow mass fluxes and the intercepted snowfall that the kernel
#     has just written.
# On host copies it then forms area sums of A, A*h and A*hs, positive, negative and net increments, the ice ledger
# (F_ice * dt = rho_i * d(A*h)) and the two-term snow ledger ((F_snow + P_intercepted) * dt = rho_s * d(A*hs)) for
# five regions, a SHA-256 hash of every sampled array, per-cell rows for the fixed loose-pack mask, and full maps at
# a declared cadence. `save_sea_ice_branch_observer` adds the final boundary state and writes the record.
#
# It does not record the kernel-local heat and branch terms (see SEA_ICE_BRANCH_OBSERVER_NOT_IN_V1). It is a state,
# mass and area ledger. It cannot attribute area change to a heat source.
#
# Why it is passive
# -----------------
# 1. Off (the default) no observer is installed. Each of the two hook calls in `rk_substep!` tests whether the
#    observer registry is empty and returns `nothing`. It reads no field, launches no kernel, allocates no device
#    memory and does not synchronise. This is the same pattern as the released momentum observer hook.
# 2. On, the hooks act in the final stage only (`clock.stage == Nstages` and `Δτ == clock.last_stage_Δt`). They
#    synchronise the device, which changes timing and not arithmetic. They copy whole interior arrays to new host
#    arrays with `Array(interior(field))`. Every sum, hash and check uses those host copies. The hooks never write a
#    model field or a timestepper field. They never call a kernel, `fill_halo_regions!`, `compute!` or `set!`. They
#    do not change the order or the arguments of the four model calls in `rk_substep!`.
# 3. The observer adds no Oceananigans callback and no model field. The simulation's callbacks, its checkpoint content
#    and its restart identity are the same with and without it.
# 4. `sea_ice_branch_observer` is a run-control field (`_DYNAMIC_RESTART_RUN_CONTROL_FIELDS`), like `scored_outputs`.
#    It does not enter the restart identity hash, so observed and unobserved runs can restart from each other's
#    checkpoints.
# test/sea_ice_branch_observer.jl checks this: observed runs end bit-identical to unobserved runs in every prognostic
# field and at every step, and with the bundle's parent reference they equal the parent source 95459f5.
#
# Cost: one synchronisation and twelve whole-array device-to-host copies per coupled step (not per substep), then
# host reductions and hashes. The observer times itself (`observer_seconds` in the receipt).

const SEA_ICE_BRANCH_OBSERVER_SCHEMA = "sea_ice_branch_observer_v1"
const SEA_ICE_BRANCH_OBSERVER_FIELD_SCHEMA_VERSION = 1
const SEA_ICE_BRANCH_OBSERVER_DIRECTORY = "sea_ice_branch_observer_v1"
const SEA_ICE_BRANCH_PASSIVITY_CALLBACK = :sea_ice_branch_passivity_hashes

# The binding measurement contract and its independent fixed 1 January loose-pack mask (Codex, 3 Nov 2026).
const SEA_ICE_BRANCH_CONTRACT_MANIFEST_SHA256 =
    "897467e3c9906104052079550e81d755ff560b8b40e003f7cc0b41a3d4ca0c11"
const SEA_ICE_BRANCH_CONTRACT_MASK_SHA256 =
    "f4cf0a8bb057e0b0b6e51bddad49a5586deb95aab6a6d0354c5f31d396571b26"
const SEA_ICE_BRANCH_CONTRACT_MASK_CELLS = 1267
const SEA_ICE_BRANCH_CONTRACT_SELECTED_CELLS = ((125, 4), (280, 20), (120, 3))

const _SEA_ICE_BRANCH_LATITUDE_MAX = -50.0
const _SEA_ICE_BRANCH_LOOSE_MIN = 0.15
const _SEA_ICE_BRANCH_LOOSE_MAX = 0.80
const _SEA_ICE_BRANCH_SECTOR = (120.0, 180.0)
# Per-cell ledger tolerance in units of the model's machine epsilon. The kernel's own rounding is at most about
# 2 eps (two products, one difference, one product and one quotient), so 8 eps leaves a factor of four.
const _SEA_ICE_BRANCH_LEDGER_FACTOR = 8.0

const _SEA_ICE_BRANCH_REGIONS = (
    :north_of_50S,
    :south_of_50S,
    :fixed_loose_pack,
    :sector_120E_180E_south_of_50S,
    :current_loose_pack,
)

const _SEA_ICE_BRANCH_REGION_DEFINITIONS = (
    north_of_50S = "wet cells with latitude > -50 degrees. The duplicate eastern half of a RightCenterFolded " *
        "tripolar fold row is excluded from every sum.",
    south_of_50S = "wet cells with latitude <= -50 degrees",
    fixed_loose_pack = "contract mask: wet cells with latitude <= -50 and 0.15 <= A <= 0.80 (Float64 comparison) " *
        "in the iteration-0 state, or a supplied mask. Fixed for the run.",
    sector_120E_180E_south_of_50S = "wet cells with latitude <= -50 and 120 <= (longitude mod 360) < 180, half-open",
    current_loose_pack = "wet cells with latitude <= -50 and 0.15 <= A <= 0.80 in the sampled state. Stage rows " *
        "use the pre-thermodynamic A and boundary rows use the boundary A. The population changes with time.",
)

# Contract items that v1 does not record. They are local to the thermodynamic kernel. Recording them needs separate
# diagnostic writes in the kernel (the v2 design in the bundle README).
const SEA_ICE_BRANCH_OBSERVER_NOT_IN_V1 = (
    "effective top heat Qui and the snow-melt correction Qs*A that enter the volume rate",
    "frazil heat, ice-ocean interface heat and their assembled bottom sum Qbi",
    "alpha, beta, K_m and K_f of the closed-form snow-ice solve",
    "the selected freeze or melt branch and the unclipped cell-net volume rate",
    "the clipped cell-net volume rate apart from snow-ice formation, which v1 sees only as their sum",
    "every wet cell at every step: v1 keeps fixed-mask cells every step and all cells at the map cadence",
)

const SEA_ICE_BRANCH_OBSERVER_INTERPRETATION =
    "v1 measures states, mass fluxes and the A*h and A*hs ledgers at the final-stage thermodynamic update. It can " *
    "show a coupled response and close the mass and snow ledgers. It cannot attribute area change to top, frazil " *
    "or interface heat. A positive January v1 screen is not promoted to July or to more members until a v2 " *
    "kernel trace passes its 96-step passivity and endpoint checks."

const SEA_ICE_BRANCH_OBSERVER_CADENCE =
    "Every coupled step: the final Runge-Kutta stage only (clock.stage == Nstages and Δτ == last_stage_Δt), " *
    "immediately before and after thermodynamic_time_step! and before consolidation. Step boundaries: the " *
    "step-start cache Ψ⁻ at each final stage and the live state at save. Full maps every map_every_n_steps " *
    "steps and at the first and last samples."

const _SEA_ICE_BRANCH_STAGE_ARRAYS = (
    (:pre_ice_concentration, "1"),
    (:pre_ice_thickness, "m"),
    (:pre_snow_thickness, "m"),
    (:post_ice_concentration, "1"),
    (:post_ice_thickness, "m"),
    (:post_snow_thickness, "m"),
    (:ice_thermodynamic_mass_flux, "kg m-2 s-1"),
    (:snow_thermodynamic_mass_flux, "kg m-2 s-1"),
    (:intercepted_snowfall, "kg m-2 s-1"),
)

const _SEA_ICE_BRANCH_BOUNDARY_ARRAYS = (
    (:ice_concentration, "1"),
    (:ice_thickness, "m"),
    (:snow_thickness, "m"),
)

# Per region and final stage. Sums are over wet cells, weighted by cell area. "pre" is after advection and before
# thermodynamics. "post" is after thermodynamics and before consolidation. Mean thicknesses are ice-covered-area
# means (volume over ice area), and 0 where the ice area is 0.
const _SEA_ICE_BRANCH_STAGE_QUANTITIES = (
    (:cells, "1"),
    (:cell_area, "m2"),
    (:pre_ice_area, "m2"),
    (:pre_ice_volume, "m3"),
    (:pre_snow_volume, "m3"),
    (:pre_mean_ice_thickness, "m"),
    (:pre_mean_snow_thickness, "m"),
    (:post_ice_area, "m2"),
    (:post_ice_volume, "m3"),
    (:post_snow_volume, "m3"),
    (:post_mean_ice_thickness, "m"),
    (:post_mean_snow_thickness, "m"),
    (:area_increase, "m2"),
    (:area_decrease, "m2"),
    (:volume_increase, "m3"),
    (:volume_decrease, "m3"),
    (:snow_volume_change, "m3"),
    (:positive_volume_tendency, "m3 s-1"),
    (:positive_concentration_tendency, "m2 s-1"),
    (:net_volume_tendency, "m3 s-1"),
    (:net_concentration_tendency, "m2 s-1"),
    (:positive_volume_tendency_from_ice_flux, "m3 s-1"),
    (:ice_thermodynamic_mass_flux, "kg s-1"),
    (:snow_thermodynamic_mass_flux, "kg s-1"),
    (:intercepted_snowfall, "kg s-1"),
    (:ice_mass_ledger_residual, "kg"),
    (:ice_mass_ledger_tolerance, "kg"),
    (:snow_mass_ledger_residual, "kg"),
    (:snow_mass_ledger_tolerance, "kg"),
    (:ice_ledger_failed_cells, "1"),
    (:snow_ledger_failed_cells, "1"),
    (:snow_rebase_volume, "m3"),
    (:snow_rebase_mean_depth_change, "m"),
    (:snow_column_mean_depth_change, "m"),
    (:volume_gain_cells, "1"),
    (:volume_loss_cells, "1"),
    (:volume_unchanged_cells, "1"),
    (:area_gain_cells, "1"),
    (:area_loss_cells, "1"),
    (:area_gain_with_volume_loss_cells, "1"),
    (:invalid_cells, "1"),
)

const _SEA_ICE_BRANCH_BOUNDARY_QUANTITIES = (
    (:cells, "1"),
    (:cell_area, "m2"),
    (:ice_area, "m2"),
    (:ice_volume, "m3"),
    (:snow_volume, "m3"),
    (:mean_ice_thickness, "m"),
    (:mean_snow_thickness, "m"),
    (:invalid_cells, "1"),
)

# Installed observers, keyed by the live ClimaSeaIce model that `rk_substep!` receives.
const _SEA_ICE_BRANCH_OBSERVERS = IdDict{Any, Any}()

mutable struct SeaIceBranchObserver
    ice_model::Any
    float_type::DataType
    Nx::Int
    Ny::Int
    coupled_timestep_seconds::Float64
    map_every_n_steps::Int
    ledger_epsilon::Float64
    longitude::Matrix{Float64}
    latitude::Matrix{Float64}
    cell_area::Matrix{Float64}
    active::BitMatrix
    unique_storage::BitMatrix
    weight::Matrix{Float64}
    ice_density::Matrix{Float64}
    snow_density::Matrix{Float64}
    region_masks::Dict{Symbol, BitMatrix}
    fixed_mask::Union{Nothing, BitMatrix}
    fixed_mask_source::String
    cell_indices::Vector{CartesianIndex{2}}
    selected_cells::Vector{Pair{String, CartesianIndex{2}}}
    provenance::Dict{String, Any}
    pending::Any
    stage_iteration::Vector{Int}
    stage_time::Vector{Float64}
    stage_dt::Vector{Float64}
    stage_values::Vector{Float64}
    stage_hashes::Vector{Vector{String}}
    cell_values::Dict{Symbol, Any}
    stage_map_iteration::Vector{Int}
    stage_maps::Dict{Symbol, Any}
    boundary_iteration::Vector{Int}
    boundary_time::Vector{Float64}
    boundary_values::Vector{Float64}
    boundary_hashes::Vector{Vector{String}}
    boundary_map_iteration::Vector{Int}
    boundary_maps::Dict{Symbol, Any}
    seconds::Float64
end

#####
##### Hooks called from rk_substep! (conservative_sea_ice_advection.jl)
#####

@inline function _observe_sea_ice_branch_before_thermodynamics!(model, Δτ)
    isempty(_SEA_ICE_BRANCH_OBSERVERS) && return nothing
    observer = get(_SEA_ICE_BRANCH_OBSERVERS, model, nothing)
    isnothing(observer) && return nothing
    return _sea_ice_branch_before!(observer, model, Δτ)
end

@inline function _observe_sea_ice_branch_after_thermodynamics!(model, Δτ)
    isempty(_SEA_ICE_BRANCH_OBSERVERS) && return nothing
    observer = get(_SEA_ICE_BRANCH_OBSERVERS, model, nothing)
    isnothing(observer) && return nothing
    return _sea_ice_branch_after!(observer, model, Δτ)
end

function _sea_ice_branch_final_stage(model, Δτ)
    model.clock.stage == model.timestepper.Nstages || return false
    Float64(Δτ) == Float64(model.clock.last_stage_Δt) || error(
        "sea-ice branch observer received a final-stage step $(Δτ) that differs from the clock's " *
        "$(model.clock.last_stage_Δt)",
    )
    return true
end

function _sea_ice_branch_before!(observer::SeaIceBranchObserver, model, Δτ)
    _sea_ice_branch_final_stage(model, Δτ) || return nothing
    started = time_ns()
    isnothing(observer.pending) || error("sea-ice branch observer: the previous final stage was not closed")
    Oceananigans.Architectures.synchronize(Oceananigans.Architectures.architecture(model.grid))
    iteration = Int(model.clock.iteration)
    time = Float64(model.clock.time)
    # Ψ⁻ is the exact copy of the state at the start of this step, made by cache_current_fields!. It is the state
    # at the previous step boundary.
    base = model.timestepper.Ψ⁻
    boundary = _sea_ice_branch_state(observer, base.ℵ, base.h, base.hs)
    _sea_ice_branch_record_boundary!(observer, iteration, time, boundary)
    pre = _sea_ice_branch_state(observer, model.ice_concentration, model.ice_thickness, model.snow_thickness)
    observer.pending = (iteration = iteration, time = time, Δτ = Float64(Δτ), A = pre.A, h = pre.h, hs = pre.hs)
    observer.seconds += (time_ns() - started) / 1e9
    return nothing
end

function _sea_ice_branch_after!(observer::SeaIceBranchObserver, model, Δτ)
    _sea_ice_branch_final_stage(model, Δτ) || return nothing
    started = time_ns()
    pending = observer.pending
    isnothing(pending) && error("sea-ice branch observer: no pre-thermodynamic sample for this final stage")
    pending.iteration == Int(model.clock.iteration) || error(
        "sea-ice branch observer: the final stage changed iteration between its two samples",
    )
    Oceananigans.Architectures.synchronize(Oceananigans.Architectures.architecture(model.grid))
    Nx, Ny = observer.Nx, observer.Ny
    fluxes = model.mass_fluxes
    sample = (
        A0 = pending.A,
        h0 = pending.h,
        hs0 = pending.hs,
        A1 = _sea_ice_branch_host(model.ice_concentration, Nx, Ny),
        h1 = _sea_ice_branch_host(model.ice_thickness, Nx, Ny),
        hs1 = _sea_ice_branch_host(model.snow_thickness, Nx, Ny),
        Fi = _sea_ice_branch_host(fluxes.thermodynamics.ice, Nx, Ny),
        Fs = _sea_ice_branch_host(fluxes.thermodynamics.snow, Nx, Ny),
        P = _sea_ice_branch_host(fluxes.intercepted_snowfall, Nx, Ny),
    )
    _sea_ice_branch_record_stage!(observer, pending, sample)
    observer.pending = nothing
    observer.seconds += (time_ns() - started) / 1e9
    return nothing
end

#####
##### Host copies, hashes and masks
#####

function _sea_ice_branch_host(field, Nx, Ny)
    values = Array(Oceananigans.interior(field))
    size(values, 1) == Nx && size(values, 2) == Ny && length(values) == Nx * Ny || error(
        "sea-ice branch observer expected a ($Nx, $Ny) surface field, got $(size(values))",
    )
    return reshape(values, Nx, Ny)
end

_sea_ice_branch_state(observer, A, h, hs) = (
    A = _sea_ice_branch_host(A, observer.Nx, observer.Ny),
    h = _sea_ice_branch_host(h, observer.Nx, observer.Ny),
    hs = _sea_ice_branch_host(hs, observer.Nx, observer.Ny),
)

_sea_ice_branch_digest(values) = bytes2hex(SHA.sha256(reinterpret(UInt8, vec(values))))

# The contract's mask convention: native (longitude, latitude) order, one byte per cell.
_sea_ice_branch_mask_sha256(mask) = bytes2hex(SHA.sha256(reinterpret(UInt8, vec(UInt8.(mask)))))

function _sea_ice_branch_file_sha256(path)
    return try
        bytes2hex(SHA.sha256(read(path)))
    catch
        "unavailable"
    end
end

_sea_ice_branch_density(density::Number, Nx, Ny) = fill(Float64(density), Nx, Ny)
_sea_ice_branch_density(density::Oceananigans.Fields.ConstantField, Nx, Ny) =
    fill(Float64(density.constant), Nx, Ny)
_sea_ice_branch_density(density, Nx, Ny) = Float64.(_sea_ice_branch_host(density, Nx, Ny))

_sea_ice_branch_ratio(numerator, denominator) = denominator > 0 ? numerator / denominator : 0.0

_sea_ice_branch_valid_state(a, h, z) =
    isfinite(a) && isfinite(h) && isfinite(z) && 0 <= a <= 1 && h >= 0 && z >= 0

_sea_ice_branch_is_loose(x) = _SEA_ICE_BRANCH_LOOSE_MIN <= Float64(x) <= _SEA_ICE_BRANCH_LOOSE_MAX

# The contract's fixed-mask rule, as in v019_column_initial_mask_prediction_codex_v1 (wet, latitude <= -50,
# 0.15 <= A <= 0.80 compared in Float64).
function _sea_ice_branch_contract_mask(observer, A)
    return observer.active .& (observer.latitude .<= _SEA_ICE_BRANCH_LATITUDE_MAX) .&
           map(_sea_ice_branch_is_loose, A)
end

_sea_ice_branch_loose_pack(observer, A) =
    observer.region_masks[:south_of_50S] .& map(_sea_ice_branch_is_loose, A)

# Representative cells by the contract's deterministic rule (first and last wet-interior cell, first and last mask
# edge cell, duplicates removed in that order), plus the first coastal cell. Longitude wraps and latitude clips.
function _sea_ice_branch_selected_cells(mask::AbstractMatrix{Bool}, active::AbstractMatrix{Bool})
    Nx, Ny = size(mask)
    selected = Pair{String, CartesianIndex{2}}[]
    candidates = findall(mask)
    isempty(candidates) && return selected
    neighbours(I) = begin
        i, j = Tuple(I)
        ((i == 1 ? Nx : i - 1, j), (i == Nx ? 1 : i + 1, j), (i, max(1, j - 1)), (i, min(Ny, j + 1)))
    end
    interior = filter(I -> all(active[ii, jj] for (ii, jj) in neighbours(I)), candidates)
    edge = filter(I -> any(!mask[ii, jj] for (ii, jj) in neighbours(I)), candidates)
    coastal = filter(I -> !all(active[ii, jj] for (ii, jj) in neighbours(I)), candidates)
    roles = Pair{String, CartesianIndex{2}}[]
    if !isempty(interior)
        push!(roles, "interior_first" => first(interior), "interior_last" => last(interior))
    end
    if !isempty(edge)
        push!(roles, "edge_first" => first(edge), "edge_last" => last(edge))
    end
    isempty(coastal) || push!(roles, "coastal_first" => first(coastal))
    for (role, I) in roles
        any(other -> last(other) == I, selected) || push!(selected, role => I)
    end
    return selected
end

function _sea_ice_branch_set_fixed_mask!(observer::SeaIceBranchObserver, mask::AbstractMatrix{Bool}, source)
    size(mask) == (observer.Nx, observer.Ny) || throw(DimensionMismatch(
        "the fixed loose-pack mask has size $(size(mask)), expected $((observer.Nx, observer.Ny))",
    ))
    fixed = BitMatrix(mask)
    observer.fixed_mask = fixed
    observer.fixed_mask_source = String(source)
    observer.region_masks[:fixed_loose_pack] = fixed .& observer.active .& observer.unique_storage
    observer.cell_indices = findall(observer.region_masks[:fixed_loose_pack])
    observer.selected_cells = _sea_ice_branch_selected_cells(fixed, observer.active)
    return observer
end

"""
    set_sea_ice_branch_fixed_mask!(observer, mask; source = "supplied")

Set the fixed loose-pack mask before the first observed step, for a run that does not start at iteration 0 (a
restart). Use the contract's 1 January mask (SHA-256 `SEA_ICE_BRANCH_CONTRACT_MASK_SHA256`).
"""
function set_sea_ice_branch_fixed_mask!(observer::SeaIceBranchObserver, mask::AbstractMatrix{Bool};
                                        source::AbstractString = "supplied")
    isempty(observer.stage_iteration) && isempty(observer.boundary_iteration) || error(
        "the fixed loose-pack mask must be set before the first observed step",
    )
    return _sea_ice_branch_set_fixed_mask!(observer, mask, source)
end

#####
##### Reductions
#####

function _sea_ice_branch_stage_values(observer, selected, s, Δτ)
    weight = observer.weight
    ice_density = observer.ice_density
    snow_density = observer.snow_density
    ε = _SEA_ICE_BRANCH_LEDGER_FACTOR * observer.ledger_epsilon
    cells = 0
    cell_area = 0.0
    pre_ice_area = 0.0
    pre_ice_volume = 0.0
    pre_snow_volume = 0.0
    post_ice_area = 0.0
    post_ice_volume = 0.0
    post_snow_volume = 0.0
    area_increase = 0.0
    area_decrease = 0.0
    volume_increase = 0.0
    volume_decrease = 0.0
    snow_volume_change = 0.0
    flux_positive_volume = 0.0
    ice_flux = 0.0
    snow_flux = 0.0
    intercepted = 0.0
    ice_residual = 0.0
    ice_tolerance = 0.0
    snow_residual = 0.0
    snow_tolerance = 0.0
    ice_failed = 0
    snow_failed = 0
    rebase_volume = 0.0
    column_volume = 0.0
    volume_gain = 0
    volume_loss = 0
    volume_unchanged = 0
    area_gain = 0
    area_loss = 0
    area_gain_volume_loss = 0
    invalid = 0
    @inbounds for I in eachindex(weight)
        selected[I] || continue
        w = weight[I]
        w > 0 || continue
        cells += 1
        cell_area += w
        a0 = Float64(s.A0[I])
        h0 = Float64(s.h0[I])
        z0 = Float64(s.hs0[I])
        a1 = Float64(s.A1[I])
        h1 = Float64(s.h1[I])
        z1 = Float64(s.hs1[I])
        fi = Float64(s.Fi[I])
        fs = Float64(s.Fs[I])
        pint = Float64(s.P[I])
        if !(_sea_ice_branch_valid_state(a0, h0, z0) && _sea_ice_branch_valid_state(a1, h1, z1) &&
             isfinite(fi) && isfinite(fs) && isfinite(pint))
            invalid += 1
            continue
        end
        ρi = ice_density[I]
        ρs = snow_density[I]
        # For Float32 states these Float64 products are exact.
        v0 = a0 * h0
        v1 = a1 * h1
        s0 = a0 * z0
        s1 = a1 * z1
        dA = a1 - a0
        dV = v1 - v0
        dS = s1 - s0
        pre_ice_area += w * a0
        pre_ice_volume += w * v0
        pre_snow_volume += w * s0
        post_ice_area += w * a1
        post_ice_volume += w * v1
        post_snow_volume += w * s1
        area_increase += w * max(dA, 0.0)
        area_decrease += w * min(dA, 0.0)
        volume_increase += w * max(dV, 0.0)
        volume_decrease += w * min(dV, 0.0)
        snow_volume_change += w * dS
        flux_positive_volume += w * max(fi, 0.0) / ρi
        ice_flux += w * fi
        snow_flux += w * fs
        intercepted += w * pint
        # Ice ledger: the kernel writes F_ice = rho_i * (h1*A1 - h0*A0) / dt.
        r_ice = fi * Δτ - ρi * dV
        t_ice = ε * ρi * (abs(v1) + abs(v0))
        ice_residual += w * r_ice
        ice_tolerance += w * t_ice
        abs(r_ice) <= t_ice || (ice_failed += 1)
        # Two-term snow ledger: the kernel writes F_snow = rho_s * (hs1*A1 - hs0*A0) / dt - P_intercepted.
        r_snow = (fs + pint) * Δτ - ρs * dS
        t_snow = ε * (ρs * (abs(s1) + abs(s0)) + abs(pint) * Δτ)
        snow_residual += w * r_snow
        snow_tolerance += w * t_snow
        abs(r_snow) <= t_snow || (snow_failed += 1)
        # Snow rebase at fixed A*hs: new ice area thins the conditional snow depth hs0 to hs0*A0/A1. The rest of the
        # depth change is column snow change (snowfall, melt and snow-ice formation).
        if a1 > 0
            rebased = z0 * a0 / a1
            rebase_volume += w * a1 * (rebased - z0)
            column_volume += w * a1 * (z1 - rebased)
        end
        dV > 0 && (volume_gain += 1)
        dV < 0 && (volume_loss += 1)
        dV == 0 && (volume_unchanged += 1)
        dA > 0 && (area_gain += 1)
        dA < 0 && (area_loss += 1)
        (dA > 0 && dV < 0) && (area_gain_volume_loss += 1)
    end
    rate = Δτ > 0 ? 1 / Δτ : 0.0
    values = (
        cells = Float64(cells),
        cell_area = cell_area,
        pre_ice_area = pre_ice_area,
        pre_ice_volume = pre_ice_volume,
        pre_snow_volume = pre_snow_volume,
        pre_mean_ice_thickness = _sea_ice_branch_ratio(pre_ice_volume, pre_ice_area),
        pre_mean_snow_thickness = _sea_ice_branch_ratio(pre_snow_volume, pre_ice_area),
        post_ice_area = post_ice_area,
        post_ice_volume = post_ice_volume,
        post_snow_volume = post_snow_volume,
        post_mean_ice_thickness = _sea_ice_branch_ratio(post_ice_volume, post_ice_area),
        post_mean_snow_thickness = _sea_ice_branch_ratio(post_snow_volume, post_ice_area),
        area_increase = area_increase,
        area_decrease = area_decrease,
        volume_increase = volume_increase,
        volume_decrease = volume_decrease,
        snow_volume_change = snow_volume_change,
        positive_volume_tendency = volume_increase * rate,
        positive_concentration_tendency = area_increase * rate,
        net_volume_tendency = (volume_increase + volume_decrease) * rate,
        net_concentration_tendency = (area_increase + area_decrease) * rate,
        positive_volume_tendency_from_ice_flux = flux_positive_volume,
        ice_thermodynamic_mass_flux = ice_flux,
        snow_thermodynamic_mass_flux = snow_flux,
        intercepted_snowfall = intercepted,
        ice_mass_ledger_residual = ice_residual,
        ice_mass_ledger_tolerance = ice_tolerance,
        snow_mass_ledger_residual = snow_residual,
        snow_mass_ledger_tolerance = snow_tolerance,
        ice_ledger_failed_cells = Float64(ice_failed),
        snow_ledger_failed_cells = Float64(snow_failed),
        snow_rebase_volume = rebase_volume,
        snow_rebase_mean_depth_change = _sea_ice_branch_ratio(rebase_volume, post_ice_area),
        snow_column_mean_depth_change = _sea_ice_branch_ratio(column_volume, post_ice_area),
        volume_gain_cells = Float64(volume_gain),
        volume_loss_cells = Float64(volume_loss),
        volume_unchanged_cells = Float64(volume_unchanged),
        area_gain_cells = Float64(area_gain),
        area_loss_cells = Float64(area_loss),
        area_gain_with_volume_loss_cells = Float64(area_gain_volume_loss),
        invalid_cells = Float64(invalid),
    )
    return Float64[getproperty(values, name) for (name, _) in _SEA_ICE_BRANCH_STAGE_QUANTITIES]
end

function _sea_ice_branch_boundary_values(observer, selected, state)
    weight = observer.weight
    cells = 0
    cell_area = 0.0
    ice_area = 0.0
    ice_volume = 0.0
    snow_volume = 0.0
    invalid = 0
    @inbounds for I in eachindex(weight)
        selected[I] || continue
        w = weight[I]
        w > 0 || continue
        cells += 1
        cell_area += w
        a = Float64(state.A[I])
        h = Float64(state.h[I])
        z = Float64(state.hs[I])
        if !_sea_ice_branch_valid_state(a, h, z)
            invalid += 1
            continue
        end
        ice_area += w * a
        ice_volume += w * (a * h)
        snow_volume += w * (a * z)
    end
    values = (
        cells = Float64(cells),
        cell_area = cell_area,
        ice_area = ice_area,
        ice_volume = ice_volume,
        snow_volume = snow_volume,
        mean_ice_thickness = _sea_ice_branch_ratio(ice_volume, ice_area),
        mean_snow_thickness = _sea_ice_branch_ratio(snow_volume, ice_area),
        invalid_cells = Float64(invalid),
    )
    return Float64[getproperty(values, name) for (name, _) in _SEA_ICE_BRANCH_BOUNDARY_QUANTITIES]
end

_sea_ice_branch_region_mask(observer, region, current) =
    region === :current_loose_pack ? current : observer.region_masks[region]

function _sea_ice_branch_record_boundary!(observer, iteration, time, state; force_map::Bool = false)
    if !isempty(observer.boundary_iteration)
        iteration == last(observer.boundary_iteration) + 1 || error(
            "sea-ice branch observer: boundary iteration $iteration does not follow " *
            "$(last(observer.boundary_iteration))",
        )
    end
    if isnothing(observer.fixed_mask)
        iteration == 0 || error(
            "sea-ice branch observer: the fixed 1 January loose-pack mask needs the iteration-0 state or a mask " *
            "set with set_sea_ice_branch_fixed_mask!, but the first observed step starts at iteration $iteration",
        )
        _sea_ice_branch_set_fixed_mask!(
            observer,
            _sea_ice_branch_contract_mask(observer, state.A),
            "iteration_0_state",
        )
    end
    current = _sea_ice_branch_loose_pack(observer, state.A)
    for region in _SEA_ICE_BRANCH_REGIONS
        selected = _sea_ice_branch_region_mask(observer, region, current)
        append!(observer.boundary_values, _sea_ice_branch_boundary_values(observer, selected, state))
    end
    push!(observer.boundary_iteration, iteration)
    push!(observer.boundary_time, time)
    push!(observer.boundary_hashes, String[
        _sea_ice_branch_digest(state.A),
        _sea_ice_branch_digest(state.h),
        _sea_ice_branch_digest(state.hs),
    ])
    if force_map || length(observer.boundary_iteration) == 1 || iteration % observer.map_every_n_steps == 0
        push!(observer.boundary_map_iteration, iteration)
        push!(observer.boundary_maps[:ice_concentration], copy(state.A))
        push!(observer.boundary_maps[:ice_thickness], copy(state.h))
        push!(observer.boundary_maps[:snow_thickness], copy(state.hs))
    end
    return nothing
end

function _sea_ice_branch_record_stage!(observer, pending, s)
    completed = pending.iteration + 1
    if !isempty(observer.stage_iteration)
        completed == last(observer.stage_iteration) + 1 || error(
            "sea-ice branch observer: step $completed does not follow $(last(observer.stage_iteration))",
        )
    end
    current = _sea_ice_branch_loose_pack(observer, s.A0)
    for region in _SEA_ICE_BRANCH_REGIONS
        selected = _sea_ice_branch_region_mask(observer, region, current)
        append!(observer.stage_values, _sea_ice_branch_stage_values(observer, selected, s, pending.Δτ))
    end
    arrays = (s.A0, s.h0, s.hs0, s.A1, s.h1, s.hs1, s.Fi, s.Fs, s.P)
    push!(observer.stage_iteration, completed)
    push!(observer.stage_time, pending.time + pending.Δτ)
    push!(observer.stage_dt, pending.Δτ)
    push!(observer.stage_hashes, String[_sea_ice_branch_digest(values) for values in arrays])
    if !isempty(observer.cell_indices)
        for ((name, _), values) in zip(_SEA_ICE_BRANCH_STAGE_ARRAYS, arrays)
            append!(observer.cell_values[name], values[observer.cell_indices])
        end
    end
    if length(observer.stage_iteration) == 1 || completed % observer.map_every_n_steps == 0
        push!(observer.stage_map_iteration, completed)
        for ((name, _), values) in zip(_SEA_ICE_BRANCH_STAGE_ARRAYS, arrays)
            push!(observer.stage_maps[name], copy(values))
        end
    end
    return nothing
end

#####
##### Installation
#####

function _sea_ice_branch_git(arguments...)
    return try
        String(strip(readchomp(Cmd(["git", "-C", PROJECT_ROOT, arguments...]))))
    catch
        "unavailable"
    end
end

function _sea_ice_branch_source_commit()
    return try
        String(_dynamic_source_commit())
    catch
        "unavailable"
    end
end

function _sea_ice_branch_provenance(every, Δt, arm, identity, lead_closing_thickness_m, FT)
    settings = join((
        "schema=$(SEA_ICE_BRANCH_OBSERVER_SCHEMA)",
        "field_schema_version=$(SEA_ICE_BRANCH_OBSERVER_FIELD_SCHEMA_VERSION)",
        "map_every_n_steps=$(every)",
        "latitude_max_deg=$(_SEA_ICE_BRANCH_LATITUDE_MAX)",
        "loose_pack=$(_SEA_ICE_BRANCH_LOOSE_MIN),$(_SEA_ICE_BRANCH_LOOSE_MAX)",
        "sector_east_deg=$(_SEA_ICE_BRANCH_SECTOR[1]),$(_SEA_ICE_BRANCH_SECTOR[2])",
        "ledger_factor=$(_SEA_ICE_BRANCH_LEDGER_FACTOR)",
        "regions=$(join(String.(_SEA_ICE_BRANCH_REGIONS), ","))",
        "stage_quantities=$(join([String(name) for (name, _) in _SEA_ICE_BRANCH_STAGE_QUANTITIES], ","))",
        "boundary_quantities=$(join([String(name) for (name, _) in _SEA_ICE_BRANCH_BOUNDARY_QUANTITIES], ","))",
    ), "\n")
    status = _sea_ice_branch_git("status", "--porcelain")
    return Dict{String, Any}(
        "schema" => SEA_ICE_BRANCH_OBSERVER_SCHEMA,
        "field_schema_version" => SEA_ICE_BRANCH_OBSERVER_FIELD_SCHEMA_VERSION,
        "contract_manifest_sha256" => SEA_ICE_BRANCH_CONTRACT_MANIFEST_SHA256,
        "contract_fixed_mask_sha256" => SEA_ICE_BRANCH_CONTRACT_MASK_SHA256,
        "source_commit" => _sea_ice_branch_source_commit(),
        "source_tree" => _sea_ice_branch_git("rev-parse", "HEAD^{tree}"),
        "source_status" => status == "unavailable" ? status : isempty(status) ? "clean" : "dirty",
        "observer_code_sha256" =>
            _sea_ice_branch_file_sha256(joinpath(PROJECT_ROOT, "src", "sea_ice_branch_observer.jl")),
        "observer_hook_code_sha256" =>
            _sea_ice_branch_file_sha256(joinpath(PROJECT_ROOT, "src", "conservative_sea_ice_advection.jl")),
        "observer_settings_sha256" => bytes2hex(SHA.sha256(codeunits(settings))),
        "arm" => String(arm),
        "config_identity_sha256" => String(identity),
        "sea_ice_lead_closing_thickness_m" => Float64(lead_closing_thickness_m),
        "coupled_timestep_seconds" => Float64(Δt),
        "map_every_n_steps" => every,
        "float_type" => string(FT),
        "cadence" => SEA_ICE_BRANCH_OBSERVER_CADENCE,
        "observer_mutates_live_state" => 0,
    )
end

"""
    install_sea_ice_branch_observer!(simulation; map_every_n_steps, fixed_loose_pack_mask, ...)
    install_sea_ice_branch_observer!(simulation, config)

Attach the read-only sea-ice branch observer v1 to a built coupled simulation. `build_dynamic_esm` calls the second
form when `config.sea_ice_branch_observer` is true. Full maps are kept every `map_every_n_steps` coupled steps
(default one day). Without `fixed_loose_pack_mask` the contract mask is taken from the iteration-0 state.
"""
function install_sea_ice_branch_observer!(
    simulation;
    map_every_n_steps::Union{Nothing, Integer} = nothing,
    fixed_loose_pack_mask::Union{Nothing, AbstractMatrix{Bool}} = nothing,
    fixed_mask_source::AbstractString = "supplied",
    arm::AbstractString = "",
    config_identity_sha256::AbstractString = "",
    lead_closing_thickness_m::Real = 0.0,
)
    earth = simulation.model
    ice = earth.sea_ice.model
    ocean = earth.ocean.model
    ice.advection isa ConservativeSeaIceAdvection || throw(ArgumentError(
        "the sea-ice branch observer needs ReadyESM's conservative sea-ice advection (sea_ice_dynamics = true)",
    ))
    ice.timestepper isa _SplitRK || throw(ArgumentError(
        "the sea-ice branch observer needs the split Runge-Kutta sea-ice timestepper",
    ))
    !isnothing(ice.snow_thickness) && !isnothing(ice.snow_thermodynamics) || throw(ArgumentError(
        "the sea-ice branch observer needs the layered snow and ice thermodynamics",
    ))
    haskey(_SEA_ICE_BRANCH_OBSERVERS, ice) && error("a sea-ice branch observer is already installed on this model")
    Δt = Float64(simulation.Δt)
    isfinite(Δt) && Δt > 0 || throw(ArgumentError("the coupled timestep must be finite and positive"))
    every = isnothing(map_every_n_steps) ? max(1, round(Int, 86_400 / Δt)) : Int(map_every_n_steps)
    every > 0 || throw(ArgumentError("map_every_n_steps must be positive"))

    grid = ocean.grid
    Nx, Ny = size(grid, 1), size(grid, 2)
    FT = eltype(ice.ice_concentration)
    longitude, latitude = _ocean_surface_coordinate_matrices(ocean.tracers.S)
    size(longitude) == size(latitude) == (Nx, Ny) || throw(DimensionMismatch(
        "ocean surface coordinates do not match the ($Nx, $Ny) grid",
    ))
    active = BitMatrix(_ocean_surface_active_mask(grid) .> 0)
    unique_storage = trues(Nx, Ny)
    if Oceananigans.topology(grid, 2) == Oceananigans.RightCenterFolded
        # The eastern half of the fold row duplicates the western half and is excluded from physical sums.
        unique_storage[(Nx ÷ 2 + 1):Nx, Ny] .= false
    end
    cell_area = _ocean_surface_cell_area(grid)
    wet = active .& unique_storage
    weight = ifelse.(wet, cell_area, 0.0)
    all(isfinite, weight) && all(>=(0), weight) || error("sea-ice branch observer: invalid wet cell areas")
    south = wet .& (latitude .<= _SEA_ICE_BRANCH_LATITUDE_MAX)
    north = wet .& (latitude .> _SEA_ICE_BRANCH_LATITUDE_MAX)
    east = mod.(longitude, 360.0)
    sector = south .& (east .>= _SEA_ICE_BRANCH_SECTOR[1]) .& (east .< _SEA_ICE_BRANCH_SECTOR[2])
    region_masks = Dict{Symbol, BitMatrix}(
        :north_of_50S => north,
        :south_of_50S => south,
        :sector_120E_180E_south_of_50S => sector,
    )
    ice_density = _sea_ice_branch_density(ice.sea_ice_density, Nx, Ny)
    snow_density = _sea_ice_branch_density(ice.snow_density, Nx, Ny)
    all(>(0), ice_density) && all(>(0), snow_density) || error("sea-ice branch observer: invalid densities")
    provenance = _sea_ice_branch_provenance(
        every,
        Δt,
        arm,
        config_identity_sha256,
        lead_closing_thickness_m,
        FT,
    )
    observer = SeaIceBranchObserver(
        ice,
        FT,
        Nx,
        Ny,
        Δt,
        every,
        Float64(eps(FT)),
        longitude,
        latitude,
        cell_area,
        active,
        unique_storage,
        weight,
        ice_density,
        snow_density,
        region_masks,
        nothing,
        "",
        CartesianIndex{2}[],
        Pair{String, CartesianIndex{2}}[],
        provenance,
        nothing,
        Int[],
        Float64[],
        Float64[],
        Float64[],
        Vector{String}[],
        Dict{Symbol, Any}(name => FT[] for (name, _) in _SEA_ICE_BRANCH_STAGE_ARRAYS),
        Int[],
        Dict{Symbol, Any}(name => Matrix{FT}[] for (name, _) in _SEA_ICE_BRANCH_STAGE_ARRAYS),
        Int[],
        Float64[],
        Float64[],
        Vector{String}[],
        Int[],
        Dict{Symbol, Any}(name => Matrix{FT}[] for (name, _) in _SEA_ICE_BRANCH_BOUNDARY_ARRAYS),
        0.0,
    )
    isnothing(fixed_loose_pack_mask) ||
        _sea_ice_branch_set_fixed_mask!(observer, fixed_loose_pack_mask, fixed_mask_source)
    _SEA_ICE_BRANCH_OBSERVERS[ice] = observer
    return observer
end

function install_sea_ice_branch_observer!(simulation, config::ExperimentConfig; kwargs...)
    return install_sea_ice_branch_observer!(
        simulation;
        arm = config.name,
        config_identity_sha256 = _dynamic_config_identity_hash(config),
        lead_closing_thickness_m = config.sea_ice_lead_closing_thickness_m,
        kwargs...,
    )
end

"Return the sea-ice branch observer installed on `simulation`, or `nothing`."
sea_ice_branch_observer(simulation) = get(_SEA_ICE_BRANCH_OBSERVERS, simulation.model.sea_ice.model, nothing)

function remove_sea_ice_branch_observer!(observer::SeaIceBranchObserver)
    get(_SEA_ICE_BRANCH_OBSERVERS, observer.ice_model, nothing) === observer &&
        delete!(_SEA_ICE_BRANCH_OBSERVERS, observer.ice_model)
    return nothing
end

#####
##### Summary and files
#####

function _sea_ice_branch_summary(observer)
    nr = length(_SEA_ICE_BRANCH_REGIONS)
    ns = length(observer.stage_iteration)
    nb = length(observer.boundary_iteration)
    stage = reshape(observer.stage_values, length(_SEA_ICE_BRANCH_STAGE_QUANTITIES), nr, ns)
    boundary = reshape(observer.boundary_values, length(_SEA_ICE_BRANCH_BOUNDARY_QUANTITIES), nr, nb)
    q(name) = findfirst(x -> first(x) === name, _SEA_ICE_BRANCH_STAGE_QUANTITIES)
    bq(name) = findfirst(x -> first(x) === name, _SEA_ICE_BRANCH_BOUNDARY_QUANTITIES)
    # Every wet cell is in exactly one hemisphere region, so these two cover each cell once.
    hemispheres = [
        findfirst(==(:north_of_50S), _SEA_ICE_BRANCH_REGIONS),
        findfirst(==(:south_of_50S), _SEA_ICE_BRANCH_REGIONS),
    ]
    invalid = Int(sum(stage[q(:invalid_cells), hemispheres, :])) +
              Int(sum(boundary[bq(:invalid_cells), hemispheres, :]))
    ice_failed = Int(sum(stage[q(:ice_ledger_failed_cells), hemispheres, :]))
    snow_failed = Int(sum(stage[q(:snow_ledger_failed_cells), hemispheres, :]))
    ratio(residual, tolerance) = tolerance > 0 ? abs(residual) / tolerance : (residual == 0 ? 0.0 : Inf)
    ice_ratio = maximum(ratio.(stage[q(:ice_mass_ledger_residual), :, :], stage[q(:ice_mass_ledger_tolerance), :, :]))
    snow_ratio = maximum(ratio.(stage[q(:snow_mass_ledger_residual), :, :], stage[q(:snow_mass_ledger_tolerance), :, :]))
    first_step = first(observer.stage_iteration)
    consecutive = observer.stage_iteration == collect(first_step:(first_step + ns - 1)) &&
                  observer.boundary_iteration == collect((first_step - 1):(first_step + ns - 1))
    fixed_cells = isnothing(observer.fixed_mask) ? 0 : count(observer.fixed_mask)
    fixed_sha = isnothing(observer.fixed_mask) ? "" : _sea_ice_branch_mask_sha256(observer.fixed_mask)
    selected = [Tuple(I) for (_, I) in observer.selected_cells]
    pass = invalid == 0 && ice_failed == 0 && snow_failed == 0 && ice_ratio <= 1 && snow_ratio <= 1 && consecutive
    return (
        pass = pass,
        stages = ns,
        boundaries = nb,
        consecutive = consecutive,
        invalid_cells = invalid,
        ice_ledger_failed_cells = ice_failed,
        snow_ledger_failed_cells = snow_failed,
        max_ice_ledger_ratio = ice_ratio,
        max_snow_ledger_ratio = snow_ratio,
        fixed_mask_cells = fixed_cells,
        fixed_mask_sha256 = fixed_sha,
        fixed_mask_source = observer.fixed_mask_source,
        fixed_mask_matches_contract = fixed_sha == SEA_ICE_BRANCH_CONTRACT_MASK_SHA256,
        selected_cells = selected,
        selected_cells_match_contract = selected == collect(SEA_ICE_BRANCH_CONTRACT_SELECTED_CELLS),
        observer_seconds = observer.seconds,
    )
end

function _sea_ice_branch_geometry_hashes(observer)
    hashes = Pair{String, String}[
        "ocean_longitude" => _sea_ice_branch_digest(observer.longitude),
        "ocean_latitude" => _sea_ice_branch_digest(observer.latitude),
        "ocean_cell_area" => _sea_ice_branch_digest(observer.cell_area),
        "ocean_surface_active_mask" => _sea_ice_branch_mask_sha256(observer.active),
        "unique_storage_mask" => _sea_ice_branch_mask_sha256(observer.unique_storage),
        "sea_ice_density" => _sea_ice_branch_digest(observer.ice_density),
        "snow_density" => _sea_ice_branch_digest(observer.snow_density),
    ]
    for region in _SEA_ICE_BRANCH_REGIONS
        haskey(observer.region_masks, region) || continue
        push!(hashes, "region_mask_" * String(region) => _sea_ice_branch_mask_sha256(observer.region_masks[region]))
    end
    isnothing(observer.fixed_mask) ||
        push!(hashes, "fixed_loose_pack_mask_contract_convention" => _sea_ice_branch_mask_sha256(observer.fixed_mask))
    return hashes
end

# NetCDF attributes: no empty strings (and the observer writes no NaN or Bool attributes).
_sea_ice_branch_attribute(value) = value isa AbstractString && isempty(value) ? "none" : value

function _sea_ice_branch_put!(dataset, name, values, dimensions, units)
    variable = defVar(dataset, name, eltype(values), dimensions)
    variable.attrib["units"] = units
    variable[ntuple(_ -> Colon(), ndims(values))...] = values
    return variable
end

function _sea_ice_branch_put_maps!(dataset, prefix, dimension, maps, arrays, FT)
    for (name, units) in arrays
        variable = defVar(dataset, prefix * String(name), FT, ("ocean_x", "ocean_y", dimension))
        variable.attrib["units"] = units
        for (k, values) in enumerate(maps[name])
            variable[:, :, k] = values
        end
    end
    return nothing
end

function _sea_ice_branch_write_netcdf(path, observer, summary)
    nr = length(_SEA_ICE_BRANCH_REGIONS)
    ns = length(observer.stage_iteration)
    nb = length(observer.boundary_iteration)
    stage = reshape(observer.stage_values, length(_SEA_ICE_BRANCH_STAGE_QUANTITIES), nr, ns)
    boundary = reshape(observer.boundary_values, length(_SEA_ICE_BRANCH_BOUNDARY_QUANTITIES), nr, nb)
    FT = observer.float_type
    NCDataset(path, "c") do dataset
        defDim(dataset, "ocean_x", observer.Nx)
        defDim(dataset, "ocean_y", observer.Ny)
        defDim(dataset, "region", nr)
        defDim(dataset, "stage", ns)
        defDim(dataset, "boundary", nb)
        surface = ("ocean_x", "ocean_y")
        _sea_ice_branch_put!(dataset, "ocean_longitude", observer.longitude, surface, "degrees_east")
        _sea_ice_branch_put!(dataset, "ocean_latitude", observer.latitude, surface, "degrees_north")
        _sea_ice_branch_put!(dataset, "ocean_cell_area", observer.cell_area, surface, "m2")
        _sea_ice_branch_put!(dataset, "ocean_surface_active_mask", Int8.(observer.active), surface, "1")
        _sea_ice_branch_put!(dataset, "unique_storage_mask", Int8.(observer.unique_storage), surface, "1")
        for region in _SEA_ICE_BRANCH_REGIONS
            haskey(observer.region_masks, region) || continue
            _sea_ice_branch_put!(
                dataset,
                "region_mask_" * String(region),
                Int8.(observer.region_masks[region]),
                surface,
                "1",
            )
        end
        _sea_ice_branch_put!(dataset, "stage_iteration", observer.stage_iteration, ("stage",), "1")
        _sea_ice_branch_put!(dataset, "stage_time_seconds", observer.stage_time, ("stage",), "s")
        _sea_ice_branch_put!(dataset, "stage_dt_seconds", observer.stage_dt, ("stage",), "s")
        for (k, (name, units)) in enumerate(_SEA_ICE_BRANCH_STAGE_QUANTITIES)
            _sea_ice_branch_put!(dataset, "stage_" * String(name), stage[k, :, :], ("region", "stage"), units)
        end
        _sea_ice_branch_put!(dataset, "boundary_iteration", observer.boundary_iteration, ("boundary",), "1")
        _sea_ice_branch_put!(dataset, "boundary_time_seconds", observer.boundary_time, ("boundary",), "s")
        for (k, (name, units)) in enumerate(_SEA_ICE_BRANCH_BOUNDARY_QUANTITIES)
            _sea_ice_branch_put!(dataset, "boundary_" * String(name), boundary[k, :, :], ("region", "boundary"), units)
        end
        nc = length(observer.cell_indices)
        if nc > 0
            defDim(dataset, "cell", nc)
            indices = observer.cell_indices
            _sea_ice_branch_put!(dataset, "cell_i", [I[1] for I in indices], ("cell",), "1")
            _sea_ice_branch_put!(dataset, "cell_j", [I[2] for I in indices], ("cell",), "1")
            _sea_ice_branch_put!(dataset, "cell_longitude", observer.longitude[indices], ("cell",), "degrees_east")
            _sea_ice_branch_put!(dataset, "cell_latitude", observer.latitude[indices], ("cell",), "degrees_north")
            _sea_ice_branch_put!(dataset, "cell_area", observer.cell_area[indices], ("cell",), "m2")
            for (name, units) in _SEA_ICE_BRANCH_STAGE_ARRAYS
                _sea_ice_branch_put!(
                    dataset,
                    "cell_" * String(name),
                    reshape(observer.cell_values[name], nc, ns),
                    ("cell", "stage"),
                    units,
                )
            end
        end
        nm = length(observer.stage_map_iteration)
        if nm > 0
            defDim(dataset, "stage_map", nm)
            _sea_ice_branch_put!(dataset, "stage_map_iteration", observer.stage_map_iteration, ("stage_map",), "1")
            _sea_ice_branch_put_maps!(dataset, "map_", "stage_map", observer.stage_maps, _SEA_ICE_BRANCH_STAGE_ARRAYS, FT)
        end
        nbm = length(observer.boundary_map_iteration)
        if nbm > 0
            defDim(dataset, "boundary_map", nbm)
            _sea_ice_branch_put!(
                dataset,
                "boundary_map_iteration",
                observer.boundary_map_iteration,
                ("boundary_map",),
                "1",
            )
            _sea_ice_branch_put_maps!(
                dataset,
                "boundary_map_",
                "boundary_map",
                observer.boundary_maps,
                _SEA_ICE_BRANCH_BOUNDARY_ARRAYS,
                FT,
            )
        end
        for (key, value) in observer.provenance
            dataset.attrib[key] = _sea_ice_branch_attribute(value)
        end
        dataset.attrib["region_names"] = join(String.(_SEA_ICE_BRANCH_REGIONS), ",")
        for region in _SEA_ICE_BRANCH_REGIONS
            dataset.attrib["region_definition_" * String(region)] =
                getproperty(_SEA_ICE_BRANCH_REGION_DEFINITIONS, region)
        end
        dataset.attrib["not_in_v1"] = join(SEA_ICE_BRANCH_OBSERVER_NOT_IN_V1, " | ")
        dataset.attrib["interpretation"] = SEA_ICE_BRANCH_OBSERVER_INTERPRETATION
        dataset.attrib["selected_cells"] = _sea_ice_branch_attribute(
            join(["$(role):$(I[1]),$(I[2])" for (role, I) in observer.selected_cells], ";"),
        )
        dataset.attrib["fixed_mask_source"] = _sea_ice_branch_attribute(observer.fixed_mask_source)
        dataset.attrib["fixed_mask_cells"] = summary.fixed_mask_cells
        dataset.attrib["fixed_mask_sha256"] = _sea_ice_branch_attribute(summary.fixed_mask_sha256)
        dataset.attrib["fixed_mask_matches_contract"] = Int8(summary.fixed_mask_matches_contract)
        dataset.attrib["selected_cells_match_contract"] = Int8(summary.selected_cells_match_contract)
        dataset.attrib["validation_pass"] = Int8(summary.pass)
        dataset.attrib["invalid_cells"] = summary.invalid_cells
        dataset.attrib["ice_ledger_failed_cells"] = summary.ice_ledger_failed_cells
        dataset.attrib["snow_ledger_failed_cells"] = summary.snow_ledger_failed_cells
    end
    return path
end

function _sea_ice_branch_write_hashes(path, observer)
    open(path, "w") do io
        println(io, "kind\titeration\tarray\tsha256")
        for (name, sha) in _sea_ice_branch_geometry_hashes(observer)
            println(io, "geometry\t-1\t", name, "\t", sha)
        end
        for (k, iteration) in enumerate(observer.stage_iteration)
            for ((name, _), sha) in zip(_SEA_ICE_BRANCH_STAGE_ARRAYS, observer.stage_hashes[k])
                println(io, "stage\t", iteration, "\t", name, "\t", sha)
            end
        end
        for (k, iteration) in enumerate(observer.boundary_iteration)
            for ((name, _), sha) in zip(_SEA_ICE_BRANCH_BOUNDARY_ARRAYS, observer.boundary_hashes[k])
                println(io, "boundary\t", iteration, "\t", name, "\t", sha)
            end
        end
    end
    return path
end

function _sea_ice_branch_write_receipt(path, observer, summary)
    entries = Dict{String, String}()
    for (key, value) in observer.provenance
        entries[key] = string(value)
    end
    for key in propertynames(summary)
        entries["result_" * String(key)] = string(getproperty(summary, key))
    end
    entries["not_in_v1"] = join(SEA_ICE_BRANCH_OBSERVER_NOT_IN_V1, " | ")
    entries["interpretation"] = SEA_ICE_BRANCH_OBSERVER_INTERPRETATION
    for region in _SEA_ICE_BRANCH_REGIONS
        entries["region_definition_" * String(region)] = getproperty(_SEA_ICE_BRANCH_REGION_DEFINITIONS, region)
    end
    open(path, "w") do io
        for key in sort!(collect(keys(entries)))
            println(io, key, " = ", entries[key])
        end
    end
    return path
end

"""
    save_sea_ice_branch_observer(simulation, config; detach = true)
    save_sea_ice_branch_observer(simulation, directory; detach = true)

After the run, add the final boundary state and write `sea_ice_branch_observer_v1.nc`, `hashes.tsv`, `receipt.txt`
and `SHA256SUMS` to `directory` (by default `<output_dir>/sea_ice_branch_observer_v1`). Refuses to overwrite. Returns
a summary with `pass`, the counts and the file hashes. With `detach` the observer is removed from the model.
"""
save_sea_ice_branch_observer(simulation, config::ExperimentConfig; kwargs...) = save_sea_ice_branch_observer(
    simulation,
    joinpath(config.output_dir, SEA_ICE_BRANCH_OBSERVER_DIRECTORY);
    kwargs...,
)

function save_sea_ice_branch_observer(simulation, directory::AbstractString; detach::Bool = true)
    ice = simulation.model.sea_ice.model
    observer = get(_SEA_ICE_BRANCH_OBSERVERS, ice, nothing)
    observer isa SeaIceBranchObserver || error("no sea-ice branch observer is installed on this simulation")
    isnothing(observer.pending) || error("sea-ice branch observer: a final stage was not closed")
    isempty(observer.stage_iteration) && error("sea-ice branch observer recorded no final-stage samples")
    iteration = Int(ice.clock.iteration)
    iteration == last(observer.stage_iteration) || error(
        "sea-ice branch observer: the model is at iteration $iteration but the last observed step is " *
        "$(last(observer.stage_iteration))",
    )
    if last(observer.boundary_iteration) != iteration
        Oceananigans.Architectures.synchronize(Oceananigans.Architectures.architecture(ice.grid))
        final = _sea_ice_branch_state(observer, ice.ice_concentration, ice.ice_thickness, ice.snow_thickness)
        _sea_ice_branch_record_boundary!(observer, iteration, Float64(ice.clock.time), final; force_map = true)
    end
    target = abspath(directory)
    ispath(target) && error("sea-ice branch observer output exists: $target")
    mkpath(dirname(target))
    staging = target * ".part"
    ispath(staging) && rm(staging; recursive = true)
    mkpath(staging)
    summary = _sea_ice_branch_summary(observer)
    _sea_ice_branch_write_netcdf(joinpath(staging, "sea_ice_branch_observer_v1.nc"), observer, summary)
    _sea_ice_branch_write_hashes(joinpath(staging, "hashes.tsv"), observer)
    _sea_ice_branch_write_receipt(joinpath(staging, "receipt.txt"), observer, summary)
    names = sort!(["hashes.tsv", "receipt.txt", "sea_ice_branch_observer_v1.nc"])
    sums = [name => _sea_ice_branch_file_sha256(joinpath(staging, name)) for name in names]
    open(joinpath(staging, "SHA256SUMS"), "w") do io
        for (name, sha) in sums
            println(io, sha, "  ", name)
        end
    end
    mv(staging, target)
    detach && remove_sea_ice_branch_observer!(observer)
    return merge(summary, (path = target, files = Dict(sums), observer = observer))
end

#####
##### Passivity wrapper: the v019 fifteen-field prognostic hashes at every step boundary
#####

_sea_ice_branch_host_data(field) = Array(Base.has_offset_axes(field.data) ? parent(field.data) : field.data)

"""
    sea_ice_branch_prognostic_hashes(simulation)

The fifteen prognostic SHA-256 hashes of the v019 column trace (atmosphere temperature, humidity, vorticity,
divergence and pressure, ocean T, S, u, v and free surface, ice thickness, concentration, snow thickness and ice
velocities), with the same array conventions.
"""
function sea_ice_branch_prognostic_hashes(simulation)
    prognostic = simulation.model.atmosphere.variables.prognostic
    result = Dict{String, String}()
    for name in (:temperature, :humidity, :vorticity, :divergence, :pressure)
        result[String(name)] = _sea_ice_branch_digest(_sea_ice_branch_host_data(getproperty(prognostic, name)))
    end
    ocean, ice = simulation.model.ocean.model, simulation.model.sea_ice.model
    for (name, field) in (
        ("ocean_T", ocean.tracers.T),
        ("ocean_S", ocean.tracers.S),
        ("ocean_u", ocean.velocities.u),
        ("ocean_v", ocean.velocities.v),
        ("ocean_eta", ocean.free_surface.displacement),
        ("ice_h", ice.ice_thickness),
        ("ice_A", ice.ice_concentration),
        ("snow_h", ice.snow_thickness),
        ("ice_u", ice.velocities.u),
        ("ice_v", ice.velocities.v),
    )
        result[name] = _sea_ice_branch_digest(Array(Oceananigans.interior(field)))
    end
    return result
end

"A read-only step-boundary callback that keeps `sea_ice_branch_prognostic_hashes` at every coupled step."
mutable struct SeaIceBranchPassivityHashes
    iterations::Vector{Int}
    hashes::Vector{Dict{String, String}}
end

function (recorder::SeaIceBranchPassivityHashes)(simulation)
    iteration = Int(simulation.model.clock.iteration)
    isempty(recorder.iterations) || iteration == last(recorder.iterations) + 1 || error(
        "sea-ice branch passivity hashes: iteration $iteration does not follow $(last(recorder.iterations))",
    )
    push!(recorder.iterations, iteration)
    push!(recorder.hashes, sea_ice_branch_prognostic_hashes(simulation))
    return nothing
end

"""
    install_sea_ice_branch_passivity_hashes!(simulation)

Add the passivity recorder as an every-step callback. Install it the same way in the unobserved control and in
the observed run, then compare the recorders' hashes step by step.
"""
function install_sea_ice_branch_passivity_hashes!(simulation)
    haskey(simulation.callbacks, SEA_ICE_BRANCH_PASSIVITY_CALLBACK) &&
        error("sea-ice branch passivity hashes are already installed")
    recorder = SeaIceBranchPassivityHashes(Int[], Dict{String, String}[])
    Oceananigans.add_callback!(
        simulation,
        recorder,
        Oceananigans.IterationInterval(1);
        name = SEA_ICE_BRANCH_PASSIVITY_CALLBACK,
    )
    return recorder
end

"Write the recorder's hashes as `iteration<TAB>field<TAB>sha256` lines. Refuses to overwrite."
function save_sea_ice_branch_passivity_hashes(recorder::SeaIceBranchPassivityHashes, path::AbstractString)
    target = abspath(path)
    ispath(target) && error("sea-ice branch passivity hash file exists: $target")
    mkpath(dirname(target))
    open(target, "w") do io
        println(io, "iteration\tfield\tsha256")
        for (iteration, hashes) in zip(recorder.iterations, recorder.hashes)
            for field in sort!(collect(keys(hashes)))
                println(io, iteration, "\t", field, "\t", hashes[field])
            end
        end
    end
    return target
end
