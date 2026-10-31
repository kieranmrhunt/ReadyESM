#=
Opt-in dissipative heating (D1): return the kinetic energy (KE) that the momentum sinks remove as heat in the same air.

Selected by `atmosphere_dissipative_heating`:
- `off` (default): the release components are built unchanged, so no code path differs.
- `diagnose_v1`: the wrappers below record each operator's KE loss and add no heat. Prognostic fields equal `off`.
- `kinetic_energy_v1`: as `diagnose_v1`, and the loss is added as a temperature tendency.

Each operator changes the wind by an increment over the leapfrog step 2Δt, applied to the older wind `u_prev`. For a
local operator, the KE that this increment removes, per unit mass and time, is exactly -ū·F with the midpoint wind
ū = u_prev + Δt F, where F is the operator's tendency. The heating is that loss divided by cp. The operators are:
- SpeedyWeather's surface wind stress (lowest layer), wrapped by `DissipativeSurfaceMomentumFlux`;
- the all-level and surface speed-limit drag, wrapped by `DissipativeSpeedLimitDrag`;
- the K-profile momentum mixing, through its `dissipative_heating` flag (pbl_k_profile_host.jl). Mixing also moves KE
  between layers, so its heat follows the non-negative face and implicit parts of its exact loss, not each layer's
  own KE change;
- the vorticity and divergence hyperdiffusion, wrapped by `DissipativeHyperDiffusion`. Its loss is exact per layer
  from the spectral coefficients. Version 1 adds it as a uniform heating in each layer (the l = 0 coefficient of the
  temperature tendency), so the vertical profile and the total are kept but the horizontal placement is not.
The Robert-Williams time filter is diagnosed (`DissipativeHyperDiffusion.filter_change`), never heated.
Start-up: SpeedyWeather's first one or two steps of a fresh run use shorter steps and no filter. The diagnostics and the
heating assume the ordinary 2Δt filtered step, so those steps are approximate. A restart restores the time stepper's
state and does not repeat them.
Design: artifacts/development_wave_20261030/dissipative_heating_v1/DESIGN.md.
=#

const DISSIPATIVE_HEATING_MODES = (:off, :diagnose_v1, :kinetic_energy_v1)
const DISSIPATIVE_HEATING_SCHEME = "midpoint_kinetic_energy_drag_speed_limit_k_profile_layer_uniform_hyperdiffusion_v1"

_dissipative_heating_on(config::ExperimentConfig) = config.atmosphere_dissipative_heating != :off
_dissipative_heating_heats(config::ExperimentConfig) = config.atmosphere_dissipative_heating == :kinetic_energy_v1

# SURFACE WIND STRESS

"""SpeedyWeather's surface momentum flux, plus the KE loss of its lowest-layer tendency (and heat, if `heat`)."""
struct DissipativeSurfaceMomentumFlux{S} <: SpeedyWeather.AbstractSurfaceMomentumFlux
    flux::S
    heat::Bool
end

Adapt.@adapt_structure DissipativeSurfaceMomentumFlux

DissipativeSurfaceMomentumFlux(spectral_grid::SpeedyWeather.SpectralGrid; heat::Bool) =
    DissipativeSurfaceMomentumFlux(SpeedyWeather.SurfaceMomentumFlux(spectral_grid), heat)

"""The plain SpeedyWeather flux inside a momentum-flux component (the component itself if it is not wrapped)."""
_plain_surface_momentum_flux(flux) = flux
_plain_surface_momentum_flux(flux::DissipativeSurfaceMomentumFlux) = flux.flux

SpeedyWeather.initialize!(scheme::DissipativeSurfaceMomentumFlux, model::SpeedyWeather.PrimitiveEquation) =
    SpeedyWeather.initialize!(scheme.flux, model)

SpeedyWeather.variables(scheme::DissipativeSurfaceMomentumFlux) = (
    SpeedyWeather.variables(scheme.flux)...,
    SpeedyWeather.ParameterizationVariable(:surface_drag_dissipation, SpeedyWeather.Grid2D();
        desc = "Kinetic-energy loss rate of the surface wind stress tendency", units = "W/m^2"),
)

Base.@propagate_inbounds function SpeedyWeather.parameterization!(ij, vars, scheme::DissipativeSurfaceMomentumFlux,
                                                                  model)
    surface = model.geometry.nlayers
    Fu = vars.tendencies.grid.u
    Fv = vars.tendencies.grid.v
    Fu0, Fv0 = Fu[ij, surface], Fv[ij, surface]
    SpeedyWeather.surface_wind_stress!(ij, vars, scheme.flux, model)
    du, dv = Fu[ij, surface] - Fu0, Fv[ij, surface] - Fv0
    NF = typeof(du)
    dt = NF(model.time_stepping.Δt_sec)
    loss = -((vars.grid.u_prev[ij, surface] + dt * du) * du + (vars.grid.v_prev[ij, surface] + dt * dv) * dv)
    mass = vars.grid.pressure_prev[ij] * model.geometry.σ_levels_thick[surface] / NF(model.planet.gravity)
    vars.parameterizations.surface_drag_dissipation[ij] = mass * loss
    if scheme.heat
        vars.tendencies.grid.temperature[ij, surface] += loss / NF(model.atmosphere.heat_capacity)
    end
    return nothing
end

# SPEED-LIMIT DRAG

"""
ReadyESM's speed-limit drag, plus the column KE loss of its tendencies (and heat, if `heat`). The drag acts on each
wind component by its sign with the same strength, so in a large step it can add KE to a small component. The rule then
returns that exact change, which can be negative.
"""
struct DissipativeSpeedLimitDrag{D <: CoupledSurfaceSpeedLimitDrag} <: SpeedyWeather.AbstractDrag
    drag::D
    heat::Bool
end

SpeedyWeather.initialize!(scheme::DissipativeSpeedLimitDrag, model::SpeedyWeather.AbstractModel) =
    SpeedyWeather.initialize!(scheme.drag, model)

SpeedyWeather.variables(scheme::DissipativeSpeedLimitDrag) = (
    SpeedyWeather.variables(scheme.drag)...,
    SpeedyWeather.ParameterizationVariable(:speed_limit_drag_dissipation, SpeedyWeather.Grid2D();
        desc = "Kinetic-energy loss rate of the all-level and surface speed-limit drag", units = "W/m^2"),
)

"""Physical (unscaled) speed-limit drag tendency of one wind component, as in SpeedyWeather and atmosphere_drag.jl."""
_speed_limit_tendency(c, speed_limit, u, v) = -c * max(0, sqrt(u^2 + v^2) - speed_limit)^2 * sign(u)

function SpeedyWeather.drag!(vars, scheme::DissipativeSpeedLimitDrag, lf::Integer, model)
    SpeedyWeather.drag!(vars, scheme.drag, lf, model)
    drag = scheme.drag
    u, v = vars.grid.u.data, vars.grid.v.data
    NF = eltype(u)
    dt = NF(model.time_stepping.Δt_sec)
    n = size(u, 2)
    # Physical tendencies of both branches (the stored tendencies are scaled by the radius in the dynamics).
    c, limit = NF(drag.all_levels.drag), NF(drag.all_levels.speed_limit)
    Fu = _speed_limit_tendency.(c, limit, u, v)
    Fv = _speed_limit_tendency.(c, limit, v, u)
    cs, limits = NF(drag.surface_drag), NF(drag.surface_speed_limit)
    us, vs = view(u, :, n), view(v, :, n)
    view(Fu, :, n) .+= _speed_limit_tendency.(cs, limits, us, vs)
    view(Fv, :, n) .+= _speed_limit_tendency.(cs, limits, vs, us)
    loss = .-((vars.grid.u_prev.data .+ dt .* Fu) .* Fu .+ (vars.grid.v_prev.data .+ dt .* Fv) .* Fv)
    thickness = reshape(model.geometry.σ_levels_thick, 1, n)
    mass = vars.grid.pressure_prev.data ./ NF(model.planet.gravity)
    vars.parameterizations.speed_limit_drag_dissipation.data .= mass .* vec(sum(loss .* thickness; dims = 2))
    if scheme.heat
        scale = vars.prognostic.scale[]
        vars.tendencies.grid.temperature.data .+= (scale / NF(model.atmosphere.heat_capacity)) .* loss
    end
    return nothing
end

# HYPERDIFFUSION

"""
SpeedyWeather's hyperdiffusion, plus the exact KE loss of its vorticity and divergence increments in each layer (and
a uniform heating of each layer, if `heat`). `layer_loss` holds the last step's area-mean loss in W/kg per layer.
The uniform placement is a stated version-1 approximation: the loss is largest near jet level in the storm tracks, and
placing it there is a later arm.
`filter_change` holds the KE change rate (W/kg per layer, signed) that the coming Robert-Williams filter makes to the
two time levels it changes, (ΔKE_now + ΔKE_new) / 2Δt. It assumes the filtered leapfrog of every step after the first.
"""
struct DissipativeHyperDiffusion{H, L, W, V} <: SpeedyWeather.AbstractHorizontalDiffusion
    diffusion::H
    heat::Bool
    vorticity_tendency::L      # tendencies before the diffusion (work copies)
    divergence_tendency::L
    weight::W                  # KE per |coefficient|^2: (m > 0 ? 2 : 1) / (l (l + 1)) / (8 pi), zero for l = 0
    layer_loss::V
    filter_change::V
end

function DissipativeHyperDiffusion(spectral_grid::SpeedyWeather.SpectralGrid, diffusion; heat::Bool)
    NF = spectral_grid.NF
    spectrum = spectral_grid.spectrum
    ls, ms = Array(spectrum.l_indices), Array(spectrum.m_indices)      # one-based degree and order
    weight = [l == 1 ? zero(NF) : NF((m == 1 ? 1 : 2) / ((l - 1) * l) / (8 * pi)) for (l, m) in zip(ls, ms)]
    architecture = spectral_grid.architecture
    work() = zeros(Complex{NF}, spectrum, spectral_grid.nlayers)
    layers() = SpeedyWeather.on_architecture(architecture, zeros(NF, spectral_grid.nlayers))
    return DissipativeHyperDiffusion(diffusion, heat, work(), work(),
        SpeedyWeather.on_architecture(architecture, weight), layers(), layers())
end

SpeedyWeather.initialize!(scheme::DissipativeHyperDiffusion, model::SpeedyWeather.AbstractModel) =
    SpeedyWeather.initialize!(scheme.diffusion, model)

"""The plain SpeedyWeather diffusion inside a diffusion component (the component itself if it is not wrapped)."""
_plain_horizontal_diffusion(diffusion) = diffusion
_plain_horizontal_diffusion(diffusion::DissipativeHyperDiffusion) = diffusion.diffusion

"""
Exact KE change rate (W/kg, area mean, per layer) of a leapfrog increment dt (after - before) applied to the older
state `old`, written into `loss` as a positive loss.
"""
function _spectral_kinetic_energy_loss!(loss, weight, old, before, after, dt, seconds)
    increment = sum(weight .* real.(conj.(old .+ (dt / 2) .* (before .+ after)) .* (after .- before)); dims = 1)
    loss .-= (2 * dt / seconds) .* vec(increment)   # KE change = 2 dt increment over the step of `seconds`
    return loss
end

"""
KE change rate of the Robert-Williams filter in SpeedyWeather's `leapfrog_kernel!`: with new = old + dt tendency and
d = old - 2 now + new, the filter sets now + w1 d and new - w2 d. Adds (ΔKE_now + ΔKE_new) / `seconds` to `change`.
"""
function _filter_kinetic_energy_change!(change, weight, old, now, tendency, dt, w1, w2, seconds)
    new = old .+ dt .* tendency
    d = old .- 2 .* now .+ new
    total = sum(weight .* (2 .* w1 .* real.(conj.(now) .* d) .- 2 .* w2 .* real.(conj.(new) .* d) .+
                           (w1^2 + w2^2) .* abs2.(d)); dims = 1)
    change .+= vec(total) ./ seconds
    return change
end

function SpeedyWeather.horizontal_diffusion!(vars::SpeedyWeather.Variables, scheme::DissipativeHyperDiffusion,
                                             model::SpeedyWeather.PrimitiveEquation, lf::Integer = 1)
    vorticity_tendency = vars.tendencies.vorticity
    divergence_tendency = vars.tendencies.divergence
    copyto!(scheme.vorticity_tendency.data, vorticity_tendency.data)
    copyto!(scheme.divergence_tendency.data, divergence_tendency.data)
    SpeedyWeather.horizontal_diffusion!(vars, scheme.diffusion, model, lf)

    stepping = model.time_stepping
    NF = eltype(scheme.layer_loss)
    dt = 2 * NF(stepping.Δt)                 # leapfrog step in the scaled time of the tendencies
    seconds = 2 * NF(stepping.Δt_sec)
    fill!(scheme.layer_loss, 0)
    vorticity = SpeedyWeather.get_step(vars.prognostic.vorticity, lf)
    divergence = SpeedyWeather.get_step(vars.prognostic.divergence, lf)
    _spectral_kinetic_energy_loss!(scheme.layer_loss, scheme.weight, vorticity.data, scheme.vorticity_tendency.data,
                                   vorticity_tendency.data, dt, seconds)
    _spectral_kinetic_energy_loss!(scheme.layer_loss, scheme.weight, divergence.data,
                                   scheme.divergence_tendency.data, divergence_tendency.data, dt, seconds)
    # The tendencies are now final: the leapfrog and its filter come next.
    robert, williams = NF(stepping.robert_filter), NF(stepping.williams_filter)
    w1, w2 = robert * williams / 2, robert * (1 - williams) / 2
    fill!(scheme.filter_change, 0)
    for (variable, tendency) in ((vars.prognostic.vorticity, vorticity_tendency),
                                 (vars.prognostic.divergence, divergence_tendency))
        _filter_kinetic_energy_change!(scheme.filter_change, scheme.weight, SpeedyWeather.get_step(variable, 1).data,
                                       SpeedyWeather.get_step(variable, 2).data, tendency.data, dt, w1, w2, seconds)
    end
    if scheme.heat
        # A uniform field c has l = 0 coefficient c sqrt(4 pi). The tendency is in the scaled time: per second times
        # Δt_sec / Δt (the radius).
        factor = sqrt(4 * NF(pi)) * NF(stepping.Δt_sec) / NF(stepping.Δt) / NF(model.atmosphere.heat_capacity)
        view(vars.tendencies.temperature.data, 1, :) .+= factor .* scheme.layer_loss
    end
    return nothing
end

# MODEL ASSEMBLY, PROVENANCE AND RESTART

"""Wrap `diffusion`, `drag` and the surface momentum flux when dissipative heating is selected (unchanged if off)."""
_dissipative_horizontal_diffusion(config::ExperimentConfig, spectral_grid, diffusion) =
    _dissipative_heating_on(config) ?
    DissipativeHyperDiffusion(spectral_grid, diffusion; heat = _dissipative_heating_heats(config)) : diffusion

_dissipative_speed_limit_drag(config::ExperimentConfig, drag) =
    _dissipative_heating_on(config) ? DissipativeSpeedLimitDrag(drag, _dissipative_heating_heats(config)) : drag

"""Keyword arguments for the model constructor: none if off, so the release default is used unchanged."""
_dissipative_surface_momentum_flux(config::ExperimentConfig, spectral_grid) =
    _dissipative_heating_on(config) ?
    (; surface_momentum_flux = DissipativeSurfaceMomentumFlux(spectral_grid;
                                                              heat = _dissipative_heating_heats(config))) : (;)

"""The selected mode of an assembled atmosphere model, from its components."""
function _dissipative_heating_mode(model)
    diffusion = model.horizontal_diffusion
    diffusion isa DissipativeHyperDiffusion || return :off
    return diffusion.heat ? :kinetic_energy_v1 : :diagnose_v1
end

"""Check that all wrappers and the K-profile flag agree with the selected mode."""
function _check_dissipative_heating_host(model, config::ExperimentConfig)
    mode = config.atmosphere_dissipative_heating
    heats = mode == :kinetic_energy_v1
    on = mode != :off
    flux, drag, diffusion = model.surface_momentum_flux, model.drag, model.horizontal_diffusion
    (flux isa DissipativeSurfaceMomentumFlux) == on && (drag isa DissipativeSpeedLimitDrag) == on &&
        (diffusion isa DissipativeHyperDiffusion) == on ||
        error("dissipative heating components do not match atmosphere_dissipative_heating=$mode")
    on && !(flux.heat == drag.heat == diffusion.heat == heats) &&
        error("dissipative heating wrappers disagree on heating")
    pbl = model.vertical_diffusion
    pbl isa KProfileBoundaryLayerMixing && pbl.dissipative_heating != (heats && pbl.mix_momentum) &&
        error("K-profile dissipative heating flag does not match atmosphere_dissipative_heating=$mode")
    return nothing
end

"""Diagnostics provenance, written only when the option is not off."""
function _dissipative_heating_provenance(config::ExperimentConfig)
    _dissipative_heating_on(config) || return (;)
    return (atmosphere_dissipative_heating = String(config.atmosphere_dissipative_heating),
            atmosphere_dissipative_heating_scheme = DISSIPATIVE_HEATING_SCHEME)
end

"""The restart record: present only when the option is not off."""
function _dissipative_heating_restart_state(simulation)
    mode = _dissipative_heating_mode(simulation.model)
    return mode == :off ? (;) : (; atmosphere_dissipative_heating = String(mode))
end

"""
Reject a checkpoint whose dissipative-heating mode differs from the reconstructed atmosphere, before any state is
changed. Checkpoints without the record (including every checkpoint written before this option existed) imply off.
"""
function _validate_dissipative_heating_restart(simulation, state)
    saved = hasproperty(state, :atmosphere_dissipative_heating) ? state.atmosphere_dissipative_heating : "off"
    saved isa AbstractString || throw(ArgumentError("checkpoint dissipative-heating mode must be a string"))
    active = String(_dissipative_heating_mode(simulation.model))
    saved == active || throw(ArgumentError(
        "checkpoint and reconstructed atmosphere use different dissipative heating " *
        "(checkpoint $saved, configuration $active)",
    ))
    return nothing
end

# NATIVE-STEP BUDGET

"""
Running totals of the dissipative-heating diagnostics, added only when the option is not off. After every step each
global-mean loss rate is multiplied by the step length Δt and added, so total / `seconds` is the time mean since
launch. Checkpoints carry the totals, so the mean between two checkpoints is their difference over the elapsed time.
Units: J/m2 for the column terms (positive = KE removed, except `filter`, where positive = KE added), J/kg for the
per-layer hyperdiffusion and filter terms (area means). The hyperdiffusion and filter columns use the global-mean
surface pressure of each step.
"""
Base.@kwdef mutable struct DissipativeHeatingBudgetCallback{W} <: SpeedyWeather.AbstractCallback
    point_weights::W
    steps::Int = 0
    seconds::Float64 = 0.0         # Float64 defaults: the generated outer constructor does not convert
    surface_drag::Float64 = 0.0
    speed_limit_drag::Float64 = 0.0
    pbl_face::Float64 = 0.0
    pbl_implicit::Float64 = 0.0
    hyperdiffusion::Float64 = 0.0
    filter::Float64 = 0.0
    hyperdiffusion_layers::Vector{Float64} = Float64[]
    filter_layers::Vector{Float64} = Float64[]
end

const _DISSIPATIVE_HEATING_BUDGET_SCALARS = (
    :steps, :seconds, :surface_drag, :speed_limit_drag, :pbl_face, :pbl_implicit, :hyperdiffusion, :filter,
)

DissipativeHeatingBudgetCallback(spectral_grid::SpectralGrid) = DissipativeHeatingBudgetCallback(;
    point_weights = _global_point_weights(spectral_grid),
    hyperdiffusion_layers = zeros(spectral_grid.nlayers), filter_layers = zeros(spectral_grid.nlayers))

"""The budget callback entry for the model's callbacks: none if off."""
_dissipative_heating_callbacks(config::ExperimentConfig, spectral_grid) =
    _dissipative_heating_on(config) ?
    (:dissipative_heating_budget => DissipativeHeatingBudgetCallback(spectral_grid),) : ()

function SpeedyWeather.initialize!(callback::DissipativeHeatingBudgetCallback, vars::SpeedyWeather.Variables,
                                   model::SpeedyWeather.AbstractModel)
    for name in _DISSIPATIVE_HEATING_BUDGET_SCALARS
        setfield!(callback, name, zero(fieldtype(typeof(callback), name)))
    end
    fill!(callback.hyperdiffusion_layers, 0)
    fill!(callback.filter_layers, 0)
    return nothing
end

function SpeedyWeather.callback!(callback::DissipativeHeatingBudgetCallback, vars::SpeedyWeather.Variables,
                                 model::SpeedyWeather.AbstractModel)
    weights = callback.point_weights
    p = vars.parameterizations
    dt = Float64(model.time_stepping.Δt_sec)
    global_mean(field) = Float64(_global_grid_mean(field.data, weights))
    callback.steps += 1
    callback.seconds += dt
    callback.surface_drag += dt * global_mean(p.surface_drag_dissipation)
    callback.speed_limit_drag += dt * global_mean(p.speed_limit_drag_dissipation)
    if hasproperty(p, :pbl_momentum_dissipation_face)
        callback.pbl_face += dt * global_mean(p.pbl_momentum_dissipation_face)
        callback.pbl_implicit += dt * global_mean(p.pbl_momentum_dissipation_implicit)
    end
    diffusion = model.horizontal_diffusion
    loss, filter = Float64.(Array(diffusion.layer_loss)), Float64.(Array(diffusion.filter_change))
    layer_mass = global_mean(vars.grid.pressure_prev) / Float64(model.planet.gravity) .*
                 Float64.(Array(model.geometry.σ_levels_thick))
    callback.hyperdiffusion += dt * sum(loss .* layer_mass)
    callback.filter += dt * sum(filter .* layer_mass)
    callback.hyperdiffusion_layers .+= dt .* loss
    callback.filter_layers .+= dt .* filter
    return nothing
end

SpeedyWeather.finalize!(::DissipativeHeatingBudgetCallback, args...) = nothing

"""Time means since launch: W/m2 for the column terms, W/kg for the per-layer terms."""
function dissipative_heating_budget_means(callback::DissipativeHeatingBudgetCallback)
    seconds = callback.seconds
    seconds > 0 || throw(ArgumentError("the dissipative-heating budget has no completed step"))
    columns = (; (name => getfield(callback, name) / seconds
                  for name in _DISSIPATIVE_HEATING_BUDGET_SCALARS if name ∉ (:steps, :seconds))...)
    return (; steps = callback.steps, seconds, columns...,
            hyperdiffusion_layers = callback.hyperdiffusion_layers ./ seconds,
            filter_layers = callback.filter_layers ./ seconds)
end

_dissipative_heating_budget_state(callback::DissipativeHeatingBudgetCallback) =
    (; (name => getfield(callback, name) for name in _DISSIPATIVE_HEATING_BUDGET_SCALARS)...,
       hyperdiffusion_layers = copy(callback.hyperdiffusion_layers), filter_layers = copy(callback.filter_layers))

function _restore_dissipative_heating_budget!(callback::DissipativeHeatingBudgetCallback, state)
    for name in _DISSIPATIVE_HEATING_BUDGET_SCALARS
        setfield!(callback, name, getproperty(state, name))
    end
    length(state.hyperdiffusion_layers) == length(callback.hyperdiffusion_layers) ||
        throw(ArgumentError("checkpoint dissipative-heating budget has a different number of layers"))
    copyto!(callback.hyperdiffusion_layers, state.hyperdiffusion_layers)
    copyto!(callback.filter_layers, state.filter_layers)
    return callback
end
