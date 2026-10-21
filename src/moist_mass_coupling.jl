# Optional moist-mass coupling (C10) for the Betts-Miller host.
#
# SpeedyWeather removes or adds water vapour in the physics without changing the air mass: the surface pressure
# does not feel evaporation or precipitation, and the dynamical core deletes the global mean of the ln(p_s)
# tendency. With `atmosphere_moist_mass_coupling = true`, the vapour sources also change the air mass.
#
# 1. `MoistMassCoupling` occupies `stochastic_physics`, the last column parameterization. In each column it runs
#    after the large-scale condensation, the convection and the surface humidity flux. Sources in kg m-2 s-1:
#      Σ   = Σ_k m_k (accumulated physics dq/dt)_k   the column vapour source, E - P (internal transport sums to 0)
#      L_k = m_k (large-scale condensation humidity tendency)_k, a per-layer phase change
#      E   = surface humidity flux, in the bottom layer
#      C   = Σ - Σ_k L_k - E, the Betts-Miller net source. Betts-Miller has no condensation profile, so C is
#            spread in sigma (S_k gets Δσ_k C). This is a stated approximation.
#      S_k = L_k + Δσ_k C + E δ_kK, so Σ_k S_k = Σ up to Float32 layer thicknesses.
#    A passive diagnostic records the moist enthalpy the sources carry, Σ_k S_k (cp T_k + Lv), in W m-2.
#    The tendencies follow the column mass-coupling reference: dilution -q_k S_k, upwind interface terms for q, T, u
#    and v from the interface mass flux F_{k+1/2} = Σ_{j<=k} S_j - σ_{k+1/2} Σ, evaporated vapour at rest in the
#    bottom layer, and tendencies.grid.pressure -= radius g Σ / p_s. The ln(p_s) tendency is not radius-scaled by
#    scale_tendencies!, while the leapfrog uses Δt / radius, so the radius factor is explicit.
# 2. Two SpeedyWeather 0.21.1 methods set the l = 0 ln(p_s) tendency to zero for mass conservation. For a model
#    whose `stochastic_physics` is a `MoistMassCoupling`, `dynamics_tendencies!` and `implicit_correction!` are
#    specialized. Each is an exact copy of the upstream method with only that reset removed (the first through
#    `_moist_mass_surface_pressure_tendency!`). Every other model type keeps the upstream methods.
#    test/moist_mass_coupling.jl checks the copies against the loaded upstream source and the method selection.

"""Column vapour-to-air-mass coupling for the Betts-Miller host. Occupies the `stochastic_physics` slot."""
struct MoistMassCoupling end

const MOIST_MASS_COUPLING_SCHEME = "betts_miller_sigma_spread_remainder_v1"

SpeedyWeather.variables(::MoistMassCoupling) = (
    SpeedyWeather.ParameterizationVariable(:mmc_column_source, SpeedyWeather.Grid2D();
        desc = "Column vapour source E - P coupled to the air mass", units = "kg/m2/s"),
    SpeedyWeather.ParameterizationVariable(:mmc_large_scale_source, SpeedyWeather.Grid2D();
        desc = "Large-scale condensation column vapour source", units = "kg/m2/s"),
    SpeedyWeather.ParameterizationVariable(:mmc_surface_source, SpeedyWeather.Grid2D();
        desc = "Surface humidity flux placed in the bottom layer", units = "kg/m2/s"),
    SpeedyWeather.ParameterizationVariable(:mmc_remainder_source, SpeedyWeather.Grid2D();
        desc = "Convective remainder source, spread in sigma", units = "kg/m2/s"),
    SpeedyWeather.ParameterizationVariable(:mmc_source_enthalpy, SpeedyWeather.Grid2D();
        desc = "Moist enthalpy carried by the vapour sources, sum S_k (cp T_k + Lv)", units = "W/m2"),
)
SpeedyWeather.initialize!(::MoistMassCoupling, ::SpeedyWeather.AbstractModel) = nothing

@inline _mmc_upwind(F, above, below) = F > 0 ? above : below

Base.@propagate_inbounds function SpeedyWeather.parameterization!(ij, vars, ::MoistMassCoupling, model)
    p = vars.parameterizations
    Δσ = model.geometry.σ_levels_thick
    K = length(Δσ)
    g = Float64(model.planet.gravity)
    a = Float64(model.planet.radius)
    ps = Float64(vars.grid.pressure_prev[ij])                  # Pa during the physics
    tq = vars.tendencies.grid.humidity
    lsc = p.large_scale_condensation_humidity_tendency
    E = Float64(p.surface_humidity_flux[ij])
    column = ps / g
    Σ = 0.0
    L = 0.0
    for k in 1:K
        m = Float64(Δσ[k]) * column
        Σ += m * Float64(tq[ij, k])
        L += m * Float64(lsc[ij, k])
    end
    C = Σ - L - E
    p.mmc_column_source[ij] = Σ
    p.mmc_large_scale_source[ij] = L
    p.mmc_surface_source[ij] = E
    p.mmc_remainder_source[ij] = C
    q, T = vars.grid.humidity_prev, vars.grid.temperature_prev
    # Passive energy diagnostic: the moist enthalpy (cp T + Lv per unit vapour mass) that the sources carry when they
    # enter or leave at the layer temperature.
    cp = Float64(model.atmosphere.heat_capacity)
    Lv = Float64(model.atmosphere.latent_heat_condensation)
    enthalpy = 0.0
    for k in 1:K
        S = Float64(Δσ[k]) * column * Float64(lsc[ij, k]) + Float64(Δσ[k]) * C + (k == K ? E : 0.0)
        enthalpy += S * (cp * Float64(T[ij, k]) + Lv)
    end
    p.mmc_source_enthalpy[ij] = enthalpy
    u, v = vars.grid.u_prev, vars.grid.v_prev
    tT, tu, tv = vars.tendencies.grid.temperature, vars.tendencies.grid.u, vars.tendencies.grid.v
    cumulative_S = 0.0
    cumulative_σ = 0.0
    Fabove = 0.0                                                # F_{k-1/2}, positive downward, F_{1/2} = 0
    for k in 1:K
        m = Float64(Δσ[k]) * column
        S = m * Float64(lsc[ij, k]) + Float64(Δσ[k]) * C + (k == K ? E : 0.0)
        cumulative_S += S
        cumulative_σ += Float64(Δσ[k])
        Fbelow = k < K ? cumulative_S - cumulative_σ * Σ : 0.0  # F_{k+1/2}, F_{K+1/2} = 0
        qk, Tk, uk, vk = Float64(q[ij, k]), Float64(T[ij, k]), Float64(u[ij, k]), Float64(v[ij, k])
        rq = -qk * S
        rT = 0.0
        ru = 0.0
        rv = 0.0
        if k > 1
            rq += Fabove * (_mmc_upwind(Fabove, Float64(q[ij, k-1]), qk) - qk)
            rT += Fabove * (_mmc_upwind(Fabove, Float64(T[ij, k-1]), Tk) - Tk)
            ru += Fabove * (_mmc_upwind(Fabove, Float64(u[ij, k-1]), uk) - uk)
            rv += Fabove * (_mmc_upwind(Fabove, Float64(v[ij, k-1]), vk) - vk)
        end
        if k < K
            rq -= Fbelow * (_mmc_upwind(Fbelow, qk, Float64(q[ij, k+1])) - qk)
            rT -= Fbelow * (_mmc_upwind(Fbelow, Tk, Float64(T[ij, k+1])) - Tk)
            ru -= Fbelow * (_mmc_upwind(Fbelow, uk, Float64(u[ij, k+1])) - uk)
            rv -= Fbelow * (_mmc_upwind(Fbelow, vk, Float64(v[ij, k+1])) - vk)
        end
        if k == K && E > 0                                       # evaporated vapour enters at rest
            ru -= E * uk
            rv -= E * vk
        end
        tq[ij, k] += rq / m
        tT[ij, k] += rT / m
        tu[ij, k] += ru / m
        tv[ij, k] += rv / m
        Fabove = Fbelow
    end
    vars.tendencies.grid.pressure[ij] -= a * g * Σ / ps
    return nothing
end

"""`PrimitiveWetModel` types whose `stochastic_physics` is a `MoistMassCoupling`. Built from the loaded
SpeedyWeather type so that the type parameter of `stochastic_physics` is found by field, not by position."""
function _moist_mass_coupled_model_type()
    body = Base.unwrap_unionall(SpeedyWeather.PrimitiveWetModel)
    tv = fieldtype(body, :stochastic_physics)
    tv isa TypeVar || error("PrimitiveWetModel.stochastic_physics is not a type parameter")
    MoistMassCoupling <: tv.ub || error("MoistMassCoupling is outside the stochastic_physics bound $(tv.ub)")
    params = collect(body.parameters)
    i = findfirst(x -> x === tv, params)
    isnothing(i) && error("stochastic_physics type parameter not found")
    restricted = TypeVar(tv.name, tv.lb, MoistMassCoupling)
    params[i] = restricted
    T = SpeedyWeather.PrimitiveWetModel{params...}
    for x in reverse(params)
        T = UnionAll(x, T)
    end
    return T
end

const MoistMassCoupledModel = _moist_mass_coupled_model_type()

# BEGIN MOIST MASS COPY: surface_pressure_tendency!
# SpeedyWeather 0.21.1 dynamics/tendencies.jl surface_pressure_tendency!, with the l = m = 0 reset removed.
function _moist_mass_surface_pressure_tendency!(
        vars::SpeedyWeather.Variables,
        S::SpeedyWeather.AbstractSpectralTransform,
    )
    pres_tend = vars.tendencies.pressure
    pres_tend_grid = vars.tendencies.grid.pressure
    (; dpres_dx, dpres_dy, u_mean_grid, v_mean_grid, div_mean) = vars.dynamics
    scratch_memory = vars.scratch.transform_memory

    # in grid-point space the the (ū, v̄)⋅∇lnpₛ term (swap sign in spectral)
    # += to allow for forcing contributions already in pres_tend_grid
    @. pres_tend_grid += u_mean_grid * dpres_dx + v_mean_grid * dpres_dy

    ūv̄∇lnpₛ = vars.scratch.a_2D             # reuse 2D work array
    SpeedyWeather.transform!(ūv̄∇lnpₛ, pres_tend_grid, scratch_memory, S)

    # for semi-implicit div_mean is calc at time step i-1 in vertical_integration!
    @. pres_tend -= ūv̄∇lnpₛ + div_mean      # add the -div_mean term in spectral, swap sign

    return nothing
end
# END MOIST MASS COPY: surface_pressure_tendency!

# BEGIN MOIST MASS COPY: dynamics_tendencies!
function SpeedyWeather.dynamics_tendencies!(
        vars::SpeedyWeather.Variables,
        lf::Integer,                # leapfrog index for tendencies
        model::MoistMassCoupledModel,
    )
    SpeedyWeather.forcing!(vars, lf, model)
    SpeedyWeather.drag!(vars, lf, model)

    (; orography, geometry, spectral_transform, geopotential, atmosphere, implicit) = model

    # for semi-implicit corrections (α >= 0.5) linear gravity-wave related tendencies are
    # evaluated at previous timestep i-1 (i.e. lf=1 leapfrog time step)
    # nonlinear terms and parameterizations are always evaluated at lf
    lf_implicit = implicit.α == 0 ? lf : 1

    # calculate ∇ln(pₛ), then (u_k, v_k)⋅∇ln(p_s)
    SpeedyWeather.pressure_gradient_flux!(vars, lf, spectral_transform)

    # calculate Tᵥ = T + Tₖμq in spectral as a approxmation to Tᵥ = T(1+μq) used for geopotential
    SpeedyWeather.linear_virtual_temperature!(vars, lf_implicit, model)

    # temperature relative to profile
    # TODO: broadcast with LTA doesn't work here becasue of a broadcast conflict (temp profile and temp_grid are different dimensions and array types)
    vars.grid.temperature.data .-= implicit.temp_profile'

    # from ∂Φ/∂ln(pₛ) = -RTᵥ for bernoulli_potential!
    SpeedyWeather.geopotential!(vars, geopotential, orography)

    # get ū, v̄, D̄ on grid; D̄ in spectral
    SpeedyWeather.vertical_integration!(vars, lf_implicit, geometry)

    # ∂ln(pₛ)/∂t = -(ū, v̄)⋅∇ln(pₛ) - D̄
    _moist_mass_surface_pressure_tendency!(vars, spectral_transform)

    # calculate vertical velocity σ̇ in sigma coordinates for the vertical mass flux M = pₛ * σ̇
    SpeedyWeather.vertical_velocity!(vars, geometry)

    # add the RTₖlnpₛ term to geopotential
    SpeedyWeather.linear_pressure_gradient!(vars, lf_implicit, atmosphere, implicit)

    # use σ̇ for the vertical advection of u, v, T, q
    SpeedyWeather.vertical_advection!(vars, model)

    # vorticity advection, pressure gradient term
    SpeedyWeather.vordiv_tendencies!(vars, model)

    # hor. advection + adiabatic term
    SpeedyWeather.temperature_tendency!(vars, model)

    # horizontal advection of humidity (nothing for wetcore)
    SpeedyWeather.humidity_tendency!(vars, model)

    # add -∇²(E + ϕ + RTₖlnpₛ) term to div tendency
    SpeedyWeather.bernoulli_potential!(vars, spectral_transform)

    # advect all tracers
    SpeedyWeather.tracer_advection!(vars, model)

    # back to absolute temperature
    vars.grid.temperature.data .+= implicit.temp_profile'

    return nothing
end
# END MOIST MASS COPY: dynamics_tendencies!

# BEGIN MOIST MASS COPY: implicit_correction!
function SpeedyWeather.implicit_correction!(
        vars::SpeedyWeather.Variables,
        implicit::SpeedyWeather.ImplicitPrimitiveEquation,
        model::MoistMassCoupledModel,
    )

    # escape immediately if explicit
    implicit.α == 0 && return nothing

    (; S⁻¹, R, U, L, W, nlayers) = implicit
    ξ = implicit.ξ[]

    temp_tend = vars.tendencies.temperature
    pres_tend = vars.tendencies.pressure
    div_tend = vars.tendencies.divergence
    div_old, div_new = SpeedyWeather.get_steps(vars.prognostic.divergence)
    G = vars.scratch.a                  # reuse work arrays, used for combined tendency G
    geopotential = vars.scratch.b       # used for geopotential
    geopotential .= 0

    # Get precomputed l_indices from the spectrum
    l_indices = temp_tend.spectrum.l_indices

    arch = SpeedyWeather.architecture(temp_tend)

    # Single kernel: All implicit correction steps for each spectral mode
    SpeedyWeather.launch!(
        arch, SpeedyWeather.LinearWorkOrder, (size(pres_tend, 1),),
        SpeedyWeather.implicit_primitive_single_kernel!,
        temp_tend, pres_tend, div_tend, G, geopotential,
        div_old, div_new, S⁻¹, R, U, L, W, l_indices, ξ, nlayers
    )

    SpeedyWeather.zero_last_degree!(div_tend)
    SpeedyWeather.zero_last_degree!(pres_tend)
    SpeedyWeather.zero_last_degree!(temp_tend)

    return nothing
end
# END MOIST MASS COPY: implicit_correction!

"""Return the `stochastic_physics` component for `config`: the coupling when the option is on, else `nothing`."""
_atmosphere_stochastic_physics(config::ExperimentConfig) =
    config.atmosphere_moist_mass_coupling ? MoistMassCoupling() : nothing

"""Check the assembled model is a supported coupling host (also run at construction, after validation)."""
function _check_moist_mass_coupling_host(model)
    model.stochastic_physics isa MoistMassCoupling || return nothing
    convection = model.convection
    convection isa NetColumnBettsMillerConvection || throw(ArgumentError(
        "atmosphere_moist_mass_coupling requires the net-column Betts-Miller convection host, got $(typeof(convection))",
    ))
    model.large_scale_condensation isa NetColumnImplicitCondensation || throw(ArgumentError(
        "atmosphere_moist_mass_coupling requires NetColumnImplicitCondensation",
    ))
    last(model.parameterizations) == :stochastic_physics || throw(ArgumentError(
        "atmosphere_moist_mass_coupling requires stochastic_physics to be the last parameterization",
    ))
    model isa MoistMassCoupledModel || error("moist-mass coupled model type not selected")
    return nothing
end

"""Diagnostics provenance, written only when the option is on."""
function _moist_mass_coupling_provenance(config::ExperimentConfig)
    config.atmosphere_moist_mass_coupling || return (;)
    return (
        requested_atmosphere_moist_mass_coupling = true,
        atmosphere_moist_mass_coupling_scheme = MOIST_MASS_COUPLING_SCHEME,
        atmosphere_pressure_zero_mode_resets = "removed_surface_and_implicit_v1",
    )
end

"""
Reject a checkpoint whose moist-mass coupling setting differs from the reconstructed atmosphere, before any state is
changed. The record is written only when the coupling is on, so checkpoints without it (including every checkpoint
written before this option existed) imply the uncoupled atmosphere.
"""
function _validate_moist_mass_restart(simulation, state)
    saved = hasproperty(state, :moist_mass_coupling) ? state.moist_mass_coupling : false
    saved isa Bool || throw(ArgumentError("checkpoint moist-mass coupling setting must be Boolean"))
    active = simulation.model.stochastic_physics isa MoistMassCoupling
    saved == active || throw(ArgumentError(
        "checkpoint and reconstructed atmosphere use different moist-mass coupling settings " *
        "(checkpoint $(saved ? "coupled" : "uncoupled"), configuration $(active ? "coupled" : "uncoupled"))",
    ))
    return nothing
end

"""The restart record: present (and true) only when the coupling is on."""
_moist_mass_restart_state(simulation) =
    simulation.model.stochastic_physics isa MoistMassCoupling ? (; moist_mass_coupling = true) : (;)
