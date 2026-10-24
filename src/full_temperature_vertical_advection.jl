#=
Opt-in port of upstream SpeedyWeather PR #1286 (commit adf0d61d, merged 2026-10-20) to SpeedyWeather 0.21.1:
"Advect the full temperature in the vertical advection".

In SpeedyWeather 0.21.1, `dynamics_tendencies!` subtracts the semi-implicit reference profile T_k from the grid
temperature before `vertical_advection!`. The default `CenteredVerticalAdvection` advects the current step, so it
advects only T' = T - T_k, and the vertical advection of the reference profile, -sigmadot dT_k/dsigma, never enters
the explicit temperature tendency. `implicit_correction!` (step 1) assumes the explicit tendency contains the full
linear term L*D, of which this vertical advection is the L1 part, and only shifts it from step i to step i-1. Without
it most of the static-stability forcing by vertical motion cancels (upstream: the Jablonowski-Williamson baroclinic
wave grows about 12 % per day too slowly).

The fix moves the subtraction of T_k to after `vertical_advection!`, so the vertical advection sees the full
temperature. Between the old and the new position nothing reads the grid temperature in 0.21.1: `geopotential!`
(spectral virtual temperature), `vertical_integration!`, `surface_pressure_tendency!`, `vertical_velocity!` and
`linear_pressure_gradient!`.

Selection: `atmosphere_vertical_advection = :centered_full_temperature` installs
`FullTemperatureCenteredVerticalAdvection`, which is the centred second-order scheme of SpeedyWeather (every face
reconstruction delegates to `CenteredVerticalAdvection`). Models whose `vertical_advection` has this type get the
reordered `dynamics_tendencies!` below, including the moist-mass copy. The default
(`:speedyweather_0_21_1_centered`) installs SpeedyWeather's own `CenteredVerticalAdvection` and keeps the upstream
0.21.1 method, so the default is bitwise unchanged. The L1 algebra of the implicit operator and a refresh of the
reference profile are not part of this change.
=#

const ATMOSPHERE_VERTICAL_ADVECTION_OPTIONS = (:speedyweather_0_21_1_centered, :centered_full_temperature)
const FULL_TEMPERATURE_VERTICAL_ADVECTION_PROVENANCE =
    "centered_order2_full_temperature_speedyweather_pr1286_adf0d61d_ported_to_0_21_1"
const DEFAULT_VERTICAL_ADVECTION_PROVENANCE = "speedyweather_0_21_1_centered_order2_anomaly_temperature"

"""SpeedyWeather's centred second-order vertical advection, marked so that the full temperature is advected."""
struct FullTemperatureCenteredVerticalAdvection{NF, B} <: SpeedyWeather.DispersiveVerticalAdvection{NF, B} end

FullTemperatureCenteredVerticalAdvection(spectral_grid; order = 2) =
    FullTemperatureCenteredVerticalAdvection{spectral_grid.NF, order ÷ 2}()

@inline SpeedyWeather.reconstructed_at_face(ξ, ij, k, u, ::FullTemperatureCenteredVerticalAdvection{NF, B}) where {NF, B} =
    SpeedyWeather.reconstructed_at_face(ξ, ij, k, u, SpeedyWeather.CenteredVerticalAdvection{NF, B}())

function _atmosphere_vertical_advection_scheme(config, spectral_grid)
    config.atmosphere_vertical_advection == :centered_full_temperature &&
        return FullTemperatureCenteredVerticalAdvection(spectral_grid)
    config.atmosphere_vertical_advection == :speedyweather_0_21_1_centered &&
        return SpeedyWeather.CenteredVerticalAdvection(spectral_grid)
    error("unsupported atmosphere_vertical_advection $(config.atmosphere_vertical_advection)")
end

_atmosphere_vertical_advection_provenance(::FullTemperatureCenteredVerticalAdvection) =
    FULL_TEMPERATURE_VERTICAL_ADVECTION_PROVENANCE
_atmosphere_vertical_advection_provenance(::SpeedyWeather.CenteredVerticalAdvection{NF, 1}) where {NF} =
    DEFAULT_VERTICAL_ADVECTION_PROVENANCE
_atmosphere_vertical_advection_provenance(scheme) = string(nameof(typeof(scheme)))

"""PrimitiveWetModel with the given type parameters (by field name) restricted to the given upper bounds."""
function _restricted_wet_model_type(restrictions::Pair{Symbol, <:Type}...)
    body = Base.unwrap_unionall(SpeedyWeather.PrimitiveWetModel)
    params = collect(body.parameters)
    names = fieldnames(body)
    for (field, bound) in restrictions
        # Raw declared field type. fieldtype(body, field) with a non-constant field name is inferred as Type, so the
        # compiler folds `isa TypeVar` to false although the value is the TypeVar (found by Codex, 24 Oct).
        index = findfirst(==(field), names)
        isnothing(index) && error("PrimitiveWetModel has no field $field")
        tv = body.types[index]
        tv isa TypeVar || error("PrimitiveWetModel.$field is not a type parameter")
        bound <: tv.ub || error("$bound is outside the $field bound $(tv.ub)")
        i = findfirst(x -> x === tv, params)
        isnothing(i) && error("$field type parameter not found")
        params[i] = TypeVar(tv.name, tv.lb, bound)
    end
    T = SpeedyWeather.PrimitiveWetModel{params...}
    for x in reverse(params)
        T = UnionAll(x, T)
    end
    return T
end

const FullTemperatureAdvectionModel =
    _restricted_wet_model_type(:vertical_advection => FullTemperatureCenteredVerticalAdvection)
const MoistMassFullTemperatureModel = _restricted_wet_model_type(
    :vertical_advection => FullTemperatureCenteredVerticalAdvection, :stochastic_physics => MoistMassCoupling)

# BEGIN PR1286 COPY: dynamics_tendencies!
# SpeedyWeather 0.21.1 dynamics/tendencies.jl dynamics_tendencies!(vars, lf, model::PrimitiveEquation), with the
# subtraction of the reference profile moved after vertical_advection! (upstream PR #1286). No other change.
function SpeedyWeather.dynamics_tendencies!(
        vars::SpeedyWeather.Variables,
        lf::Integer,                # leapfrog index for tendencies
        model::FullTemperatureAdvectionModel,
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

    # from ∂Φ/∂ln(pₛ) = -RTᵥ for bernoulli_potential!
    SpeedyWeather.geopotential!(vars, geopotential, orography)

    # get ū, v̄, D̄ on grid; D̄ in spectral
    SpeedyWeather.vertical_integration!(vars, lf_implicit, geometry)

    # ∂ln(pₛ)/∂t = -(ū, v̄)⋅∇ln(pₛ) - D̄
    SpeedyWeather.surface_pressure_tendency!(vars, spectral_transform)

    # calculate vertical velocity σ̇ in sigma coordinates for the vertical mass flux M = pₛ * σ̇
    SpeedyWeather.vertical_velocity!(vars, geometry)

    # add the RTₖlnpₛ term to geopotential
    SpeedyWeather.linear_pressure_gradient!(vars, lf_implicit, atmosphere, implicit)

    # use σ̇ for the vertical advection of u, v, T, q. PR #1286: this advects the full temperature T = T' + Tₖ, so it
    # includes the vertical advection of the reference profile Tₖ, which implicit_correction! expects here.
    SpeedyWeather.vertical_advection!(vars, model)

    # temperature relative to profile
    vars.grid.temperature.data .-= implicit.temp_profile'

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
# END PR1286 COPY: dynamics_tendencies!

# BEGIN PR1286 MOIST MASS COPY: dynamics_tendencies!
# The moist-mass copy (moist_mass_coupling.jl) with the same single reordering. More specific than both the
# moist-mass method and the method above, so a model with both options dispatches here without ambiguity.
function SpeedyWeather.dynamics_tendencies!(
        vars::SpeedyWeather.Variables,
        lf::Integer,                # leapfrog index for tendencies
        model::MoistMassFullTemperatureModel,
    )
    SpeedyWeather.forcing!(vars, lf, model)
    SpeedyWeather.drag!(vars, lf, model)

    (; orography, geometry, spectral_transform, geopotential, atmosphere, implicit) = model

    lf_implicit = implicit.α == 0 ? lf : 1

    SpeedyWeather.pressure_gradient_flux!(vars, lf, spectral_transform)
    SpeedyWeather.linear_virtual_temperature!(vars, lf_implicit, model)
    SpeedyWeather.geopotential!(vars, geopotential, orography)
    SpeedyWeather.vertical_integration!(vars, lf_implicit, geometry)
    _moist_mass_surface_pressure_tendency!(vars, spectral_transform)
    SpeedyWeather.vertical_velocity!(vars, geometry)
    SpeedyWeather.linear_pressure_gradient!(vars, lf_implicit, atmosphere, implicit)

    # PR #1286: advect the full temperature, then work with the anomaly
    SpeedyWeather.vertical_advection!(vars, model)
    vars.grid.temperature.data .-= implicit.temp_profile'

    SpeedyWeather.vordiv_tendencies!(vars, model)
    SpeedyWeather.temperature_tendency!(vars, model)
    SpeedyWeather.humidity_tendency!(vars, model)
    SpeedyWeather.bernoulli_potential!(vars, spectral_transform)
    SpeedyWeather.tracer_advection!(vars, model)

    vars.grid.temperature.data .+= implicit.temp_profile'

    return nothing
end
# END PR1286 MOIST MASS COPY: dynamics_tendencies!
