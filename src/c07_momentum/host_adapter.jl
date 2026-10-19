module C07HostAdapter


include(joinpath(@__DIR__, "flux_form_cmt.jl"))
using .FluxFormCMT

const R = parentmodule(@__MODULE__)
const SW = R.SpeedyWeather
const S = SW
const C08 = R.FiniteC08Core

export C07FluxFormCMTConvection,
    subcycled_transport_momentum

function _zero_host_exchange(::Type{T}, n, raw_step_number; skipped) where {T}
    zeros_n = zeros(T, n)
    return (
        du = copy(zeros_n),
        dv = copy(zeros_n),
        final_u = copy(zeros_n),
        final_v = copy(zeros_n),
        raw_step_number = T(raw_step_number),
        maximum_applied_step_number = zero(T),
        subcycles = 0,
        skipped,
        momentum_residual_u = zero(T),
        momentum_residual_v = zero(T),
        momentum_gross_u = zero(T),
        momentum_gross_v = zero(T),
        kinetic_energy_dissipation_W_m2 = zero(T),
    )
end

"""Apply C07 over equal convex substeps and return one host-step tendency."""
function subcycled_transport_momentum(
    u,
    v,
    layer_mass,
    Mhalf,
    entrainment,
    detrainment,
    mass_flux,
    effective_step_seconds;
    strength = one(mass_flux),
    maximum_subcycles::Int = 16,
)
    n = length(u)
    maximum_subcycles > 0 || throw(ArgumentError(
        "C07 maximum_subcycles must be positive",
    ))
    T = promote_type(
        eltype(u),
        eltype(v),
        eltype(layer_mass),
        typeof(mass_flux),
        typeof(effective_step_seconds),
        typeof(strength),
    )
    u_initial = T.(u)
    v_initial = T.(v)
    mass = T.(layer_mass)
    entrainment_typed = T.(entrainment)
    effective_mass_flux = T(mass_flux) * T(strength)
    full_step = T(effective_step_seconds)
    raw_step_number = maximum(
        effective_mass_flux .* entrainment_typed .* full_step ./ mass,
    )

    if iszero(effective_mass_flux)
        exchange = _zero_host_exchange(T, n, raw_step_number; skipped = false)
        return merge(exchange, (; final_u = u_initial, final_v = v_initial))
    end
    if !isfinite(raw_step_number) || raw_step_number > T(maximum_subcycles)
        exchange = _zero_host_exchange(T, n, raw_step_number; skipped = true)
        return merge(exchange, (; final_u = u_initial, final_v = v_initial))
    end

    subcycles = max(1, ceil(Int, raw_step_number))
    substep = full_step / T(subcycles)
    if subcycles == 1
        exchange = FluxFormCMT.transport_momentum(
            u_initial,
            v_initial,
            mass,
            Mhalf,
            entrainment_typed,
            detrainment,
            mass_flux,
            full_step;
            strength,
        )
        final_u = u_initial .+ full_step .* exchange.du
        final_v = v_initial .+ full_step .* exchange.dv
        closure = FluxFormCMT.check_momentum(exchange, mass)
        old_energy = sum(
            mass .* (u_initial .^ 2 .+ v_initial .^ 2),
        ) / 2
        new_energy = sum(mass .* (final_u .^ 2 .+ final_v .^ 2)) / 2
        return (
            du = exchange.du,
            dv = exchange.dv,
            final_u,
            final_v,
            raw_step_number,
            maximum_applied_step_number =
                exchange.max_entrainment_step_number,
            subcycles,
            skipped = false,
            momentum_residual_u = closure.residual_u,
            momentum_residual_v = closure.residual_v,
            momentum_gross_u = closure.gross_u,
            momentum_gross_v = closure.gross_v,
            kinetic_energy_dissipation_W_m2 =
                (old_energy - new_energy) / full_step,
        )
    end
    u_work = copy(u_initial)
    v_work = copy(v_initial)
    maximum_applied_step_number = zero(T)
    for _ in 1:subcycles
        exchange = FluxFormCMT.transport_momentum(
            u_work,
            v_work,
            mass,
            Mhalf,
            entrainment_typed,
            detrainment,
            mass_flux,
            substep;
            strength,
        )
        maximum_applied_step_number = max(
            maximum_applied_step_number,
            exchange.max_entrainment_step_number,
        )
        u_work .+= substep .* exchange.du
        v_work .+= substep .* exchange.dv
    end

    du = (u_work .- u_initial) ./ full_step
    dv = (v_work .- v_initial) ./ full_step
    closure = FluxFormCMT.check_momentum((; du, dv), mass)
    old_energy = sum(mass .* (u_initial .^ 2 .+ v_initial .^ 2)) / 2
    new_energy = sum(mass .* (u_work .^ 2 .+ v_work .^ 2)) / 2
    return (
        du,
        dv,
        final_u = u_work,
        final_v = v_work,
        raw_step_number,
        maximum_applied_step_number,
        subcycles,
        skipped = false,
        momentum_residual_u = closure.residual_u,
        momentum_residual_v = closure.residual_v,
        momentum_gross_u = closure.gross_u,
        momentum_gross_v = closure.gross_v,
        kinetic_energy_dissipation_W_m2 =
            (old_energy - new_energy) / full_step,
    )
end

struct C07FluxFormCMTConvection{B, F} <: SW.AbstractConvection
    base::B
    strength::F
    maximum_subcycles::Int
end

function C07FluxFormCMTConvection(
    base::R.FiniteMultiPlumeConvection;
    strength::Real = 1,
    maximum_subcycles::Int = 16,
)
    isfinite(strength) && strength >= 0 || throw(ArgumentError(
        "C07 strength must be finite and non-negative",
    ))
    maximum_subcycles > 0 || throw(ArgumentError(
        "C07 maximum_subcycles must be positive",
    ))
    FT = eltype(base.spectrum.weights)
    isfinite(FT(strength)) || throw(ArgumentError(
        "C07 strength must remain finite in the atmosphere storage precision",
    ))
    return C07FluxFormCMTConvection(
        base,
        FT(strength),
        maximum_subcycles,
    )
end

SW.initialize!(scheme::C07FluxFormCMTConvection, model::SW.AbstractModel) =
    SW.initialize!(scheme.base, model)

SW.variables(scheme::C07FluxFormCMTConvection) = (
    SW.variables(scheme.base)...,
    SW.ParameterizationVariable(
        :c07_triggered,
        SW.Grid2D();
        desc = "One when C07 applies momentum transport",
        units = "1",
    ),
    SW.ParameterizationVariable(
        :c07_subcycles,
        SW.Grid2D();
        desc = "Number of C07 momentum subcycles",
        units = "1",
    ),
    SW.ParameterizationVariable(
        :c07_subcycled,
        SW.Grid2D();
        desc = "One when C07 uses more than one subcycle",
        units = "1",
    ),
    SW.ParameterizationVariable(
        :c07_skipped,
        SW.Grid2D();
        desc = "One when C07 refuses a pathological subcycle request",
        units = "1",
    ),
    SW.ParameterizationVariable(
        :c07_raw_step_number,
        SW.Grid2D();
        desc = "Full-step explicit entrainment number",
        units = "1",
    ),
    SW.ParameterizationVariable(
        :c07_applied_step_number,
        SW.Grid2D();
        desc = "Largest applied substep entrainment number",
        units = "1",
    ),
    SW.ParameterizationVariable(
        :c07_kinetic_energy_dissipation,
        SW.Grid2D();
        desc = "C07 kinetic-energy dissipation",
        units = "W/m2",
    ),
    SW.ParameterizationVariable(
        :c07_momentum_residual_u,
        SW.Grid2D();
        desc = "C07 zonal column-momentum tendency residual",
        units = "N/m2",
    ),
    SW.ParameterizationVariable(
        :c07_momentum_residual_v,
        SW.Grid2D();
        desc = "C07 meridional column-momentum tendency residual",
        units = "N/m2",
    ),
    SW.ParameterizationVariable(
        :c07_u_tendency,
        SW.Grid3D();
        desc = "C07 zonal wind tendency",
        units = "m/s2",
    ),
    SW.ParameterizationVariable(
        :c07_v_tendency,
        SW.Grid3D();
        desc = "C07 meridional wind tendency",
        units = "m/s2",
    ),
)

@inline function _reset_c07_column!(parameterizations, ij, number_layers, zero_value)
    parameterizations.c07_triggered[ij] = zero_value
    parameterizations.c07_subcycles[ij] = zero_value
    parameterizations.c07_subcycled[ij] = zero_value
    parameterizations.c07_skipped[ij] = zero_value
    parameterizations.c07_raw_step_number[ij] = zero_value
    parameterizations.c07_applied_step_number[ij] = zero_value
    parameterizations.c07_kinetic_energy_dissipation[ij] = zero_value
    parameterizations.c07_momentum_residual_u[ij] = zero_value
    parameterizations.c07_momentum_residual_v[ij] = zero_value
    for layer in 1:number_layers
        parameterizations.c07_u_tendency[ij, layer] = zero_value
        parameterizations.c07_v_tendency[ij, layer] = zero_value
    end
    return nothing
end

function _c07_column_and_composite(ij, vars, scheme, model, ::Type{NF}) where {NF}
    (; geometry, planet) = model
    temperature = vars.grid.temperature_prev
    humidity = vars.grid.humidity_prev
    geopotential = vars.grid.geopotential
    surface_pressure = NF(vars.grid.pressure_prev[ij])
    sigma = geometry.σ_levels_full
    sigma_thickness = geometry.σ_levels_thick
    number_layers = length(sigma)
    column = C08.Column(
        NF[temperature[ij, layer] for layer in 1:number_layers],
        NF[max(humidity[ij, layer], zero(NF)) for layer in 1:number_layers],
        zeros(NF, number_layers),
        NF[surface_pressure * sigma[layer] for layer in 1:number_layers],
        NF[
            surface_pressure * sigma_thickness[layer]
            for layer in 1:number_layers
        ],
        NF[geopotential[ij, layer] for layer in 1:number_layers],
    )
    thermo = R._finite_multiplume_thermo(scheme.base, model, NF)
    col, th, spectrum = R.C08FiniteConvection.MixedPrecisionC08.wide_inputs(
        column, thermo, scheme.base.spectrum)
    parameters = [C08.bin_parameters(spectrum, b) for b in eachindex(spectrum.weights)]
    profiles = [C08.EP.plume_profile(col, th, p) for p in parameters]
    active = [C08.bin_triggers(profiles[b], parameters[b]) && spectrum.weights[b] > 0
        for b in eachindex(parameters)]
    wide_composite = C08.composite_profile(profiles, spectrum.weights, active)
    cast(x::AbstractArray{<:AbstractFloat}) = NF.(x)
    cast(x::AbstractFloat) = NF(x)
    cast(x) = x
    names = propertynames(wide_composite)
    composite = NamedTuple{names}(Tuple(cast(getproperty(wide_composite, k)) for k in names))
    return (; column, thermo, active, composite)
end

Base.@propagate_inbounds function SW.parameterization!(
    ij,
    vars,
    scheme::C07FluxFormCMTConvection,
    model,
)
    number_layers = length(model.geometry.σ_levels_full)
    NF = eltype(vars.grid.temperature_prev)
    zero_value = zero(NF)
    parameterizations = vars.parameterizations
    _reset_c07_column!(parameterizations, ij, number_layers, zero_value)

    SW.parameterization!(ij, vars, scheme.base, model)
    iszero(scheme.strength) && return nothing
    parameterizations.c08_triggered[ij] > 0 || return nothing

    diagnosed = _c07_column_and_composite(ij, vars, scheme, model, NF)
    expected_active = round(Int, parameterizations.c08_active_plumes[ij])
    count(diagnosed.active) == expected_active || error(
        "C07 recomputed active-plume count differs from C08 host diagnosis",
    )
    u = NF[
        vars.grid.u_prev[ij, layer] for layer in 1:number_layers
    ]
    v = NF[
        vars.grid.v_prev[ij, layer] for layer in 1:number_layers
    ]
    layer_mass = diagnosed.column.dp ./ diagnosed.thermo.g
    timestep = NF(model.time_stepping.Δt_sec)
    effective_step =
        scheme.base.spectrum.shared.effective_step_factor * timestep
    exchange = subcycled_transport_momentum(
        u,
        v,
        layer_mass,
        diagnosed.composite.Mhalf,
        diagnosed.composite.E,
        diagnosed.composite.D,
        NF(parameterizations.c08_mass_flux[ij]),
        effective_step;
        strength = scheme.strength,
        maximum_subcycles = scheme.maximum_subcycles,
    )

    parameterizations.c07_triggered[ij] = one(NF)
    parameterizations.c07_subcycles[ij] = NF(exchange.subcycles)
    parameterizations.c07_subcycled[ij] = NF(exchange.subcycles > 1)
    parameterizations.c07_skipped[ij] = NF(exchange.skipped)
    parameterizations.c07_raw_step_number[ij] = exchange.raw_step_number
    parameterizations.c07_applied_step_number[ij] =
        exchange.maximum_applied_step_number
    parameterizations.c07_kinetic_energy_dissipation[ij] =
        exchange.kinetic_energy_dissipation_W_m2
    parameterizations.c07_momentum_residual_u[ij] =
        exchange.momentum_residual_u
    parameterizations.c07_momentum_residual_v[ij] =
        exchange.momentum_residual_v
    exchange.skipped && return nothing

    for layer in 1:number_layers
        vars.tendencies.grid.u[ij, layer] += exchange.du[layer]
        vars.tendencies.grid.v[ij, layer] += exchange.dv[layer]
        parameterizations.c07_u_tendency[ij, layer] = exchange.du[layer]
        parameterizations.c07_v_tendency[ij, layer] = exchange.dv[layer]
    end
    return nothing
end

include("thermal_return.jl")

end
