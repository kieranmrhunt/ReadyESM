"""
Column-local flux-form convective momentum transport for the C07 Stage 0
component gate.

The operator receives the entrainment, detrainment and interface mass-flux
profile of an independently diagnosed plume. It mixes horizontal momentum
into that plume, returns detrained momentum to the environment and applies the
same top-down implicit compensating subsidence used by the qualified C08
thermodynamic exchange. There is no relaxation to a layer-mean wind.
"""
module FluxFormCMT

export check_mass_continuity,
    check_momentum,
    kinetic_energy_change,
    transport_momentum

function _require_vector_lengths(n, vectors...)
    all(length(vector) == n for vector in vectors) ||
        throw(ArgumentError("C07 column-vector lengths differ"))
    return nothing
end

function _require_finite(name, values)
    all(isfinite, values) || throw(ArgumentError("non-finite C07 $name"))
    return nothing
end

"""Return layer residuals of `Mbelow + E = Mabove + D`."""
function check_mass_continuity(Mhalf, entrainment, detrainment)
    n = length(entrainment)
    length(detrainment) == n ||
        throw(ArgumentError("C07 entrainment and detrainment lengths differ"))
    length(Mhalf) == n + 1 ||
        throw(ArgumentError("C07 Mhalf must have n + 1 interfaces"))
    return [
        Mhalf[layer + 1] + entrainment[layer] -
        Mhalf[layer] - detrainment[layer]
        for layer in 1:n
    ]
end

function _validate_profile(layer_mass, Mhalf, entrainment, detrainment)
    n = length(layer_mass)
    _require_vector_lengths(n, entrainment, detrainment)
    length(Mhalf) == n + 1 ||
        throw(ArgumentError("C07 Mhalf must have n + 1 interfaces"))
    _require_finite("layer mass", layer_mass)
    _require_finite("interface mass flux", Mhalf)
    _require_finite("entrainment", entrainment)
    _require_finite("detrainment", detrainment)
    all(>(0), layer_mass) ||
        throw(ArgumentError("C07 layer masses must be positive"))
    all(>=(0), Mhalf) ||
        throw(ArgumentError("C07 interface mass flux must be non-negative"))
    all(>=(0), entrainment) ||
        throw(ArgumentError("C07 entrainment must be non-negative"))
    all(>=(0), detrainment) ||
        throw(ArgumentError("C07 detrainment must be non-negative"))

    residual = check_mass_continuity(Mhalf, entrainment, detrainment)
    scale = max(
        maximum(abs, Mhalf),
        maximum(abs, entrainment),
        maximum(abs, detrainment),
        eps(eltype(residual)),
    )
    tolerance = 64 * eps(eltype(residual)) * scale
    maximum(abs, residual) <= tolerance || throw(ArgumentError(
        "C07 plume mass continuity failed: " *
        "maximum residual $(maximum(abs, residual)), tolerance $tolerance",
    ))
    abs(Mhalf[1]) <= tolerance ||
        throw(ArgumentError("C07 plume carries mass through the model top"))
    abs(Mhalf[end]) <= tolerance ||
        throw(ArgumentError("C07 plume carries mass through the surface"))
    return (; residual, tolerance)
end

function _plume_momentum(environment, Mhalf, entrainment, detrainment)
    n = length(environment)
    updraft = similar(environment)
    fill!(updraft, zero(eltype(environment)))
    for layer in n:-1:1
        incoming_mass = Mhalf[layer + 1]
        mixed_mass = incoming_mass + entrainment[layer]
        if mixed_mass > 0
            incoming_momentum = layer == n ? zero(mixed_mass) :
                incoming_mass * updraft[layer + 1]
            updraft[layer] = (
                incoming_momentum +
                entrainment[layer] * environment[layer]
            ) / mixed_mass
        elseif detrainment[layer] > 0 || Mhalf[layer] > 0
            throw(ArgumentError(
                "C07 profile has outgoing plume mass without incoming mass",
            ))
        else
            updraft[layer] = environment[layer]
        end
    end
    return updraft
end

function _environment_tendency(
    environment,
    updraft,
    layer_mass,
    Mhalf,
    entrainment,
    detrainment,
    mass_flux,
    effective_step_seconds,
)
    n = length(environment)
    tendency = similar(environment)
    previous_new = zero(eltype(environment))
    for layer in 1:n
        source = mass_flux * (
            detrainment[layer] * updraft[layer] -
            entrainment[layer] * environment[layer]
        )
        flux_from_above = mass_flux * Mhalf[layer]
        flux_to_below = mass_flux * Mhalf[layer + 1]
        # Algebraically identical implicit exchange, evaluated directly as
        # a tendency. Avoid losing weak transport in (new_value - old_value).
        numerator = source + flux_from_above * previous_new -
                    flux_to_below * environment[layer]
        tendency[layer] = numerator /
            (layer_mass[layer] + effective_step_seconds * flux_to_below)
        previous_new = environment[layer] +
                       effective_step_seconds * tendency[layer]
    end
    return tendency
end

"""
Apply flux-form convective momentum transport to one column.

`Mhalf`, `entrainment` and `detrainment` are the dimensionless plume profile
for unit cloud-base mass flux. `mass_flux` has units kg m-2 s-1. The returned
`du` and `dv` have units m s-2. `strength = 0` returns exact zero tendencies.
"""
function transport_momentum(
    u,
    v,
    layer_mass,
    Mhalf,
    entrainment,
    detrainment,
    mass_flux,
    effective_step_seconds;
    strength = one(mass_flux),
)
    n = length(u)
    _require_vector_lengths(n, v, layer_mass, entrainment, detrainment)
    length(Mhalf) == n + 1 ||
        throw(ArgumentError("C07 Mhalf must have n + 1 interfaces"))
    _require_finite("zonal wind", u)
    _require_finite("meridional wind", v)
    isfinite(mass_flux) && mass_flux >= 0 ||
        throw(ArgumentError("C07 cloud-base mass flux must be finite and non-negative"))
    isfinite(strength) && strength >= 0 ||
        throw(ArgumentError("C07 strength must be finite and non-negative"))
    isfinite(effective_step_seconds) && effective_step_seconds > 0 ||
        throw(ArgumentError("C07 effective step must be finite and positive"))
    continuity = _validate_profile(
        layer_mass,
        Mhalf,
        entrainment,
        detrainment,
    )
    T = promote_type(
        eltype(u),
        eltype(v),
        eltype(layer_mass),
        typeof(mass_flux),
        typeof(effective_step_seconds),
        typeof(strength),
    )
    effective_mass_flux = T(mass_flux) * T(strength)
    step = T(effective_step_seconds)
    mass_typed = T.(layer_mass)
    entrainment_typed = T.(entrainment)
    max_entrainment_step_number = maximum(
        effective_mass_flux .* entrainment_typed .* step ./ mass_typed,
    )
    if !isfinite(max_entrainment_step_number) ||
            max_entrainment_step_number > one(T)
        throw(ArgumentError(
            "C07 explicit entrainment step number exceeds one: " *
            "$max_entrainment_step_number",
        ))
    end
    if iszero(effective_mass_flux)
        return (
            du = fill(zero(T), n),
            dv = fill(zero(T), n),
            updraft_u = T.(u),
            updraft_v = T.(v),
            mass_continuity = continuity,
            effective_mass_flux,
            max_entrainment_step_number,
        )
    end

    # Work in a column-local translating frame. Constant wind must produce
    # zero exchange, including when the supplied mass profile is rounded.
    # Return plume winds in the original frame after evaluating tendencies.
    reference_u, reference_v = T(u[end]), T(v[end])
    u_typed = T.(u) .- reference_u
    v_typed = T.(v) .- reference_v
    Mhalf_typed = T.(Mhalf)
    detrainment_typed = T.(detrainment)
    updraft_u = _plume_momentum(
        u_typed,
        Mhalf_typed,
        entrainment_typed,
        detrainment_typed,
    )
    updraft_v = _plume_momentum(
        v_typed,
        Mhalf_typed,
        entrainment_typed,
        detrainment_typed,
    )
    du = _environment_tendency(
        u_typed,
        updraft_u,
        mass_typed,
        Mhalf_typed,
        entrainment_typed,
        detrainment_typed,
        effective_mass_flux,
        step,
    )
    dv = _environment_tendency(
        v_typed,
        updraft_v,
        mass_typed,
        Mhalf_typed,
        entrainment_typed,
        detrainment_typed,
        effective_mass_flux,
        step,
    )
    return (
        du,
        dv,
        updraft_u = updraft_u .+ reference_u,
        updraft_v = updraft_v .+ reference_v,
        mass_continuity = continuity,
        effective_mass_flux,
        max_entrainment_step_number,
    )
end

"""Column momentum residuals and gross exchange in tendency units."""
function check_momentum(exchange, layer_mass)
    residual_u = sum(layer_mass .* exchange.du)
    residual_v = sum(layer_mass .* exchange.dv)
    gross_u = sum(layer_mass .* abs.(exchange.du))
    gross_v = sum(layer_mass .* abs.(exchange.dv))
    return (; residual_u, residual_v, gross_u, gross_v)
end

"""Finite-step change of column kinetic energy in J m-2."""
function kinetic_energy_change(exchange, u, v, layer_mass, step_seconds)
    u_new = u .+ step_seconds .* exchange.du
    v_new = v .+ step_seconds .* exchange.dv
    old_energy = sum(
        layer_mass .* (u .^ 2 .+ v .^ 2),
    ) / 2
    new_energy = sum(
        layer_mass .* (u_new .^ 2 .+ v_new .^ 2),
    ) / 2
    return new_energy - old_energy
end

end
