"""
BS01/BS02 shared convective detrainment interface, v1 (prototype, CPU, Base only).

One column, one convective call. Convection (BS01 plume, or the current
Betts-Miller scheme through an identity adapter) reports what it did to the
environment as layer TENDENCIES over the call interval `dt`, plus surface
precipitation FLUXES. BS02's transported cloud reads only `dql` and `dqi`.
Nothing else crosses the interface, so condensate is defined exactly once.

Layers run k = 1 (top) ... n (surface), as in SpeedyWeather. The layer air mass
is `m[k] = p_s Δσ[k] / g` [kg m-2].

Fields, all per layer in kg kg-1 s-1 or K s-1 unless stated:
  dqv  vapour tendency (any sign)
  dql  environment liquid condensate tendency (detrainment, minus any
       entrainment/subsidence of existing condensate); post-call value >= 0
  dqi  the same for ice
  dT   temperature tendency, including every latent heat release
  rain, snow  surface precipitation fluxes [kg m-2 s-1] (>= 0)

Conservation contracts (column integrals, no surface exchange inside the call):
  W:  Σ m (dqv + dql + dqi) + rain + snow = 0
  E:  Σ m (cp dT + Lv dqv - Lf dqi) = Lf snow
E is the frozen moist enthalpy with liquid water as the reference phase. Rain
leaving the column carries no E. Snow leaving removes -Lf per kg of ice, so the
column gains Lf*snow. Sensible heat carried by precipitation is neglected, as
in SpeedyWeather and ReadyESM.

Phase: one function for every component that phase-partitions condensate,
`liquid_fraction(T)`, equal to the mixed_phase saturation blend weight
(1 at and above 273.15 K, 0 at and below 250.16 K, quadratic between). Today
the stack has four different rules:
  - saturation blend 250.16-273.15 K (src/mixed_phase_saturation.jl);
  - large-scale snow below 263 K (flagged_condensation.jl:63,227);
  - convective precipitation always rain (convection.jl:136);
  - radiation cloud phase linear between 253 and 273 K
    (rrtmgp_radiation.jl `_set_diagnostic_cloud_column!`).
The interface fixes the rule only for new code. Changing an existing scheme
is a separate, reviewed decision.

Positivity: q + dt*dq >= 0 for vapour, liquid and ice after the call.
Each `check_*` function returns a NamedTuple of residuals.
`assert_contracts` applies tolerances scaled by the column's gross exchange,
not by its residual.
"""
module ConvectiveDetrainmentInterface

export ColumnExchange, liquid_fraction, dliquid_fraction_dT, latent_heat_released, layer_mass,
       check_water, check_enthalpy, check_positivity, assert_contracts,
       identity_adapter

const T_FREEZE = 273.15
const T_ALL_ICE = 250.16

"""Liquid fraction shared with the mixed_phase saturation law."""
@inline function liquid_fraction(T::Real)
    T >= T_FREEZE && return one(T)
    a = clamp((T - T_ALL_ICE) / (T_FREEZE - T_ALL_ICE), zero(T), one(T))
    return a * a
end

"""Analytic derivative of `liquid_fraction` with respect to temperature (zero outside the blend)."""
@inline function dliquid_fraction_dT(T::Real)
    x = (T - T_ALL_ICE) / (T_FREEZE - T_ALL_ICE)
    return (zero(T) < x < one(T)) ? 2x / (T_FREEZE - T_ALL_ICE) : zero(T)
end

"""Heat released per kg of vapour condensed at temperature T, using the shared phase rule."""
@inline latent_heat_released(T, Lv, Lf) = Lv + (1 - liquid_fraction(T)) * Lf

struct ColumnExchange{V <: AbstractVector, S <: Real}
    dqv::V
    dql::V
    dqi::V
    dT::V
    rain::S
    snow::S
end

function ColumnExchange(n::Integer, ::Type{NF} = Float64) where {NF}
    z() = zeros(NF, n)
    return ColumnExchange(z(), z(), z(), z(), zero(NF), zero(NF))
end

layer_mass(ps, dσ, g) = ps .* dσ ./ g

"""Column water residual [kg m-2 s-1] and the gross exchange used to scale it."""
function check_water(x::ColumnExchange, m)
    net = sum(m .* (x.dqv .+ x.dql .+ x.dqi)) + x.rain + x.snow
    gross = sum(m .* (abs.(x.dqv) .+ abs.(x.dql) .+ abs.(x.dqi))) + x.rain + x.snow
    return (; residual = net, gross)
end

"""Column frozen moist enthalpy residual [W m-2] and its gross scale."""
function check_enthalpy(x::ColumnExchange, m; cp, Lv, Lf)
    lhs = sum(m .* (cp .* x.dT .+ Lv .* x.dqv .- Lf .* x.dqi))
    gross = sum(m .* (cp .* abs.(x.dT) .+ Lv .* abs.(x.dqv) .+ Lf .* abs.(x.dqi))) + Lf * x.snow
    return (; residual = lhs - Lf * x.snow, gross)
end

"""
Most negative post-call vapour and condensate, and the sign of the surface fluxes.
`ql`, `qi` are the environment condensate before the call (zero if not carried):
convection may entrain or subside existing condensate, so dql/dqi can be negative,
but the post-call values may not.
"""
function check_positivity(x::ColumnExchange, qv, dt; ql = zero(qv), qi = zero(qv))
    min_qv = minimum(qv .+ dt .* x.dqv)
    min_qc = min(minimum(ql .+ dt .* x.dql), minimum(qi .+ dt .* x.dqi))
    return (; min_qv, min_qc, condensate_nonnegative = min_qc >= 0,
            fluxes_nonnegative = x.rain >= 0 && x.snow >= 0)
end

"""Throw if any contract fails. `rtol` multiplies the gross exchange."""
function assert_contracts(x::ColumnExchange, qv, m, dt; cp, Lv, Lf, rtol = 1e-6,
                          ql = zero(qv), qi = zero(qv))
    w = check_water(x, m)
    e = check_enthalpy(x, m; cp, Lv, Lf)
    p = check_positivity(x, qv, dt; ql, qi)
    abs(w.residual) <= rtol * max(w.gross, eps()) || error("water contract: $(w)")
    abs(e.residual) <= rtol * max(e.gross, eps()) || error("enthalpy contract: $(e)")
    p.min_qv >= -rtol * maximum(abs, qv) || error("vapour positivity: $(p.min_qv)")
    p.min_qc >= -rtol * max(maximum(abs, qv), eps()) || error("condensate positivity: $(p.min_qc)")
    p.fluxes_nonnegative || error("negative surface precipitation")
    return (; water = w, enthalpy = e, positivity = p)
end

"""
Identity adapter for an existing scheme that detrains nothing (Betts-Miller as
run in ReadyESM). Given its temperature and vapour tendencies and its NET
column rain (NetColumnBettsMillerConvection), report them unchanged with
dql = dqi = 0 and snow = 0. With latent heat Lv only, the E contract then holds
exactly when cp ΣmdT = -Lv Σm dqv, which is the identity a conserving
Betts-Miller column satisfies.
"""
function identity_adapter(dT, dqv, rain)
    n = length(dT)
    z = zero(eltype(dT))
    return ColumnExchange(collect(dqv), fill(z, n), fill(z, n), collect(dT), rain, z)
end

end # module
