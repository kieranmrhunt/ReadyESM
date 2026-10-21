# Copied unchanged from artifacts/development_wave_20261020/c10_moist_mass_coupling_v2/column.jl (sha256 7e931969...),
# the accepted column reference for the moist-mass coupling. Used by test/moist_mass_coupling.jl only.
# Moist-mass coupling in sigma coordinates: one-column reference implementation (Claude, prototype).
# Layers k = 1..K run top to bottom. Δσ sums to 1. Layer mass m_k = Δσ_k p_s / g. Specific humidity q_k is vapour mass
# per moist-air mass. Condensate leaves the air where it forms.
#
# Physics supplies the net vapour source S_k (kg m-2 s-1) in each layer: surface evaporation E in the bottom layer,
# minus condensation, plus re-evaporation. The column source is Σ = sum(S) = E - P. With mass coupling:
#   d p_s/dt = g Σ                                          (the surface pressure carries the vapour mass)
#   d m_k/dt = Δσ_k Σ                                        (sigma keeps each layer's share of the column)
#   F_{k+1/2} = Σ_{j<=k} S_j - σ_{k+1/2} Σ                  (interface mass flux, positive downward)
# so that d m_k/dt = S_k + F_{k-1/2} - F_{k+1/2}, with F_{1/2} = F_{K+1/2} = 0. F is a coordinate flux: air crosses
# the sigma surfaces because those surfaces are fixed fractions of a changing p_s.
# Each carried quantity c (q, T, u, v) obeys d(m c)_k/dt = S⁺_k c_src + S⁻_k c_k + F_{k-1/2} ĉ_{k-1/2} - F_{k+1/2} ĉ_{k+1/2}
# with upwind interface values ĉ. S⁺ are additions and S⁻ removals. For vapour both carry c = 1 (pure water), so
# d(m q)/dt = S⁺ - S⁻ + interface terms and m dq/dt = (S⁺ - S⁻)(1 - q) + interface terms: the dilution term.
# Temperature additions carry the source temperature (surface evaporation: T_evap, re-evaporation: the layer value).
# Momentum additions carry zero for surface evaporation and the layer value for re-evaporation. Removed condensate
# leaves with the layer's temperature and momentum, so removals change neither.
# True mass sources and internal transport are separate inputs. A native scheme can also redistribute vapour
# between layers by exchanging air (convective or diffusive transport). Such a tendency T_k has sum(m .* T) = 0, it
# changes no layer's air mass, and the native q tendency already describes it. It is NOT a mass source, and the
# coupling correction must not depend on it (Codex review). Only E, condensate removal and re-evaporation are sources.
module MoistMassColumn

export ColumnState, coupled_step, uncoupled_step, coupled_tendencies, mass_coupling_correction, budgets

struct ColumnState{F}
    ps::F
    q::Vector{F}
    T::Vector{F}
    u::Vector{F}
end

layer_mass(ps, Δσ, g) = Δσ .* (ps / g)

"Interface fluxes F_{k+1/2}, k = 0..K (F[1] is the top, F[K+1] the surface; both zero by construction)."
function interface_flux(S, Δσ)
    K = length(S)
    Σ = sum(S)
    F = zeros(eltype(S), K + 1)
    cumulative_S, cumulative_σ = zero(eltype(S)), zero(eltype(S))
    for k in 1:K-1
        cumulative_S += S[k]
        cumulative_σ += Δσ[k]
        F[k+1] = cumulative_S - cumulative_σ * Σ
    end
    return F
end

upwind(F, above, below) = F > 0 ? above : below

"""
Exact finite-volume update over dt. `S_add` are additions (re-evaporation, plus surface evaporation `E` in the
bottom layer) and `S_remove` removals (condensate, which falls at once in the production configuration).
Assumptions, all explicit:
- `u_evap`: horizontal velocity of the surface exchange that supplies evaporated vapour (required, no default).
- `T_evap`: temperature of evaporated vapour (`nothing` = the bottom-layer temperature, so no sensible term).
- Re-evaporated vapour enters with the layer's own temperature and momentum. Native latent cooling already holds.
- Removed condensate leaves with the layer's temperature and momentum. Native latent heating already holds.
- No latent heat appears here. Native schemes have already applied it.
Zero sources return the input state unchanged (bitwise).
"""
function coupled_step(s::ColumnState, S_add, S_remove, dt, Δσ, g; u_evap, T_evap = nothing, E = 0.0,
                      transport_q = nothing)
    K = length(Δσ)
    if transport_q !== nothing
        m0 = layer_mass(s.ps, Δσ, g)
        abs(sum(m0 .* transport_q)) <= 1e3 * eps(eltype(s.q)) * sum(m0 .* abs.(transport_q)) ||
            throw(ArgumentError("internal transport must conserve column vapour"))
    end
    if all(iszero, S_add) && all(iszero, S_remove) && iszero(E)
        q = transport_q === nothing ? copy(s.q) : s.q .+ dt .* transport_q   # the native transport step, unchanged
        return ColumnState(s.ps, q, copy(s.T), copy(s.u)), (; F = zeros(eltype(S_add), K + 1),
            m = layer_mass(s.ps, Δσ, g), m_new = layer_mass(s.ps, Δσ, g), T_add = copy(s.T), u_add = copy(s.u))
    end
    S = S_add .- S_remove
    m = layer_mass(s.ps, Δσ, g)
    F = interface_flux(S, Δσ)
    ps_new = s.ps + g * sum(S) * dt
    m_new = layer_mass(ps_new, Δσ, g)
    function carry(c, c_add, c_remove)
        mc = m .* c .+ dt .* (S_add .* c_add .- S_remove .* c_remove)
        for k in 1:K
            Ftop, Fbot = F[k], F[k+1]
            k > 1 && (mc[k] += dt * Ftop * upwind(Ftop, c[k-1], c[k]))
            k < K && (mc[k] -= dt * Fbot * upwind(Fbot, c[k], c[k+1]))
        end
        return mc ./ m_new
    end
    # Additions: re-evaporation brings the layer's own T and u, surface evaporation (the E part of the bottom
    # layer's additions) brings T_evap (default: the layer temperature) and zero momentum.
    T_add = copy(s.T); u_add = copy(s.u)
    if E > 0
        Kb = K
        Tsrc = T_evap === nothing ? s.T[Kb] : T_evap
        T_add[Kb] = ((S_add[Kb] - E) * s.T[Kb] + E * Tsrc) / S_add[Kb]
        u_add[Kb] = ((S_add[Kb] - E) * s.u[Kb] + E * u_evap) / S_add[Kb]
    end
    water = ones(eltype(s.q), K)
    q_new = carry(s.q, water, water)
    transport_q === nothing || (q_new .+= dt .* transport_q .* m ./ m_new)   # transport moves vapour, not air
    T_new = carry(s.T, T_add, s.T)
    u_new = carry(s.u, u_add, s.u)
    return ColumnState(ps_new, q_new, T_new, u_new), (; F, m, m_new, T_add, u_add)
end

"The current model: q changes by S/m, p_s, T and u are untouched by the vapour mass."
function uncoupled_step(s::ColumnState, S_add, S_remove, dt, Δσ, g)
    m = layer_mass(s.ps, Δσ, g)
    return ColumnState(s.ps, s.q .+ dt .* (S_add .- S_remove) ./ m, copy(s.T), copy(s.u))
end

"Continuous-time tendencies of the coupled form (what a model parameterization would add)."
function coupled_tendencies(s::ColumnState, S_add, S_remove, Δσ, g; T_add, u_add)
    K = length(Δσ)
    S = S_add .- S_remove
    m = layer_mass(s.ps, Δσ, g)
    F = interface_flux(S, Δσ)
    function rate(c, c_add, c_remove)
        r = S_add .* (c_add .- c) .- S_remove .* (c_remove .- c)        # sources relative to the layer value
        for k in 1:K
            Ftop, Fbot = F[k], F[k+1]
            k > 1 && (r[k] += Ftop * (upwind(Ftop, c[k-1], c[k]) - c[k]))
            k < K && (r[k] -= Fbot * (upwind(Fbot, c[k], c[k+1]) - c[k]))
        end
        return r ./ m
    end
    return (; dlnps = g * sum(S) / s.ps,
            dq = rate(s.q, ones(eltype(s.q), K), ones(eltype(s.q), K)),
            dT = rate(s.T, T_add, s.T),
            du = rate(s.u, u_add, s.u))
end

"""
What a model parameterization adds to the native tendencies for the true sources only. The native q tendency
already contains S/m (and any internal transport), so the q correction is the coupled rate minus S/m. The native
schemes give no mass-related T, u or p_s tendencies, so those corrections are the coupled rates. Internal transport
is not an argument: the correction cannot depend on it.
"""
function mass_coupling_correction(s::ColumnState, S_add, S_remove, Δσ, g; T_add, u_add)
    t = coupled_tendencies(s, S_add, S_remove, Δσ, g; T_add, u_add)
    m = layer_mass(s.ps, Δσ, g)
    return (; dlnps = t.dlnps, dq = t.dq .- (S_add .- S_remove) ./ m, dT = t.dT, du = t.du)
end

"Column budgets: dry mass, vapour, cp-enthalpy (with surface geopotential), momentum."
function budgets(s::ColumnState, Δσ, g; cp = 1004.0, Φs = 0.0)
    m = layer_mass(s.ps, Δσ, g)
    return (; dry = sum(m .* (1 .- s.q)), vapour = sum(m .* s.q), mass = s.ps / g,
            enthalpy = sum(m .* cp .* s.T) + Φs * s.ps / g, momentum = sum(m .* s.u))
end

end
