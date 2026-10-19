# Plume thermodynamics, transport and precipitation from the qualified C08 core.
module EntrainingPlume

include(joinpath(@__DIR__, "convective_detrainment_interface.jl"))
using .ConvectiveDetrainmentInterface
export ConvectiveDetrainmentInterface, PlumeParameters, Thermo,
       plume_profile, Column, saturation_humidity, esat, des_dT

# ---------------------------------------------------------------- thermodynamics
Base.@kwdef struct Thermo{F}
    cp::F = 1004.0
    Lv::F = 2501.0e3
    Lf::F = 3.3e5
    Rd::F = 1000 * 8.3145 / 28.9649
    Rv::F = 1000 * 8.3145 / 18.0153
    g::F = 9.81
    e0::F = 610.78
    T0::F = 273.15
    T_ice::F = 250.16
    law::Symbol = :mixed_phase      # :liquid or :mixed_phase (as ReadyESM)
end

@inline eps_ratio(th) = th.Rd / th.Rv
@inline mu_virt(th) = (1 - eps_ratio(th)) / eps_ratio(th)

@inline e_liquid(T, th) = th.e0 * exp(th.Lv / th.Rv * (1 / th.T0 - 1 / T))
@inline e_ice(T, th) = th.e0 * exp((th.Lv + th.Lf) / th.Rv * (1 / th.T0 - 1 / T))
@inline alpha_blend(T, th) = T >= th.T0 ? one(T) :
    clamp((T - th.T_ice) / (th.T0 - th.T_ice), zero(T), one(T))^2

"""Saturation vapour pressure of the host law (identical to SpeedyWeather / ReadyESM forms)."""
@inline function esat(T, th)
    th.law === :liquid && return e_liquid(T, th)
    T >= th.T0 && return e_liquid(T, th)
    a = alpha_blend(T, th)
    return a * e_liquid(T, th) + (1 - a) * e_ice(T, th)
end

"""Analytic derivative of `esat` with respect to temperature."""
@inline function des_dT(T, th)
    el = e_liquid(T, th); del = el * th.Lv / (th.Rv * T^2)
    (th.law === :liquid || T >= th.T0) && return del
    ei = e_ice(T, th); dei = ei * (th.Lv + th.Lf) / (th.Rv * T^2)
    a = alpha_blend(T, th)
    x = (T - th.T_ice) / (th.T0 - th.T_ice)
    da = (0 < x < 1) ? 2x / (th.T0 - th.T_ice) : zero(T)
    return da * (el - ei) + a * del + (1 - a) * dei
end

"""SpeedyWeather form q* = ε e*/p."""
@inline saturation_humidity(T, p, th) = eps_ratio(th) * esat(T, th) / p
@inline dqs_dT(T, p, th) = eps_ratio(th) * des_dT(T, th) / p

"""Frozen moist static energy of saturated-or-not air with total water qt at T."""
@inline function mse(T, qt, p, Φ, th)
    qs = saturation_humidity(T, p, th)
    qv = min(qt, qs)
    c = qt - qv
    qi = (1 - liquid_fraction(T)) * c
    return th.cp * T + Φ + th.Lv * qv - th.Lf * qi
end

"""
Invert h = mse(T, qt, p, Φ) for T (saturation adjustment). Returns (T, qv, ql, qi).
Newton on a smooth function; at most 30 iterations, |ΔT| < 1e-10 K.
"""
function adjust(h, qt, p, Φ, th)
    # unsaturated guess
    T = (h - Φ - th.Lv * qt) / th.cp
    if qt <= saturation_humidity(T, p, th)
        return T, qt, zero(T), zero(T)
    end
    for _ in 1:30
        qs = saturation_humidity(T, p, th)
        f = liquid_fraction(T)
        c = max(qt - qs, zero(T))
        F = th.cp * T + Φ + th.Lv * min(qt, qs) - th.Lf * (1 - f) * c - h
        # derivative (saturated branch), analytic in every term
        dq = dqs_dT(T, p, th)
        df = dliquid_fraction_dT(T)
        dF = th.cp + th.Lv * dq + th.Lf * (df * c + (1 - f) * dq)
        ΔT = F / dF
        T -= ΔT
        abs(ΔT) < 1e-10 && break
    end
    qs = saturation_humidity(T, p, th)
    qv = min(qt, qs); c = qt - qv; f = liquid_fraction(T)
    return T, qv, f * c, (1 - f) * c
end

@inline virtual_density_temperature(T, qv, c, th) = T * (1 + mu_virt(th) * qv - c)

# ---------------------------------------------------------------- parameters
Base.@kwdef struct PlumeParameters{F}
    entrainment_m::F = 1.0e-3        # ε0 [m-1]
    rh_dependent_entrainment::Bool = true   # ε = ε0 (1.3 - RH) (qs/qs_base)^3, Bechtold et al. 2008
    min_entrainment_m::F = 1.0e-4
    subcloud_entrainment_m::F = 0.0  # entrainment between launch and cloud base
    detrainment_m::F = 0.75e-4       # turbulent δ0 [m-1]
    precip_conversion_m::F = 1.5e-3  # c0 [m-1]
    condensate_threshold::F = 5.0e-4 # kg/kg retained before conversion
    time_scale_s::F = 7200.0         # CAPE relaxation τ
    cape_trigger::F = 25.0           # J/kg (dilute)
    cin_max::F = 150.0               # J/kg allowed below the LFC
    source_depth_Pa::F = 5000.0     # well-mixed source block for the plume
    rain_evaporation::F = 0.5        # fraction of the subsaturation removable per layer
    snow_sublimation::F = 0.5
    melt_below_freezing::Bool = true
    effective_step_factor::F = 2.0   # leapfrog applies tendencies over 2Δt
    test_temperature_change_K::F = 0.02   # size of the closure's linearisation step
    closure::Symbol = :cape           # :cape (v1) | :cape_shallow_blq | :max_cape_blq
    min_parcel_temperature::F = 150.0 # K: below this (or non-finite) the parcel cannot reach the level (v3)
    shallow_depth_Pa::F = 20000.0     # cloud base -> top pressure depth below which a plume is shallow
    blq_min_excess::F = 1000.0        # J/kg floor on (h_u - h_env) at cloud base in the BLQ mass flux
    max_layer_entrainment::F = Inf    # v3c: bound on eps dz per layer (Inf = v3b)
end

"""One column of environment state. Layers 1 (top) ... n (surface)."""
struct Column{V}
    T::V
    qv::V
    qi::V        # existing environmental ice (BS02 tracer), zero if absent
    p::V         # full-level pressure [Pa]
    dp::V        # layer pressure thickness [Pa]
    Φ::V         # full-level geopotential [m2 s-2]
end

# ---------------------------------------------------------------- plume profile
"""
Plume properties for unit cloud-base mass flux. Returns a NamedTuple with
per-layer plume values (Tu, qvu, qlu, qiu, hu, qtu), mass fluxes at interfaces
(Mhalf, length n+1: Mhalf[k] is the upward flux through the interface ABOVE
layer k, Mhalf[n+1] = 0 at the surface), entrainment E, detrainment D and
precipitation formed (Pl, Pi) as mass rates per unit M_b, plus the source-top,
cloud base, LFC and top layer indices, CAPE and CIN.

The plume draws its mass from a well-mixed source block, the lowest layers
within `source_depth_Pa`. Each source layer is entrained in proportion to its
mass, so the plume leaving the block carries the mass-weighted block mean and
nothing depends on how finely the block is resolved. Every layer satisfies the
plume budget exactly:
  (M_in + E) φ_mix = M_in φ_in + E φ_env;   (M_in + E) φ_u = (M_in + E) φ_mix - S;
  M_out + D = M_in + E.
"""
function plume_profile(col::Column, th::Thermo, pp::PlumeParameters)
    n = length(col.T)
    F = eltype(col.T)
    m = col.dp ./ th.g
    hbar = [th.cp * col.T[k] + col.Φ[k] + th.Lv * col.qv[k] - th.Lf * col.qi[k] for k in 1:n]
    qtbar = [col.qv[k] + col.qi[k] for k in 1:n]
    Tu = fill(F(NaN), n); qvu = zeros(F, n); qlu = zeros(F, n); qiu = zeros(F, n)
    hu = zeros(F, n); qtu = zeros(F, n)
    Mhalf = zeros(F, n + 1); E = zeros(F, n); D = zeros(F, n)
    Pl = zeros(F, n); Pi = zeros(F, n); B = zeros(F, n)
    # source block: exactly source_depth_Pa of air from the bottom up. Layers wholly
    # inside it are entrained completely; the top source layer s only by the fraction
    # `wtop` of its mass needed to complete the depth (the whole lowest layer if it
    # is already deeper than source_depth_Pa).
    s = n; depth = col.dp[n]
    while s > 1 && depth < pp.source_depth_Pa
        s -= 1; depth += col.dp[s]
    end
    wtop = s == n ? one(F) : clamp((pp.source_depth_Pa - (depth - col.dp[s])) / col.dp[s], zero(F), one(F))
    weight(j) = j == s ? wtop : one(F)
    msrc = sum(weight(j) * m[j] for j in s:n)
    Min = zero(F); hin = zero(F); qin = zero(F)
    for j in n:-1:s
        Ej = weight(j) * m[j] / msrc
        hu[j] = (Min * hin + Ej * hbar[j]) / (Min + Ej)
        qtu[j] = (Min * qin + Ej * qtbar[j]) / (Min + Ej)
        Tu[j], qvu[j], qlu[j], qiu[j] = adjust(hu[j], qtu[j], col.p[j], col.Φ[j], th)
        E[j] = Ej
        Min += Ej; hin = hu[j]; qin = qtu[j]
        Mhalf[j] = Min
    end
    base = 0; lfc = 0; top = 0; guarded = 0
    cape = zero(F); cin = zero(F)
    qs_base = zero(F)
    # If the mixed source air is already saturated at the top of the block, cloud base
    # is the source top. Its condensate is carried up (no precipitation or detrainment
    # inside the source block), so every layer budget still closes.
    if qlu[s] + qiu[s] > 0
        base = s
        qs_base = saturation_humidity(Tu[s], col.p[s], th)
    end
    k = s
    while k > 1
        k -= 1
        dz = (col.Φ[k] - col.Φ[k + 1]) / th.g
        Min = Mhalf[k + 1]
        ε = if base == 0
            pp.subcloud_entrainment_m
        elseif pp.rh_dependent_entrainment
            qs_env = saturation_humidity(col.T[k], col.p[k], th)
            rh = clamp(col.qv[k] / qs_env, zero(F), one(F))
            max(pp.entrainment_m * (F(1.3) - rh) * (qs_env / qs_base)^3, pp.min_entrainment_m)
        else
            pp.entrainment_m
        end
        δ = base == 0 ? zero(F) : pp.detrainment_m
        # exact within-layer solution of dM/dz = εM and dφ/dz = -ε(φ - φ_env) for a
        # piecewise-constant environment: mixing weight exp(-ε dz) on the incoming plume
        Ek = Min * expm1(min(ε * dz, pp.max_layer_entrainment))   # v3c bound (Inf: v3b bitwise)
        Mmix = Min + Ek
        hmix = (Min * hu[k + 1] + Ek * hbar[k]) / Mmix
        qmix = (Min * qtu[k + 1] + Ek * qtbar[k]) / Mmix
        T, qv, ql, qi = adjust(hmix, qmix, col.p[k], col.Φ[k], th)
        # v3 guard: the mixture cannot exist at this level at a physical temperature
        if !(isfinite(T) && pp.min_parcel_temperature < T < 400 && isfinite(qv) && isfinite(ql) && isfinite(qi))
            if lfc > 0
                D[k + 1] += Mhalf[k + 1]; Mhalf[k + 1] = zero(F); top = k + 1   # terminate at the last physical level
                guarded = 1
            else
                top = 0; lfc = 0                                                # never buoyant: no convection
                guarded = 2
            end
            break
        end
        # precipitation from condensate above the threshold, per unit mass of the mixture
        c = ql + qi
        Plk = Pik = zero(F)
        if c > pp.condensate_threshold
            P = (c - pp.condensate_threshold) * (1 - exp(-pp.precip_conversion_m * dz))
            fl = ql / c
            Plk, Pik = P * fl, P * (1 - fl)
            ql -= Plk; qi -= Pik
            qmix -= P
            hmix += th.Lf * Pik          # removing ice (h content -Lf) at constant T
        end
        Dk = -Mmix * expm1(-δ * dz)            # dM/dz = -δM over the layer
        Tρu = virtual_density_temperature(T, qv, ql + qi, th)
        Tρe = virtual_density_temperature(col.T[k], col.qv[k], col.qi[k], th)
        b = th.g * (Tρu - Tρe) / Tρe
        if base == 0 && c > 0
            base = k
            qs_base = saturation_humidity(T, col.p[k], th)
        end
        if b > 0 && base > 0 && lfc == 0
            lfc = k
        end
        if b <= 0 && lfc > 0
            Dk = Mmix; top = k                      # neutral buoyancy: detrain everything
        elseif b <= 0
            cin -= b * dz
            cin > pp.cin_max && (Dk = Mmix; top = k)
        else
            cape += b * dz
        end
        k == 1 && top == 0 && (Dk = Mmix; top = 1)  # model top reached while buoyant
        Tu[k], qvu[k], qlu[k], qiu[k] = T, qv, ql, qi
        hu[k], qtu[k] = hmix, qmix
        E[k], D[k], B[k] = Ek, Dk, b
        Pl[k], Pi[k] = Plk * Mmix, Pik * Mmix
        Mhalf[k] = Mmix - Dk
        top > 0 && break
    end
    return (; Tu, qvu, qlu, qiu, hu, qtu, Mhalf, E, D, Pl, Pi, B, launch = s, base, lfc, top,
            cape, cin, hbar, qtbar, guarded)
end

# ---------------------------------------------------------------- environment exchange
"""
Environment exchange in detrainment/entrainment form, for every carried variable φ
in (h, qv, ql, qi):
    m_k (φ'_k - φ_k)/τ = D_k φ_u,k - E_k φ_k + M_{k-1/2} φ*_{k-1} - M_{k+1/2} φ*_k,
with φ* = φ (explicit, `τ = nothing`) or φ* = φ' (implicit compensating
subsidence, solved top-down, unconditionally stable). Entrainment is explicit so
the environment loses exactly what the plume took. The column sums telescope, so
water closes to -(precipitation formed) and h to +Lf (ice formed) only if every
layer's plume budget closes: this is a real test of the plume, not a tautology.
`Mb` scales all mass fluxes. Environmental liquid is zero on entry (not carried
by `Column` v1). Returns (x::ColumnExchange before fall-out, m).
"""
function environment_exchange(col::Column, pr, th::Thermo, Mb; τ = nothing)
    n = length(col.T); F = eltype(col.T)
    m = col.dp ./ th.g
    ql0 = zeros(F, n)
    fields = ((pr.hbar, pr.hu), (col.qv, pr.qvu), (ql0, pr.qlu), (col.qi, pr.qiu))
    out = map(fields) do (φ, φu)
        dφ = zeros(F, n)
        prev = zero(F)                      # φ* of the layer above (none above layer 1)
        for k in 1:n
            src = Mb * (pr.D[k] * φu[k] - pr.E[k] * φ[k])
            Mabove = Mb * pr.Mhalf[k]       # subsidence into k from above
            Mbelow = Mb * pr.Mhalf[k + 1]   # subsidence out of k through its base
            if τ === nothing
                dφ[k] = (src + Mabove * prev - Mbelow * φ[k]) / m[k]
                prev = φ[k]
            else
                # Solve directly for the tendency. The equivalent
                # `(new - φ[k]) / τ` form subtracts two nearly equal values.
                # For enthalpy in Float32 that loses most of the small
                # convective increment before the column ledger is formed.
                numerator = src + Mabove * prev - Mbelow * φ[k]
                dφ[k] = numerator / (m[k] + τ * Mbelow)
                prev = φ[k] + τ * dφ[k]
            end
        end
        dφ
    end
    dh, dqv, dql, dqi = out
    dT = (dh .- th.Lv .* dqv .+ th.Lf .* dqi) ./ th.cp     # cp dT = dh - Lv dqv + Lf dqi at fixed Φ
    return ColumnExchange(dqv, dql, dqi, dT, Mb * sum(pr.Pl), Mb * sum(pr.Pi)), m
end

"""Explicit exchange for unit M_b (closure linearisation and diagnostics)."""
unit_exchange(col::Column, pr, th::Thermo) = environment_exchange(col, pr, th, one(eltype(col.T)))

"""
Precipitation fall-out for the scaled exchange `x`. Precipitation formed in layer k
falls into k+1, k+2, ... Below cloud base it may evaporate (rain, Lv) or sublimate
(snow, Lv+Lf): the fraction removed in a layer is 1 - exp(-β (1-RH) Δp / 100 hPa),
capped so the layer cannot pass saturation over the effective step. Snow entering
a layer warmer than freezing melts there (Lf). Water and frozen enthalpy are
conserved exactly.
"""
function fallout(x::ColumnExchange, col::Column, pr, m, Mb, teff, th::Thermo, pp::PlumeParameters)
    n = length(col.T); F = eltype(col.T)
    dqv = copy(x.dqv); dT = copy(x.dT)
    rain = zero(F); snow = zero(F)
    for k in 1:n
        if (rain + snow) > 0
            # layer temperature after the tendencies already applied over the step
            T = col.T[k] + teff * dT[k]
            if pp.melt_below_freezing && T > th.T0 && snow > 0
                # melt at most what the layer's heat above freezing can supply over the step
                melt = min(snow, th.cp * m[k] * (T - th.T0) / (th.Lf * teff))
                dT[k] -= th.Lf * melt / (th.cp * m[k])
                T -= teff * th.Lf * melt / (th.cp * m[k])
                rain += melt; snow -= melt
            end
            if pr.base > 0 && k > pr.base
                qs = saturation_humidity(T, col.p[k], th)
                rh = clamp(col.qv[k] / qs, zero(F), one(F))
                capacity = max(qs - col.qv[k], zero(F)) /
                    (1 + th.Lv / th.cp * dqs_dT(T, col.p[k], th)) * m[k] / teff
                er = min(rain * (1 - exp(-pp.rain_evaporation * (1 - rh) * col.dp[k] / 1.0e4)), capacity)
                es = max(min(snow * (1 - exp(-pp.snow_sublimation * (1 - rh) * col.dp[k] / 1.0e4)), capacity - er), zero(F))
                rain -= er; snow -= es
                dqv[k] += (er + es) / m[k]
                dT[k] -= (th.Lv * er + (th.Lv + th.Lf) * es) / (th.cp * m[k])
            end
        end
        rain += Mb * pr.Pl[k]; snow += Mb * pr.Pi[k]
    end
    return ColumnExchange(dqv, x.dql, x.dqi, dT, rain, snow)
end


end # module
