module DeterministicMultiPlumeC08

include(joinpath(@__DIR__, "entraining_plume.jl"))
using .EntrainingPlume
using .EntrainingPlume.ConvectiveDetrainmentInterface
const EP = EntrainingPlume
const MIN_MASS_FLUX = 1.0e-12
const DEFAULT_CLOSURE_ITERATIONS = 8      # selected by the frozen cost gate (cost_gate_contract.toml, cost_gate_result.toml)

struct PlumeSpectrum{F}
    entrainment_m::Vector{F}
    weights::Vector{F}
    shared::PlumeParameters{F}
    closure_iterations::Int
end
PlumeSpectrum{F}(e, w, p) where {F} = PlumeSpectrum{F}(e, w, p, DEFAULT_CLOSURE_ITERATIONS)
function PlumeSpectrum(; entrainment_m = [0.3e-3, 1.0e-3, 3.0e-3], weights = fill(1 / 3, 3), max_layer_entrainment = 2.0,
                       closure_iterations = DEFAULT_CLOSURE_ITERATIONS, kwargs...)
    length(entrainment_m) == length(weights) || throw(ArgumentError("bins and weights differ in length"))
    all(>=(0), weights) && sum(weights) <= 1 + 1e-12 || throw(ArgumentError("weights must be non-negative, sum <= 1"))
    shared = PlumeParameters(; max_layer_entrainment, kwargs...)
    F = typeof(shared.entrainment_m)
    return PlumeSpectrum{F}(F.(entrainment_m), F.(weights), shared, closure_iterations)
end
function bin_parameters(s::PlumeSpectrum, i)
    names = fieldnames(PlumeParameters)
    nt = NamedTuple{names}(Tuple(getfield(s.shared, f) for f in names))
    return PlumeParameters(; merge(nt, (; entrainment_m = s.entrainment_m[i]))...)
end
bin_triggers(pr, pp) = pr.lfc > 0 && pr.top > 0 && isfinite(pr.cape) && isfinite(pr.cin) &&
                       pr.cape >= pp.cape_trigger && pr.cin <= pp.cin_max

function composite_profile(prs, w, active)
    if count(active) == 1 && w[findfirst(active)] == 1      # one bin: its own profile, unchanged
        return prs[findfirst(active)]
    end
    n = length(prs[1].E); F = eltype(prs[1].E)
    E = zeros(F, n); D = zeros(F, n); Mhalf = zeros(F, n + 1); Pl = zeros(F, n); Pi = zeros(F, n)
    Dh = zeros(F, n); Dqv = zeros(F, n); Dql = zeros(F, n); Dqi = zeros(F, n)
    base = 0; top = n + 1
    for (i, pr) in enumerate(prs)
        active[i] || continue
        wi = w[i]
        E .+= wi .* pr.E; D .+= wi .* pr.D; Mhalf .+= wi .* pr.Mhalf; Pl .+= wi .* pr.Pl; Pi .+= wi .* pr.Pi
        Dh .+= wi .* pr.D .* pr.hu; Dqv .+= wi .* pr.D .* pr.qvu; Dql .+= wi .* pr.D .* pr.qlu; Dqi .+= wi .* pr.D .* pr.qiu
        base = max(base, pr.base); top = min(top, pr.top)
    end
    div(a) = [D[k] > 0 ? a[k] / D[k] : zero(F) for k in 1:n]
    return (; hbar = prs[1].hbar, hu = div(Dh), qvu = div(Dqv), qlu = div(Dql), qiu = div(Dqi), E, D, Mhalf, Pl, Pi, base, top)
end

include("scalar_cape_workspace.jl")
include("typed_constructors.jl")
include("float32_repairs.jl")
include("source_order_repair.jl")

export Column, Thermo, PlumeSpectrum, typed_column, typed_thermo,
    typed_plume_parameters, typed_plume_spectrum
end
