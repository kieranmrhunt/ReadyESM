"""Construct a `Thermo` whose floating fields all use `FT`."""
function typed_thermo(::Type{FT}; kwargs...) where {FT <: AbstractFloat}
    defaults = Thermo()
    overrides = (; kwargs...)
    names = fieldnames(typeof(defaults))
    values = ntuple(length(names)) do index
        name = names[index]
        value = hasproperty(overrides, name) ? getproperty(overrides, name) :
            getfield(defaults, name)
        fieldtype(Thermo{FT}, index) <: AbstractFloat ? FT(value) : value
    end
    return Thermo{FT}(values...)
end

"""Construct plume parameters without mixed Float32 and Float64 defaults."""
function typed_plume_parameters(::Type{FT}; kwargs...) where {
    FT <: AbstractFloat,
}
    defaults = PlumeParameters()
    overrides = (; kwargs...)
    names = fieldnames(typeof(defaults))
    values = ntuple(length(names)) do index
        name = names[index]
        value = hasproperty(overrides, name) ? getproperty(overrides, name) :
            getfield(defaults, name)
        fieldtype(PlumeParameters{FT}, index) <: AbstractFloat ? FT(value) :
            value
    end
    return PlumeParameters{FT}(values...)
end

"""Construct a plume spectrum whose arrays and shared parameters use `FT`."""
function typed_plume_spectrum(
    ::Type{FT};
    entrainment_m = (0.3e-3, 1.0e-3, 3.0e-3),
    weights = (1 / 3, 1 / 3, 1 / 3),
    closure_iterations = DEFAULT_CLOSURE_ITERATIONS,
    kwargs...,
) where {FT <: AbstractFloat}
    length(entrainment_m) == length(weights) ||
        throw(ArgumentError("bins and weights differ in length"))
    typed_entrainment = FT.(collect(entrainment_m))
    typed_weights = FT.(collect(weights))
    all(>=(zero(FT)), typed_weights) ||
        throw(ArgumentError("weights must be non-negative"))
    sum(typed_weights) <= one(FT) + FT(16) * eps(FT) ||
        throw(ArgumentError("weights must sum to at most one"))
    shared = typed_plume_parameters(FT; kwargs...)
    return PlumeSpectrum{FT}(
        typed_entrainment,
        typed_weights,
        shared,
        Int(closure_iterations),
    )
end

"""Convert a plume test column to one floating-point precision."""
function typed_column(::Type{FT}, column::Column) where {FT <: AbstractFloat}
    return Column(
        FT.(column.T),
        FT.(column.qv),
        FT.(column.qi),
        FT.(column.p),
        FT.(column.dp),
        FT.(column.Φ),
    )
end
