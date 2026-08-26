"""
    DeterministicMaxRandomOverlap(counters, seed)

Maximum-random McICA cloud overlap with an explicit per-column sample counter.

RRTMGP's `MaxRandomOverlap` calls the implicit device `Random.rand()` from
inside its one-thread-per-column CUDA kernels. CUDA.jl stores that RNG counter
once per warp, so columns whose cloud profiles require different numbers of
draws race on shared state. This sampler instead hashes the explicit
`(seed, column, spectral sample, draw)` tuple. It is independent of warp
scheduling, allocation-free in the radiation kernel, and resettable before
each solve.
"""
struct DeterministicMaxRandomOverlap{C, S} <:
       RRTMGP.AtmosphericStates.AbstractCloudMask
    counters::C
    seed::S
end

Adapt.@adapt_structure DeterministicMaxRandomOverlap

function DeterministicMaxRandomOverlap(DA, ncol::Integer)
    ncol > 0 || throw(ArgumentError("ncol must be positive"))
    return DeterministicMaxRandomOverlap(
        DA{UInt32}(zeros(UInt32, ncol)),
        DA{UInt64}(zeros(UInt64, 1)),
    )
end

"""Reset every per-column McICA counter and install the solve-specific seed."""
function _reset_deterministic_cloud_sampler!(
    sampler::DeterministicMaxRandomOverlap,
    seed::Integer,
)
    fill!(sampler.counters, zero(eltype(sampler.counters)))
    fill!(sampler.seed, UInt64(seed))
    return nothing
end

_reset_deterministic_cloud_sampler!(::Any, seed::Integer) = nothing

function _reset_solver_cloud_sampler!(solver, seed::Integer)
    cloud_state = solver.as.cloud_state
    isnothing(cloud_state) && return nothing
    _reset_deterministic_cloud_sampler!(cloud_state.mask_type, seed)
    return nothing
end

@inline function _splitmix64(value::UInt64)
    value += UInt64(0x9e3779b97f4a7c15)
    value = (value ⊻ (value >> 30)) * UInt64(0xbf58476d1ce4e5b9)
    value = (value ⊻ (value >> 27)) * UInt64(0x94d049bb133111eb)
    return value ⊻ (value >> 31)
end

@inline function _cloud_uniform(
    ::Type{FT},
    seed::UInt64,
    gcol::Integer,
    sample::UInt32,
    draw::UInt32,
) where {FT}
    key = seed ⊻
          UInt64(gcol) * UInt64(0xd2b74407b1ce6e93) ⊻
          UInt64(sample) * UInt64(0xca5a826395121157) ⊻
          UInt64(draw) * UInt64(0x9e3779b97f4a7c15)
    # Twenty-four high-quality mantissa bits are ample for cloud-fraction
    # thresholding and give identical sampling decisions for Float32/Float64.
    return FT(_splitmix64(key) >> 40) * FT(0x1p-24)
end

@inline _cloud_column_index(values::SubArray{<:Any, 1}) =
    Int(last(parentindices(values)))
@inline _cloud_column_index(::AbstractVector) = 1

@inline function _first_cloudy_layer(cloud_fraction)
    @inbounds for layer in eachindex(cloud_fraction)
        cloud_fraction[layer] > 0 && return layer
    end
    return 0
end

@inline function _last_cloudy_layer(cloud_fraction)
    @inbounds for layer in reverse(eachindex(cloud_fraction))
        cloud_fraction[layer] > 0 && return layer
    end
    return 0
end

"""
Deterministic, GPU-safe implementation of RRTMGP's maximum-random overlap.

The sampled distribution and vertical-overlap rule match RRTMGP's
`MaxRandomOverlap`; only the source of uniform variates changes.
"""
function RRTMGP.Optics.build_cloud_mask!(
    cloud_mask::AbstractVector{Bool},
    cloud_fraction::AbstractVector{FT},
    sampler::DeterministicMaxRandomOverlap,
) where {FT}
    gcol = _cloud_column_index(cloud_fraction)
    @inbounds sample = sampler.counters[gcol]
    @inbounds sampler.counters[gcol] = sample + one(sample)
    @inbounds seed = sampler.seed[1]

    first_cloudy = _first_cloudy_layer(cloud_fraction)
    if first_cloudy == 0
        @inbounds for layer in eachindex(cloud_mask)
            cloud_mask[layer] = false
        end
        return nothing
    end

    last_cloudy = _last_cloudy_layer(cloud_fraction)
    @inbounds for layer in 1:(first_cloudy - 1)
        cloud_mask[layer] = false
    end
    @inbounds for layer in (last_cloudy + 1):length(cloud_mask)
        cloud_mask[layer] = false
    end

    draw = UInt32(0)
    @inbounds fraction_above = cloud_fraction[last_cloudy]
    random_above = _cloud_uniform(FT, seed, gcol, sample, draw)
    draw += one(draw)
    @inbounds cloud_mask[last_cloudy] =
        cloudy_above = random_above >= (one(FT) - fraction_above)

    layer = last_cloudy - 1
    while layer >= first_cloudy
        @inbounds fraction = cloud_fraction[layer]
        if fraction > zero(FT)
            random = if cloudy_above
                random_above
            else
                value = _cloud_uniform(FT, seed, gcol, sample, draw) *
                        (one(FT) - fraction_above)
                draw += one(draw)
                value
            end
            cloudy = random >= (one(FT) - fraction)
            random_above = random
        else
            cloudy = false
        end
        @inbounds cloud_mask[layer] = cloudy
        fraction_above = fraction
        cloudy_above = cloudy
        layer -= 1
    end
    return nothing
end
