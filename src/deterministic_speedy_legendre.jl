"""
Deterministic CUDA implementation of SpeedyTransforms' forward Legendre
transform.

SpeedyTransforms 0.1.5 assigns one work item to every `(layer, latitude,
order)` and atomically accumulates its contribution into the packed spectral
coefficients. The atomics are data-race free, but their arrival order changes
between executions and therefore changes Float32 roundoff. ReadyESM instead
assigns one work item to each `(layer, coefficient)`, sums northern/southern
latitude-pair contributions in ascending latitude order, and writes the result
once. This is the grid-to-spectral direction only; the inverse Legendre
transform already gives each output element a single writer.
"""

const _DETERMINISTIC_LEGENDRE_MMAX_CACHE = IdDict{Any, Any}()

function _deterministic_legendre_mmax(S)
    return get!(_DETERMINISTIC_LEGENDRE_MMAX_CACHE, S.kjm_indices) do
        CUDA.CuArray(S.mmax_truncation)
    end
end

KernelAbstractions.@kernel function _deterministic_forward_legendre_kernel!(
    spectral_coefficients,
    legendre_polynomials,
    fourier_north,
    fourier_south,
    spectral_l_indices,
    spectral_m_indices,
    longitude_offsets,
    solid_angles,
    maximum_order_per_latitude,
    nharmonics,
    nlat_half,
)
    linear_index = @index(Global, Linear)
    harmonic = mod1(linear_index, nharmonics)
    layer = (linear_index - 1) ÷ nharmonics + 1

    degree = spectral_l_indices[harmonic]
    order = spectral_m_indices[harmonic]
    symmetric_about_equator = iseven(degree - order)
    coefficient = zero(eltype(spectral_coefficients))

    @inbounds for latitude in 1:nlat_half
        # `maximum_order_per_latitude` is zero based, while `order` is the
        # one-based Fourier/spectral column index.
        if order <= maximum_order_per_latitude[latitude] + 1
            northern = fourier_north[order, layer, latitude]
            southern = fourier_south[order, layer, latitude]
            paired = ifelse(
                symmetric_about_equator,
                northern + southern,
                northern - southern,
            )
            rotated_solid_angle =
                solid_angles[latitude] * conj(longitude_offsets[order, latitude])
            coefficient += legendre_polynomials[harmonic, latitude] *
                           (rotated_solid_angle * paired)
        end
    end

    @inbounds spectral_coefficients[harmonic, layer] = coefficient
end

function SpeedyWeather.SpeedyTransforms._legendre!(
    coefficients::SpeedyWeather.LowerTriangularArray,
    fourier_north::CUDA.CuArray{<:Complex, 3},
    fourier_south::CUDA.CuArray{<:Complex, 3},
    ::SpeedyWeather.SpeedyTransforms.ColumnScratchMemory,
    transform::SpeedyWeather.SpeedyTransforms.SpectralTransform{
        NF,
        <:SpeedyWeather.Architectures.GPU,
    },
) where {NF}
    SpeedyWeather.SpeedyTransforms.ismatching(transform, coefficients) ||
        throw(DimensionMismatch(transform, coefficients))
    nlayers = size(coefficients, 2)
    expected_scratch_size = (
        transform.nfreq_max,
        transform.nlayers,
        transform.grid.nlat_half,
    )
    size(fourier_north) == size(fourier_south) == expected_scratch_size ||
        throw(DimensionMismatch(transform, coefficients))
    nlayers <= transform.nlayers ||
        throw(DimensionMismatch(transform, coefficients))

    nharmonics = size(coefficients.data, 1)
    maximum_order_per_latitude = _deterministic_legendre_mmax(transform)
    backend = KernelAbstractions.get_backend(fourier_north)
    _deterministic_forward_legendre_kernel!(backend)(
        coefficients.data,
        transform.legendre_polynomials.data,
        fourier_north,
        fourier_south,
        transform.spectrum.l_indices,
        transform.spectrum.m_indices,
        transform.lon_offsets,
        transform.solid_angles,
        maximum_order_per_latitude,
        nharmonics,
        transform.grid.nlat_half;
        ndrange = nharmonics * nlayers,
    )
    return nothing
end
