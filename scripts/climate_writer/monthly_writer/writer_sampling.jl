module WriterSampling

export weighted_mean, unique_surface_mask

"""Host reduction mask. Exclude the redundant half-row of RightCenterFolded storage."""
function unique_surface_mask(active::AbstractMatrix; right_center_folded=false)
    all(v -> v == 0 || v == 1, active) || error("surface mask must be binary")
    mask=Float64.(active)
    if right_center_folded
        nx,ny=size(mask)
        iseven(nx) && nx > 0 && ny > 0 || error("invalid center-folded mask shape")
        mask[nx÷2+1:nx,ny] .= 0
    end
    return mask
end

"""Area-weighted mean over positive-area cells. Ignore absent cells, including NaN dry cells."""
function weighted_mean(values, weights)
    size(values) == size(weights) || throw(DimensionMismatch("field and area weights differ"))
    numerator = 0.0
    denominator = 0.0
    for i in eachindex(values, weights)
        w = Float64(weights[i])
        isfinite(w) && w >= 0 || error("invalid area weight at $i")
        iszero(w) && continue
        x = Float64(values[i])
        isfinite(x) || error("nonfinite value on an active cell at $i")
        numerator += w * x
        denominator += w
    end
    denominator > 0 || error("no active area in weighted mean")
    return numerator / denominator
end

end
