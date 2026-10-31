module SeasonalVegetation

using Dates

"Twelve prescribed monthly nodes on the model's native land columns."
struct MonthlyLandCycle{A}
    cover::A
    area_leaf_index::A
    snow_free_albedo::A
end

function MonthlyLandCycle(cover::AbstractMatrix, leaf_index::AbstractMatrix,
                         snow_free_albedo::AbstractMatrix)
    size(cover) == size(leaf_index) == size(snow_free_albedo) ||
        throw(DimensionMismatch("monthly land fields must share column/month dimensions"))
    size(cover,1) > 0 && size(cover,2) == 12 ||
        throw(ArgumentError("monthly land cycle requires nonempty columns and twelve months"))
    all(x->isfinite(x) && 0 <= x <= 1,cover) || throw(ArgumentError("invalid vegetation cover"))
    all(x->isfinite(x) && 0 <= x <= 10,leaf_index) || throw(ArgumentError("invalid leaf area index"))
    all(x->isfinite(x) && 0 <= x <= 1,snow_free_albedo) || throw(ArgumentError("invalid snow-free albedo"))
    # Own the input data. Preserve leaf area per land-column area through
    # interpolation rather than averaging a conditional LAI independently.
    c,l,a = Matrix{Float64}(cover),Matrix{Float64}(leaf_index),Matrix{Float64}(snow_free_albedo)
    return MonthlyLandCycle{Matrix{Float64}}(c,c.*l,a)
end

month_centre(date::DateTime) = begin
    start = DateTime(year(date),month(date),1)
    stop = start + Month(1)
    start + Millisecond(div(Dates.value(stop-start),2))
end

"Gregorian month-centred interpolation, including December/January and leap years."
function bracket(date::DateTime)
    centre = month_centre(date)
    if date < centre
        left,right = month_centre(date-Month(1)),centre
    else
        left,right = centre,month_centre(date+Month(1))
    end
    weight = Dates.value(date-left)/Dates.value(right-left)
    return (; lower=month(left),upper=month(right),weight,left,right)
end

function sample!(cover,leaf_index,snow_free_albedo,cycle::MonthlyLandCycle,date::DateTime)
    n=size(cycle.cover,1)
    all(x->size(x)==(n,), (cover,leaf_index,snow_free_albedo)) ||
        throw(DimensionMismatch("destination fields must be native land-column vectors"))
    all(x->eltype(x)<:AbstractFloat,(cover,leaf_index,snow_free_albedo)) ||
        throw(ArgumentError("destination fields must use floating-point storage"))
    b=bracket(date)
    for i in 1:n
        f=(1-b.weight)*cycle.cover[i,b.lower]+b.weight*cycle.cover[i,b.upper]
        area_lai=(1-b.weight)*cycle.area_leaf_index[i,b.lower]+b.weight*cycle.area_leaf_index[i,b.upper]
        cover[i]=f
        leaf_index[i]=iszero(f) ? 0 : area_lai/f
        snow_free_albedo[i]=(1-b.weight)*cycle.snow_free_albedo[i,b.lower]+b.weight*cycle.snow_free_albedo[i,b.upper]
    end
    return b
end

function sample_elapsed!(cover,leaf_index,albedo,cycle,start::DateTime,seconds::Real)
    isfinite(seconds) && seconds >= 0 || throw(ArgumentError("invalid elapsed model time"))
    return sample!(cover,leaf_index,albedo,cycle,start+Millisecond(round(Int,1000*seconds)))
end

"One shared exchange value for land and radiation, held until their next radiation update."
function shared_albedo(::Type{FT},snow_free,snow_cover,snow_albedo) where {FT<:AbstractFloat}
    all(x->isfinite(x) && 0 <= x <= 1,(snow_free,snow_cover,snow_albedo)) ||
        throw(ArgumentError("albedo inputs must be finite fractions"))
    return FT((1-snow_cover)*snow_free+snow_cover*snow_albedo)
end

end
