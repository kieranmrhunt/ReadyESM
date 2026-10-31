module SeasonalLandRemapping

const FIELD_NAMES = (:cvh,:cvl,:lai_hv,:lai_lv,:aluvp,:aluvd,:alnip,:alnid,:lsm)
const ALBEDO_NAMES = (:aluvp,:aluvd,:alnip,:alnid)

"Validated global regular-grid monthly nodes, stored latitude × longitude × month."
struct Source
    latitude::Vector{Float64}
    longitude::Vector{Float64}
    fields::Dict{Symbol,Array{Float64,3}}
    identity::String
end

function source(latitude,longitude,fields;identity="synthetic")
    lat,lon=Float64.(latitude),Float64.(longitude)
    length(lat)>=3 && length(lon)>=4 || throw(ArgumentError("source grid too small"))
    all(isfinite,lat) && all(isfinite,lon) || throw(ArgumentError("nonfinite coordinates"))
    lat[1]==-90 && lat[end]==90 || throw(ArgumentError("latitude must increase pole to pole"))
    lon[1]==0 || throw(ArgumentError("longitude must start at zero"))
    dy=180/(length(lat)-1); dx=360/length(lon)
    all(isapprox.(diff(lat),dy;atol=1e-10,rtol=0)) || throw(ArgumentError("irregular latitude"))
    all(isapprox.(diff(lon),dx;atol=1e-10,rtol=0)) || throw(ArgumentError("irregular or duplicate-seam longitude"))
    Set(keys(fields))==Set(FIELD_NAMES) || throw(ArgumentError("source fields differ"))
    owned=Dict{Symbol,Array{Float64,3}}()
    for name in FIELD_NAMES
        values=Float64.(fields[name])
        size(values)==(length(lat),length(lon),12) || throw(DimensionMismatch("wrong dimensions for $name"))
        upper=name in (:lai_hv,:lai_lv) ? 10 : 1
        all(x->!isfinite(x) || 0<=x<=upper,values) || throw(ArgumentError("out of range $name"))
        owned[name]=copy(values)
    end
    all(isfinite,owned[:lsm]) || throw(ArgumentError("missing source land mask"))
    land=owned[:lsm].>0
    for name in FIELD_NAMES
        all(isfinite,owned[name][land]) || throw(ArgumentError("missing $name over source land"))
    end
    # High/low cover partitions must describe one land area. Do not silently
    # clamp their sum and leave the associated leaf area unscaled.
    total=owned[:cvh].+owned[:cvl]
    all(x->x<=1+8eps(Float32),total[land]) || throw(ArgumentError("overlapping high/low vegetation cover"))
    return Source(lat,lon,owned,String(identity))
end

function stencil(s::Source,longitude,latitude)
    isfinite(longitude) && isfinite(latitude) && -90<=latitude<=90 ||
        throw(ArgumentError("invalid target location"))
    x=mod(Float64(longitude),360)/(360/length(s.longitude))
    y=(Float64(latitude)+90)/(180/(length(s.latitude)-1))
    i0=floor(Int,x)+1; i1=mod1(i0+1,length(s.longitude))
    j0=min(floor(Int,y)+1,length(s.latitude)-1); j1=j0+1
    wx=x-floor(x); wy=y-(j0-1)
    return ((j0,i0),(j0,i1),(j1,i0),(j1,i1)),
           ((1-wy)*(1-wx),(1-wy)*wx,wy*(1-wx),wy*wx)
end

"Land-conditioned weights. A bounded nearest-land fallback is reported explicitly."
function weights(s::Source,longitude,latitude,month;max_distance=12.0)
    1<=month<=12 || throw(ArgumentError("invalid month"))
    isfinite(max_distance) && max_distance>=0 || throw(ArgumentError("invalid search radius"))
    indices,geometric=stencil(s,longitude,latitude)
    w=[g*s.fields[:lsm][j,i,month] for ((j,i),g) in zip(indices,geometric)]
    total=sum(w)
    if total>0
        return (;indices=collect(indices),weights=w./total,fallback=false,distance=0.0)
    end
    # Match the accepted convention's bounded angular search, including the
    # periodic seam. Require land, not a convenient finite ocean placeholder.
    best=nothing; best_distance=Inf
    for j in eachindex(s.latitude)
        dy=abs(s.latitude[j]-latitude)
        dy<=max_distance || continue
        for i in eachindex(s.longitude)
            s.fields[:lsm][j,i,month]>1e-6 || continue
            dx=abs(mod(s.longitude[i]-longitude+180,360)-180)
            distance=hypot(dx,dy)
            if distance<=max_distance && distance<best_distance
                best=(j,i);best_distance=distance
            end
        end
    end
    isnothing(best) && error("no source land within $max_distance angular degrees of ($longitude,$latitude)")
    return (;indices=[best],weights=[1.0],fallback=true,distance=best_distance)
end

"Keep roundoff from taking a physical cover fraction outside [0, 1]."
function bounded_cover(c::Real)
    tolerance=16eps(Float64)
    isfinite(c) && -tolerance<=c<=1+tolerance ||
        throw(ArgumentError("remapped vegetation cover is outside [0, 1] beyond roundoff"))
    return clamp(c,0.0,1.0)
end

function remap(s::Source,longitudes,latitudes;max_distance=12.0)
    length(longitudes)==length(latitudes)>0 || throw(DimensionMismatch("target coordinate lengths"))
    n=length(longitudes)
    cover=zeros(n,12); area_leaf_index=zeros(n,12)
    albedo=Dict(name=>zeros(n,12) for name in ALBEDO_NAMES)
    fallback=falses(n,12); distance=zeros(n,12)
    for m in 1:12,k in 1:n
        mapping=weights(s,longitudes[k],latitudes[k],m;max_distance)
        fallback[k,m]=mapping.fallback;distance[k,m]=mapping.distance
        for ((j,i),w) in zip(mapping.indices,mapping.weights)
            w>0 || continue # avoid 0*NaN in absent ocean source components
            high,low=s.fields[:cvh][j,i,m],s.fields[:cvl][j,i,m]
            total=high+low
            # Only roundoff-level overlap was allowed at input. Normalize
            # both cover and leaf area consistently in that case.
            scale=total>1 ? 1/total : 1.0
            cover[k,m]+=w*scale*total
            area_leaf_index[k,m]+=w*scale*(high*s.fields[:lai_hv][j,i,m]+low*s.fields[:lai_lv][j,i,m])
            for name in ALBEDO_NAMES
                albedo[name][k,m]+=w*s.fields[name][j,i,m]
            end
        end
        cover[k,m]=bounded_cover(cover[k,m])
    end
    leaf_index=map((c,a)->iszero(c) ? 0.0 : a/c,cover,area_leaf_index)
    return (;cover,area_leaf_index,leaf_index,albedo,fallback,distance,source_identity=s.identity)
end

end
