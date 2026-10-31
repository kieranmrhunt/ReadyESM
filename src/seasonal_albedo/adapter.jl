module SeasonalAlbedoAdapter

export visible_weights, pack_albedos!, equivalent_land_albedo

unit(x)=isfinite(x)&&0<=x<=1

"""Visible weights by band. The crossing-band choice must be explicit."""
function visible_weights(bounds;crossing,threshold=14286.)
    size(bounds,1)==2 || throw(ArgumentError("expected two band edges"))
    isfinite(threshold)&&threshold>0 || throw(ArgumentError("invalid visible threshold"))
    crossing in (:equal_halves,:visible,:near_infrared) || throw(ArgumentError("unknown crossing-band rule"))
    weights=zeros(Float64,size(bounds,2))
    for band in eachindex(weights)
        lo,hi=bounds[:,band]
        isfinite(lo)&&isfinite(hi)&&0<=lo<hi || throw(ArgumentError("invalid band edges"))
        weights[band]=lo>=threshold ? 1. : hi<=threshold ? 0. :
            crossing==:equal_halves ? .5 : crossing==:visible ? 1. : 0.
    end
    return weights
end

"""Pack separate direct/diffuse spectral land and grey ocean albedos."""
function pack_albedos!(direct,diffuse,weights,fraction,ocean,vis_direct,nir_direct,vis_diffuse,nir_diffuse)
    nband,ncol=size(direct)
    size(diffuse)==size(direct) || throw(DimensionMismatch("direct/diffuse shape"))
    length(weights)==nband || throw(DimensionMismatch("band weights"))
    all(x->length(x)==ncol,(fraction,ocean,vis_direct,nir_direct,vis_diffuse,nir_diffuse)) ||
        throw(DimensionMismatch("surface columns"))
    all(unit,weights) || throw(ArgumentError("invalid band weights"))
    for col in 1:ncol
        f=fraction[col]
        unit(f) || throw(ArgumentError("invalid land fraction"))
        o=f==1 ? 0. : ocean[col]
        unit(o) || throw(ArgumentError("invalid active ocean albedo"))
        if f==0
            direct[:,col].=o
            diffuse[:,col].=o
            continue
        end
        vd,nd,vf,nf=vis_direct[col],nir_direct[col],vis_diffuse[col],nir_diffuse[col]
        all(unit,(vd,nd,vf,nf)) || throw(ArgumentError("invalid active land albedo"))
        for band in 1:nband
            w=weights[band]
            direct[band,col]=(1-f)*o+f*(w*vd+(1-w)*nd)
            diffuse[band,col]=(1-f)*o+f*(w*vf+(1-w)*nf)
        end
    end
    return nothing
end

"""
Recover the land albedo matching a solved mixed-cell broadband reflection.

Assumes grey ocean albedo, a shared downwelling flux and linear area mixing.
Return a value rounded once to exchange precision plus the weighted residual.
The tolerance bounds arithmetic error, not a climate-energy tuning allowance.
"""
function equivalent_land_albedo(::Type{FT},fraction,ocean,down,up;previous) where {FT<:AbstractFloat}
    FT in (Float32,Float64) || throw(ArgumentError("unsupported exchange precision"))
    unit(fraction)&&unit(previous) || throw(ArgumentError("invalid fraction or held albedo"))
    isfinite(down)&&down>=0&&isfinite(up)&&up>=0 || throw(ArgumentError("invalid shortwave flux"))
    f,d,u=Float64(fraction),Float64(down),Float64(up)
    o=f==1 ? 0. : Float64(ocean)
    unit(o) || throw(ArgumentError("invalid active ocean albedo"))
    tolerance=64*Float64(eps(FT))*max(d,u,1.)
    ocean_up=(1-f)*o*d
    available=f*d
    target=u-ocean_up
    -tolerance<=target<=available+tolerance || throw(ArgumentError("mixed reflection is inconsistent with bounded land albedo"))
    old=FT(previous)
    weak_signal=available<=tolerance
    use_previous=iszero(available)||(weak_signal&&abs(ocean_up+available*old-u)<=tolerance)
    albedo=use_previous ? old : FT(clamp(target/available,0.,1.))
    residual=ocean_up+available*Float64(albedo)-u
    abs(residual)<=tolerance || throw(ArgumentError("land shortwave closure exceeds arithmetic bound"))
    return (;albedo,residual,tolerance,weak_signal)
end

end
