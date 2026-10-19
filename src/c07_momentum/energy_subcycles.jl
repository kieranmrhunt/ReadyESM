module C07EnergySubcycles
include("energy_return.jl")
const E=C07EnergyReturn

"""Keep the qualified momentum result and transport dissipated heat between substeps."""
function exchange(reference, u,v,mass,M,en,det,mass_flux,dt;
        strength=one(mass_flux),maximum_subcycles=16,cp=1004)
    isfinite(cp) && cp>0 || throw(ArgumentError("positive heat capacity required"))
    old=reference(u,v,mass,M,en,det,mass_flux,dt;strength,maximum_subcycles)
    F=eltype(old.du)
    n=length(u)
    if old.skipped || old.subcycles==0
        return merge(old,(;dT=zeros(F,n)))
    end
    full_step=F(dt)
    if old.subcycles==1
        x=E.exchange(u,v,mass,M,en,det,mass_flux,full_step;strength,cp)
        x.du==old.du && x.dv==old.dv || error("single-step momentum mismatch")
        return merge(old,(;dT=x.dT))
    end
    step=full_step/F(old.subcycles)
    uw,vw=F.(u),F.(v)
    heat=zeros(F,n)
    for _ in 1:old.subcycles
        x=E.exchange(uw,vw,mass,M,en,det,mass_flux,step;strength,cp)
        anomaly=heat .- heat[end]
        plume=E.C._plume_momentum(anomaly,F.(M),F.(en),F.(det))
        dh=E.C._environment_tendency(anomaly,plume,F.(mass),F.(M),F.(en),F.(det),
            F(mass_flux)*F(strength),step)
        heat .+= step .* (dh .+ F(cp).*x.dT)
        uw .+= step.*x.du
        vw .+= step.*x.dv
    end
    uw==old.final_u && vw==old.final_v || error("subcycle trajectory mismatch")
    return merge(old,(;dT=heat./(F(cp)*full_step)))
end
end
