module FiniteC08Response

const C08=parentmodule(@__MODULE__).C08MultiPlume
const EP=C08.EP
const CDI=EP.ConvectiveDetrainmentInterface

"""Fixed initial plume spectrum and the inherited explicit-entrainment bound."""
function prepare(column,dt,thermo,spectrum)
    F=eltype(column.T)
    isfinite(dt) && dt>0 || throw(ArgumentError("Positive timestep required"))
    all(iszero,column.qi) || throw(ArgumentError("This host carries vapor only"))
    all(isfinite,column.T) && all(>(0),column.T) || throw(ArgumentError("Invalid temperature"))
    all(x->isfinite(x)&&0<=x<=1,column.qv) || throw(ArgumentError("Invalid vapor"))
    parameters=[C08.bin_parameters(spectrum,b) for b in eachindex(spectrum.weights)]
    profiles=[EP.plume_profile(column,thermo,p) for p in parameters]
    active=[C08.bin_triggers(profiles[b],parameters[b]) && spectrum.weights[b]>0
        for b in eachindex(parameters)]
    cape=zero(F)
    for b in eachindex(active)
        active[b] && (cape+=spectrum.weights[b]*profiles[b].cape)
    end
    mass=column.dp./thermo.g
    step=F(spectrum.shared.effective_step_factor*dt)
    composite=any(active) ? C08.composite_profile(profiles,spectrum.weights,active) : nothing
    bound=zero(F)
    if any(active)
        bound=F(Inf)
        for k in eachindex(mass)
            composite.E[k]>0 && (bound=min(bound,F(.9)*mass[k]/(composite.E[k]*step)))
        end
        isfinite(bound) && bound>=0 || error("Invalid entrainment bound")
    end
    return (;column,thermo,spectrum,parameters,profiles,active,cape,mass,step,composite,bound)
end

"""
Evaluate the actual finite implicit exchange, fallout and host re-evaporation.
Reject invalid trial states before CAPE evaluation. No water or temperature
clipping repairs a rejected trial. The caller receives the rejection reason.
"""
function probe(c,mass_flux)
    (;column,thermo,spectrum,step,mass)=c
    F=eltype(column.T)
    M=F(mass_flux)
    isfinite(M) && 0<=M<=c.bound || throw(ArgumentError("Flux outside inherited bound"))
    n=length(column.T)
    empty=CDI.ColumnExchange(n,F)
    invalid(reason,x,T,q)=(valid=false,reason,exchange=x,host_exchange=empty,
        cape=F(NaN),temperature=T,vapor=q,water_relative=F(NaN),enthalpy_relative=F(NaN))
    if isnothing(c.composite)
        return (valid=true,reason=:inactive,exchange=empty,host_exchange=empty,
            cape=zero(F),temperature=copy(column.T),vapor=copy(column.qv),
            water_relative=zero(F),enthalpy_relative=zero(F))
    end
    x,_=EP.environment_exchange(column,c.composite,thermo,M;τ=step)
    T=column.T.+step.*x.dT
    q=column.qv.+step.*x.dqv
    ql=step.*x.dql
    qi=column.qi.+step.*x.dqi
    all(isfinite,T) && all(>(0),T) || return invalid(:exchange_temperature,x,T,q)
    all(a->all(v->isfinite(v)&&v>=0,a),(q,ql,qi)) || return invalid(:exchange_water,x,T,q)
    x=EP.fallout(x,column,c.composite,mass,M,step,thermo,spectrum.shared)
    hostq=x.dqv.+x.dql.+x.dqi
    hostt=x.dT.-(thermo.Lv.*x.dql.+(thermo.Lv+thermo.Lf).*x.dqi)./thermo.cp
    T=column.T.+step.*hostt
    q=column.qv.+step.*hostq
    all(isfinite,T) && all(>(0),T) || return invalid(:host_temperature,x,T,q)
    all(v->isfinite(v)&&0<=v<=1,q) || return invalid(:host_water,x,T,q)
    all(v->isfinite(v)&&v>=0,(x.rain,x.snow)) || return invalid(:precipitation,x,T,q)
    host=CDI.ColumnExchange(hostq,zeros(F,n),zeros(F,n),hostt,x.rain,x.snow)
    cape=zero(F)
    for b in eachindex(c.active)
        c.active[b] || continue
        cape+=spectrum.weights[b]*C08.dilute_cape_perturbed_scalar(
            column,thermo,c.parameters[b],host,step)
    end
    isfinite(cape) && cape>=0 || return invalid(:cape,x,T,q)
    w=CDI.check_water(host,mass)
    e=CDI.check_enthalpy(host,mass;cp=thermo.cp,Lv=thermo.Lv,Lf=thermo.Lf)
    return (valid=true,reason=:valid,exchange=x,host_exchange=host,cape,temperature=T,vapor=q,
        water_relative=abs(w.residual)/max(w.gross,floatmin(F)),
        enthalpy_relative=abs(e.residual)/max(e.gross,floatmin(F)))
end

end
