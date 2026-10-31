module MonthlyCore
using Dates

mutable struct Accumulator
    start::DateTime
    stop::DateTime
    step_ms::Int
    sigma::Vector{Float64}
    dsigma::Vector{Float64}
    levels_Pa::Vector{Float64}
    gravity::Float64
    water_density::Float64
    state_count::Int
    flux_count::Int
    water_sum::Vector{Float64}
    ps_sum::Vector{Float64}
    temperature_sum::Matrix{Float64}
    temperature_count::Matrix{Int}
    precipitation_kg_m2::Vector{Float64}
end

function Accumulator(start,stop,npoints,sigma,dsigma,gravity,water_density;
                     step_ms=900000,levels_Pa=[20000.,50000.,85000.])
    span=Dates.value(stop-start)
    span>0 && span%3600000==0 || error("Window must contain whole hours")
    0<step_ms<=3600000 && 3600000%step_ms==0 || error("Step must divide one hour")
    npoints>0 || error("Empty grid")
    length(sigma)==length(dsigma)>=2 || error("Layer count mismatch")
    all(isfinite,sigma) && all(>(0),sigma) && all(<=(1),sigma) &&
        all(>(0),diff(sigma)) || error("Sigma must increase toward the surface")
    all(isfinite,dsigma) && all(>(0),dsigma) || error("Invalid layer thickness")
    isfinite(gravity) && gravity>0 && isfinite(water_density) && water_density>0 || error("Invalid constants")
    all(isfinite,levels_Pa) && all(>(0),levels_Pa) && all(>(0),diff(levels_Pa)) || error("Invalid target levels")
    return Accumulator(start,stop,step_ms,Float64.(sigma),Float64.(dsigma),Float64.(levels_Pa),
        Float64(gravity),Float64(water_density),0,0,zeros(npoints),zeros(npoints),
        zeros(npoints,length(levels_Pa)),zeros(Int,npoints,length(levels_Pa)),zeros(npoints))
end

# No extrapolation below the lowest full level or above the highest full level.
function pressure_temperature(temperature,sigma,ps,target)
    target<sigma[1]*ps && return nothing
    target>sigma[end]*ps && return nothing
    x=target/ps
    upper=searchsortedfirst(sigma,x)
    upper>length(sigma) && return nothing
    upper==1 && return Float64(temperature[1])
    sigma[upper]==x && return Float64(temperature[upper])
    lower=upper-1
    fraction=log(x/sigma[lower])/log(sigma[upper]/sigma[lower])
    return (1-fraction)*Float64(temperature[lower])+fraction*Float64(temperature[upper])
end

function state!(a::Accumulator,time,temperature,humidity,ps)
    time==a.start+Hour(a.state_count) && time<a.stop || error("Missing, duplicate or out-of-window hourly state")
    n=length(a.ps_sum)
    size(temperature)==size(humidity)==(n,length(a.sigma)) && length(ps)==n || error("State shape mismatch")
    all(isfinite,temperature) && all(isfinite,humidity) && all(isfinite,ps) && all(>(0),ps) || error("Invalid state")
    for i in 1:n
        a.ps_sum[i]+=Float64(ps[i])
        a.water_sum[i]+=sum(Float64(humidity[i,k])*a.dsigma[k] for k in eachindex(a.sigma))*Float64(ps[i])/a.gravity
        for j in eachindex(a.levels_Pa)
            value=pressure_temperature(view(temperature,i,:),a.sigma,Float64(ps[i]),a.levels_Pa[j])
            isnothing(value) && continue
            a.temperature_sum[i,j]+=value
            a.temperature_count[i,j]+=1
        end
    end
    a.state_count+=1
    return nothing
end

function flux!(a::Accumulator,time,rain,snow)
    time==a.start+Millisecond((a.flux_count+1)*a.step_ms) && time<=a.stop || error("Missing, duplicate or out-of-window flux")
    length(rain)==length(snow)==length(a.ps_sum) || error("Flux shape mismatch")
    all(isfinite,rain) && all(isfinite,snow) || error("Nonfinite precipitation")
    dt=a.step_ms/1000
    for i in eachindex(rain)
        a.precipitation_kg_m2[i]+=(Float64(rain[i])+Float64(snow[i]))*a.water_density*dt
    end
    a.flux_count+=1
    return nothing
end

function finish(a::Accumulator)
    duration_ms=Dates.value(a.stop-a.start)
    a.state_count==duration_ms÷3600000 || error("Incomplete hourly states")
    a.flux_count==duration_ms÷a.step_ms || error("Incomplete precipitation intervals")
    valid=a.temperature_count.==a.state_count
    temperature=fill(NaN,size(a.temperature_sum))
    temperature[valid].=a.temperature_sum[valid]./a.state_count
    seconds=duration_ms/1000
    return (;temperature_K=temperature,temperature_valid=valid,
        temperature_sample_count=copy(a.temperature_count),
        temperature_sum_K=copy(a.temperature_sum),
        column_water_kg_m2=a.water_sum./a.state_count,
        surface_pressure_Pa=a.ps_sum./a.state_count,
        precipitation_kg_m2=copy(a.precipitation_kg_m2),
        precipitation_kg_m2_s=a.precipitation_kg_m2./seconds,
        precipitation_mm_day=a.precipitation_kg_m2.*(86400/seconds),
        state_samples=a.state_count,flux_samples=a.flux_count,duration_seconds=seconds)
end
end
