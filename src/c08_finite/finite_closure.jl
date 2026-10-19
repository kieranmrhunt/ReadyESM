module FiniteC08Closure
const FR=parentmodule(@__MODULE__).FiniteC08Response
const C08=FR.C08
const CDI=FR.CDI

"""
Choose mass flux from the finite host-state CAPE response over the effective
step. The target is exp(-step/relaxation_time) times initial active-bin CAPE.

Scan the inherited entrainment domain in increasing flux order. Refine the
first sampled target crossing. If no crossing exists, retain the strongest
valid stabilising response found and report partial adjustment. A column
whose sampled exchange never reduces CAPE receives zero exchange with an
explicit no_stabilising_response flag. No invalid state enters a CAPE call.
"""
function solve(column,dt,thermo,spectrum;intervals=16,iterations=12)
    intervals>=1 && iterations>=0 || throw(ArgumentError("Invalid search resolution"))
    isfinite(spectrum.shared.time_scale_s) && spectrum.shared.time_scale_s>0 ||
        throw(ArgumentError("Positive finite relaxation time required"))
    c=FR.prepare(column,dt,thermo,spectrum)
    F=eltype(column.T)
    target=c.cape*exp(-c.step/spectrum.shared.time_scale_s)
    evaluations=0
    invalid_probes=0
    function checked_probe(M)
        p=FR.probe(c,M)
        evaluations+=1
        invalid_probes+=!p.valid
        return p
    end
    chosen=checked_probe(zero(F))
    chosen.valid || error("Initial state has no valid zero exchange")
    chosen_M=zero(F)
    status=:inactive
    bracket_width=zero(F)
    if any(c.active) && c.cape>0 && target<c.cape
        status=:no_stabilising_response
        previous_M=zero(F)
        bracket_found=false
        stopped_invalid=false
        for j in 1:intervals
            M=F(j)/F(intervals)*c.bound
            p=checked_probe(M)
            if !p.valid
                stopped_invalid=true
                break
            end
            if p.cape<chosen.cape
                chosen,chosen_M=p,M
                status=:partial_adjustment
            end
            if p.cape<=target
                lower,upper=previous_M,M
                upper_probe=p
                for _ in 1:iterations
                    middle=(lower+upper)/F(2)
                    middle==lower || middle==upper ? break : nothing
                    candidate=checked_probe(middle)
                    if !candidate.valid
                        stopped_invalid=true
                        break
                    elseif candidate.cape<=target
                        upper,upper_probe=middle,candidate
                    else
                        lower=middle
                    end
                end
                chosen,chosen_M=upper_probe,upper
                bracket_width=upper-lower
                status=stopped_invalid ? :valid_crossing_with_invalid_midpoint : :target_crossing
                bracket_found=true
                break
            end
            previous_M=M
        end
        if !bracket_found && chosen_M>0
            status=stopped_invalid ? :admissibility_limited :
                chosen_M==c.bound ? :entrainment_limited : :partial_adjustment
        elseif !bracket_found && stopped_invalid && iszero(chosen_M)
            status=:no_valid_stabilising_response
        end
    end
    triggered=chosen_M>=F(C08.MIN_MASS_FLUX)
    if !triggered && chosen_M!=0
        chosen=checked_probe(zero(F))
        chosen_M=zero(F)
        status=:rounding_zero
    end
    response=triggered ? max(c.cape-chosen.cape,zero(F))/(chosen_M*c.step) : zero(F)
    diagnostics=(triggered,Mb=chosen_M,active=c.active,
        triggers=[C08.bin_triggers(c.profiles[b],c.parameters[b]) for b in eachindex(c.active)],
        cape=[p.cape for p in c.profiles],cin=[p.cin for p in c.profiles],
        base=[p.base for p in c.profiles],top=[p.top for p in c.profiles],
        bin_Mb=[c.active[b] ? spectrum.weights[b]*chosen_M : zero(F) for b in eachindex(c.active)],
        R=response,R_steps=F[maximum(abs,c.step.*chosen.host_exchange.dT)],
        limited=status in (:entrainment_limited,:admissibility_limited),degenerate=false,
        finite_closure_status=status,finite_closure_evaluations=evaluations,
        finite_closure_invalid_probes=invalid_probes,finite_closure_target=target,
        finite_closure_initial_cape=c.cape,finite_closure_final_cape=chosen.cape,
        finite_closure_mass_bound=c.bound,finite_closure_bracket_width=bracket_width,
        finite_closure_temperature_minimum=minimum(chosen.temperature),
        finite_closure_temperature_maximum=maximum(chosen.temperature),
        finite_closure_vapor_minimum=minimum(chosen.vapor))
    return chosen.exchange,diagnostics
end
end
