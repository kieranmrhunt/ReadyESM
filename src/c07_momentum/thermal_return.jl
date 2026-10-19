include("energy_subcycles.jl")

struct C07ThermalReturn{C} <: S.AbstractConvection
    parent::C
    enabled::Bool
end
const ENERGY_FIELDS=(:c07_total_energy_residual,:c07_total_energy_gross,
    :c07_thermal_return,:c07_heat_rounding_units,:c07_heat_rounding_W_m2,
    :c07_heat_rounding_bound_W_m2,:c07_energy_checked)
S.variables(c::C07ThermalReturn)=(S.variables(c.parent)...,
    S.ParameterizationVariable(:c07_temperature_tendency,S.Grid3D();units="K/s",desc="C07 mixing heat"),
    (S.ParameterizationVariable(k,S.Grid2D();units=k in (:c07_heat_rounding_units,:c07_energy_checked) ? "1" : "W/m2",
        desc="C07 operator energy ledger") for k in ENERGY_FIELDS)...)
function S.initialize!(c::C07ThermalReturn,model::S.AbstractModel)
    model.architecture isa S.Architectures.AbstractCPU || error("CPU-only energy reference")
    return S.initialize!(c.parent,model)
end
function S.parameterization!(ij,vars,c::C07ThermalReturn,model)
    p=vars.parameterizations
    n=length(model.geometry.σ_levels_full)
    for name in ENERGY_FIELDS
        getproperty(p,name)[ij]=0
    end
    for k in 1:n
        p.c07_temperature_tendency[ij,k]=0
    end
    S.parameterization!(ij,vars,c.parent,model)
    c.enabled || return nothing
    iszero(c.parent.strength) && return nothing
    p.c07_triggered[ij]>0 || return nothing
    p.c07_skipped[ij]==0 || return nothing
    F=eltype(vars.grid.temperature_prev)
    d=C07HostAdapter._c07_column_and_composite(ij,vars,c.parent,model,F)
    u=F[vars.grid.u_prev[ij,k] for k in 1:n]
    v=F[vars.grid.v_prev[ij,k] for k in 1:n]
    mass=d.column.dp./d.thermo.g
    dt=c.parent.base.spectrum.shared.effective_step_factor*F(model.time_stepping.Δt_sec)
    cp=model.atmosphere.heat_capacity
    x=C07EnergySubcycles.exchange(C07HostAdapter.subcycled_transport_momentum,
        u,v,mass,d.composite.Mhalf,d.composite.E,d.composite.D,F(p.c08_mass_flux[ij]),dt;
        strength=c.parent.strength,maximum_subcycles=c.parent.maximum_subcycles,cp)
    all(x.du[k]==p.c07_u_tendency[ij,k] && x.dv[k]==p.c07_v_tendency[ij,k] for k in 1:n) ||
        error("Energy path changed momentum")
    ledger=C07EnergySubcycles.E.ledger(x,u,v,mass,dt,cp)
    p.c07_total_energy_residual[ij]=ledger.residual
    p.c07_total_energy_gross[ij]=ledger.gross
    p.c07_thermal_return[ij]=ledger.heating_W_m2
    p.c07_energy_checked[ij]=1
    # Replace the cancellation-prone diagnostic in this active prototype only.
    p.c07_kinetic_energy_dissipation[ij]=ledger.kinetic_loss_W_m2
    rounding_units=0.0
    rounding_energy=0.0
    rounding_bound=0.0
    for k in 1:n
        before=vars.tendencies.grid.temperature[ij,k]
        vars.tendencies.grid.temperature[ij,k]+=x.dT[k]
        after=vars.tendencies.grid.temperature[ij,k]
        p.c07_temperature_tendency[ij,k]=x.dT[k]
        actual=Float64(after)-Float64(before)
        difference=actual-Float64(x.dT[k])
        bound=Float64(eps(F))*max(abs(Float64(before)),abs(Float64(after)),Float64(floatmin(F)))
        rounding_units=max(rounding_units,abs(difference)/bound)
        rounding_energy+=Float64(mass[k])*Float64(cp)*difference
        rounding_bound+=Float64(mass[k])*Float64(cp)*bound
    end
    p.c07_heat_rounding_units[ij]=rounding_units
    p.c07_heat_rounding_W_m2[ij]=rounding_energy
    p.c07_heat_rounding_bound_W_m2[ij]=rounding_bound
    return nothing
end
