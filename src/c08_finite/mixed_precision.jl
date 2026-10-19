# Evaluate the complete closure in one precision, then cast its physical exchange.
module MixedPrecisionC08
const R=parentmodule(@__MODULE__)
const C08=R.C08MultiPlume
const CDI=C08.EP.ConvectiveDetrainmentInterface

field_values(x)=NamedTuple{fieldnames(typeof(x))}(
    Tuple(getfield(x,n) for n in fieldnames(typeof(x))))

function wide_inputs(column,thermo,spectrum)
    col=C08.typed_column(Float64,column)
    th=C08.typed_thermo(Float64;field_values(thermo)...)
    pp=C08.typed_plume_parameters(Float64;field_values(spectrum.shared)...)
    sp=C08.PlumeSpectrum{Float64}(Float64.(spectrum.entrainment_m),
        Float64.(spectrum.weights),pp,spectrum.closure_iterations)
    return col,th,sp
end

function cast_exchange(x,::Type{F}) where F
    CDI.ColumnExchange(F.(x.dqv),F.(x.dql),F.(x.dqi),F.(x.dT),F(x.rain),F(x.snow))
end

function solve(column,timestep,thermo,spectrum;enabled=true)
    enabled || return R.FiniteC08Closure.solve(column,timestep,thermo,spectrum)
    col,th,sp=wide_inputs(column,thermo,spectrum)
    exchange,diagnostics=R.FiniteC08Closure.solve(col,Float64(timestep),th,sp)
    # Diagnostics describe the wide solve. Cast only the returned physical
    # tendencies and precipitation to the atmospheric storage precision.
    return cast_exchange(exchange,eltype(column.T)),diagnostics
end
end
