module C08FiniteConvection

include("deterministic_multiplume_core.jl")
const C08MultiPlume = DeterministicMultiPlumeC08
include("finite_response.jl")
include("finite_closure.jl")
include("mixed_precision.jl")

const solve = MixedPrecisionC08.solve
const Column = C08MultiPlume.Column
const typed_thermo = C08MultiPlume.typed_thermo
const typed_plume_spectrum = C08MultiPlume.typed_plume_spectrum

export Column, typed_thermo, typed_plume_spectrum, solve
end
