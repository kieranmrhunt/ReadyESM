include(joinpath(@__DIR__, "c07_momentum", "host_adapter.jl"))

"""Add conservative momentum transport and local mixing heat to finite C08."""
function ConvectiveMomentumTransport(
    base::FiniteMultiPlumeConvection;
    strength::Real = 0.5,
    maximum_subcycles::Int = 16,
)
    momentum = C07HostAdapter.C07FluxFormCMTConvection(
        base; strength, maximum_subcycles)
    return C07HostAdapter.C07ThermalReturn(momentum, true)
end

function _convective_momentum_metadata(config::ExperimentConfig)
    config.atmosphere_convection == :finite_multiplume_c07_v1 || return (;)
    return (
        atmosphere_c07_momentum_strength = config.atmosphere_c07_momentum_strength,
        atmosphere_c07_maximum_subcycles = config.atmosphere_c07_maximum_subcycles,
        atmosphere_c07_transport = "conservative_flux_form_convex_subcycles",
        atmosphere_c07_energy_return = "local_kinetic_to_thermal_with_subcycle_transport",
        atmosphere_c07_profile_precision = "Float64_selection_Float32_transport",
    )
end
