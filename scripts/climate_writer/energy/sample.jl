# Read-only native component sampler. No live regridding buffers are touched.
include("energy.jl")

function climate_energy_samples(R, earth, unique_ocean_mask)
    atmosphere = earth.atmosphere
    land_model = atmosphere.model.land.model
    isnothing(land_model.snow) || error("energy observer requires snow-free Terrarium")
    land_fraction = vec(R._host_array(atmosphere.model.land_sea_mask.mask))
    if earth.clock.iteration == 0
        # An initial build has not completed a land coupling interval.
        # Initial files are outside monthly accumulators. Mark unavailable
        # fluxes rather than interpreting constructor zeros as observations.
        return Dict{String,Any}(name => fill(NaN, kind == :surface ? size(unique_ocean_mask) : size(land_fraction))
            for (name, (kind, _, _)) in ClimateEnergy.SPECS)
    end
    ocean_flux = earth.interfaces.atmosphere_ocean_interface.fluxes.latent_heat
    ice_flux = earth.interfaces.atmosphere_sea_ice_interface.fluxes.latent_heat
    concentration = earth.sea_ice.model.ice_concentration
    result = Dict{String,Any}(ClimateEnergy.ocean_samples(R._interior_matrix(ocean_flux), R._interior_matrix(ice_flux),
        R._interior_matrix(concentration), unique_ocean_mask))
    merge!(result, ClimateEnergy.vapor_samples(
        R._interior_matrix(earth.interfaces.atmosphere_ocean_interface.fluxes.water_vapor),
        R._interior_matrix(earth.interfaces.atmosphere_sea_ice_interface.fluxes.water_vapor),
        R._interior_matrix(concentration), unique_ocean_mask))
    humidity = vec(R._host_array(atmosphere.variables.prognostic.land.surface_humidity_flux))
    latent_heat = land_model.constants.thermodynamics.latent_heat_vaporization
    merge!(result, ClimateEnergy.land_samples(humidity, land_fraction, latent_heat))
    return result
end
