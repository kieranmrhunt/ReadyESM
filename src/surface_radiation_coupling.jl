"""
Live NumericalEarth surface-radiation component driven by the prognostic
SpeedyWeather/RRTMGP atmosphere.

SpeedyWeather's NumericalEarth exchanger already conservatively remaps
downwelling shortwave and longwave radiation onto the shared exchange grid.
This bridge copies those live fields into NumericalEarth's radiation state at
the net-flux assembly phase, then reuses NumericalEarth's native CPU/GPU
surface-radiation kernels to heat the ocean and sea ice.
"""
mutable struct AtmosphereDrivenRadiation{S, FT}
    surface_properties::S
    stefan_boltzmann_constant::FT
    interface_fluxes
end

struct AtmosphereDrivenRadiationProperties{FT, S}
    σ::FT
    surface_properties::S
end

function AtmosphereDrivenRadiation(
    ::Type{FT} = Float32;
    ocean_albedo = 0.06,
    ocean_emissivity = 0.98,
    # Match SpeedyWeather.OceanSeaIceAlbedo and the single emissivity used by
    # RRTMGP so atmosphere and surface close the same radiative boundary.
    sea_ice_albedo = 0.60,
    sea_ice_emissivity = 0.98,
) where {FT}
    surface_properties = (
        ocean = NumericalEarth.SurfaceRadiationProperties(
            FT(ocean_albedo),
            FT(ocean_emissivity),
        ),
        sea_ice = NumericalEarth.SurfaceRadiationProperties(
            FT(sea_ice_albedo),
            FT(sea_ice_emissivity),
        ),
    )
    return AtmosphereDrivenRadiation(
        surface_properties,
        FT(NumericalEarth.Radiations.default_stefan_boltzmann_constant),
        nothing,
    )
end

Oceananigans.TimeSteppers.time_step!(::AtmosphereDrivenRadiation, Δt) = nothing
NumericalEarth.EarthSystemModels.adopt_clock(r::AtmosphereDrivenRadiation, clock) = r

function NumericalEarth.EarthSystemModels.InterfaceComputations.ComponentExchanger(
    ::AtmosphereDrivenRadiation,
    grid,
)
    state = (
        ℐꜜˢʷ = Oceananigans.Field{
            Oceananigans.Center,
            Oceananigans.Center,
            Nothing,
        }(grid),
        ℐꜜˡʷ = Oceananigans.Field{
            Oceananigans.Center,
            Oceananigans.Center,
            Nothing,
        }(grid),
    )
    return NumericalEarth.EarthSystemModels.InterfaceComputations.ComponentExchanger(
        state,
        nothing,
    )
end

NumericalEarth.EarthSystemModels.interpolate_state!(
    exchanger,
    grid,
    ::AtmosphereDrivenRadiation,
    coupled_model,
) = nothing

function NumericalEarth.EarthSystemModels.update_net_fluxes!(
    coupled_model,
    ::AtmosphereDrivenRadiation,
)
    atmosphere_state = coupled_model.interfaces.exchanger.atmosphere.state
    radiation_state = coupled_model.interfaces.exchanger.radiation.state
    radiation_state.ℐꜜˢʷ .= atmosphere_state.ℐꜜˢʷ
    radiation_state.ℐꜜˡʷ .= atmosphere_state.ℐꜜˡʷ
    return nothing
end

function NumericalEarth.EarthSystemModels.allocate_interface_fluxes!(
    radiation::AtmosphereDrivenRadiation,
    exchange_grid,
    surfaces,
)
    pairs = (
        surface => NumericalEarth.InterfaceRadiationFlux(exchange_grid)
        for surface in surfaces
    )
    radiation.interface_fluxes = NamedTuple(pairs)
    return nothing
end

@inline function NumericalEarth.EarthSystemModels.InterfaceComputations.kernel_radiation_properties(
    radiation::AtmosphereDrivenRadiation,
)
    return AtmosphereDrivenRadiationProperties(
        radiation.stefan_boltzmann_constant,
        radiation.surface_properties,
    )
end

@inline function NumericalEarth.EarthSystemModels.InterfaceComputations.air_sea_interface_radiation_state(
    radiation_properties::AtmosphereDrivenRadiationProperties,
    exchanger_state,
    i,
    j,
    k,
    grid,
    time,
)
    return NumericalEarth.Radiations._surface_radiation_state(
        radiation_properties.surface_properties.ocean,
        radiation_properties,
        exchanger_state,
        i,
        j,
        k,
        grid,
        time,
    )
end

@inline function NumericalEarth.EarthSystemModels.InterfaceComputations.air_sea_ice_interface_radiation_state(
    radiation_properties::AtmosphereDrivenRadiationProperties,
    exchanger_state,
    i,
    j,
    k,
    grid,
    time,
)
    return NumericalEarth.Radiations._surface_radiation_state(
        radiation_properties.surface_properties.sea_ice,
        radiation_properties,
        exchanger_state,
        i,
        j,
        k,
        grid,
        time,
    )
end
