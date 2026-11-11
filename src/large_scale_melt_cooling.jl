# Opt-in snowmelt cooling for SpeedyWeather 0.21.1's implicit large-scale condensation (Claude, 1 Nov 2026).
#
# SpeedyWeather.jl is licensed under the European Union Public Licence v1.2 (EUPL-1.2), the same licence as ReadyESM:
# Copyright (c) 2020 Milan Kloewer and (c) 2021 The SpeedyWeather.jl Contributors for SpeedyWeather.jl, and (c) 2022
# Fred Kucharski and Franco Molteni for the SPEEDY parametrization schemes. The column routine below is a copy of
# `large_scale_condensation!` in SpeedyWeather 0.21.1 src/parameterizations/large_scale_condensation.jl with two
# changes, the same as Codex's melt-energy repair (b6a9dd4) for the unreleased process-conservative candidate:
#   1. the melt cooling is a tendency, -Lᵢ/cₚ δq_melt / Δt, so the temperature in the implicit derivative is the
#      temperature after one step of melting (upstream adds Δt times a per-step change);
#   2. the condensation heating is added to the melt cooling instead of replacing it.
# Upstream computes the melt cooling and then overwrites it, while freezing still releases Lᵢ. Snow that forms aloft
# and melts before the surface therefore added Lᵢ per kilogram with no sink: about 0.6 to 0.8 W/m² in the v0.1.19
# July and January months (artifacts/development_wave_20261101/lsc_melt_energy_claude_v1). With the fix, the column
# heat input of the scheme, sum of (cₚ δT + Lᵥ δq) times the layer mass, equals Lᵢ times the snow leaving the column.
# Every other operation, parameter and output is unchanged. Where no snow melts, the tendencies are bit-identical to
# upstream. One passive output is added: large_scale_snowmelt_rate, the column total of melted snow in this step
# (m/s of water), so the melt cooling can be checked directly, as Lᵢ ρ_w times this rate.

"""
SpeedyWeather 0.21.1 `ImplicitCondensation` with the cooling of melting snow kept. It holds an upstream
`ImplicitCondensation` for its parameters and output variables.
"""
struct MeltCoolingImplicitCondensation{C} <: SpeedyWeather.AbstractCondensation
    condensation::C
end

Adapt.@adapt_structure MeltCoolingImplicitCondensation

function MeltCoolingImplicitCondensation(spectral_grid::SpeedyWeather.SpectralGrid; kwargs...)
    condensation = SpeedyWeather.ImplicitCondensation(spectral_grid; kwargs...)
    return MeltCoolingImplicitCondensation{typeof(condensation)}(condensation)
end

SpeedyWeather.variables(condensation::MeltCoolingImplicitCondensation) = (
    SpeedyWeather.variables(condensation.condensation)...,
    SpeedyWeather.ParameterizationVariable(
        :large_scale_snowmelt_rate,
        SpeedyWeather.Grid2D();
        desc = "Large-scale snow melted while falling through the column this step, as water",
        units = "m/s",
    ),
)

SpeedyWeather.initialize!(condensation::MeltCoolingImplicitCondensation, model::SpeedyWeather.AbstractModel) =
    SpeedyWeather.initialize!(condensation.condensation, model)

Base.@propagate_inbounds function SpeedyWeather.parameterization!(
    ij,
    vars,
    condensation::MeltCoolingImplicitCondensation,
    model,
)
    (; geometry, planet, atmosphere, time_stepping) = model
    _melt_cooling_large_scale_condensation!(
        ij, vars, condensation.condensation, geometry, planet, atmosphere, time_stepping)
    return nothing
end

Base.@propagate_inbounds function _melt_cooling_large_scale_condensation!(
    ij,
    vars,
    condensation::SpeedyWeather.ImplicitCondensation,
    geometry,
    planet,
    atmosphere,
    time_stepping,
)
    # use previous time step for more stable Euler forward step of the parameterizations
    temp = vars.grid.temperature_prev                  # temperature [K]
    humid = vars.grid.humidity_prev                # specific humidity [kg/kg]
    temp_tend = vars.tendencies.grid.temperature       # temperature tendency [K/s]
    humid_tend = vars.tendencies.grid.humidity     # specific humidity tendency [kg/kg/s]
    nlayers = size(temp, 2)

    # precompute scaling constants to minimize divisions (used to convert between humidity [kg/kg] and precipitation [m])
    pₛ = vars.grid.pressure_prev[ij]               # surface pressure [Pa]
    (; Δt_sec) = time_stepping                      # time step [s]
    σ = geometry.σ_levels_full                      # vertical sigma coordinate [1]
    Δσ = geometry.σ_levels_thick                    # layer thickness in sigma coordinates
    ρ = atmosphere.water_density                    # air density [kg/m³]
    pₛ_gρ = pₛ / (planet.gravity * ρ)                 # [Pa / (m/s² * kg/m³)] = [Pa m² * s² / kg] = [m]

    # thermodynamics
    Lᵥ = atmosphere.latent_heat_condensation        # latent heat of vaporization
    Lᵢ = atmosphere.latent_heat_fusion              # latent heat of fusion
    cₚ = atmosphere.heat_capacity                   # heat capacity
    Rᵥ = atmosphere.R_vapor                         # gas constant for water vapor
    Lᵥ_cₚ = Lᵥ / cₚ
    Lᵢ_cₚ = Lᵢ / cₚ
    cₚ_Lᵢ = cₚ / Lᵢ

    (; time_scale, relative_humidity_threshold, freezing_threshold, melting_threshold) = condensation
    let_it_snow = condensation.snow                     # flag to switch snow on/off
    rain_flux_down = zero(eltype(temp))                 # start with zero rain/snow flux at top of atmosphere
    snowmelt_rate = zero(eltype(temp))                  # ReadyESM: column total of melted snow [m/s]
    snow_flux_down = zero(eltype(temp))                 # both have units of [m/s] (precipitation)

    for k in 1:nlayers

        # Condensation from humidity in this layer (for a negative humidity tendency)
        # relative to threshold that can be <100%, e.g. 95%
        sat_humid_k = SpeedyWeather.saturation_humidity(temp[ij, k], σ[k] * pₛ, atmosphere)
        δq_cond = sat_humid_k * relative_humidity_threshold - humid[ij, k]

        # skip if no condensation has occurred yet in this layer or above
        if δq_cond < 0 || snow_flux_down > 0 || rain_flux_down > 0

            # 0. convert between humidity tendency [kg/kg/s] and precipitation amount [m] or rate [m/s]
            Δp_gρ = Δσ[k] * pₛ_gρ                           # pressure thickness of layer Δp times 1/g/ρ [m]
            Δp_Δtgρ = Δp_gρ / Δt_sec                          # pressure thickness of layer Δp times 1/Δt/g/ρ [m/s]
            Δtgρ_Δp = inv(Δp_Δtgρ)                          # [s/m]

            # 1. Melting of snow from layer above
            δT_melt = max(temp[ij, k] - melting_threshold, 0)   # only if temperature above melting threshold
            E_avail = cₚ_Lᵢ * δT_melt                       # available energy [1] in multiples of latent heat [J/kg of air] for melting
            melt_rate_max = E_avail * Δp_Δtgρ               # max rate of snow melt [m/s]
            melt_rate = min(snow_flux_down, melt_rate_max)  # cap to snow flux [m/s], don't melt more than available
            snow_flux_down -= melt_rate                     # move melted snow to rain
            snowmelt_rate += melt_rate
            rain_flux_down += melt_rate                     # this can evaporate now too
            δq_melt = melt_rate * Δtgρ_Δp                   # convert back to humidity increase over timestep [kg/kg]
            # ReadyESM: the melt cooling as a tendency [K/s]. Upstream stores the per-step change here and then
            # overwrites it below, so melting snow released no cooling.
            δT = -Lᵢ_cₚ * δq_melt / Δt_sec
            # don't use δq_melt as humidity tendency as it's only snow->rain

            # 2. Reevaporation and condensation
            dq = sat_humid_k - humid[ij, k]                 # reevaporation relative to 100% humidity
            r = condensation.reevaporation * max(dq, 0)     # reevaporation efficiency [1], scale linearly with dryness
            rain_evaporated = min(r, 1) * rain_flux_down    # [m/s], min(r, 1) to not evaporate more than available
            rain_flux_down -= rain_evaporated               # remove reevaporated rain
            δq_evap = rain_evaporated * Δtgρ_Δp             # convert to humidity tendency over timestep [kg/kg]
            δq = min(0, δq_cond) + δq_evap                  # [kg/kg] sum with condensation (negative) and evaporation (positive)
            # division by timestep for tendency below

            # Solve for melting of snow, condensation, reevaporation (and possibly sublimation) implicitly in time
            # implicit correction, Frierson et al. 2006 eq. (21)
            # derivative of qsat wrt to temp
            T = temp[ij, k] + Δt_sec * δT                     # use temperature minus latent heat for melting in gradient
            dqsat_dT = sat_humid_k * relative_humidity_threshold * Lᵥ_cₚ / (Rᵥ * T^2)
            δq /= ((1 + Lᵥ_cₚ * dqsat_dT) * time_scale * Δt_sec)
            # ReadyESM: add the condensation heating to the melt cooling instead of replacing it.
            δT -= Lᵥ_cₚ * δq

            # If there is large-scale condensation at a level higher (i.e. smaller k) than
            # the cloud-top previously diagnosed due to convection, then increase the cloud-top
            # Fortran SPEEDY documentation Page 7 (last sentence)
            vars.parameterizations.cloud_top[ij] = min(vars.parameterizations.cloud_top[ij], k)

            # 2. Precipitation (rain) due to large-scale condensation [kg/m²/s] /ρ for [m/s]
            δq_rain = min(0, δq)                # precipitation only for negative humidity tendency
            rain = Δp_gρ * -δq_rain             # precipitation rain [m/s] on this layer k, Formula 4
            snow = zero(rain)                   # start with zero snow but potentially swap below

            # decide whether to turn precip into snow (all rain freezes to snow)
            rain, snow = let_it_snow && temp[ij, k] < freezing_threshold ? (snow, rain) : (rain, snow)
            rain_flux_down += rain              # accumulate into downward fluxes [m/s] (used in layer below)
            snow_flux_down += snow              # accumulate into downward fluxes [m/s] (used in layer below)

            # latent heat release when freezing for enthalpy conservation
            δT -= (snow > 0) * Lᵢ_cₚ * δq_rain  # snow for negative δq (condensation then freezing)

            # only accumulate into humid_tend now to allow humid_tend != 0 before this scheme is called
            humid_tend[ij, k] += δq             # [kg/kg/s]
            temp_tend[ij, k] += δT              # [K/s]
        end
    end

    vars.parameterizations.large_scale_snowmelt_rate[ij] = snowmelt_rate

    # avoid negative precipitation from rounding errors
    rain_flux_down = max(rain_flux_down, 0)
    snow_flux_down = max(snow_flux_down, 0)

    # precipitation from rain/snow [m] during time step whatever is fluxed out
    vars.parameterizations.rain_large_scale[ij] += Δt_sec * rain_flux_down
    vars.parameterizations.snow_large_scale[ij] += Δt_sec * snow_flux_down

    # and the rain/snow fall rate [m/s]
    vars.parameterizations.rain_rate_large_scale[ij] = rain_flux_down
    vars.parameterizations.snow_rate_large_scale[ij] = snow_flux_down

    # accumulate into total rain/snow rate including convection [m/s]
    vars.parameterizations.rain_rate[ij] += rain_flux_down
    vars.parameterizations.snow_rate[ij] += snow_flux_down

    return nothing
end
