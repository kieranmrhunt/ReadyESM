"""
Earth atmosphere constants with saturation over liquid water and ice.

SpeedyWeather 0.21.1 applies liquid-water Clausius--Clapeyron at every
temperature, so its large-scale condensation, Betts--Miller convection and
ReadyESM's diagnostic cloud search never see ice saturation. At 210 K that law
is 2.2 times ice saturation, and the cold upper troposphere stays strongly
ice-supersaturated.

`MixedPhaseEarthAtmosphere` keeps every `EarthAtmosphere` constant and changes
only the saturation vapour pressure: SpeedyWeather's own liquid expression at
and above freezing, ice Clausius--Clapeyron below 250.16 K and an IFS-style
quadratic blend between. The ice curve uses vaporisation plus fusion heat, the
energy SpeedyWeather's condensation already releases when it forms snow; it is
within 2% of Murphy and Koop (2005) at 210 K, whereas the separate 2801 kJ/kg
sublimation constant would overestimate ice saturation there by about 9%.

Latent heating and the derivative approximations in Betts--Miller's moist
adiabat and implicit large-scale condensation are unchanged. Those inherited
approximations are not fully consistent with the mixed-phase saturation curve.
This is the tested saturation-only change, not a complete phase-aware
thermodynamic closure; derivative and parcel-heating repairs are separate
experiments.
"""
struct MixedPhaseEarthAtmosphere{A <: SpeedyWeather.EarthAtmosphere} <:
       SpeedyWeather.AbstractWetAtmosphere
    earth::A
end

Adapt.@adapt_structure MixedPhaseEarthAtmosphere

MixedPhaseEarthAtmosphere(spectral_grid::SpeedyWeather.SpectralGrid; kwargs...) =
    MixedPhaseEarthAtmosphere(SpeedyWeather.EarthAtmosphere(spectral_grid; kwargs...))

# SpeedyWeather reads constants as fields; serve all of them from EarthAtmosphere.
@inline Base.getproperty(A::MixedPhaseEarthAtmosphere, name::Symbol) =
    name === :earth ? getfield(A, :earth) : getproperty(getfield(A, :earth), name)
Base.propertynames(A::MixedPhaseEarthAtmosphere) = propertynames(getfield(A, :earth))
Base.eltype(::Type{MixedPhaseEarthAtmosphere{SpeedyWeather.EarthAtmosphere{NF}}}) where {NF} = NF
Base.eltype(A::MixedPhaseEarthAtmosphere) = eltype(typeof(A))

const MIXED_PHASE_ALL_ICE_TEMPERATURE_K = 250.16

@inline function SpeedyWeather.saturation_vapor_pressure(T, A::MixedPhaseEarthAtmosphere)
    earth = getfield(A, :earth)
    T₀ = earth.temperature_freezing
    NF = typeof(T₀)
    T = convert(NF, T)
    liquid = SpeedyWeather.saturation_vapor_pressure(T, earth)
    sublimation = earth.latent_heat_condensation + earth.latent_heat_fusion
    ice = earth.saturation_vapor_pressure *
        exp(sublimation / earth.R_vapor * (inv(T₀) - inv(T)))
    T_ice = NF(MIXED_PHASE_ALL_ICE_TEMPERATURE_K)
    α = clamp((T - T_ice) / (T₀ - T_ice), zero(NF), one(NF))^2
    return ifelse(T >= T₀, liquid, α * liquid + (1 - α) * ice)
end

const _ATMOSPHERE_SATURATION_LAWS = (:liquid, :mixed_phase)

"""Thermodynamic constants for the configured saturation law."""
function _atmosphere_thermodynamics(config::ExperimentConfig, spectral_grid)
    config.atmosphere_saturation_law == :liquid &&
        return SpeedyWeather.EarthAtmosphere(spectral_grid)
    config.atmosphere_saturation_law == :mixed_phase &&
        return MixedPhaseEarthAtmosphere(spectral_grid)
    throw(ArgumentError("unsupported atmosphere_saturation_law $(config.atmosphere_saturation_law)"))
end

_atmosphere_saturation_provenance(config::ExperimentConfig) = (
    atmosphere_saturation_law = String(config.atmosphere_saturation_law),
    atmosphere_saturation_all_ice_temperature_K =
        config.atmosphere_saturation_law == :mixed_phase ?
        MIXED_PHASE_ALL_ICE_TEMPERATURE_K : NaN,
)
