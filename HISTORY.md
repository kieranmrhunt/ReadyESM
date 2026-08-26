# History

The first tagged repository version is `0.1.0`. Earlier entries below are
development milestones reconstructed from the internal labbook, not releases.

## 0.1.0 — 2026-08-26

First GitHub version. It contains the accepted T31/L27 tripolar configuration,
Terrarium 0.1.6, Oceananigans 0.110.15, conservative atmospheric water handling,
direct sea-ice/snow volume transport, restart support and final-state diagnostics.

## 24–25 August 2026

Closed the diagnosed Betts-Miller rainfall water source, made rejected metric
vertical diffusion opt-in, and repaired sea-ice transport to advect represented
ice and snow volume directly. A strong sigma-convection change was rejected on
temperature drift despite improving radiation and column water.

## 17–23 August 2026

Established the T31/L27, 1° tripolar production identity and exact coupled
restart path. A 730-day control completed. It exposed the remaining positive
TOA imbalance, local ocean salinity and heat failures, and coastal Arctic ice
pile-up. Terrarium 0.1.2 was used at this stage.

## 12–16 August 2026

Stabilised the full SpeedyWeather–RRTMGP–Oceananigans–ClimaSeaIce–Terrarium
stack, moved to ERA5/ECCO initial conditions and a high-top L27 atmosphere, and
fixed several nondeterministic GPU operations. Earlier low-latitude ice, static
soil moisture and invalid surface-temperature diagnostics were superseded.

## Initial prototype

SpeedyWeather was first coupled to RRTMGP, a slab ocean, thermodynamic sea ice
and simple land. This proved the component interfaces and prescribed CO2/AOD
forcing, but was not the dynamic Earth-system configuration retained here.

## Retired estuary experiments

The old `v1`–`v5` labels referred to river-mouth mixing experiments, not ReadyESM
releases: fixed broad mixing, narrowed mixing, salinity-only mixing, active-runoff
gating, then freshwater-loading scaling. None is enabled in `production.yml`;
the salinity problem remains open.
