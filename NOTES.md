# Notes

## Current model

- SpeedyWeather 0.21.1, T31/L27, ERA5 initial atmosphere.
- RRTMGP all-sky radiation with 420 ppm prescribed CO2, diagnostic clouds and
  zero sulfate AOD in the production configuration.
- Oceananigans on a 360 x 180 x 60 tripolar grid, initialized from ECCO V4r4.
- ClimaSeaIce dynamics and thermodynamics. Ice and snow volume are advected
  directly and conservatively; velocity and tracer halos are refreshed first.
- Terrarium 0.1.6 `LandModel` with 16 soil layers and prescribed vegetation.
  Its new snow component is disabled for the current compatibility baseline,
  so snowfall follows the previously qualified liquid-input treatment.
- Terrarium runoff is transferred conservatively to nearby ocean cells.

ReadyESM carries small patches to Terrarium and SpeedyTransforms. Their source
and licences are included under `vendor/`.

## Main known limitations

- Recent production-path tests reproduced an illegal GPU memory access during
  the first sea-ice momentum step. A synchronization-only candidate did not
  resolve it. This remains under investigation; v0.1.3 does not fix it.
- Initial ice-surface/radiation reconciliation can give excessively cold
  surface temperatures. A checked startup solver is being tested separately
  and is not included in this maintenance release.
- The late model energy balance is still too positive (about +20 W m-2 all-sky
  and +45 W m-2 clear-sky in the completed control). The model is not spun up
  or calibrated for projections.
- Clouds are diagnosed from humidity; cloud liquid, ice and overlap are not a
  prognostic cloud system.
- Arctic ice can still converge into coastal cells. Direct volume transport
  fixes a real conservation error, but the remaining hotspot needs a coupled
  qualification and probably better coastal/ice-thickness physics.
- One-degree river mouths can retain very fresh surface lenses. Broad enhanced
  mixing was too intrusive; a narrowly loading-scaled closure is still being
  tested and is not selected in `production.yml`.
- Terrarium 0.1.6 has passed CPU and GPU water-conservation gates plus a coupled
  RRTMGP/Terrarium GPU step. Restart and long-run qualification remain.
- Forcing is concentration/AOD driven. Carbon-cycle, SO2 and aerosol-emissions
  modules are not present.

## Decisions retained from development

- The upstream-observed large-scale precipitation route is used after repairing
  the Betts-Miller rain diagnostic so column water removal matches rainfall.
- SpeedyWeather's zero-adjacent-transport vertical-diffusion control is used.
  The tested metric bulk-Richardson replacement worsened the 30-day climate.
- A strong sigma-dependent convection profile improved radiation and water
  drift but worsened temperature drift, so it is not in the production config.
- No polar sponge, global depth correction, arbitrary ice-thickness cap or
  broad river-mouth diffusivity is applied.
- Initial conditions are inserted directly. A balanced incremental launch is a
  useful future experiment, not part of this configuration.

The internal development labbook and rejected experiment matrix are intentionally
not shipped in this repository.
