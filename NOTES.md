# Notes

## Current model

- SpeedyWeather 0.21.1, T31/L27, ERA5 initial atmosphere.
- RRTMGP all-sky radiation with 420 ppm prescribed CO2, diagnostic clouds and
  zero sulfate AOD. Liquid/ice cloud paths are 75/25 g/m².
- Oceananigans on a 360 x 180 x 60 tripolar grid, initialized from ECCO V4r4.
- ClimaSeaIce dynamics and thermodynamics. Ice and snow volume are advected
  directly and conservatively; velocity and tracer halos are refreshed first.
- Terrarium 0.1.6 `LandModel` with 16 soil layers and prescribed vegetation.
  Land state and water accounting use Float64; atmosphere and ocean use Float32.
  Its new snow component is disabled for the current compatibility baseline,
  so snowfall follows the previously qualified liquid-input treatment.
- Terrarium runoff is accumulated at every native land step and transferred
  to unique physical ocean cells. Pending water and the held discharge are
  included in checkpoints.

ReadyESM carries small patches to Terrarium, SpeedyTransforms and ClimaSeaIce.
Their source and licences are included under `vendor/`.

## Main known limitations

- Earlier versions have an intermittent GPU illegal-access failure during
  sea-ice stepping. The first v0.1.5 candidate passes two short restart runs,
  including normal package loading, but the subsequent 48-hour attempt fails
  in its first resumed step. Exact saved-state restoration passes. The
  initiating operation and change in failure frequency remain unproven.
  With the added snow type change, an ordinary restart from hour 2 to hour 48
  now passes full diagnostics and exact restored-state comparison. The GPU
  and driver also changed, so patch causality and general restart reliability
  remain unqualified. A fresh restart from the new hour-48 state also passes,
  continuing through hour 49 with full diagnostics and exact restored state.
  See the [validation record](docs/validation-v0.1.5.md).
- Initial ice-surface/radiation reconciliation can give excessively cold
  surface temperatures. A checked startup solver is being tested separately
  and is not included in v0.1.6. Its later ice-interface non-convergence
  also remains under investigation.
- The model energy balance remains too positive: +16.6 W/m² all-sky over
  days 20–30, versus +21.7 in the matched cloud60 control. Clear-sky net flux
  is almost unchanged (+43.6 versus +43.5 W/m²); zero clear-sky flux is not
  an independent equilibrium target. The model is not spun up or calibrated.
- Clouds are diagnosed from humidity; cloud liquid, ice and overlap are not a
  prognostic cloud system.
- Arctic ice can still converge into coastal cells. Direct volume transport
  fixes a real conservation error, but the remaining hotspot needs a coupled
  qualification and probably better coastal/ice-thickness physics.
- The pinned NumericalEarth air–ice flux kernel uses zero surface ice velocity
  when calculating relative wind, despite the evolving ice dynamics. This
  approximation needs review, including consistent vector mapping on the
  tripolar grid. Its role in weak-wind solver failures is not established.
- One-degree river mouths can retain very fresh surface lenses. Broad enhanced
  mixing was too intrusive; a narrowly loading-scaled closure is still being
  tested and is not selected in `production.yml`.
- Double-precision land closes its water budget to 2.8e-10 mm through the
  coupled 30-day run, including exact native-precision checkpoint restoration.
  This does not qualify century-scale closure or conversion of old checkpoints.
  Atmospheric water accounting is separate: its day-30 residual is
  -0.178 kg/m², versus -0.099 in the cloud60 control. The normal diagnostic
  checks require finite later residuals but do not impose a quantitative
  atmospheric closure bound.
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
