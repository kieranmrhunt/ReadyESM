# Notes

## Current model

- SpeedyWeather 0.21.1, T31/L27, ERA5 initial atmosphere.
- Production removes the extra `l = lmax + 1` row from the ERA5 spectral
  launch before the first grid transform. The generic option remains off.
- Two-hour Betts–Miller convection, reference RH0.7, net-column rain diagnosis.
- Mixed liquid/ice saturation and net-column large-scale rain diagnosis:
  surface precipitation equals that scheme's net column water removal,
  while transport, filtering and coupling remain in separate ledgers. The
  inherited moist-adiabat and condensation derivative approximations remain.
  A fully phase-consistent thermodynamic closure is future work.
- RRTMGP all-sky radiation with 420 ppm prescribed CO2, diagnostic clouds and
  zero sulfate AOD. Liquid/ice cloud paths are 75/25 g/m².
- Oceananigans on a 360 x 180 x 60 tripolar grid. v0.1.10 initializes
  ocean fields halfway between the December 1992 and January 1993 ECCO V4r4
  monthly means. Sea ice uses the January state.
  Main now separates the surface-height and grid-scale updates to fix stale
  neighbouring reads. See the [ordering tests](docs/validation-ocean-grid-ordering.md).
- ClimaSeaIce dynamics and thermodynamics. Ice and snow volume are advected
  directly and conservatively. Velocity and tracer halos are refreshed first.
- Terrarium 0.1.6 `LandModel` with 16 soil layers and prescribed vegetation.
  Land state and water accounting use Float64. Atmosphere and ocean use Float32.
  Its new snow component is disabled for the current compatibility baseline,
  so snowfall follows the previously qualified liquid-input treatment.
- Terrarium runoff is accumulated at every native land step and transferred
  to unique physical ocean cells. Pending water and the held discharge are
  included in checkpoints.

ReadyESM carries small patches to Terrarium, SpeedyTransforms and ClimaSeaIce.
Their source and licences are included under `vendor/`.
The Oceananigans grid-ordering override is in `src/ordered_ocean_grid.jl`.

## Main known limitations

- Earlier versions intermittently failed with GPU illegal access during
  sea-ice stepping. v0.1.6 and the v0.1.7 convection branch subsequently passed
  through day120 with native Float64 land and exact restored
  states. The GPU and driver also changed between failed and successful
  attempts. The cause of those illegal-access failures remains unresolved.
  See the [v0.1.5 record](docs/validation-v0.1.5.md) and
  [v0.1.7 record](docs/validation-v0.1.7.md).
- Initial ice-surface/radiation reconciliation can give excessively cold
  surface temperatures. Production keeps the qualified startup while a new
  solver and its later ice-interface non-convergence are investigated.
- The pre-fix v0.1.8 D experiment reaches day390, with +3.49 W/m² all-sky net
  TOA uptake over days360–384 versus +8.75 for the matched liquid/upstream control.
  D is also 1.36 K cooler and has 4.07 kg/m² less column water in that window.
  Temperature/water gaps grow between day 120 and day 384, though not
  monotonically. Ocean heat uptake remains 9.25 W/m² per global area,
  exceeding TOA uptake by 5.76 W/m² (control: 3.72). The 5.76 W/m² gap remains
  unexplained. A component-storage ledger is needed. These are transient
  single trajectories without equilibrium or calibration. Zero clear-sky
  flux is no independent equilibrium target.
- The analytic ozone profile has no latitude/season dependence. A preliminary
  fixed-column comparison gives about 460 DU versus 245–378 DU in the sampled
  ERA5 columns. Its upper-tropospheric radiative effect needs validation.
  v0.1.10 keeps the analytic ozone profile.
- Long-run GPU trajectory spread is unmeasured. The long D run has one
  trajectory. Short full-grid repeats with the ordered ocean update match
  exactly.
- Clouds are diagnosed from humidity. A prognostic system would carry cloud
  liquid, ice and overlap.
- Arctic ice can still converge into coastal cells. Direct volume transport
  fixes a real conservation error, but the remaining hotspot needs a coupled
  qualification and probably better coastal/ice-thickness physics.
- The pinned NumericalEarth air–ice flux kernel uses zero surface ice velocity
  when calculating relative wind, despite the evolving ice dynamics. This
  approximation needs review, including consistent vector mapping on the
  tripolar grid. Its role in weak-wind solver failures remains unknown.
- One-degree river mouths can retain very fresh surface lenses. Broad enhanced
  mixing was too intrusive. Production keeps native mixing while a narrowly
  loading-scaled closure is tested.
- Double-precision land's maximum cumulative water residual is 2.97e-9 mm
  through D's day 390, with native-precision checkpoint restoration.
  Century-scale closure and old-checkpoint conversion remain unqualified.
- Atmospheric water accounting is separate. v0.1.9 removes a frozen launch
  row that was a major source: in one fresh matched pair the days 5–30 residual
  rate changes from +0.01560 to -0.00081 kg/m²/day. This is a 30-day result.
  The older v0.1.8 D run retains the imprint and grows by 5.43 kg/m² over days
  90–384. Net-column process rain closes the scheme's column removal.
  Transport, filtering and coupling remain separate open terms.
- Forcing is concentration/AOD driven. The current scope excludes carbon-cycle,
  SO2 and aerosol-emissions modules. RRTMGP receives the configured 420 ppm CO2.
  Terrarium carries a default 380 ppm ambient input, but it is physically
  dormant in production: no photosynthesis model is constructed and the
  prescribed evapotranspiration has no CO2 term. A shared driver and assertion
  are still required before CO2-responsive vegetation is enabled.

## Decisions retained from development

- Both convection and large-scale precipitation use net column water removal
  to diagnose surface precipitation. Large-scale precipitation's inherited
  latent-heating and derivative approximations remain unchanged.
- SpeedyWeather's zero-adjacent-transport vertical-diffusion control is used.
  The tested metric bulk-Richardson replacement worsened the 30-day climate.
- A strong sigma-dependent convection profile improved radiation and water
  drift but worsened temperature drift. Production rejects this profile.
- No polar sponge, global depth correction, arbitrary ice-thickness cap or
  broad river-mouth diffusivity is applied.
- Atmosphere and ocean initial conditions are inserted directly. The atmosphere
  then receives the spectral consistency cleanup above. A multi-day incremental
  launch remains future work.

The internal development labbook and rejected experiment matrix are intentionally
not shipped in this repository.
