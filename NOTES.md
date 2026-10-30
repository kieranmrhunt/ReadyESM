# Notes

## Current model

- SpeedyWeather 0.21.1, T31/L27, ERA5 initial atmosphere.
- Production removes the extra `l = lmax + 1` row from the ERA5 spectral
  launch before the first grid transform. The generic option remains off.
- Full-temperature vertical advection, ported from SpeedyWeather PR #1286.
- K-profile boundary-layer momentum mixing at scale 1. Heat and humidity
  mixing are off in production.
- Two-hour Betts–Miller convection with reference RH 0.5 above sigma 0.3,
  increasing to 0.7 at sigma 0.7, then constant below. Net-column rain diagnosis.
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

- The land model uses a wrong exponent in the law for how fast water moves
  through soil that is not full of water. In drier soil the predicted ability
  to carry water is 8.6 times too high at half saturation and about 520 times
  too high at 10 percent saturation. The effect on the climate has not been
  measured. A correction is being tested and is not in this release.
- Water vapour that leaves sea ice goes into the air, but the ice and snow do
  not lose that water. The energy budget also misses the heat of fusion for
  each kilogram. In a four-year run, the recorded vapour from sea ice was
  equivalent to about 2.8 mm of water per year spread across Earth. This
  measures the sampled flow, not how a corrected model would change the
  ice or climate. The correction is not in this release.
- In the sea-ice motion code, ice that is not moving still feels half of its
  full internal pressure. The correct value is zero. Where ice strength
  changes from place to place, this can push ice that should stay still. A
  correction has passed its own tests but is not in this release.
- When a run starts from ECCO ocean data, the starting ocean currents are not
  put exactly where the ocean model keeps them, including along the fold of
  the grid near the North Pole. A correction has passed its own tests, but
  its effect on a coupled run has not been tested for release.
- When falling snow melts in the air, the rain and snow code does not take
  the heat of melting from the air. This comes from the SpeedyWeather version
  used. Its effect on the full model has not been measured.
- The coupler turns the wind into the directions of the ocean grid. Near the
  North Pole this turn is less accurate. In a test of the same rotation
  with a flow of 1 m/s, the error was up to 1.4 percent of that speed at
  85 to 88 degrees N and up to 5.4 percent at 88 to 89 degrees N. Its
  effect on the climate has not been measured.
- In the standard output, ocean and sea-ice surface velocities are labelled
  east and north, but they follow the directions of the ocean grid. Near the
  Arctic, where those grid lines turn, the values are not east and north.
  A correction is being tested and is not in this release.
- Some atmospheric diagnostics are labelled with a time about one model
  step (15 minutes) away from the state they describe. This is small for
  monthly means but matters for six-hourly output and the timing of the
  water budget. A correction is being tested and is not in this release.
- Two smaller issues remain. In the sea-ice heat code, energy is not conserved
  in one case: when the ice area grows and the snow on it then runs out.
  When the model removes negative humidity, it also removes very small
  positive humidity, below about one ten-millionth of a kilogram of water
  per kilogram of air. The size of either effect on the climate has not been
  measured. Corrections for both are being tested.
- In one 30-day test for each season with v0.1.17 physics, upper-level
  storm energy at 250 hPa was 0.65 to 0.79 of ERA5 at model resolution.
  The sampling uncertainty has not been measured.
- Land albedo is fixed at 0.4. Over Antarctic land-bearing grid cells,
  the ERA5 forecast value is about 0.85 in January and July. The effect
  of this difference on the coupled climate has not been measured.
- Sea-ice albedo is fixed at 0.60 and does not change with snow cover.
  In May of the v0.1.15 long run, the model blended surface albedo
  north of 80 degrees N was 0.57 to 0.58, against 0.77 to 0.79 in ERA5.
- In the v0.1.15 long run, global snowfall was about half the ERA5
  value, while total precipitation was close. Large-scale condensate
  becomes snow only below minus 10 degrees C, and convective
  precipitation is always rain. The rule's share of the shortfall
  has not been measured.

- v0.1.15 improves January/July climate scores and low-level
  winds, but increases monthly global ice-area bias by 0.19 and 0.33 million
  km². Two January members exceed the day-31 ice-area change allowance.
  Local ice accumulation remains, with mixed changes in July maxima.
  A consistent heat return for momentum mixing remains to be implemented.
  Its energy effect at production resolution needs further testing.
  See the [comparison](docs/validation-v0.1.15.md).

- v0.1.14 improves the initialized January/July climate scores and SST, but
  strengthens low-level winds. Tropical January SST becomes slightly too cold.
  July ice thickness increases locally near 72.4°S, 170.5°E, with daily peaks
  of 5.7–6.0 m area-equivalent thickness across three runs. Global ice area
  and volume decrease. See the [release comparison](docs/validation-v0.1.14.md).
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
- Optional moist-mass coupling reduces dry-air drift in one January month.
  It retains a roughly 0.2 kg/m² launch offset. The enthalpy carried by
  the coupled water source is recorded but not yet applied to the energy
  equations. The option stays off by default pending longer coupled tests.
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
