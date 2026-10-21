# Version history

Latest tagged release: [v0.1.11](https://github.com/kieranmrhunt/ReadyESM/releases/tag/v0.1.11),
21 October 2026.

See [known limitations](NOTES.md), [current experiments](docs/experiments.md)
and the [detailed development record](https://github.com/kieranmrhunt/ReadyESM/blob/9c6b388a18d0480d5e0145768fb7f9b72385811b/HISTORY.md)
for test results and earlier work.

## 0.1.11 — 2026-10-21

Production settings remain those of v0.1.10.

- Average land sensible and latent heat fluxes over each coupling interval.
  This makes exported fluxes match those used during the land substeps.
  The month test removes the regular-step land-water exchange mismatch to
  rounding precision. Moisture error improves slightly, while temperature
  and rainfall errors increase slightly. [Tests](docs/validation-v0.1.11.md).
- Fixed the order of ocean surface-height and vertical grid-scale updates.
  [Tests](docs/validation-ocean-grid-ordering.md).
- Fixed standalone atmosphere restart to restore both transform work arrays.
  Added restart checks for solar, methane and nitrous-oxide settings.
  [Atmosphere tests](docs/validation-atmosphere-restart.md).
- Added separate convective snowfall diagnostics and support for atmosphere
  models without land.
- Added optional CPU multi-plume convection and convective momentum transport.
  [Settings and limitations](docs/finite-convection.md).
- Added optional [carbon exchange and emissions-driven CO2](docs/carbon.md),
  prescribed CH4/N2O, and [seasonal Earth–Sun distance forcing](docs/solar-forcing.md).
- Added optional [runoff heat transport](docs/bs08-runoff-enthalpy.md).
  Keep this option off. Month-long tests found excessive soil heating.
- Added [speed measurements](docs/performance.md) from T15 to T127 and
  simplified the setup instructions.

## 0.1.10 — 2026-10-17

- Centred the ocean initial state on 1 January using December and January ECCO
  monthly means. Sea ice uses January.
- Added input preparation, download and verification scripts.
- Three matched pairs reduced days 1–10 SST error from 1.237 to 1.029 K.
  Antarctic ice-area error improved, while ice-extent error increased.
  Model physics stayed the same as v0.1.9.

## 0.1.9 — 2026-10-14

- Removed an extra spectral row from the ERA5 atmosphere at initialization.
  Its pressure, temperature and humidity coefficients had remained fixed
  throughout the run.
- Reduced the atmospheric-water residual rate by 94.8% in a matched month.
  Removed a spurious polar surface-pressure depression and exposed a cold bias.
- Existing checkpoints keep their original spectra unless explicitly migrated.
  [Validation](docs/validation-v0.1.9.md).

## 0.1.8 — 2026-10-13

- Selected mixed liquid/ice saturation and net-column large-scale precipitation.
- Reduced late net TOA heat uptake from 8.75 to 3.49 W/m² in the matched
  day-360–384 comparison. Cooling and drying increased. The candidate reached
  day 390, while its control stopped near day 385 with negative ocean salinity.
- Added checks of temperature and salinity throughout saved output.
  Updated the README figures to day 390.
  [Validation](docs/validation-v0.1.8.md).

## 0.1.7 — 2026-10-07

- Shortened production Betts–Miller convection adjustment from four hours to two.
- Reduced days 90–120 net TOA heat uptake from 9.34 to 7.27 W/m² in a matched
  continuation. Rain increased and snow decreased. Coastal ice and ocean
  temperature problems remained.
  [Validation](docs/validation-v0.1.7.md).

## 0.1.6 — 2026-10-05

- Changed Terrarium state and water accounting to Float64. Atmosphere and ocean
  storage stayed Float32. Land-water residuals fell to 2.8e-10 mm in the month test.
- Increased diagnostic liquid cloud path from 60 to 75 g/m².
  Days 20–30 net TOA heat uptake fell from 21.70 to 16.56 W/m².
  Atmospheric-water accounting and late peak ice thickness worsened.
  [Validation](docs/validation-v0.1.6.md).

## 0.1.5 — 2026-09-30

- Fixed mixed-type branches in sea-ice momentum, stress and snow thermodynamics.
  Physical equations and production settings stayed unchanged.
- Completed GPU integration and fresh-process restart tests. Follow-up runs
  reached days 30 and 120. The cause of earlier GPU illegal-access failures
  remained unresolved.
  [Validation](docs/validation-v0.1.5.md).

## 0.1.4 — 2026-09-24

- Saved only recorded output after interrupted runs, removing spurious 0 K samples.
- Accumulated runoff at each land step and routed it to unique tripolar ocean cells.
  Checkpoints gained pending runoff and held discharge.
- Fixed GPU transform cache keys and replaced sparse exchange multiplication
  with a fixed-order kernel.
- Added CPU continuous integration and GPU regression tests.
  Ordinary coupled GPU restart still failed at this release.
  [Validation](docs/validation-v0.1.4.md).

## 0.1.3 — 2026-09-08

- Removed the redundant Darcy boundary override reported in
  [issue #1](https://github.com/kieranmrhunt/ReadyESM/issues/1).
  Prescribed pressure gradients now reach the soil solver.
- Kept the porosity, evaporation and runoff-conservation repairs.
  Production settings stayed unchanged.

## 0.1.2 — 2026-08-27

- Added the surface albedo used by RRTMGP to coupled NetCDF output.
- Added checks of its bounds and direct/diffuse agreement.
  Physics and parameters matched v0.1.1.

## 0.1.1 — 2026-08-27

- Fixed Terrarium method selection for production runs with prescribed vegetation.
  The intended evapotranspiration method is now selected explicitly.
- Added a regression test and release-matched 30-day README figure.

## 0.1.0 — 2026-08-26

First public release: T31/L27 atmosphere, RRTMGP radiation, 1° tripolar ocean
with 60 levels, dynamic sea ice and 16-layer Terrarium land. ERA5 and ECCO
supply initial conditions. Production uses prescribed 420 ppm CO2 and zero
sulfate aerosol optical depth.

Included net-column convective rainfall, conservative ice and snow volume
transport, land-water repairs, checkpoints and coupled diagnostics.

## Before 0.1.0

August development added the coupled components, ERA5/ECCO initialization,
GPU reproducibility repairs and conservation checks. Early V3/V4 labels refer
to experimental ocean setups. The [archived chronology](https://github.com/kieranmrhunt/ReadyESM/blob/9c6b388a18d0480d5e0145768fb7f9b72385811b/HISTORY.md#development-chronology)
records these tests and the rejected candidates.
