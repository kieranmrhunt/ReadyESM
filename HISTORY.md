# Version history

Latest tagged version: [v0.1.18](https://github.com/kieranmrhunt/ReadyESM/tree/v0.1.18),
30 October 2026.

See [known limitations](NOTES.md), [current experiments](docs/experiments.md)
and the [detailed development record](https://github.com/kieranmrhunt/ReadyESM/blob/9c6b388a18d0480d5e0145768fb7f9b72385811b/HISTORY.md)
for test results and earlier work.

## 0.1.18 — 2026-10-30

- Add the energy from the wind and from surface cooling that drives ocean
  mixing near the surface. The standard ocean combines four or five mixing
  schemes. With four or more, the energy input to surface turbulence was lost.
  The model now adds it.
- Add optional output files for the release scores: three files of means over
  the whole run, and one file of six-hourly values. The option needs RRTMGP
  radiation and a run length that is a whole number of six-hour periods. A run
  with this option on cannot continue from a saved state.
- Restart files now record the model settings, the exact source commit and the
  package versions. A restart stops if any of them differ. The run name,
  output folder, run length, analysis start day and output options do not
  count. Any new commit, even one that changes only the documentation, stops a
  restart. Saved states from earlier versions lack this record. To continue
  from one, use `allow_identity_mismatch=true` when restoring it, or wrap the
  `run!` call in `with_dynamic_restart_identity_override`.
- After each step the model checks its atmosphere, ocean, sea-ice and land
  state. It stops and names the first field that is not a finite number or is
  outside a fixed limit, for example an ocean temperature above 40 °C.
- The README example now runs for 31 days, and the figure data are named by
  release files and hashes, not by local paths.
- GPU runs now ask for version 13.2 of the CUDA runtime. Use an NVIDIA driver
  that supports it.
- Over four pairs run on the same GPU host, v0.1.18 took 4 to 16 percent longer
  per model day than the matched control. The new step check accounted for
  less than 1 percent in a separate timing test.
- Six 31-day runs, three in January and three in July, were compared with
  matched runs using the v0.1.17 settings from the same observed starts. The
  controls used an earlier source commit than the published v0.1.17 tag. In
  January the change in error was within the spread of the matched runs. In
  July the error grew by slightly more than that spread, so the original
  seven-field score rejected this version. The nine-field score did not reach a
  decision. The separate release rule also rejected the change. The paired
  radiation and model-state checks passed. These results are reported without
  changing the scores or the rule outcome. Runs that start from observations
  do not show long-term climate skill.
- The new ocean mixing energy ran correctly on the GPU. In a simple six-hour
  test with wind over one ocean column, it deepened the mixed layer by 0.06 m.
  The target was 0.1 m. The [validation record](docs/validation-v0.1.18.md)
  describes this, and [NOTES.md](NOTES.md) lists other known limitations.

## 0.1.17 — 2026-10-28

- Make one part of the atmosphere's numerical smoothing weaker. The model
  smooths the divergent part of the wind, which is air that spreads out or
  flows together, to remove small-scale noise. This smoothing now acts over
  four hours instead of one hour. The smoothing of the other wind,
  temperature and humidity fields stays at four hours.
- The tests used three pairs of month-long January runs and three pairs of
  July runs. The two runs in each pair start from the same state and differ
  only in this setting. The release score compares seven weather fields with
  observations, relative to the CMIP6 climate models, an international set
  of comparison models. A lower score is better. In July the mean score
  falls from 1.00332 to 0.99175, and on average ReadyESM is better than 36
  of 60 CMIP6 models instead of 33. In January the mean score changes from
  0.82418 to 0.81848, which is within the spread of the repeat runs. One
  January pair became worse, from 0.79986 to 0.83093. In one July pair, the
  small error in the atmosphere's water budget rose from 0.009 to 0.016 kg
  per square metre.
- A run that restarts from a saved state now checks that the smoothing
  settings are the same as when the state was saved. It stops if they
  differ. Saved states from earlier versions do not contain these settings,
  so this version cannot restart from them.
- The change does not reach a separate target for stronger storms high in
  the atmosphere. The [test results](docs/validation-v0.1.17.md) give the
  details.

## 0.1.16 — 2026-10-27

- Enable seasonal Earth–Sun distance in the standard configuration.
  Other production settings and model equations stay unchanged.
- In three paired runs per season, January normalised error changes from
  0.89050 to 0.88447 and mean CMIP6 position from 53.67 to 54.67 of 59.
  July changes from 1.03108 to 1.03048 and from 31.00 to 31.67.
  January improves beyond the repeat spread. July stays within it.
- The two-season radiation checks, CPU suite and GPU restart checks pass.
  Test fixtures now match the production defaults and cover the ERA5 launch.
  Local ice accumulation remains open. See the [test results](docs/validation-v0.1.16.md).

## 0.1.15 — 2026-10-25

- Added vertical momentum mixing in the atmospheric boundary layer.
  Production uses the K-profile scheme at scale 1, with heat and humidity
  mixing off. Convection and cloud settings are unchanged.
- Three paired January and July runs reduce the v4 normalised error by
  11.9% and 9.5%. Both CMIP6 positions and all radiation checks improve.
  The large tropical low-level wind jump is much smaller.
- Global ice-area bias rises in both months. Local ice accumulation and
  the heat lost through momentum mixing remain open.
  [Tests](docs/validation-v0.1.15.md).

## 0.1.14 — 2026-10-25

- Fixed vertical temperature transport by porting SpeedyWeather PR #1286.
  Production now transports the background temperature profile as well as
  local variations. Convection and cloud settings are unchanged.
- Three paired runs in each of January and July improve both v4 climate
  scores, radiation, SST and global ice-area bias. Stronger low-level winds
  and local July ice thickening remain. [Tests](docs/validation-v0.1.14.md).
- Older production checkpoints need
  `atmosphere_vertical_advection: speedyweather_0_21_1_centered`.
  Start a fresh run to use the corrected transport.

## 0.1.13 — 2026-10-23

- Changed the production Betts–Miller reference humidity by height. It is
  0.5 above sigma 0.3, rises to 0.7 at sigma 0.7, and stays 0.7 below.
  Model code and the two-hour adjustment time are unchanged.
- Three fresh paired runs in each of January and July improved the v4 climate
  score and TOA radiation. The paired-input, conservation, restart, CPU and
  physical-safety checks passed. The high sea-ice-area bias remains.
  [Tests](docs/validation-v0.1.13.md).

## 0.1.12 — 2026-10-21

Production settings remain those of v0.1.11.

- Added an optional atmospheric moist-mass coupling for the Betts–Miller
  configuration. In one January month, dry-air change fell from
  −6.51 to −0.24 kg/m². The option defaults to off.
- Checked exact option-off output, fresh restarts, checkpoint compatibility,
  the full CPU suite and an A100 run. The launch offset and energy carried
  by water mass remain open. [Tests](docs/validation-v0.1.12.md).

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
