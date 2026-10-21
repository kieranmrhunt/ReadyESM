# History

Release changes and earlier development results. The
[experiment notes](docs/experiments.md) track current research. Detailed run
records are in the development labbook.

## Main branch after 0.1.10

### 2026-10-21

The monthly runoff-heat tests found soil overheating in both enabled runs.
The disabled controls passed. Keep `terrarium_runoff_enthalpy` off while the
soil heat-transport repair is tested. See the [feature guide](docs/bs08-runoff-enthalpy.md).

### 2026-10-20

Added optional runoff-heat routing. Water now carries its land-derived heat
through the runoff store to the ocean when `terrarium_runoff_enthalpy` is on.
The model records the transferred heat, storage and energy in transit, and
restores these quantities from checkpoints. Production keeps the option off.

Added prescribed methane and nitrous oxide concentrations for RRTMGP through
`forcing.ch4_ppb` and `forcing.n2o_ppb`. Omitted values retain 1900 and 340 ppb.
The two extensions have joint CPU restart and full-grid A100 tests. See the
[validation record](docs/validation-runoff-trace-gases.md).

Added optional ocean carbon exchange and emissions-driven atmospheric CO2.
The model writes DIC, alkalinity and carbon-budget histories and restores
them from checkpoints. Four full-grid GPU days, CPU and GPU output/restart
comparisons, and the combined CPU suite pass. Production keeps carbon off.
See the [carbon guide](docs/carbon.md) and
[validation record](docs/validation-carbon.md).

Completed the T15 and T31 speed table on one A100. The standard T31 setup
with the 1° ocean runs at 5.74 simulated years per day, after 46–53 minutes
of start-up. See the [timing details](docs/performance.md).

Added T63 and T127 measurements with the same 1° ocean. Rates are 2.99 and
1.07 simulated years per day. These timing runs use a wider land-data search
for one T63 island column and three T127 columns. Production is unchanged.

Completed four January C06 comparisons. Rain evaporation and downdrafts gave
small precipitation and upper-temperature gains, with moisture and wind
penalties. A 20-member screen now separates their parameter effects.
The runoff-heat extension passed its earlier whole-model cost comparison,
while its land-only timing exceeded the local limit. A later short coupled
pair measured 7.6% overhead. It is available as an option, with its monthly
climate comparison still running. The convection settings remain under assessment.
See the [experiment notes](docs/experiments.md).

Fixed standalone atmosphere checkpoints so they restore both transform-column
work arrays. The coupled restart path already copied these arrays. The new
regression, full CPU suite and two ERA5 restore probes pass. The same
regression exposes the omission on the old source. See the
[restart checks](docs/validation-atmosphere-restart.md).

Fixed an ordering error in the Oceananigans z-star grid update. Surface heights
are now copied before neighbouring grid scales are calculated. The scaling
formula and launch bounds are unchanged. Four full-grid GPU restarts and
their independent repeats match all 1680 state entries at five saved points.
The package regression test and full CPU suite pass. See the
[validation record](docs/validation-ocean-grid-ordering.md).

Added optional CPU multi-plume convection and convective momentum transport.
The schemes use Float64 plume calculations with Float32 atmospheric storage.
Momentum mixing returns dissipated kinetic energy as heat. Rain and snow are
accounted for separately. Both schemes are selected through the configuration
file. The full CPU suite, six coupled pilots and six fresh-process restarts
passed. Restored states and four further steps matched the uninterrupted runs
exactly. See the [validation record](docs/validation-finite-convection.md).

Fixed the timestep wrapper and Terrarium initialization callback for atmosphere
models without land. Production continues to use GPU Betts–Miller convection
and the v0.1.10 configuration.

Research during the same period tested the new convection settings against
January and July ERA5. The preferred combination improved moisture and wind
scores against the initial multi-plume setting. Against native convection,
moisture and upper- and mid-level temperature errors remained larger.
Regional wind results were mixed. Two 30-day rain-evaporation and downdraft comparisons
reduced late net downward radiation by 0.63–0.97 W/m².

The ocean fix first passed compact comparisons across 425 arrays before
the full-grid tests above. A cloud-transport prototype increased
sharp-cloud peak retention from 56–57% to 87–88% after a T15 rotation.
The cloud prototype remained in a research branch. Its results are recorded
in the [experiment notes](docs/experiments.md).

The uniform convection-depth cutoff was closed after all four tested settings
worsened temperature scores at each of the three tested pressure levels.

### 2026-10-19

Added convective snowfall to output and diagnostic history, with separate
liquid and frozen precipitation checks. Older checkpoints reconstruct the
new history on load. The full CPU suite passed. Matched eight-step runs
reproduced all 1,643 existing state fields and 208 existing output variables
exactly. See the [validation record](docs/validation-snow-accounting.md).

## Releases

### 0.1.10 — 2026-10-17

- Used the v0.1.9 model source and physics. Changed the production ocean launch
  state to the midpoint of the December 1992 and January 1993 ECCO monthly
  means. Sea-ice area and effective thickness use January.
- Added a preparation recipe and a public January downloader. The
  scripts verify source receipts, physical bounds and exact output hashes. The
  production runner verifies the prepared bundle again before construction.
- Replaced the legacy v0.1.9 EVEL inpainting cache with a fresh cache from the
  same source. This changes 478 raw or masked EVEL cells by at most
  `1.8e-5 m s-1`. The matched control used fresh caches, so the measured SST
  gain comes from calendar centring rather than this cache cleanup.
- Three matched 30-day pairs improve days 1–10 SST RMSE by 0.208 K, from
  1.237 K to 1.029 K. All paired gains lie between 0.207 and 0.210 K. The gain
  is about 76 times the frozen 0.0027 K native-control repeat envelope. Of the
  day-zero gain, 95% remains on days 8–10.
- Over days 1–10, Southern Ocean SST RMSE falls from 2.043 K to 1.568 K.
  Antarctic ice-area error improves by 24%. Antarctic 15% extent error worsens
  by 15%, and IIEE worsens by 1.1%. All pairs pass the preregistered Antarctic
  regrowth guard on ice area and 15% extent. IIEE is report-only.
- Ocean freshwater relative closure peaks at 0.00331 and stays below the 0.01
  stop. The launch transient decays to control levels by about day 8.
- The full CPU suite, the real ERA5 launch tests and the fresh-cache A100 gate
  pass. The GPU gate rebuilds all five ocean caches from the eight base files,
  reproduces the qualified construction hashes and completes a one-day
  production run.

### 0.1.9 — 2026-10-14

- Removed the extra `l = lmax + 1` row from the ERA5 atmosphere immediately
  before the first grid transform. ReadyESM previously left non-zero pressure,
  temperature and humidity coefficients in this row. SpeedyWeather truncates
  tendencies but never updates those coefficients, so the launch imprint stayed
  frozen for the whole integration. Production enables the repair. The generic
  configuration default remains off for compatibility.
- In a matched fresh 30-day pair, the days 5–30 atmospheric-water residual rate
  changes from +0.01560 to -0.00081 kg/m²/day, a 94.8% reduction in magnitude.
  The one-time launch change is -0.134 kg/m² total atmospheric mass and
  -0.00097 kg/m² water. A separate diagnostic cleanup of an old day-120
  checkpoint reduces its subsequent residual by 99.8%.
- The fix removes a spurious 80–90 hPa polar surface-pressure depression
  and exposes an existing cold polar bias.
  The change alters restart trajectories. Existing checkpoints keep their
  original spectra until the migration helper is applied explicitly.
- The fresh 30-day run passes the production-output safety gates. It has no
  tropical ice and lowers maximum area-equivalent ice thickness from 8.36 to
  4.90 m in this pair. Coastal pile-up, ocean hot cells, whole-system energy
  closure and long-run climate drift are the main follow-up problems. See the
  [validation record](docs/validation-v0.1.9.md).

### 0.1.8 — 2026-10-13

- Selected mixed liquid and ice saturation and net-column-conservative large-scale
  precipitation in production (the tested D combination). The generic defaults
  remain liquid saturation and upstream-observed precipitation for existing
  configurations. Convection stays at two hours. Liquid and ice cloud paths
  remain 75 and 25 g/m².
- In the matched day-360–384 window, net TOA uptake falls from 8.75 to
  3.49 W/m². D reaches day 390. Its liquid-saturation and upstream-observed
  control stops near day 385 after ocean salinity becomes negative. Each arm
  contains one trajectory.
- Cooling and drying increase. Atmospheric water closure remains positive:
  over days 90–384 the residual grows by 5.43 kg/m² in D versus 2.17 in the
  control (0.0185 versus 0.00739 kg/m²/day). It remains positive late in D.
  Ocean warming and coastal ice pile-up remain. The release keeps the existing
  derivatives, analytic ozone profile, clouds and convection.
- The day-390 maps and timelines replaced the README figures. This trajectory
  used v0.1.8 physics from day 90 and was the longest accepted run at release.
- Local and GitHub CPU tests pass, alongside a production-resolution GPU
  launch, independent-process restart and exact restored-state comparison.
  See the [validation record](docs/validation-v0.1.8.md).
- Checked ocean temperature and salinity bounds throughout the saved history.
  A recovered endpoint no longer hides earlier negative salinity or
  overheating. Physical limits and model physics are unchanged. Added 46
  regression checks. The full CPU suite passes.

### 0.1.7 — 2026-10-07

- Shortened production Betts–Miller convection adjustment from four hours to
  two. Added `atmosphere_convection_timescale_seconds`. Its generic default
  remains 14400 for existing configurations, while production selects 7200.
  Cloud paths, reference humidity and the net-column rainfall fix are unchanged.
- Matched native-Float64-land continuations from day 30 to day 120 reduce
  days 90–120 net TOA uptake from 9.34 to 7.27 W/m² (22.2%). Both runs pass
  exact restoration and retain land-water residuals near 1.1e-9 mm.
- Trade-offs: rain +10.4%, snow -32.2%, surface-air temperature -0.13 K.
  The smaller cumulative atmospheric water residual partly cancels inherited
  error. Ice pile-up and ocean hot cells remain.
- Updated README maps and timelines to the day-120 candidate. The first 30 days
  used four-hour convection. The subsequent 90 days used the release setting.

### 0.1.6 — 2026-10-05

- Used Float64 for the production Terrarium state, geometry and cumulative
  water accounting. Atmosphere and ocean remain Float32. In the coupled
  30-day test, maximum land-water residual is 2.8e-10 mm. Native Float64
  fresh-process restarts restore exactly. Old Float32 checkpoints retain their
  original precision. Added 45 land-precision and restart-history checks.
- Increased diagnostic liquid cloud path from 60 to 75 g/m². Ice path stayed
  at 25 g/m², with four-hour convection. Against a matched Float64-land control,
  days 20–30 net TOA uptake falls from 21.70 to 16.56 W/m² (23.7%).
- Clear-sky uptake changed little. Atmospheric water residual worsened, and
  the late mean of worst-cell ice thickness increased despite a slightly
  lower final value.
- Replaced the README maps and timelines with the tested 30-day trajectory:
  fresh Float64 land, cloud75 enabled at day 2, another restart at day 10.
  Production starts with cloud75 immediately. The configuration change resets
  the previous 120-day validation horizon.
- Faster convection and reduced ice-cloud path remain separate experiments.
  See the [validation record](docs/validation-v0.1.6.md).

### 0.1.5 — 2026-09-30

- Used consistent branch types in the sea-ice momentum and EVP stress kernels.
  The change removes five literal pointer expressions from the affected GPU
  specializations. Physical equations, parameters and stepping order are unchanged.
- Vendored the patched ClimaSeaIce 0.5.8 runtime and added CPU tests for active,
  marginal and absent ice with single and mixed numerical precision.
- Used the computed snow-rebasing type for its zero alternative in layered
  thermodynamics. This removes another literal pointer from offline compiled
  code. All 3,374,400 tested state and flux values matched. Added 456 thermal
  checks. The ice-volume helper kept its original arithmetic because broader
  type promotion introduced an artificial mass flux.
- The updated public CPU suite passes, including the existing land, runoff,
  precipitation and output regressions. GitHub CPU CI also passes on the
  release candidate. The final release adds documentation to that candidate.
- Production-state CPU comparisons are exact. Two earlier candidate GPU
  restarts pass, one using ordinary package loading. Both pass the unchanged
  diagnostics and independent exact restored-state checks.
- The earlier 48-hour attempt fails with CUDA700 in its first resumed step, before
  loading the ice-consolidation kernel. The initiating operation is unknown.
  The GPU fault remained after the momentum and stress change.
- The snow-updated candidate completes 184 resumed steps from model hour 2 to
  hour 48, with full diagnostics and exact restored-state comparison, in 59m03s.
  It uses ordinary package loading, without a runtime patch or diagnostic
  observer. Its A100 driver differs from the failed run, leaving patch causality
  and general GPU reliability unresolved. See the
  [validation record](docs/validation-v0.1.5.md).
- Replaced the README snapshot with hour-48 maps and added timeline plots.
- A separate process restores that hour-48 checkpoint and reaches hour 49,
  again passing full diagnostics and exact restored-state comparison (56m24s).
- Production keeps the existing land and air–ice surface-solver settings.
- Follow-up on 1 October: the same runtime reaches day 30 and then day 120
  through fresh-process restarts, passing normal diagnostics and exact
  restored-state checks. The README maps and timelines were updated to day 120.
  Follow-up work shifts to ocean drift, radiative imbalance and local ice
  pile-up.

### 0.1.4 — 2026-09-24

- Saved only recorded temperature and radiation samples after an interrupted
  run. Unused capacity can no longer appear as spurious zero-kelvin output.
- Kept recorded non-finite samples and clarified that the legacy
  `global_surface_temperature` field is lowest-layer air temperature.
- Those output changes leave model physics and production settings unchanged.
- Integrated runoff at each native Terrarium step before passing it to the
  ocean. The previous endpoint-rate approximation could lose runoff water.
- Routed runoff only to unique ocean cells at the tripolar join, and filled
  the copied boundary cells from their physical partners before applying the flux.
- Checkpoints now retain pending runoff and held discharge. Restoring an older
  checkpoint fails with an explicit incompatibility.
- Checked rainfall-component sums at native Float32 precision. Exporting rates
  as Float64 retains their original rounding error.
- Used the pressure from convection and condensation to check precipitation
  against their humidity tendencies. The previous check used later pressure.
- Added both scratch buffers and array layouts to the Fourier GPU graph
  cache key, and kept the captured arrays until the cache is cleared. A
  focused GPU test reproduced incorrect forward and inverse results with
  changed scratch storage. The repaired cache passes 35 regression checks,
  including cleanup while another CUDA stream still has work queued.
- Allowed configuration inspection before input downloads. Model construction
  still requires the input files. CPU regression tests can run without ERA5
  credentials, and GitHub Actions runs them from the pinned environment.
- Replaced GPU sparse exchange multiplication with a fixed-order kernel that
  assigns each destination row to one work item. The production grid weights
  reproduce shared-memory sanitizer reports in cuSPARSE 12.7.10. The replacement
  passes the same CPU-reference and repeatability checks with zero reports.
  A standalone regression passes 176 assertions under both initcheck and
  memcheck with zero errors, covering both directions and Float32/Float64
  edge cases.
- CPU CI, focused GPU regressions, a fresh production smoke and short coupled
  sanitizer checks pass. This maintenance release retains the first-step GPU
  illegal-access limitation already documented in v0.1.3. An ordinary restart
  of its fresh smoke still fails. The next release qualifies GPU restart.
  The [validation record](docs/validation-v0.1.4.md) includes that failure,
  passing checks and source hashes.
- Replaced the README's v0.1.1 figure with v0.1.4's surface
  fields at model hour two, with source hashes and an explicit short-run label.

### 0.1.3 — 2026-09-08

- Removed the redundant Darcy top-face override reported in [#1](https://github.com/kieranmrhunt/ReadyESM/issues/1).
  Terrarium already fills a zero-gradient pressure-head halo for flux boundary
  conditions. Prescribed pressure gradients are no longer forced to zero.
- Kept the separate porosity, evaporation and runoff-conservation fixes.
  Production physics settings are unchanged.
- Checked plain and masked soil and land models on CPU and GPU, including
  prescribed-gradient boundaries and water conservation. Stored land states
  match the previous default-boundary implementation in these tests.
- This maintenance release keeps the new snow, land-energy and coupled-startup
  work in experimental branches. The README continues to show the labelled
  v0.1.1 run. Version 0.1.3 reuses that integration.

### 0.1.2 — 2026-08-27

- Added `atmosphere_radiative_surface_albedo` to coupled NetCDF output. It is
  copied from the live RRTMGP solver state, so it records the value used by the
  radiation calculation rather than a later reconstruction.
- Added the solver-state source and checks for direct and diffuse
  agreement, spectral-band invariance and physical bounds.
- Passed the complete development regression and a real coupled
  RRTMGP-to-NetCDF smoke test.
- Model physics, configuration and parameters match `v0.1.1`. The new field
  helps diagnose the top-of-atmosphere radiation imbalance.

### 0.1.1 — 2026-08-27

- Fixed an ambiguous Terrarium 0.1.6 land-process dispatch in the exact T31/L27
  production path with prescribed ERA5 vegetation.
- The ambiguity could select Terrarium's vegetation-free forwarding method in
  place of ReadyESM's prescribed-vegetation evapotranspiration method. The fix
  makes the intended production method unambiguous.
- Added a regression that requires the production call to resolve to
  ReadyESM's evapotranspiration method.
- Added a release-matched 30-day output overview to the README.

### 0.1.0 — 2026-08-26

First published version. It contains the accepted T31/L27 atmosphere, 1°
tripolar 60-level ocean, dynamic sea ice and 16-layer Terrarium land model.
ERA5 and ECCO provide the initial state. RRTMGP uses prescribed 420 ppm CO2
and zero sulfate AOD in the production configuration.

The release also includes repaired atmospheric water accounting, direct ice
and snow volume transport, restart machinery, strict configuration contracts
and final-state diagnostics.

## Development chronology

### 4 August — first coupled stack

ReadyESM began around SpeedyWeather's moist primitive-equation atmosphere. The
first T15/L8 baseline used a slab ocean, thermodynamic ice and simplified
radiation, then gained:

- an RRTMGP adapter with explicit pressure, ordering and unit conversions.
- radiatively active prescribed CO2 and humidity-aware prescribed sulfate AOD.
- diagnostic all-sky cloud optics.
- Terrarium soil energy and water.
- Oceananigans and ClimaSeaIce through NumericalEarth and ClimaOcean.
- NetCDF diagnostics, Makie figures and Orchid GPU launchers.

A fixed-profile CO2 sweep gave a monotonic longwave response from 280 to
1120 ppm. AOD varied with humidity even when sulfate mass stayed fixed.
Recalibrating column mass against realised AOD closed that test.

Terrarium 0.1.4 was initially incompatible with SpeedyWeather 0.21.1's land
extension, so the early stack used 0.1.2. Renamed initialization keywords,
changed return types, an ambiguous `run!`, an unsynchronised callback period
and masked-land diagnostics were repaired before the first simultaneous
atmosphere–radiation–land–ocean–ice step passed.

### 5–6 August — GPU year and the first unacceptable climate

The initial one-degree ocean was unstable with a WENO momentum choice that
differed from ClimaOcean's long-running examples. Vector-invariant momentum
and WENO7 tracers gave stable ten-day ocean-and-ice-only and fully coupled runs.
A 365-day full-A100 integration then completed 35,040 coupled steps in about
two hours. The simulated climate had severe atmospheric cooling, tropical
sea ice, extreme coastal salinity and implausible ice growth.

Two definite coupling defects were fixed:

- SpeedyWeather precipitation was passed as depth flux to a NumericalEarth
  interface expecting mass flux, suppressing ocean precipitation by 1000.
- The atmosphere lacked live ClimaSeaIce concentration, so atmospheric
  reflection and surface absorption used inconsistent ice albedos.

The uniform soil-moisture output was also traced beyond plotting. Terrarium's
state and clock needed to advance through the coupled path, and diagnostics
needed the evolving internal state rather than the initializer.

### 7–10 August — deterministic GPU execution and code audit

Step-zero and step-one bisection isolated four independent reproducibility
defects:

- `RingGrids.similar(::Field)` could round-trip undefined memory when the
  element type already matched.
- an incompletely written RRTMGP net-flux work buffer.
- cuSPARSE `CSR_ALG1` admitted documented run-to-run variation.
- an atomic Float32 Legendre accumulation depended on arrival order despite a
  clean racecheck.

The Legendre transform moved to one work item per coefficient, eliminating the
unordered atomic sum without measurable cost. Runs then exposed reproducible
local climate biases.

A code audit added checks on configuration, dates and forcing sources. It
separated diagnostics from restart state, corrected the label of the old
"surface temperature" series, and excluded dry-cell fill values from ocean
extrema. The production clock was aligned with the January
1993 ERA5/ECCO state.

### 11–12 August — high-top tripolar model and conservation work

The atmosphere moved from L8 to L27 so radiation included a useful upper
troposphere and lower stratosphere. ERA5 was mapped into the T31 atmosphere.
ECCO supplied the 360 × 180 × 60 tripolar ocean, currents, free surface and
sea-ice state.

Two 730-day latitude–longitude development candidates, called V3 and V4 in the
labbook, informed later work. V3 was rejected when a shallow Indonesian
column exceeded 40 °C. V4 imposed a 40 m minimum on wet cells and reduced the
bias, but still crossed the unchanged gate for five days by up to 0.09 °C. The
global depth correction was later retired from the tripolar production route as
a coarse-resolution workaround rather than a general physical treatment.

The land gained ERA5 soil temperature and moisture, prescribed vegetation
fraction, LAI and roots, moisture-limited transpiration and conservative runoff
transfer.
A runoff and storage defect was repaired with porosity- and layer-thickness-aware
reconciliation. Daily land storage and flux ledgers replaced endpoint-only soil
moisture checks.

Atmospheric water accounting became prognostic and conservative. Betts–Miller
rainfall now diagnoses net column humidity removal, instead of adding every
drying layer while ignoring compensating moistening. This removed a genuine
water source. An inactive-cell mask also removed ECCO ice retained in dry
cells, which had contaminated coastal mechanics and early low-latitude maps.

Import latency fell substantially after Makie became lazy. Persistent batch
processes and shared precompiled depots were kept. Constructed-model caching
stays disabled because the GPU and package state is mutable.

### 13–16 August — restart exactness and Arctic ice attribution

ECCO `SIheff` was corrected to effective thickness. Geographic EVEL and NVEL
currents were rotated together into tripolar axes. Both were real defects, but
the Arctic pile-up persisted in matched 30-day tests and after the necessary
fold-vector halo corrections.

Force and transport diagnostics showed conservative dynamic convergence into
poorly connected coastal cells, followed by a pressure-supported near-stall.
The measured growth came from local convergence. Global ice volume and the
conditional-thickness display ceiling did not explain it. Increased EVP substeps, generic
immersed coastal drag, current rotation alone and broad mask changes all failed
their declared gates.

Exact fresh-process restart required more than the obvious prognostic arrays:

- lagged radiation and exchange fields.
- SpeedyWeather implicit operators and leapfrog state.
- transform scratch arrays read during continuation.
- Oceananigans closure halo storage.
- current Runge–Kutta tendencies.
- adaptive vertically implicit advection timesteps.

With these restored, a two-step checkpoint plus two-step fresh-process resume
was bitwise equal to its causally linked uninterrupted CPU control. Frozen-source
Slurm chains were added for GPU and segmented production use.

A one-column bathymetry alignment error at an original Arctic hotspot was
repaired, but the maximum relocated. A majority-area coastline, doubled EVP
iterations, generic coastal drag and replacement pressure were rejected. The
remaining issue was classified as missing coarse-coast, landfast-ice and
thickness-distribution physics at coarse coastal cells.

SpeedyWeather's original vertical-diffusion neighbour indices also made its
mixing effectively zero. A metric-aware adjacent-level operator was
mechanically sound but worsened the matched 30-day climate, so production kept
the zero-operator control.

### 17–18 August — atmospheric drift attribution

The high-top atmosphere reduced the radiative error, but positive TOA imbalance
and first-month cooling and moistening remained. Crossed offline replays separated
initial-state, cloud and clear-sky contributions. The adjustment was vertically
broad. Latent moistening offset much of the dry-temperature enthalpy loss but
left a residual in the full energy budget.

Several cloud changes failed their comparison tests: a physical-sigma diagnostic rule, a SimCloud transplant, altered fixed
water paths and the first conservative prognostic-condensate prototype. Some
improved individual radiation terms but worsened the coupled climate or failed
to reproduce the required cloud-radiation effects. Diagnostic clouds therefore
remained the baseline, with prognostic clouds retained as future work.

Restart testing found and repaired a CATKE closure-history omission.
Repeat-aware comparison envelopes separate tiny GPU reduction-order variation
from climate responses.

### 19–21 August — long tripolar evidence and river mouths

The strict tripolar chain completed two years, but produced an unacceptable
climate. It reproduced late coastal subsurface heating, a widening two-sided
salinity envelope and Arctic ice above 20 m area-equivalent thickness. Its
final atlas is archived with the rejected result.

The first salinity failure was localised to shallow, poorly connected runoff
receivers. Equal areal freshwater flux causes much greater fractional dilution
in a 10 m column than in a deep one. WOA23 showed why a global 20 psu floor was
invalid, so the audit moved to local climatological anomalies, non-negative
physical bounds and explicit receiver histories.

Column-volume routing, strong vertical and horizontal estuary mixing, weaker
receiver-centred mixing and salinity-only mixing were compared. Strong mixing
removed fresh pits but overmixed broadly. Weaker variants improved salinity
while worsening SST or dormant-mouth locality. Production kept its existing
closure. The tests established that salinity and trapped heat overlap spatially
without sharing one acceptable coefficient.

A coupled incremental-analysis launch was also built for atmosphere, land,
ocean and ice, including hidden component memory and conservation sources.

### 22–24 August — launch, process ledgers and fail-closed jobs

The five-day balanced launch crossed the old early failure horizon and reduced
insertion shock. The later climate bias persisted. The work exposed separate
land-state and skin-energy defects, which were fixed.

Layer-resolved humidity and temperature tendencies were added for convection,
large-scale condensation, radiation and surface exchange. This observer found
the Betts–Miller water source. A net-column correction to SpeedyWeather's
large-scale precipitation was then rejected at 30 days because it overcorrected
the climate. Production kept the upstream scheme and recorded it diagnostically.

Slurm orchestration became fail closed. Invalid dependency jobs were removed,
submitters gained `--kill-on-invalid-dep=yes`, source and configuration digests
became mandatory, and the scheduler cancels downstream analyses after failed
branches.
Driver stalls and zero-runtime submission errors are recorded separately from
scientific failures.

Other prototype tests ran alongside this work. Finite cloud memory smoothed
weather while the climate bias persisted. Five ice categories still needed a mechanical
redistribution scheme. A subgrid coastal graph
recovered channels, while the ice peak persisted with connectivity alone.

### 25 August — constrained successors

A loading-scaled river-mouth closure removed the 60-day fresh pit but failed its
global-locality guardrail. A sigma-dependent Betts–Miller humidity profile
improved radiation and column-water drift at 30 days but worsened temperature.
Both failed their declared gates.

Sea-ice transport violated represented-volume conservation under some
concentration changes. ReadyESM changed to direct conservative transport of
ice and snow volume. This removes the mass error. Coastal convergence remains
the next target.

The first full-GPU prognostic-cloud prototype ran with cloud state coupled
to radiation, but worsened the tested climate scores.

### 26 August — Terrarium 0.1.6 and the public repository

Terrarium moved from the patched 0.1.2 baseline to 0.1.6, with Oceananigans
0.110.15. The API and physics migration required the new state layout and
model-owned timestepper, separate boundary and auxiliary phases, changed
evapotranspiration and hydrology signatures, new stratigraphy fields and
physical water-flux units.

The compatibility baseline disables the new snow store and passes snowfall
once as liquid water equivalent. This avoids the temporary double count found
during migration. Existing albedo, emissivity, skin conductivity, prescribed
vegetation and porosity-aware water repairs were retained.

Fresh GPU compilation exposed a Terrarium generated-function helper defined too
late for Julia 1.12. Moving the unchanged helper before its generated method
allowed CPU and CUDA water ledgers to pass. A real GPU step that couples RRTMGP
and Terrarium also passed.
The clean accepted model was then published as `v0.1.0`.

The first release-exact T31/L27 documentation run exposed an ambiguity between
Terrarium's vegetation-free forwarding method and ReadyESM's prescribed-
vegetation method. The earlier small GPU test selected a different dispatch
path. The exact bridge and regression form `v0.1.1`.

## Naming note

Only `0.1.x` identifiers are GitHub releases. Capital V3/V4 were pre-release
ocean and bathymetry configurations. `localized_estuary_v1` through `v5` were
river-mouth experiments. Production uses none of them.

Current scientific limitations and active work are summarised in
[`NOTES.md`](NOTES.md).
