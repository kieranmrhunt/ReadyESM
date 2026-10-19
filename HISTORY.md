# History

This is a concise history of the accepted model releases and the main
experiments that shaped them. The internal labbook retains the job-by-job
record, including failed commands, frozen-source hashes and diagnostic
artefacts.

## Main branch after 0.1.10

### 2026-10-19

- Track convective snowfall separately in output and diagnostic history.
- Check liquid and frozen precipitation budgets separately. Reconstruct the
  added history when loading older checkpoints.
- Keep dynamics, forcing and surface fluxes unchanged. Full CPU tests pass.
  Matched eight-step production-grid runs preserve all 1,643 existing state
  fields and all 208 existing output variables exactly. See the
  [validation record](docs/validation-snow-accounting.md).

## Releases

### 0.1.10 — 2026-10-17

- Keep the v0.1.9 model source and physics. Change the production ocean launch
  state to the midpoint of the December 1992 and January 1993 ECCO monthly
  means. Keep January sea-ice area and effective thickness.
- Add a deterministic preparation recipe and a public January downloader. The
  scripts verify source receipts, physical bounds and exact output hashes. The
  production runner verifies the prepared bundle again before construction.
- Replace the legacy v0.1.9 EVEL inpainting cache with a fresh cache from the
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

- Remove the extra `l = lmax + 1` row from the ERA5 atmosphere immediately
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
- The fix also removes a spurious 80–90 hPa polar surface-pressure depression.
  This materially changes the climate and exposes an existing cold polar bias.
  The change alters restart trajectories. Existing checkpoints keep their
  original spectra until the migration helper is applied explicitly.
- The fresh 30-day run passes the production-output safety gates. It has no
  tropical ice and lowers maximum area-equivalent ice thickness from 8.36 to
  4.90 m in this pair. Coastal pile-up, ocean hot cells, whole-system energy
  closure and long-run climate drift are the main follow-up problems. See the
  [validation record](docs/validation-v0.1.9.md).

### 0.1.8 — 2026-10-13

- Select mixed liquid and ice saturation and net-column-conservative large-scale
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
- Follow-up work covers full-system energy closure and a repeat envelope for D.
- Local and GitHub CPU tests pass, alongside a production-resolution GPU
  launch, independent-process restart and exact restored-state comparison.
  See the [validation record](docs/validation-v0.1.8.md).
- Check ocean temperature and salinity bounds throughout the saved history.
  A recovered endpoint no longer hides earlier negative salinity or
  overheating. Physical limits and model physics are unchanged. Add 46
  regression checks. The full CPU suite passes.

### 0.1.7 — 2026-10-07

- Shorten production Betts–Miller convection adjustment from four hours to
  two. Add `atmosphere_convection_timescale_seconds`. Its generic default
  remains 14400 for existing configurations, while production selects 7200.
  Cloud paths, reference humidity and the net-column rainfall fix are unchanged.
- Matched native-Float64-land continuations from day 30 to day 120 reduce
  days 90–120 net TOA uptake from 9.34 to 7.27 W/m² (22.2%). Both runs pass
  exact restoration and retain land-water residuals near 1.1e-9 mm.
- Trade-offs: rain +10.4%, snow -32.2%, surface-air temperature -0.13 K.
  The smaller cumulative atmospheric water residual partly cancels inherited
  error. Ice pile-up and ocean hot cells remain.
- Update README maps and timelines to the day-120 candidate. The first 30 days
  used four-hour convection. The subsequent 90 days used the release setting.

### 0.1.6 — 2026-10-05

- Use Float64 for the production Terrarium state, geometry and cumulative
  water accounting. Atmosphere and ocean remain Float32. In the coupled
  30-day test, maximum land-water residual is 2.8e-10 mm. Native Float64
  fresh-process restarts restore exactly. Old Float32 checkpoints retain their
  original precision. Add 45 land-precision and restart-history checks.
- Increase diagnostic liquid cloud path from 60 to 75 g/m². Retain ice path
  at 25 g/m² and four-hour convection. Against a matched Float64-land control,
  days 20–30 net TOA uptake falls from 21.70 to 16.56 W/m² (23.7%).
- This release gives an incremental improvement. Clear-sky uptake is nearly
  unchanged, atmospheric water residual worsens, and the late mean of
  worst-cell ice thickness increases despite a slightly lower final value.
  Follow-up work covers ocean drift, coastal pile-up and long-run balance.
- Replace the README maps and timelines with the tested 30-day trajectory:
  fresh Float64 land, cloud75 enabled at day 2, another restart at day 10.
  Production starts with cloud75 immediately. The configuration change resets
  the previous 120-day validation horizon.
- Faster convection and reduced ice-cloud path remain separate experiments.
  See the [validation record](docs/validation-v0.1.6.md).

### 0.1.5 — 2026-09-30

- Use consistent branch types in the sea-ice momentum and EVP stress kernels.
  The change removes five literal pointer expressions from the affected GPU
  specializations. Physical equations, parameters and stepping order are unchanged.
- Vendor the patched ClimaSeaIce 0.5.8 runtime and add CPU tests for active,
  marginal and absent ice with single and mixed numerical precision.
- Use the computed snow-rebasing type for its zero alternative in layered
  thermodynamics. This removes another literal pointer from offline compiled
  code while preserving 3,374,400 tested state and flux values exactly. Add 456
  portable thermal checks. Keep the ice-volume helper unchanged to avoid
  artificial mass flux from broader type promotion.
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
- Replace the README snapshot with actual hour-48 maps and add timeline plots.
- A separate process restores that hour-48 checkpoint and reaches hour 49,
  again passing full diagnostics and exact restored-state comparison (56m24s).
- Production keeps the existing land and air–ice surface-solver settings.
- Follow-up on 1 October: the same runtime reaches day 30 and then day 120
  through fresh-process restarts, passing normal diagnostics and exact
  restored-state checks. Update the README maps and timelines to day 120.
  Follow-up work shifts to ocean drift, radiative imbalance and local ice
  pile-up.

### 0.1.4 — 2026-09-24

- Save only recorded temperature and radiation samples after an interrupted
  run. Unused capacity can no longer appear as spurious zero-kelvin output.
- Preserve genuinely non-finite samples and clarify that the legacy
  `global_surface_temperature` field is lowest-layer air temperature.
- Those output changes leave model physics and production settings unchanged.
- Integrate runoff at each native Terrarium step before passing it to the
  ocean. The previous endpoint-rate approximation could lose runoff water.
- Route only to unique ocean cells at the tripolar join, and fill the copied
  boundary cells from their physical partners before applying the flux.
- Checkpoints now retain pending runoff and held discharge. Restoring an older
  checkpoint fails with an explicit incompatibility.
- Check rainfall-component sums at native Float32 precision. Exporting rates
  as Float64 retains their original rounding error.
- Retain the pressure used by convection and condensation when checking
  precipitation against humidity tendencies, rather than using later pressure.
- Include both scratch buffers and array layouts in the Fourier GPU graph
  cache key, and retain the captured arrays until the cache is cleared. A
  focused GPU test reproduced incorrect forward and inverse results with
  changed scratch storage. The repaired cache passes 35 regression checks,
  including cleanup while another CUDA stream still has work queued.
- Allow configuration inspection before input downloads. Model construction
  still requires the input files. CPU regression tests can run without ERA5
  credentials, and GitHub Actions runs them from the pinned environment.
- Replace GPU sparse exchange multiplication with a fixed-order kernel that
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
  passing checks, source hashes and the scope of this incremental release.
- Replace the README's v0.1.1 figure with v0.1.4's actual surface
  fields at model hour two, with source hashes and an explicit short-run label.

### 0.1.3 — 2026-09-08

- Remove the redundant Darcy top-face override reported in [#1](https://github.com/kieranmrhunt/ReadyESM/issues/1).
  Terrarium already fills a zero-gradient pressure-head halo for flux boundary
  conditions. Prescribed pressure gradients are no longer forced to zero.
- Retain the separate porosity, evaporation and runoff-conservation fixes.
  Production physics settings are unchanged.
- Checked plain and masked soil and land models on CPU and GPU, including
  prescribed-gradient boundaries and water conservation. Stored land states
  match the previous default-boundary implementation in these tests.
- This maintenance release keeps the new snow, land-energy and coupled-startup
  work in experimental branches. The README continues to show the labelled
  v0.1.1 run. Version 0.1.3 reuses that integration.

### 0.1.1 to 0.1.2 at a glance

| | 0.1.1 | 0.1.2 |
|---|---|---|
| Model physics and parameters | Accepted 0.1.1 model | Unchanged |
| Land coupling | Correct Terrarium evapotranspiration method | Unchanged |
| Saved shortwave surface albedo | Unavailable | Exact albedo used by RRTMGP |
| Additional checks | — | Albedo provenance, consistency and bounds |
| Expected model trajectory from the same state | Baseline | Same as 0.1.1 |

Version `v0.1.1` was a model fix. Version `v0.1.2` changed output and
verification only. The extra albedo field diagnoses the outstanding
top-of-atmosphere radiation imbalance.

### 0.1.2 — 2026-08-27

- Added `atmosphere_radiative_surface_albedo` to coupled NetCDF output. It is
  copied from the live RRTMGP solver state, so it records the value used by the
  radiation calculation rather than a later reconstruction.
- Added solver-state provenance and fail-closed checks for direct and diffuse
  agreement, spectral-band invariance and physical bounds.
- Passed the complete development regression and a real coupled
  RRTMGP-to-NetCDF smoke test.
- Model physics, configuration and parameters match `v0.1.1`. This release adds
  observability for the outstanding TOA imbalance.

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
and final-state diagnostics. It is a research baseline for later calibration
and projection work.

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
two hours, proving the engineering path but producing an unacceptable climate:
severe atmospheric cooling, tropical sea ice, extreme coastal salinity and
implausible ice growth.

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
unordered atomic sum without measurable cost. The failure mode changed from
unpredictable NaNs to reproducible local physical biases.

A whole-code audit then tightened configuration validation, dates and forcing
provenance, distinguished diagnostics from restart state, corrected the meaning
of the old “surface temperature” series, and separated dry-fill contamination
from active-ocean extrema. The production clock was aligned with the January
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
The evidence pointed away from global ice-volume creation and the marginal
conditional-thickness display ceiling. Increased EVP substeps, generic
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
thickness-distribution physics rather than one bad cell.

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

Several plausible cloud changes were rejected rather than tuned into the
control: a physical-sigma diagnostic rule, a SimCloud transplant, altered fixed
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

Small blue-sky screens were also informative. Finite cloud memory smoothed
weather while the climate bias persisted. Five ice categories required real
redistribution physics rather than a new container. A subgrid coastal graph
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

The first full-GPU prognostic-cloud prototype was scientifically rejected. It
established a working cloud-state and radiation interface and showed that
condensate memory alone was insufficient.

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
