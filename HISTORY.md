# History

This is a curated history of the accepted model and the main experiments that
shaped it. The internal labbook contains the full job-by-job record, including
failed commands, frozen-source hashes and diagnostic artefacts.

## Releases

In plain terms, `v0.1.1` fixed a land-coupling dispatch error. `v0.1.2` uses
the same model equations, configuration and parameters as `v0.1.1`; it only
makes the shortwave surface albedo actually used by RRTMGP visible in the
saved diagnostics. A run started from the same state should therefore follow
the same model trajectory in both versions, while `v0.1.2` writes one extra
field and its provenance.

### 0.1.2 — 2026-08-27

- Added `atmosphere_radiative_surface_albedo` to coupled NetCDF output. This is
  read from the live RRTMGP solver state rather than reconstructed later.
- Added solver-state provenance and fail-closed checks for direct/diffuse
  agreement, spectral-band invariance and physical bounds.
- Passed the complete development regression and a real coupled
  RRTMGP-to-NetCDF smoke test.
- No model physics, configuration or parameter changed from `v0.1.1`. This is
  an observability release and does not repair the outstanding TOA imbalance.

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
ERA5 and ECCO provide the initial state; RRTMGP uses prescribed 420 ppm CO2
and zero sulfate AOD in the production configuration.

The release also includes repaired atmospheric water accounting, direct
ice/snow-volume transport, restart machinery, strict configuration contracts
and final-state diagnostics. It is a research baseline rather than a calibrated
projection model.

## Development chronology

### 4 August — first coupled stack

ReadyESM began around SpeedyWeather's moist primitive-equation atmosphere. The
first T15/L8 baseline used a slab ocean, thermodynamic ice and simplified
radiation, then gained:

- an RRTMGP adapter with explicit pressure, ordering and unit conversions;
- radiatively active prescribed CO2 and humidity-aware prescribed sulfate AOD;
- diagnostic all-sky cloud optics;
- Terrarium soil energy and water;
- Oceananigans and ClimaSeaIce through NumericalEarth/ClimaOcean; and
- NetCDF diagnostics, Makie figures and Orchid GPU launchers.

A fixed-profile CO2 sweep gave a monotonic longwave response from 280 to
1120 ppm. The aerosol test showed that fixed sulfate mass does not imply fixed
AOD as humidity changes; recalibrating column mass against realised AOD closed
that test.

Terrarium 0.1.4 was initially incompatible with SpeedyWeather 0.21.1's land
extension, so the early stack used 0.1.2. Renamed initialization keywords,
changed return types, an ambiguous `run!`, an unsynchronised callback period
and masked-land diagnostics were repaired before the first simultaneous
atmosphere–radiation–land–ocean–ice step passed.

### 5–6 August — GPU year and the first bad climate

The initial one-degree ocean was unstable with a WENO momentum choice that did
not match ClimaOcean's long-running examples. Vector-invariant momentum and
WENO7 tracers gave stable ten-day ocean/ice-only and fully coupled runs. A
365-day full-A100 integration then completed 35,040 coupled steps in about two
hours, proving the engineering path but producing an unacceptable climate:
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

Step-zero/step-one bisection isolated four independent reproducibility defects:

- `RingGrids.similar(::Field)` could round-trip undefined memory when the
  element type already matched;
- an RRTMGP net-flux work buffer was not fully written;
- cuSPARSE `CSR_ALG1` admitted documented run-to-run variation; and
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
troposphere and lower stratosphere. ERA5 was mapped into the T31 atmosphere;
ECCO supplied the 360 × 180 × 60 tripolar ocean, currents, free surface and
sea-ice state.

Two 730-day latitude–longitude candidates, called V3 and V4 in the labbook,
were useful but were not releases. V3 was rejected when a shallow Indonesian
column exceeded 40 °C. V4 imposed a 40 m minimum on wet cells and reduced the
bias, but still crossed the unchanged gate for five days by up to 0.09 °C. The
global depth correction was later retired from the tripolar production route as
a coarse-resolution workaround rather than a general physical treatment.

The land gained ERA5 soil temperature/moisture, prescribed vegetation fraction,
LAI and roots, moisture-limited transpiration and conservative runoff transfer.
A runoff/storage defect was repaired with porosity- and layer-thickness-aware
reconciliation. Daily land storage and flux ledgers replaced endpoint-only soil
moisture checks.

Atmospheric water accounting became prognostic and conservative. Betts–Miller
rainfall now diagnoses net column humidity removal, instead of adding every
drying layer while ignoring compensating moistening. This removed a genuine
water source. An inactive-cell mask also removed ECCO ice retained in dry
cells, which had contaminated coastal mechanics and early low-latitude maps.

Import latency fell substantially after Makie became lazy. Persistent batch
processes and shared precompiled depots were kept; caching complete constructed
models was rejected because too much mutable GPU/package state was unsafe.

### 13–16 August — restart exactness and Arctic ice attribution

ECCO `SIheff` was corrected to effective thickness, and geographic EVEL/NVEL
currents were rotated together into tripolar axes. Both were real defects, but
matched 30-day tests showed that neither removed the Arctic pile-up. Fold-vector
halo corrections were also necessary but did not cure it.

Force and transport diagnostics showed conservative dynamic convergence into
poorly connected coastal cells, followed by a pressure-supported near-stall.
It was not global ice-volume creation or the display-only marginal conditional-
thickness ceiling. Increased EVP substeps, generic immersed coastal drag,
current rotation alone and broad mask changes all failed their declared gates.

Exact fresh-process restart required more than the obvious prognostic arrays:

- lagged radiation and exchange fields;
- SpeedyWeather implicit operators and leapfrog state;
- transform scratch arrays read during continuation;
- Oceananigans closure halo storage;
- current Runge–Kutta tendencies; and
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
and first-month cooling/moistening remained. Crossed offline replays separated
initial-state, cloud and clear-sky contributions. The adjustment was vertically
broad; latent moistening offset much of the dry-temperature enthalpy loss but
did not close the full energy budget.

Several plausible cloud changes were rejected rather than tuned into the
control: a physical-sigma diagnostic rule, a SimCloud transplant, altered fixed
water paths and the first conservative prognostic-condensate prototype. Some
improved individual radiation terms but worsened the coupled climate or failed
to reproduce the required cloud-radiation effects. Diagnostic clouds therefore
remained the baseline, with prognostic clouds retained as future work.

Restart testing found and repaired a CATKE closure-history omission.
Repeat-aware comparison envelopes were introduced so tiny GPU reduction-order
variation could not be mistaken for a climate response.

### 19–21 August — long tripolar evidence and river mouths

The strict tripolar two-year chain reached its horizon but failed scientific
acceptance. It reproduced late coastal subsurface heating, a widening two-sided
salinity envelope and Arctic ice above 20 m area-equivalent thickness. Its
final atlas is negative evidence, not current model output.

The first salinity failure was localised to shallow, poorly connected runoff
receivers. Equal areal freshwater flux causes much greater fractional dilution
in a 10 m column than in a deep one. WOA23 showed why a global 20 psu floor was
invalid, so the audit moved to local climatological anomalies, non-negative
physical bounds and explicit receiver histories.

Column-volume routing, strong vertical/horizontal estuary mixing, weaker
receiver-centred mixing and salinity-only mixing were compared. Strong mixing
removed fresh pits but overmixed broadly; weaker variants improved salinity
while worsening SST or dormant-mouth locality. None entered production. The
tests established that salinity and trapped heat overlap spatially without
sharing one acceptable coefficient.

A coupled incremental-analysis launch was also built for atmosphere, land,
ocean and ice, including hidden component memory and conservation sources.

### 22–24 August — launch, process ledgers and fail-closed jobs

The five-day balanced launch crossed the old early failure horizon and reduced
insertion shock, but did not repair the later climate. It remained experimental.
The work did expose separate land-state and skin-energy defects, which were
fixed.

Layer-resolved humidity and temperature tendencies were added for convection,
large-scale condensation, radiation and surface exchange. This observer found
the Betts–Miller water source. A net-column correction to SpeedyWeather's
large-scale precipitation was then rejected at 30 days because it overcorrected
the climate; production observes the upstream scheme without applying it.

Slurm orchestration became fail closed. Invalid dependency jobs were removed,
submitters gained `--kill-on-invalid-dep=yes`, source/configuration digests
became mandatory, and downstream analyses could not run after failed branches.
Driver stalls and zero-runtime submission errors are recorded separately from
scientific failures.

Small blue-sky screens were also informative: finite cloud memory smoothed
weather but not the climate bias; five ice categories required real
redistribution physics rather than a new container; and a subgrid coastal graph
recovered channels but connectivity alone could not remove the ice peak.

### 25 August — constrained successors

A loading-scaled river-mouth closure removed the 60-day fresh pit but failed its
global-locality guardrail. A sigma-dependent Betts–Miller humidity profile
improved radiation and column-water drift at 30 days but worsened temperature.
Both remained experiments; thresholds were not relaxed after seeing results.

Sea-ice transport was found to advect concentration and conditional thickness
in a way that did not conserve represented volume under all concentration
changes. ReadyESM changed to direct conservative transport of ice and snow
volume. This removes that mass error but is not claimed to solve coastal
convergence.

The first full-GPU prognostic-cloud prototype was scientifically rejected. It
established a working cloud-state/radiation interface while showing that simply
adding condensate memory is not enough.

### 26 August — Terrarium 0.1.6 and the public repository

Terrarium moved from the patched 0.1.2 baseline to 0.1.6, with Oceananigans
0.110.15. This was an API and physics migration, not a drop-in dependency bump.
It required the new state layout and model-owned timestepper, separate boundary
and auxiliary phases, changed evapotranspiration/hydrology signatures, new
stratigraphy fields and physical water-flux units.

The new snow store remains disabled in the compatibility baseline; snowfall is
passed once as liquid water equivalent. This avoids the temporary double count
found during migration. Existing albedo, emissivity, skin conductivity,
prescribed vegetation and porosity-aware water repairs were retained.

Fresh GPU compilation exposed a Terrarium generated-function helper defined too
late for Julia 1.12. Moving the unchanged helper before its generated method
allowed CPU and CUDA water ledgers and a real RRTMGP/Terrarium GPU step to pass.
The clean accepted model was then published as `v0.1.0`.

The first release-exact T31/L27 documentation run exposed an ambiguity between
Terrarium's vegetation-free forwarding method and ReadyESM's prescribed-
vegetation method. The earlier small GPU test had not selected that exact full
production dispatch. The exact bridge and regression form `v0.1.1`.

## Naming note

Only `0.1.x` identifiers are GitHub releases. Capital V3/V4 were pre-release
ocean/bathymetry configurations. `localized_estuary_v1` through `v5` were
river-mouth experiments; none is active in `config/production.yml`.

Current scientific limitations and active work are summarised in
[`NOTES.md`](NOTES.md).
