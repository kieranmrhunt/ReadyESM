# Experiment notes

Updated 20 October 2026. Production is v0.1.10. Research candidates stay in
separate source snapshots until their results support a model update.
This page tracks the convection and blue-sky experiments.
Main now includes the accepted [snow-accounting update](validation-snow-accounting.md).
It also includes [optional CPU C08 and C07 convection](validation-finite-convection.md).
The [optional-land fix](validation-optional-land.md) passes the full CPU suite.
Main also includes the [ordered ocean grid update](validation-ocean-grid-ordering.md).
Production keeps its current Betts–Miller scheme.

## Convection portfolio

| Track | Status | Result or next test |
| --- | --- | --- |
| C01: mixed-layer parcel launch | Closed | The stronger launch reduces tropical deep-convection area too much. |
| C02: depth-limited shallow adjustment | Closed, uniform cutoff | Four cutoffs improve moisture but worsen temperature at every tested level and seed. Reopen for a regime-dependent or resolution-based mechanism. |
| C03: surface-forced shallow closure | Closed | Moisture improves, but the intended storm-track response and spatial mechanism fail their tests. |
| C04: elevated and slantwise convection | Closed at T31 | Eligible columns are rare at this resolution. Reopen for a higher-resolution or event-based test. |
| C05: non-equilibrium deep closure | Closed | The tested memory produces little change in convective timing. |
| C06: rain evaporation and cold pools | Active | Both 30-day pairs retain moistening and reduced net downward radiation. Coupled source and restart comparisons pass. Four January skill runs are submitted. Integration with the latest ocean fix is under test. |
| C07: convective momentum transport | Integrated option, assessment active | Coupled and exact CPU restart tests pass. July confirms global wind gains against C08, with regional trade-offs. |
| C08: deterministic multi-plume convection | Integrated option, assessment active | Coupled and exact CPU restart tests pass. Case 12 improves moisture against its original setting in January and July. |
| C09: stochastic plume population | Closed | Eight seeds show no resolved change in the tested mean heating. |
| C10: shallow eddy-diffusivity/mass-flux transport | Active | Positive grid transport includes phase exchange and heating. Wider owned storage reduces heat-budget error. Centered stepping and wider pressure coefficients reduce the largest six-hour mass error to 0.00309 kg/m². L27 and production-timestep tests follow. |
| C11: saturated-parcel trigger census | Closed | The affected columns contribute too little rain to justify this route. |

## Latest results

### C02

The 0.7 shallow-depth setting has entered a date-matched assessment. All 45
current-model column checks pass. A zero limit reproduces all 4,608 columns
exactly in both precisions. Three two-day ERA5 host pilots also pass, with
exact native/disabled reports and maps and zero water-ledger failures.

The monthly sampler passes 42 unit checks and three passive host pilots.
Every original report and saved map remains unchanged. It collects hourly
states, integrates rain plus snow over every step, and excludes unbracketed
pressure-level temperatures. Complete January coverage requires 744 hourly
states and 2,976 precipitation intervals.

All nine 31-day runs pass: native, zero limit and 0.7 for three matched seeds.
Their complete-run reader verifies every day and exact native/zero agreement.
The active setting reduces mean column water by 1.79–1.82 kg/m² and raises
precipitation by 0.054–0.069 mm/day. January ERA5 comparisons reduce column-water
RMSE by 0.820–0.870 kg/m² and precipitation RMSE by 0.036–0.078 mm/day.
Temperature RMSE increases by 0.274–0.296 K at 200 hPa, 1.197–1.255 K at 500 hPa
and 0.548–0.623 K at 850 hPa. Keep this setting separate from production.
The weaker-cutoff sweep has now completed. All fifteen new runs pass, and
the repeated off/0.7 controls match their archived fields. Independent checks
reproduce 432 scores, 360 paired differences and every common mask.

Cutoffs 0.4, 0.5, 0.6 and 0.7 worsen global temperature RMSE at 200, 500 and
850 hPa in all three seeds. All twelve moisture comparisons improve. Rain
error rises at 0.4 and 0.5, has mixed signs at 0.6, and falls at 0.7. The weaker
settings reduce the mid-level temperature penalty but lose the consistent rain
gain. Close this uniform-cutoff form for the current T31 setup. Keep the source
and results. Reopen for a physical change to the depth rule or new evidence at
another resolution.

### C06

All five one-day coupled pilots completed. Native and disabled controls pass.
The three active cases pass their internal conservation checks but fail the
rainfall diagnostic and the source-coverage screen.

The rainfall diagnostic omitted C06's separate humidity tendency. Adding it
reduces the worst saved-column discrepancy from 4–8e-5 to about 3e-9 kg/m²/s.
Every final column then meets the unchanged 3e-6 tolerance. The corrected
collector exports base, C06 and total tendencies. Thirteen component tests
pass. Native, disabled and strongest active GPU pilots all pass the collector,
file output and full validator. The active case has zero internal physical
failures. Its dry-parcel fraction is 3.16%, so the inherited 1.25% screening
flag remains failed.

The source review finds 14 cases with heating or drying below the diagnosed
cloud interval. They carry at most 0.01485% of sampled rain. C06 leaves all
unmapped rain unchanged. The larger dry-parcel class is tracked separately.
The earlier 1.25% coverage limit came from a project screening choice.

The failed screen remains archived. A separate ten-day assessment retains
the physical checks, reports both coverage flags, and checks that unmapped
rain passes through without evaporation. It compares native, middle and
strongest settings with two repeats each. The diagnostic pilots have passed
and the ten-day runs have started. The first two attempts reached day ten but
their report builder rejected missing launch diagnostics. A separate repair
keeps those values in NetCDF and checks every later sample. The replacement
assessment now passes in all six members, including reports and restart files.
Independent checks cover 84 output histories and 28 paired metrics. Over days
5–10, the stronger setting reduces all-sky net downward radiation by 0.684 and
0.708 W/m² in the two repeats. Column water rises by 0.0665 and 0.0703 kg/m².
Clear-sky net downward radiation rises by 0.0972 and 0.0805 W/m². Source-coverage
flags remain failed. Matched 30-day runs will test persistence.

Both 30-day pairs pass physical checks and production diagnostics. The reader
verifies all 56 output histories and the paired initial atmosphere/ocean-tracer
hashes. Over days 20–30, column water rises by 0.0462 and 0.1177 kg/m².
Net downward TOA radiation falls by 0.6271 and 0.9743 W/m² against native.
Clear-sky net downward radiation rises by 0.0672 and 0.0284 W/m².
Both aggregate dry-parcel fractions are below 1.25%. Source placement still
fails, and unmapped rain passes through unchanged.

The current-main candidate passes 141 focused assertions and the full CPU
suite. New checks reject unsupported C06/C07/C08 combinations and invalid
parameters. Coupled source comparisons and fresh-process restarts now pass.
Both two-day sampler pilots complete. Their first reader rejects the
Float32 model weights against a Float64 sum tolerance. An independent
quadrature reconstruction matches every saved weight exactly. Explicit
reader-only normalization then passes all original coverage, units, masks,
precipitation and input-passivity checks.

The current-main native reference finishes eight steps, then its report
builder fails because a snow-history variable hides a Julia function.
A corrected runner reuses both verified reference integrations. The native
report is explicitly recovered from execution evidence. The strong reference
report is unchanged. Twelve report-fragment tests pass. Both candidate
integrations and fresh-process restarts complete.

The native comparison reader then rejects four new diagnostic variables.
The two base-convection fields match the existing total tendencies exactly.
The two C06 fields are zero. A separate reader checks these identities,
units, dimensions and data types. All 1,645 native state entries and all
210 shared output variables match. Strong-case source checks match
1,994 shared entries, with two added snow-history entries. Restarts match
all state entries and all 214 output variables in both cases.
The original failures remain in the run records.

Four replacement 31-day native/strong runs will compare January ERA5 skill
and execution-repeat spread. They retain the frozen C06 candidate.
A separate candidate applies the same C06 code to the latest main branch,
including the accepted ocean-ordering fix. Its 518 focused CPU assertions
and full CPU suite pass. Coupled compatibility and restart tests are running.

### C08 and C07

The revised C08 closure passes eight coupled steps with the ocean, sea ice,
land and radiation components. All snow and conservation checks pass. The
disabled version reproduces the original C08 endpoint exactly. The active
case completes 101,396 closure evaluations with zero rejected probes.

Packed spectral arrays required a repair to the restart checker. Its 22
regression checks pass. Four full-state restart tests pass, with exact
checkpoint and continued states. The active one-day coupled run also passes.
Day-one restart extensions remain running. Four full-ESM pilots pass for finite C08,
zero-strength C07, momentum transport, and momentum transport with heat return.

The clean native, C08 centre, case 3 and case 12 runs now complete one coupled
day. All four reproduce their qualified short-run prefixes. State bounds,
closure checks, snow histories and production diagnostics pass. A fresh
reader run reproduces every saved report value.

Both day-one restart cases reproduce the step-96 checkpoint and steps 97–100
exactly. C08 also matches all 1,816 entries in the older reference.
The native reference comparison initially fails because its older wrapper
contains nine entries for an unused snow diagnostic. A separate audit verifies
the zero values, shape and grid metadata. All 1,645 shared native entries
match exactly. The archive retains the original failed readout.
Six current-main C07/C08 parameter comparisons continue.

The atmosphere screen starts from ERA5 on 1 January 1993. It uses T31/L27,
prescribed seasonal SST and simple radiation. Seven configurations cross
three perturbed initial states for 30 days. The original-C08 controls finish
their trajectories but fail rare column enthalpy checks. Their maximum
relative residuals are 0.023–0.031%, above the existing 0.02% limit. The
comparison retains these failure flags.

All nine revised-closure cases have now finished. The centre setting at seed
707 passes every check. The other eight have rare enthalpy failures. The
faster setting at seed 808 also has one host-water failure. All closure-search
checks pass, with zero rejected probes. These residuals concern the returned
exchange and its host representation. All 435 saved near-limit exchanges
replay exactly. Wider ledger reductions retain the residual. Fixed-mass-flux
tests locate error in the rounded plume profiles, their averaging and the
environmental exchange. Fallout changes the sampled residuals by at most
4.1e-9 W/m². The separate host-water failure is about 1.3e-12 kg/m²/s.

A mixed-precision prototype calculates the full plume and closure in Float64,
then returns Float32 tendencies. All 435 captured columns pass, with maximum
host energy ratio 3.58e-8 against the unchanged 2e-4 limit. All six two-day
host pilots also pass. Disabled reports and daily maps match exactly.
The centre case changes second-day column water by +0.000430 kg/m² and
temperature by -0.000116 K. It takes 233 seconds versus 259 seconds for the
existing solver in this single matched host test. The earlier selected-column
benchmark cost 3.42 times as much, so retain both timing results.

All nine mixed-precision thirty-day runs pass, with zero captured water or
energy failures. Three fresh Float32 controls reproduce the earlier reports,
maps and failure flags exactly. The precision change has little effect on
global means. Large cold and dry differences from native convection remain.
The full 30-day parameter grid now passes all 84 members, including three native
controls. Lower rain conversion and greater condensate retention improve water
and temperature relative to the centre case. Four selected settings also pass
two-hour coupled tests. Their maximum relative host-water residual is 1.94e-7.
All 12 selected January runs pass, with exact 30-day prefixes and centre controls.
Their January ERA5 scores now favour case 12 as the next coupled parameter test.
Against the original C08 setting, it reduces column-water RMSE by about
0.927 kg/m², 500 hPa temperature RMSE by 0.119–0.131 K and 850 hPa temperature
RMSE by 0.087–0.090 K. Rainfall RMSE improves by 0.001–0.004 mm/day.
The 200 hPa temperature RMSE increases slightly, by 0.002–0.009 K.
Case 12 halves precipitation conversion and doubles the condensate threshold,
while retaining the two-hour relaxation. Keep the one-hour case 3 as a warmer
contrast. The four-hour case 21 has larger temperature errors.

Against native convection, case 12 still increases column-water RMSE by
0.375–0.383 kg/m² and 500 hPa temperature RMSE by 1.295–1.312 K. Rainfall RMSE
falls by 0.793–0.804 mm/day. Longer coupled tests must assess these trade-offs.

The clean C08 package removes the experimental source wrappers. It passes 4,350
solver comparisons and 1,740 host comparisons exactly. This CPU-only option uses
the normal model configuration and leaves native defaults unchanged. The full
CPU suite passes. All six clean coupled pilots now pass for native, C08,
zero-strength C07 and active C07 at three settings. Native and C08 reproduce
their reference states and all 210 output variables. Zero-strength C07 matches
all 1,816 common state entries at six saved stages. Active C07 passes its
water, energy, momentum and snow checks. All six fresh-process clean restarts
now pass. Each restored state and four following steps match the uninterrupted
run exactly across 1,645–1,987 entries. Both schemes are available on main as
optional CPU configurations. Production defaults remain unchanged.

All four research day-one restarts also pass. Their complete states match
at the checkpoint and four following steps. These compare 1,654–1,816 entries
per member. The clean implementations have now passed their own restart tests.

All six combined C07/C08 cases pass their physical checks. Days 10–30 mean
column water is 16.535 and 16.440 kg/m² for the centre and faster settings,
versus 22.996 for native convection. Mean atmospheric temperature is 249.240
and 249.362 K, versus 251.256 K. These large differences need observation-based
assessment. The completed comparison retains every earlier failure flag.

The first mixed-precision C07 pilots exposed a profile mismatch. C08 selected
plumes in Float64, while C07 rebuilt them in Float32. Two exact failure replays
show one-versus-two and two-versus-three active plumes. Both the momentum and
heat paths now use C08's selection precision. The component test passes 6,961
checks and all eight paired two-day pilots pass, including exact controls.
All 15 January runs compare native, C08 and C07 across three seeds and pass
the complete-run checks. Their reader initially loaded a shared helper into
the wrong Julia module. The repaired reader passes without rerunning the models.
Wind output at 200, 500 and 850 hPa and the ERA5 wind reference are ready.
Both the 15-member wind scorer and 33-member parameter scorer have completed.
Independent interpolation and weighted sums reproduce all 1,692 score rows and
1,872 paired differences, including every common mask.

C07 centre reduces global wind-vector RMSE relative to C08 centre by
0.539–0.544 m/s at 200 hPa, 0.227–0.236 m/s at 500 hPa and 0.220–0.234 m/s at
850 hPa. Temperature RMSE improves at all three levels. Column-water RMSE rises
by 0.072–0.074 kg/m² and rainfall RMSE by 0.037–0.039 mm/day.
Tropical 850 hPa wind improves by 1.037–1.046 m/s, while extratropical wind
worsens by 0.150–0.163 m/s. Tropical 500 hPa temperature RMSE rises by 0.326 K.
Test momentum transport with the improved moisture parameters next.

All observed-skill ranges above cover three perturbed initial states for January
1993, with prescribed seasonal SST and simple radiation. Selected parameters
need a separate validation period. Production defaults remain unchanged.

July validation now passes all 21 members and the physical and ERA5 readers.
Independent interpolation and weighted sums verify 1,260 score rows and
2,340 paired differences. Keep these results separate from January.

Adding C07 to case 12 reduces global July RMSE against matched C08:

| Variable | RMSE reduction across three seeds |
| --- | ---: |
| Column water | 0.136–0.142 kg/m² |
| Rainfall | 0.052–0.054 mm/day |
| Temperature at 500 hPa | 0.412–0.447 K |
| Wind at 200 hPa | 0.697–0.752 m/s |
| Wind at 500 hPa | 0.339–0.406 m/s |
| Wind at 850 hPa | 0.378–0.439 m/s |

Temperature at 200 and 850 hPa also improves. Against native convection,
combined case 12 still has worse moisture and upper/mid-level temperature
errors, despite better rainfall and wind scores. Extratropical 200 hPa wind
worsens against C08. Case 12 remains the preferred coupled test, with case 3
as contrast. All six one-day coupled cases and their reader now pass,
including native, plain C08 and zero-strength C07 controls. Every case
reproduces its qualified eight-step prefix. Zero-strength C07 matches
plain C08 in all 1816 state entries at three saved boundaries and all
210 output variables. The preferred case 12 has no failed closure,
momentum or energy-return checks.

### C10

The stronger shallow-transport settings improve some global means but cause
local temperatures near 358 K in two days. Saved-column tests locate direct
convection increments as large as 77.6 K in one effective step.

Recomputing transport over shorter internal steps reduces those increments.
The host receives the divergence of the time-averaged water and energy fluxes.
All nine host cases and their comparison pass the original state and
conservation checks. The native and single-step controls reproduce earlier runs.

| Updraft area fraction | One-step maximum | Eight-substep maximum | Sixteen-substep maximum | Thirty-two-substep maximum |
| --- | ---: | ---: | ---: | ---: |
| 0.04 | 358.50 K | 324.16 K | 324.18 K | 324.22 K |
| 0.08 | 358.43 K | 324.28 K | 324.33 K | 324.36 K |

Both cases use 1000 m entrainment length and buoyancy-defined plume height.
The maxima cover all grid points, layers and completed steps. Their locations
and times can differ. Global means use 24 hourly samples from the second day.
The native control reaches 327.43 K. These tests measure numerical behaviour
and model sensitivity. Observational climate skill remains to assess.

Eight substeps change second-day global mean temperature by −0.0018 K at
area fraction 0.04 and −0.0069 K at 0.08 relative to one step. Rainfall shows
more internal-step sensitivity than the temperature maximum.

All six ten-day follow-ups pass and reproduce their two-day prefixes. The
temperature peak reaches 338.62–339.05 K. Within each updraft-area setting,
the range across 8, 16 and 32 substeps is less than 0.26 K. Over days 5–10,
mean temperature is 1.55–1.60 K colder than native and column water is
0.88–1.06 kg/m² lower. These are differences from the native model. The peak
is at 131.25°E, 20.41°S in the lowest layer. Passive replays will separate the
local convection, surface-flux and other atmospheric contributions.
All six eight-step full-model pilots pass. These retain the actual ocean and
land surface fluxes. Disabled prognostic fields, surface memories and scientific
arrays match native exactly. All six fresh-process restart tests pass.
The complete coupled state and all three forcing memories match at step four
and each continued step through step eight. Native has 1,681 compared entries,
and each experimental arm has 1,771. Active cases have zero physical failures.

A cloud-bearing prototype now carries liquid, ice and their latent energy.
Stored-column tests pass, followed by 432 live phase tests and 16 fresh-process
restart checks. Grid water still becomes negative after spectral projection.

A global smoothing repair was rejected. It cannot preserve each phase's
inventory when that inventory is already negative. In the strongest admissible
case it also removes over 99.9% of a cloud layer's horizontal variation.

A passive test of the native time filter passes 552 assertions and reproduces
the original default-filter states exactly. At step three, the filter turns
a near-zero unfiltered liquid coefficient into a negative value in both
precisions. With the filter off, this effect disappears but negative grid
water remains. The next prototype must handle phase exchange and time history
together, while conserving total water and energy.

A joint spectral-water projection passes 12 of 16 saved-state cases.
Four dry Float64 cases reach the iteration limit with tiny negative phases.
All cases pass water and energy checks. The 12 passing states also pass
physical phase equilibrium, which evaporates their warm, undersaturated cloud.

A separate balanced-cloud screen tests both precisions, three temperatures
and smooth or cloud-edge fields. Smooth Float32 clouds retain 99.9877–99.9974%
of condensate. At cloud edges, the final spectral transform creates negative
phases again. Water and energy budgets still close.

The diagnostic audit accounts for temperature and pressure roundoff in
saturation humidity. All six smooth-cloud cases now pass. Every physical
summary and equilibrium residual matches the original test exactly.
All six cloud-edge cases still fail. The archive keeps both sets of verdicts.
Cloud-edge positivity and consistent live time histories remain the next
physical tasks.

A conserved-variable prototype now passes projection inventories and grid
phase diagnosis in all 24 cases. It covers both precisions, three temperatures,
smooth and edge clouds, and uniform or variable pressure. The required native
temperature/vapour handoff passes nine cases. All 12 cloud-edge handoffs fail,
with local total-water changes up to 2.15e-6 kg/kg. The next adapter must
preserve local water and energy through that handoff and both stored time levels.

A separate grid-storage test now passes 28 of 32 cases. Keeping vapour,
liquid and ice on the grid passes all 16 cases for 32 cycles, with positive
owned water and conserved phase-exchange and layer budgets. Saturated Float32
clouds retain at least 99.857% of their initial mass.

Cloud-only grid storage passes 12 of 16 cases. Four dry cloud-edge cases
develop negative vapour after the native spectral transform. The projected
vapour also remains negative in the all-water arm, while its owned water
stays positive.

Shared-mass grid transport now passes 32 prescribed-flow cases and exact
fresh-process restarts. It moves air, vapour, liquid, ice and thermal content
through the same faces. The first-order method spreads sharp clouds too much.

Limited flux correction improves return error against donor Euler and donor
RK2 in all 16 full-rotation cases. Both new 32-case arms preserve nonnegative
water, conserved inventories and temperature bounds. All 64 restarts match
exactly. At T15, sharp-cloud peak retention rises from 56–57% to 87–88%,
and relative return error falls from 63–64% to 14–15%. Smooth-cloud peaks
retain about 98.9%.

The moving-cloud test adds phase equilibrium and prescribed heating or cooling.
All 128 cases and their fresh-process restarts pass the numerical checks.
Float32 storage still leaves a heat residual of 6498.5 J/m² in the worst
case. A matched repeat promotes the original extensive values to Float64
before stepping. Its largest heat residual among Float32-launched cases is
0.0000444 J/m². All 64 Float64 controls remain exact, as do all 128 restarts.

Read-only Float32 atmospheric views have up to 126.34 J/m² instantaneous
energy quantization error. These views never overwrite the owned state in
the test.

The native-dynamics sweep completes 24 cases and exact halfway restarts.
It checks pressure-derived air mass and both stored time levels. The diagnostic
leaves the selected state arrays unchanged. In the Float64 T7 rotating-flow
case, mean air mass rises by 0.989901 kg/m² over 3840 seconds while mean log
pressure changes by only 3.20e-14. Float32 exponentiation rounding is too
small to explain the gain. Resting Float64 controls remain steady to rounding.

The pinned core zeros the mean log-pressure tendency in two places.
The eight-member attribution test completes all 96 trajectories and exact
restarts. Removing either assignment alone leaves the diagnostic reports
unchanged. Removing both reduces absolute mass change in every moving-flow
case: by 99.607–99.985% in Float64 and 14.658–99.687% in Float32. In the
Float64 T7 rotating case, +0.989901 kg/m² becomes −0.000483958 kg/m².
Float32 retains a larger residual.

Four members stopped before construction because Julia dedented the method
lookup string. The repaired lookup passes seven regression assertions.
Only those four members were resubmitted.

A 36-member follow-up completes six hours at T7/T15/T31 and three timestep
settings in both precisions. All 216 trajectories and exact halfway restarts
pass. Compatible cases reproduce the earlier 64-step records. Removing
both assignments reduces mass error in all 72 moving-flow comparisons:
by 93.08–99.98% in Float32 and 96.78–99.83% in Float64. At T31 and a
15-second timestep, the Float64 rotating case changes from +16.0725 to
−0.0421167 kg/m².

The Float32 readouts initially failed on constant-field pressure means
that differed from identical extrema by up to 1.31e-10 Pa. A separate
reader uses the existing arithmetic tolerance for constant fields only.
All other checks remain. It recovers those 18 members and reproduces all
18 passing Float64 summaries exactly. An independent raw-record calculation
reproduces all 108 paired mass comparisons. No models were rerun.

Float64 residuals are roughly proportional to timestep. Float32 adds a
larger, less regular error.

The 36-member time-filter test completes all 216 trajectories and exact
restarts. Turning off the Robert–Williams filter changes the retained-mean
Float64 endpoint errors by less than 0.122%. The T31 15-second rotating
case changes from −0.0421167 to −0.0420756 kg/m². The nearly linear
timestep dependence remains. Keep the production filter.

The implicit-time sweep now passes all 216 trajectories and exact
restarts. Centered stepping reduces mass error by 99.14–99.84% in all
36 moving-flow Float64 cases. The worst error falls from 0.3313 to
0.00283 kg/m². Explicit stepping also improves all 36 cases.

Float32 behaves differently: 26 of 36 moving cases improve, ten worsen,
and the worst error rises from 1.11394 to 1.14887 kg/m². Keep the
production setting.

A pressure-storage probe passes 123 ownership, copy and short-evolution checks.
It widens only the pressure coefficients, or pressure coefficients and their
tendencies, while retaining the Float32 3-D atmosphere and transform buffers.
A 27-member sweep passes all 162 six-hour trajectories and exact halfway
restarts at T7/T15/T31 and three timesteps. An independent audit reproduces
every summary. Native controls match the earlier centered Float32 records.
Wider pressure coefficients improve all 36 moving cases. The largest
endpoint mass error falls from 1.148871 to 0.00308867 kg/m². Wider pressure
tendencies leave this largest error unchanged. Ten endpoints improve and three
worsen, with a maximum change of 0.000317 kg/m².

The wider type also changes spectral arithmetic, with a maximum initial
mean-pressure difference of 0.04833 Pa. A matched control keeps that wider
launch and arithmetic, then rounds the stored pressure through Float32
after each leapfrog update. All 54 trajectories and exact restarts pass.
Rounding worsens all 36 moving cases and raises the largest mass error to
1.144323 kg/m². Initial records match exactly, which supports retained
pressure precision as the cause of the improvement.

A 30-member follow-up tests 8/27 layers and 300/600/900-second timesteps.
It separates pressure precision, implicit centering and the mean-mode change,
with native controls at every setting. All 332 source assertions and five
reader test groups pass. Full-physics and coupled tests follow this dynamics test.

### BS08 runoff and heat

The twelve-process test on one A100 passes its registered compatibility rule.
Construction and initialization match exactly across 443 stored state fields. At
step two, the repeat-variation statistic is 0.22222 against a limit of 1,
with exact permutation p=395/462. Independent arithmetic reproduces the result.
This qualifies option-off compatibility for the tested configuration. The
earlier failed five-plus-five comparison remains archived.

The device-kernel test passes its physical checks and exact repeats, but its
34.0% cost increase exceeds the 25% limit. The dependent coupled test cancels.
A separate unchanged-model profile tests the suspected synchronization cost.
The coupled test retains its 3% limit.

The completed profile gives three ABBA ratios of 1.235, 1.346 and 1.382.
All 3,600 raw substep samples and repeat hashes verify. Reducing the partition
flag through its field interface takes 91.433 μs, versus 38.909 μs on the
stored interior. The one-line candidate now passes physical checks, exact
old-source state comparisons and invalid-flag tests. Its three ABBA ratios
are 1.162, 1.288 and 1.270. Two still exceed the retained 1.25 local cost limit.

A separate full-coupled diagnostic runs off/on/on/off on one A100, then
checks a fresh-process restart. It retains the routing, heat, freshwater,
correction-count and state checks. All cost ratios will be reported against
the existing 1.03 whole-coupled criterion. The earlier local cost failure
remains recorded.

The first full-coupled attempt completed its off member, then failed in a
test snapshot. The TKE boundary flux is a function, while the snapshot
assumed every flux entry owned an array. The repaired harness passes 15
regression checks. A fresh four-member comparison is running with unchanged
model code and thresholds.

### BS10 carbon inventory

Both compact CPU fresh-process restarts pass. Live-carbon and zero-exchange
runs reproduce all 1,548 state entries at the checkpoint and four continued
steps. Maximum relative inventory errors are 2.93e-8 and 2.90e-8.

The apparent GPU carbon loss came from its inventory calculation. Copying the
grid to the CPU reconstructed reference coordinates and lost changing cell
volumes. A state-preserving copy reduces maximum relative error from 2.21e-5
to 3.99e-8 in the matched 24-step GPU test. All 25 DIC and alkalinity snapshots
and air–sea uptake match the original run exactly. Both CPU controls also pass.

The zero-exchange GPU test also passes at 3.96e-8. Its DIC snapshots match
throughout, while alkalinity hashes match at 8 of 25 samples with equal extrema.

The committed repair has tests for CPU/GPU, both precisions, three grid types
and immersed boundaries. The first gate passed 552 assertions but hit 144
offset-array conversion errors in its test comparisons. The corrected test
compares offset axes and parent storage. All 1,416 geometry assertions now
pass. The compact GPU host test also passes, with relative inventory error
6.58e-9 and exact saved carbon-controller state. Both production-tripolar pilots
and their reader now pass. They retain real land, ice, radiation and ERA5/ECCO
launch states. Maximum relative carbon-inventory errors are 1.37e-9 with live
exchange and 1.38e-9 with zero exchange, below the unchanged 2e-6 limit.
The live-exchange production-grid GPU run now passes a full day. All 97
inventory samples verify, with maximum relative error 5.51e-9 against the 2e-6
limit. A nine-case sweep crosses
three prescribed CO2 levels with three alkalinity settings. It resumes the test
stopped by the earlier inventory bug and checks every step. The emissions-driven
candidate passes its compact CPU restart checks.

The full-grid zero-exchange day and all nine CO2/alkalinity cases now pass.
GPU emissions host tests pass all 27 saved samples. All 56 carbon-history
rows from four GPU restart cases pass the ledger checks, with maximum relative
inventory error 2.956e-8. Two cases reproduce the full state exactly after
restart. Two first differ in ocean fields after exact restored boundaries.
The dependent full-grid emissions batch cancelled before allocation.

The compact numeric repeats now show that the second live run and both
restores match exactly. Only the first live run differs after their common
boundary. Ocean T, S, DIC and TA and atmospheric prognostic arrays remain
exact. The mismatch can therefore occur without a restart.

An isolated test of pinned Oceananigans 0.110.15 finds a grid-update defect.
The native kernel writes sea-surface height and immediately interpolates
neighbouring heights. It can mix old and new values. Separate copy and
scaling kernels pass 472 checks, including repeatability and a CPU reference.
Both coupled restarts match all 1,550 state entries at the restored boundary
and four following steps. Carbon relative error stays below 2.957e-8.
The numeric reader verifies all 20 restart/repeat comparisons across 425
arrays each, including halos. Every difference is zero.

Four full-grid tests pass, with live and zero carbon exchange repeated twice.
All 1680 state entries match at the restored boundary and four following steps.
All 20 cross-repeat comparisons also match exactly. An independent check
reproduces all 40 restart and repeat comparisons.
A current-main candidate adds the fix as a package include. The first CPU
tests stopped at an outdated constructor call in the test. With that call
repaired, all 377 focused assertions pass. The full CPU suite also passes,
with the existing expected-broken truncation test. Both receipts verify.
The ordered update is now included on main. The carbon extension remains
under development.

Candidate `68988e3` combines the carbon extension with current main. Its three
carbon modules and five dedicated test files match the earlier emissions
candidate. All 6,710 focused CPU assertions pass, including the ordered ocean
update and changing-volume inventory checks. The full CPU suite also passes,
with the existing expected-broken initial truncation test.

Four compact GPU tests will repeat the emissions restarts on this source.
A separate CPU comparison checks that carbon-off matches main and that passive
carbon tracers leave the common model state unchanged. Four full-grid runs
follow the compact restart reader: prescribed 420 ppm, interactive CO2 with
zero emissions, and 10 PgC/year with live or zero ocean exchange. Each runs
for one model day and checks all 97 carbon-ledger boundaries. These tests
retain the earlier tolerances and exact restart comparisons.

## Evidence record

The development labbook retains source hashes, reports and failed attempts.
These job identifiers locate the evidence in that archive.

| Evidence | Jobs |
| --- | --- |
| Carbon on ordered main, CPU tests and option-off comparison | 60869035, 60874069 |
| Ordered-main emissions GPU restarts and reader | 60870291, 60870410 |
| Ordered-main full-grid emissions tests | 60874070 |
| Ordered ocean repeats and numeric comparison | 60817815, 60817817 |
| Full-grid ocean follow-up and current-main CPU tests | 60821733, 60830346, 60835380 |
| Grid-water baseline and limited-flux comparison | 60826008, 60830217 |
| Moving-cloud phase/source test and owned precision | 60831748, 60835591 |
| Native pressure continuity and zero-mode attribution | 60840424, 60841330, 60841567 |
| Six-hour pressure sweep, original reader and arithmetic recovery | 60845375, 60846103, 60846859 |
| Pressure time-filter isolation | 60850042 |
| Implicit time-centering sweep | 60861271 |
| Pressure-storage probe, sweep and audit | 60877770, 60878718, 60878998 |
| Matched pressure-rounding probe, sweep and audit | 60884700, 60884701, 60884702 |
| Pressure layer/timestep follow-up and audit | 60886934, 60886935 |
| Eight-thread convection on ordered ocean main | 60876558, 60876559 |
| July monthly convection, physical reader and ERA5 scores | 60783806, 60783807, 60784800 |
| Combined one-day coupled cases and reader | 60804346, 60804368 |
| C06 30-day persistence test | 60781768, 60781800 |
| C06 current-main CPU tests, coupled compatibility and skill samplers | 60813131, 60815837, 60809445, 60809446 |
| C06 repaired coupled comparisons, January skill and latest-main tests | 60824167, 60858832, 60859107 |
| Clean C08 one-day parameter runs and reader | 60782349, 60782351 |
| C10 joint projection, physical phases and balanced clouds | 60809034, 60811701, 60818123 |
| C10 diagnostic roundoff audit | 60818202 |
| C10 conserved-variable representation | 60819225 |
| BS08 full-coupled cost diagnostic and repaired harness | 60818157, 60835381 |
| Ocean grid-ordering witness and coupled replay | 60815972, 60817815, 60817817 |
| C10 live cloud phase, rejected spectral bound and filter witness | 60800258, 60801160, 60804150 |
| BS08 completed cost profile and reduction candidate | 60783921, 60803170 |
| BS10 GPU emissions restart and numeric repeats | 60775069, 60775070, 60800419, 60800421 |
| C02 current-model columns and ERA5 host pilots | 60681564, 60683567, 60683670 |
| C02 passive monthly capture and reader | 60691659, 60691662 |
| C02 nine-member January assessment and complete-run reader | 60692616, 60694629 |
| C02 weaker-cutoff sweep and ERA5 scores | 60772893, 60772921 |
| Clean C07/C08 coupled pilots and reader | 60776941, 60776959 |
| Research day-one restarts and reader | 60710185, 60710248 |
| BS10 live-carbon production day | 60720035_0 |
| Finite-C08 coupled pilot and readout | 60575928, 60575955 |
| Full-state coupled restart | 60652124, 60652400, 60652401 |
| Active finite-C08 coupled day | 60610464 |
| Combined C07/C08 coupled pilots | 60614913 |
| Thirty-day ERA5-start screen | 60585172 |
| Complete-run comparison with safety flags | 60616516 |
| C08 passive trace pilot and full replays | 60626405, 60626411 |
| C08 exact-column and high-precision ledger readout | 60627161 |
| C08 fixed-flux stage attribution | 60668881 |
| C08 mixed-precision component and host pilots | 60671261, 60673654, 60673724 |
| C08 capture pilot and thirty-day comparison | 60678217, 60678269, 60679415 |
| C10 exact hot-spot replays and readout | 60603150, 60605096 |
| C10 column substeps and stored-flux tests | 60606681, 60607862 |
| C10 two-day substep screen | 60607957 |
| C10 ten-day follow-up | 60618862 |
| C10 full-model pilots and readout | 60632808, 60632809 |
| C10 coupled restart tests | 60676074, 60676200 |
| C06 saved rainfall reconciliation | 60651213 |
| C06 corrected diagnostic pilots | 60654516 |
| C06 ten-day reporting repair and replacement assessment | 60656531, 60709864, 60709903 |
| C08 thirty-day parameter grid and selected coupled tests | 60710059, 60710060, 60716600, 60716619 |
| C07 and selected C08 January runs and repaired readers | 60712609, 60718700, 60746567 |
| C07 and C08 parameter scoring input checks | 60748518, 60748519 |
| Clean C08 solver and host checks | 60731314, 60738239 |
| Clean C08 coupled and restart tests | 60739236, 60752549 |
| BS08 six-plus-six GPU comparison and next tests | 60506397, 60508170, 60508197 |
| BS10 matched inventory diagnosis | 60542721, 60542722, 60649076, 60649077 |
| BS10 original grid regression and cancelled pilots | 60661629, 60662199 |
| BS10 corrected grid regression, production pilots and reader | 60668131, 60668383, 60690215 |
| BS10 full-grid restarts, day-one budgets and sensitivity sweep | 60717159, 60720035, 60756040 |

Model releases use the combined climate response, conservation, stability,
restart behaviour and cost. Each experiment keeps its stated comparison
window and limits. A useful improvement can enter a release while other
model problems remain open.
