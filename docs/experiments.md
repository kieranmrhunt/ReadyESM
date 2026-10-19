# Experiment notes

Updated 19 October 2026. Production is v0.1.10. Research candidates stay in
separate source snapshots until their results support a model update.
The [README](../README.md#experimental-roadmap) tracks the ten blue-sky ideas.
Main now includes the accepted [snow-accounting update](validation-snow-accounting.md).

## Convection portfolio

| Track | Status | Result or next test |
| --- | --- | --- |
| C01: mixed-layer parcel launch | Closed | The stronger launch reduces tropical deep-convection area too much. |
| C02: depth-limited shallow adjustment | Active | All nine January runs and the complete-run reader pass. The ERA5 scorer is tested and awaits its queued reference files. The 0.8 setting is closed. |
| C03: surface-forced shallow closure | Closed | Moisture improves, but the intended storm-track response and spatial mechanism fail their tests. |
| C04: elevated and slantwise convection | Closed at T31 | Eligible columns are rare at this resolution. Reopen for a higher-resolution or event-based test. |
| C05: non-equilibrium deep closure | Closed | The tested memory produces little change in convective timing. |
| C06: rain evaporation and cold pools | Active | The rainfall-diagnostic correction passes 13 component tests and all three coupled pilots. Ten-day climate runs have started. The active dry-parcel screening flag remains failed. |
| C07: convective momentum transport | Active | Eight repaired mixed-precision and control pilots pass. Fifteen January comparisons are submitted for wind scoring. |
| C08: deterministic multi-plume convection | Active | All nine mixed-precision 30-day runs pass. January scoring and a three-seed parameter grid are in progress. |
| C09: stochastic plume population | Closed | Eight seeds show no resolved change in the tested mean heating. |
| C10: shallow eddy-diffusivity/mass-flux transport | Active | All nine two-day tests, six full-model pilots and six coupled restart tests pass. Ten-day climate tests are running. |
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
precipitation by 0.054–0.069 mm/day. These changes need observational scoring.
The scorer passes eight metric tests and independent checks of all nine files.
ERA5's server confirms the reference request is queued. Water, precipitation
and temperature will use common masks and area weights. These
atmospheric tests use prescribed seasonal SST and simple radiation. Coupled
assessment follows if the combined observed response supports it.

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
and the ten-day runs have started.

### C08 and C07

The revised C08 closure passes eight coupled steps with the ocean, sea ice,
land and radiation components. All snow and conservation checks pass. The
disabled version reproduces the original C08 endpoint exactly. The active
case completes 101,396 closure evaluations with zero rejected probes.

Packed spectral arrays required a repair to the restart checker. Its 22
regression checks pass. Four full-state restart tests pass, with exact
checkpoint and continued states. A day-one restart extension and an active
one-day coupled run are running. Four full-ESM pilots pass for finite C08,
zero-strength C07, momentum transport, and momentum transport with heat return.

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
January comparisons and a 27-setting, three-seed parameter grid are running.

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
Fifteen January runs compare native, C08 and C07 across three seeds. Passive
wind output at 200, 500 and 850 hPa will support comparison with ERA5.

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

Six ten-day follow-ups are running to check
persistence, global moisture and rain, and sensitivity to the internal step.
Each longer run must reproduce its two-day history and daily maps.
All six eight-step full-model pilots pass. These retain the actual ocean and
land surface fluxes. Disabled prognostic fields, surface memories and scientific
arrays match native exactly. All six fresh-process restart tests pass.
The complete coupled state and all three forcing memories match at step four
and each continued step through step eight. Native has 1,681 compared entries,
and each experimental arm has 1,771. Active cases have zero physical failures.

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

### BS10 carbon inventory

Both compact CPU fresh-process restarts pass. Live-carbon and zero-exchange
runs reproduce all 1,548 state entries at the checkpoint and four continued
steps. Maximum relative inventory errors are 2.93e-8 and 2.90e-8. The first
production-grid GPU member is running, with the second waiting for priority.

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
6.58e-9 and exact saved carbon-controller state. Two production-tripolar pilots
await GPU priority. They retain the real land, ice, radiation and ERA5/ECCO launch.
A queued reader checks all inventory samples and preserves failed outcomes.
Longer response and restart tests follow.

## Evidence record

The development labbook retains source hashes, reports and failed attempts.
These job identifiers locate the evidence in that archive.

| Evidence | Jobs |
| --- | --- |
| C02 current-model columns and ERA5 host pilots | 60681564, 60683567, 60683670 |
| C02 passive monthly capture and reader | 60691659, 60691662 |
| C02 nine-member January assessment and complete-run reader | 60692616, 60694629 |
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
| C06 separate ten-day coupled assessment | 60656531, 60656533 |
| BS08 six-plus-six GPU comparison and next tests | 60506397, 60508170, 60508197 |
| BS10 matched inventory diagnosis | 60542721, 60542722, 60649076, 60649077 |
| BS10 original grid regression and cancelled pilots | 60661629, 60662199 |
| BS10 corrected grid regression, production pilots and reader | 60668131, 60668383, 60690215 |

Model releases use the combined climate response, conservation, stability,
restart behaviour and cost. Each experiment keeps its stated comparison
window and limits. A useful improvement can enter a release while other
model problems remain open.
