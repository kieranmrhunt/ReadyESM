# Experiment notes

Updated 19 October 2026. Production is v0.1.10. Research candidates stay in
separate source snapshots until their results support a model update.
The [README](../README.md#experimental-roadmap) tracks the ten blue-sky ideas.

## Convection portfolio

| Track | Status | Result or next test |
| --- | --- | --- |
| C01: mixed-layer parcel launch | Closed | The stronger launch reduces tropical deep-convection area too much. |
| C02: depth-limited shallow adjustment | Reserve | Retain the 0.7-depth setting for an observational comparison. The 0.8 setting is closed. |
| C03: surface-forced shallow closure | Closed | Moisture improves, but the intended storm-track response and spatial mechanism fail their tests. |
| C04: elevated and slantwise convection | Closed at T31 | Eligible columns are rare at this resolution. Reopen for a higher-resolution or event-based test. |
| C05: non-equilibrium deep closure | Closed | The tested memory produces little change in convective timing. |
| C06: rain evaporation and cold pools | Active | Saved data confirm a rainfall-diagnostic error. The correction passes 13 component tests. Three coupled pilots precede six ten-day climate runs. Source coverage remains under review. |
| C07: convective momentum transport | Active | Four coupled pilots pass. Both momentum-plus-heat settings pass all three 30-day repeats, but retain large cold and dry differences from native convection. |
| C08: deterministic multi-plume convection | Active | A finite-response closure replaces invalid trial columns in the original search. Coupled and climate tests are in progress. |
| C09: stochastic plume population | Closed | Eight seeds show no resolved change in the tested mean heating. |
| C10: shallow eddy-diffusivity/mass-flux transport | Active | All nine two-day tests pass. Six ten-day tests and six full-model pilots are running. |
| C11: saturated-parcel trigger census | Closed | The affected columns contribute too little rain to justify this route. |

## Latest results

### C06

All five one-day coupled pilots completed. Native and disabled controls pass.
The three active cases pass their internal conservation checks but fail the
rainfall diagnostic and the source-coverage screen.

The rainfall diagnostic omitted C06's separate humidity tendency. Adding it
reduces the worst saved-column discrepancy from 4–8e-5 to about 3e-9 kg/m²/s.
Every final column then meets the unchanged 3e-6 tolerance. The corrected
collector exports base, C06 and total tendencies. Thirteen component tests
pass. Three GPU pilots test the collector, file output and full validator.

The source review finds 14 cases with heating or drying below the diagnosed
cloud interval. They carry at most 0.01485% of sampled rain. C06 leaves all
unmapped rain unchanged. The larger dry-parcel class is tracked separately.
The earlier 1.25% coverage limit came from a project screening choice.

The failed screen remains archived. A separate ten-day assessment retains
the physical checks, reports both coverage flags, and checks that unmapped
rain passes through without evaporation. It compares native, middle and
strongest settings with two repeats each, after the diagnostic pilots pass.

### C08 and C07

The revised C08 closure passes eight coupled steps with the ocean, sea ice,
land and radiation components. All snow and conservation checks pass. The
disabled version reproduces the original C08 endpoint exactly. The active
case completes 101,396 closure evaluations with zero rejected probes.

Packed spectral arrays required a repair to the restart checker. Its 22
regression checks pass. Four full-state restart tests and an active one-day
coupled run are running. Four full-ESM pilots pass for finite C08,
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
exchange and its host representation. Their absolute budget effect is under
investigation. Passive column traces require exact reproduction of each
reference report and daily map. The acceptance limits remain unchanged.

All six combined C07/C08 cases pass their physical checks. Days 10–30 mean
column water is 16.535 and 16.440 kg/m² for the centre and faster settings,
versus 22.996 for native convection. Mean atmospheric temperature is 249.240
and 249.362 K, versus 251.256 K. These large differences need observation-based
assessment. The completed comparison retains every earlier failure flag.

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
Six eight-step full-model pilots are also running. These retain the actual
ocean and land surface fluxes and include native and disabled controls.

### BS10 carbon inventory

The apparent GPU carbon loss came from its inventory calculation. Copying the
grid to the CPU reconstructed reference coordinates and lost changing cell
volumes. A state-preserving copy reduces maximum relative error from 2.21e-5
to 3.99e-8 in the matched 24-step GPU test. All 25 DIC and alkalinity snapshots
and air–sea uptake match the original run exactly. Both CPU controls also pass.

The committed repair now has tests for CPU/GPU, both precisions, three grid
types and immersed boundaries. Two short pilots on the full production
tripolar mesh depend on these tests. They retain the real land, ice, radiation
and ERA5/ECCO launch. Longer response and restart tests follow.

## Evidence record

The development labbook retains source hashes, reports and failed attempts.
These job identifiers locate the evidence in that archive.

| Evidence | Jobs |
| --- | --- |
| Finite-C08 coupled pilot and readout | 60575928, 60575955 |
| Full-state coupled restart | 60652124, 60652400, 60652401 |
| Active finite-C08 coupled day | 60610464 |
| Combined C07/C08 coupled pilots | 60614913 |
| Thirty-day ERA5-start screen | 60585172 |
| Complete-run comparison with safety flags | 60616516 |
| C08 passive trace pilot and full replays | 60626405, 60626411 |
| C08 exact-column and high-precision ledger readout | 60627161 |
| C10 exact hot-spot replays and readout | 60603150, 60605096 |
| C10 column substeps and stored-flux tests | 60606681, 60607862 |
| C10 two-day substep screen | 60607957 |
| C10 ten-day follow-up | 60618862 |
| C10 full-model pilots | 60632808 |
| C06 saved rainfall reconciliation | 60651213 |
| C06 corrected diagnostic pilots | 60654516 |
| C06 separate ten-day coupled assessment | 60656531, 60656533 |
| BS10 matched inventory diagnosis | 60542721, 60542722, 60649076, 60649077 |
| BS10 committed grid regression and production pilots | 60661629, 60662199 |

Model releases use the combined climate response, conservation, stability,
restart behaviour and cost. Each experiment keeps its stated comparison
window and limits. A useful improvement can enter a release while other
model problems remain open.
