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
| C06: rain evaporation and cold pools | Active | All 33 ten-day screen members pass safety checks. Five full-ESM GPU pilots are queued. |
| C07: convective momentum transport | Active | A translating-frame calculation removes a small spurious wind force. A conservative heat-return option now runs in full-ESM pilots. |
| C08: deterministic multi-plume convection | Active | A finite-response closure replaces invalid trial columns in the original search. Coupled and climate tests are in progress. |
| C09: stochastic plume population | Closed | Eight seeds show no resolved change in the tested mean heating. |
| C10: shallow eddy-diffusivity/mass-flux transport | Active | Shorter internal transport steps remove the severe hot spots in the two-day tests. Longer tests are queued. |
| C11: saturated-parcel trigger census | Closed | The affected columns contribute too little rain to justify this route. |

## Latest results

### C08 and C07

The revised C08 closure passes eight coupled steps with the ocean, sea ice,
land and radiation components. All snow and conservation checks pass. The
disabled version reproduces the original C08 endpoint exactly. The active
case completes 101,396 closure evaluations with zero rejected probes.

Four full-state restart tests and an active one-day coupled run are running.
Four further full-ESM pilots compare finite C08, zero-strength C07, momentum
transport alone, and momentum transport with heat return.

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

The six combined C07/C08 cases are running. Wind profiles will be compared
with ERA5 January climatology.

### C10

The stronger shallow-transport settings improve some global means but cause
local temperatures near 358 K in two days. Saved-column tests locate direct
convection increments as large as 77.6 K in one effective step.

Recomputing transport over shorter internal steps reduces those increments.
The host receives the divergence of the time-averaged water and energy fluxes.
Both eight-substep host cases pass the original state and conservation checks.

| Updraft area fraction | One-step maximum | Eight-substep maximum | Change in second-day global mean temperature |
| --- | ---: | ---: | ---: |
| 0.04 | 358.50 K | 324.16 K | −0.0018 K |
| 0.08 | 358.43 K | 324.28 K | −0.0069 K |

Both cases use 1000 m entrainment length and buoyancy-defined plume height.
The maxima cover all grid points, layers and completed steps. Their locations
and times can differ. Global means use 24 hourly samples from the second day.
The native control reaches 327.43 K. These tests measure numerical behaviour
and model sensitivity. Observational climate skill remains to assess.

Both 16-substep cases also pass. Their maxima are 324.18 K at area fraction
0.04 and 324.33 K at 0.08, close to the eight-substep results.
The 32-substep tests are running. Six ten-day follow-ups will check
persistence, global moisture and rain, and sensitivity to the internal step.
Each longer run must reproduce its two-day history and daily maps.

## Evidence record

The development labbook retains source hashes, reports and failed attempts.
These job identifiers locate the evidence in that archive.

| Evidence | Jobs |
| --- | --- |
| Finite-C08 coupled pilot and readout | 60575928, 60575955 |
| Full-state coupled restart | 60595966 |
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

Model releases use the combined climate response, conservation, stability,
restart behaviour and cost. Each experiment keeps its stated comparison
window and limits. A useful improvement can enter a release while other
model problems remain open.
