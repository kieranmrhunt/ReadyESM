# v0.1.8 validation

Candidate; release-package qualification is pending.

## Change

Production selects `atmosphere_saturation_law: mixed_phase` and
`atmosphere_large_scale_precipitation: net_column_conservative` (D).
The implementation is based on experiment source `8fd44289b8e8e213bad235f0da35ea66b536444b`;
the saturation equations are unchanged. Other production physics and parameters
remain those of v0.1.7. Generic configuration defaults retain the old choices.
The saved-history ocean-bound checks already on main are included.

## Coupled experiment

D changes the two settings at model day 90, then continues through native
checkpoints to day 390. Production enables them from initialization; the long
trajectory is a continuation, not a fresh-start 390-day qualification.
The matched control keeps liquid saturation and upstream-observed large-scale
precipitation. It stops near day 384.948 after
negative salinity; no comparison uses its missing final days.

Over the common day 360–384 window:

| Quantity | Control | D |
| --- | ---: | ---: |
| Net TOA uptake (W/m²) | 8.75 | 3.49 |
| Surface-air temperature (K) | 288.57 | 287.21 |
| Column water vapour (kg/m²) | 32.49 | 28.41 |
| Ocean heat uptake per global area (W/m²) | 12.47 | 9.25 |
| Atmospheric water residual change (kg/m²) | -0.167 | 0.407 |

Flux and state entries are time-window means; water residual is an endpoint
difference, not a mean of cumulative values. This is one paired trajectory,
without a long-run repeat envelope. Smaller net TOA uptake and about five more
days of survival
do not establish equilibrium, eliminate regional bias or repair total water
conservation. D's diagnostic checks pass through day 390; those checks do not
impose a tight quantitative atmospheric water-closure bound.

D's maximum saved ocean temperature is 39.76°C, final minimum salinity 2.56 psu
and final maximum area-equivalent ice thickness 22.19 m. These extrema can move
between cells. Maximum absolute land-water residual is 2.97e-9 mm; it does not
measure atmospheric or whole-system closure.

An independent netCDF4/NumPy recomputation verifies the table's units, clocks,
exact window endpoints and rounded values. D diagnostics SHA256:
`acd18f951d340a061c35e21e13427b9cd44c7719dd9ee7ff91fef87eaa39f20e`.
The control ends early and is not labelled a completed 390-day run.

## Trade-offs

- From day 90 to 384, the atmospheric water residual grows by 5.431 kg/m²
  in D versus 2.172 in the control: 0.01847 versus 0.00739 kg/m²/day. D's
  late source persists at 0.01694 kg/m²/day over days 360–384. Surface rain
  matching an individual scheme's water removal does not close this budget.
  Separate standalone tests find this remainder sensitive to consistency
  between initial surface pressure and model orography; coupled verification
  is pending, and v0.1.8 does not change the launch.
- D is 1.36 K cooler and 4.07 kg/m² drier in the matched late mean. The
  endpoint gaps increase from 0.42 K/1.39 kg/m² at day 120 to 1.59/4.21 at
  day 384, but do not grow monotonically. No equilibrium claim is supported.
- Ocean heat uptake exceeds net TOA uptake by 5.76 W/m² in D and 3.72 in
  the control. This difference alone is not an energy-conservation residual:
  atmosphere, land and ice storage and coupling transfers must also be
  included. Whole-system energy closure remains unqualified.
- The latitude/season-independent analytic ozone profile is unchanged.
  A preliminary fixed-column comparison gives about 460 DU versus 245–378
  in sampled ERA5 columns; this is not a global observed range.
- GPU trajectories are not guaranteed bitwise reproducible, and D has no
  independent long repeat. Exact checkpoint restoration is a narrower test.

The saturation law also affects diagnostic cloud detection. The net-radiation
response is not attributed solely to one process, and the inherited parcel
and condensation derivative approximations remain unchanged.

## Remaining checks and limitations

- Full CPU suite passes (job 59313164, 19m55), including 419 mixed-phase assertions,
  legacy configuration defaults and the production saturation/precipitation
  combination. Runtime hashes match before and after the tests.
- Exercise a fresh production-resolution GPU launch and an independent-process
  restore, checking the restored state exactly and subsequent diagnostics.
- README maps and timelines now show D through day 390. Figure job 59314808
  passes all 11 sample-contract assertions and both images have been visually
  checked. The day 90 switch is labelled; missing initial samples stay missing.
  [Figure provenance](readiesm-v0.1.8-snapshot.json) records source/image hashes.

The experimental moist-adiabat derivative repair, ozone replacement, prognostic
clouds, new convection, water fixers and salinity advection changes are excluded.
Ocean hot cells, positive ocean heat uptake and coastal ice pile-up remain.
See [NOTES.md](../NOTES.md) for the retained limitations.
