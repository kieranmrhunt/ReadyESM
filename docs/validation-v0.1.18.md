# v0.1.18 validation

This record describes the exact-source gates completed by 30 October 2026.
The failed CATKE
wind-size gate, inherited limitations and neutral skill-rule REJECT are
disclosed. The final GPU restart passed and the figures are installed.
The release tag and remote assets are checked separately by the post-tag
audit. This source-scoped record does not establish their status.

## Candidate and changes

The tested source is commit
`48c232cc6d313370c9d5c8fd2a52745af4e6c2ac`, tree
`eb1ea01c088cbe39e5e22f29e2d377d3e27f1fe4`. The candidate includes
surface-energy input to native CATKE mixing, optional scored outputs,
restart-identity checks, a per-step prognostic-state guard, and a tracked
CUDA 13.2 runtime preference. The standard production configuration uses
RRTMGP radiation. A scored-output run cannot currently resume from a saved
state. The guard stops a run if a prognostic value is non-finite or outside
its fixed physical limit; it changes run behaviour and is not just a report.

The first attempt at the six paired month runs on earlier commit `381c1492`
stopped before its first step because the new guard read SpeedyWeather's
spectral and grid arrays on the GPU one element at a time, which CUDA does
not allow. The wrapper test reproduced that fault on the failed commit.
A first repair, `ad90eeff`, passed 24 of 25 wrapper cases but still failed
on the error path of an Oceananigans `Field`. Commit `48c232cc` repaired
that case. The six replacement paired runs completed on exact source
`48c232cc`; no result from the failed array is counted.

## Completed checks

| Check | Evidence | Scope |
| --- | --- | --- |
| GPU prognostic guard wrapper | `v018_guard_gpu_test_v1/results/b3_48c232cc_run1`, 27/27 | Real A100 wrapper coverage on the exact candidate, including Oceananigans Field and array wrappers. The same test reproduced the scalar-indexing fault on `381c1492`. This is not a production integration run. |
| Six paired configuration preflights | Slurm array `63369713_[0-5]`, all `COMPLETED 0:0`, and `v018_paired_preflight_review_v1/review.json` SHA-256 `cd6f755f61d0efb37528ca4bf327f70093b731468e724e8779f1ecd0cc626648` | JAN/JUL members 1–3, matched control inputs from a source equivalent to v0.1.17, and resolved configurations. |
| Focused CPU tests | `63365948/focused_tests.log`, `V018_COMBINED_FOCUSED_CPU_PASS` | Restart identity 18/18, prognostic guard 22/22, native CATKE surface flux 133/133, scored-output configuration 4/4, accumulators 11/11 and coupled run 18/18. These tests preceded the full one-day qualification. |
| One-day scored-output CPU qualification | Job `63365948`, `COMPLETED 0:0`, verdict `SCORED_OUTPUT_WRITER_CPU_V5_PASS`, result-manifest SHA-256 `f520ce5b8593d3b35ee09ce8e25d0cfcf5a2fae293fde58e2419be18639f4eea`, frozen input seal `b3f28ae174e4c321f5b7459c8abd5e8c5f4bfdabb8b42a04cde6656e9dbf3502` | Exact source `48c232cc`/tree `eb1ea01c`, January one-day coupled run. All four writer/observer files are byte-identical. The final saved model states match exactly, with zero excluded identity fields. The result manifest and full job inventory were independently verified. |
| Restart-identity CPU preflight | Job `63368988`, `COMPLETED 0:0`, verdict `V018_RESTART_PREFLIGHT_PASS`, 16-file result-manifest SHA-256 `7797df0942db49326e22837d24f1089796703c625fe6c2b95820695193b0a913` | Bound to the audited CPU qualification receipt `f520ce5b8593d3b35ee09ce8e25d0cfcf5a2fae293fde58e2419be18639f4eea`. The same-identity checkpoint round trip passes. Commit, manifest and config mutations each produce the intended one-field refusal. Verifier tests 16/16 pass. This is a CPU preflight, not the A100 restart gate. |
| Pre-tag metadata audit | Manifest-aware checker SHA-256 `d565d1a2b7ad1130e537a8bc32cbd915af3d89c6c3636347398be4081d08cced`, clean-b3 check log SHA-256 `53d5431e6cbdcb8018fee50761fab6e1b9965abc38a07137273fc75765ed63a6` | Project and local Manifest ReadyESM versions both read 0.1.18. Manifest path/UUID, changelog, citation date, vendor records and tracked path checks passed on the earlier draft. CITATION lists Kieran Hunt and the ReadyESM contributors, with no DOI. This pre-tag audit did not check a release tag; the tagged tree is checked separately. |
| Release-doc review wrapper fixtures | Isolated wrapper SHA-256 `91630c5309f08cc85112a6d8bf3e318cec3afed78c1e429849e9cf5372ad9107`, 24/24 tests, frozen-input manifest SHA-256 `41e0ce4b4fab5f4f81d3819f0f5ab0943353a0b7c001785809177548432db90f` | Checker fixtures cover both lightweight and annotated tags. One fixture deliberately uses the earlier heading `## 0.1.18 — 2026-10-29`; the release heading is 30 October. Matching metadata passes, a stale Manifest fails, and an annotated tag that no longer points to HEAD is refused. These fixtures do not establish a release PASS; a separate substantive review and actual tagged-tree check are required. |
| Production-shaped GPU step | `63365988`, `COMPLETED 0:0` | One A100 step on exact b3 source with scalar indexing disabled and the prognostic guard active. This establishes startup and one-step integration only. |
| Six paired GPU months | `63369995_[0-5]`, all `COMPLETED 0:0` | Three January and three July 1993 pairs, 31 days each, against matched v0.1.17-equivalent D4 controls. The frozen scorer verified all six candidate outputs and matched-source inputs. |
| Frozen skill and neutral readout | `63449460`, `COMPLETED 0:0`, result-manifest SHA-256 `3a08fd5990d445be02f417215f35470fcf737bfea4e1f4951cc941769afbd555` | Independently verified with the v3 readout-input seal `31b8f22b532e59ac815831c0d022623dcb00d48195e37eea0b376a63c57fad1c`. Frozen v3 REJECT, v4 UNDECIDED and separate neutral-rule REJECT are preserved below. |
| Native CATKE GPU | `63374110`, `COMPLETED 0:0`, result-manifest SHA-256 `d5119315152f41c01d09549f5a4534c42ca04aaf58f4b6d1a09b0911bf154539` | Exact-source A100 native surface-flux test, 133/133; verdict `V018_CATKE_NATIVE_GPU_PASS`. This does not change the failed idealized six-hour wind-size gate below. |
| Exact-source A100 restart gate | `63495868`, `COMPLETED 0:0`, verdict `V018_RESTART_GPU_PASS`, 55-entry result-manifest SHA-256 `2ef4cac0bf6a3421d101b0361eb8606bb46ba2c92610aca224593fefd1a8cf47` | Exact source `48c232cc`/tree `eb1ea01c`. Six exact saved-state comparisons and six nonvacuity difference controls passed. The frozen input seal `25b71ad7413cabd7161dd51fe09881e643be76acb7661710e7834a2347d329b4`, result manifest and frozen verifier were independently checked. The prior `63368989` harness failure is superseded, as explained below. |
| Canonical-start version history | `63378103_[0-1]` and CPU staging job `63496364`, all `COMPLETED 0:0` | January and July 31-day starts with matching 15-field initial hashes. The original publication figure has result-manifest SHA-256 `04aaf8ffceb67d558366edcf807962e53e65e66901a9ea4c96f608d8ece05e75`; its reviewed plain-label v3 successor has manifest `a9ede187a87ca900dd70c1229e9408a6381106e70605d15bb6ed0d16479a3c1e` and uses the same scored inputs. These single starts are descriptive and are not the six-pair release estimate. |

The full one-day CPU qualification passed. Its sealed report records identical
writer and observer outputs for all four files and exact equality of the
writer-on and writer-off final states. The independent audit checked the
complete result manifest, input seal, report, verdict, source tree and result
inventory.

## Skill results and release disposition

The frozen original seven-field v3 scorer returned **REJECT**: January was
NEUTRAL, but July was WORSE. July mean seven-field normalised error (ONR7)
rose by `0.01154899` against its noise bound of `0.01026302`, exceeding it by
`0.00128596`. That bound is the
standard deviation of the three matched control members under the frozen v3
rule, with a floor of `0.01`.
The nine-field v4 scorer returned **UNDECIDED** in both seasons and overall.
Its January mean nine-field normalised error (ONR9) changed by `-0.00023490`;
July rose by `0.00801674`, within the v4 `0.01000000` bound. All frozen v4
global radiation, net and
physical-breakage guards passed. The separate neutral adjudication applied
its written rule and returned **REJECT**, because July v3 ONR7 was clearly
worse. This result remains REJECT in the release record; it is not converted
to PASS.

The six members are observation-initialised January and July months. In a
separate six-field baseline test, climatology and persistence beat all 59
CMIP6 models. CMIP6 position here measures agreement with the observed month;
it does not establish free-running climate skill.

Using each sealed `report.toml`'s elapsed time from day 1 to day 31, b3 took
4.23, 4.82, 10.75 and 15.64 percent longer per model day over the four pairs
whose control and candidate ran on the same GPU host; the mean is 8.86
percent. JAN member 1 and JUL member 3 ran on different hosts and are
excluded. A separate per-step guard timing test took under 1 percent; the
remaining cost is not attributed to a single change.

## Restart gate and remaining release work

The replacement exact-source A100 restart job `63495868` passed. Its CPU
preflight `63449982` passed. The earlier GPU job `63368989` matched saved
states exactly at steps 4 and 8 and passed the refusal cases, then failed at
`control_boundary` because the test script registered a callback after the
refusal branch. It is a superseded test-harness failure, not evidence of a
model restart mismatch. The replacement's strict saved-state comparisons,
refusal and non-mutation controls passed with a sealed result, closing O-088
for the b3 default path. The prepared scorecards are installed in this tree.
Documentation review, tag and remote audit are separate release checks.

The native-stage six-hour wind-column gate requiring at least 0.1 m of
mixed-layer-depth change **failed**: the repaired minus zero-source depth was
+0.0624 m at six hours. The fixed arm had 1.244 times the top-cell turbulent
kinetic energy (TKE) and 1.043 times the top-4-m tracer diffusivity at that
point. The native-stage ledger and ordinary-step trajectory checks passed.
In that ledger, the added surface source was offset mostly by extra dissipation
(75–81%) and partly by less shear production (19–23%), so little of it remained
as TKE. The depth is interpolated across 2 m cells using a 0.03 kg/m³ density
threshold, so the 0.0624 m difference is well below the grid spacing.
These are idealized CPU-column mechanism results, not production-grid effect
sizes. The failed wind-size result remains a limitation for the CATKE release
decision and must not be relabeled as a pass. A separate cooling test initially
used an incorrect positive-potential-energy criterion. That criterion was
withdrawn, and the refined 0.5 m cooling-centre and conservation sign checks
passed 8/8. Their effect size is sensitive to the grid. The CATKE change ships
with the failed wind-size gate disclosed. This does not make that gate pass.

## Inherited dependency limitations

The selected production path retains two dependency limitations also present
in the v0.1.17-equivalent control. Both are disclosed in NOTES.md.

- **Snowmelt latent heat:** `config/production.yml` selects
  `net_column_conservative`. Its `NetColumnImplicitCondensation` wrapper calls
  SpeedyWeather 0.21.1 `ImplicitCondensation` and adjusts diagnosed rain and
  snow rates after that call without changing its temperature tendency. The
  upstream omission of cooling when snow melts aloft therefore remains in the
  production path. The full-model effect has not been measured. The isolated
  repair at `b6a9dd4` is in the unselected opt-in
  `process_conservative_v1` scheme, whose 28 October screen was UNDECIDED.
  The b3 and control manifests both pin
  SpeedyWeather 0.21.1 at tree
  `3c1a4dbb5558c16639fe6f74af9b10d89ddfdce0`. The selected wrapper is
  implemented in `src/conservative_atmosphere_humidity.jl` and selected in
  `src/baseline.jl`.
- **Near-polar wind rotation:** the coupler turns the atmosphere's east and
  north winds into the ocean grid's directions through NumericalEarth
  `_rotate_winds!`, which uses Oceananigans `rotation_angle`. That angle is
  estimated from corner latitude differences and is less accurate near the
  geographic pole. The candidate and matched control pin NumericalEarth
  0.6.0 at tree `7f07b8657ac5d29bd6b71b5191c1529f505fd6ea` and
  Oceananigans 0.110.15 at tree `0b9f838634eb636b493cf473b35be62d3184179a`.
  In an offline 1 m/s cross-polar test flow on the production 360 x 180 x 60
  tripolar grid, the largest east/north velocity error on active ocean cells
  is below 0.0025 m/s south of 85 N, 0.014 m/s at 85-88 N and 0.054 m/s
  at 88-89 N. Its result record is
  `climatebench_cmor_v1/rotate_selftest_record_63355687.toml`, SHA-256
  `756bf4f546f69e33b24337e0f768fd44a4ac295bcccae3ecd42af53601c84af0`.
  The first production-grid probe `63375518` stopped at config load because
  b3 rejects the gas-file and ozone-file keys. Its rotation job `63377880`
  was cancelled. A later offline CMIP writer chain on a clone of `48c232cc`,
  with those forcing files removed, passed one-step output, rotation and
  conversion checks. It is not release-line production qualification.
  The effect of the near-pole direction error on coupled climate remains
  unmeasured.

Both limitations remain open.

## Other inherited production limitations

The exact b3 candidate still has the following selected-path issues from
`ORACLE_ISSUES.md`. A plain-language summary appears in `NOTES.md`. None of
the separate candidate fixes for these issues is part of v0.1.18.

| Issue | Production effect | Evidence and release boundary |
| --- | --- | --- |
| O-005 | Sea-ice vapour enters the air without an ice or snow mass debit; the heat-of-fusion term is also missing. | In the separate v0.1.15 long run, 51 complete monthly files through March 1997 were verified against the committed output-prefix hashes. Integrating the sampled, phase-weighted flux over the four full years 1993–1996 gives 11.127 mm global water equivalent, or 2.782 mm/year. See `v015_ice_vapor_budget_v1/receipt.json` (SHA-256 `118e0ea4d13b19e3deea2945f7d4fe492d60daa53f69736c5105e8ad4f7af194`). This measures the writer's latest post-step cached flux, not a prognostic mass-closure residual or the response to a fix. The opt-in correction is unmerged. |
| O-007 | The selected Terrarium van Genuchten–Mualem conductivity has the wrong exponent. | A native reference comparison predicts conductivity 8.6 times too high at half saturation and about 520 times too high at 10% saturation. The coupled soil-water and climate effect is unmeasured; a corrected candidate is unmerged. |
| O-009 | Resting sea ice feels half its full internal pressure instead of zero, which can drive motion where ice strength varies. | A native GPU constitutive repair passed its focused checks, but coupled dynamics and restart qualification remain open. The fix is not in b3. |
| O-010 | The layered ice thermodynamics does not conserve energy when ice area grows and snow is then exhausted. | The affected case is known; its coupled size is unmeasured and the correction is unmerged. |
| O-002 | The negative-humidity repair also removes very small positive humidity below about 1.2e-7 kg/kg. | CPU/GPU focused corrections exist, but coupled and restart qualification remain open. The fix is not in b3. |
| O-086 | ECCO initial ocean velocities are misregistered to model faces and the tripolar fold. | A correction has focused evidence, but its coupled release effect is unmeasured and the fix is not in b3. |

O-023 is the CATKE change being qualified for this release and is covered by
the gates above. O-087 belongs to an unmerged combined integration, so it is
not described as a b3 limitation. These inherited issues are disclosed in
NOTES.md.

Two further open diagnostic issues need public disclosure. In O-012, some
atmospheric diagnostic states carry a timestamp about one 900-second model
step away from the state they describe, and the water residual is staggered.
This is small for monthly means but matters for six-hourly output and water
budget timing. In O-015, the standard output labels tripolar ocean and
sea-ice surface velocities east/north even though they are grid-axis values
at tracer centres; the direction is wrong where the grid turns near the
Arctic. Both have isolated repair candidates, but full writer, restart and
GPU qualification is pending. Neither correction is claimed as shipped.
Both limitations are described in `NOTES.md`.
The offline CMIP converter is a separate workflow and is not claimed here as
part of the release.

O-088 was closed in the issue tracker on 30 October after the A100 gate. A
prior candidate did not bind the precipitation scheme in its checkpoint
identity. The exact b3 candidate hashes every `ExperimentConfig` field except
named run controls,
including `atmosphere_large_scale_precipitation`. A targeted Julia check on
b3 produced distinct config hashes for `:upstream_observed` and
`:net_column_conservative`. Checkpoints without the new identity are refused
by default, while the explicit mismatch override remains available. The
production precipitation setting is `net_column_conservative` in both the
`a2f2ac4` control and b3 (`config/production.yml:17`). The exact-source A100
restart gate `63495868` passed, including refusal before state copy and exact
uninterrupted/resumed states. This identity evidence does not settle
legacy-checkpoint compatibility under the override.

Four inherited model biases are also listed in `NOTES.md`.
They are observations from earlier versions, not measured v0.1.18 changes:

- The v0.1.17 eddy screen used one 30-day member per season. Its 250-hPa
  eddy-energy ratios to resolution-matched ERA5 are 0.65–0.79 across the two
  hemispheres and seasons; its predeclared 0.10 improvement target was missed.
  Sampling uncertainty is not measured. Source:
  `d4_publication_v017_v1/repo/docs/validation-v0.1.17-eddy-screen.json`,
  SHA-256 `58687d4a064df83a3af2acaefad5bbf4204c0698ef3d84ae9ffdeb4d25fd414f`.
- The standard land albedo is 0.4. An offline comparison over 509
  land-bearing T31 cells south of 60°S gives ERA5 forecast-albedo means
  0.8455 in January and 0.8484 in July. This mask is not a map of permanent
  ice, and its coupled effect is unmeasured. Source:
  `antarctic_fal_attribution_v1/summary.json`, SHA-256
  `7dc8101a5fa61c46e7843cb35dfba5aba65de3e159a7e33cf087530953216964`.
- The sea-ice albedo setting is 0.60 with no snow response. In the v0.1.15
  long run, May blended surface albedo north of 80°N was 0.57–0.58 versus
  0.77–0.79 in ERA5. This is a grid-box comparison, not a pure ice-tile
  albedo measurement. Source: `arctic_sea_ice_check_v1/arctic_sea_ice_check.json`,
  SHA-256 `4c2e904468ff443b2742565be8a707d69b32742390c1779e9864ce337588b297`.
- The 1993–1996 v0.1.15 run produced about half of ERA5 global snowfall
  while total precipitation was close. The selected phase rule makes
  large-scale snow only below minus 10°C and always treats convective
  precipitation as rain. No experiment has isolated that rule's share of
  the snowfall gap. Source: `arctic_sea_ice_check_v1/snowfall_vs_era5.json`,
  SHA-256 `9358801b2675f925cb05b94d09d3f2ecfbf6e00f042d37f28217cc5a41e74ea2`.

The long forced climate continuation from the v0.1.15 line is separate
research evidence. Its final segment `63270822` and scorer `63270887`
completed with exit code 0, and the scorer's sealed result verifies with
`CLIMATE_RUN_SCORE_PASS`. For the 1994–1997 period, one trajectory has annual
ONR9 `2.9629` and beats 0 of 54 complete CMIP6 models; January and July
ONR9 are `2.4411` and `2.7225`. The atmosphere is cold and dry. The five-year
result's annual net top-of-atmosphere gain rose from `-0.29` W/m² in 1993 to
`+4.18` W/m² in 1997. This pattern is not a tested causal chain. The run
cannot replace paired v0.1.18 skill and safety checks, and no
long-term climate-skill claim follows from the 31-day tests.

## Figures and release record

The public README maps and timelines have been rendered from sealed b3
January member 1 (`v018_readme_maps_b3_v1`, result-manifest SHA-256
`92e0ddd80b1ee3dbc297a70f3595bd5847f84a7ffea43576bd79d98097542dde`).
The six-pair paired skill scorecard is prepared in
`v018_postgate_evidence_codex_v1/scorecard_prepared_v4`, with
`staged_outputs.sha256` SHA-256
`aded5f53baf212bbefe9ffcde8b63eaf4f32a4f12e685915595c0f794a73f7e4`
and PNG SHA-256
`1678b4a442de6a52f17dd3e6b0740456ec25c6869c1d3b4a54ee215349e71a91`.
The reviewed development scorecard is prepared in
`development_scorecard_v018_b3_v6_historybound/release_stage_63495868_v3`, with
`result.sha256` SHA-256
`3bdf1565bb72cee2e3d5b1fd5d138b5434c83686ab19a479c43bc9ca0265e681`
and PNG SHA-256
`bb50a0ee4240d4cedb9798bb0701853a9c89673ac232f9b51bfa99c5e6784aee`.
Both stages include the exact restart and native CATKE GPU gate receipts.
The development stage uses paired v4 and the reviewed release HISTORY text.
Its scientific rows and portable score JSON are unchanged from the previous
development stage; only the render timestamp and provenance changed.
The plain-label README comparison is staged in `v018_comparison_figure_plain_v2`
with SHA-256 manifest
`286f757095b4d22aae6c161d1d1517394d899b1c1dc3a5b1377abe0347ec180c`.
Its scored fields match the earlier sealed record. The reviewed canonical-start
version-history PNG and JSON are a matched pair in
`version_history_v018_plain_v3/canonical_v018_publication_stage_63378103`.
Its `result.sha256` has SHA-256
`a9ede187a87ca900dd70c1229e9408a6381106e70605d15bb6ed0d16479a3c1e`;
the PNG has SHA-256
`226071786d0f444e59b0ac36e8011717182e61bf399c5366f8c6328a77444270`
and the matching JSON has SHA-256
`257227b407f0928d95ce06b7c414914c313575c16c8cc89370bc60de247f7e94`.
It preserves the scored inputs from the original publication stage, manifest
`04aaf8ffceb67d558366edcf807962e53e65e66901a9ea4c96f608d8ece05e75`.
The `release_status` and `publication_status` fields in the figure records give
the state when each figure was staged, before the tag. The paired and
development scorecards record the completed GPU gates. The installed PNGs
match their sealed stages byte for byte.
The maps, timelines, comparison, reviewed version-history v3 pair, and both
root scorecards are installed in this tree. The original v3/v4 metric and
season verdicts remain visible beside the neutral
adjudication.
Every README caption, validation statement, figure alt text and scorecard
that gives an ONR/CMIP6 rank must carry this scope note: these runs start
from observations. In a separate six-field test without radiation, a
climatology and an anomaly-persistence baseline beat all 59 CMIP6 models in
both seasons. The ranks measure agreement with the observed month, not
free-running climate skill. The six-field baseline record is
`initialisation_advantage_v1/baselines.json`, SHA-256
`417fe8872b678b397b1ae253d91960d313fb0295c813fdaceaff8f6836163aa0`;
it is not the full v4 release score.
The release audit checks the final documentation commit, sealed source data,
figure scope, and post-tag source and asset hashes. The paired v0.1.18 skill
outcome is the REJECT result reported above.
