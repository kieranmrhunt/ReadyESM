# v0.1.19 candidate validation record

## Paired score decision, 1 November

All six fresh exact-source January/July candidate months in array `63776410`
ended Slurm COMPLETED 0:0 with strict `V019_PAIRED_MONTH_PASS` receipts. The
frozen treatment gate passed. The final scorer job `63794836` ended Slurm
COMPLETED 0:0, and its independent terminal verifier returned
`V019_PAIRED_SCORE_VERIFIED REJECT`. Its result-manifest SHA-256 is
`4110bafb73c4a8a7517c6dfad48edc001edf99382f2d2d02ddcbcbb760667c64`;
its frozen input-manifest SHA-256 is
`0a691e8e6674e83d216378ba5da0969a6caec342c6ba2a738e7f4f598a92161f`.
The two intervening scorer failures were output-receipt bookkeeping errors;
their traces remain in the sealed successor bundles. The final scorer uses
the original frozen v3/v4 metrics, treatment thresholds, safety gate and
neutral adjudication rule.

The v3 metric gate was `UNDECIDED`; the v4 metric gate and neutral release
adjudication were `REJECT`. Both seasons improved on the nine-field score,
and the albedo and outgoing-longwave guards passed. The two-season absolute
net top-of-atmosphere bias increased from 0.556 to 6.974 W m⁻², a
6.418 W m⁻² deterioration against the frozen 0.5 W m⁻² allowance. Signed
net biases were positive in all three candidate January months (about
+6.0 W m⁻²) and all three July months (about +8.0 W m⁻²). The control's
small two-season number includes cancellation between January positive
(about +2.9 W m⁻²) and July negative (about −1.8 W m⁻²) biases.

The physical breakage gate also failed: January day-31 sea-ice area exceeded
the paired allowance in members 1 and 3. January and July day-31 ocean
surface temperature, and July sea-ice area, passed. The separate NSIDC
monthly-area readout showed January observed 17.314, control 18.381 and
candidate 18.910 million km²; total monthly bias therefore worsened from
1.067 to 1.595 million km². July total monthly bias improved slightly, from
2.582 to 2.557 million km². The observation comparison is reported
separately and does not replace the predeclared day-31 gate.

The day-31 area increase is concentrated in Southern Hemisphere summer ice:
control areas of 1.61–1.64 million km² become 2.68–2.70 million km² in all
three pairs, an increase of 1.04–1.08 million km². Northern day-31 changes
have mixed signs (+0.94, −0.25 and +0.08 million km²). This pattern is
consistent with less Antarctic summer melt, but the paired experiment changes
several processes together and does not attribute the cause. NSIDC's January
Southern Hemisphere **monthly mean** is 3.12 million km²; comparing that
with model **day 31** is only directional, not a matched observational score.
The monthly global total comparison above uses matching time windows and
shows a larger candidate bias. Both views belong in a release discussion;
neither changes the frozen physical-breakage verdict.

**Release status: HOLD.** The frozen net-radiation and January sea-ice
guards both fail. LONG1 drift evidence and a scientific disposition are
still pending. No v0.1.19 tag or GitHub release has been made; scored figures
are staged as recipes but are not installed in the release checkout.

## Superseding restart finding, 31 October

The `8ad3df7` candidate and its frozen GPU restart, paired launcher and scorer
seals are on **HOLD**. A real step-12 CPU checkpoint from the production
spectral-albedo profile failed to restore. Claude's sealed probe, job
`63738389` (Slurm COMPLETED 0:0, `SHA256SUMS` SHA-256
`d41acdf4dcb63418b691efa6bd7082948551a880357991cc92bc8db503551428`),
loaded the saved checkpoint through the production loader and found the
seasonal-albedo record fields in the order
`(:has_solved, :direct, :diffuse, :identity)`. The old validator required the
same four fields in a different order. Its identity values and the vegetation
identity survived the save and load.

The narrow successor is commit `4caab6896ee4780f9fe283d03d955407bdfa3748`,
tree `b6e94f088272a99868d333111118f284ce420b59`, descended from `8ad3df7`.
It compares exact seasonal-albedo fields and identity values by name, while
rejecting missing, extra or changed fields. Its focused tests passed 80/80.
Claude and Barry independently read the fix. Real Float32 and production-Float64
spectral-albedo CPU save/restore jobs `63738587` and `63738488` both passed
with bit-identical continuations. Exact-source T31 January/July jobs
`63738594_0/_1` both passed. The full CPU suite `63738471` ended Slurm
COMPLETED 0:0 with a strict-sealed `READYESM_CPU_SUITE_PASS` receipt.
The one-A100 restart test `63764225` ended Slurm COMPLETED 0:0. Its 55-file
strict result manifest has SHA-256
`f59c79ed0fd097ca4c1c243c3db6b4518d3a11c3687fdfefa07d5e4663634a47`;
the bundled independent verifier returned `V019_RESTART_GPU_VERIFIED`. The
first exact-source paired candidate array `63769557` began on 31 October.
Its three January members passed. July member 1 failed during the input
check, before the model started, because an earlier test run had written
working files derived from the July ocean start files into the shared input
folder. Those five files were moved, with hashes, to
`v019_july_ecco_cache_quarantine_codex_v1`. The seven raw inputs and receipt
were unchanged, and the input check passed again. Pending
July tasks 2 and 3 from the first array were cancelled. A fresh six-task array
`63776410` uses the unchanged frozen launcher, source, inputs and scorer. Its
dependency cleared after the first array's January member 3 finished Slurm
COMPLETED 0:0 with a strict-verifying `V019_PAIRED_MONTH_PASS` receipt and a
passing tropical-land treatment check. Fresh January tasks 0 and 1 ended
Slurm COMPLETED 0:0 on 1 November with strict-verifying
`V019_PAIRED_MONTH_PASS` receipts; their result-manifest SHA-256 values are
`eb046c073f931a9a088e6858c1d49df9e68879e0a87498111c378aae9b901351`
and `75394848e73078de786b273d7f61b0854b1cc77ec7e2685978441f5c6155dcdf`.
Fresh January task 2 also passed on 1 November, with result-manifest SHA-256
`80905a2536d0eba5d42007d791e10f948adde88eced734317d7ebb8d5091b69d`.
Fresh July task 3 (member 1) passed with result-manifest SHA-256
`7cf3b18d45baa1176af59fc904bf24aa475776a66a29ef5b4aed4117a1f7daf9`.
Only the fresh array entered the six-pair score. The failed first attempt
stays in the audit record; the terminal scored result and release hold are
recorded above.

This is a working record for the release candidate. It is not a release verdict.
The earlier local candidate is commit `8ad3df733949038fda6cef68065c30cadd90bffb`,
tree `4163431e7fb6a84d760fbee01738b0771f45647b`, descended from
the published v0.1.18 tag through the earlier candidate `93db53a`. The final tagged source and its exact inputs must
be checked again before publication.

## Production selection

The versioned `config/production.yml` profile selects atmospheric
`kinetic_energy_v1` dissipative heating, `era5_spectral` land albedo and
`cover_lai` seasonal vegetation. The generic `ExperimentConfig` defaults
remain `off`, `constant` and `off` so ad hoc configurations do not silently
require the data-backed options. The albedo input SHA-256 is
`2f8e9193e2fea340e81dacc0e486c1d7de7f3b3a1b1e44dfc16ae3eef76f5c0c`.
The vegetation input SHA-256 is
`b80d7e01a4390340794030f1f924e5583d813fb163a1f9d875ae749091c13f60`.
Both are packaged in `data/` and the albedo derivation, its erratum and limits
are documented there.

## Evidence required before release

| Gate | Current state |
| --- | --- |
| Exact packaged-input hashes and effective production config | The six exact-source `4caab68` paired configuration tasks `63738777_0..5` ended Slurm COMPLETED 0:0 with strict-sealed PASS receipts. Their aggregate receipt SHA-256 is `aee5a21ee925961d42faa1a7c543e9ca017864fbe5b696e18969a68045709634`. The actual paired runner's January/July configuration checks `63738895_0/_3` also passed. The final tag must retain these packaged data and production settings. |
| Component tests for albedo, vegetation and dissipative heating | Albedo passes 69/69, vegetation 40/40 and corrected D1 148/148 on the equivalent substantive source. The first D1 test used an outdated fixed input. Its replacement passed. These focused results do not substitute for the exact-source full suite. |
| Full CPU suite | Exact-source job `63738471` ended Slurm COMPLETED 0:0. Its strict result-manifest SHA-256 is `a75d0a2f16dcc107d97172836e92c90e6260877a0289fc280ac1b02fda94b125`, with status complete, exit code 0, source `4caab68`/tree `b6e94f08` and verdict `READYESM_CPU_SUITE_PASS`. The frozen launcher runs the complete `test/runtests.jl`, including D1, seasonal albedo and vegetation tests. |
| Six-hour CPU smoke with all three options active | The 90-day full-set run's preflight `63667747` and smoke `63667748` both completed 0:0 on source `d325fc3`. Their verified result manifests have SHA-256 `fcf8a958ab9749e93ae6a859245b7330be3f77f86ed42e2d0af4904056b1b191` and `bc182579ae9623e6ccf4d80dfa91201fd2559730389e19d372eb157dbe37365f`. The 24-step run and its independent output-file check passed. A separate six-step check measured positive D1 drag and hyperdiffusion and 183 mapped columns for each land treatment. An independent review cleared the 90-day launch. |
| Same-start 90-day July full-set run versus the v0.1.18 baseline | The one-member 90-day run with all three production changes, job `63690775`, completed 0:0 with a strict-verifying 38-entry result-manifest SHA-256 `9d3521208650ee9c6b134986e638800215340c949e958bed14d62c938d561c0f`, 15-field same start and output-file checks. The independently reviewed run reader reported `treatment_applied`: 40–60 N land reflected shortwave fell 49.92 W m⁻² all sky and 62.61 W m⁻² clear sky against the matched run with v0.1.18 settings. Its predeclared 30–60 N 500 hPa weeks-4–12 ERA5 error moved from −3.31 K to +0.05 K, a +3.36 K change above the 0.6 K detection threshold. This is one combined-treatment start, not an attribution or six-member skill verdict. Its tested source `d325fc3` differs from `93db53a` under `src/` only by a trailing blank line, and the selected production physics and two packaged data hashes match. D1 is selected in the loaded GPU config but its budget was not measured in this run. The CPU trio checks activity, and later paired months require live budget data. |
| Energy, vertical-profile, land/ocean and regional guardrails | The independent readout reports global net top-of-atmosphere changes of +9.71 W m⁻² in July, +7.46 in August and +6.65 in September against the matched run with v0.1.18 settings. Direct July comparison with the frozen v4 ERA5 TOA reference finds the combined run's reflected shortwave bias −0.97 W m⁻² and outgoing-longwave bias −6.25 W m⁻². The matched v0.1.18 run's respective biases were +10.60 and −8.11 W m⁻². The matched v0.1.18 run's near net-TOA agreement combined two large compensating errors, while the combined run's net bias is +7.95 W m⁻². Global 250 hPa ERA5 errors are +1.29, +1.76 and +1.64 K in the combined run, versus +0.49, +0.04 and +0.38 K in the matched v0.1.18 run. Other monthly regional, land/ocean, precipitation and surface-flux rows are in the sealed readout `v019_v5b_independent_readout_codex_v1` (SHA256SUMS SHA-256 `089b72ba4959843c0516c44043752d5083b73d610aa96bf2e95a5ef0a4571b19`). These are reported guardrails, not an energy-closure pass or a repeatability result. Independent scientific disposition remains open. |
| January and July paired release months, scores and safety checks | All six exact-source fresh months in `63776410` passed with strict receipts. The treatment gate passed. Scorer `63794836` completed and its frozen terminal verifier returned `V019_PAIRED_SCORE_VERIFIED REJECT` (result-manifest SHA-256 `4110bafb73c4a8a7517c6dfad48edc001edf99382f2d2d02ddcbcbb760667c64`). The v3 gate was UNDECIDED; the v4 and neutral gates rejected on the two-season net TOA guard and January sea-ice physical breakage. Full values and the retained bookkeeping-failure audit are above. |
| GPU restart/output checks, final figures and two independent scoped reviews | Real `4caab68` CPU spectral-albedo round trips passed: Float32 `63738587` (result SHA-256 `5d3d1ba2147c0159b32e11a13f2431b41893b1a6d889aac2dc13631647b33775`) and production Float64 `63738488` (`92f99f79ccdf0e7c7d57f0f489c0ccd71e76ec6d678eaf2715f3062b7f237576`). T31 January/July `63738594_0/_1` passed with strict result SHA-256 values `82f6770b27721304bf384af9eccfb52a3bd877b8c7f17ff04c4bdea1533ba503` and `5a7d410ca8daee608a2d58a6d99c6abc4c5647e0e67f8f50d1bd7a4d3469e5f5`. Both reports show eight CPU T31 steps, Float64 land, 2,123 mapped columns, cover in [0,1] and production diagnostics PASS. A corrected task-local ECCO rerun, `63777419_0/_1`, also ended Slurm COMPLETED 0:0 for both seasons. Its prewritten independent reader strict-verified both result manifests, eight-step Float64 T31 reports, source identity, local raw input/cache inventory and the clean read-only shared July raw folder; it returned `V019_T31_TASK_LOCAL_ECCO_MODEL_PASS` with sealed result artifact `v019_t31_tasklocal_readout_codex_v1` (SHA256SUMS SHA-256 `1221a5fcb86e0d8be3d35f6a37f983efa8cf78e1dc30728b955ce34615f1b4aa`). The frozen five-receipt pre-submit verifier returned PASS with JSON SHA-256 `7c10e0bc8e720ab9de0575e2c0c2af3cf701dc935f57908e0a12f38b8e09a18e`. The exact-source A100 restart input seal is `2df0ff60e4ba267b0b30aac2ff5c36efc5020218ef6ff07859a02c6c940f4772`. Job `63764225` passed with a strict 55-file result-manifest SHA-256 `f59c79ed0fd097ca4c1c243c3db6b4518d3a11c3687fdfefa07d5e4663634a47` and independent verifier PASS. The paired output is verified and rejects under the frozen rule. README figures, all three root scorecards and release disposition remain pending. |
| Matched one-year LONG1 segment A | Same-start T1a and T1b jobs `63776468_0/_1`, both initialized on 1 January 1993 and covering January–June in segment A, were RUNNING at 00:33 UTC on 1 November. Their terminal strict reports, global energy and temperature drift, and next-segment decision remain pending. This long-run readout is distinct from the six-pair release score. |
| Drift-screen calculation unit tests | Nine Python `unittest` cases pass on the tag-descended candidate. |
| Final clean tree, tag, remote readback and post-tag figure audit | The separate preliminary docs integration checkout descends from tested `4caab68` and currently changes only documentation. It is not the final release commit. After scored figures and validation are integrated, re-run path/subtree equivalence and a real CPU save/restore on the exact final-tag commit. Then verify clean tree, tag, remote refs and all figure provenance. |

## Supplemental January sea-ice transient, exactly repeated in fresh member 1

The completed first-array January member 1 had a **0.711936 m** rise in the
daily maximum area-equivalent sea-ice thickness, `max(h × A)`, from day 22 to
day 23. The sealed supplemental reader
`v019_ice_transient_readout_codex_v1` (SHA256SUMS SHA-256
`713e29f7aa9224ab44a006076980d658703e270a1529fb4281a8e7d93729fda5`)
found that global ice area and volume changed on their ordinary daily scale.
On day 31, Claude's sealed locator `jan_m1_ice_extreme_claude_v1`
(SHA256SUMS SHA-256
`98dd19c67df08ad274b0d16f8ef6470e3cd2bcad9c3f69a3afcb7cde02ddb12e`)
placed the maximum `h × A` of **4.916815 m** at 80.212° N, 92.022° E on the
Severnaya Zemlya coast; the saved day-31 scalar agrees exactly with that map
maximum. The cell had thickness 5.218 m and concentration 0.942. The day-22
and day-23 outputs save only the daily global-maximum scalar, so the jump
cannot be assigned to the day-31 cell or to a specific local mechanism. The
fresh January member-1 task `63776410_0` had the same 15-field initial state
and **exactly identical saved daily values** for maximum `h × A`, extent-maximum
thickness, global ice area and volume. The strict-sealed repeat reader
`v019_ice_repeat_jan_m1_codex_v1` has SHA256SUMS SHA-256
`b52700ea01ab971903914eab4085550f62435840d44fd745782c432e3f6b5e01`.
The +0.711936 m step is therefore repeatable in this configuration; its
location and physical cause remain unresolved. This is a post-hoc observation,
separate from the frozen paired score and its predeclared physical-breakage
gate.

The prewritten full three-member January repeat check completed after fresh
member 3 (`63776410_2`) reached a strict-sealed PASS receipt (result-manifest
SHA-256 `80905a2536d0eba5d42007d791e10f948adde88eced734317d7ebb8d5091b69d`).
The final report `v019_repeat_jan_final_codex_v1` (SHA256SUMS SHA-256
`c6888a92166805e2d292b239392f20d53bcbbd635fa772d6fa758bab466a47ee`)
compared fresh array `63776410` with first array `63769557`. All 15 initial
state hashes match for each member. All four scored NetCDF files are
byte-identical for all three pairs. Every one of 210 diagnostic variables and
its metadata match by shape, dtype and data bit pattern; the sole differing
global attribute is the wall-clock `created_at_utc`. The report returns
`all_reproducible=true`. This is a repeatability result, separate from the
frozen six-member skill score and treatment gate.

The staged treatment gate has also been checked against all six sealed
v0.1.18 control months used as unchanged candidates. Each gave zero reflected
shortwave change and failed the predeclared treatment thresholds. This is a
negative control for the gate, not a v0.1.19 climate result. A code read
confirmed that the staged scorer evaluates this gate before the frozen v3/v4
skill scores and keeps a scientific `REJECT` distinct from a failed treatment
gate. The candidate treatment passed in all six fresh months; the final scorer result is recorded above.

The v0.1.18 comparison figure binds a separate GPU guard-wrapper receipt to
its exact candidate source. The v0.1.19 figure must bind any corresponding
runtime claim to an exact-source v0.1.19 receipt or describe the verified
restart and production-diagnostic checks instead. No v0.1.18 guard-wrapper
receipt can substantiate a v0.1.19 claim. The v0.1.19 figure recipe and final
figure remain pending.

The [v0.1.18 record](validation-v0.1.18.md) is the previous release baseline.
The one-member combined-treatment screen removed its measured 30–60 N 500 hPa cold bias.
The three January members repeated exactly in their saved scientific fields.
The six-pair release adjudication is REJECT under the frozen rule; the long-run energy budget remains open.

## Known large-scale snow energy limitation

The production large-scale precipitation path inherits a missing snow-melt
cooling term from the pinned SpeedyWeather `ImplicitCondensation` routine. The
routine computes a melt cooling amount but overwrites it before adding its
temperature tendency, while retaining the heat released when snow forms. Its
intermediate temperature estimate also multiplies that per-step cooling amount
by the timestep as though it were a tendency. The released-source candidate
therefore has an energy-budget limitation in snow-forming and melting columns.

Claude's corrected sealed saved-profile audit `lsc_melt_energy_claude_v2`
(SHA256SUMS SHA-256 `861b18e6617e028ec66cfafb44b8cffd486cef8ffe1943b79dd9dce7a5a5a989`)
estimated the global-mean large-scale column heating in excess of latent heat
carried by surface snow as 0.803 W m⁻² for exact-source January member 1 and
0.628 W m⁻² for the A1 July member 1. An independent sealed comparison
(`v019_a1_july_equivalence_codex_v1`, SHA256SUMS SHA-256
`10b6c41ff868c8b363d1e55b775e2b6c82effcdc7a38a18526bfec20d4792450`)
found its four scored files and all 210 diagnostic variables identical to the
fresh exact-source July member 1; only job metadata differ. These are
approximate upper bounds on missing melt cooling: the
calculation uses global-mean surface pressure for layer masses, and the
net-column wrapper may reduce diagnosed snowfall without changing the
thermodynamic tendencies. The corrected readout uses the pinned model's
330 kJ kg⁻¹ fusion constant; its v1 predecessor used 334 kJ kg⁻¹. The
convection column provided a zero-residual control under the same calculation.
A direct melt diagnostic and a controlled correction are being developed for a
later version. The v0.1.19 paired score and LONG1 readout must not be presented
as a complete atmospheric energy closure result.

## Experimental snow and canopy work

A land snowpack and a store for rain held on leaves are being developed on
separate branches. A combined snow-energy prototype passed a 24-step real-snow
CPU coupled smoke and restart integration test on merge commit `3cc12d2`
(sealed result-manifest SHA-256 `90bbed357f6863a51bb9d7fa35e8b8737fbbbc2a4fe2cec8b97dbde6c67070de`).
These arms are not part of v0.1.19 and were not used in any release result in
this record. Month-scale energy closure remains untested.
