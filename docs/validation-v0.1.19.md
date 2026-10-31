# v0.1.19 candidate validation record

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
rejecting missing, extra or changed fields. Its focused tests passed 80/80;
Claude and Barry independently read the fix. Real Float32 and production-Float64
spectral-albedo CPU save/restore jobs `63738587` and `63738488` both passed
with bit-identical continuations. Exact-source T31 January/July jobs
`63738594_0/_1` both passed. The full CPU suite `63738471` is running.
No GPU restart or paired candidate array has been submitted for `4caab68`.

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
`2f8e9193e2fea340e81dacc0e486c1d7de7f3b3a1b1e44dfc16ae3eef76f5c0c`;
the vegetation input SHA-256 is
`b80d7e01a4390340794030f1f924e5583d813fb163a1f9d875ae749091c13f60`.
Both are packaged in `data/` and the albedo derivation, its erratum and limits
are documented there.

## Evidence required before release

| Gate | Current state |
| --- | --- |
| Exact packaged-input hashes and effective production config | The six exact-source `4caab68` paired configuration tasks `63738777_0..5` ended Slurm COMPLETED 0:0 with strict-sealed PASS receipts. Their aggregate receipt SHA-256 is `aee5a21ee925961d42faa1a7c543e9ca017864fbe5b696e18969a68045709634`; the actual paired runner's January/July configuration checks `63738895_0/_3` also passed. The final tag must retain these packaged data and production settings. |
| Component tests for albedo, vegetation and dissipative heating | Albedo passes 69/69, vegetation 40/40 and corrected D1 148/148 on the equivalent substantive source. The first D1 test used an outdated fixed input. Its replacement passed. These focused results do not substitute for the exact-source full suite. |
| Full CPU suite | The exact-source job `63738471` is running on `4caab68`. Slurm COMPLETED 0:0 and a strict-sealed `READYESM_CPU_SUITE_PASS` are required before A100 restart. Earlier-source receipts do not qualify this commit. |
| Six-hour CPU smoke with all three options active | The 90-day full-set run's preflight `63667747` and smoke `63667748` both completed 0:0 on source `d325fc3`. Their verified result manifests have SHA-256 `fcf8a958ab9749e93ae6a859245b7330be3f77f86ed42e2d0af4904056b1b191` and `bc182579ae9623e6ccf4d80dfa91201fd2559730389e19d372eb157dbe37365f`. The 24-step run and its independent output-file check passed. A separate six-step check measured positive D1 drag and hyperdiffusion and 183 mapped columns for each land treatment. An independent review cleared the 90-day launch. |
| Same-start 90-day July full-set run versus the v0.1.18 baseline | The one-member 90-day run (internal label V5b), job `63690775`, completed 0:0 with a strict-verifying 38-entry result-manifest SHA-256 `9d3521208650ee9c6b134986e638800215340c949e958bed14d62c938d561c0f`, 15-field same start and output-file checks. The reviewed V5b reader v2 reported `treatment_applied`: 40–60 N land reflected shortwave fell 49.92 W m⁻² all sky and 62.61 W m⁻² clear sky against V2. Its predeclared 30–60 N 500 hPa weeks-4–12 ERA5 error moved from −3.31 K to +0.05 K, a +3.36 K change above the 0.6 K detection threshold. This is one combined-treatment start, not an attribution or six-member skill verdict. Its tested source `d325fc3` differs from `93db53a` under `src/` only by a trailing blank line, and the selected production physics and two packaged data hashes match. D1 is selected in the loaded GPU config but its budget was not measured in this run; the CPU trio checks activity, and later paired months require live budget data. |
| Energy, vertical-profile, land/ocean and regional guardrails | The V5b reader reports global net top-of-atmosphere changes of +9.71 W m⁻² in July, +7.46 in August and +6.65 in September against V2. Direct July comparison with the frozen v4 ERA5 TOA reference finds V5b reflected shortwave bias −0.97 W m⁻² and outgoing-longwave bias −6.25 W m⁻²; V2's respective biases were +10.60 and −8.11 W m⁻². Thus V2's near net-TOA agreement combined two large compensating errors, while V5b's net bias is +7.95 W m⁻². Global 250 hPa ERA5 errors are +1.29, +1.76 and +1.64 K in V5b, versus +0.49, +0.04 and +0.38 K in V2. Other monthly regional, land/ocean, precipitation and surface-flux rows are in the sealed readout `v019_v5b_independent_readout_codex_v1` (SHA256SUMS SHA-256 `089b72ba4959843c0516c44043752d5083b73d610aa96bf2e95a5ef0a4571b19`). These are reported guardrails, not an energy-closure pass or a repeatability result; independent scientific disposition remains open. |
| January and July paired release months, scores and safety checks | Six matched v0.1.18 controls have strict-sealed PASS receipts. Exact-source six-member configuration array `63738777_0..5` and actual runner configuration jobs `63738895_0/_3` passed. Candidate launcher input seal `f94f01b9537f99ff532bf00da2e994c1dff83434f40ec2d03b774a71ffb5950d`; scorer seal `44215f42b7a7da81304ab19d733dbc26e27b7e2340128548c9a4e28aea59fb7e`. Six candidate GPU months, treatment gate, frozen v3/v4 scores, neutral rule and safety checks remain pending until the restart PASS and a fresh vote. |
| GPU restart/output checks, final figures and two independent scoped reviews | Real `4caab68` CPU spectral-albedo round trips passed: Float32 `63738587` (result SHA-256 `5d3d1ba2147c0159b32e11a13f2431b41893b1a6d889aac2dc13631647b33775`) and production Float64 `63738488` (`92f99f79ccdf0e7c7d57f0f489c0ccd71e76ec6d678eaf2715f3062b7f237576`). T31 January/July `63738594_0/_1` passed with strict result SHA-256 values `82f6770b27721304bf384af9eccfb52a3bd877b8c7f17ff04c4bdea1533ba503` and `5a7d410ca8daee608a2d58a6d99c6abc4c5647e0e67f8f50d1bd7a4d3469e5f5`; both reports show eight CPU T31 steps, Float64 land, 2,123 mapped columns, cover in [0,1] and production diagnostics PASS. A100 restart input seal `2df0ff60e4ba267b0b30aac2ff5c36efc5020218ef6ff07859a02c6c940f4772` and pre-submit verifier seal `0961ef9b3085a7cdab90c1fb2de6ed78538ec09c6a881d24e36c0c8181d79bbc` are staged. No A100 restart job has been submitted; full CPU and pre-submit PASS remain prerequisites. Final figures and reviews are pending. |
| Drift-screen calculation unit tests | Nine Python `unittest` cases pass on the tag-descended candidate. |
| Final clean tree, tag, remote readback and post-tag figure audit | The separate preliminary docs integration checkout is clean at `18c0c46`, descended from tested `4caab68`, and changes only four prose paths. It is not the final release commit. After scored figures and validation are integrated, re-run path/subtree equivalence and a real CPU save/restore on the exact final-tag commit; then verify clean tree, tag, remote refs and all figure provenance. |

The staged treatment gate has also been checked against all six sealed
v0.1.18 control months used as unchanged candidates. Each gave zero reflected
shortwave change and failed the predeclared treatment thresholds. This is a
negative control for the gate, not a v0.1.19 climate result. A code read
confirmed that the staged scorer evaluates this gate before the frozen v3/v4
skill scores and keeps a scientific `REJECT` distinct from a failed treatment
gate; execution with candidate months remains pending.

The v0.1.18 comparison figure binds a separate GPU guard-wrapper receipt to
its exact candidate source. The v0.1.19 figure must bind any corresponding
runtime claim to an exact-source v0.1.19 receipt or describe the verified
restart and production-diagnostic checks instead. No v0.1.18 guard-wrapper
receipt can substantiate a v0.1.19 claim. The v0.1.19 figure recipe and final
figure remain pending.

The [v0.1.18 record](validation-v0.1.18.md) is the previous release baseline.
The one-member V5b screen removed its measured 30–60 N 500 hPa cold bias.
Repeatability, whole-model skill and the long-run energy budget remain open.

## Experimental snow and canopy work

A land snowpack and a store for rain held on leaves are being developed on
separate branches. A combined snow-energy prototype passed a 24-step real-snow
CPU coupled smoke and restart integration test on merge commit `3cc12d2`
(sealed result-manifest SHA-256 `90bbed357f6863a51bb9d7fa35e8b8737fbbbc2a4fe2cec8b97dbde6c67070de`).
These arms are not part of v0.1.19 and were not used in any release result in
this record. Month-scale energy closure remains untested.
