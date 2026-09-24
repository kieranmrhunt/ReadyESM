# v0.1.4 validation

Version 0.1.4 is an incremental maintenance release of the research model.
The fixes below have component checks, CPU CI and short production evidence.
**Ordinary GPU restart reliability remains unqualified:** current-source
restart `54787077` fails with CUDA700 in its first resumed step. The earlier
[v0.1.3 notes](https://github.com/kieranmrhunt/ReadyESM/blob/v0.1.3/NOTES.md#main-known-limitations) already document the same
first-sea-ice-step illegal-access symptom. Its initiating operation and any
change in failure frequency remain unknown.

This release follows the incremental scope: publish demonstrated improvements
and retain known limitations explicitly. It does not claim all checks passed
or that the GPU fault was repaired. The [machine-readable evidence](validation-v0.1.4.json)
retains the failed restart alongside the passing tests. Ice and ocean
construction, core ice code, restart machinery and production settings checked
against v0.1.3 are unchanged; the runoff integration and GPU exchange repairs
are described below.

## Changes checked

The candidate integrates runoff removed during each native land step, retains
pending water across coupling intervals and checkpoints, and delivers water
only to unique physical ocean cells at the tripolar join. In a controlled
20 mm land-water pulse, native integration transfers 0.0044558364 m of runoff;
the previous endpoint-rate approximation represents 0.0038864274 m, about
12.8% less in that test. This is a test-specific difference, not an estimated
global bias.

Precipitation checks now respect their source Float32 precision and use the
pressure seen by convection and condensation. These diagnostic changes do not
alter the physical moisture tendencies. Output also excludes unrecorded
diagnostic capacity after early termination while preserving invalid samples
that were actually recorded.

The Fourier CUDA graph cache previously keyed only on the field pointer,
although a captured graph also contains the two scratch-buffer pointers. A
GPU test that changes scratch storage reproduced silent incorrect forward
and inverse results. The fix keys every captured buffer and its layout,
retains its Julia array owner, and waits for graph work before clearing those
owners. Initial job `54302667` passed 32 checks at T7/L2 and T31/L27, including
scratch switching, repeated view wrappers, changing values and garbage
collection. Follow-up `54308967` passed all 35 checks, including work queued
on another CUDA stream. Its control version correctly failed the new
stream-completion assertion. This identifies real transform defects; it does not establish
the cause of the separate coupled restart error.

## Qualification record

The following completed coupled checks use the Fourier graph candidate
`49d32fb`, before the sparse exchange change. Documentation-only commit
`7b050ed` also passes [CPU workflow 35843941778](https://github.com/kieranmrhunt/ReadyESM/actions/runs/35843941778).
Results are distinguished by the source and checkpoint each run used.

| Check | Result |
| --- | --- |
| Fresh GitHub environment | [Workflow 35719234028](https://github.com/kieranmrhunt/ReadyESM/actions/runs/35719234028) passed on `49d32fb`; it installs the pinned environment and runs the complete CPU suite without private ERA5 inputs. |
| Fourier GPU cache regression | Job `54308967` passed all 35 assertions. The old implementation fails the scratch-switch tests; the first cleanup correction fails the separate-stream control. |
| Repaired-source production GPU smoke | Job `54303131` passed four 900 s steps, the full diagnostic validator, NetCDF output and checkpoint creation. T31/L27 atmosphere; 360 × 180 × 60 ocean. |
| Earlier repaired-source GPU restart | Job `54303130` passed four resumed steps through hour two, full diagnostics and exact restored-boundary comparison, starting from the earlier smoke checkpoint `54135528`. |
| Final restart from the repaired smoke | Job `54683480` failed after 55 minutes on a full A100. Construction, checkpoint loading and owned runoff/alias checks passed. The first resumed step reported CUDA error 700 before the first u-momentum kernel was loaded. The initiating device operation is unidentified. |
| First allocation of the final restart | Job `54362946` failed CUDA availability before model construction. It supplies no model continuation result. |

The latest model failure occurred on the same physical GPU as the earlier
clean memcheck. This demonstrates that the graph repair has not eliminated
the coupled fault. The two new diagnostic paths check uninitialized device
memory and enable capture of GPU exception state. They retain the frozen
model and input checkpoint; their results are diagnostic, not release gates.

The earlier CPU suite passed 372 assertions in named test sets plus top-level
checks. Production CPU smoke `54174157` and CPU restart `54269734` passed
full diagnostics; the CPU restored boundary compared exactly.

All seven raw January 1993 ECCO files match the SHA-256 digests published by
the [download mirror](https://github.com/NumericalEarth/NumericalEarthArtifacts/releases/tag/data-v1).
Derived inpainting and bathymetry caches were hashed retrospectively on
22 September, then checked before the final restart. They were not
independently regenerated.

Original pre-repair restart `54166085` timed out after four hours, and
`54269500` failed with CUDA error 700 during its first resumed step.
The unchanged pre-repair source subsequently passed the full restart under
compute-sanitizer (`54318236`, zero reported errors). Instrumentation changes
execution timing; this does not establish the cause of the ordinary failure.
Original attempts, statuses and source manifests are retained.

The smoke and earlier repaired restart use frozen `qualification_graph_v5`,
with source-manifest SHA-256
`484306fcb2d09ce466f2ff7c1e757cab7ca91507f91a2e947724be13ffaef6f1`.
The final restart and focused GPU regression use `qualification_graph_v6`:
`a24c6bd74a5fc545275154b84df080a0d480fa282f156765ac7617b44d0042fc`.
Its only runtime difference is `clear_fourier_graph_cache!`, which production
construction and stepping do not call. Transform execution, cache keys and
buffer retention match v5 byte for byte. Dependency versions and physical
settings are unchanged.

## Sparse exchange follow-up

Independent CPU comparison `54754887` confirms exact equality of checkpointed
model state at the restored boundary of failed ordinary restart `54683480`.
The failure occurs during subsequent evolution.

Full restart initcheck `54722750` completes all four resumed steps and the
full diagnostic validator, but fails the sanitizer gate with 432 reports.
The first 100 printed reports are uninitialized shared-memory reads in
`cusparse::csrmv_v3_kernel`; the remaining 332 were not printed and are not
classified from that log. This is not a clean sanitizer result.

The exact production horizontal geometry provides a small independent
reproduction: a 64800 × 4608 intersection matrix with 106158 nonzeros and
its reverse operator. With initialized synthetic vectors, cuSPARSE CSR_ALG2
passes the numerical/reference checks but reports 160 shared-memory reads
under initcheck (`54766994`). Generic sparse matrices with the same dimensions
do not reproduce the reports (`54763192`, 236 assertions and zero reports).

The candidate now computes each CSR destination row in fixed storage order,
with one work item owning its accumulation and output. The isolated
production-geometry test passes its CPU-reference, repeatability and dimension
checks with zero initcheck reports (`54769921`). The public standalone test is
`test/csr_regridding_gpu.jl`; its additional Float32/Float64 edge cases include
empty rows, signed weights, cancellation and invalid vector dimensions.
The frozen integrated source has manifest SHA-256
`0f76e3eec61a85f83709225404eea279e5b98b8e24e911e68d06ed3b32a9bb7d`.
Integrated GPU regression `54775695` passes all 176 assertions separately
under initcheck (global and shared memory) and memcheck, both with zero
reported errors. Current coupled results are:

| Current runtime check | Result |
| --- | --- |
| CPU suite | [Workflow 35849068979](https://github.com/kieranmrhunt/ReadyESM/actions/runs/35849068979) passed on exact code commit `9843415`. |
| Full coupled initcheck | `54777253` completed in 57:00, with four resumed steps, full diagnostics, exact restored boundary and zero reported errors. |
| Ordinary GPU restart | `54780489` completed in 54:09, with four resumed steps, full diagnostics and exact restored boundary. |
| Fresh production smoke | `54787075` passed in 51:00: four 900 s steps, full diagnostics, NetCDF output and a one-hour checkpoint. |
| Restart of the new smoke | `54787077` failed in 2:41:14. Construction, checkpoint loading and owned-runoff checks passed; step 5 reported CUDA700 before the first u-momentum module loaded. No resumed step completed. |

The passing coupled initcheck and ordinary restart restore smoke `54303131`'s
one-hour checkpoint and advance through model hour two. The failed final
restart instead restores current-source smoke `54787075`. It used the same
physical A100 as passing ordinary restart `54780489` and the new smoke.
Its failure is detected during context synchronization before loading the
first u-momentum kernel; that kernel has not launched, and the initiating
operation is unidentified. The original log is retained with SHA-256
`5258b5ecca1804b52e9fab5a5c04aed61dc0be2c23320a51d0689b8cc42787bd`.

The earlier exception-capture run `54741841`, using unchanged v6 source,
completed without reproducing the fault or creating a GPU dump. Its repeat
`54789905` timed out after four hours while still in the first resumed step,
with no GPU dump. The timeout supplies no continuation pass or identified
fault. Independent comparison `55269611` passed exact equality of all
checkpointed model state between smoke `54787075` and failed restart
`54787077`'s saved restored boundary, excluding wall-clock bookkeeping. A full
memcheck of that actual failed checkpoint (`55270909`) completes all four
resumed steps, full diagnostics and exact restored-boundary comparison with
zero reported errors (58:07). The preceding allocation `55269692` failed CUDA initialization before model or
sanitizer execution. The passing memcheck retains frozen v7 source and physical
settings; it does not erase the ordinary failure or identify its cause.

An isolated CPU replay (`55282841`) restores the same checkpoint's ocean/ice
model state and net ice fluxes exactly, then completes one native 900 s ice
step with finite velocities and thickness, and concentration in [0, 1]. It changes
component execution, compilation and memory layout and is a diagnostic result,
not whole-model qualification. Its initial harness attempt (`55275430`) stopped
at a direct payload comparison before stepping. A second harness (`55279930`)
identified offset-array/checkpoint representation differences. The corrected
harness round-trips the restored state through the native serializer before
exact comparison and has positive and negative controls. Those earlier
harness failures did not execute a model step. The GPU replay remains a
separate investigation.

[CPU workflow 35973457753](https://github.com/kieranmrhunt/ReadyESM/actions/runs/35973457753) passes on `d2dd4a9`.
Runtime, configuration, dependencies and tests match the frozen v7 source;
this release finalization changes documentation only.

[CPU workflow 35854890920](https://github.com/kieranmrhunt/ReadyESM/actions/runs/35854890920)
passed on the README figure commit `7dfc85b`. This CPU result does not qualify
the failed GPU continuation.

The README now shows the actual hour-two surface fields from ordinary restart
`54780489`. Its [figure provenance](readiesm-v0.1.4-snapshot.json) records the
source commit, NetCDF and PNG hashes, units and plotting method. The image is
a short validation snapshot; it replaces the older v0.1.1 illustration.

## Scope and reproduction

The [README test commands](../README.md#test) reproduce the CPU tests, coupled
smoke, fresh-process continuation and restored-boundary comparison. Production
checks require the ERA5/ECCO input files and a CUDA device. Model construction
and the first GPU step can take hours.

These short runs do not establish a stationary or calibrated climate, seasonal
river-mouth behaviour, or equality between continuous and restarted
trajectories. The separate experimental ice-startup/fallback changes are not
included. The [known limitations](../NOTES.md#main-known-limitations) remain
applicable; the README snapshot does not provide a new long integration.

Checkpoints from versions before this runoff change omit pending exchange
water and are rejected. Start a new run with v0.1.4.
