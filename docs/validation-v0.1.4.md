# v0.1.4 validation

Status on 23 September 2026: candidate. The final ordinary GPU restart failed
with CUDA illegal memory access after the Fourier graph repair. A separate
sparse-exchange sanitizer finding now has a focused, tested replacement;
coupled qualification of that change is pending. The initiating cause of the
ordinary failure remains unidentified, so v0.1.4 is not tagged.

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
reported errors. Coupled checks of this source are pending. This evidence
does not establish that sparse exchange caused the ordinary CUDA700 failure.

## Scope and reproduction

The [README test commands](../README.md#test) reproduce the CPU tests, coupled
smoke, fresh-process continuation and restored-boundary comparison. Production
checks require the ERA5/ECCO input files and a CUDA device. Model construction
and the first GPU step can take hours.

These short runs do not establish a stationary or calibrated climate, seasonal
river-mouth behaviour, or equality between continuous and restarted
trajectories. The separate experimental ice-startup/fallback changes are not
included. The [known limitations](../NOTES.md#main-known-limitations) and the
labelled v0.1.1 thirty-day figure remain applicable.

Checkpoints from versions before this runoff change omit pending exchange
water and are rejected. Start a new run with v0.1.4.
