# v0.1.4 validation

Status on 22 September 2026: candidate; GPU restart continuation failed with
CUDA illegal memory access. A subsequently identified Fourier graph cache
defect is repaired and passes its focused GPU regression. Fresh coupled
qualification is pending; no v0.1.4 tag has been created.

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

The CPU and coupled results below apply to the candidate **before** the
Fourier graph repair. New ordinary GPU smoke `54303131` and restart
`54303130` use the repaired frozen source and remain pending qualification.

| Check | Result |
| --- | --- |
| CPU suite without a private input-data directory | Passed: 372 assertions in named test sets, plus top-level checks. Covers runoff amounts and owned restart state, fold topology, Darcy boundaries, diagnostic sample counts, rain precision and physics-time pressure. |
| Fresh GitHub environment | [Workflow 35703695263](https://github.com/kieranmrhunt/ReadyESM/actions/runs/35703695263) passed on the current code commit `e184d64`; it installs the pinned environment and runs the complete CPU suite. Earlier candidates also passed workflows `35604989343` and `35615665722`. |
| Production GPU smoke | Job `54135528` passed four 900 s steps, the full diagnostic validator, NetCDF output and checkpoint creation. T31/L27 atmosphere; 360 × 180 × 60 ocean. |
| Independent production CPU smoke | Job `54174157` passed the same one-hour diagnostic path using `device=cpu`. |
| Production CPU restart continuation | Job `54269734` passed all four resumed steps, the two-hour endpoint, full diagnostics and exact restored-boundary comparison. |
| Raw ECCO input provenance | All seven January 1993 files match the SHA-256 digests published by the [download mirror](https://github.com/NumericalEarth/NumericalEarthArtifacts/releases/tag/data-v1). Derived inpainting and bathymetry caches were hashed retrospectively on 22 September; they were not independently regenerated. |
| Fresh GPU checkpoint restoration | Owned runoff stores and their aliasing passed. A separate process also found exact equality of checkpointed model state between the original checkpoint and the restored boundary from job `54166085`. Wall-clock bookkeeping is excluded. |
| GPU continuation from that checkpoint | Failed. MIG job `54166085` timed out after four hours. Full-A100 job `54269500` restored the owned state but failed in its first resumed step with CUDA error 700, detected before loading the first sea-ice u-momentum kernel. The preceding device operation responsible is not yet identified. |

Model source, dependencies, vendor code, production settings and CPU tests are
identical across the pre-repair candidates. Commits `25757df` and `e184d64` change only
the restart validation driver and documentation relative to the smoke-tested
`9f9d3b2`: the driver extends the resumed diagnostic horizon, reports progress
and checks every resumed clock increment. Frozen source manifests and original
failed attempts are retained in the development evidence.

The new `qualification_graph_v5` snapshot contains the graph extension fix,
its output provenance and GPU regression. Its source-manifest SHA-256 is
`484306fcb2d09ce466f2ff7c1e757cab7ca91507f91a2e947724be13ffaef6f1`.
Dependency versions and physical settings remain unchanged. Earlier coupled
passes are not counted as qualification of this changed GPU code.

The final cache-cleanup correction is frozen in `qualification_graph_v6`, with
source-manifest SHA-256
`a24c6bd74a5fc545275154b84df080a0d480fa282f156765ac7617b44d0042fc`.
Its only runtime difference from v5 is `clear_fourier_graph_cache!`, which
production construction and stepping do not call. All transform execution,
cache keys and buffer retention match v5 byte for byte. The separate stream
regression checks this utility change; the ongoing ordinary coupled gates
use v5. CPU [workflow 35714709883](https://github.com/kieranmrhunt/ReadyESM/actions/runs/35714709883)
passed on `7e6a1d9`, before the cleanup follow-up.

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
