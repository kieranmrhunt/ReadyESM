# Package precompilation on the v0.2 development line

The public v0.1.19 tag has a proposed code-only bounded PrecompileTools port at
`aacf3f02fd8e7d7f4c46d88e1165a0e6d1baf29b`. Antigravity developed the
workload; this public-tag port applies the same code change to the published
release lineage. These measurements concern the proposed port, not a new
published release or the production-resolution timings in the README.

## Choosing a mode

`READYESM_PRECOMPILE_WORKLOAD_MODE` supports `bounded` (the default), `off`,
and experimental `full`. The older `READYESM_PRECOMPILE_WORKLOAD` variable is
used only when `_MODE` is unset.

- `bounded` loads and validates small analytic configurations, constructs a
  T7 spectral grid, exercises selected atmosphere helpers, and requests four
  constructor/runner precompile signatures. It does not build or step a
  coupled model during package precompilation.
- `off` omits this ReadyESM workload from a newly built package cache.
- `full` additionally builds and runs a one-step analytic coupled model and
  collects diagnostics. It is experimental and outside the qualification
  reported below.

**Set the mode before ReadyESM's compiled cache is created.** An ordinary
`using ReadyESM` can reuse an existing cache, so changing the environment
variable afterward does not change the workload stored in that image. To
change modes, rebuild ReadyESM's compiled image or use a fresh depot with no
fallback compiled ReadyESM image. Prepending an empty depot while leaving an
older compiled-cache depot in the search path does not guarantee a mode change.
Package sources and artifacts can still be supplied from a separate prepared
store.

## Measured analytic CPU path, 5 November 2026

The sealed public-tag pair compared clean control `816fc17` with the code-only
port `aacf3f0`. Each arm created a fresh compiled image, loaded ReadyESM,
constructed a T7 atmosphere with a 32×16×4 idealized ocean, took two coupled
steps, saved diagnostics, and built a second model. The arms used different
AMD EPYC 9654 hosts. All 17 initial and 17 final prognostic hashes, 203 saved
NetCDF arrays with dimensions and dtypes, and four diagnostic scalars matched
bit for bit. All listed module-global caches were empty after import. The pair
input SHA-256 was `ea49f6f6887da069a2ae0939c250ee1d581c245d092fec6fe4d41f5a8c1811fd`;
the control and port terminal result-manifest SHA-256 values were
`f852bdeae122b9562bbf4512b16a602ebb131ef98e662ba72c8658acb51f4528`
and `edba819054583b1b7a7a1bf2b6e702a3824358ec9f9e93dce1037b0ee259a0d7`.

| Observed phase | Public control | Bounded port |
| --- | ---: | ---: |
| Package cache creation | 180.93 s | 479.04 s |
| Module load | 10.26 s | 24.45 s |
| First model build | 709.71 s | 383.87 s |
| Two-step run | 362.88 s | 365.94 s |
| Diagnostic collection and save | 91.54 s | 88.91 s |
| Full first analytic path after cache creation | 1180.57 s | 863.30 s |
| Second build after run and diagnostics | 235.35 s | 376.69 s |
| Native ReadyESM cache image (`.so`) | 6.47 MB | 172.88 MB |

The table reports what happened on those two hosts. It does not establish a
causal speedup or a project-wide default for compute time. Cache creation and
image storage cost more in the port. The second build after a run and
diagnostic save was slower in both this pair and an earlier candidate pair;
a same-node build-only test did not reproduce that penalty. Antigravity's
staged run/diagnostic profile remains needed before changing the workload on
performance grounds.

For repeated jobs on one immutable checkout path, a verified cache built for
that exact checkout can be reused. A candidate-specific test found no
additional package compilation on same-path reuse and verified that its
20-file copied depot remained byte-identical after ReadyESM import. That
result does not establish cache portability to a relocated checkout, a
different source commit, or a GPU node.

The code-only public port has passed the two-step analytic CPU parity gate.
Its complete public CPU suite, exact-port A100 smoke, `off`-mode fallback and
combination with other v0.2 physics changes require separate terminal results
before release. This guide will be updated from those receipts before
publication.
