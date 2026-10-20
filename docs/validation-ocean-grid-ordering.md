# Ocean grid-update ordering

## Change

Oceananigans 0.110.15 updates the ocean surface height and vertical grid scales
in one kernel. A grid scale can read a neighbouring height before that height
has been updated. The resulting differences follow the kernel block boundaries.

ReadyESM now copies the new heights in one kernel, then calculates the scales
in a second kernel. It keeps the upstream scaling formula, launch bounds and
previous-scale assignments. The change covers Adams–Bashforth, Runge–Kutta
and vertical-coordinate reconciliation.

The override applies to hydrostatic models with mutable z-star grids.
Other model types retain the upstream methods. A version check requires
review of this override before changing the Oceananigans dependency.
The checkpoint format and configuration stay the same.

## Tests

| Test | Result |
| --- | --- |
| Isolated GPU ordering test | All 472 assertions pass |
| Package regression test | All 377 assertions pass |
| Direct packaged GPU methods | All 385 assertions pass, with 48 device-method checks |
| Full CPU suite | Passes, with the existing expected-broken spectral test |
| Compact coupled restart and repeat tests | All 20 comparisons match across 425 numeric arrays |
| Full-grid coupled restart tests | Four members match all 1680 state entries at five saved points |
| Full-grid independent repeats | All 20 comparisons match exactly |

The full-grid tests use T31/L27, the 360×180×60 tripolar ocean, land and
dynamic sea ice. Live and zero carbon exchange each have two independent
runs. Each run takes eight coupled steps and restarts from step four in
a fresh process. The largest carbon-inventory relative error is 1.31e-9.
Both uninterrupted and restored repeats match exactly.

Independent checks reproduce all 40 full-grid comparisons, covering
67,200 state-entry comparisons. Input and result checksums verify.
The accepted package source and tests match the qualified CPU candidate.
Direct GPU job 60896206 also verifies dispatch to the packaged ordered update
for both precisions and two grid patterns. Device repeats match exactly.

The first CPU test attempt stopped at an outdated model-constructor call.
Passing `grid` positionally fixed the test. The kernel code was unchanged.

The development records are GPU test 60815972, compact runs 60817815 and
reader 60817817, CPU tests 60835380, and full-grid runs 60821733 with
reader 60821734, plus packaged GPU test 60896206.

## Reproduce the package tests

```bash
julia --startup-file=no --project=. \
    -e 'using Test, ReadyESM' \
    -e 'include("test/ordered_ocean_grid.jl")'
julia --startup-file=no --project=. test/runtests.jl
```

The ordered update can change trajectories from the earlier implementation.
Use the same source revision for a continued run and its reference.
