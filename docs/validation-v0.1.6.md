# v0.1.6 validation

Two accepted changes: Float64 land state/water accounting, and diagnostic
liquid cloud path 60 → 75 g/m². Ice path stays at 25 g/m²; convection stays at
four hours. Runtime code is identical to the frozen land-precision candidate;
only production defaults, version metadata, tests and documentation differ.

## Coupled tests

All runs use T31/L27, the 360 × 180 × 60 tripolar ocean, dynamic sea ice and
16-layer Terrarium land. Atmosphere and ocean remain Float32.

- CPU job 56735190 passes the full candidate suite, including 45 new land
  precision/configuration and restart-history assertions.
- GPU job 56735382 starts fresh, reaches day 1, restores in a new process and
  reaches day 2. Maximum land-water residual: 2.16e-11 mm. Normal diagnostics
  and independent exact restored-state comparison pass.
- Jobs 56738876_0 and 56738876_6 branch from that same native Float64 day-2
  checkpoint with cloud60 and cloud75 respectively. Each reaches day 10,
  restores in a separate process, and reaches day 30. Both pass normal
  diagnostics and exact restored-state comparisons at both boundaries.
- An independent review recomputes 23 metrics from raw NetCDF, checks units,
  clocks and hashes, and confirms identical inherited day-0–2 histories.

The clean release checkout also passes the full CPU suite (job 57299834,
17m51s), with the new production defaults. Source identity and runtime hashes
pass before and after testing. Figure rendering passes 11 sample-contract
assertions; both images were visually checked and their hashes verified.
GitHub's independent [CPU test run](https://github.com/kieranmrhunt/ReadyESM/actions/runs/37284855061)
also passes on candidate `7f2c4fc`. The final release changes only documentation
from that tested commit.

The table gives time-weighted means over days 20–30, except labelled endpoints.
Both columns use Float64 land: this isolates the cloud change, not the combined
effect against an earlier Float32 trajectory.

| Quantity | Cloud60 control | Cloud75 |
|---|---:|---:|
| Net TOA uptake, W/m² | 21.699 | 16.562 |
| Reflected shortwave, W/m² | 101.875 | 108.240 |
| Outgoing longwave, W/m² | 217.677 | 216.449 |
| Clear-sky net uptake, W/m² | 43.504 | 43.627 |
| Surface-air temperature, K | 285.902 | 285.826 |
| Maximum ocean temperature, °C | 32.989 | 32.773 |
| Worst-cell area-equivalent ice thickness, m | 7.375 | 8.050 |
| Day-30 atmospheric water residual, kg/m² | -0.0988 | -0.1784 |
| Maximum absolute land-water residual through day 30, mm | 2.76e-10 | 2.80e-10 |

Net uptake falls 23.7%, mainly through greater reflected sunlight. Rain falls
2.1%. The atmospheric water residual and late mean ice extreme worsen; the
slightly lower day-30 ice endpoint (8.33 → 8.22 m) is not evidence of repaired
ice dynamics. These are model-to-model differences, not observational scores.

## Scope and compatibility

This is a 30-day initialized trajectory, not a stationary climate or long-run
reliability claim. Cloud75 was enabled at day 2; production enables it from
initialization. v0.1.5's 120-day validation does not qualify this changed setup.
Faster convection, reduced ice-cloud path and startup-solver experiments are
not included. No diagnostic tolerances were relaxed.

Native Float64 checkpoints are qualified. Old Float32 land checkpoints cannot
be silently converted: start fresh, use a matching Float64 checkpoint, or
explicitly retain Float32 and the original forcing to continue an old run.
The generic configuration default remains Float32 for compatibility; the
production YAML explicitly selects Float64.

The current atmospheric-water check detects non-finite residuals but has no
later-time numerical closure threshold. Passing it does not establish water
conservation of the atmosphere. Land closure is assessed separately above.

## Provenance

- Tested source manifest SHA256:
  `cac7e4fa1889d05cc954fe2322d8d20609fac93d67e15c0f08464dae55b9af28`.
- Raw-output review SHA256:
  `5b1870d28d9c7db98745adedc415392a246dd7643304969733a30fca56db3cd8`.
- [Figure provenance](readiesm-v0.1.6-snapshot.json) records the actual output
  and image hashes. Grey marks inactive cells; exactly zero sea ice is white.
