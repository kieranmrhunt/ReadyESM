# Optional CPU convection validation

Accepted on main, 20 October 2026. Production remains v0.1.10 with GPU
Betts–Miller convection. C08 and C07 are optional CPU schemes with Float32
atmospheric storage and diagnostic clouds. See the
[configuration guide](finite-convection.md).

## Source and tests

The tested source is `5f31ed2d5ae1fa53da587e929741f97b28dc0c06`.
The integration preserves its source, tests, dependencies, vendor code,
configuration files and runners exactly. Production YAML and dependencies
also match the previous public main, `3b012af`.

- C08 matches 4,350 solver comparisons and 1,740 host comparisons.
- C07 matches 2,610 host comparisons across three parameter settings.
- The full CPU regression suite passes. Its existing spectral-truncation
  test retains one expected broken assertion.
- Six coupled pilots pass. They cover native, C08, zero-strength C07 and
  C07 at parameter cases 14, 3 and 12.
- Six fresh-process restarts pass. Every saved state matches the uninterrupted
  run at restoration and the four following steps.

The restart comparison covers 1,645 entries for native, 1,816 for C08 and
1,987 for each C07 case. It includes atmosphere, ocean, sea ice, land,
radiation, callbacks and coupling interfaces. Independent checks compare all
36 saved dictionaries, including launch and checkpoint states.

Native and C08 reproduce their qualified reference states and all 210 output
variables. Zero-strength C07 matches all 1,816 common C08 state entries at
six saved stages and all 210 output variables. Active C07 exercises momentum
transport and local heat return in every tested setting. Restart relative
energy residuals stay below 2.68e-7 and momentum residuals below 4.05e-7.

## Climate assessment

January atmosphere-only comparisons favour C08 case 12 for moisture and
lower-level temperature relative to the original C08 setting. C07 improves
global wind errors relative to matched C08, with regional trade-offs.
Production parameters stay unchanged. July and longer coupled comparisons
continue. [Experiment notes](experiments.md) retain the numerical results.

Eight-thread pilots pass the physical checks but differ between repeated
processes. Their maximum final surface-ocean temperature differences range
from 0.00331 to 0.04139 K after eight steps. This integration uses the exact
one-thread restart evidence. Threaded restart qualification continues.

## Evidence record

| Check | Slurm jobs |
|---|---|
| C07 host comparisons and full CPU suite | 60775991 tasks 0–1 |
| Clean coupled pilots and reader | 60776941 tasks 0–5, 60776959 |
| Fresh-process restarts and reader | 60782476 tasks 0–5, 60782673 |
| Eight-thread pilots and reader | 60786590 tasks 0–5, 60786591 |

The labbook retains frozen inputs, file receipts, detailed checks and failed
development attempts. Experimental runners and outputs stay outside this repo.
