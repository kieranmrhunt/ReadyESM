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

January and July atmosphere-only comparisons favour C08 case 12 for moisture and
lower-level temperature relative to the original C08 setting. C07 improves
global wind errors relative to matched C08, with regional trade-offs.
Production parameters stay unchanged. Longer coupled comparisons continue.
[Experiment notes](experiments.md) retain the numerical results.

## Eight-thread restarts

After the ocean grid-ordering fix, all six eight-thread members pass on
source `e9e50b4`. Native, C08 and C07 case 12 each run twice for eight coupled
steps, with a fresh-process restart at step four.

All saved state entries match at restoration and each following step.
Independent repeats match at all six saved stages and across all 210 output
variables for each case. Both reviewers verified the reader receipt.

Earlier eight-thread failures remain in the labbook. They first diverged
in ocean and coupling fields at step five. The ordering repair removes the
observed restart and repeat differences in this test.

## Thread-count comparison

One-thread runs of native convection and C07 case 12 match both eight-thread
repeats on `e9e50b4`. Each comparison covers the initial state, steps four
through seven and the final state after eight full-grid coupled steps.
All 1,645 native entries, 1,987 C07 entries and 210 output variables match
exactly. The reader and result receipt pass. Thread count remains fixed
within each longer experiment.

## Evidence record

| Check | Slurm jobs |
|---|---|
| C07 host comparisons and full CPU suite | 60775991 tasks 0–1 |
| Clean coupled pilots and reader | 60776941 tasks 0–5, 60776959 |
| Fresh-process restarts and reader | 60782476 tasks 0–5, 60782673 |
| Eight-thread pilots and reader | 60786590 tasks 0–5, 60786591 |
| Eight-thread restarts and repeats after ocean fix | 60876558 tasks 0–5, 60876559 |
| One-thread against eight-thread comparison | 60896182 tasks 0–1, 60896186 |

The labbook retains frozen inputs, file receipts, detailed checks and failed
development attempts. Experimental runners and outputs stay outside this repo.
