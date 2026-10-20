# Standalone atmosphere restart

Accepted on main, 20 October 2026.

The standalone SpeedyWeather restore now copies both transform-column work
arrays from the checkpoint. SpeedyWeather's generic Variables copy skipped
these nested arrays. ReadyESM's coupled restore already copied them.
The fix adds that same copy to the standalone path.

Tested source: `1d7188a`. Main preserves its source, tests, configuration,
dependencies and runners exactly.

- The focused regression passes 58 assertions. It checks direct and JLD2
  restore, repeated restore, array ownership and shared work storage.
- The same regression on the old source fails eight assertions.
- The full CPU suite passes, with the existing expected-broken
  spectral-truncation assertion.
- Two ERA5 T31/L27 probes restore every saved state entry exactly. They cover
  native pressure and the experimental centered/wider-pressure setting.

The probes also distinguish model state from observer work. An observer that
uses the model's transform storage changes four saved work arrays. The same
calculation with private storage changes none and returns identical pressure
records.

| Check | Slurm jobs |
| --- | --- |
| Focused regression | 60906292 task 0 |
| Full CPU suite | 60906292 task 1 |
| Old-source negative control | 60906292 task 2 |
| ERA5 restore and observer probes | 60922856 tasks 0 and 3 |

Result manifests verify. Failed wrappers and earlier probes remain in the
development labbook. The pressure experiments remain separate from main.
