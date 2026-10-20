# Optional land support

ReadyESM now supports wet-atmosphere runs with `land=nothing`. Two wrapper
paths previously assumed that the land namespace existed. Both now check
for it before accessing Terrarium state. Coupled land physics is unchanged.

The regression compares the ReadyESM wrapper with the native atmospheric
timestep. It covers Float32 and Float64, absent and prescribed land, and
internal synchronization on and off. Each case takes four steps, including
startup and filtered steps. Both stored spectral levels and current and
previous grid fields match exactly.

| Test | Job | Result |
| --- | --- | --- |
| Focused native-equivalence tests | 60802879_0 | 426 assertions pass |
| Full CPU suite | 60802879_1 | Pass, with one existing expected-broken truncation test |

Both jobs tested source `85a854adb2f19af25ef85f6f7b95a61639c9bf5f`.
Output receipts and the frozen source checks pass. The full suite took
26 minutes 7 seconds. Production configuration and physics remain unchanged.
