# v0.1.17 tests

The standard atmosphere uses a four-hour divergence-damping time. The
previous value was one hour. Other production physics settings stay the
same. A checkpoint records the damping setting and refuses a changed or
missing setting on restart.

## Initialised runs

Three paired 31-day runs in January and July 1993 compare v0.1.17
with v0.1.16. Each pair starts from the same state. The setup uses a
T31/L27 atmosphere and a 1° tripolar ocean.

| Three-member mean | January, v0.1.16 → v0.1.17 | July, v0.1.16 → v0.1.17 |
| --- | ---: | ---: |
| Original 42-cell normalised error | 0.82418 → 0.81848 | 1.00332 → 0.99175 |
| Better than N of 60 CMIP6 models | 58.67 → 59.00 | 33.00 → 36.00 |
| Radiation-inclusive normalised error | 0.88447 → 0.87559 | 1.03048 → 1.01896 |

The original release score averages seven fields across 42 regional
cells. Its CMIP6 position uses six fields. The frozen gate classifies
January as neutral and July as improved. It accepts the two-season
result. The separate radiation-inclusive score gives the same season
verdicts. Its outgoing-longwave and albedo checks pass. The two-season
absolute net top-of-atmosphere bias rises from 0.456 to 0.586 W/m²,
within the 0.5 W/m² change allowance. Paired sea-surface temperature
and sea-ice-area safety checks pass in all six runs.

The January mean is neutral despite member 3's original normalised error
rising from 0.7999 to 0.8309. Conservation diagnostics are similar across
the matched runs, except the July member-1 atmospheric water residual
rises from 0.0094 to 0.0161 kg/m². This is still a result from short runs
that start from observations. The local ice pile-up and the ocean
freshening in longer runs need separate work.

The separate upper-level storm-strength screen missed its predeclared
improvement target. Mean absolute log EKE250 error fell by 0.034 in
January and 0.059 in July. The target was a fall of at least 0.10 in
both seasons. This is a missed improvement threshold, with no observed
EKE regression in that screen.
The [saved screen decision](validation-v0.1.17-eddy-screen.json) records
the original rule and values for all tested settings. It came from scorer
`62473933`, whose result-manifest SHA-256 is
`c1451804a1ee8488567b87f17b8c48b979834e39502293d4f5f84d8538806ba5`.
The screen used one matched member per season. The release comparison
above uses three.

The [README maps](readiesm-v0.1.17-maps.png) and
[timelines](readiesm-v0.1.17-timelines.png) use January member 1.
The [figure record](readiesm-v0.1.17-figures.json) gives the exact
diagnostic hash and plotting method. The run used the scored source
`a2f2ac4` with the four-hour setting.
The [README comparison](readiesm-v0.1.17-comparison.png) shows all six
paired original-score changes and the separate one-member-per-season eddy
screen. Its [numeric record](readiesm-v0.1.17-comparison.json) pins the
six scorecard hashes, the score job's result manifest and the saved eddy
decision. The [plot recipe](plot-v0.1.17-comparison.py) renders that record.

## Software checks

- The exact-source full CPU suite passed in job `62784235`.
- Coupled CPU restart tests `62786547` and `62786626` passed the
  changed-setting refusal and same-setting restore.
- Six GPU climate runs `62719568_[0-5]` and scorer `62719569`
  completed with verified result seals. Independent original-score
  checker `check_v3_d4_v2.py` returned `ACCEPT`.
- Fresh-process GPU source and restart test `62797763` passed with verdict
  `D4_RELEASE_GPU_PASS` and verifier result `D4_RELEASE_GPU_VERIFIED`.
  Its result-manifest SHA-256 is
  `635d97a9527439585ad6695d296d24cec6945ab4e64b5073729ad06f333eea9e`.

Scorer result-manifest SHA-256:
`068baa8b7cce44599bc6e38cf068be393ffb7dd0b6337f02ca01c4c87dee3e53`.

The native restart guard refuses a checkpoint whose damping setting
differs from the receiver. It also refuses a legacy checkpoint without
the new damping identity record. Such checkpoints cannot be restarted
directly with v0.1.17. Continue them with an earlier release, or start
a new run.
