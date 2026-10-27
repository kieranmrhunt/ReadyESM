# v0.1.16 tests

The standard configuration now includes seasonal changes in Earth–Sun
distance. Other production settings and model equations stay unchanged.
Use the same solar-distance setting when restoring a checkpoint.

## Initialised runs

Three paired 31-day runs in January and July 1993 compare the change with
v0.1.15. Each pair shares its initial state. Both use T31/L27 and the
1° tripolar ocean.

| Three-member mean | January, control → candidate | July, control → candidate |
| --- | ---: | ---: |
| Normalised error, lower is better | 0.89050 → 0.88447 | 1.03108 → 1.03048 |
| Better than N of 59 CMIP6 models | 53.67 → 54.67 | 31.00 → 31.67 |
| Monthly ice-area bias, million km² | +1.3145 → +1.1723 | +2.4356 → +2.4366 |

The v4 error score averages nine fields after regional normalisation and
area weighting. CMIP6 position uses eight fields, excluding column water.
These values describe the initialised months. January improves beyond the
repeat spread. July stays within it, with neither score worse.

The two-season radiation, sea-surface temperature and ice-area checks pass.
July's outgoing-longwave bias increases by 0.82 W/m², while January's falls
by 0.71 W/m². Wind errors rise slightly in both months. Local ice accumulation
and longer-run drift remain under investigation.

## Software checks

- Targeted CPU tests pass in job 62168797.
- The full CPU suite passes in job 62199749. Its private-input ERA5 spectral
  branch was skipped because that checkout had no private data files.
  Separate data-enabled job 62388615 passes the complete unchanged module:
  97 assertions, including all seven ERA5 assertions, with no skips.
- Three GPU smoke runs pass in job 62172806.
- Both opposite-solar-setting rejection checks pass in job 62334655.
- GPU job 62383233 passes exact full-checkpoint comparisons at the restored
  boundary and final step. A different-time comparison fails as expected.

The first GPU attempt stopped at a script-scope error in its rejection test.
The follow-up passed both rejection checks, then found that the test's
uninterrupted run saved an extra callback. The resumed run lacked it.
CPU diagnostic 62382468 isolated the callback difference. The corrected GPU
test registers the same callback in both runs and passes the full comparison.
Each comparison group uses one allocated A100. The three smoke runs, two
rejection runs and new restart pair retain their separate allocation records.

GPU qualification job: 62383233. Result SHA-256:
`0d5d6749844aac17d495d8d57750405165a45558c1a38efe2f91ccd6c979139c`.

Scoring job: 62159035. Result SHA-256:
`3626764aaf31af8ffe27946d8af9e1236152d6995ccd489c51021bfd07f355ba`.

Qualified solar source: `cc588b2958acd5caf39f8d53beef800889e6cf2e`.
CPU-suite source: `0792497ec407d381ae2a43cf4d692e7e3516b031`, which adds
test-fixture repairs without changing the model.

The README figure uses the January run on source `a2f2ac4` with solar-distance
correction enabled. This matches the v0.1.16 model and configuration. Its
figure metadata retains the original run and render records, including the
candidate status at the time of rendering.
