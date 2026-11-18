# v0.1.21 validation record

## Release scope and decision, 18 November 2026

The release source is commit `a6cbb2648c8e5c908d058e8358a4a2e95ae7a272`.
The January–June 1993 science run used its physics parent
`396474675f78358401fcea445c4b805e1432a4d0`. The two later commits
change restart identity and the monthly resume runner; they do not change
fresh-run radiation or climate settings. The standard profile uses the
packaged 1993 greenhouse-gas values, monthly historical ozone, a sampled
monthly upper radiative cap and a diagnostic-cloud minimum relative humidity
of 0.4. The earlier 420 ppm CO₂ setting is an inactive fallback while file
forcing is selected.

The declared radiation target improved against the same ERA5 monthly
reference. The calendar-day-weighted January–June global net top-of-atmosphere
bias changed from **+2.30** W m⁻² in v0.1.20 to **+0.55** W m⁻² in v0.1.21.
Clear-sky outgoing longwave bias changed from **−7.07** to **−3.17** W m⁻².
These are one-start, six-month radiation comparisons, not an equilibrium
energy balance. The v0.1.21 A100 run `53279` completed 181 days and 17,376
steps; source-bound terminal audit `56952` passed. Its `result.sha256` file
has SHA-256 `f0ba2b6a6257bfff746488371df2933162a0351ae43cb666e2c32d56fdcccfe0`.

The frozen Part C V3 joint score `2584091` compared v0.1.19, v0.1.20,
v0.1.21 and a no-cap diagnostic on a single common-valid mask with 55 CMIP6
models. Its `RESULT_SHA256SUMS` has SHA-256
`42ad50448892d02bde93352d5dce80550a97701495a34475b35df57742734bb4`;
the scored JSON has SHA-256
`69449672e24f2ae4fb4f23ce75d5baf0be6b147a2c28b08111f2ec5b32ba370a`.
The eight-field overall normalised error ratios are **1.2005** for v0.1.19,
**1.1603** for v0.1.20 and **1.1788** for v0.1.21. Lower is better.
The v0.1.21 increase of **0.0185** against v0.1.20 is a regression, within
the observed single-start month-choice spread, and was accepted as a
recorded release limit. The March–June monthly fields are scored equally;
January–February influence the trajectory but do not enter ONR8 directly.
The [scorecard](../scorecard.png), [figure record](scorecard-v0.1.21-partc.json)
and [plot recipe](plot-v0.1.21-scorecard.py) retain this method and the
regional values. The [version history](../readyesm_scorecard_versions.png)
keeps the older January/July 31-day method in a separate panel.

## Restart qualification

The January–February calendar restart pair has an exact physical comparison
on the qualified source. Comparison job `2606587` completed successfully.
The first terminal audit, `2606736`, failed because it treated
`diagnostics.nc/cmip_aerosol_validated_writes` as continuous model state.
That int32 attribute counts validated aerosol writes **within each process**:
1417 in the uninterrupted run and 672 in the resumed February process.
The isolated corrected audit `2937501` requires those exact values and
matches every other diagnostic value, mask, schema and attribute across 13
NetCDF pairs. It also passes seven exact physical endpoint groups and the
radiation restart state. Its audit JSON SHA-256 is
`f2e6500c51ed44915f16a152dfd7657404b5b7908d1e80edc2d0cb5b9a7ae36a`.
The original failed audit remains part of the record. This qualifies the
tested calendar restart, including scored outputs, under the narrow
process-local metadata exception; it does not establish arbitrary-date or
multi-year restart identity.

## Known physical limits

The northern 30–60°N 500-hPa temperature score is worse than in v0.1.20:
its RMSE to the CMIP6 median rises from about **1.04** to **1.56** on the
joint mask. The model's winter top is roughly **20–35 K too cold**.
Factor tests indicate that v0.1.20's warmer free troposphere partly depended
on incorrect 1993 forcing: fixed 420 ppm CO₂ and analytic ozone. Restoring
those inputs would hide this bias rather than correct its physics. A drier
upper convection reference is being tested separately for v0.1.22.

The common-valid mask excludes missing cells. At 850 hPa, the frozen scorer
uses ERA5 surface pressure above 875 hPa and shared finite data; it does
not screen each CMIP6 model against its own surface pressure. Some model
below-surface values may remain, which can affect 850-hPa errors and ranks.
One observed start and four scored months do not demonstrate long-term
climate skill, seasonal robustness, equilibrium or repeatability. Ocean,
land, snow and sea-ice physical limits described in [the notes](../NOTES.md)
remain open unless separately qualified.
