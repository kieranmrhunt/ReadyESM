# v0.1.9 validation

v0.1.9 fixes the ERA5 atmosphere launch. ReadyESM previously populated
SpeedyWeather's extra `l = lmax + 1` spectral row in logarithmic surface
pressure, temperature and humidity. Tendencies are truncated, so those
coefficients remained frozen rather than evolving with the atmosphere.

Production now removes that row before SpeedyWeather constructs the first grid
state. Coefficients at `l <= lmax` are unchanged bitwise. The generic option is
off for compatibility. Existing checkpoints are not changed automatically;
the explicit restored-state helper changes the trajectory and preserves the
previous-step grids by default.

## Fresh coupled comparison

Matched production-resolution runs start from ERA5/ECCO on 1 January 1993 and
differ only in launch truncation. Each reaches day 30 and passes its built-in
diagnostic validation and a source-compatible production safety audit.

| Quantity | Unfixed launch | v0.1.9 launch |
| --- | ---: | ---: |
| Atmospheric-water residual rate, days 5–30 (kg/m²/day) | +0.01560 | -0.00081 |
| Late net TOA uptake (W/m²) | 10.77 | 11.64 |
| Final atmospheric-water residual (kg/m²) | +0.475 | -0.0156 |
| Maximum area-equivalent sea-ice thickness (m) | 8.36 | 4.90 |
| Maximum ocean temperature (°C) | 32.62 | 32.39 |
| Minimum ocean salinity (psu) | 23.30 | 23.24 |
| Maximum land-water residual (kg/m²) | 2.79e-10 | 2.82e-10 |

The magnitude of water drift falls by 94.8%. The one-time launch difference is
-0.134 kg/m² total atmospheric mass, -0.00097 kg/m² water and -40.9 kJ/m²
`cpT`. No tropical ice occurs in either run. The sea-ice difference is
encouraging but one 30-day pair does not establish that coastal pile-up is
fixed.

The unfixed row also creates a persistent polar pressure error. At 87.2°N the
day-30 ring mean is 927 hPa with the imprint, 1017 hPa after removal and about
1007 hPa in ERA5. Removing it exposes an existing cold polar bias, including
roughly 15–18 K cooling in the upper-tropospheric comparison at that ring. This
is a material climate change caused by correcting the state, not parameter
tuning.

An explicit cleanup of a v0.1.8 day-120 checkpoint reduces its days 121–150
water-residual rate from +0.02317 to -0.000049 kg/m²/day. That experiment
demonstrates causality but is not a release trajectory or a recommended silent
checkpoint migration.

Result record SHA256:
`bab8514b3625a6d49139f4e0428df38b871506eeae0c27afd609ba9255983eb1`.
Physical-gate record SHA256:
`6ced66e37477550e463b2dca351ff9b8a61de6b917f3788dfa77a4b8abd94b8b`.

## Release qualification

The focused regression checks exact row removal, idempotence, retained
coefficients, first-grid reconstruction, previous-step grids, restored-state
semantics and finite stepping; all 97 assertions pass locally. Full CPU and production-resolution GPU
launch/restart qualification are recorded before the release tag is created.

This validation does not establish equilibrium, centennial stability or a
complete energy budget. Ocean hot cells and coastal ice pile-up remain open.
