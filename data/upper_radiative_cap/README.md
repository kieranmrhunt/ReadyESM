# Provisional ERA5 January 1993 upper-radiative-cap temperatures

`era5_19930101_zonal_5_1hpa.csv` is the longitude mean of ERA5 temperature
at each integer latitude, at 5 and 1 hPa on 1993-01-01 00Z. It is derived by
`prepare_era5_jan1993.py` from the local ERA5 stratosphere NetCDF with SHA-256
`fb3a51a140a665fddbba724079ef8b2fdf9587ebcd89c81027e53872fea354c7`.
The CSV SHA-256 is
`748910cd885d0d52f57f368373a444bb8aba54bbf05cd17505986dd2cbd2f9c3`.
The raw NetCDF remains outside Git because it is much larger.

This is a **single start-state analysis**, not an ERA5 monthly zonal mean or a
climatology. Its use is restricted by configuration validation to January 1993.
The two-layer cap prototype interpolates in log pressure between 5 and 1 hPa
and provisionally holds the 1 hPa temperature above 1 hPa. A general
production option needs a sourced monthly upper-atmosphere temperature field
and review of the extension above 1 hPa.

## Experimental sampled annual cycle

`era5_1993_sampled_monthly_upper_cap.csv` has twelve months and 181 integer
latitudes. It is derived from eight 1993 ERA5 model-level temperature
snapshots per month, using pressure-overlap mass means in the 500–100 Pa and
100–1.005183574463 Pa layers. It also contains 500-Pa, 100-Pa and near-1-Pa
boundary temperatures. The table SHA-256 is
`5f07536991bd4a6a59536f7e5d86b58af803813c043d154d7a7d0d9184f5ca28`.
The source-file list and 96 file hashes are in
`era5_1993_sampled_monthly_provenance.json`. The exact preparation code and
ECMWF hybrid coefficients are included here as
`prepare_era5_1993_sampled_monthly.py` and
`ecmwf_l137_half_levels_0_24.csv`. The independent audit is preserved in
`artifacts/development_wave_20261111/v021_era5_monthly_cap_codex_v1/`
and `artifacts/development_wave_20261111/v021_era5_monthly_cap_audit_codex_v1/`
in the development workspace.

The experimental adapter repeats this 1993 sampled annual cycle in any model
year and interpolates continuously between Gregorian monthly centres,
including across December and January. It is not a complete hourly monthly
mean, multi-year climatology or general production forcing. Its checkpoint
identity binds the table SHA and interpolation method.
