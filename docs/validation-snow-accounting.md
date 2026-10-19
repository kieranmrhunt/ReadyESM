# Convective snow accounting

Accepted on 19 October 2026. Model commit: `92daafd`.

The update adds `atmosphere_convective_snowfall_flux` and
`global_convective_snowfall_flux`. It checks liquid and frozen precipitation
separately and restores the added history from older checkpoint fields.
Dynamics, forcing and surface fluxes stay unchanged.

The comparison uses public main `bb2eb09` and the candidate, each for eight
CPU steps at T31/L27 with the 360 × 180 × 60 ocean grid, dynamic sea ice,
Terrarium and RRTMGP. Both start from the same ERA5/ECCO configuration.

| Check | Result |
| --- | --- |
| Existing model-state fields | All 1,643 exact at initialization and the final step |
| Existing output variables | All 208 exact, including masks, types and units |
| New native-convection snow output | Zero, with the expected missing initial history value |
| Older diagnostic-history reconstruction | Exact |
| Full CPU regression suite | Pass, with one pre-existing broken spectral-truncation test |

The comparison jobs are 60701556, with reader 60701557. Their file hashes and
independent state comparisons were checked before publication. The existing
v0.1.10 tag stays unchanged.
