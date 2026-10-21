# Land surface flux averaging

The CPU and A100 fresh-restart checks pass.

Terrarium takes several land steps during each coupling interval. ReadyESM
previously exported the flux evaluated at the end of that interval. Land
storage had changed under the fluxes used during the individual steps.

The repair accumulates sensible and latent heat flux multiplied by each
land step's duration. It exports the interval mean. The scratch integrals
reset each interval. Production parameters and the atmospheric time filter
stay unchanged. Runoff heat transport stays off by default.

## January comparison

Two initialized 31-day runs use T31/L27 atmosphere, the 1° tripolar ocean
with 60 levels, and 16-layer Terrarium land. ERA5 scores use common spatial
masks and identical averaging periods.

| Measure | Previous exchange | Averaged exchange |
| --- | ---: | ---: |
| Cumulative regular-step water exchange mismatch, kg/m² | −0.0143931 | 2.635e−10 |
| Net downward TOA flux, W/m² | 9.695 | 9.401 |
| Column-water RMSE, kg/m² | 5.989 | 5.870 |
| 500 hPa temperature RMSE, K | 4.770 | 4.822 |
| Precipitation RMSE, mm/day | 2.492 | 2.495 |

The water ledger compares land-applied water with the exported flux.
Startup retains a separate −6.33e−5 kg/m² mismatch. The atmospheric water
budget and cold bias remain open. The small climate differences come from
one paired trajectory.

## Checks

The full CPU suite and Float32/Float64 land-adapter checks pass. Fresh CPU
restarts after coupled steps 1 and 3 reproduce all subsequent saved states
and fluxes through step 8. An independent comparison verifies 26 records.

The A100 full-grid run and fresh restarts after coupled steps 1 and 3 pass.
An independent comparison verifies the restored boundaries and all later
saved states and fluxes through step 8, for 26 exact record comparisons.
The release tree matches the tested candidate byte for byte across 333
model files. These restart checks cover eight coupled steps.

The monthly reader verifies common masks and saved grid weights. A separate
calculation reproduces 60 raw water, radiation and ocean/ice summary values.
The detailed run records remain in the development labbook.
