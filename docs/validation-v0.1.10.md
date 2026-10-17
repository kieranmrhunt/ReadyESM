# v0.1.10 validation

v0.1.10 keeps the v0.1.9 model physics and changes the production ocean launch
state. Ocean temperature, salinity, velocity and sea-surface height start from
the midpoint of the December 1992 and January 1993 ECCO V4r4 monthly means.
Sea-ice area and effective thickness keep the January state.

The preparation tools verify the source receipts, physical bounds and output
hashes. The production runner verifies the eight base files and five generated
inpainting caches before the first time step.

## Matched coupled comparison

Three 30-day calendar-centred and control pairs use the same model source and
physics. The calendar-centred arm uses the new ocean state. The control uses
January ECCO. Each seeded perturbation survives coupled initialization.

| Quantity | January control | Calendar-centred launch | Change |
| --- | ---: | ---: | ---: |
| Days 1–10 OISST RMSE (K) | 1.237 | 1.029 | -0.208 |
| Southern Ocean days 1–10 RMSE (K) | 2.043 | 1.568 | -0.475 |
| Antarctic area absolute error (10-day mean, 10¹² m²) | 0.104 | 0.079 | -24% |
| Antarctic 15% extent absolute error (10-day mean, 10¹² m²) | 0.548 | 0.630 | +15% |
| Antarctic IIEE (10-day mean, 10¹² m²) | 2.810 | 2.841 | +1.1% |

Each pair improves OISST RMSE by 0.207–0.210 K. The frozen native-control
repeat envelope is 0.0027 K. The mean gain on days 8–10 is 0.193 K, which
retains 95.3% of the day-zero gain. All pairs pass the preregistered Antarctic
regrowth guard on ice area and 15% extent. IIEE is report-only.

Calendar-centred ocean freshwater relative closure peaks at 0.00331. The
frozen stop is 0.01. Ocean temperature stays between -1.97 and 32.47 °C.
Salinity stays between 21.92 and 40.64 g/kg. The atmosphere-water residual
rate is between -0.00044 and -0.00038 mm/day across the calendar-centred
members.

The frozen score returns `PROCEED`. Its SHA256 is
`63a97bd02c6d97534f1cd252c14e4c7ffc2dbb1f836c51ae4d2d35c472987def`.

## Release qualification

The full CPU suite passes in Slurm job `60222743`. The real ERA5 spectral-launch
test passes 97 of 97 checks in job `60224166`, with no skipped or broken result.
Its log SHA256 is
`90cb438aad037aec27dc7ad3dca5b1761f98bb0da0d653236bf5feb2607040eb`.

The first A100 release gate, job `60217860`, reproduces the qualified ocean and
sea-ice construction and initialized states. It also completes the 96-step
production smoke and writes its diagnostics and summary figure.

Fresh-cache A100 job `60222791` starts from the eight base ocean files. It
rebuilds all five caches, verifies their release hashes, completes a 96-step
production run and reproduces the qualified calendar and monthly-control state
hashes. The qualified source commit is
`8a5545ecab3b29b55df1032febe931a22a504443`.

This release improves the initial ocean fit. Antarctic extent, coastal ice
pile-up, ocean hot cells, whole-system energy closure and long-run climate
drift remain active model-development targets.
