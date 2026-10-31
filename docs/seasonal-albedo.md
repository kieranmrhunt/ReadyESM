# Seasonal ERA5 surface albedo

The v0.1.19 production profile selects `terrarium_surface_albedo: era5_spectral`.
The generic experiment default `constant` retains the existing 0.4 land
albedo. The option requires `land_model: terrarium` and
`ocean_model: oceananigans`. The production input is the pinned realistic
candidate under `data/era5_realistic_albedo_1993/`; the original ERA5
snow-free bands remain under `data/era5_seasonal_land_cycle_1993/`.

Configure the source with its recorded SHA-256 and choose how the RRTMGP band
that crosses the visible/near-infrared boundary is assigned:

```yaml
terrarium_surface_albedo: era5_spectral
terrarium_albedo_cycle_path: data/era5_realistic_albedo_1993/era5_monthly_land_cycle_1993_realistic_v2.nc
terrarium_albedo_cycle_sha256: 2f8e9193e2fea340e81dacc0e486c1d7de7f3b3a1b1e44dfc16ae3eef76f5c0c
terrarium_albedo_crossing_band: equal_halves
```

The source's visible and near-infrared direct and diffuse albedos remain
separate. They are mapped onto native RRTMGP shortwave bands. The crossing-band
rule is `equal_halves`, `visible`, or `near_infrared`; it is part of the restart
identity. After each radiation solve, ReadyESM diagnoses the scalar land
albedo that reproduces the mixed-cell reflected shortwave flux and holds that
same value for Terrarium until the next radiation solve. The source, mapping,
calendar and convention are checked when restoring a restart.

This is a production boundary choice under qualification, not a validated
climate correction. The source adjusts four snow-free ERA5 bands toward a
held-out ERA5 forecast-albedo climatology; the [data provenance](../data/era5_realistic_albedo_1993/README.md)
explains its spectral and direct/diffuse limits. The combined 90-day run tests
the three selected options together and cannot attribute its response to one
of them. Separate paired arms are needed for attribution.

**Implementation status:** an earlier development version replaced the
seasonal land albedo with the old single value just before each radiation
calculation, so the option had no effect in coupled runs. The fix now refreshes
the RRTMGP state first and then packs the seasonal bands. The CPU and A100
smoke checks verify that the solver receives those bands and that reflected
shortwave responds to the albedo change; both pass. A one-member July month
sensitivity test found an approximately 2 K improvement in 30–60°N 500-hPa
temperature bias, but global net top-of-atmosphere flux shifted by about
+9.65 W m⁻². The test is short and has vertical and regional overshoots.
The one-member 90-day combined run has completed on development source
`d325fc3`. Its 30–60°N 500-hPa weeks-4–12 error against ERA5 moved from
−3.31 K in the matched v0.1.18-settings run to +0.05 K. Global net
top-of-atmosphere flux increased by 9.71 W m⁻² in July, 7.46 in August and
6.65 over the 28 September days in the file; the upper-air warm bias worsened.
This combined run does not isolate the albedo effect or establish
repeatability or long-run energy closure. The six paired release months
remain to be scored.
