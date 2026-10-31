# Realistic seasonal land albedo candidate

`era5_monthly_land_cycle_1993_realistic_v2.nc` is the albedo input selected by
the v0.1.19 production profile for qualification. It retains the twelve 1993
monthly ERA5 cover, leaf-area and land-mask fields from the source land-cycle
file. Its four snow-free albedo bands have been adjusted toward the 1994–2013
ERA5 forecast-albedo (`fal`) monthly climatology, including the effects of snow
and land ice in that target. The ocean and sea ice do not use these land bands.

For each land cell and month, `make_realistic.py` adds one offset to the four
bands and clips each to [0, 1]. It chooses the offset so their **arithmetic
mean** equals the area-weighted `fal` target. The model's realised broadband
albedo is weighted by its radiation state and need not equal that target.
This construction does not recover snow's true spectral split and uses an
all-sky diffuse target for both direct and diffuse bands. See `ERRATUM.md`
for the correction to the unchanged derivation script's docstring.

`derivation_receipt.json` records the source land-cycle, held-out `fal`, script
and output hashes and summary checks. The original land-cycle source is
packaged in `../era5_seasonal_land_cycle_1993/`; the `fal` source is an external
ERA5 input and is identified by SHA-256
`65db76e869a0605dede1dd1ee6a2a5ac64cea38a759f63c8a9845eca6167f65a`.
The output and the files copied from the derivation stage are pinned in
`SHA256SUMS`. The source NetCDF's historical attribute calls this a test-arm
input; selection in the v0.1.19 production profile does not establish a
validated climate correction. The combined 90-day qualification remains to
be read before release.
