# Historical radiation inputs for the 1993–1994 production run

`production.yml` starts on 1993-01-01 and requests 730 days. These two files
make the v0.1.21 gas and ozone settings independent of a developer's absolute
artifact paths.

| File | SHA-256 | Use |
| --- | --- | --- |
| `wmghg_cmip6_ssp245_1850-2100.nc` | `d399d2d530ba946f4a53ec01a25caba956a8e2ac21646f5e02c2f45ef7f1b1b9` | CMIP6 annual CO₂, CH₄, N₂O, CFC-12 and CFC-11-equivalent; production holds 1993 values. |
| `o3_cmip6_historical_199201-199512.nc` | `fc687e57baba2df2e83182453827f08e1543131b5296602c298c72998b021d15` | CMIP6 historical monthly zonal ozone for 1993–1994, with adjacent years for exact interpolation at both boundaries. |

The gas file is byte-identical to the prepared `ssp_forcing_v1` SSP245 file.
Its NetCDF attributes retain the UoM input4MIPs source identifiers and source
checksums. The ozone file is a 48-month extract of the prepared 1850–2099
SSP245 ozone file, SHA-256
`873d52845442e8e5b833d9fbd29281f3f81337b54e8adc811c8c21fabcd66945`.
Its attributes retain the UReading CCMI historical source identifiers, source
checksums, zonal reduction method and calendar details. Every selected NetCDF
variable was compared with its source. The model's 15-minute interpolation
brackets, source fields and Float64 weights match at all 70,081 instants from
1993-01-01 through 1995-01-01, inclusive. The portable [subset builder](prepare_1993_ozone.py),
[interpolation checker](verify_interpolation.py) and
[receipt](ozone_subset_receipt.json) are tracked here. The receipt preserves
the source and subset content hashes while replacing private workspace paths
with portable names. With the prepared full ozone source, regenerate the
1992–1995 subset using `--first-year 1992 --last-year 1995`, then compare
its SHA-256 with the packaged file. The checker accepts the full source,
packaged subset and this receipt; use `--first-model-year 1993
--last-model-year 1994` to check all 70,081 production timesteps. The full
source file remains outside Git.

This compact ozone file covers the released 1993–1994 configuration. Runs for
other years require an ozone file with appropriate monthly coverage. The gas
schedule's 1993 fixed-year setting is intentional for the matched development
screen; a transient historical run should instead set `well_mixed_gas_fixed_year`
to zero. A coupled model and CPU/GPU comparison remain required before this
candidate is promoted to a release.
