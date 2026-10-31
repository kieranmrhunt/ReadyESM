# File-based gas and ozone forcing

File forcing is optional and applies only with RRTMGP radiation. Leave both
paths empty to use the standard constant gases and analytic ozone profile.
Each input file is checked when the configuration is loaded. Its content hash
and selected year are part of restart identity. Use the same file contents and
settings when resuming a run.

## Well-mixed greenhouse gases

Set `forcing.well_mixed_gas_path` to an annual file with schema
`readiesm_well_mixed_gases_v1`. It contains global annual mean mole fractions
for `co2`, `ch4`, `n2o`, `cfc12` and `cfc11eq`, each in `mol mol-1`. The file
must have at least two consecutive years, finite nonnegative values, and
positive CO2, CH4 and N2O values. Annual values are placed at the middle of
each labelled year. The model linearly interpolates between those midpoints
using its Gregorian clock and stops if a time-varying run would need values
outside the file.

The `cfc12` field is CFC-12. `cfc11eq` is the CFC-11 equivalent of the other
halocarbons and is passed to RRTMGP as `cfc11`. The file replaces the
`forcing.co2_ppm` value and the optional `forcing.ch4_ppb` and
`forcing.n2o_ppb` settings. Leave the CH4 and N2O settings unset when using
the file. It cannot be combined with an active carbon mode, which also sets
CO2.

For a time-slice or pre-industrial control, set `forcing.well_mixed_gas_fixed_year`
to a year included in the file. The selected annual values are then held for
the whole run. The default `0` follows the model date.

```yaml
forcing:
  radiation: rrtmgp_all_sky
  well_mixed_gas_path: data/forcing/well_mixed_gases.nc
  well_mixed_gas_fixed_year: 0
```

## Zonal-mean ozone

Set `forcing.ozone_path` to a monthly file with schema
`readiesm_zonal_ozone_v1`. The file contains `vmro3` in `mole mole-1`, with
latitude, pressure in Pa, and month coordinates and bounds. Latitude and
pressure coordinates must be strictly ordered and in range. The model
interpolates between monthly midpoints in time, then linearly in latitude and
in log pressure to each radiation column and layer. It holds the end month
between the file boundary and its first or last midpoint, and stops if the
run goes outside the file's time bounds. Outside the file's pressure range,
the nearest ozone level is held.

For a repeating seasonal cycle, set `forcing.ozone_fixed_year` to a complete
year in the file. The model repeats its twelve months and interpolates across
December to January. The default `0` follows the model date.

```yaml
forcing:
  radiation: rrtmgp_all_sky
  ozone_path: data/forcing/zonal_ozone.nc
  ozone_fixed_year: 0
```

The input preparation scripts and their source records are maintained in the
development archive. The model reads the prepared files. It does not fetch
forcing data from the internet.
