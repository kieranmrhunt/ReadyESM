# CMIP6 total solar irradiance values

`tsi_values.toml` holds the two values the opt-in `forcing.solar_constant_w_m2` setting can use with a citable source:

- the CMIP6 pre-industrial control value, for runs with constant forcing.
- the 1993 annual mean of the CMIP6 historical series, for runs that start from observed 1993 states. This is a constant
  value for the whole run, not the time-varying CMIP6 historical forcing.

The values are read from the CMIP6 SOLARIS-HEPPA v3.2 input4MIPs files in the JASMIN CEDA archive. The file records
their paths, SHA-256 sums, the variable and the reduction. Reference: Matthes, K., et al. (2017), Solar forcing for
CMIP6 (v3.2), Geoscientific Model Development, 10, 2247-2302.
