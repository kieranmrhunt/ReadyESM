# v0.1.15 tests

GPU and fresh-process restart tests passed in job 62105211 on 25 October 2026.

## Change

Add K-profile vertical momentum mixing in the atmospheric boundary layer.
Production uses scale 1. Heat and humidity mixing stay off. The atmosphere
retains full-temperature vertical transport and two-hour Betts–Miller
convection. Cloud and forcing settings are unchanged.

## Initialised-month comparison

Three paired 31-day runs in each of January and July 1993 compare the
candidate with the v0.1.14 configuration. Each pair has the same initial
state. Both use T31/L27 and the 1° tripolar ocean.

The v4 score includes radiation. Regional errors are divided by the CMIP6
median, weighted by region area, then averaged over nine fields. Lower
normalised error is better. CMIP6 position uses eight fields, omitting
column water, and counts models with a larger error among 59 coupled
historical models. Positions below are means over three members.

ERA5 and ECCO constrain ReadyESM's starting state. CMIP6 historical runs
branch from long control simulations and evolve freely, so their weather
and ocean variability differ from the observed 1993 sequence. The rank is
a reference for this short-run error measure. The paired ReadyESM changes
provide the release evidence. Climate skill needs longer runs and comparison
of multi-year means. See the [CMIP6 experiment design, section A2](https://gmd.copernicus.org/articles/9/1937/2016/).

| Measure | January control → candidate | July control → candidate |
| --- | ---: | ---: |
| Normalised error | 1.011 → 0.891 | 1.139 → 1.031 |
| CMIP6 position, out of 59 | 31.0 → 53.7 | 14.3 → 31.0 |
| 850-hPa zonal-wind bias, 20°S–20°N, m/s | −5.89 → −1.67 | −3.33 → −0.75 |
| Monthly global ice-area bias, million km² | 1.124 → 1.315 | 2.101 → 2.436 |

Both skill measures improve in all six pairs and exceed the repeat-noise
allowances. The two-season absolute radiation biases decrease: net flux
0.294 → 0.253 W/m², outgoing longwave 7.80 → 7.42 W/m², and albedo
0.02373 → 0.02276. Each radiation check first averages the signed January
and July biases for a member, then takes the absolute value, then averages
over the three members.

The day-31 tropical zonal-wind jump between the lowest two atmospheric
levels falls from 9.5–10.1 to 0.76–0.93 m/s in January, and from 4.2–5.2
to 0.48–0.52 m/s in July. This diagnostic covers 20°S–20°N.

## Physical review

All six surface-temperature checks pass. The January day-31 ice-area
check fails in members 1 and 2. The increases are 1.05 and 1.23 million
km², against allowances near 0.91 million km². January member 3 and all
July members pass.

A physical-review override supports an incremental release because the
atmospheric skill gains are large and repeatable, and radiation improves.
The automatic physical result remains failed. Monthly global ice-area
bias increases in both months, including all three July members.

January maximum area-equivalent ice thickness falls in all three members.
July maxima change from 5.61 to 6.52 m, 6.02 to 4.30 m, and 4.26 to
4.30 m. The original Ross Sea hotspot thins in all three. Local accumulation
remains under investigation.

A consistent heat return for momentum mixing remains to be implemented.
Its energy effect at production resolution and long-term drift remain
to be tested.

## Source and software checks

- Scientific branch: `f5572aab934f87cdea0d55199f04ae4a460d3e52`.
- Release port: `73aef734d2c381faf63934ac37eaef67a62e93e5`.
- Repeats: `62078754`, with member-1 runs from `62068262`.
- v4 scoring: `62078821`. Physical review: `62105048`.
- Port CPU tests: `62080639`, 2,206 assertions passed. Default-off and
  option-on CPU reference states match their tested sources.
- Port GPU and fresh-process restart: `62105211`, passed. With mixing off,
  all ten initial and final field hashes match v0.1.14. With mixing on,
  they match the scientific branch. Restored states match exactly at the
  step-4 boundary and step-8 endpoint in a fresh process. The comparison
  rejects states from different times. All runs used the same A100 GPU.
- GPU result-manifest SHA-256:
  `455eb1cea5929816d389ece7ee60647f0d4428afad015b3679cdaff7713efa3f`.

Production configuration matches the tested momentum-only settings.
Start a fresh run to use the new
production settings. Keep the original diffusion settings for existing
checkpoints.

The README image shows the saved final state of January member 2. Rainfall
is an instantaneous flux expressed in mm/day. Figure metadata records
the source run, data checksum and renderer checksum.
