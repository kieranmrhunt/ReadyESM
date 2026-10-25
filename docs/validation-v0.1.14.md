# Full-temperature vertical transport

This change ports [SpeedyWeather PR #1286](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/1286),
commit `adf0d61`, to SpeedyWeather 0.21.1. Vertical motion now transports
the background temperature gradient as well as local temperature variations.
Production enables `atmosphere_vertical_advection: centered_full_temperature`.
Convection, clouds and the implicit solver retain their v0.1.13 settings.

## Climate comparison

Three fresh paired runs start in each of January and July 1993 and last
31 days. Each pair has matching initial states and input files. All twelve
runs pass production diagnostics. The agreed v4 comparison gives:

| Start month | v0.1.13 ONR9 | Candidate ONR9 | CMIP6 models worse, before → after |
| --- | ---: | ---: | ---: |
| January | 1.1342 | 1.0105 | 16 → 31 of 59 |
| July | 1.3102 | 1.1394 | 6.33 → 14.33 of 59 |

Values are three-member means. Each pair improves both scores. ONR9 is
regional error normalized by the CMIP6 median, averaged over nine fields.
CMIP6 position uses eight fields, omitting column water. Fractional positions
come from averaging the three integer ranks.

The absolute two-season mean OLR bias falls from 12.12 to 7.80 W/m²,
planetary-albedo bias from 0.02504 to 0.02373, and net downward TOA bias
from 4.17 to 0.29 W/m². All frozen radiation guards pass.

## Ocean and sea ice

Endpoint SST is 0.61–0.73 K cooler than the controls. All six pairs improve
global open-water bias and RMSE against ERA5. This holds against either
bracketing monthly mean or their midpoint. The references are monthly
fields, used to approximate conditions at the endpoint.
Against their midpoint, the three-member mean global open-water bias falls
from +0.86 to +0.08 K in January and from +0.98 to +0.33 K in July.

The monthly excess global ice area against NSIDC falls from 1.760 to
1.124 million km² in January and from 3.617 to 2.101 million km² in July.
Every pair improves. The changes exceed the existing review bounds for
SST in six pairs and ice area in five. The observational comparisons
support accepting these changes. Codex and Claude agreed to a documented
physical-review override, with the failed original bounds retained.

Low-level wind errors increase. Tropical 850-hPa zonal-wind bias changes
from −0.49…−0.42 to −5.95…−5.85 m/s in January and from −0.19…−0.16
to −3.53…−3.11 m/s in July. The normalized u850 field error rises from
0.90–0.91 to 1.52–1.63 and from 1.28–1.29 to 1.59–1.71, respectively.
January tropical SST has a small cold
overshoot. July ice becomes thicker locally near 72.4°S, 170.5°E, with
daily maximum area-equivalent thickness reaching 5.7–6.0 m across the
three candidates. Total July ice volume falls in every pair. Wind mixing
and coastal ice dynamics remain development priorities.

## Software checks and restart

The largest absolute final atmospheric-water residual across the six
candidates is 0.01433 kg/m², and the land-water residual is
2.94×10⁻¹⁰ kg/m². The largest final relative ocean heat, freshwater and
salt closure residuals are 4.93×10⁻⁵, 7.69×10⁻⁴ and 3.65×10⁻⁷.
The July freshwater relative residual reaches 0.0165 during the month.
The day-1 absolute residual also increases in that pair, so the relative
increase is partly due to the residual and partly to its changing denominator.

The release-source CPU suite passes 614 checks. Eight saved arrays from
each short legacy and corrected run match their reference runs exactly.
Production-grid fresh-process restarts on both the scored source and the
final release source reproduce the restored state and continued endpoint
exactly. Each test resumes four coupled steps in a separate process.
Their comparators also reject checkpoints from different model times.

To continue a v0.1.13 production checkpoint, set:

```yaml
atmosphere_vertical_advection: speedyweather_0_21_1_centered
```

The checkpoint and configuration must use the same setting. Start a fresh
run to use the correction. Configurations that omit the option retain the
legacy behaviour.
