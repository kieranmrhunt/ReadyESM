# Sigma-RH convection profile

The v0.1.13 candidate changes the production Betts–Miller reference humidity.
It is 0.5 above sigma 0.3, rises linearly to 0.7 at sigma 0.7, and stays
0.7 below. The two-hour adjustment time and model code are unchanged.

In the matched 31-day January and July screens, the frozen v4 score was:

| Start month | v0.1.12 ONR9 | Candidate ONR9 | CMIP6 models worse, before → after |
| --- | ---: | ---: | ---: |
| January 1993 | 1.168365 | 1.132393 | 12 → 16 of 59 |
| July 1993 | 1.330989 | 1.311963 | 5 → 6 of 59 |

The absolute value of the two-month mean OLR bias fell from 17.24 to
12.09 W/m². The corresponding planetary-albedo bias fell from 0.0328 to
0.0253, and net downward TOA bias from 6.65 to 4.05 W/m².

The final production config passed a Julia parser check. The full CPU
suite on the final release base passed on 23 October 2026, with an
unchanged-file receipt. A compact T7/L8 idealized
fresh-process coupled restart reproduced all compared components exactly
at 0.125 day. This restart check does not cover the full production grid.

Three fresh control/candidate pairs in each month completed 31 days. The
frozen pair checker verified their initial states and input files. Every
run passed the production diagnostics check. The formal v4 gate accepted
the metric change in both months:

| 1993 start | Control ONR9 | Candidate ONR9 | CMIP6 models worse, control → candidate |
| --- | ---: | ---: | ---: |
| January | 1.16538 | 1.13424 | 13 → 16 of 59 |
| July | 1.33117 | 1.31016 | 5 → 6.33 of 59 |

All six pairs improve ONR9 and CMIP6 position. The two-season mean absolute
OLR bias is 17.21 → 12.12 W/m², planetary-albedo bias is 0.03282 →
0.02504, and net TOA bias is 6.61 → 4.17 W/m². The frozen radiation guards
pass.

The day-31 gross-breakage check passes for every paired ocean mean surface
temperature and sea-ice area. The monthly global sea-ice-area bias against
NSIDC remains high. It changes from +1.775 to +1.760 million km² in January
and from +3.626 to +3.617 million km² in July. Ice bias is a report-only
metric under the declared physical-safety addendum.

Conservation remains within the model's production diagnostics checks in
all 12 runs. Across the candidate members, the largest absolute final
atmospheric-water residual is 0.0140 kg/m² and the largest land-water
residual is 2.96×10⁻¹⁰ kg/m². The largest final relative ocean heat,
freshwater and salt closure values are 1.40×10⁻⁴, 4.38×10⁻⁴ and
6.12×10⁻⁷. The July heat and freshwater residuals rose in some pairs,
though their final values remain small. The release preflight requires the
saved conservation readouts and a passing production diagnostics report for
every control and candidate member.
