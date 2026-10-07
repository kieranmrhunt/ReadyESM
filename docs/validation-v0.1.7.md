# v0.1.7 validation

Two-hour convection is the only physical change from v0.1.6. Reference RH
remains 0.7, liquid/ice cloud paths 75/25 g/m² and land precision Float64.
The generic timescale default remains four hours; production selects two.

## Coupled comparison

Jobs 57355307_0/1 start from the same native day-30 checkpoint and reach day120.
Both pass exact restored-state comparison, native diagnostic checks and the
additional 1e-6 mm land-water residual bound. Maximum land residuals are
1.112e-9 and 1.101e-9 mm. No tolerances were relaxed.

Time-weighted days90–120 means:

| Quantity | Four hours | Two hours |
| --- | ---: | ---: |
| Net TOA uptake, W/m² | 9.33946 | 7.27024 |
| Rain, mm/day | 2.11967 | 2.33982 |
| Snow, mm/day water equivalent | 0.14104 | 0.09556 |
| Surface-air temperature, K | 288.09539 | 287.96684 |
| Column water vapour, kg/m² | 28.22759 | 27.50861 |

The net-flux reduction is 22.16%. An earlier four-hour continuation gives
9.68699 W/m², so the observed control-repeat difference is smaller than the
parameter contrast. One candidate and two controls are not an ensemble.

## Limitations

The atmospheric residual inherits -0.19139 kg/m² at day30. At day120 it is
-0.41507 (control) and +0.08170 (candidate). Branch-period residual increments
are therefore -0.22368 and +0.27309 kg/m²: the smaller cumulative endpoint
partly cancels inherited error and is not a conservation fix.

Day120 maximum ocean temperature is 36.725 versus 36.889°C; maximum ice
area-equivalent thickness is 18.282 versus 18.382 m. Neither problem is fixed.
The small late-window ice difference is below observed control-repeat variation.
Passing native atmospheric checks establishes finiteness, not tight closure.

The tested branch switches convection at day30; production uses two hours
from initialization. No claim of fresh 120-day two-hour, equilibrium or
century-scale qualification is made. Native Float64 checkpoints are retained;
old Float32 land checkpoints are not silently converted.

## Evidence

The complete CPU regression suite passes in job57573912 (7 October,
16m33s), including 21 new convection-timescale assertions and the native
land-precision checks. CPU log SHA256:
`83a5e16a99a8da46d05acac51c7c09d32fb4a0db60a27a65dbc3c610f03f007c`.
Figure job57574158 passes 11 sample-contract assertions; both images have
been visually checked.

Raw review SHA256:
`e5af1b66a69373b9e34504cbc4edd76210c616497b9de78894edd628b705177a`.
Candidate source manifest SHA256:
`85d3a2f7b24d82d826ac1cfb3aff1399ef677ad8cf83eca933a7fcf5fae9edd1`.
The review checks 23 metrics, hashes, units/clocks and exact inherited day0–30
history. [Figure provenance](readiesm-v0.1.7-snapshot.json) records the output
and image hashes. Grey marks inactive cells; zero sea ice is white.
