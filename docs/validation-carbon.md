# Carbon extension checks

Accepted as an optional model component on 20 October 2026.
Production keeps carbon off. See [configuration and output](carbon.md).

The extension adds ocean DIC and alkalinity, air–sea CO2 exchange and a
well-mixed atmospheric carbon reservoir. Constant emissions can change
atmospheric CO2 and its radiative forcing. Land carbon and aerosols remain
outside this implementation.

## Evidence

| Check | Result |
| --- | --- |
| Carbon physics, configuration and changing-volume inventories | 6,710 focused assertions and the full CPU suite pass. |
| Carbon disabled | All 431 common state entries match the previous model at construction, initialization and step two. |
| Four compact GPU emissions restarts | All 1,552 saved state entries match at the restart boundary and four subsequent steps. |
| Four full-grid GPU days | All 97 ledger boundaries pass in each run, with maximum relative inventory error below 5.7 × 10⁻⁹. |
| Output and history | 1,992 focused assertions and 64 fresh-process history checks pass. |
| Five compact CPU output cases | Off, passive tracers, prescribed CO2, zero emissions and 10 PgC/year pass. |
| Four compact GPU output cases | Off, prescribed CO2, zero emissions and 10 PgC/year pass. Passive tracer mode remains CPU-only. |
| Combined main-branch source | Full CPU suite passes. Carbon code and tests match the qualified source. |

The full-grid runs use T31/L27 with a 360×180×60 ocean. They cover fixed
420 ppm, zero emissions with ocean uptake, and 10 PgC/year with live or
zero ocean exchange.

Output checks compare saved physical states at construction and steps 1–8.
Fresh-process restarts reproduce steps 4–8 and their diagnostic histories.
The 193 shared NetCDF variables retain their values and metadata.
Exchange modes add 13 variables and 13 global attributes. File creation
times are checked separately. Collecting and writing output leaves the
model state unchanged.

## Provenance

Physics source: `68988e3`. Output source: `6094781`. Combined CPU-qualified
source: `75c35f9`. Current-main integration: `df834a8`.
The integration retains main's standalone atmosphere restart repair.
Its model source, tests, dependency pins and run scripts match `75c35f9`.

Slurm records: `60869035`, `60870291`, `60870410`, `60874069`,
`60874070`, `60925029`, `60928226`, `60942296`, `60945298` and `60956342`.
Independent audits recheck source mappings, receipts, raw saved states,
NetCDF values and attributes.

These runs cover conservation, coupling, output and restart behaviour.
The ocean starts from prescribed carbonate equilibrium. Long-term carbon
uptake and climate sensitivity remain to be evaluated.
