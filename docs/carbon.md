# Carbon

The standard configuration keeps carbon off and prescribes CO2 for radiation.
The optional carbon model adds a well-mixed atmospheric reservoir, ocean
dissolved inorganic carbon (DIC), alkalinity and air–sea CO2 exchange.

Copy `config/production.yml` and add:

```yaml
carbon:
  mode: emissions_driven
  atmospheric_co2_ppm: 420.0
  ocean_equilibrium_co2_ppm: 280.0
  emissions_pgc_per_year: 10.0
```

Keep the existing `forcing.co2_ppm` equal to `carbon.atmospheric_co2_ppm`
at launch. Run the copied configuration with `scripts/run_dynamic_esm.jl`.
Oceananigans and RRTMGP must be active.

| Mode | Behaviour |
| --- | --- |
| `off` | Standard prescribed forcing, with no carbon tracers or output. |
| `tracer_only` | Passive DIC and alkalinity transport on CPU. |
| `prescribed` | Air–sea exchange with fixed atmospheric CO2. Implied emissions balance ocean uptake. |
| `emissions_driven` | Constant emissions and ocean uptake change atmospheric CO2 and radiation. |

Set emissions to zero for modes other than `emissions_driven`. Emissions
are nonnegative, in petagrams of carbon per 365.25-day year. The ocean starts
from carbonate equilibrium at `ocean_equilibrium_co2_ppm`. Default alkalinity
is 2300 µmol/kg at salinity 35. Land vegetation remains prescribed.

## Output

`diagnostics.nc` contains final physical DIC and alkalinity maps. Numerical
storage offsets have been removed. Dry cells are NaN. File attributes record
the storage convention, offsets, reference density and final tracer inventory.

Exchange modes also write these histories on `carbon_time`:

- Atmospheric CO2 and atmospheric and ocean carbon ledgers.
- Cumulative emissions, ocean uptake and implied emissions.
- Total ledger residual and uptake rates.

The sample schedule matches the coupled budgets, including the final
endpoint. The interval-mean uptake rate integrates to the change in cumulative
uptake. The next-step rate describes the boundary supplied to the next coupling
step. Positive uptake moves carbon into the ocean.

Carbon histories are included in restart files. Use checkpoints from the
same carbon-output version. Output reports the final difference between the
tracer inventory and ocean ledger, using live cell volumes and the model's
reference density.
