# ReadyESM

ReadyESM couples SpeedyWeather, RRTMGP, Oceananigans, ClimaSeaIce and Terrarium.
The current model is a T31/L27 atmosphere over a 1° tripolar, 60-level ocean,
with dynamic sea ice and 16-layer land. ERA5 and ECCO provide the initial state.

This is a research model under development, not a calibrated projection model.
The production configuration is [`config/production.yml`](config/production.yml).

## Run

ReadyESM requires Julia 1.12 and an NVIDIA GPU with working CUDA drivers.

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
python scripts/download_era5_initial_state.py --include-stratosphere
julia --project=. scripts/run_dynamic_esm.jl config/production.yml
```

The ERA5 download requires CDS credentials. ECCO V4r4 data are loaded through
the ClimaOcean/NumericalEarth data path. `scripts/prepare_ecco_monthly_bundle.jl`
can verify and arrange a local monthly ECCO bundle.

The runner validates its diagnostics and writes results under
`artifacts/production`. To make a final-day atlas:

```bash
julia --project=. scripts/plot_final_day_summary.jl \
    artifacts/production/diagnostics.nc artifacts/production/final_day.png
```

CO2 is an active prescribed concentration in RRTMGP. Sulfate aerosol can be
applied as prescribed optical depth. There is no interactive carbon cycle,
emissions-driven CO2 or aerosol chemistry yet.

See [`NOTES.md`](NOTES.md) for the present scientific limitations and
[`HISTORY.md`](HISTORY.md) for the short model history.
