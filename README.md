# ReadyESM

ReadyESM couples SpeedyWeather, RRTMGP, Oceananigans, ClimaSeaIce and Terrarium.
The current model has a T31/L27 atmosphere and a 1° tripolar ocean with 60
levels, together with dynamic sea ice and a 16-layer land model. ERA5 and ECCO
provide its initial state.

ReadyESM is a research model under active development.
The production configuration is [`config/production.yml`](config/production.yml).

`v0.1.9` removes an invalid, frozen spectral row from the ERA5 launch state.
In a matched fresh 30-day pair, this change reduces the magnitude of
atmospheric-water drift by 94.8% and removes a spurious 80–90 hPa polar pressure
depression. The change has a material climate effect. Old checkpoints keep
their original launch state until explicit migration. The model retains
two-hour convection, mixed-phase saturation, double-precision land and liquid
and ice cloud paths of 75 and 25 g/m².
See the [release history](HISTORY.md) and
[validation record](docs/validation-v0.1.9.md).

The longest accepted trajectory remains the pre-fix v0.1.8 run to day 390.
The figure below shows this run. Validation of the repaired v0.1.9 launch
currently reaches day 30. Current work targets ocean hot cells, coastal ice
pile-up and whole-system energy closure.

![ReadyESM D: ocean, sea ice, atmosphere and land at model day 390](docs/readiesm-v0.1.8-snapshot.png)

*A 390-day initialised run. It uses v0.1.8 physics from day 90.*
[Timelines](docs/readiesm-v0.1.8-timelines.png) ·
[Figure provenance](docs/readiesm-v0.1.8-snapshot.json)

## Experimental roadmap

Snapshot: **15 October 2026**. Current blue-sky work:

| Track | Direction | Progress | Next step |
| --- | --- | --- | --- |
| BS01 | Conservative convection, then a learned closure | An idealised shallow-convection latitude taper recovers much of the missing Southern Hemisphere eddy activity. Its GPU lifecycle test passes, and a full day-30-to-120 comparison is running. | Measure circulation, water and radiation together in the coupled run. |
| BS02 | Prognostic clouds | The first local liquid and ice store is stable, conservative and restartable, but its coupled run has too little tropical and mid-level cloud and worsens net TOA by 19.8 W/m². | Test the relative-humidity threshold and conservative convective detrainment separately in short paired runs. |
| BS03 | Virtual vertical layers | Bottom-layer refinement improves an offline surface-longwave column, but the first coupled option worsens TOA balance. | Use conservative heating and a clean-launch control for the next test. |
| BS04 | Balanced initialisation | A targeted spectral launch cleanup removes a frozen pressure, temperature and humidity row. Version 0.1.9 uses this cleanup. | Extend the clean launch beyond 30 days. Test gradual analysis increments separately. |
| BS05 | Parameter calibration | Two 33-member real-column screens and the nine-arm convection screen are complete. Calibration remains at the screening stage. | Use the heating-budget diagnosis to choose identifiable parameters before a bounded coupled sweep. |
| BS06 | Sea-ice thickness distribution | Seven-day attribution shows that the Rothrock strength law caused almost all of the prototype's added pile-up. Retaining native aggregate pressure closes 94% of that excess. A repeated day-30 check is running. | Compare the original target, its 3-by-3 neighbourhood and relocated maxima with native ice. |
| BS07 | Coastal connectivity and landfast ice | Fine-bathymetry connectivity and metric-interface benchmarks pass around the problem coast. | Translate resolved channel geometry into a local dynamics benchmark. |
| BS08 | Runoff storage and heat | Standalone routing closes water and heat budgets and restarts exactly. A receiver heat-replacement GPU test is submitted. | Connect the land heat debit and source temperature without double-counting ocean heat. |
| BS09 | Automated experiment branching | Bounded planning and archived-failure replay tested. The live evidence adapter is under review. | Drive one reviewed shadow experiment with verified inputs and no duplicate submissions. |
| BS10 | Carbon, aerosols and emissions | Fixed-CO₂ radiation identity and closed carbon-box tests pass. Dated gases give the expected lower forcing. Production land uses prescribed evapotranspiration. | Add a shared CO₂ driver for interactive vegetation, then replay spatial air–sea exchange with a closed ledger. |

## Run

The production configuration requires Julia 1.12 and an NVIDIA GPU with working
CUDA drivers.

```bash
julia --project=. -e 'using Pkg' -e 'Pkg.instantiate()'
python -m pip install -r scripts/requirements.txt
python scripts/download_era5_initial_state.py --include-stratosphere
julia --project=. scripts/run_dynamic_esm.jl config/production.yml
```

The ERA5 download requires CDS credentials. The January 1993 ECCO V4r4 files
download through ClimaOcean and NumericalEarth. A public
[NumericalEarthArtifacts mirror](https://github.com/NumericalEarth/NumericalEarthArtifacts/releases/tag/data-v1)
is available when the primary ECCO download fails. Initial downloads require
network access. Later runs reuse the cached inputs.
`scripts/prepare_ecco_monthly_bundle.jl`
can verify and arrange a local monthly bundle for an explicit
`ecco_initial_conditions_directory`.

The runner validates its diagnostics and writes results under
`artifacts/production`. To make a final-day atlas:

```bash
julia --project=. scripts/plot_final_day_summary.jl \
    artifacts/production/diagnostics.nc artifacts/production/final_day.png
```

CO2 is an active prescribed concentration in RRTMGP. Production land uses
prescribed evapotranspiration without a CO2 response. Terrarium's ambient-CO2
input stays inactive in this configuration. A shared driver will connect CO2
to interactive vegetation. Sulfate aerosol can be applied as prescribed
optical depth. Interactive carbon, emissions-driven CO2 and aerosol chemistry
are planned extensions.

See [`NOTES.md`](NOTES.md) for the present scientific limitations and
[`HISTORY.md`](HISTORY.md) for the short model history.

## Test

The CPU regression suite uses synthetic initial conditions and requires no
ERA5 credentials or GPU. After instantiating the pinned environment, run:

```bash
julia --project=. --startup-file=no test/runtests.jl
```

The Fourier graph cache and conservative regridding regressions need a CUDA
GPU but no input downloads:

```bash
julia --project=. --startup-file=no test/cuda_graph_cache.jl
julia --project=. --startup-file=no test/csr_regridding_gpu.jl
```

On a CUDA machine with the production inputs, the short coupled and restart
checks run in separate Julia processes:

```bash
julia --project=. scripts/validate_coupled_smoke.jl artifacts/release-smoke
julia --project=. scripts/validate_coupled_restart.jl \
    artifacts/release-smoke/hour1_restart.jld2 artifacts/release-restart
julia --project=. scripts/compare_dynamic_restart_boundaries.jl \
    artifacts/release-restart artifacts/release-smoke/hour1_restart.jld2 \
    artifacts/release-restart/restored_boundary.jld2
```

These checks exercise one simulated hour followed by an independently restored
hour at production resolution. They verify exact restoration of checkpointed
model state and finite continuation. Longer experiments provide the climate
tests. Model construction and the first GPU step can take hours despite the
short simulated duration.

Developed by Kieran Hunt with AI coding assistance from Claude and OpenAI Codex.
