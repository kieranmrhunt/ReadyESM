# ReadyESM

ReadyESM couples SpeedyWeather, RRTMGP, Oceananigans, ClimaSeaIce and Terrarium.
The current model has a T31/L27 atmosphere and a 1° tripolar ocean with 60
levels. It also has dynamic sea ice and a 16-layer land model. ERA5 and ECCO
provide the initial state.

ReadyESM is under development. Do not use it for projections.
The production configuration is [`config/production.yml`](config/production.yml).

`v0.1.9` removes an invalid, frozen spectral row from the ERA5 launch state.
In a matched fresh 30-day pair, this change reduces the magnitude of
atmospheric-water drift by 94.8%. It also removes a spurious 80–90 hPa polar
pressure depression. The change has a material climate effect. ReadyESM does
not apply it silently to old checkpoints. The model retains two-hour
convection, mixed-phase saturation and double-precision land. It also retains
liquid and ice cloud paths of 75 and 25 g/m².
See the [release history](HISTORY.md) and
[validation record](docs/validation-v0.1.9.md).

The longest accepted trajectory remains the pre-fix v0.1.8 run to day 390.
The figure below shows this run. The repaired v0.1.9 launch has been tested for
30 days. It is not a centennial or equilibrium qualification. Ocean hot cells,
coastal ice pile-up and whole-system energy closure remain unresolved.

![ReadyESM D: ocean, sea ice, atmosphere and land at model day 390](docs/readiesm-v0.1.8-snapshot.png)

*A 390-day initialised run. It uses v0.1.8 physics from day 90.*
[Timelines](docs/readiesm-v0.1.8-timelines.png) ·
[Figure provenance](docs/readiesm-v0.1.8-snapshot.json)

## Experimental roadmap

Snapshot: **14 October 2026**. These ten blue-sky tracks are experimental, not
features of the released model unless stated otherwise. A successful component
test does not establish coupled climate skill. Releases do not have to wait for
all tracks.

| Track | Direction | Progress | Next step |
| --- | --- | --- | --- |
| BS01 | Conservative convection, then a learned closure | A nine-arm 60-day prescribed-SST screen found no material trigger, humidity or plume candidate. The plume arms weaken overturning, dry column water and run about five times slower. | Decompose latent and radiative heating in ascent and subtropical subsidence before another closure sweep. |
| BS02 | Prognostic clouds | The first local liquid and ice store is stable, conservative and restartable. Its coupled run has far too little tropical and mid-level cloud. It also worsens net TOA by 19.8 W/m². | Test the relative-humidity threshold and conservative convective detrainment separately in short paired runs. |
| BS03 | Virtual vertical layers | Bottom-layer refinement improves an offline surface-longwave column, but the first coupled option worsens TOA balance. | Revisit only with conservative heating and a clean-launch control. |
| BS04 | Balanced initialisation | A targeted spectral launch cleanup removes a frozen pressure, temperature and humidity row. Version 0.1.9 uses this cleanup. | Extend the clean launch beyond 30 days. Test gradual analysis increments separately. |
| BS05 | Parameter calibration | Two 33-member real-column screens and the nine-arm convection screen are complete. No fitted calibration has been accepted. | Use the heating-budget diagnosis to choose identifiable parameters before a bounded coupled sweep. |
| BS06 | Sea-ice thickness distribution | Five-category GPU evolution and bitwise fresh-process restart pass on a synthetic immersed island. A coupled five-arm Arctic factorial is queued. | Test material hotspot reduction and displacement. The pile-up is not fixed yet. |
| BS07 | Coastal connectivity and landfast ice | Fine-bathymetry connectivity and metric-interface benchmarks pass around the problem coast. No live coastal stress scheme yet. | Translate resolved channel geometry into a local dynamics benchmark. |
| BS08 | Runoff storage and heat | Standalone routing closes water and heat budgets and restarts exactly. A receiver heat-replacement GPU test is submitted. | Connect the land heat debit and source temperature without double-counting ocean heat. |
| BS09 | Automated experiment branching | Bounded planning and archived-failure replay tested. The live evidence adapter is under review. | Drive one reviewed shadow experiment with verified inputs and no duplicate submissions. |
| BS10 | Carbon, aerosols and emissions | Fixed-CO₂ radiation identity and closed carbon-box tests pass. Dated gases give the expected lower forcing. Production land currently has no CO₂-responsive vegetation. | Add a shared driver before enabling interactive vegetation. Then replay spatial air–sea exchange with a closed ledger. |

## Run

The production configuration requires Julia 1.12. It also requires an NVIDIA
GPU with working CUDA drivers.

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
input is dormant because production does not construct a photosynthesis model.
ReadyESM requires a shared driver before it can enable CO2-responsive
vegetation. Sulfate aerosol can be applied as prescribed optical depth. The
model does not yet have an interactive carbon cycle, emissions-driven CO2 or
aerosol chemistry.

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
model state and finite continuation. They do not establish a stationary climate
or exact equality between continuous and restarted trajectories. Model
construction and the first GPU step can take hours despite the short simulated
duration.

Developed by Kieran Hunt with AI coding assistance from Claude and OpenAI Codex.
