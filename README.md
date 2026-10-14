# ReadyESM

ReadyESM couples SpeedyWeather, RRTMGP, Oceananigans, ClimaSeaIce and Terrarium.
The current model is a T31/L27 atmosphere over a 1° tripolar, 60-level ocean,
with dynamic sea ice and 16-layer land. Initial state comes from ERA5 and ECCO.

This is still under development so please don't use it for projections.
The production configuration is [`config/production.yml`](config/production.yml).

`v0.1.8` adds mixed liquid/ice saturation and net-column
large-scale rain diagnosis: surface rain matches that scheme's column water
removal, but total atmospheric water is still(!) not conserved. Two-hour convection,
double-precision land and cloud paths of 75/25 g/m² are retained.
See the [release history](HISTORY.md) and
[validation record](docs/validation-v0.1.8.md).

The current version reaches at least day 390. Compared with 0.17, over days 360–384, net
top-of-atmosphere energy uptake falls from 8.75 to 3.49 W/m². The atmosphere
is 1.36 K cooler and 4.07 kg/m² drier, its unclosed water
budget grows by 5.43 kg/m² over days 90–384, versus 2.17 in the control.
Ocean hot cells, coastal ice pile-up and whole-system energy closure remain
unresolved. CPU tests and a fresh-process GPU launch/restart check pass.

![ReadyESM D: ocean, sea ice, atmosphere and land at model day 390](docs/readiesm-v0.1.8-snapshot.png)

*A 390-day initialised run; v0.1.8 physics from day 90.*
[Timelines](docs/readiesm-v0.1.8-timelines.png) ·
[Figure provenance](docs/readiesm-v0.1.8-snapshot.json)

## Experimental roadmap

Snapshot: **14 October 2026**. These ten blue-sky tracks are experimental, not
features of the released v0.1.8 model. Passing a component test does not establish
coupled climate skill; releases need not wait for every track.

| Track | Direction | Progress | Next step |
| --- | --- | --- | --- |
| BS01 | Conservative convection, then a learned closure | Standalone plume runs reach 20 days, but upper-air moisture is biased high. No learned closure yet. | Port the plume to GPU after the cloud microphysics; test water and energy budgets. |
| BS02 | Prognostic clouds | A local liquid/ice store avoids the negative-condensate repair needed by the transported prototype. GPU port underway. | Couple stored condensate and cloud fraction to radiation; test restart and conservation. |
| BS03 | Virtual vertical layers | Bottom-layer refinement removes about 69% of the offline surface longwave error against the high-resolution reference. Coupled tests submitted. | Compare 30-day runs with matched controls, including heating and surface fluxes. |
| BS04 | Balanced initialisation | Remapping strongly reduces water drift in a 10-day test, but changes the initial inventory. More targeted spectral-cleanup tests submitted. | Qualify the GPU cleanup and compare fresh/restarted runs; gradual analysis increments later. |
| BS05 | Parameter calibration | Two 33-member real-column sensitivity screens completed; an expanded plume screen now separates moistening and drying terms. No fitted calibration yet. | Apply the diagnostics to the newer plume, then select a bounded coupled parameter sweep. |
| BS06 | Sea-ice thickness distribution | Five-category native transport and exact restart pass on a synthetic masked CPU grid. GPU qualification submitted. | Test mechanical feedback at the Arctic hotspot; the pile-up is not fixed yet. |
| BS07 | Coastal connectivity and landfast ice | Fine-bathymetry connectivity and metric-interface benchmarks pass around the problem coast. No live coastal stress scheme yet. | Translate resolved channel geometry into a local dynamics benchmark. |
| BS08 | Runoff storage and heat | Standalone routing closes water/heat budgets and restarts exactly. Receiver heat-replacement GPU test submitted. | Connect the land heat debit and source temperature without double-counting ocean heat. |
| BS09 | Automated experiment branching | Bounded planning and archived-failure replay tested. The live evidence adapter is under review. | Drive one reviewed shadow experiment with verified inputs and no duplicate submissions. |
| BS10 | Carbon, aerosols and emissions | Fixed-CO₂ radiation identity and closed carbon-box tests pass; dated greenhouse-gas coupled tests submitted. | Replay spatial air–sea carbon exchange with a closed ledger, before interactive coupling. |

## Run

The production configuration requires Julia 1.12 and an NVIDIA GPU with
working CUDA drivers.

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
python -m pip install -r scripts/requirements.txt
python scripts/download_era5_initial_state.py --include-stratosphere
julia --project=. scripts/run_dynamic_esm.jl config/production.yml
```

The ERA5 download requires CDS credentials. The January 1993 ECCO V4r4 files
download through ClimaOcean/NumericalEarth, with a public
[NumericalEarthArtifacts mirror](https://github.com/NumericalEarth/NumericalEarthArtifacts/releases/tag/data-v1)
when the primary ECCO download fails. Initial downloads require network access;
later runs reuse the cached inputs. `scripts/prepare_ecco_monthly_bundle.jl`
can verify and arrange a local monthly bundle for an explicit
`ecco_initial_conditions_directory`.

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
