# ReadyESM

ReadyESM couples SpeedyWeather, RRTMGP, Oceananigans, ClimaSeaIce and Terrarium.
The current model has a T31/L27 atmosphere and a 1° tripolar ocean with 60
levels, together with dynamic sea ice and a 16-layer land model. ERA5 and ECCO
provide its initial state.

ReadyESM is a research model under active development.
The production configuration is [`config/production.yml`](config/production.yml).

`v0.1.10` keeps the v0.1.9 model physics and changes the ocean launch state to
the midpoint of the December 1992 and January 1993 ECCO monthly means. Sea ice
stays on the January state. In three matched 30-day pairs, days 1–10 SST RMSE
falls by 0.208 K. Over the same window, Antarctic ice-area error improves by
24% and Antarctic extent error worsens by 15%. All pairs pass the Antarctic
regrowth guard. Ocean freshwater closure stays within its 0.01 bound.

`v0.1.9` removes an invalid frozen spectral row from the ERA5 launch state.
In a matched fresh 30-day pair, this reduces the magnitude of atmospheric-water
drift by 94.8% and removes a spurious 80–90 hPa polar pressure depression. Old
checkpoints keep their original launch state until explicit migration. The
model retains two-hour convection, mixed-phase saturation, double-precision
land and liquid and ice cloud paths of 75 and 25 g/m². See the
[release history](HISTORY.md) and
[v0.1.10 validation record](docs/validation-v0.1.10.md).

The longest accepted trajectory remains the pre-fix v0.1.8 run to day 390.
The figure below shows this run. v0.1.10 has three matched 30-day tests and a
passing day-one production gate. Current work targets ocean hot cells, coastal
ice pile-up and whole-system energy closure.

![ReadyESM D: ocean, sea ice, atmosphere and land at model day 390](docs/readiesm-v0.1.8-snapshot.png)

*A 390-day initialised run. It uses v0.1.8 physics from day 90.*
[Timelines](docs/readiesm-v0.1.8-timelines.png) ·
[Figure provenance](docs/readiesm-v0.1.8-snapshot.json)

## Experimental roadmap

Snapshot: **19 October 2026**. Production remains at v0.1.10.
The tables below track research candidates. See the
[experiment notes](docs/experiments.md) for the convection portfolio and test scope.

| Track | Direction | Progress | Next step |
| --- | --- | --- | --- |
| BS01 | Conservative convection, then a learned closure | The tested single-plume form is closed. Its conservative core supports the newer convection experiments. | Continue C08 and C10. |
| BS02 | Prognostic clouds | Transport, water storage and restart components pass. The tested cloud closure worsens radiation. | Replace cloud fraction and overlap before another screen. |
| BS03 | Virtual vertical layers | The partial-layer follow-up failed its registered effect tests. This form is closed. | Reopen with a clear observational benefit. |
| BS04 | Balanced initialisation | Accepted in v0.1.10. Three matched pairs improve early SST error, with mixed Antarctic ice results. | Extend the accepted launch beyond 30 days. |
| BS05 | Parameter calibration | Column screens establish parameter sensitivity. | Select the coupled parameter set after convection and cloud structure tests. |
| BS06 | Sea-ice thickness distribution | The tested pressure hybrid is closed. Category transport and restart tools are retained. | Reopen after a defined coastal-growth or thickness trigger. |
| BS07 | Coastal connectivity and landfast ice | The tested opening and generic coastal stress reduce target growth but worsen or relocate total growth. This form is closed. | Test geometry-specific landfast physics after a growth trigger. |
| BS08 | Runoff storage and heat | Component and restart tests pass. A new six-release/six-candidate GPU repeat test is running. | Complete the repeatability, cost and coupled tests. |
| BS09 | Automated experiment branching | The tested incident-recovery tool is complete and archived. | Reopen for a general tool with resource and environment checks. |
| BS10 | Carbon, aerosols and emissions | The live carbon adapter works on CPU and GPU. A longer GPU case fails its inventory check. CPU diagnostic repeats pass. | Diagnose the running GPU repeats with gas exchange on and off. |

Current convection work:

| Track | Progress | Next step |
| --- | --- | --- |
| C06 | Rain evaporation and cold-pool memory pass the ten-day atmosphere screen. All five full-ESM GPU pilots are running. | Assess the coupled results, then run the queued ten-day repeats. |
| C07 | Conservative momentum transport and heat return pass component and atmosphere tests. | Complete four full-ESM pilots and the 30-day ERA5-start comparisons. |
| C08 | The revised closure passes its eight-step full-ESM pilot. Thirty-day runs expose rare energy residuals and one host-water residual. | Trace the returned exchange and complete coupled restart and day-one tests. |
| C10 | All nine two-day tests pass. Eight to thirty-two internal steps reduce peak temperatures from about 358 K to 324 K. | Complete the six running ten-day tests and six full-model pilots. |

## Run

The production configuration requires Julia 1.12 and an NVIDIA GPU with working
CUDA drivers.

```bash
julia --project=. -e 'using Pkg' -e 'Pkg.instantiate()'
python -m pip install -r scripts/requirements.txt
python scripts/download_era5_initial_state.py --include-stratosphere
julia --project=. scripts/prepare_ecco_monthly_bundle.jl \
    /path/to/ecco-december-raw data/ecco_v4r4_1992_12 1992-12
python -I scripts/download_ecco_january_mirror.py \
    data/ecco_january_mirror_1993_01
python -I scripts/prepare_calendar_centred_ecco.py \
    data/ecco_v4r4_1992_12 data/ecco_january_mirror_1993_01 \
    data/ecco_calendar_centred_1993_01
julia --project=. scripts/run_dynamic_esm.jl config/production.yml
```

The ERA5 download requires CDS credentials. The December ECCO step requires
the four PO.DAAC V4r4 monthly files and their published `.sha512` files. The
January step uses the public
[NumericalEarthArtifacts mirror](https://github.com/NumericalEarth/NumericalEarthArtifacts/releases/tag/data-v1).
The preparation scripts verify every source and output hash. The production
entry point is `scripts/run_dynamic_esm.jl`. It checks the NetCDF files before
construction and checks the five generated inpainting caches before the first
time step. Direct API calls check file presence but skip these release hashes.
The cache bytes were qualified on A100. Another architecture that produces
different cache bytes stops before the first time step.

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
