# Setup and tests

Run the commands below from the repository directory after installing the
dependencies listed in the [README](../README.md#run).

## Input data

The standard configuration starts on 1 January 1993. It uses ERA5 for the
atmosphere and land, and ECCO V4r4 for the ocean and sea ice.

### ERA5

Configure your Copernicus Climate Data Store credentials, then run:

```bash
python scripts/download_era5_initial_state.py --include-stratosphere
```

This writes pressure-level, single-level and land data under
`data/era5_initial_state`. The stratospheric levels are required by the
27-level atmosphere.

### ECCO

Download the December 1992 monthly files from NASA PO.DAAC. The preparation
script requires these four files and a published `.sha512` file alongside
each one:

| Collection | File |
| --- | --- |
| `ECCO_L4_TEMP_SALINITY_05DEG_MONTHLY_V4R4` | `OCEAN_TEMPERATURE_SALINITY_mon_mean_1992-12_ECCO_V4r4_latlon_0p50deg.nc` |
| `ECCO_L4_OCEAN_VEL_05DEG_MONTHLY_V4R4` | `OCEAN_VELOCITY_mon_mean_1992-12_ECCO_V4r4_latlon_0p50deg.nc` |
| `ECCO_L4_SSH_05DEG_MONTHLY_V4R4` | `SEA_SURFACE_HEIGHT_mon_mean_1992-12_ECCO_V4r4_latlon_0p50deg.nc` |
| `ECCO_L4_SEA_ICE_CONC_THICKNESS_05DEG_MONTHLY_V4R4` | `SEA_ICE_CONC_THICKNESS_mon_mean_1992-12_ECCO_V4r4_latlon_0p50deg.nc` |

The [preparation script](../scripts/prepare_ecco_monthly_bundle.jl) records
the PO.DAAC source URLs. Replace `/path/to/ecco-december-raw` below with the
directory containing the downloaded files.

```bash
julia --project=. scripts/prepare_ecco_monthly_bundle.jl \
    /path/to/ecco-december-raw data/ecco_v4r4_1992_12 1992-12
python -I scripts/download_ecco_january_mirror.py \
    data/ecco_january_mirror_1993_01
python -I scripts/prepare_calendar_centred_ecco.py \
    data/ecco_v4r4_1992_12 data/ecco_january_mirror_1993_01 \
    data/ecco_calendar_centred_1993_01
```

The January downloader uses the public
[NumericalEarthArtifacts mirror](https://github.com/NumericalEarth/NumericalEarthArtifacts/releases/tag/data-v1).
The prepared ocean state is the midpoint of the December and January
monthly means. Sea ice uses the January state.

The run script checks input hashes before model construction. It also checks
five generated ocean-fill caches before stepping. Those cache hashes were
tested on an A100. If another architecture produces different cache bytes,
the run stops at that check. See the
[v0.1.10 validation record](validation-v0.1.10.md).

## Run and output

```bash
julia --project=. scripts/run_dynamic_esm.jl config/production.yml
```

The supplied configuration requests 730 days, with diagnostics averaged from
day 365. For a shorter run, copy the file and set `analysis_start_day` below
`duration_days`. Also set a separate `name` and `output_dir`.

The runner writes NetCDF diagnostics, checks the simulated fields, then plots
the results. The [README](../README.md#configuration) lists common
switches. The complete standard settings are in
[`config/production.yml`](../config/production.yml).

## Tests

The CPU regression suite uses synthetic initial conditions. It needs neither
ERA5 credentials nor a GPU:

```bash
julia --project=. --startup-file=no test/runtests.jl
```

The Fourier graph-cache and conservative regridding tests need a CUDA GPU,
with no input downloads:

```bash
julia --project=. --startup-file=no test/cuda_graph_cache.jl
julia --project=. --startup-file=no test/csr_regridding_gpu.jl
```

With the production inputs and a CUDA GPU, run the coupled and restart checks
in separate Julia processes:

```bash
julia --project=. scripts/validate_coupled_smoke.jl artifacts/release-smoke
julia --project=. scripts/validate_coupled_restart.jl \
    artifacts/release-smoke/hour1_restart.jld2 artifacts/release-restart
julia --project=. scripts/compare_dynamic_restart_boundaries.jl \
    artifacts/release-restart artifacts/release-smoke/hour1_restart.jld2 \
    artifacts/release-restart/restored_boundary.jld2
```

These checks run one simulated hour, then restore the checkpoint and run a
second hour. They compare the restored state with the checkpoint and check
that the continuation remains finite.
