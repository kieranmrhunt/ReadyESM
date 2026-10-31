# Climate output

ReadyESM includes an opt-in streamed climate-output runner at
`scripts/run_streamed_esm.jl`. It records monthly and daily climate fields,
annual ocean fields, profile and energy diagnostics, and CMIP-style atmosphere
and ocean products while the coupled model runs. The ordinary
`scripts/run_dynamic_esm.jl` workflow is unchanged.

## Run modes

The command line is:

```sh
scripts/run_streamed_esm.jl CONFIG.yml off|full|segment|resume OUTPUT_DIR [STOP_ITERATION] [COMMITTED_CHECKPOINT_DIR]
```

The checkout must be clean. Set `READYESM_EXPECTED_MODEL_SOURCE` to the
40-character Git commit hash for that checkout before launching; the runner
checks both the commit and clean working tree. `OUTPUT_DIR` must not already
exist.

- `off` runs the coupled model without constructing or writing climate output.
- `full` runs the configured duration and writes the complete output set.
- `segment` runs through `STOP_ITERATION`, before the configured final step,
  and writes a committed restart pair.
- `resume` continues from a committed checkpoint directory through the
  requested `STOP_ITERATION`.

For example, after preparing a clean checkout and environment:

```sh
export READYESM_EXPECTED_MODEL_SOURCE="$(git rev-parse HEAD)"
julia --project=. scripts/run_streamed_esm.jl config/my-run.yml full output/my-run
```

GPU configurations require a functional CUDA device. Production ECCO launches
are checked against the qualified calendar-centred January 1993 input bundle.

## Output conventions

The output includes source and configuration hashes, field schemas, sampling
timing, and aggregation conventions in its metadata. Atmospheric means use the
model's Gaussian ring quadrature. Ocean means use unique active surface-cell
areas, excluding duplicate center-fold storage. Profile means retain coverage
and should be scored only when coverage is complete. Launch-time values that
are not yet available are recorded as missing and excluded from period means.
CMIP three-dimensional ocean fields are annual; the listed surface and
overturning fields are monthly. The writer's checkpoint snapshot is an
output-only endpoint sample and does not alter persistent model restart
history.

The bundled converter and native-velocity rotation helper are under
`scripts/climate_writer/cmip/`; run `python3 -I
scripts/climate_writer/cmip/cmorize_v2.py --help` for conversion arguments.
Writer and CMIP component tests are included from `test/runtests.jl`.
