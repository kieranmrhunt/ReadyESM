# Finite multi-plume convection

ReadyESM has two optional convection schemes for CPU runs:

- C08 represents three rising plumes with different entrainment rates.
- C07 adds vertical momentum transport and returns dissipated kinetic energy
  as local heat. It uses C08 for moisture and thermal convection.

Both use the normal model factory. They require Float32 atmospheric storage
and diagnostic clouds. The plume calculation and closure use Float64.
Production keeps Betts–Miller convection on the GPU.

## Configuration

Copy `config/production.yml` to a separate file. For a C08 experiment, change
these entries and add the three `atmosphere_c08_*` entries:

```yaml
name: finite_convection
device: cpu
duration_days: 1
analysis_start_day: 0
output_dir: artifacts/finite_convection
atmosphere_convection: finite_multiplume_c08_v1
atmosphere_convection_timescale_seconds: 7200
atmosphere_c08_entrainment_scale: 0.75
atmosphere_c08_precipitation_conversion_m_inv: 0.0015
atmosphere_c08_condensate_threshold_kg_kg: 0.0005
```

Keep the other production settings and input paths. Run the copied file with
`scripts/run_dynamic_esm.jl`. The initial-state files are the same as production.

For C07, change the scheme and add its momentum settings:

```yaml
atmosphere_convection: finite_multiplume_c07_v1
atmosphere_c07_momentum_strength: 0.5
atmosphere_c07_maximum_subcycles: 16
```

A momentum strength of zero retains the underlying C08 solution. Each
substep uses conservative momentum fluxes. The heat calculation follows
the same momentum path. Momentum transport is skipped when the required
number of substeps exceeds the configured maximum.

The current parameter experiments also test lower rain conversion and a larger
condensate threshold:

```yaml
atmosphere_c08_precipitation_conversion_m_inv: 0.00075
atmosphere_c08_condensate_threshold_kg_kg: 0.001
```

This is case 12 with a 7,200-second relaxation time. Case 3 uses 3,600 seconds.
Keep these choices explicit in each experiment configuration.

## Water, energy and output

C08 evaluates the response of the actual finite-step atmospheric state.
It tracks rain and snow separately. Detrained liquid and ice evaporate locally
with the corresponding phase-energy change. The scheme stores no persistent
convective condensate.

Output metadata records the scheme, parameter values, arithmetic precision
and heat-return policy. Convective rainfall and snowfall remain separate
fields. C07 records its momentum and energy checks in the model state.

## Current evidence

The clean C08 implementation matches 4,350 solver comparisons and 1,740 host
comparisons. The combined C07 implementation matches 2,610 host comparisons
at three parameter settings. The full CPU regression suite passes.

January 1993 atmosphere-only tests favour case 12 over the original C08
setting for moisture and lower-tropospheric temperature. Relative to native
convection, case 12 improves rainfall error but retains larger moisture and
500 hPa temperature errors. C07 improves global wind error relative to C08,
with mixed regional effects. See [experiment results](experiments.md).

July tests provide a separate seasonal check. Coupled integration, restart
checks and longer persistence tests are tracked in the same experiment notes.
