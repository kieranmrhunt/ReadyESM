# Seasonal solar forcing

Add this setting under the existing `forcing` section:

```yaml
forcing:
  solar_distance_correction: true
```

RRTMGP then varies incoming sunlight with seasonal Earth–Sun distance.
The calculation uses the model date and a present-day orbital approximation.
It follows the fraction of each Gregorian year, including leap years.
The solar constant and radiation interval stay unchanged.

The default is `false`, which retains constant-distance forcing. The
option requires RRTMGP. Long-term orbital changes require a separate scheme.

## Restart

Keep the setting unchanged when you resume a checkpoint. The model checks
it before copying saved state. Older checkpoints restore with the option off.
CH4 and N2O have separate checks. See the
[trace-gas notes](validation-runoff-trace-gases.md#scope).

## Checks

The full CPU suite passes. Three coupled CPU tests check unchanged results
with the option off, agreement with the earlier solar calculation when on,
fresh-process restarts and rejection of changed forcing settings.

The A100 test uses T31 with 27 atmosphere levels and a 360×180×60 tripolar
ocean. It runs eight coupled steps and restarts from step three, between
radiation calls. All continued states and diagnostic fields agree exactly.
Device checks confirm GPU storage and disable scalar indexing.

These tests use source `6b9313b`. Their job IDs are CPU suite `61093863`,
coupled array `61094976` and GPU test `61110025`.
Monthly climate comparisons are tracked in the [experiment notes](experiments.md).
