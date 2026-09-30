# v0.1.5 candidate validation

This candidate changes numerical branch types in ClimaSeaIce, not the model
physics or production configuration. It is not published yet. The 48-hour
attempt has failed; no new long-run figure was generated.

## Change

The velocity kernels promote the computed dynamic and free-drift alternatives
to a common type before selecting between them; the zero alternative uses
that type. EVP stress updates likewise use the increment's zero type.
This removes five literal pointer expressions from the captured GPU kernel
specializations. The computed tendencies and native stepping order are unchanged.
ClimaSeaIce 0.5.8 is vendored with changes to two runtime source files.

## Completed checks

- Production-state CPU replay: 2,952,600 stored values agree exactly with
  the original kernels, including matching auxiliary NaNs.
- Portable CPU comparison: nine precision/ice cases, 792 assertions and
  12,096 exactly matching values.
- Packaged public CPU suite: passes, including 298 new branch checks. Native
  package precompilation/loading and offline GPU specialization checks pass.
- GPU candidate job 56380978: isolated 900-second ice step, four resumed
  coupled steps, full diagnostics and exact restored checkpoint state pass.
  This first run used a startup code overlay.
- Native-package job 56416366: four resumed coupled steps through model hour
  two, full diagnostics, output/checkpoint saving and independent exact
  restored-state comparison pass, without an overlay or diagnostic observer.
  This ran on a different physical A100 from the first candidate pass.
- Clean release checkout, job 56435805: all public CPU tests pass, including
  298 branch checks, plus 489 source-identity and five package-loading checks.
  It needs no private forcing files. The tested state is local commit
  `dca067d`; subsequent changes here only record the qualification results.

The frozen native candidate has source-manifest SHA-256
`474bdf0c1b66f930ee862f0d97f147b3ed2375864103a75a0e957dfc6452374a`.
Release staging changes only ReadyESM version metadata and documentation
beyond that tested source. Its separate CPU recheck passes in 19m12s.

## Remaining qualification and limits

Job 56422865 attempted to continue the native-package checkpoint from hour two
to hour 48, using unchanged resolution, parameters and diagnostic thresholds.
It failed with CUDA700 after 1h04m13s in its first resumed step (iteration 9);
no resumed step completed. The error was detected during context synchronization
before loading the ice-consolidation kernel, after the host had submitted the
dynamic and thermodynamic updates. The initiating operation is not identified.

The restored checkpoint bytes exactly match the input. This verifies saved
state, not internal work buffers or execution ordering. The dependent figure
job cancelled without running, and the README retains its labelled v0.1.4 image.
The failure demonstrates that the momentum/stress type change is not sufficient
to eliminate ordinary coupled illegal accesses.

Two repeat allocations failed CUDA initialization before loading the model;
another native attempt lacked its input-data link. These are retained as
setup failures, not successful model runs. The v0.1.4 CUDA700 failure remains
in its [validation record](validation-v0.1.4.md).

Successful short runs do not prove the initiating cause of the earlier fault,
quantify its frequency, establish continuous-versus-restarted trajectory
equality or demonstrate a stationary climate. Experimental land coupling and
air–ice surface-solver changes are excluded from this candidate.
