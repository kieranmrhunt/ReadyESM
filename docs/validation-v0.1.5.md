# v0.1.5 candidate validation

This candidate changes numerical branch types in ClimaSeaIce, not the model
physics or production configuration. It is not published yet. The earlier
48-hour attempt failed; a newer snow change passes native CPU comparisons
and awaits coupled GPU testing. No new long-run figure was generated.

## Change

The velocity kernels promote the computed dynamic and free-drift alternatives
to a common type before selecting between them; the zero alternative uses
that type. EVP stress updates likewise use the increment's zero type.
This removes five literal pointer expressions from the captured GPU kernel
specializations. Layered thermodynamics additionally uses the computed snow
rebasing type for its zero alternative, removing one literal pointer from
the reconstructed offline thermal specialization. The ice-volume helper is
unchanged. ClimaSeaIce 0.5.8 is vendored with changes to three runtime files;
physical parameters and native stepping order are unchanged.

## Earlier momentum/stress candidate

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
  `dca067d`, before the snow change was added.

That frozen native candidate has source-manifest SHA-256
`474bdf0c1b66f930ee862f0d97f147b3ed2375864103a75a0e957dfc6452374a`.
Its separate staged-release CPU recheck passes in 19m12s.

## Added snow change

Native-package job 56457539 passes in 6m59s, without runtime overlays:

- All 3,374,400 values across eight thermal state/flux fields and six timestep
  cases match the original CPU reference exactly, with no nonfinite values.
- All 456 portable thermal checks pass, including zero-tendency mass flux.
- The full offline thermal kernel contains no literal pointer expressions.

A broader promotion of ice-volume helper returns was rejected: it creates
artificial mass flux in 48 of 144 zero-tendency cases despite identical stored
thickness/concentration. The snow-only change leaves that helper untouched.

The updated 245-file source-manifest SHA-256 is
`b9c8c228703f860882999583f3d2c1989768a99848f525624509c5415f3f6111`.
Release staging differs only in version metadata and documentation; all 490
staging identity checks pass. Public CPU suite 56457566 passes in 21m05s,
including 1,126 named regression assertions and 486 package source-identity
checks, without private forcing data. Ordinary GPU continuation 56457732 is
waiting for resources, not yet qualified. Offline thermal argument types remain
reconstructed, not independently matched to a live thermal GPU capture.

## Remaining qualification and limits

Before the snow change, job 56422865 attempted to continue the native-package checkpoint from hour two
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
