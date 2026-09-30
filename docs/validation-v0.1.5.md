# v0.1.5 candidate validation

This candidate changes numerical branch types in ClimaSeaIce, not the model
physics or production configuration. It is not published yet. The earlier
48-hour attempt failed; the newer snow change passes native CPU comparisons
and an ordinary GPU continuation to hour 48. This is not climate or multi-year
qualification.

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
checks using synthetic inputs. Offline thermal argument types remain
reconstructed, not independently matched to a live thermal GPU capture.

## Ordinary GPU continuation to hour 48

Job 56457732 completes in 59m03s on 30 September. It restores the earlier
hour-two checkpoint and advances 184 native 900-second steps to hour 48 at
T31/L27 and 360 x 180 x 60 ocean resolution. Native package loading, unchanged
full diagnostic validation, output/checkpoint saving and an independent exact
restored-state comparison all pass. There is no runtime overlay or diagnostic
observer. The saved restored checkpoint is byte-identical to its input.

Model construction takes 30m39s; the first resumed step takes 13m19s and the
remaining 183 steps take 2m17s in total. These are measured phases of this run,
not a general runtime forecast.

The successful run uses an A100 on gpuhost011 with driver 610.57.04. The failed
earlier attempt used a different A100 on gpuhost007 with driver 610.43.02.
The source and hardware/software environment both changed: success does not
isolate the patch's causal effect or quantify failure frequency. The separate
fresh-process restart described below also passes.

SHA-256 records:

| File | SHA-256 |
| --- | --- |
| Input and saved restored boundary | `403bfdf8eff0fd980ce2dd11e436b17987834cd8a78521e5b4ecaf16b63e51cb` |
| Hour-48 checkpoint | `d744ca05fb63c9d586d24a3b5161803f35c6b95074279ccff975b6d5e89da06b` |
| Diagnostics | `37a7bf9e41ede7908a866a71863471942afe0c3e546e3ba74a25f3954c07c524` |
| Coupled log | `ec024dc97c61635bd06afc0d0d1bc2db9a2808449bce2b070dcf0ec987f56b61` |
| Independent comparison log | `23b9dc6e80ad202eea8801538f890e9420555a727c400e626bfc0a7b56685cc6` |

Figure job 56469061 passes in 1m04s. The README maps and linked timelines use
this actual output, with [source and rendering provenance](readiesm-v0.1.5-snapshot.json).
Missing initial diagnostic samples remain missing, not zero. Both exported
figures were visually checked for units, labels, masks and clipping.

## Fresh restart of the hour-48 checkpoint

Job 56463210 completes in 56m24s on 30 September. A new native model process
restores the hour-48 checkpoint, advances four 900-second steps to hour 49,
passes full diagnostics and saves its output. An independent process verifies
exact restored model-state equality; the saved restored checkpoint is also
byte-identical to the input. The owned runoff buffers are poisoned before the
normal restore, as in the preceding check.

This uses the same A100 and driver 610.57.04 as the successful 48-hour run.
It verifies restoration and finite continuation, not equality with an
independent uninterrupted trajectory or general GPU failure frequency.

| File | SHA-256 |
| --- | --- |
| Input and restored boundary | `d744ca05fb63c9d586d24a3b5161803f35c6b95074279ccff975b6d5e89da06b` |
| Hour-49 checkpoint | `eab13cd0511ca0ac58b5f679f8f25b67756ee481753fff49991e2cb32e1a61c3` |
| Diagnostics | `42a5848302c24d9d271cee154de0323870dd7a7c2535a98babda09845d3ccfc8` |
| Coupled log | `7b20cbc399d13264d39ca9d3ae2dd78f009a32eb89458b0ae0033346f3e342f2` |
| Independent comparison log | `23b9dc6e80ad202eea8801538f890e9420555a727c400e626bfc0a7b56685cc6` |

## Retained failures and limits

GitHub CPU run [36777929865](https://github.com/kieranmrhunt/ReadyESM/actions/runs/36777929865)
ended with a hosted-runner communication failure after dependency installation
passed. GitHub retained no downloadable job log. This is not a passing CI run
or an identified test failure; the cause of the runner loss is unknown. The
local CPU regression passes above are separate evidence.

Before the snow change, job 56422865 attempted to continue the native-package checkpoint from hour two
to hour 48, using unchanged resolution, parameters and diagnostic thresholds.
It failed with CUDA700 after 1h04m13s in its first resumed step (iteration 9);
no resumed step completed. The error was detected during context synchronization
before loading the ice-consolidation kernel, after the host had submitted the
dynamic and thermodynamic updates. The initiating operation is not identified.

The restored checkpoint bytes exactly match the input. This verifies saved
state, not internal work buffers or execution ordering. The dependent figure
job cancelled without running. The README now uses the later successful run.
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
