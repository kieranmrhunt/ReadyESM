# Runoff heat and prescribed trace gases

Tested model source: `5ae5d603ff9291379a94ba7c0954661de5075fca`.
Base main: `1cd7d46db68007a8fc37a2b0f09a967eaab1a5bc`.
The production configuration and dependency manifest are unchanged.

## Scope

`terrarium_runoff_enthalpy: true` routes land-derived water enthalpy through
the runoff stores and into the mutable z-star ocean. It requires Terrarium.
The option adds store and transfer histories to output and checkpoints.

`forcing.ch4_ppb` and `forcing.n2o_ppb` set fixed, well-mixed concentrations
for RRTMGP. Values must be finite, positive and representable in the model's
floating-point type. Omitted values leave RRTMGP's 1900 ppb methane and
340 ppb nitrous oxide unchanged. The model records requested and applied
values. These controls prescribe concentrations, without adding gas chemistry
or an emissions cycle for either gas.

## Checks

| Test | Result | Job |
| --- | --- | --- |
| Combined off versus main | Exact stored state and NetCDF output | 61056251_0 |
| Combined runoff versus runoff-only source | Exact state, output and ledger | 61056251_0 |
| Combined gases versus gases-only source | Exact state, output and gas values | 61056251_1 |
| Both options and fresh-process restart | Exact restart at step 3, then continuation through step 8 | 61056251_1 |
| Full CPU suite with ERA5 inputs | 44 testsets, including 97 spectral-initialisation checks, no failures or skips | 61065689 |
| Combined readout | 39 checks pass with corrected state-derived attribute comparison | 61087310 |
| Full-grid A100, both options | Eight steps, valid gas values and runoff ledger, arrays on GPU | 61070683 |
| Device-probe positive and negative controls | Pending | 61073955 |

The compact comparison uses T7/L8 and a 32×16×4 ocean. The GPU test uses
the production T31/L27 atmosphere and 360×180×60 ocean. Its routed-enthalpy
ledger closes to a relative residual of 6.3e-9, with 3.51e17 J released.
The tests include output collection without changes to model state.

The original combined reader failed because it required the in-transit runoff
heat attribute to stay fixed when greenhouse-gas forcing changed. The corrected
reader ties that attribute exactly to each run's endpoint receiver rate and
checks it against released minus received heat. All other comparisons remain
in force. The original failure and ten negative tests are retained.

An earlier GPU checker mistook a view of a GPU array for host storage. The
corrected checker follows the underlying storage and records its type.

## Runtime and climate follow-up

The local kernel speed test exceeded its 1.25 ratio limit. The earlier
four-block whole-model comparison gave a mean ratio of 1.0044 and passed its
1.03 criterion. A later short pair reported 1.076, with no cost verdict.
These measurements are retained separately because their harnesses and source
revisions differ. The optional publication decision accepts this runtime
trade-off. The failed local criterion remains recorded as failed.

Production leaves runoff heat disabled and keeps the previous gas values.
Four matched monthly runs are measuring runoff heat's effects on coastal
temperature, salinity, radiation and water and energy budgets. Their results
will inform any proposed default change.
