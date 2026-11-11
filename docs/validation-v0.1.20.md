# v0.1.20 validation record

## Release decision, 11 November 2026

v0.1.20 uses the production physics at source commit
`fba77cd3ffc42f0f9227f9cd2ca24840bb8ad396`, tree
`b735b45610c67e24e0d99cf9d31b4e2f894f4d2c`. The release branch
contains commit `3605b2f856da3f964416afc424249024283e9a2b`, which
updates six test files. It leaves the model source and
production configurations used by the GPU run unchanged. The release merge
also retains the published main branch's CI change. Release metadata and
figures were added after the tests. Neither changes the physical model.

The accepted scientific comparison is one continuous January–June 1993 run
per version from the same observed starting state. It uses only locations
with valid data for every model. The v0.1.20 run improves the eight-field
overall normalised error ratio from 1.20 (published v0.1.19) to 1.16.
Lower is better. The typical model among 55 CMIP6 models scored 1.04.
Those models supplied 1,650 model-years of comparison data. The earlier
v0.1.20 candidate with several new atmosphere physics options scored 1.60
and was rejected. The released run improves the average error across all
nine reported fields. Some regions worsened. The Southern midlatitude
500-hPa temperature ratio rises from 1.85 to 2.10, and the subtropical
850-hPa zonal wind ratio rises from 1.18 to 1.40. Those ratios are against
the same observational fields as the overall score.

Locations with any missing data are dropped from the comparison. This can
make every model look slightly better. Values at 850 hPa are omitted where
that pressure level is below ground or lacks a valid value. The retained
fraction of the ERA5 850-hPa area ranges from 0.77 in the southern polar
region to almost 1.00 in the southern midlatitudes. The other seven fields
retain their ERA5 masks. This score describes agreement with six observed months
following one observed initial state. It does not establish repeatability,
equilibrium, or long-term climate skill. The [scorecard](../scorecard.png),
[machine-readable figure record](scorecard-v0.1.20-partc.json), and
[plot recipe](plot-v0.1.20-scorecard.py) show the compared values.

## Exact run and score receipts

The A100 six-month production run, Slurm job `66841309`, completed with exit
code 0, 181 days and 17,376 model steps. Its result-manifest SHA-256
is `5be8e9032614f7a2785bac0fd04785ad93750402b2ee2e90a757f22dc40cd118`.
The run recorded six January–June monthly windows and 18 scored NetCDF
outputs. The run-input manifest SHA-256 is
`5e5ca394b0f4437d966665bf0bdf8c662db5b8f6470d96814891d30c1f71b244`.
The production six-month config has the same physical settings as
[`config/production.yml`](../config/production.yml), changing only the run length, analysis start day,
name and output path. The January [monthly figure](readiesm-v0.1.20-january-monthly.png)
comes from that run. Its [figure record](readiesm-v0.1.20-january-monthly.json)
and [plot recipe](plot-v0.1.20-january-monthly.py) record the sample window
and hashes. The figure shows atmosphere and top-of-atmosphere radiation,
not ocean, ice or land.

The joint six-month scorer, job `66904500`, completed with exit code 0.
Its result-manifest SHA-256 is
`802ec3617e06ee7b0e6b60aaee26c7aea33c0bcfbc76c49f0fcffa1bcb414c53`.
The scored JSON SHA-256 is
`2bfcdf0756a2a6359a0ede6ac333f80cb8f33330412412c2e5bf547476cbc856`.
The comparison uses the published v0.1.19 run `66742450`, released v0.1.20
run `66841309`, and rejected earlier candidate run `66632474`. Its input
manifest SHA-256 is
`22ff02c2808dd9031abbc37226bbe9ce8757eecec7638ebd8e8c65cfa1d40417`.

The full CPU suite, job `66856523`, passed on the released model code. Its audit
JSON SHA-256 is
`3ba4348b478974857d979e0a97ac5ffa903c2f92ba2e3d0f41960a035ee918a9`.
The result-manifest SHA-256 is
`18dd99e0e9fa41b528ed43d5f762329cc957cd2fb60520be6a161e62e93d9e13`.
Focused tests, job `66840968`, also passed. Their audit JSON SHA-256 is
`bb92232a3d075389c7d9dddd5556e5e78bbf8ec2b170872505306dac0fc97f09`.
The provenance check linked both test results to the production source,
its two configuration hashes and the exact six-month GPU run. The binding
JSON SHA-256 is
`e5d7c7ff07f37e9bd71e365898483ddc2571eb1f12e33b48258d01338bebfeab`.
These receipts establish provenance and execution. Scientific acceptance is
the separate comparison and review above.

## Production settings and unresolved checks

The production profile keeps v0.1.19 atmosphere physics: two-hour
Betts–Miller convection, mixed-phase saturation, momentum-only K-profile
boundary-layer mixing, no large-scale snow-melt cooling, diagnostic cloud
paths of 75/25 g m⁻², and fixed 420 ppm CO₂. It uses 0.5 m sea-ice lead
closing, prepared MACv2/MACv2-SP/UOEXETER prescribed aerosols, and the
1360.7 W m⁻² CMIP6 piControl solar constant. Ozone is analytic.
The prepared aerosol input files are required and are not distributed in
the repository. See the [input guide](../data/README_aerosol.md).

The paired January and July month test asked for more pairs before a
decision. They were not run for this release. The v0.1.19
[six-pair score](validation-v0.1.19.md) remains
`REJECT`, including its energy and Antarctic sea-ice failures. This release
does not retrospectively change that result. A clear-sky outgoing-longwave
gap is present at the January start. A provisional upper-atmosphere
correction is under separate evaluation for a later release. A one-step
CPU-versus-GPU check found a larger difference than allowed in one land
field. A production-sized GPU restart has not been
qualified. Restart evidence is from a smaller grid, and runs with
`scored_outputs: true` cannot restart. Multi-year energy balance, ocean and
sea-ice drift, and repeatability from other starts still need testing.
