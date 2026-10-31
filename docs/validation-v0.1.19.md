# v0.1.19 candidate validation record

This is a working record for the release candidate. It is not a release verdict.
The final tagged source and its exact inputs must replace the candidate
identities below before publication.

## Production selection

The versioned `config/production.yml` profile selects atmospheric
`kinetic_energy_v1` dissipative heating, `era5_spectral` land albedo and
`cover_lai` seasonal vegetation. The generic `ExperimentConfig` defaults
remain `off`, `constant` and `off` so ad hoc configurations do not silently
require the data-backed options. The albedo input SHA-256 is
`2f8e9193e2fea340e81dacc0e486c1d7de7f3b3a1b1e44dfc16ae3eef76f5c0c`;
the vegetation input SHA-256 is
`b80d7e01a4390340794030f1f924e5583d813fb163a1f9d875ae749091c13f60`.
Both are packaged in `data/` and the albedo derivation, its erratum and limits
are documented there.

## Evidence required before release

| Gate | Current state |
| --- | --- |
| Exact packaged-input hashes and effective production config | Tag-descended candidate passes 52/52 assertions and both packaged `SHA256SUMS` sets; recheck on final source. |
| Component tests for albedo, vegetation and dissipative heating | Albedo and vegetation pass. The first D1 run used a stale compact fixture that inherited the new production options; the corrected fixture is under retest. |
| Six-hour CPU smoke with all three options active | V5b successor preflight passes on source d325fc3; trio smoke running, with receipt and independent review pending. |
| Same-start 90-day July full-set run versus the v0.1.18 baseline | V5b launch and readout pending. This tests the combination; it cannot attribute any change to one option. |
| Energy, vertical-profile, land/ocean and regional guardrails | Pending V5b readout. The prior albedo-only month shifted global net top-of-atmosphere flux by about +9.65 W m⁻², so temperature improvement alone is insufficient. |
| Exact-final-candidate January and July paired release months, scores and safety checks | Pending. The v0.1.18 release record retains its rejected seven-field and neutral-rule results. |
| Full CPU suite, GPU restart/output checks, final figures and two independent scoped reviews | Pending on the final candidate. |
| Drift-screen calculation unit tests | Nine Python `unittest` cases pass on the tag-descended candidate. |
| Final clean tree, tag, remote readback and post-tag figure audit | Pending. |

The [v0.1.18 record](validation-v0.1.18.md) is the previous release baseline.
The v0.1.19 production choices have not been shown to cure temperature drift
or close the long-run energy budget.
