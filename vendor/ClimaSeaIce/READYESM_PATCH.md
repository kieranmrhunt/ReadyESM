# ReadyESM patch

Based on ClimaSeaIce 0.5.8, tree
`220b3720187532df0ed33bc48fb262c342e31bb8`. Runtime source, package metadata,
README and license are retained; upstream examples and tests are omitted.

The split-explicit velocity kernels promote the already-computed dynamic and
free-drift alternatives before selecting between them. Their zero branch uses
that same type. EVP stress increments likewise use a zero of the increment
type. This avoids mixed Float32/Float64 branch results and the literal pointer
expressions observed in the affected GPU specializations. Computed tendencies,
physical parameters, native launches and momentum ordering are unchanged.

Layered thermodynamics uses the computed snow-rebasing type for its zero
alternative. The ice-volume helper is deliberately unchanged: promoting its
returns can create an ice mass flux even when the stored state is unchanged.
The narrow snow change passes exact stored-state/flux and zero-tendency checks.

The original and candidate CPU kernels agree exactly in production-state and
small-grid comparisons. Coupled GPU reliability remains under qualification;
this patch alone is not proof that CUDA700 is resolved.

## Opt-in lead-closing thickness (3 Nov 2026)

`SeaIceThermodynamics.LeadClosingEvolution(h0)` is a `concentration_evolution` rule that is the same as
`ProportionalEvolution`, except that its freezing term uses the lead-closing thickness `h0` and not the consolidation
thickness `hc`. `freezing_concentration_thickness` returns `hc` for `ProportionalEvolution`, so the closed-form snow
step and the concentration step keep their released arithmetic exactly. ReadyESM selects the new rule only with
`sea_ice_lead_closing_thickness_m > 0`.

## Opt-in hemispheric lead-closing thickness (4 Nov 2026)

`SeaIceThermodynamics.HemisphericLeadClosingEvolution(north, south)` picks the lead-closing thickness per cell from
the sign of the cell-centre latitude (`φnode`, equator counted as north). A zero value keeps `hc` in that hemisphere.
`local_concentration_evolution(evolution, i, j, grid, hc)` returns the cell's rule. It returns the same object for
every other rule, and it returns `LeadClosingEvolution(h0)` (with `h0 = hc` for a zero value) for the hemispheric rule.
The slab and layered kernels now call it once per cell and pass the result to `freezing_concentration_thickness` and,
as `(; concentration_evolution = evolution)`, to `ice_volume_update`. For the existing rules this gives the same
arithmetic. ReadyESM selects the new rule only when `sea_ice_lead_closing_thickness_south_m` is set.
