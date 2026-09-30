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

The original and candidate CPU kernels agree exactly in production-state and
small-grid comparisons. Coupled GPU reliability remains under qualification;
this patch alone is not proof that CUDA700 is resolved.
