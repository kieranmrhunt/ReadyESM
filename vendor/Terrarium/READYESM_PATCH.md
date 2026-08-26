# ReadyESM Terrarium 0.1.6 integration patch

This directory starts from the official Terrarium `v0.1.6` source archive.
The General registry identifies that release with
`git-tree-sha1 = 7fef90df5d9e676d1b67941ee66fd66750ee5d26`.

ReadyESM retains five focused correctness changes that are not present in the
release tag:

1. Coupled soil initialization diagnoses biogeochemical composition and soil
   thermodynamics before hydraulic conductivity, whose calculation depends on
   both. This prevents conductivity from remaining at its allocation value
   until the first auxiliary update.
2. Saturation-profile correction transfers physical water volume, including
   both source and destination porosity as well as layer thickness. Excess
   moved into `surface_excess_water` is likewise converted to physical water
   depth with top-layer porosity. The generic three-argument upstream method is
   retained; coupled soil closures select the porosity-aware four-argument
   method.
3. The nonlinear skin-temperature solve uses a Kelvin-shifted numerical
   coordinate while retaining Celsius in the state and physical equations.
   This avoids RootSolvers' relative Newton limiter treating an ordinary
   `0 °C` starting point as zero-scale, and rejects trial temperatures outside
   the positive-absolute-temperature domain so line search can backtrack.
4. The generic evapotranspiration forcing unwraps `ColumnGrid` and
   `ColumnRingGrid` before asking Oceananigans for layer thickness and the top
   index. The release implementation passes the wrapper directly, which has no
   `z` field and fails on the first soil-water tendency evaluation.
5. The generation-time helper for typed `get_fields` is defined before the
   generated method that calls it. Julia 1.12 can compile that method while the
   module is still loading; the release order then raises `UndefVarError` on
   the GPU path even though the helper appears later in the file.

ReadyESM separately specializes the exact coupled Richards/runoff configuration
in `src/terrarium_gpu_coupling.jl` to conserve infiltration, drainage,
evapotranspiration and GPU accumulation. Those host-model changes are kept out
of this vendored package because they are specific to ReadyESM's coupling and
ledger definitions.
