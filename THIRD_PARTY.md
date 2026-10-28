# Third-party licences

ReadyESM is released under the European Union Public Licence v1.2
(EUPL-1.2). See [LICENSE](LICENSE).

The `vendor/` folder holds modified copies of three packages. The table
records each copy's location, upstream licence and the local changes. Each
copy keeps its own licence file. ReadyESM's EUPL-1.2 notice does not replace
those package licences.

| Package | Location | Licence | What ReadyESM changed | Licence file |
| --- | --- | --- | --- | --- |
| ClimaSeaIce 0.5.8 | `vendor/ClimaSeaIce/` | Apache-2.0 | Used computed result types for zero alternatives in split-explicit velocity and EVP stress increments, and for the snow-rebasing zero alternative. | [vendor/ClimaSeaIce/LICENSE](vendor/ClimaSeaIce/LICENSE) |
| SpeedyTransforms 0.1.5 | `vendor/SpeedyTransforms/` | EUPL-1.2 | The GPU Fourier transform cache now uses the scratch arrays, their sizes and their strides as its key. It keeps the captured input and output arrays alive. It waits for every recorded CUDA context before it is cleared. | [vendor/SpeedyTransforms/LICENSE](vendor/SpeedyTransforms/LICENSE) |
| Terrarium 0.1.6 | `vendor/Terrarium/` | EUPL-1.2 | Soil processes now start in a fixed order. The correction for oversaturated soil moves water by porosity and layer thickness. The skin-temperature solver works in Kelvin. Evapotranspiration reads the layer thicknesses from the grid directly. A helper function is defined before the generated method that calls it. | [vendor/Terrarium/LICENSE](vendor/Terrarium/LICENSE) |

Some files in `src/` adapt code from SpeedyWeather, which uses EUPL-1.2.
For example, the full-temperature vertical advection follows SpeedyWeather
[pull request 1286](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/1286),
merged as commit `adf0d61da03bc3a9bd6698edb9fc0a7ddfac002d`.

Julia installs the other packages listed in `Manifest.toml`. Each package
keeps its own licence. Most use the MIT licence. Others use Apache-2.0,
BSD, MPL-2.0, LGPL, EUPL-1.2 or other licences.
