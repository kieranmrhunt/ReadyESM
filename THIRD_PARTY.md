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
The optional snowmelt cooling in `src/large_scale_melt_cooling.jl` copies the
`large_scale_condensation!` routine of SpeedyWeather 0.21.1
(`src/parameterizations/large_scale_condensation.jl`) and changes two lines,
so that melting snow cools the air. The file header names the copyright
holders and the licence.

Julia installs the other packages listed in `Manifest.toml`. Each package
keeps its own licence. Most use the MIT licence. Others use Apache-2.0,
BSD, MPL-2.0, LGPL, EUPL-1.2 or other licences.

## Bundled CMIP6 forcing data

The two NetCDF files in [`data/forcing_1993/`](data/forcing_1993/README.md)
are derived from CMIP6 input4MIPs data and retain their source attribution.
The source collections below are listed as **CC BY-SA 4.0** by the World Data
Center for Climate. ReadyESM's EUPL-1.2 code licence does not replace the
source-data licence. The NetCDF attributes and the forcing README record the
preparation method, original file identifiers and checksums.

| Bundled file | Source collection and version | Original creators and citation | Preparation |
| --- | --- | --- | --- |
| `wmghg_cmip6_ssp245_1850-2100.nc` | [UoM-CMIP-1-2-0, v20160830](https://doi.org/10.22033/ESGF/input4MIPs.1118) and [UoM-MESSAGE-GLOBIOM-ssp245-1-2-1, v20181127](https://doi.org/10.22033/ESGF/input4MIPs.9866) | Malte Meinshausen and Elisabeth Vogel (2016); Malte Meinshausen and Zebedee R. J. Nicholls (2018), respectively. | Annual global and hemispheric greenhouse-gas mole fractions assembled into the ReadyESM gas file. The production configuration holds the 1993 values. |
| `o3_cmip6_historical_199201-199512.nc` | [UReading-CCMI-1-0, v20160711](https://doi.org/10.22033/ESGF/input4MIPs.1115); latitude and pressure bounds were borrowed from [UReading-CCMI-ssp119-1-1, v20190316](https://doi.org/10.22033/ESGF/input4MIPs.10381) during preparation of the full ozone file. | Michaela Hegglin, Douglas Kinnison, Jean-Francois Lamarque and David Plummer (2016; 2019 for the bounds donor). | Historical ozone was reduced to a zonal mean and subset to January 1992–December 1995. The prepared full file also contains SSP245 ozone after 2014; those scenario months are absent from this bundled subset. |

The source catalogues state the [CC BY-SA 4.0 licence](https://creativecommons.org/licenses/by-sa/4.0/)
and ask users to cite the collection with its version. See the linked catalogue
records for the source terms and recommended citations. This section covers
the bundled gas and ozone files; other input datasets have separate provenance.
