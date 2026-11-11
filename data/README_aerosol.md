# Prescribed aerosol input data

ReadyESM's prescribed aerosol option (`aerosol_scheme: cmip_prescribed`) reads three files. They are not in the
repository. This note tells you how to download the source data and make the three files. All scripts are in
`scripts/aerosol/`. Run them with Python 3.11 and the packages numpy, scipy, netCDF4 and mpmath.

| File the model reads | What it holds | How you get it |
|---|---|---|
| `data/macv2_2018_background_prepared_v1/macv2_background_v1.nc` | MACv2 natural background aerosol | Download, then prepare |
| `data/macv2sp_cmip6_v1_0_raw/MACv2.0-SP_v1.nc` | MACv2-SP anthropogenic plumes | Download only |
| `data/uoexeter_cmip_2_2_1_prepared_v7/uoexeter_strat_bands_199212_199501_v7.nc` | Stratospheric (volcanic) aerosol | Download, then prepare |

The prepared stratospheric file covers December 1992 to January 1995. A run that goes outside these months stops
with an error. It does not extrapolate.

## 1. MACv2 background

Licence: MACv2 has no published data licence. Do not redistribute the raw or prepared MACv2 files. Cite Kinne (2019).

Download these files from `ftp://ftp-projects.mpimet.mpg.de/aerocom/climatology/MACv2_2018/`.
The path column below identifies each file on the server. Put all five `.nc`
files directly in `data/macv2_2018_background_raw/` (without the `for_ICON/`
subdirectories), as required by the preparer.

| File (path under MACv2_2018/) | SHA-256 |
|---|---|
| for_ICON/spec_2005/aeropt_kinne_sw_b14_fin_pre.nc | 8a0bf8871b75c6f4d27b45c839af54418863fa32533aad073061ecee0f884f26 |
| for_ICON/spec_2005/aeropt_kinne_sw_b14_coa.nc | 047f9c79fb8028d70a7641c94c12d006d07d163c21b8a7647ac80cf0b2dd4aa3 |
| for_ICON/spec_2005/aeropt_kinne_lw_b16_coa.nc | 134cf36ede1241e5c6c5c5c0317b02abcd1e38287b4b796f92994e8abdbc9c82 |
| for_ICON/altitude/aeropt_fin_alt_km20.nc | 4d01453fe3b2b0b57af17e88cf9b3d55b45c5a5ce31ee6a81653f2303cf987f8 |
| for_ICON/altitude/aeropt_coa_alt_km20.nc | 787983ab57ff1f4f15a32de23e3bcc654fe5ef7a6aab959f8c6690c87f6811e2 |

Check the complete raw-file hashes, then run:

```bash
(cd data/macv2_2018_background_raw && sha256sum --strict -c ../../scripts/aerosol/macv2_raw.sha256)
mkdir -p data/macv2_2018_background_prepared_v1
python3 -I scripts/aerosol/prepare_macv2_background_v1.py data/macv2_2018_background_raw \
    scripts/aerosol/band_map_macv2_rrtmgp_v1.json \
    data/macv2_2018_background_prepared_v1/macv2_background_v1.nc
```

The script checks the first 16 hex characters of each raw file's SHA-256 before
it uses the file. The `sha256sum` command above checks all 64. The band map
`band_map_macv2_rrtmgp_v1.json` assigns each MACv2 spectral record to an RRTMGP band. It was made by
`band_map_macv2_rrtmgp_v1.py` and is shipped as is.

## 2. MACv2-SP plumes

Licence: CMIP6 input4MIPs terms. Cite Stevens et al. (2017).

Download `MACv2.0-SP_v1.nc` into `data/macv2sp_cmip6_v1_0_raw/` from ESGF:
`input4MIPs/CMIP6/RFMIP/MPI-M/MPI-M-MACv2-SP-1-0/atmos/yr/multiple/none/v20170201/`.
Its SHA-256 is `b99f13727d42ed24ba71713064cbb40b3a0211ced0e725c56c794bce394bc5f8`. The model reads it directly.

## 3. Stratospheric aerosol (UOEXETER-CMIP-2-2-1)

Licence: CC BY 4.0. Contact t.aubry@exeter.ac.uk.

Download the monthly extinction, single-scattering albedo and asymmetry files from ESGF:
`input4MIPs/CMIP7/CMIP/uoexeter/UOEXETER-CMIP-2-2-1/atmos/mon/{ext,ssa,asy}/gnz/v20250521/`.
Each file is about 2.7 GB. Put them in one directory.

| File | SHA-256 |
|---|---|
| ext_input4MIPs_aerosolProperties_CMIP_UOEXETER-CMIP-2-2-1_gnz_175001-202312.nc | 2a3980ff420487f94816f53f3fa2b6518bd28104d30ad66d55dd5f24baabc422 |
| ssa_input4MIPs_aerosolProperties_CMIP_UOEXETER-CMIP-2-2-1_gnz_175001-202312.nc | c7406b8651a629eaec6cf9d84614e6b838ba5ff0d0b948dc60bdc7000b6f5abe |
| asy_input4MIPs_aerosolProperties_CMIP_UOEXETER-CMIP-2-2-1_gnz_175001-202312.nc | 1f3f3f984b5c788cabe8df682f7ff31205688c4733575ee198339976e9154f60 |

The preparation also needs two more inputs:
- The CMIP6 solar spectrum `solarforcing-ref-mon_input4MIPs_solar_CMIP_SOLARIS-HEPPA-3-2_gn_185001-229912.nc`, from
  `input4MIPs/CMIP6/CMIP/SOLARIS-HEPPA/SOLARIS-HEPPA-3-2/atmos/mon/solarforcing-ref-mon/gn/v20170103/`.
  SHA-256 `4b844a3324b74be16faf7272b64ed7853ea2c606c4db841a318366ad223ebcac`.
- The RRTMGP data directory `rrtmgp-data-1.9`, which the Julia RRTMGP package installs as an artifact. The script
  uses `rrtmgp-gas-sw-g224.nc` (SHA-256 `584f1dd41ea9fc07d4ee3754eb1dafbd46ad3161cd6fd20fa06b6922b6f0702e`) and
  `rrtmgp-gas-lw-g256.nc` (SHA-256 `4048360199d1917ed8f2ccaae2ec097d0f990da3bbad9830337b739b4fa01be7`).

Make the output directory first (`mkdir -p data/uoexeter_cmip_2_2_1_prepared_v7`), then run:

```bash
python3 -I scripts/aerosol/prepare_uoexeter_v7.py UOEXETER_DIR SOLAR_FILE RRTMGP_DATA_DIR \
    data/uoexeter_cmip_2_2_1_prepared_v7/uoexeter_strat_bands_199212_199501_v7.nc \
    data/uoexeter_cmip_2_2_1_prepared_v7/prepare_receipt_v7.json
```

The script checks every input hash and writes a receipt. To check a prepared file, run `verify_prepared_v8.py`. It
prepares the file again in a new directory and compares the two:

```bash
python3 -I scripts/aerosol/verify_prepared_v8.py UOEXETER_DIR SOLAR_FILE RRTMGP_DATA_DIR \
    data/uoexeter_cmip_2_2_1_prepared_v7/uoexeter_strat_bands_199212_199501_v7.nc \
    data/uoexeter_cmip_2_2_1_prepared_v7/prepare_receipt_v7.json scripts/aerosol/prepare_uoexeter_v7.py WORK_DIR
```

The verifier also requires the exact Python and library versions that made the reference file. These are Python
3.11.15, numpy 2.4.3, scipy 1.17.1, netCDF4 1.7.4 (netCDF-C 4.10.0, HDF5 1.14.6) and mpmath 1.4.1. With other
versions the preparation still runs, but the verifier reports a runtime mismatch.

Reference result (job 66237764, 9 Nov 2026): the prepared file has SHA-256
`30e4e7383d6511c85c9aecdf84505bdb9f56d1aba717055aa0097c160daf607a`, and a fresh preparation reproduced it exactly.

### Frozen six-month evaluation input

`config/production_six_month.yml` retains the UOEXETER v6 file used by the
frozen v0.1.20 six-month evaluation. Prepare it from the same raw, solar and
RRTMGP inputs listed above with the included v6 script:

```bash
mkdir -p data/uoexeter_cmip_2_2_1_prepared_v6
python3 -I scripts/aerosol/prepare_uoexeter_v6.py UOEXETER_DIR SOLAR_FILE RRTMGP_DATA_DIR \
    data/uoexeter_cmip_2_2_1_prepared_v6/uoexeter_strat_bands_199212_199401_v6.nc \
    data/uoexeter_cmip_2_2_1_prepared_v6/prepare_receipt_v6.json
python3 -I scripts/aerosol/verify_prepared_v7.py UOEXETER_DIR SOLAR_FILE RRTMGP_DATA_DIR \
    data/uoexeter_cmip_2_2_1_prepared_v6/uoexeter_strat_bands_199212_199401_v6.nc \
    data/uoexeter_cmip_2_2_1_prepared_v6/prepare_receipt_v6.json \
    scripts/aerosol/prepare_uoexeter_v6.py WORK_DIR
```

The frozen file SHA-256 is
`83b0a9d647f41a540331b8ab2c0f64c80fb46fd4aff1ec82457866909366af37`.
The v6 receipt SHA-256 is
`da6e7ffc7a70ebf24c7263cb073f4030efd8ae47032bf961b949c56d4b7715b1`.
The v7 file is the production default for the 730-day horizon.

## References

- Kinne, S. (2019). The MACv2 aerosol climatology. Tellus B 71. doi:10.1080/16000889.2019.1623639
- Stevens, B. et al. (2017). MACv2-SP: a parameterization of anthropogenic aerosol optical properties and an
  associated Twomey effect for use in CMIP6. GMD 10, 433-452. doi:10.5194/gmd-10-433-2017
- UOEXETER-CMIP-2-2-1 stratospheric aerosol optical properties, CMIP7 input4MIPs, version v20250521.
