"""Band-mapping table for the MACv2 representative wavelengths onto RRTMGP bands (Claude Code, 6 Nov 2026; Paul's
PA2 v2 condition 2). Reads the representative wavelengths from each MACv2 file's own NCO history (the ncecat order of
gt_?_NNNNNnm.nc, which is the record order), the RRTMGP band limits from the RRTMGP-data v1.9 gas files that RRTMGP
0.22.2 loads (get_artifact_path), and writes, for every source record: record index, wavelength (nm), wavenumber
(cm-1), the RRTMGP band containing it and that band's limits. It also lists the RRTMG band limits (AER RRTMG_SW bands
16-29, RRTMG_LW bands 1-16) as published in the RRTMG documentation (stated, not read from a file) and checks they
equal the RRTMGP limits. PASS requires: every record in exactly one RRTMGP band, every band used exactly once.
usage: python3 -I band_map_macv2_rrtmgp_v1.py MACV2_DIR RRTMGP_DATA_DIR OUTPUT.json"""
import json
import re
import sys

import netCDF4

macv2, rrtmgp, output = sys.argv[1:4]
# RRTMG band limits (cm-1), as in the AER RRTMG_SW/LW documentation: SW bands 16-29 (band 29 is 820-2600 in RRTMG_SW),
# LW bands 1-16. Listed for the record; the check below compares them with RRTMGP's.
RRTMG_SW = [(2600, 3250), (3250, 4000), (4000, 4650), (4650, 5150), (5150, 6150), (6150, 7700), (7700, 8050),
            (8050, 12850), (12850, 16000), (16000, 22650), (22650, 29000), (29000, 38000), (38000, 50000), (820, 2600)]
RRTMG_LW = [(10, 350), (350, 500), (500, 630), (630, 700), (700, 820), (820, 980), (980, 1080), (1080, 1180),
            (1180, 1390), (1390, 1480), (1480, 1800), (1800, 2080), (2080, 2250), (2250, 2380), (2380, 2600),
            (2600, 3250)]
report = {"rrtmgp_data": rrtmgp, "tables": {}}
ok = True
for spectrum, mfile, gfile, rrtmg in (("sw", "aeropt_kinne_sw_b14_fin_pre.nc", "rrtmgp-gas-sw-g224.nc", RRTMG_SW),
                                      ("sw_coarse", "aeropt_kinne_sw_b14_coa.nc", "rrtmgp-gas-sw-g224.nc", RRTMG_SW),
                                      ("lw", "aeropt_kinne_lw_b16_coa.nc", "rrtmgp-gas-lw-g256.nc", RRTMG_LW)):
    with netCDF4.Dataset(f"{macv2}/{mfile}") as ds:
        history = ds.history
        nrec = len(ds.dimensions["record"])
    waves = [int(w) for w in re.findall(r"gt_[a-z]_(\d{5})nm\.nc", history.split("\n")[0])]
    with netCDF4.Dataset(f"{rrtmgp}/{gfile}") as ds:
        limits = [(float(a), float(b)) for a, b in ds["bnd_limits_wavenumber"][:]]
    rows, used = [], []
    for i, w in enumerate(waves):
        k = 1.0e7 / w
        bands = [j for j, (a, b) in enumerate(limits) if a <= k < b]
        rows.append({"record": i + 1, "wavelength_nm": w, "wavenumber_cm1": round(k, 3),
                     "rrtmgp_band": [j + 1 for j in bands], "rrtmgp_band_limits_cm1": [limits[j] for j in bands]})
        used += bands
    one_each = len(waves) == nrec == len(limits) and all(len(r["rrtmgp_band"]) == 1 for r in rows) and \
        sorted(used) == list(range(len(limits)))
    rrtmg_sorted = sorted(rrtmg)
    same = [tuple(map(float, x)) for x in rrtmg_sorted] == sorted(limits)
    differ = [(a, b) for a, b in zip(rrtmg_sorted, sorted(limits)) if tuple(map(float, a)) != b]
    report["tables"][spectrum] = {"macv2_file": mfile, "records": nrec, "rrtmgp_file": gfile, "rrtmgp_limits": limits,
                                  "rows": rows, "one_to_one": one_each, "rrtmg_limits_equal_rrtmgp": same,
                                  "rrtmg_vs_rrtmgp_differences": differ}
    ok &= one_each
    print(f"== {spectrum} {mfile}: records={nrec} one_to_one={one_each} rrtmg_limits_equal_rrtmgp={same} differ={differ}")
    for r in rows:
        print(f"   rec {r['record']:2d} {r['wavelength_nm']:6d} nm {r['wavenumber_cm1']:10.3f} cm-1 -> band "
              f"{r['rrtmgp_band']} {r['rrtmgp_band_limits_cm1']}")
report["pass"] = ok
json.dump(report, open(output, "w"), indent=1)
print("BAND_MAP_V1_PASS" if ok else "BAND_MAP_V1_FAIL")
