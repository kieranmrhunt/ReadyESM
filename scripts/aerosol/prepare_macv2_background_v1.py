"""Prepare the MACv2 pre-industrial background (PA2 component A) for ReadyESM (Claude Code, 6 Nov 2026). Design v3
(0e5660da) section A; Kieran 13:5x chose the MACv2 pre-industrial background, development use only (not redistributed).

Inputs (pinned by SHA-256 below, retrieved 6 Nov 2026, manifest source_manifest_v1.json 8fa48bd0):
  for_ICON/spec_2005/aeropt_kinne_sw_b14_fin_pre.nc  SW pre-industrial fine mode (gt_p, 1850), 14 records
  for_ICON/spec_2005/aeropt_kinne_sw_b14_coa.nc      SW coarse mode (gt_c, no annual changes), 14 records
  for_ICON/spec_2005/aeropt_kinne_lw_b16_coa.nc      LW: the MACv2 IR total (gt_i, 2005; see design v3), 16 records
  for_ICON/altitude/aeropt_fin_alt_km20.nc, aeropt_coa_alt_km20.nc  fraction of fine/coarse AOD in 20 layers
and the sealed band table band_map_macv2_rrtmgp_v1.json (4af3094c).
Output (schema readiesm_macv2_background_v1): optical properties in RRTMGP band order (record r goes to the RRTMGP
band the table assigns it), latitude ascending, for fine and coarse SW separately (they have different vertical
profiles) and the LW total; the fine and coarse layer fractions with their layer bounds (m above sea level) and the
source surface altitude. The LW uses the coarse profile (declared: the IR is dominated by the coarse mode, and the file
is labelled coarse by ICON). Checks: no NaN, aod >= 0, 0 <= ssa <= 1, 0 <= g < 1, fractions in [0, 1] summing to 1
within 1e-5 in every column, no fraction in layers below the source surface; the reorder is checked by re-reading.
usage: python3 -I prepare_macv2_background_v1.py RAW_DIR BAND_MAP.json OUTPUT.nc
"""
import hashlib
import json
import os
import sys

import netCDF4
import numpy as np

PINNED = {
    "aeropt_kinne_sw_b14_fin_pre.nc": "8a0bf8871b75c6f4",
    "aeropt_kinne_sw_b14_coa.nc": "047f9c79fb8028d7",
    "aeropt_kinne_lw_b16_coa.nc": "134cf36ede1241e5",
    "aeropt_fin_alt_km20.nc": "4d01453fe3b2b0b5",
    "aeropt_coa_alt_km20.nc": "787983ab57ff1f4f",
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def read_optics(path):
    with netCDF4.Dataset(path) as ds:
        ds.set_auto_mask(False)
        lat = ds["lat"][:].astype(float)
        out = {v: ds[v][:].astype(np.float64) for v in ("aod", "ssa", "asy")}   # (record, month, lat, lon)
        lon = ds["lon"][:].astype(float)
    order = np.argsort(lat)
    return {k: v[:, :, order, :] for k, v in out.items()}, lat[order], lon


def check_optics(name, d):
    for v, x in d.items():
        if not np.all(np.isfinite(x)):
            raise ValueError(f"{name} {v} has nonfinite values")
    if np.any(d["aod"] < 0) or np.any(d["ssa"] < 0) or np.any(d["ssa"] > 1) or np.any(d["asy"] < 0) or \
            np.any(d["asy"] >= 1):
        raise ValueError(f"{name} has values outside aod >= 0, 0 <= ssa <= 1, 0 <= g < 1")


def main(raw, band_map, output):
    if os.path.exists(output):
        sys.exit("refusing to overwrite " + output)
    inputs = {}
    for f, prefix in PINNED.items():
        h = sha256(os.path.join(raw, f))
        if not h.startswith(prefix):
            raise ValueError(f"{f} sha256 {h} does not match the retrieval receipt prefix {prefix}")
        inputs[f] = h
    bm = json.load(open(band_map))
    assert bm["pass"], "band map did not pass"

    def order_for(table):
        rows = bm["tables"][table]["rows"]
        target = [r["rrtmgp_band"][0] - 1 for r in rows]                   # record -> RRTMGP band (0-based)
        perm = np.empty(len(target), int)
        for rec, band in enumerate(target):
            perm[band] = rec                                                  # band -> source record
        return perm, [r["wavelength_nm"] for r in rows]

    sw_perm, sw_wave = order_for("sw")
    lw_perm, lw_wave = order_for("lw")
    fine, lat, lon = read_optics(os.path.join(raw, "aeropt_kinne_sw_b14_fin_pre.nc"))
    coarse, lat2, lon2 = read_optics(os.path.join(raw, "aeropt_kinne_sw_b14_coa.nc"))
    lw, lat3, lon3 = read_optics(os.path.join(raw, "aeropt_kinne_lw_b16_coa.nc"))
    assert np.array_equal(lat, lat2) and np.array_equal(lat, lat3) and np.array_equal(lon, lon2) and \
        np.array_equal(lon, lon3)
    for name, d in (("fine", fine), ("coarse", coarse), ("lw", lw)):
        check_optics(name, d)
    fine = {k: v[sw_perm] for k, v in fine.items()}
    coarse = {k: v[sw_perm] for k, v in coarse.items()}
    lw = {k: v[lw_perm] for k, v in lw.items()}
    fractions = {}
    for mode, f, v in (("fine", "aeropt_fin_alt_km20.nc", "zkm_aer_fine_mo"),
                       ("coarse", "aeropt_coa_alt_km20.nc", "zkm_aer_coarse_mo")):
        with netCDF4.Dataset(os.path.join(raw, f)) as ds:
            ds.set_auto_mask(False)
            flat = ds["lat"][:].astype(float)
            x = ds[v][:].astype(np.float64)                                   # (month, layer, lat, lon)
            ztop, zbot = ds["ztop"][:].astype(float), ds["zbot"][:].astype(float)
            asl = ds["asl"][:].astype(np.float64)
        o = np.argsort(flat)
        assert np.array_equal(flat[o], lat)
        x, asl = x[:, :, o, :], asl[o, :]
        s = x.sum(axis=1)
        if not (np.all(np.isfinite(x)) and np.all(x >= 0) and np.all(x <= 1) and np.max(np.abs(s - 1)) <= 1e-5):
            raise ValueError(f"{mode} layer fractions are not finite fractions summing to 1")
        if np.any((x > 0) & (ztop[None, :, None, None] <= asl[None, None, :, :])):
            raise ValueError(f"{mode} has aerosol in layers below the source surface")
        fractions[mode] = x
    with netCDF4.Dataset(output, "w", format="NETCDF4") as ds:
        ds.schema = "readiesm_macv2_background_v1"
        ds.source = ("MACv2 (Kinne 2019, Tellus B 71, doi:10.1080/16000889.2019.1623639), "
                     "ftp://ftp-projects.mpimet.mpg.de/aerocom/climatology/MACv2_2018/for_ICON, retrieved 2026-11-06")
        ds.licence = ("no data licence published; development use only, not to be redistributed in release or GMD "
                      "archives (Kieran Hunt, 2026-11-06)")
        ds.input_sha256 = json.dumps(inputs)
        ds.band_map = "band_map_macv2_rrtmgp_v1.json sha256 " + sha256(band_map)
        ds.band_rule = ("MACv2 band-representative monochromatic values; record r is placed in the RRTMGP band whose "
                        "limits contain its wavelength (one to one, see the band map); not band integrals")
        ds.composition = ("SW: pre-industrial fine (gt_p, 1850) and coarse (gt_c), kept separate; LW: the MACv2 IR "
                          "total (gt_i, 2005), which includes a small 2005 anthropogenic fine-mode term")
        ds.vertical = ("fine and coarse AOD fractions in 20 layers (m above sea level); the LW uses the coarse profile")
        ds.preparation_script_sha256 = sha256(os.path.abspath(__file__))
        for name, n in (("month", 12), ("lat", len(lat)), ("lon", len(lon)), ("sw_band", 14), ("lw_band", 16),
                        ("layer", 20)):
            ds.createDimension(name, n)
        ds.createVariable("lat", "f8", ("lat",), fill_value=False)[:] = lat
        ds.createVariable("lon", "f8", ("lon",), fill_value=False)[:] = lon
        v = ds.createVariable("sw_wavelength_nm", "i4", ("sw_band",), fill_value=False)
        v[:] = [sw_wave[r] for r in sw_perm]
        v = ds.createVariable("lw_wavelength_nm", "i4", ("lw_band",), fill_value=False)
        v[:] = [lw_wave[r] for r in lw_perm]
        ds.createVariable("layer_bottom_m", "f8", ("layer",), fill_value=False)[:] = zbot
        ds.createVariable("layer_top_m", "f8", ("layer",), fill_value=False)[:] = ztop
        ds.createVariable("surface_altitude_m", "f8", ("lat", "lon"), fill_value=False)[:] = asl
        for mode, d in (("fine_sw", fine), ("coarse_sw", coarse), ("lw", lw)):
            band = "lw_band" if mode == "lw" else "sw_band"
            for k, src in (("tau", "aod"), ("ssa", "ssa"), ("g", "asy")):
                ds.createVariable(f"{k}_{mode}", "f4", (band, "month", "lat", "lon"), fill_value=False)[:] = d[src]
        for mode in ("fine", "coarse"):
            ds.createVariable(f"layer_fraction_{mode}", "f4", ("month", "layer", "lat", "lon"),
                              fill_value=False)[:] = fractions[mode]
    # Re-read check of the band reorder: RRTMGP band b of the output equals source record perm[b].
    with netCDF4.Dataset(output) as ds:
        ds.set_auto_mask(False)
        src, _, _ = read_optics(os.path.join(raw, "aeropt_kinne_sw_b14_fin_pre.nc"))
        for b in range(14):
            assert np.array_equal(ds["tau_fine_sw"][b].astype(np.float32), src["aod"][sw_perm[b]].astype(np.float32))
    i550 = [sw_wave[r] for r in sw_perm].index(550)
    receipt = {"inputs": inputs, "output_sha256": sha256(output),
               "sw_band_order_wavelength_nm": [sw_wave[r] for r in sw_perm],
               "lw_band_order_wavelength_nm": [lw_wave[r] for r in lw_perm],
               "global_mean_aod550_annual_fine_pre": float((fine["aod"][i550].mean(axis=0) * np.cos(np.deg2rad(lat))[:, None]).sum()
                                                           / (np.cos(np.deg2rad(lat)).sum() * len(lon))),
               "global_mean_aod550_annual_coarse": float((coarse["aod"][i550].mean(axis=0) * np.cos(np.deg2rad(lat))[:, None]).sum()
                                                         / (np.cos(np.deg2rad(lat)).sum() * len(lon)))}
    print(json.dumps(receipt, indent=1))
    print("MACV2_BACKGROUND_PREPARE_V1_DONE")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
