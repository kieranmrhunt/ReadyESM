"""Fully realistic copy of the ERA5 seasonal land albedo cycle, for the V3 test arm (Claude, 2026-10-30).

usage: python3 -I -B make_realistic.py ALBEDO_CYCLE ERA5_FAL_1994_2013 OUTPUT_DIR

Kieran, 30 Oct 19:14 UTC: "use fully realistic albedo". ALBEDO_CYCLE is the opt-in option's input
`era5_monthly_land_cycle_1993.nc` (sha256 b80d7e01...): ERA5 snow-free background albedos in four bands (visible and
near-infrared, direct and diffuse) on a 1-degree grid. Ice sheets and snow are missing from it, so the Greenland interior
is about 0.19 there. ERA5_FAL_1994_2013 is the held-out ERA5 monthly forecast albedo `fal` (sha256 65db76e8..., 0.5
degree, 1994 to 2013), ERA5's realised surface albedo with snow and ice included. The 1994-2013 mean by calendar
month is used, not 1993, so the 1993 test is not fed its own year (the held-out rule of the 22 October land map).

For each month and 1-degree cell:
- F is the mean of the 1994-2013 `fal` climatology over the land pixels (ERA5 land-sea mask 0.5 or more) of the
  0.5-degree box centred on the cell, with half weight for pixels on the box edge and cos(latitude) weights.
- B is the mean of the four background bands.
- Each band becomes clip(band + d, 0, 1), with one offset d per cell and month chosen by bisection so that the mean
  of the four clipped bands equals F (Codex, 19:17: a plain shift F - B undershoots where a band clips). The broadband
  value then equals ERA5's realised albedo, and the background's spectral and direct/diffuse differences are kept where
  no band clips. A cell with no land pixel in its box keeps the background.
v2 of albedo_realistic_claude, which used the plain shift.
Limits: a uniform shift does not give snow its real visible/near-infrared split, and `fal` is a diffuse all-sky value.
Sea ice is not in this file (the model's ocean scheme sets it). This is a labelled test-arm input.
"""
import hashlib
import json
import shutil
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, num2date

FIELDS = ("aluvp", "aluvd", "alnip", "alnid")
ALBEDO_SHA = "b80d7e01a4390340794030f1f924e5583d813fb163a1f9d875ae749091c13f60"
FAL_SHA = "65db76e869a0605dede1dd1ee6a2a5ac64cea38a759f63c8a9845eca6167f65a"
NAME = "era5_monthly_land_cycle_1993_realistic_v2.nc"


def sha256(path):
    with open(path, "rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def edge_weights(fine, centre, half, periodic):
    d = fine - centre
    if periodic:
        d = (d + 180.0) % 360.0 - 180.0
    w = np.where(np.abs(d) < half - 1e-6, 1.0, 0.0)
    return np.where(np.abs(np.abs(d) - half) <= 1e-6, 0.5, w)


def climatology(path):
    with Dataset(path) as nc:
        clock = nc["valid_time"]
        stamps = num2date(clock[:], clock.units, calendar=getattr(clock, "calendar", "standard"))
        if [(s.year, s.month) for s in stamps] != [(y, m) for y in range(1994, 2014) for m in range(1, 13)]:
            raise ValueError("fal source is not the complete 1994-2013 monthly sequence")
        lat = np.asarray(nc["latitude"][:], float)
        lon = np.asarray(nc["longitude"][:], float) % 360
        fal = np.ma.filled(np.ma.asarray(nc["fal"][:], float), np.nan)
        lsm = np.ma.filled(np.ma.asarray(nc["lsm"][:], float), np.nan)
    if not (np.allclose(np.diff(lat), -0.5) and np.allclose(np.diff(lon), 0.5)):
        raise ValueError("fal source is not the 0.5-degree grid")
    clim = np.stack([np.nanmean(fal[m::12], axis=0) for m in range(12)])
    land = np.stack([np.nanmean(lsm[m::12], axis=0) >= 0.5 for m in range(12)])
    valid = land & np.isfinite(clim)
    if np.any(valid & ((clim < 0) | (clim > 1))):
        raise ValueError("fal climatology outside 0..1 on land")
    return lat, lon, clim, valid


def box_mean(flat, flon, values, valid, lat, lon):
    out = np.full((lat.size, lon.size), np.nan)
    coslat = np.cos(np.deg2rad(flat))
    lon_w = [edge_weights(flon, x, 0.5, True) for x in lon]
    lon_i = [np.flatnonzero(w) for w in lon_w]
    for i, y in enumerate(lat):
        wr = edge_weights(flat, y, 0.5, False) * coslat
        rows = np.flatnonzero(wr)
        for j in range(lon.size):
            cols = lon_i[j]
            w = wr[rows][:, None] * lon_w[j][cols][None, :] * valid[np.ix_(rows, cols)]
            if w.sum() > 0:
                out[i, j] = float((w * np.nan_to_num(values[np.ix_(rows, cols)])).sum() / w.sum())
    return out


def main(albedo_path, fal_path, output_dir):
    albedo_path, fal_path, output_dir = Path(albedo_path), Path(fal_path), Path(output_dir)
    if sha256(albedo_path) != ALBEDO_SHA or sha256(fal_path) != FAL_SHA:
        raise ValueError("input hash differs from the recorded source")
    flat, flon, clim, valid = climatology(fal_path)
    with Dataset(albedo_path) as nc:
        lat = np.asarray(nc["latitude"][:], float)
        lon = np.asarray(nc["longitude"][:], float)
        old = {name: np.ma.filled(np.ma.asarray(nc[name][:], float), np.nan) for name in FIELDS}
    if not (np.allclose(np.diff(lat), -1.0) and np.allclose(np.diff(lon), 1.0)):
        raise ValueError("albedo file is not the 1-degree grid")
    target = np.stack([box_mean(flat, flon, clim[m], valid[m], lat, lon) for m in range(12)])
    background = sum(old[name] for name in FIELDS) / len(FIELDS)
    bands = np.stack([old[name] for name in FIELDS])
    if not np.all(np.isfinite(bands)) or np.any((bands < 0) | (bands > 1)):
        raise ValueError("background albedo outside 0..1")
    # mean(clip(bands + d)) rises monotonically from 0 at d = -1 to 1 at d = 1, so bisection finds d for any F in 0..1.
    lo, hi = np.full(background.shape, -1.0), np.full(background.shape, 1.0)
    goal = np.where(np.isfinite(target), target, background)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        too_low = np.clip(bands + mid[None], 0.0, 1.0).mean(axis=0) < goal
        lo, hi = np.where(too_low, mid, lo), np.where(too_low, hi, mid)
    shift = np.where(np.isfinite(target), 0.5 * (lo + hi), 0.0)
    residual = np.abs(np.clip(bands + shift[None], 0.0, 1.0).mean(axis=0) - goal)
    if np.nanmax(np.where(np.isfinite(target), residual, 0.0)) > 1e-6:
        raise ValueError("broadband target not met")
    output_dir.mkdir(parents=False, exist_ok=False)
    out = output_dir / NAME
    shutil.copyfile(albedo_path, out)
    L, O = np.meshgrid(lat, lon % 360, indexing="ij")
    boxes = {"greenland_65_80N_305_325E": (L > 65) & (L < 80) & (O > 305) & (O < 325),
             "land_60_90N": L >= 60, "land_50_60N": (L >= 50) & (L < 60), "antarctica_south_of_70S": L < -70}
    stats = {}
    with Dataset(out, "a") as nc:
        lsm = np.ma.filled(np.ma.asarray(nc["lsm"][:], float), np.nan)
        new = {}
        for name in FIELDS:
            values = np.where(np.isfinite(old[name]), np.clip(old[name] + shift, 0.0, 1.0), old[name])
            nc[name][:] = values.astype(np.float32)
            new[name] = values
        broadband = sum(new[name] for name in FIELDS) / len(FIELDS)
        for label, box in boxes.items():
            for m in (0, 6):
                w = np.cos(np.deg2rad(L)) * box * (lsm[m] > 0.5)
                stats[f"{label}_month{m + 1}"] = {
                    "before": float(np.nansum(background[m] * w) / w.sum()),
                    "after": float(np.nansum(broadband[m] * w) / w.sum()),
                    "era5_fal_target": float(np.nansum(np.nan_to_num(target[m]) * w * np.isfinite(target[m])) /
                                             (w * np.isfinite(target[m])).sum())}
        clipped = int(sum(((old[n] + shift) < 0).sum() + ((old[n] + shift) > 1).sum() for n in FIELDS))
        nc.setncattr("readiesm_realistic_albedo",
                     "Test-arm input (Claude, 30 Oct 2026, Kieran: use fully realistic albedo). Land only, sea ice not "
                     "included. Each of aluvp, aluvd, alnip, alnid becomes clip(band + d, 0, 1), with d per cell and "
                     "month such that the four-band mean equals F, the 1-degree land box mean of the ERA5 fal 1994-2013 "
                     "monthly climatology. Not a production ancillary.")
        nc.setncattr("readiesm_realistic_albedo_source_sha256", ALBEDO_SHA)
        nc.setncattr("readiesm_realistic_albedo_fal_sha256", FAL_SHA)
        nc.setncattr("readiesm_realistic_albedo_script_sha256", sha256(__file__))
    receipt = {"schema": "readiesm_albedo_realistic_v2", "path": str(out), "sha256": sha256(out),
               "albedo_source_sha256": ALBEDO_SHA, "fal_source_sha256": FAL_SHA, "script_sha256": sha256(__file__),
               "cells_without_land_target": int((~np.isfinite(target)).sum()), "band_values_clipped": clipped,
               "max_broadband_residual": float(np.nanmax(np.where(np.isfinite(target), residual, 0.0))),
               "broadband_land_means": stats}
    with open(output_dir / "receipt.json", "x") as stream:
        json.dump(receipt, stream, indent=1)
    print("REALISTIC_ALBEDO_WRITTEN", out, receipt["sha256"], "clipped", clipped)
    for key, value in stats.items():
        print(f"  {key:40s} before {value['before']:.3f} after {value['after']:.3f} ERA5 fal {value['era5_fal_target']:.3f}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
