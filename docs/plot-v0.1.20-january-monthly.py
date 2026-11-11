#!/usr/bin/env python3
"""Render the reverted v0.1.20 January monthly atmosphere/radiation README figure."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm
import netCDF4
import numpy as np
from scipy.spatial import cKDTree

COMMIT = "fba77cd3ffc42f0f9227f9cd2ca24840bb8ad396"
TREE = "b735b45610c67e24e0d99cf9d31b4e2f894f4d2c"
START = "1993-01-01T00:00:00"
STOP = "1993-02-01T00:00:00"

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def read(ds, name, units=None):
    variable = ds.variables[name]
    if units is not None and getattr(variable, "units", None) != units:
        raise ValueError(f"{name}: unexpected units {getattr(variable, 'units', None)}")
    return np.ma.filled(variable[:], np.nan).astype(float)

def check_window(ds, start, stop):
    if getattr(ds, "start_time_UTC", getattr(ds, "window_start", None)) != start:
        raise ValueError("January start time mismatch")
    if getattr(ds, "stop_time_UTC", getattr(ds, "window_end", None)) != stop:
        raise ValueError("January stop time mismatch")

def xyz(lon, lat):
    lon, lat = np.deg2rad(lon), np.deg2rad(lat)
    return np.column_stack((np.cos(lat) * np.cos(lon),
                            np.cos(lat) * np.sin(lon), np.sin(lat)))

def raster_indices(lon, lat):
    x, y = np.meshgrid(np.arange(-179.75, 180, .5),
                       np.arange(-89.75, 90, .5))
    indices = cKDTree(xyz(lon, lat)).query(xyz(x.ravel(), y.ravel()))[1]
    return indices.reshape(x.shape)

def render(run_dir, output):
    run_dir, output = map(Path, (run_dir, output))
    month_dir = run_dir / "model_output" / "months" / "1993-01"
    receipt = run_dir / "receipt"
    status = run_dir / "status"
    verdict = run_dir / "verdict"
    result_manifest = run_dir / "result.sha256"
    state = dict(line.split("=", 1) for line in status.read_text().splitlines() if "=" in line)
    if state.get("phase") != "complete" or state.get("exit_code") != "0":
        raise ValueError("six-month run is not terminal and successful")
    if verdict.read_text().strip() != "V020R_PARTC_SIX_MONTH_PASS":
        raise ValueError("six-month verdict missing")
    subprocess.run(["sha256sum", "--strict", "--quiet", "-c", result_manifest.name],
                   cwd=run_dir, check=True)
    paths = {name: month_dir / name for name in ("monthly.nc", "wind.nc", "toa_monthly.nc")}
    png, metadata_path = Path(f"{output}.png"), Path(f"{output}.json")
    if png.exists() or metadata_path.exists():
        raise FileExistsError("Choose a new output prefix")
    lines = dict(line.split("=", 1) for line in receipt.read_text().splitlines() if "=" in line)
    if lines.get("mode") != "full" or lines.get("source_commit") != COMMIT or lines.get("source_tree") != TREE:
        raise ValueError("figure must use the reverted fba six-month source receipt")
    with (netCDF4.Dataset(paths["monthly.nc"]) as monthly,
          netCDF4.Dataset(paths["wind.nc"]) as wind,
          netCDF4.Dataset(paths["toa_monthly.nc"]) as toa):
        if (monthly.schema, wind.schema, toa.schema) != (
            "c02_monthly_capture_v1", "monthly_wind_capture_v1", "toa_monthly_v1"):
            raise ValueError("unexpected monthly NetCDF schemas")
        for ds in (monthly, wind, toa):
            check_window(ds, START, STOP)
        if (monthly.state_samples, monthly.flux_samples, wind.state_samples, toa.samples) != (
            744, 2976, 744, 744):
            raise ValueError("incomplete January monthly sampling")
        if (monthly.pressure_interpolation != "Linear in log pressure, no extrapolation. Complete hourly coverage required."
            or wind.pressure_interpolation != monthly.pressure_interpolation):
            raise ValueError("850-hPa interpolation contract changed")
        lon, lat = read(monthly, "longitude", "degrees_east"), read(monthly, "latitude", "degrees_north")
        weights = read(monthly, "area_weight", "1")
        if lon.shape != (4608,) or lat.shape != lon.shape or weights.shape != lon.shape:
            raise ValueError("unexpected atmosphere grid")
        if not (np.all(np.isfinite(lon)) and np.all(np.isfinite(lat)) and
                np.all(np.isfinite(weights)) and np.all(weights > 0) and
                np.isclose(weights.sum(), 1, atol=1e-5)):
            raise ValueError("invalid atmosphere coordinates or weights")
        for ds in (wind, toa):
            if not (np.array_equal(read(ds, "longitude"), lon) and
                    np.array_equal(read(ds, "latitude"), lat) and
                    np.allclose(read(ds, "area_weight"), weights, rtol=0, atol=1e-8)):
                raise ValueError("monthly grids differ")
        pressure = read(monthly, "pressure", "Pa")
        if not np.array_equal(pressure, read(wind, "pressure", "Pa")):
            raise ValueError("temperature and wind pressure levels differ")
        where850 = np.flatnonzero(pressure == 85000)
        if where850.size != 1:
            raise ValueError("missing unique 850-hPa level")
        p = int(where850[0])
        tvalid = read(monthly, "temperature_valid")[p]
        wvalid = read(wind, "valid")[p]
        tcount = read(monthly, "temperature_sample_count")[p]
        wcount = read(wind, "sample_count")[p]
        if not (np.isin(tvalid, (0, 1)).all() and np.isin(wvalid, (0, 1)).all() and
                np.array_equal(tvalid == 1, tcount == 744) and
                np.array_equal(wvalid == 1, wcount == 744)):
            raise ValueError("850-hPa valid/sample-count contract changed")
        t850 = read(monthly, "temperature", "K")[p] - 273.15
        u850 = read(wind, "u", "m s-1")[p]
        v850 = read(wind, "v", "m s-1")[p]
        speed850 = np.hypot(u850, v850)
        t850[tvalid != 1] = np.nan
        speed850[wvalid != 1] = np.nan
        cw = read(monthly, "column_water", "kg m-2")
        rain = read(monthly, "precipitation_mm_day", "mm day-1")
        olr = read(toa, "outgoing_longwave", "W m-2")
        lwcre = read(toa, "clear_outgoing_longwave", "W m-2") - olr
        for name, field in (("850-hPa temperature", t850[tvalid == 1]),
                            ("850-hPa wind", speed850[wvalid == 1]),
                            ("column water", cw), ("precipitation", rain),
                            ("outgoing longwave", olr), ("longwave cloud effect", lwcre)):
            if field.size == 0 or not np.all(np.isfinite(field)):
                raise ValueError(f"non-finite or empty {name}")
        if np.any(rain < 0) or np.any(cw < 0) or np.any(speed850[wvalid == 1] < 0):
            raise ValueError("negative water or wind speed")

        panels = [
            ("850-hPa air temperature", t850, "°C", "coolwarm", Normalize(-35, 35)),
            ("850-hPa wind speed", speed850, "m s⁻¹", "viridis", Normalize(0, 40)),
            ("Column water vapour", cw, "kg m⁻²", "YlGnBu", Normalize(0, 70)),
            ("Precipitation", rain, "mm day⁻¹", "Blues", LogNorm(.05, 40)),
            ("Outgoing longwave radiation", olr, "W m⁻²", "magma", Normalize(120, 350)),
            ("Longwave cloud effect", lwcre, "W m⁻²; clear − all", "RdBu_r", TwoSlopeNorm(vmin=-30, vcenter=0, vmax=100)),
        ]
        indices = raster_indices(lon, lat)
        plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                             "axes.spines.top": False, "axes.spines.right": False})
        fig, axes = plt.subplots(3, 2, figsize=(14, 11), layout="constrained")
        fig.suptitle("ReadyESM v0.1.20 · January 1993 monthly mean", fontsize=18, fontweight="bold")
        metadata = {
            "schema": "v020_reverted_readme_monthly_figure_v1",
            "source_commit": COMMIT, "source_tree": TREE,
            "window_start": START, "window_end": STOP,
            "duration_days": 31, "observed_start_count": 1,
            "receipt_sha256": sha256(receipt),
            "status_sha256": sha256(status),
            "verdict_sha256": sha256(verdict),
            "run_result_manifest_sha256": sha256(result_manifest),
            "input_sha256": {name: sha256(path) for name, path in paths.items()},
            "regrid": "nearest native atmosphere point on 0.5-degree equirectangular display",
            "pressure_rule": "850 hPa shown only with complete valid hourly coverage; no below-ground extrapolation",
            "panels": {},
        }
        for i, (ax, (title, field, units, cmap_name, norm)) in enumerate(zip(axes.flat, panels)):
            cmap = plt.get_cmap(cmap_name).copy()
            cmap.set_bad("#c5c8cc")
            display = field[indices]
            if title == "Precipitation":
                display = np.where(display > 0, display, np.nan)
            image = ax.imshow(display, origin="lower", extent=(-180, 180, -90, 90),
                              cmap=cmap, norm=norm, interpolation="nearest", aspect="equal")
            ax.set_title(f"{chr(65 + i)}  {title}", loc="left", fontweight="bold")
            ax.set_xticks([-180, -90, 0, 90, 180], ["180°", "90°W", "0°", "90°E", "180°"])
            ax.set_yticks([-60, 0, 60], ["60°S", "0°", "60°N"])
            cb = fig.colorbar(image, ax=ax, orientation="horizontal", pad=.035, fraction=.08)
            cb.set_label(units)
            finite = field[np.isfinite(field)]
            metadata["panels"][title] = {"units": units, "native_valid_points": int(finite.size),
                                         "native_min": float(finite.min()), "native_max": float(finite.max())}
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(png, dpi=160, facecolor="white")
        plt.close(fig)
        metadata["figure_sha256"] = sha256(png)
        metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(f"V020_REVERTED_MONTHLY_FIGURE_PASS {png}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir")
    parser.add_argument("output_prefix")
    args = parser.parse_args()
    render(args.run_dir, args.output_prefix)
