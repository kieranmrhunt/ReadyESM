#!/usr/bin/env python3
"""Make sampled 1993 ERA5 zonal upper-cap temperatures from L137 model levels."""

import argparse
import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path

import netCDF4
import numpy as np


HERE = Path(__file__).resolve().parent
ARCHIVE = Path("/badc/ecmwf-era5/data/oper/an_ml")
COEFFICIENTS = HERE / "ecmwf_l137_half_levels_0_24.csv"
COEFFICIENTS_SHA = "44d92619fd206f9258e1a2ee131285ac8542031ff2170964da656976cb573728"
ECMWF_TABLE_URL = "https://confluence.ecmwf.int/spaces/UDOC/pages/108117123/L137%2Bmodel%2Blevel%2Bdefinitions"
ECMWF_PAGE_SHA = "db6787b4339832df937f77b518eb1a2b73c2013b08b9548fce1e09d61f3ea34e"
TOP_PA = 1.005183574463
SPLIT_PA = 100.0
BOTTOM_PA = 500.0
DAYS = (1, 8, 15, 22)
HOURS = (0, 12)
FIELD_NAMES = (
    "temperature_5_1hpa_k", "temperature_1hpa_1pa_k",
    "temperature_5hpa_boundary_k", "temperature_1hpa_boundary_k",
    "temperature_top_1pa_k",
)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def half_levels():
    assert sha256(COEFFICIENTS) == COEFFICIENTS_SHA
    with COEFFICIENTS.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 25
    assert [int(row["half_level"]) for row in rows] == list(range(25))
    assert all(float(row["b"]) == 0.0 for row in rows)
    a = np.array([float(row["a_pa"]) for row in rows])
    assert a[0] == 0.0 and np.all(np.diff(a) > 0)
    assert a[1] > TOP_PA and a[-1] > BOTTOM_PA
    return a


def weights(a, low, high):
    overlap = np.maximum(0.0, np.minimum(a[1:], high) - np.maximum(a[:-1], low))
    assert np.isclose(overlap.sum(), high - low, atol=1e-9, rtol=0)
    return overlap / overlap.sum()


def boundary(zonal, full_pressure, target):
    upper = int(np.searchsorted(full_pressure, target))
    assert 0 < upper < len(full_pressure)
    lower = upper - 1
    fraction = (np.log(target) - np.log(full_pressure[lower])) / (
        np.log(full_pressure[upper]) - np.log(full_pressure[lower])
    )
    return (1 - fraction) * zonal[lower] + fraction * zonal[upper]


def paths():
    for month in range(1, 13):
        for day in DAYS:
            for hour in HOURS:
                timestamp = datetime(1993, month, day, hour)
                stamp = timestamp.strftime("%Y%m%d%H%M")
                path = ARCHIVE / "1993" / f"{month:02d}" / f"{day:02d}" / (
                    f"ecmwf-era5_oper_an_ml_{stamp}.t.nc"
                )
                yield timestamp, path


def read_zonal(path, timestamp):
    assert path.is_file() and path.stat().st_size > 1_000_000
    with netCDF4.Dataset(path) as ds:
        assert set(ds.dimensions) == {"time", "level", "latitude", "longitude"}
        assert tuple(len(ds.dimensions[name]) for name in
                     ("time", "level", "latitude", "longitude")) == (1, 137, 721, 1440)
        assert np.array_equal(ds.variables["level"][:24], np.arange(1, 25))
        lat = np.asarray(ds.variables["latitude"][::4], dtype=float)
        lon = np.asarray(ds.variables["longitude"][:], dtype=float)
        assert lat.shape == (181,) and np.allclose(lat, np.arange(90, -91, -1), atol=1e-6)
        assert lon.shape == (1440,) and np.allclose(lon, np.arange(1440) / 4, atol=1e-6)
        stamp = netCDF4.num2date(ds.variables["time"][0], ds.variables["time"].units,
                                 calendar=getattr(ds.variables["time"], "calendar", "standard"))
        assert stamp.strftime("%Y%m%d%H%M") == timestamp.strftime("%Y%m%d%H%M")
        field = ds.variables["t"]
        assert field.dimensions == ("time", "level", "latitude", "longitude")
        temperature = np.asarray(np.ma.filled(field[0, :24, ::4, :], np.nan), dtype=float)
        assert temperature.shape == (24, 181, 1440)
        assert np.all(np.isfinite(temperature)) and np.all((temperature > 150) & (temperature < 350))
    return temperature.mean(axis=2)[:, ::-1]


def build(output):
    assert not output.exists()
    output.mkdir(parents=True)
    a = half_levels()
    full_pressure = (a[:-1] + a[1:]) / 2
    lower_weights = weights(a, SPLIT_PA, BOTTOM_PA)
    upper_weights = weights(a, TOP_PA, SPLIT_PA)
    assert np.isclose(lower_weights.sum(), 1) and np.isclose(upper_weights.sum(), 1)
    months = {month: [] for month in range(1, 13)}
    inputs = []
    for timestamp, path in paths():
        zonal = read_zonal(path, timestamp)
        values = np.stack((
            lower_weights @ zonal,
            upper_weights @ zonal,
            boundary(zonal, full_pressure, BOTTOM_PA),
            boundary(zonal, full_pressure, SPLIT_PA),
            zonal[0],
        ), axis=1)
        assert values.shape == (181, 5) and np.all(np.isfinite(values))
        months[timestamp.month].append(values)
        inputs.append({"timestamp_utc": timestamp.strftime("%Y-%m-%dT%H:%M:%S"),
                       "path": str(path), "bytes": path.stat().st_size,
                       "sha256": sha256(path)})
        print("ERA5_CAP_SAMPLE_PASS", timestamp.strftime("%Y-%m-%dT%H:%M"), flush=True)
    assert len(inputs) == 96 and all(len(months[month]) == 8 for month in months)
    csv_path = output / "era5_1993_sampled_monthly_upper_cap.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("month", "latitude_degrees") + FIELD_NAMES)
        for month in range(1, 13):
            mean = np.mean(months[month], axis=0)
            assert mean.shape == (181, 5) and np.all((mean > 150) & (mean < 350))
            for lat, row in zip(range(-90, 91), mean):
                writer.writerow((month, lat) + tuple(f"{value:.9f}" for value in row))
    provenance = {
        "schema": "readyesm_era5_1993_sampled_monthly_upper_cap_v1",
        "ecmwf_l137_table_url": ECMWF_TABLE_URL,
        "ecmwf_l137_page_sha256_at_preparation": ECMWF_PAGE_SHA,
        "coefficient_csv_sha256": COEFFICIENTS_SHA,
        "archive": str(ARCHIVE),
        "temperature_field": "ERA5 an_ml t, model levels 1–24, 0.25-degree grid",
        "pressure_method": "ECMWF L137 half-level a+b*surface_pressure, with b=0 for all used half-levels",
        "layer_mean_method": "Pressure-overlap mass weighted mean of model-level temperatures",
        "boundary_method": "Log-pressure interpolation between neighboring model-level midpoint pressures",
        "longitude_method": "Unweighted mean of all 1440 longitudes at each integer latitude",
        "temporal_method": "Eight snapshots per month: days 1,8,15,22 at 00 and 12 UTC; sampled monthly estimate, not full monthly mean",
        "cap_top_pressure_pa": TOP_PA,
        "cap_split_pressure_pa": SPLIT_PA,
        "cap_bottom_pressure_pa": BOTTOM_PA,
        "latitude_degrees": list(range(-90, 91)),
        "sample_count": 96,
        "source_files": inputs,
        "output_csv_sha256": sha256(csv_path),
    }
    provenance_path = output / "provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    manifest = output / "result.sha256"
    manifest.write_text("".join(f"{sha256(output / name)}  {name}\n" for name in
                                ("era5_1993_sampled_monthly_upper_cap.csv", "provenance.json")))
    print("ERA5_1993_SAMPLED_MONTHLY_CAP_PASS", sha256(manifest), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    build(parser.parse_args().output)
