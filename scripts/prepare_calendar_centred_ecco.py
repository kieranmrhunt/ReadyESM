#!/usr/bin/env python3
"""Build the qualified 1993-01-01 ReadyESM ocean and sea-ice launch bundle.

Usage:
    python -I scripts/prepare_calendar_centred_ecco.py \
        DECEMBER_1992_DIRECTORY JANUARY_1993_DIRECTORY OUTPUT_DIRECTORY
    python -I scripts/prepare_calendar_centred_ecco.py --self-test

The five ocean fields are evaluated halfway between the December 1992 and
January 1993 ECCO monthly means. The two sea-ice fields are copied from
January. Source receipts and every output are pinned to the qualified BS04
experiment.
"""

import datetime as dt
import hashlib
import json
import os
import shutil
import sys
import tempfile

import netCDF4 as nc
import numpy as np


OCEAN_VARIABLES = ("THETA", "SALT", "EVEL", "NVEL", "SSH")
ICE_VARIABLES = ("SIarea", "SIheff")
ALL_VARIABLES = OCEAN_VARIABLES + ICE_VARIABLES
LAUNCH = dt.datetime(1993, 1, 1, 0, 0)
FILL = np.float32(9.96921e36)
DECEMBER_DIMS = {
    4: ("time", "Z", "latitude", "longitude"),
    3: ("time", "latitude", "longitude"),
}
JANUARY_DIMS = {
    4: ("time", "k", "j", "i"),
    3: ("time", "j", "i"),
}

EXPECTED_RECEIPTS = {
    "december": "a9d8c07d16cc0c1ee4354547ff0d89cf7946f448ffc18d964f3a5e87e60f9010",
    "january": "7b4966c32c1bc649be79a3640b67c2d81a6c326c0b6876abce282059f5fd02c1",
}
EXPECTED_DECEMBER = {
    "THETA": "9630835fe35a2f032abf401ceb3d2fc6e31ad87dc9b2ce7b75192b599e7a1732",
    "SALT": "9630835fe35a2f032abf401ceb3d2fc6e31ad87dc9b2ce7b75192b599e7a1732",
    "EVEL": "502649bae0079a8418e6e2112c96ecd5695c9682e0b06546dfa31db2ba5acd14",
    "NVEL": "502649bae0079a8418e6e2112c96ecd5695c9682e0b06546dfa31db2ba5acd14",
    "SSH": "de11060277c46a4e919eb68c1f39b9b7ed1651febfd7b22a590cf2a2925a5cfb",
}
EXPECTED_JANUARY = {
    "THETA": "0b5e88bc1e865a8bdc566d858cda24333a967705eb731af2171e8a941fd318db",
    "SALT": "8e8bdcfcefed4ef0e81ab4be3d6b1e2e9c6b99b3036dfcd0e4922ac78f596f99",
    "EVEL": "d3d2927d3af6d5860ab821447bf2aef44e303e3d7e5a1988d098c0d6b5b7f68d",
    "NVEL": "2ea1731bd8f5c653e0c89f9936d79865d47360df7a2ff2d044b077464301badf",
    "SSH": "95b48394993f4dd8e87e383b35c70253a58c2ae1e81e3f1c0196c464b2c85488",
    "SIarea": "04518f541120851fc16484b065f5be520612ee17f6ec620cc3645bbd93308717",
    "SIheff": "36699cf11a99402782473c5fc051aa3e23dabadc46381cc7c03837c205353440",
}
EXPECTED_OUTPUT = {
    "THETA": "305e40f0a318c003e172ef30524e42cfc26ff27185d96d9e134334307dadbf1f",
    "SALT": "e7bc983f2be6cc917e47736947b19a027b300fd79b4b248798746f5078936aec",
    "EVEL": "8f9f925a60c0daa10ceb4d4c3f2aae913b3f8921add5d3ac2f5f9df7688a40d0",
    "NVEL": "f19e9d239135b7d9e14c768b82596d8dab139cdf62581c31583cbc715e8a489a",
    "SSH": "71a4e330563957d2ef74c262bb4afa1c60c2ee92b0d09857ecbf47026a7be9e8",
    "SIarea": "04518f541120851fc16484b065f5be520612ee17f6ec620cc3645bbd93308717",
    "SIheff": "36699cf11a99402782473c5fc051aa3e23dabadc46381cc7c03837c205353440",
}


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_sha256(path, expected, label):
    if not os.path.isfile(path):
        raise RuntimeError(f"missing {label}: {path}")
    actual = sha256(path)
    if actual != expected:
        raise RuntimeError(
            f"{label} SHA256 mismatch: expected {expected}, found {actual}"
        )
    return actual


def raw(dataset, variable_name):
    variable = dataset.variables[variable_name]
    variable.set_auto_maskandscale(False)
    return np.array(variable[:])


def centre(dataset):
    time = dataset.variables["time"]
    return nc.num2date(
        time[:],
        time.units,
        only_use_cftime_datetimes=False,
        only_use_python_datetimes=True,
    )[0]


def check_physical(variable_name, values, location):
    if not np.all(np.isfinite(values)):
        raise RuntimeError(f"{location} {variable_name} contains a non-finite value")
    if variable_name == "SALT" and not np.all(values > 0):
        raise RuntimeError(f"{location} SALT contains ocean salinity at or below zero")
    if variable_name == "THETA" and not np.all((values >= -3) & (values <= 40)):
        raise RuntimeError(f"{location} THETA lies outside [-3, 40] degC")
    if variable_name == "SIarea" and not np.all((values >= 0) & (values <= 1)):
        raise RuntimeError(f"{location} SIarea lies outside [0, 1]")
    if variable_name == "SIheff" and not np.all(values >= 0):
        raise RuntimeError(f"{location} SIheff is negative")


def source_paths(december_directory, january_directory, variable_name):
    december = os.path.join(december_directory, f"{variable_name}_1992_12.nc")
    january = os.path.join(january_directory, f"{variable_name}_1993_01.nc")
    return december, january


def verify_pinned_sources(december_directory, january_directory):
    receipt_paths = {
        "december": os.path.join(
            december_directory, "ECCO4Monthly_1992_12_RECEIPT.txt"
        ),
        "january": os.path.join(
            january_directory, "ECCO4Monthly_1993_01_RECEIPT.txt"
        ),
    }
    receipt_hashes = {
        name: require_sha256(path, EXPECTED_RECEIPTS[name], f"{name} receipt")
        for name, path in receipt_paths.items()
    }
    source_hashes = {"december": {}, "january": {}}
    for variable_name in OCEAN_VARIABLES:
        december, january = source_paths(
            december_directory, january_directory, variable_name
        )
        source_hashes["december"][variable_name] = require_sha256(
            december,
            EXPECTED_DECEMBER[variable_name],
            f"December {variable_name}",
        )
        source_hashes["january"][variable_name] = require_sha256(
            january,
            EXPECTED_JANUARY[variable_name],
            f"January {variable_name}",
        )
    for variable_name in ICE_VARIABLES:
        january = os.path.join(
            january_directory, f"{variable_name}_1993_01.nc"
        )
        source_hashes["january"][variable_name] = require_sha256(
            january,
            EXPECTED_JANUARY[variable_name],
            f"January {variable_name}",
        )
    return receipt_hashes, source_hashes


def prepare(december_directory, january_directory, output_directory, enforce_pins=True):
    december_directory = os.path.abspath(december_directory)
    january_directory = os.path.abspath(january_directory)
    output_directory = os.path.abspath(output_directory)
    working_directory = output_directory + ".tmp"
    if os.path.exists(output_directory):
        raise RuntimeError(f"refusing to replace existing output: {output_directory}")
    if os.path.exists(working_directory):
        raise RuntimeError(f"temporary output already exists: {working_directory}")
    if not os.path.isdir(december_directory):
        raise RuntimeError(f"December directory does not exist: {december_directory}")
    if not os.path.isdir(january_directory):
        raise RuntimeError(f"January directory does not exist: {january_directory}")

    if enforce_pins:
        receipt_hashes, source_hashes = verify_pinned_sources(
            december_directory, january_directory
        )
    else:
        receipt_hashes = {}
        source_hashes = {"december": {}, "january": {}}

    january_salt = os.path.join(january_directory, "SALT_1993_01.nc")
    january_ssh = os.path.join(january_directory, "SSH_1993_01.nc")
    with nc.Dataset(january_salt) as salt_dataset, nc.Dataset(january_ssh) as ssh_dataset:
        land_3d = raw(salt_dataset, "SALT") == 0
        land_2d = raw(ssh_dataset, "SSH") == 0
        if np.isnan(raw(salt_dataset, "SALT")).any():
            raise RuntimeError("January SALT contains NaN")
        if np.isnan(raw(ssh_dataset, "SSH")).any():
            raise RuntimeError("January SSH contains NaN")
    if not np.array_equal(land_3d[:, 0], land_2d):
        raise RuntimeError("January SALT and SSH surface land masks differ")

    os.makedirs(working_directory)
    report = {
        "schema": "readiesm_calendar_centred_ecco_receipt_v1",
        "launch": LAUNCH.isoformat(),
        "ice_rule": "January 1993 monthly state copied without interpolation",
        "land_cells_3d": int(land_3d.sum()),
        "land_cells_2d": int(land_2d.sum()),
        "source_receipt_sha256": receipt_hashes,
        "source_sha256": source_hashes,
        "variables": {},
    }

    for variable_name in OCEAN_VARIABLES:
        december_path, january_path = source_paths(
            december_directory, january_directory, variable_name
        )
        output_path = os.path.join(
            working_directory, f"{variable_name}_1993_01.nc"
        )
        with nc.Dataset(december_path) as december_dataset, nc.Dataset(
            january_path
        ) as january_dataset:
            december_centre = centre(december_dataset)
            january_centre = centre(january_dataset)
            if not december_centre < LAUNCH < january_centre:
                raise RuntimeError(
                    f"{variable_name} source times do not bracket the launch"
                )
            january_weight = (LAUNCH - december_centre) / (
                january_centre - december_centre
            )
            december_weight = 1 - january_weight
            for coordinate in ("latitude", "longitude", "Z"):
                if coordinate in january_dataset.variables:
                    december_coordinate = np.asarray(
                        december_dataset.variables[coordinate][:], np.float64
                    )
                    january_coordinate = np.asarray(
                        january_dataset.variables[coordinate][:], np.float64
                    )
                    if not np.array_equal(december_coordinate, january_coordinate):
                        raise RuntimeError(
                            f"{variable_name} {coordinate} coordinates differ"
                        )

            fill_value = december_dataset.variables[variable_name].getncattr(
                "_FillValue"
            )
            if np.float32(fill_value) != FILL:
                raise RuntimeError(
                    f"{variable_name} has unexpected December fill value"
                )
            december_values = raw(december_dataset, variable_name)
            january_values = raw(january_dataset, variable_name)
            if december_values.shape != january_values.shape:
                raise RuntimeError(f"{variable_name} source shapes differ")
            if december_dataset.variables[variable_name].dimensions != DECEMBER_DIMS[
                december_values.ndim
            ]:
                raise RuntimeError(f"{variable_name} has unexpected December dimensions")
            if january_dataset.variables[variable_name].dimensions != JANUARY_DIMS[
                january_values.ndim
            ]:
                raise RuntimeError(f"{variable_name} has unexpected January dimensions")
            if december_values.dtype != np.float32 or january_values.dtype != np.float64:
                raise RuntimeError(f"{variable_name} has unexpected source precision")

            land = land_3d if december_values.ndim == 4 else land_2d
            if land.shape != december_values.shape:
                raise RuntimeError(f"{variable_name} land mask shape differs")
            if not np.array_equal(december_values == FILL, land):
                raise RuntimeError(
                    f"{variable_name} December fill mask differs from January land"
                )
            if not np.all(january_values[land] == 0):
                raise RuntimeError(f"{variable_name} January land is not zero")
            check_physical(
                variable_name,
                december_values[~land].astype(np.float64),
                "December",
            )
            check_physical(variable_name, january_values[~land], "January")
            output_values = january_values.copy()
            output_values[~land] = (
                december_weight * december_values[~land].astype(np.float64)
                + january_weight * january_values[~land]
            )

        check_physical(variable_name, output_values[~land], "output")
        lower = np.minimum(
            december_values[~land].astype(np.float64), january_values[~land]
        )
        upper = np.maximum(
            december_values[~land].astype(np.float64), january_values[~land]
        )
        if not np.all((output_values[~land] >= lower) & (output_values[~land] <= upper)):
            raise RuntimeError(f"{variable_name} output lies outside its sources")

        shutil.copyfile(january_path, output_path)
        with nc.Dataset(output_path, "r+") as output_dataset:
            output_variable = output_dataset.variables[variable_name]
            output_variable.set_auto_maskandscale(False)
            if output_variable.dimensions != JANUARY_DIMS[output_values.ndim]:
                raise RuntimeError(f"{variable_name} output dimensions changed")
            output_variable[:] = output_values
            output_dataset.setncattr(
                "readiesm_bs04_l1",
                f"{december_weight:.6f} * {os.path.basename(december_path)} "
                f"(PO.DAAC, fill -> land) + {january_weight:.6f} * "
                f"{os.path.basename(january_path)}, valid {LAUNCH.isoformat()}",
            )

        output_hash = sha256(output_path)
        if enforce_pins and output_hash != EXPECTED_OUTPUT[variable_name]:
            raise RuntimeError(
                f"{variable_name} output SHA256 mismatch: expected "
                f"{EXPECTED_OUTPUT[variable_name]}, found {output_hash}"
            )
        report["variables"][variable_name] = {
            "mode": "linear_ocean_interpolation",
            "december_centre": december_centre.isoformat(),
            "january_centre": january_centre.isoformat(),
            "december_weight": december_weight,
            "january_weight": january_weight,
            "output_sha256": output_hash,
            "ocean_cells": int((~land).sum()),
            "mean_abs_change_from_january": float(
                np.mean(np.abs(output_values[~land] - january_values[~land]))
            ),
            "max_abs_change_from_january": float(
                np.max(np.abs(output_values[~land] - january_values[~land]))
            ),
        }

    for variable_name in ICE_VARIABLES:
        january_path = os.path.join(
            january_directory, f"{variable_name}_1993_01.nc"
        )
        output_path = os.path.join(
            working_directory, f"{variable_name}_1993_01.nc"
        )
        with nc.Dataset(january_path) as january_dataset:
            january_values = raw(january_dataset, variable_name)
        if january_values.shape != land_2d.shape:
            raise RuntimeError(f"{variable_name} January mask shape differs")
        if not np.all(january_values[land_2d] == 0):
            raise RuntimeError(f"{variable_name} January land is not zero")
        check_physical(variable_name, january_values[~land_2d], "January")
        shutil.copyfile(january_path, output_path)
        output_hash = sha256(output_path)
        if enforce_pins and output_hash != EXPECTED_OUTPUT[variable_name]:
            raise RuntimeError(
                f"{variable_name} output SHA256 mismatch: expected "
                f"{EXPECTED_OUTPUT[variable_name]}, found {output_hash}"
            )
        report["variables"][variable_name] = {
            "mode": "january_ice_copy",
            "output_sha256": output_hash,
            "ocean_cells": int((~land_2d).sum()),
        }

    receipt_path = os.path.join(
        working_directory, "CALENDAR_CENTRED_ECCO_1993_01_RECEIPT.json"
    )
    with open(receipt_path, "x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    os.rename(working_directory, output_directory)
    print(
        "CALENDAR_CENTRED_ECCO_PASS "
        f"ocean_fields={len(OCEAN_VARIABLES)} ice_fields={len(ICE_VARIABLES)} "
        f"output={output_directory}"
    )
    return report


def write_synthetic(path, variable_name, values, source_format):
    nz, ny, nx = values.shape[-3:]
    with nc.Dataset(path, "w") as dataset:
        three_dimensional = values.ndim == 4
        dataset.createDimension("time", 1)
        if source_format == "december":
            dimensions = (
                ("time", "Z", "latitude", "longitude")
                if three_dimensional
                else ("time", "latitude", "longitude")
            )
            dataset.createDimension("latitude", ny)
            dataset.createDimension("longitude", nx)
            if three_dimensional:
                dataset.createDimension("Z", nz)
            variable = dataset.createVariable(
                variable_name, "f4", dimensions, fill_value=FILL
            )
            time = dataset.createVariable("time", "f8", ("time",))
            time.units = "hours since 1992-01-01T12:00:00"
            time[:] = [8400]
            z_name, y_name, x_name = "Z", "latitude", "longitude"
        else:
            dimensions = (
                ("time", "k", "j", "i")
                if three_dimensional
                else ("time", "j", "i")
            )
            dataset.createDimension("j", ny)
            dataset.createDimension("i", nx)
            if three_dimensional:
                dataset.createDimension("k", nz)
            variable = dataset.createVariable(variable_name, "f8", dimensions)
            time = dataset.createVariable("time", "f8", ("time",))
            time.units = "days since 1993-01-16 12:00:00"
            time[:] = [0]
            z_name, y_name, x_name = "k", "j", "i"
        dataset.createVariable("latitude", "f8", (y_name,))[:] = np.linspace(
            -60, 60, ny
        )
        dataset.createVariable("longitude", "f8", (x_name,))[:] = np.linspace(
            -150, 150, nx
        )
        if three_dimensional:
            dataset.createVariable("Z", "f8", (z_name,))[:] = [-5.0, -15.0, -25.0]
        variable.set_auto_maskandscale(False)
        variable[:] = values


def self_test():
    rng = np.random.default_rng(1)
    with tempfile.TemporaryDirectory() as root:
        december = os.path.join(root, "december")
        january = os.path.join(root, "january")
        os.makedirs(december)
        os.makedirs(january)
        nz, ny, nx = 3, 4, 6
        land_3d = np.zeros((1, nz, ny, nx), dtype=bool)
        land_3d[0, :, 0, :2] = True
        land_3d[0, 2, :, 5] = True
        land_2d = land_3d[:, 0]
        fields = {}
        for variable_name in ALL_VARIABLES:
            three_dimensional = variable_name in ("THETA", "SALT", "EVEL", "NVEL")
            land = land_3d if three_dimensional else land_2d
            january_values = rng.normal(10, 3, land.shape)
            december_values = rng.normal(10, 3, land.shape).astype(np.float32)
            if variable_name == "SIarea":
                january_values = rng.uniform(0, 1, land.shape)
                december_values = rng.uniform(0, 1, land.shape).astype(np.float32)
            if variable_name == "SIheff":
                january_values = np.abs(january_values)
                december_values = np.abs(december_values)
            if variable_name in ("EVEL", "SIarea"):
                january_values.flat[7] = 0
                december_values.flat[9] = 0
            january_values[land] = 0
            december_values[land] = FILL
            if variable_name == "SALT":
                january_values[~land] = np.abs(january_values[~land]) + 1
                december_values[~land] = np.abs(december_values[~land]) + 1
            write_synthetic(
                os.path.join(december, f"{variable_name}_1992_12.nc"),
                variable_name,
                december_values,
                "december",
            )
            write_synthetic(
                os.path.join(january, f"{variable_name}_1993_01.nc"),
                variable_name,
                january_values,
                "january",
            )
            fields[variable_name] = (december_values, january_values, land)

        output = os.path.join(root, "output")
        report = prepare(december, january, output, enforce_pins=False)
        for variable_name, (december_values, january_values, land) in fields.items():
            with nc.Dataset(
                os.path.join(output, f"{variable_name}_1993_01.nc")
            ) as dataset:
                output_values = raw(dataset, variable_name)
            if variable_name in OCEAN_VARIABLES:
                expected = (
                    0.5 * december_values[~land].astype(np.float64)
                    + 0.5 * january_values[~land]
                )
                if not np.array_equal(output_values[~land], expected):
                    raise AssertionError(f"{variable_name} interpolation differs")
                if report["variables"][variable_name]["january_weight"] != 0.5:
                    raise AssertionError(f"{variable_name} weight differs")
            elif not np.array_equal(output_values, january_values):
                raise AssertionError(f"{variable_name} is not an exact January copy")
            if not np.all(output_values[land] == 0):
                raise AssertionError(f"{variable_name} land changed")

        bad_mask = os.path.join(root, "bad_mask")
        shutil.copytree(december, bad_mask)
        december_values, _, land = fields["THETA"]
        changed = december_values.copy()
        changed[tuple(np.argwhere(~land)[0])] = FILL
        write_synthetic(
            os.path.join(bad_mask, "THETA_1992_12.nc"),
            "THETA",
            changed,
            "december",
        )
        try:
            prepare(bad_mask, january, os.path.join(root, "bad_mask_output"), False)
            raise AssertionError("mask mismatch was accepted")
        except RuntimeError as error:
            if "mask" not in str(error):
                raise

        bad_salinity = os.path.join(root, "bad_salinity")
        shutil.copytree(december, bad_salinity)
        december_values, _, land = fields["SALT"]
        changed = december_values.copy()
        changed[tuple(np.argwhere(~land)[0])] = -1
        write_synthetic(
            os.path.join(bad_salinity, "SALT_1992_12.nc"),
            "SALT",
            changed,
            "december",
        )
        try:
            prepare(
                bad_salinity,
                january,
                os.path.join(root, "bad_salinity_output"),
                False,
            )
            raise AssertionError("negative salinity was accepted")
        except RuntimeError as error:
            if "salinity" not in str(error):
                raise

        try:
            prepare(december, january, output, False)
            raise AssertionError("existing output was accepted")
        except RuntimeError as error:
            if "existing output" not in str(error):
                raise

        source = os.path.join(root, "hash_source")
        with open(source, "wb") as stream:
            stream.write(b"qualified")
        try:
            require_sha256(source, "0" * 64, "synthetic source")
            raise AssertionError("bad source hash was accepted")
        except RuntimeError as error:
            if "SHA256 mismatch" not in str(error):
                raise

    print("CALENDAR_CENTRED_ECCO_SELF_TEST_PASS")


if __name__ == "__main__":
    if sys.argv[1:] == ["--self-test"]:
        self_test()
    elif len(sys.argv) == 4:
        prepare(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        raise SystemExit(
            "usage: prepare_calendar_centred_ecco.py "
            "DECEMBER_1992_DIRECTORY JANUARY_1993_DIRECTORY OUTPUT_DIRECTORY"
        )
