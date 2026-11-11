#!/usr/bin/env python3
"""Make a compact historical ozone input with adjacent months for interpolation."""

import argparse
import hashlib
import json
from pathlib import Path

import netCDF4
import numpy as np

SOURCE_SHA256 = "873d52845442e8e5b833d9fbd29281f3f81337b54e8adc811c8c21fabcd66945"


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--first-year", type=int, default=1992)
    parser.add_argument("--last-year", type=int, default=1994)
    args = parser.parse_args()
    if args.last_year < args.first_year + 1:
        raise SystemExit("subset must span at least two years")
    years = tuple(range(args.first_year, args.last_year + 1))
    months_expected = 12 * len(years)
    if args.output.exists() or args.receipt.exists():
        raise SystemExit("refusing existing output or receipt")
    if sha256(args.source) != SOURCE_SHA256:
        raise SystemExit("historical ozone source SHA-256 differs")

    with netCDF4.Dataset(args.source) as original:
        variable_names = sorted(original.variables)
        year = np.asarray(original.variables["year"][:], dtype=np.int32)
        month = np.asarray(original.variables["month"][:], dtype=np.int32)
        selected = np.flatnonzero((year >= years[0]) & (year <= years[-1]))
        if len(selected) != months_expected or not np.array_equal(
            selected, np.arange(selected[0], selected[0] + months_expected)
        ):
            raise SystemExit("requested ozone interval is missing or discontinuous")
        if not np.array_equal(year[selected], np.repeat(years, 12)) or not np.array_equal(
            month[selected], np.tile(np.arange(1, 13), len(years))
        ):
            raise SystemExit("selected month labels differ")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with netCDF4.Dataset(args.output, "w", format="NETCDF4") as target:
            for key in original.ncattrs():
                target.setncattr(key, original.getncattr(key))
            target.setncattr("subset_source_sha256", SOURCE_SHA256)
            target.setncattr("subset_selection",
                f"{years[0]}-01 through {years[-1]}-12 inclusive; adjacent-month interpolation support")
            for name, dimension in original.dimensions.items():
                target.createDimension(name, months_expected if name == "time" else len(dimension))
            for name, variable in original.variables.items():
                fill = variable.getncattr("_FillValue") if "_FillValue" in variable.ncattrs() else None
                dest = target.createVariable(name, variable.datatype, variable.dimensions,
                    fill_value=fill, zlib=name == "vmro3", complevel=4)
                dest.setncatts({key: variable.getncattr(key) for key in variable.ncattrs()
                    if key != "_FillValue"})
                index = tuple(selected if dim == "time" else slice(None) for dim in variable.dimensions)
                dest[:] = variable[index]
        with netCDF4.Dataset(args.output) as target:
            for name, variable in original.variables.items():
                index = tuple(selected if dim == "time" else slice(None) for dim in variable.dimensions)
                if not np.array_equal(np.asarray(variable[index]), np.asarray(target.variables[name][:])):
                    raise SystemExit(f"subset variable {name} differs from source")

    receipt = {
        "schema": "readyesm_1993_ozone_subset_v1",
        "source_path": str(args.source.resolve()),
        "source_sha256": SOURCE_SHA256,
        "subset_path": str(args.output.resolve()),
        "subset_sha256": sha256(args.output),
        "first_month": f"{years[0]}-01",
        "last_month": f"{years[-1]}-12",
        "months": months_expected,
        "verified_variables": variable_names,
        "purpose": "historical runs with unchanged adjacent-month ozone interpolation",
    }
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"V021_OZONE_SUBSET_VALUES_EXACT years={years[0]}-{years[-1]}")
    print(f"subset_sha256={receipt['subset_sha256']}")


if __name__ == "__main__":
    main()
