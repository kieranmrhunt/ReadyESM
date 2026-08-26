#!/usr/bin/env python3
"""Download one dense, instantaneous global ERA5 atmosphere initial state."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import cdsapi
from netCDF4 import Dataset, num2date


PRESSURE_LEVELS = [
    "1000", "975", "950", "925", "900", "875", "850", "825", "800",
    "775", "750", "700", "650", "600", "550", "500", "450", "400",
    "350", "300", "250", "225", "200", "175", "150", "125", "100",
]
STRATOSPHERE_PRESSURE_LEVELS = [
    *PRESSURE_LEVELS,
    "70", "50", "30", "20", "10", "7", "5", "3", "2", "1",
]
PRESSURE_VARIABLES = [
    "temperature",
    "u_component_of_wind",
    "v_component_of_wind",
    "specific_humidity",
]
SINGLE_VARIABLES = [
    "surface_pressure",
    "2m_temperature",
    "total_column_water_vapour",
]
LAND_VARIABLES = [
    "land_sea_mask",
    "soil_type",
    "skin_temperature",
    "volumetric_soil_water_layer_1",
    "volumetric_soil_water_layer_2",
    "volumetric_soil_water_layer_3",
    "volumetric_soil_water_layer_4",
    "soil_temperature_level_1",
    "soil_temperature_level_2",
    "soil_temperature_level_3",
    "soil_temperature_level_4",
    "high_vegetation_cover",
    "low_vegetation_cover",
    "leaf_area_index_high_vegetation",
    "leaf_area_index_low_vegetation",
]
LAND_SHORT_NAMES = [
    "lsm", "slt", "skt",
    "swvl1", "swvl2", "swvl3", "swvl4",
    "stl1", "stl2", "stl3", "stl4",
    "cvh", "cvl", "lai_hv", "lai_lv",
]
RADIATION_VARIABLES = [
    "top_net_solar_radiation",
    "top_net_thermal_radiation",
    "top_net_solar_radiation_clear_sky",
    "top_net_thermal_radiation_clear_sky",
    "toa_incident_solar_radiation",
]
RADIATION_SHORT_NAMES = ["tsr", "ttr", "tsrc", "ttrc", "tisr"]


def parse_timestamp(value: str) -> datetime:
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    parsed = parsed.astimezone(timezone.utc)
    if parsed.minute or parsed.second or parsed.microsecond:
        raise argparse.ArgumentTypeError("ERA5 timestamp must be on an exact UTC hour")
    return parsed


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def time_values(dataset: Dataset) -> list[datetime]:
    name = "valid_time" if "valid_time" in dataset.variables else "time"
    variable = dataset.variables[name]
    values = num2date(
        variable[:],
        variable.units,
        calendar=getattr(variable, "calendar", "standard"),
        only_use_cftime_datetimes=False,
        only_use_python_datetimes=True,
    )
    return [value.replace(tzinfo=timezone.utc) for value in values]


def validate_file(
    path: Path,
    expected_variables: list[str],
    timestamp: datetime,
    pressure_file: bool,
) -> None:
    if not path.is_file() or path.stat().st_size < 100_000:
        raise RuntimeError(f"ERA5 response is missing or unexpectedly small: {path}")
    with Dataset(path, "r") as dataset:
        missing = [name for name in expected_variables if name not in dataset.variables]
        if missing:
            raise RuntimeError(f"ERA5 response {path} lacks variables {missing}")
        times = time_values(dataset)
        if times != [timestamp]:
            raise RuntimeError(
                f"ERA5 response {path} has timestamps {times}, expected {[timestamp]}"
            )
        for name in expected_variables:
            variable = dataset.variables[name]
            step_type = getattr(variable, "GRIB_stepType", None)
            if step_type != "instant":
                raise RuntimeError(
                    f"ERA5 variable {name} is GRIB_stepType={step_type!r}, not 'instant'"
                )
        if pressure_file:
            coordinate = dataset.variables.get("pressure_level")
            if coordinate is None:
                raise RuntimeError("ERA5 pressure response lacks pressure_level")
            levels = [int(round(value)) for value in coordinate[:].tolist()]
            expected = [int(value) for value in PRESSURE_LEVELS]
            if levels != expected:
                raise RuntimeError(
                    f"ERA5 pressure levels are {levels}, expected {expected}"
                )


def validate_radiation_file(path: Path, timestamp: datetime) -> None:
    if not path.is_file() or path.stat().st_size < 100_000:
        raise RuntimeError(
            f"ERA5 radiation response is missing or unexpectedly small: {path}"
        )
    with Dataset(path, "r") as dataset:
        missing = [
            name for name in RADIATION_SHORT_NAMES if name not in dataset.variables
        ]
        if missing:
            raise RuntimeError(f"ERA5 radiation response {path} lacks {missing}")
        times = time_values(dataset)
        if times != [timestamp]:
            raise RuntimeError(
                f"ERA5 radiation response {path} has timestamps {times}, "
                f"expected {[timestamp]}"
            )
        for name in RADIATION_SHORT_NAMES:
            variable = dataset.variables[name]
            step_type = getattr(variable, "GRIB_stepType", None)
            if step_type not in {"accum", "avg"}:
                raise RuntimeError(
                    f"ERA5 radiation variable {name} is "
                    f"GRIB_stepType={step_type!r}, not accumulated/averaged"
                )
            units = str(getattr(variable, "units", "")).replace(" ", "").lower()
            if units not in {"jm**-2", "jm^-2", "wm**-2", "wm^-2"}:
                raise RuntimeError(
                    f"ERA5 radiation variable {name} has unexpected units "
                    f"{getattr(variable, 'units', None)!r}"
                )


def validate_stratosphere_file(path: Path, timestamp: datetime) -> None:
    validate_file(path, ["t", "u", "v", "q"], timestamp, False)
    with Dataset(path, "r") as dataset:
        coordinate = dataset.variables.get("pressure_level")
        if coordinate is None:
            raise RuntimeError("ERA5 stratosphere response lacks pressure_level")
        levels = [int(round(value)) for value in coordinate[:].tolist()]
        expected = [int(value) for value in STRATOSPHERE_PRESSURE_LEVELS]
        if levels != expected:
            raise RuntimeError(
                f"ERA5 stratosphere levels are {levels}, expected {expected}"
            )


def retrieve(
    client: cdsapi.Client,
    dataset: str,
    request: dict[str, object],
    target: Path,
    expected_variables: list[str],
    timestamp: datetime,
    pressure_file: bool,
    force: bool,
) -> str:
    if target.exists() and not force:
        validate_file(target, expected_variables, timestamp, pressure_file)
        return "cached"
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + ".part")
    if temporary.exists():
        temporary.unlink()
    client.retrieve(dataset, request, str(temporary))
    validate_file(temporary, expected_variables, timestamp, pressure_file)
    temporary.replace(target)
    return "downloaded"


def retrieve_radiation(
    client: cdsapi.Client,
    request: dict[str, object],
    target: Path,
    timestamp: datetime,
    force: bool,
) -> str:
    if target.exists() and not force:
        validate_radiation_file(target, timestamp)
        return "cached"
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + ".part")
    if temporary.exists():
        temporary.unlink()
    client.retrieve("reanalysis-era5-single-levels", request, str(temporary))
    validate_radiation_file(temporary, timestamp)
    temporary.replace(target)
    return "downloaded"


def retrieve_stratosphere(
    client: cdsapi.Client,
    request: dict[str, object],
    target: Path,
    timestamp: datetime,
    force: bool,
) -> str:
    if target.exists() and not force:
        validate_stratosphere_file(target, timestamp)
        return "cached"
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + ".part")
    if temporary.exists():
        temporary.unlink()
    client.retrieve("reanalysis-era5-pressure-levels", request, str(temporary))
    validate_stratosphere_file(temporary, timestamp)
    temporary.replace(target)
    return "downloaded"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--timestamp",
        type=parse_timestamp,
        default=parse_timestamp("1993-01-01T00:00:00Z"),
        help="Exact UTC reanalysis time (default: 1993-01-01T00:00:00Z).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/era5_initial_state"),
    )
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--force-land",
        action="store_true",
        help="Refresh only the land-state file while reusing valid atmosphere files.",
    )
    parser.add_argument(
        "--include-radiation",
        action="store_true",
        help=(
            "Also download the five accumulated ERA5 TOA radiation fields for "
            "an all-sky/clear-sky RRTMGP benchmark."
        ),
    )
    parser.add_argument(
        "--include-stratosphere",
        action="store_true",
        help=(
            "Also download a separate 1000--1 hPa atmosphere profile so the "
            "model top is not extrapolated above the standard 100 hPa file."
        ),
    )
    args = parser.parse_args()

    timestamp: datetime = args.timestamp
    stamp = timestamp.strftime("%Y%m%dT%H%MZ")
    pressure_path = args.output_dir / f"era5_pressure_levels_{stamp}.nc"
    single_path = args.output_dir / f"era5_single_levels_{stamp}.nc"
    land_path = args.output_dir / f"era5_land_state_{stamp}.nc"
    radiation_path = args.output_dir / f"era5_toa_radiation_{stamp}.nc"
    stratosphere_path = (
        args.output_dir / f"era5_pressure_levels_stratosphere_{stamp}.nc"
    )
    common = {
        "product_type": ["reanalysis"],
        "year": [timestamp.strftime("%Y")],
        "month": [timestamp.strftime("%m")],
        "day": [timestamp.strftime("%d")],
        "time": [timestamp.strftime("%H:00")],
        "grid": [1.0, 1.0],
        "data_format": "netcdf",
        "download_format": "unarchived",
    }
    pressure_request = {
        **common,
        "variable": PRESSURE_VARIABLES,
        "pressure_level": PRESSURE_LEVELS,
    }
    single_request = {**common, "variable": SINGLE_VARIABLES}
    land_request = {**common, "variable": LAND_VARIABLES}
    radiation_request = {**common, "variable": RADIATION_VARIABLES}
    stratosphere_request = {
        **common,
        "variable": PRESSURE_VARIABLES,
        "pressure_level": STRATOSPHERE_PRESSURE_LEVELS,
    }

    client = cdsapi.Client(quiet=False)
    pressure_status = retrieve(
        client,
        "reanalysis-era5-pressure-levels",
        pressure_request,
        pressure_path,
        ["t", "u", "v", "q"],
        timestamp,
        True,
        args.force,
    )
    single_status = retrieve(
        client,
        "reanalysis-era5-single-levels",
        single_request,
        single_path,
        ["sp", "t2m", "tcwv"],
        timestamp,
        False,
        args.force,
    )
    land_status = retrieve(
        client,
        "reanalysis-era5-single-levels",
        land_request,
        land_path,
        LAND_SHORT_NAMES,
        timestamp,
        False,
        args.force or args.force_land,
    )
    radiation_status = None
    if args.include_radiation:
        radiation_status = retrieve_radiation(
            client,
            radiation_request,
            radiation_path,
            timestamp,
            args.force,
        )
    stratosphere_status = None
    if args.include_stratosphere:
        stratosphere_status = retrieve_stratosphere(
            client,
            stratosphere_request,
            stratosphere_path,
            timestamp,
            args.force,
        )
    receipt = {
        "timestamp_utc": timestamp.isoformat(),
        "pressure_levels_hpa": [int(value) for value in PRESSURE_LEVELS],
        "pressure_levels": {
            "status": pressure_status,
            "dataset": "reanalysis-era5-pressure-levels",
            "request": pressure_request,
            "path": str(pressure_path.resolve()),
            "bytes": pressure_path.stat().st_size,
            "sha256": sha256(pressure_path),
        },
        "single_levels": {
            "status": single_status,
            "dataset": "reanalysis-era5-single-levels",
            "request": single_request,
            "path": str(single_path.resolve()),
            "bytes": single_path.stat().st_size,
            "sha256": sha256(single_path),
        },
        "land_state": {
            "status": land_status,
            "dataset": "reanalysis-era5-single-levels",
            "request": land_request,
            "path": str(land_path.resolve()),
            "bytes": land_path.stat().st_size,
            "sha256": sha256(land_path),
        },
    }
    if args.include_radiation:
        receipt["toa_radiation"] = {
            "status": radiation_status,
            "dataset": "reanalysis-era5-single-levels",
            "request": radiation_request,
            "path": str(radiation_path.resolve()),
            "bytes": radiation_path.stat().st_size,
            "sha256": sha256(radiation_path),
            "processing_period": "one hour ending at valid time (ERA5 reanalysis)",
        }
    if args.include_stratosphere:
        receipt["pressure_levels_stratosphere"] = {
            "status": stratosphere_status,
            "dataset": "reanalysis-era5-pressure-levels",
            "request": stratosphere_request,
            "path": str(stratosphere_path.resolve()),
            "bytes": stratosphere_path.stat().st_size,
            "sha256": sha256(stratosphere_path),
        }
    receipt_path = args.output_dir / f"era5_initial_state_{stamp}.json"
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    main()
