#!/usr/bin/env python3
"""Check a compact ozone asset gives the same modeled-year interpolation brackets."""

import argparse
import hashlib
import json
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4
import numpy as np


def decimal_year(date):
    start = datetime(date.year, 1, 1)
    stop = datetime(date.year + 1, 1, 1)
    return date.year + (date - start) / (stop - start)


def bracket(times, value):
    if value <= times[0]:
        return 0, 0, 0.0
    if value >= times[-1]:
        return len(times) - 1, len(times) - 1, 0.0
    first = int(np.searchsorted(times, value, side="right") - 1)
    second = first + 1
    return first, second, float((value - times[first]) / (times[second] - times[first]))


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("subset", type=Path)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--first-model-year", type=int, default=1993)
    parser.add_argument("--last-model-year", type=int, default=1993)
    args = parser.parse_args()
    if args.last_model_year < args.first_model_year:
        raise SystemExit("model year range is empty")
    receipt = json.loads(args.receipt.read_text())
    if receipt["source_sha256"] != sha256(args.source) or receipt["subset_sha256"] != sha256(args.subset):
        raise SystemExit("source or subset digest differs from receipt")
    with netCDF4.Dataset(args.source) as full, netCDF4.Dataset(args.subset) as compact:
        original_time = np.asarray(full.variables["time"][:], dtype=np.float64)
        subset_time = np.asarray(compact.variables["time"][:], dtype=np.float64)
        original_year = np.asarray(full.variables["year"][:], dtype=np.int32)
        subset_year = np.asarray(compact.variables["year"][:], dtype=np.int32)
        original_month = np.asarray(full.variables["month"][:], dtype=np.int32)
        subset_month = np.asarray(compact.variables["month"][:], dtype=np.int32)
        first_subset_year = int(receipt["first_month"].split("-")[0])
        last_subset_year = int(receipt["last_month"].split("-")[0])
        indices = np.flatnonzero((original_year >= first_subset_year) &
                                 (original_year <= last_subset_year))
        if not np.array_equal(full.variables["vmro3"][indices, :, :], compact.variables["vmro3"][:, :, :]):
            raise SystemExit("subset interpolation fields differ")
        start = datetime(args.first_model_year, 1, 1)
        stop = datetime(args.last_model_year + 1, 1, 1)
        step = timedelta(minutes=15)
        now, count = start, 0
        while now <= stop:
            left = bracket(original_time, decimal_year(now))
            right = bracket(subset_time, decimal_year(now))
            for original_idx, subset_idx in zip(left[:2], right[:2]):
                if (original_year[original_idx], original_month[original_idx]) != (
                    subset_year[subset_idx], subset_month[subset_idx]
                ):
                    raise SystemExit(f"interpolation months differ at {now}")
            if np.float64(left[2]).view(np.uint64) != np.float64(right[2]).view(np.uint64):
                raise SystemExit(f"interpolation weight bits differ at {now}")
            now += step
            count += 1
    print(f"V021_OZONE_INTERPOLATION_EXACT years={args.first_model_year}-{args.last_model_year} samples={count}")


if __name__ == "__main__":
    main()
