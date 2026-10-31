"""Week-by-week and month-by-month temperature error of ReadyESM runs against ERA5 at matched times (Claude, 2026-10-30).

usage: python3 -I -B excess_series.py ERA5_6H_DIR OUTPUT.json LABEL=RUN_DIR [LABEL=RUN_DIR ...]

Each RUN_DIR holds `eddy_6hourly.nc` (6-hourly T at 250, 500 and 850 hPa, `time_seconds` since 1993-01-01). The run may
be 31 or 90 days long. ERA5_6H_DIR holds `era5_6h_pressure_1993_MM.nc` with t on the same three levels.

Each run sample is compared with ERA5 at the same time and at the nearest ERA5 grid point. This is the error ("excess")
of the run. Points that the run marks invalid (below ground) are left out of both. Samples with no ERA5 sample at the
same time are dropped and counted. The script reports, globally and in seven latitude bands:
- the mean excess in each full week from the run start (28 samples),
- the mean excess in each calendar month with at least 90 % of its samples,
- the drift: the last full week minus the first, and the least-squares trend of the global excess in K per 30 days.
Labels sharing a prefix before ':' (for example control:m1, control:m2) are grouped, with a member mean and sample
standard deviation.
"""
import calendar
import json
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, num2date

BANDS = ((-90, -60), (-60, -30), (-30, -10), (-10, 10), (10, 30), (30, 60), (60, 90))
LEVELS = (250.0, 500.0, 850.0)
WEEK = 28  # 7 days of 6-hourly samples
EPOCH = datetime(1993, 1, 1)
MIN_MONTH_FRACTION = 0.9


def band_mask(lat, lo, hi):
    return (lat >= lo) & ((lat < hi) if hi < 90 else (lat >= lo))


def band_names():
    return ["global"] + [f"{lo}..{hi}" for lo, hi in BANDS]


def hours_since_epoch(variable):
    """Whole hours since 1993-01-01 for a CF time variable."""
    times = num2date(np.asarray(variable[:]), variable.units, only_use_cftime_datetimes=False,
                     only_use_python_datetimes=True)
    out = np.array([(t - EPOCH).total_seconds() / 3600.0 for t in times])
    if not np.allclose(out, np.round(out)):
        raise ValueError("sample times are not on whole hours")
    return np.round(out).astype(np.int64)


def era5_at_points(directory, months, lat, lon):
    """ERA5 T at the nearest grid point to each run point, keyed by hour: {hour: array (levels, points)}."""
    table = {}
    for month in months:
        path = Path(directory) / f"era5_6h_pressure_1993_{month:02d}.nc"
        if not path.exists():
            continue
        with Dataset(path) as nc:
            hours = hours_since_epoch(nc["valid_time"])
            grid_lat = np.asarray(nc["latitude"][:], float)
            grid_lon = np.asarray(nc["longitude"][:], float)
            p = [float(x) for x in np.asarray(nc["pressure_level"][:], float)]
            T = np.asarray(nc["t"][:], float)[:, [p.index(level) for level in LEVELS]]  # time, level, lat, lon
        i = np.abs(grid_lat[:, None] - lat[None, :]).argmin(axis=0)
        gap = np.abs(((grid_lon[:, None] - lon[None, :]) + 180.0) % 360.0 - 180.0)
        j = gap.argmin(axis=0)
        for h, value in zip(hours, T[:, :, i, j]):
            table[int(h)] = value
    return table


def read_run(run_dir):
    with Dataset(Path(run_dir) / "eddy_6hourly.nc") as nc:
        nc.set_auto_mask(False)
        hours = hours_since_epoch(nc["time_seconds"])
        w = np.asarray(nc["area_weight"][:], float)
        lat = np.asarray(nc["latitude"][:], float)
        lon = np.asarray(nc["longitude"][:], float)
        p = [float(x) / 100.0 for x in np.asarray(nc["pressure"][:], float)]
        order = [p.index(level) for level in LEVELS]
        T = np.asarray(nc["temperature"][:], float)[:, order]
        ok = np.asarray(nc["valid"][:])[:, order] > 0
    return hours, w, lat, lon, T, ok


def area_means(difference, ok, w, lat):
    """Weighted means over valid points: (time, levels, points) to (time, bands, levels)."""
    masks = [np.ones_like(lat, bool)] + [band_mask(lat, lo, hi) for lo, hi in BANDS]
    out = np.empty((difference.shape[0], len(masks), len(LEVELS)))
    for b, mask in enumerate(masks):
        use = ok & mask[None, None, :]
        out[:, b, :] = np.where(use, difference * w, 0.0).sum(axis=2) / (w * use).sum(axis=2)
    return out


def to_bands(array):
    """(bands, levels) array to {band: [levels]}."""
    return {name: array[b].tolist() for b, name in enumerate(band_names())}


def series(era5_dir, run_dir):
    hours, w, lat, lon, T, ok = read_run(run_dir)
    if hours.size < 2 or np.any(np.diff(hours) != 6):
        raise ValueError(f"{run_dir}: samples must be 6-hourly and in order")
    first, last = (EPOCH + timedelta(hours=int(x)) for x in hours[[0, -1]])
    if not first.year == last.year == 1993:
        raise ValueError(f"{run_dir}: ERA5 references cover 1993 only")
    era = era5_at_points(era5_dir, range(first.month, last.month + 1), lat, lon)
    matched = np.array([int(h) in era for h in hours])
    if matched.sum() < 2 * WEEK:
        raise ValueError(f"{run_dir}: only {int(matched.sum())} samples match ERA5 times, need {2 * WEEK}")
    if not matched[: 2 * WEEK].all():
        raise ValueError(f"{run_dir}: the first two weeks must all match ERA5 times")
    start = int(hours[0]) - 6  # the run starts one sample before its first output
    h = hours[matched]
    # Run minus ERA5 at the same points and times, over the run's valid points only (below-ground points are
    # left out of both)
    excess = area_means(T[matched] - np.stack([era[int(x)] for x in h]), ok[matched], w, lat)  # time, band, level
    week = (h - hours[0]) // (7 * 24)  # week 0 holds the first 28 samples
    weekly = []
    for k in range(int(week.max()) + 1):
        rows = excess[week == k]
        if rows.shape[0] < WEEK:
            break
        weekly.append(rows.mean(axis=0))
    monthly = {}
    stamp = [EPOCH + timedelta(hours=int(x)) for x in h]
    for key in sorted({(s.year, s.month) for s in stamp}):
        rows = excess[[(s.year, s.month) == key for s in stamp]]
        expected = 4 * calendar.monthrange(*key)[1]
        if rows.shape[0] >= MIN_MONTH_FRACTION * expected:
            monthly[f"{key[0]}-{key[1]:02d}"] = {"samples": int(rows.shape[0]), **to_bands(rows.mean(axis=0))}
    days = (h - h[0]) / 24.0
    trend = np.polyfit(days, excess[:, 0, :], 1)[0] * 30.0  # global, K per 30 days, per level
    return {"run_dir": str(run_dir), "start_hour": start, "samples": int(hours.size),
            "unmatched_samples": int((~matched).sum()), "weeks": len(weekly),
            "weekly": [to_bands(x) for x in weekly], "monthly": monthly,
            "drift_last_minus_first_week": to_bands(weekly[-1] - weekly[0]),
            "global_trend_k_per_30_days": trend.tolist()}


def main(era5_dir, output, runs):
    result = {"schema": "readiesm_excess_series_v1", "levels_hpa": list(LEVELS), "week_samples": WEEK,
              "min_month_fraction": MIN_MONTH_FRACTION, "runs": {}, "groups": {}}
    groups = {}
    for spec in runs:
        label, path = spec.split("=", 1)
        result["runs"][label] = series(era5_dir, path)
        groups.setdefault(label.split(":")[0], []).append(result["runs"][label])
    for name, members in groups.items():
        weeks = min(m["weeks"] for m in members)
        weekly = np.array([[m["weekly"][k]["global"] for k in range(weeks)] for m in members])
        drift = np.array([m["drift_last_minus_first_week"]["global"] for m in members])
        trend = np.array([m["global_trend_k_per_30_days"] for m in members])
        months = sorted(set.intersection(*(set(m["monthly"]) for m in members)))
        sd = (lambda a: a.std(axis=0, ddof=1) if len(members) > 1 else np.zeros(a.shape[1:]))
        result["groups"][name] = {
            "members": len(members), "weeks": weeks,
            "weekly_global_mean": weekly.mean(axis=0).tolist(), "weekly_global_sd": sd(weekly).tolist(),
            "drift_global_mean": drift.mean(axis=0).tolist(), "drift_global_sd": sd(drift).tolist(),
            "trend_global_mean": trend.mean(axis=0).tolist(), "trend_global_sd": sd(trend).tolist(),
            "monthly_global_mean": {k: np.mean([m["monthly"][k]["global"] for m in members], axis=0).tolist()
                                    for k in months}}
    with open(output, "x") as stream:
        json.dump(result, stream, indent=1)
    for name, g in result["groups"].items():
        print(f"{name}: n={g['members']} weeks={g['weeks']} global drift (250/500/850) "
              + " ".join(f"{x:+.2f}" for x in g["drift_global_mean"]) + "  sd "
              + " ".join(f"{x:.2f}" for x in g["drift_global_sd"]) + "  trend K/30d "
              + " ".join(f"{x:+.2f}" for x in g["trend_global_mean"]))
        for month, value in g["monthly_global_mean"].items():
            print(f"  {month} global error " + " ".join(f"{x:+.2f}" for x in value))
    return result


if __name__ == "__main__":
    if len(sys.argv) < 4:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
