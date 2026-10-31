"""Tests for excess_series.py with synthetic ERA5 and run files of known error (Claude, 2026-10-30)."""
import calendar
import importlib.util
import tempfile
import unittest
from datetime import datetime
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "excess_series", HERE.parent.parent / "scripts" / "drift_screen" / "excess_series.py"
)
X = importlib.util.module_from_spec(spec)
spec.loader.exec_module(X)

LAT = np.linspace(90.0, -90.0, 73)
LON = np.arange(0.0, 360.0, 45.0)
SINCE_1970 = (datetime(1993, 1, 1) - datetime(1970, 1, 1)).total_seconds()


def era5_signal(hours):
    """A strong 6-hourly signal, so any time misalignment shows as error."""
    return 250.0 + 3.0 * np.sin(2 * np.pi * hours / 30.0) + 0.01 * hours / 24.0


def write_era5(root, months=(1, 2, 3, 4)):
    for month in months:
        start = (datetime(1993, month, 1) - datetime(1993, 1, 1)).total_seconds() / 3600
        hours = start + 6 * np.arange(4 * calendar.monthrange(1993, month)[1])
        with Dataset(root / f"era5_6h_pressure_1993_{month:02d}.nc", "w") as nc:
            nc.createDimension("valid_time", hours.size)
            nc.createDimension("pressure_level", 3)
            nc.createDimension("latitude", LAT.size)
            nc.createDimension("longitude", LON.size)
            t = nc.createVariable("valid_time", "i8", ("valid_time",))
            t.units = "seconds since 1970-01-01"
            t[:] = (SINCE_1970 + 3600 * hours).astype(np.int64)
            nc.createVariable("pressure_level", "f8", ("pressure_level",))[:] = [850.0, 500.0, 250.0]
            nc.createVariable("latitude", "f8", ("latitude",))[:] = LAT
            nc.createVariable("longitude", "f8", ("longitude",))[:] = LON
            field = era5_signal(hours)[:, None, None, None] * np.ones((hours.size, 3, LAT.size, LON.size))
            nc.createVariable("t", "f4", ("valid_time", "pressure_level", "latitude", "longitude"))[:] = field


def write_run(path, days=90, rate=0.0, offset=0.0, start_hour=0, points=300):
    """Run file whose T is ERA5 plus offset plus rate per day since the run start."""
    rng = np.random.default_rng(1)
    lat = rng.uniform(-89, 89, points)
    hours = start_hour + 6 * np.arange(1, 4 * days + 1)
    path.mkdir(parents=True, exist_ok=True)
    with Dataset(path / "eddy_6hourly.nc", "w") as nc:
        nc.createDimension("time", hours.size)
        nc.createDimension("pressure", 3)
        nc.createDimension("point", points)
        t = nc.createVariable("time_seconds", "f8", ("time",))
        t.units = "seconds since 1993-01-01 00:00:00"
        t[:] = 3600.0 * hours
        nc.createVariable("area_weight", "f8", ("point",))[:] = np.cos(np.deg2rad(lat))
        nc.createVariable("latitude", "f8", ("point",))[:] = lat
        nc.createVariable("longitude", "f8", ("point",))[:] = rng.uniform(0, 360, points)
        nc.createVariable("pressure", "f8", ("pressure",))[:] = [25000.0, 50000.0, 85000.0]
        value = era5_signal(hours) + offset + rate * (hours - start_hour) / 24.0
        nc.createVariable("temperature", "f8", ("time", "pressure", "point"))[:] = (
            value[:, None, None] * np.ones((hours.size, 3, points)))
        nc.createVariable("valid", "i1", ("time", "pressure", "point"))[:] = 1


class ExcessSeriesTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        write_era5(self.root)

    def tearDown(self):
        self.tmp.cleanup()

    def run_main(self, *runs):
        return X.main(self.root, self.root / "out.json", [f"{label}={self.root / label.replace(':', '_')}"
                                                         for label in runs])

    def test_matched_times_give_zero_error(self):
        # The run equals ERA5 at the same times. A one-sample shift would give errors of order 1 K here.
        write_run(self.root / "x_m1")
        r = self.run_main("x:m1")["runs"]["x:m1"]
        self.assertEqual(r["unmatched_samples"], 0)
        self.assertEqual(r["weeks"], 12)
        for week in r["weekly"]:
            for band in week.values():
                np.testing.assert_allclose(band, 0.0, atol=1e-4)

    def test_points_and_mask(self):
        # ERA5 varies with latitude and longitude, and the run sits on ERA5 grid points. A fifth of the 850 hPa points
        # are below ground: NaN and invalid in the run. The error must still be zero.
        for path in sorted(self.root.glob("era5_*.nc")):
            with Dataset(path, "a") as nc:
                pattern = 0.1 * LAT[:, None] + np.cos(np.deg2rad(LON))[None, :]
                nc["t"][:] = nc["t"][:] + pattern[None, None]
        write_run(self.root / "x_m1")
        rng = np.random.default_rng(2)
        with Dataset(self.root / "x_m1" / "eddy_6hourly.nc", "a") as nc:
            i = rng.integers(0, LAT.size, 300)
            j = rng.integers(0, LON.size, 300)
            nc["latitude"][:] = LAT[i]
            nc["longitude"][:] = LON[j]
            T = nc["temperature"][:] + (0.1 * LAT[i] + np.cos(np.deg2rad(LON[j])))[None, None, :]
            below = rng.random(300) < 0.2
            T[:, 2, below] = np.nan
            nc["temperature"][:] = T
            valid = np.ones(T.shape, np.int8)
            valid[:, 2, below] = 0
            nc["valid"][:] = valid
        r = self.run_main("x:m1")["runs"]["x:m1"]
        for week in r["weekly"]:
            for band in week.values():
                np.testing.assert_allclose(band, 0.0, atol=2e-4)

    def test_known_cooling(self):
        # -0.05 K/day: weeks 12 and 1 are 77 days apart, and the trend is -1.5 K per 30 days.
        write_run(self.root / "x_m1", rate=-0.05)
        r = self.run_main("x:m1")["runs"]["x:m1"]
        np.testing.assert_allclose(r["drift_last_minus_first_week"]["global"], -0.05 * 77, atol=1e-4)
        np.testing.assert_allclose(r["drift_last_minus_first_week"]["30..60"], -0.05 * 77, atol=1e-4)
        np.testing.assert_allclose(r["global_trend_k_per_30_days"], -1.5, atol=1e-4)

    def test_calendar_months(self):
        # A 90-day January start fills January, February and March exactly. The January mean error of a
        # constant offset is the offset.
        write_run(self.root / "x_m1", offset=0.7)
        r = self.run_main("x:m1")["runs"]["x:m1"]
        self.assertEqual(sorted(r["monthly"]), ["1993-01", "1993-02", "1993-03"])
        self.assertEqual(r["monthly"]["1993-02"]["samples"], 112)
        np.testing.assert_allclose(r["monthly"]["1993-01"]["-10..10"], 0.7, atol=1e-4)

    def test_partial_month_dropped(self):
        # A 40-day run has all of January and 9 days of February: February is below 90 % and is left out.
        write_run(self.root / "x_m1", days=40)
        r = self.run_main("x:m1")["runs"]["x:m1"]
        self.assertEqual(sorted(r["monthly"]), ["1993-01"])

    def test_samples_beyond_era5_are_counted(self):
        # Without the April ERA5 file, a 40-day run from 1 March keeps its 123 March samples and loses 37.
        (self.root / "era5_6h_pressure_1993_04.nc").unlink()
        start = int((datetime(1993, 3, 1) - datetime(1993, 1, 1)).total_seconds() // 3600)
        write_run(self.root / "x_m1", days=40, start_hour=start)
        r = self.run_main("x:m1")["runs"]["x:m1"]
        self.assertEqual(r["unmatched_samples"], 37)
        self.assertEqual(r["weeks"], 4)

    def test_group_mean_and_spread(self):
        write_run(self.root / "c_m1", rate=-0.05)
        write_run(self.root / "c_m2", rate=-0.10)
        g = self.run_main("c:m1", "c:m2")["groups"]["c"]
        self.assertEqual(g["members"], 2)
        np.testing.assert_allclose(g["drift_global_mean"], -0.075 * 77, atol=1e-4)
        np.testing.assert_allclose(g["drift_global_sd"], np.std([-0.05 * 77, -0.10 * 77], ddof=1), atol=1e-4)

    def test_no_overlap_refused(self):
        start = int((datetime(1993, 6, 1) - datetime(1993, 1, 1)).total_seconds() // 3600)
        write_run(self.root / "x_m1", days=31, start_hour=start)
        with self.assertRaises(ValueError):
            self.run_main("x:m1")

    def test_existing_output_refused(self):
        write_run(self.root / "x_m1")
        (self.root / "out.json").write_text("{}")
        with self.assertRaises(FileExistsError):
            self.run_main("x:m1")


if __name__ == "__main__":
    unittest.main()
