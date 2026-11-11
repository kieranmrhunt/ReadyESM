#!/usr/bin/env python3
"""Same-date zonal 5/1 hPa temperatures for a provisional January cap experiment."""
import csv
import hashlib
import sys
from pathlib import Path
import numpy as np
import xarray as xr

SOURCE_SHA256 = 'fb3a51a140a665fddbba724079ef8b2fdf9587ebcd89c81027e53872fea354c7'
source, output = map(Path, sys.argv[1:])
assert hashlib.sha256(source.read_bytes()).hexdigest() == SOURCE_SHA256
with xr.open_dataset(source) as ds:
    assert ds['t'].dims == ('valid_time', 'pressure_level', 'latitude', 'longitude')
    assert dict(ds.sizes) == {'valid_time': 1, 'pressure_level': 37, 'latitude': 181, 'longitude': 360}
    lat = ds.latitude.values.astype(float)
    assert np.array_equal(lat, np.arange(90, -91, -1))
    five = ds['t'].sel(pressure_level=5).isel(valid_time=0).values.mean(axis=1)
    one = ds['t'].sel(pressure_level=1).isel(valid_time=0).values.mean(axis=1)
    assert np.isfinite(five).all() and np.isfinite(one).all()
    assert 180 < five.min() < five.max() < 330
    assert 180 < one.min() < one.max() < 330
with output.open('w', newline='') as f:
    writer = csv.writer(f, lineterminator='\n')
    writer.writerow(('latitude_degrees', 'temperature_5hpa_k', 'temperature_1hpa_k'))
    for la, t5, t1 in zip(lat[::-1], five[::-1], one[::-1]):
        writer.writerow((f'{la:.0f}', f'{t5:.9f}', f'{t1:.9f}'))
print('UPPER_CAP_ERA5_19930101_PREPARE_PASS')
