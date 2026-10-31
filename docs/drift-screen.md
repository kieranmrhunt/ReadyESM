# Short-run drift screen

`scripts/drift_screen/excess_series.py` compares ReadyESM output with ERA5 at
matched 6-hour times and the nearest ERA5 grid point. It uses the run's own
below-ground mask for both model and reference values at 250, 500 and 850 hPa.
For each run it reports global and seven-band temperature error by full week,
by calendar month when at least 90% of the samples are present, the change
from the first to last full week, and the global error trend in K per 30 days.
Runs whose labels share the text before `:` are summarized with a member mean
and spread. The tool supports 31-day and 90-day runs.

Install the Python dependencies from `scripts/requirements.txt`. Provide the
directory containing the required 1993 ERA5 six-hourly pressure files. The
available reference set currently covers January through September, enough
for January and July 90-day runs. Then run:

```sh
python3 -I -B scripts/drift_screen/excess_series.py ERA5_6H_DIR output.json control:m1=RUN_DIR control:m2=RUN2_DIR
```

Each run directory must contain `eddy_6hourly.nc` with `time_seconds`,
temperature, pressure levels, latitude, longitude, area weights, and a
validity mask. Samples without matching ERA5 times are counted and excluded.
The output is created exclusively, so choose a new output path for each
analysis. The v0.1.18 January and July paired-month baselines are stored in
`docs/drift_screen/`.

The synthetic tests check time matching, below-ground masking, known drift,
calendar-month coverage, partial months, missing ERA5 coverage, group means,
non-overlap handling, and refusal to overwrite an output:

```sh
python3 -I -B -m unittest discover -s test/drift_screen -p 'test_*.py'
```

This screen does not replace the release's paired-month score. Treat regional
and global results as descriptive diagnostics, especially with small member
counts; it does not establish a cause for a drift change.
