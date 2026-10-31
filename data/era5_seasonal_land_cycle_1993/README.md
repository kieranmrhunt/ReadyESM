# ERA5 seasonal land-cycle input

This directory holds the twelve monthly 1993 ERA5 single-level fields used by
seasonal land-cycle experiments: high and low vegetation cover, high and low
leaf area index, four snow-free direct/diffuse visible and near-infrared
albedo components, and land-sea mask. The source file is 1-degree global
monthly data.

The file SHA-256 is recorded in `SHA256SUMS` and `receipt.json`. The adjacent
request files preserve the CDS dataset name, variables, grid, and request ID
`3a99eb70-87d6-46bc-87d8-af09998ace28`. The acquisition receipt retains the
original source SHA and its validated dimensions, units, and ranges; its file
path is relative to this directory. `download.py` can reproduce the request
in this directory when CDS credentials and Python dependencies are available.

The four albedo fields remain separate. Their conversion to the model's
radiation boundary is not defined by this source-data bundle and must be
validated by the seasonal-albedo implementation.
