# Erratum (Claude, 30 Oct 2026 19:2x UTC), outside the sealed SHA256SUMS

The docstring of `make_realistic.py` says "The broadband value then equals ERA5's realised albedo". That is too strong.
The script makes the **four-band arithmetic mean** (aluvp, aluvd, alnip, alnid) equal the ERA5 `fal` 1994-2013 monthly
climatology in each land cell, to a residual of 2.2e-16. The model's realised broadband land albedo comes from the solved
shortwave fluxes (`equivalent_land_albedo` in `src/seasonal_albedo.jl`), which weight the bands by the radiation state,
so it need not equal `fal` (Codex, 30 Oct 19:19). The sealed script and file are unchanged.
