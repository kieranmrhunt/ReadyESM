# Model speed

Measured on 20 October 2026 using model source `e9e50b4`, Julia 1.12.5,
one NVIDIA A100-SXM4-40GB and four CPU threads on JASMIN.

The test retains production physics, including radiation, dynamic sea ice
and Terrarium. Each fresh process runs eight warm-up steps, then times three
model days. Daily budget diagnostics are included. GPU work is synchronized
at the timing boundaries.

| Atmosphere | Ocean cells | Coupled step | Radiation interval | Years per day |
| --- | --- | ---: | ---: | ---: |
| T15, 27 levels | 180×90×60 | 30 min | 2 h | 15.100 |
| T15, 27 levels | 360×180×60 | 30 min | 2 h | 11.134 |
| T31, 27 levels | 180×90×60 | 15 min | 1 h | 8.175 |
| T31, 27 levels | 720×360×60 | 15 min | 1 h | 2.380 |

The rate uses total simulated time divided by total wall time across both
processes for each grid, with 365 days per year. Individual daily windows span
14.507–16.103 simulated years per day for T15 with the 180×90 ocean,
10.695–11.579 for T15 with the 360×180 ocean, and 8.154–8.193 for T31
with the 180×90 ocean. T31 with the 720×360 ocean spans 2.364–2.394.

The following phase times are for T15.

| Phase | 180×90 first | 180×90 second | 360×180 first | 360×180 second |
| --- | ---: | ---: | ---: | ---: |
| Package loading | 372.0 s | 14.2 s | 357.1 s | 13.4 s |
| Input verification | 7.4 s | 6.9 s | 6.9 s | 6.9 s |
| Model construction | 1841.2 s | 1802.2 s | 1792.1 s | 1804.7 s |
| Initialization | 105.7 s | 106.7 s | 106.3 s | 106.5 s |
| Eight warm-up steps | 795.6 s | 782.5 s | 787.1 s | 784.5 s |
| Three timed model days | 48.3 s | 45.8 s | 61.4 s | 66.2 s |
| Final output | 62.3 s | 62.7 s | 61.1 s | 62.7 s |

Total start-up was 52.0 and 45.2 minutes for the 180×90 ocean, and 50.8 and
45.3 minutes for the 360×180 ocean. Final output is timed separately.
Each second process reuses the package cache. All processes still construct
and compile their own model.

For T31 with the 180×90 ocean, start-up takes 52.3 and 45.5 minutes.
With the 720×360 ocean, it takes 53.6 and 48.2 minutes. The phase times are:

| Phase | 180×90 first | 180×90 second | 720×360 first | 720×360 second |
| --- | ---: | ---: | ---: | ---: |
| Package loading | 364.5 s | 14.3 s | 365.0 s | 14.6 s |
| Input verification | 6.9 s | 7.1 s | 6.9 s | 7.0 s |
| Model construction | 1872.1 s | 1808.6 s | 1932.8 s | 1955.3 s |
| Initialization | 106.2 s | 107.1 s | 108.4 s | 108.4 s |
| Eight warm-up steps | 787.4 s | 791.1 s | 803.5 s | 808.0 s |
| Three timed model days | 86.9 s | 86.8 s | 299.9 s | 296.9 s |
| Final output | 62.4 s | 61.2 s | 61.2 s | 62.4 s |

Slurm jobs `60854483_0`, `_1`, `_3` and `_5` completed with exit code 0. All report
checks, finite-field checks and result receipts pass. Each pair used one GPU.
The benchmark covers short January integrations.
