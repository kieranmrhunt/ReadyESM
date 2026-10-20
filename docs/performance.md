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

The rate uses total simulated time divided by total wall time across both
processes for each grid, with 365 days per year. Individual daily windows span
14.507–16.103 simulated years per day for the 180×90 ocean and
10.695–11.579 for the 360×180 ocean.

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

Slurm jobs `60854483_0` and `60854483_1` completed with exit code 0. All report
checks, finite-field checks and result receipts pass. Each pair used one GPU.
The benchmark covers short January integrations.
