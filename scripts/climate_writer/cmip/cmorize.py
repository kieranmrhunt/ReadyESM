"""CMOR-style NetCDF from the sealed monthly files of the free-running climate run (Claude, 2026-10-29).

Offline analysis code, not model code. It writes one CF-1.8 file per variable with CMIP short names, so that
ClimateEval (commit d2559b2f) and ClimateBench2 (commit 3ff9a3cc) can read ReadyESM output directly. The files are
"CMOR-style": they follow the CMIP names, units and coordinates, but the CMOR library did not write them, and they
are not a CMIP submission.

Input: the 48 monthly files of segment B (job 62507908), January 1993 to December 1996, each checked against that
segment's result manifest. v0.1.15 physics, historical gases and ozone, observed start in January 1993.

Grids:
  - Atmosphere (grid label gn): the native T31 full Gaussian grid, 48 rings of 96 points, written as a rectilinear
    latitude-longitude grid with no interpolation. Latitude bounds come from the model's own ring weights, so an
    area-weighted mean on the output grid equals the model's quadrature.
  - Ocean and sea ice (grid label gr): the native 1 degree tripolar grid has no regular rows, so the fields are
    remapped to a regular 1 degree grid. Each target cell is sampled at 5 x 5 points. Each sample takes the value of
    the nearest native cell centre (a Voronoi approximation of the native cells, which handles the north fold by
    geometry). The target value is the cos(latitude)-weighted mean over the samples that land on wet cells. sftof is
    the wet share of the target cell. tos and siconc are masked where sftof is below 50 percent, which is the usual
    land-sea mask choice. The duplicate fold cells that the writer excludes from area sums are also excluded here.

Variables and sources (monthly means):
  rsdt rsut rlut rsutcs rlutcs   toa_incoming_shortwave, toa_outgoing_shortwave, toa_outgoing_longwave,
                                 toa_clear_outgoing_shortwave, toa_clear_outgoing_longwave
  rsds rsus rlds rlus            surface_shortwave_down/up, surface_longwave_down/up (grid-box blend)
  hfss                           surface_sensible_heat_upward
  evspsbl                        surface_water_vapor_upward (kg m-2 s-1)
  hfls                           derived: L_v x surface_water_vapor_upward, L_v = 2.5008e6 J kg-1 (the ClimateBench
                                 E = hfls / L_v convention, so the two stay consistent)
  pr, prsn                       total_precipitation_rate and snow_rate x 1000 kg m-3
  ps, prw                        surface_pressure, column_water
  ta (200, 500, 850 hPa)         temperature_*hpa_contribution, masked unless *_valid_fraction == 1
  ua (200, 850 hPa)              zonal_wind_*hpa_contribution, masked unless *_valid_fraction == 1
  tos, siconc                    sst (degC), ice_concentration x 100 (percent of the ocean part of the cell)
  sftof (fixed)                  wet share of each target cell, percent
Not available from this run: tas, ts, clt, clwvi, clivi, hus, zg, daily fields, and ocean heat content. ClimateEval's
atmosphere suites interpolate ta and ua linearly to 200, 400 and 800 hPa. With only 200, 500 and 850 hPa written here,
its 400 and 800 hPa values would be interpolated across wide gaps.

Compatibility helpers for cmorize_v2.py. Use cmorize_v2.py as the supported command-line interface.
"""
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from netCDF4 import Dataset, date2num
from scipy.spatial import cKDTree

JOB = None
YEARS = (1993, 1994, 1995, 1996)
SOURCE_ID = "ReadyESM-v0-1-15"
EXPERIMENT = "historical"
VARIANT = "r1i1p1f1"
L_V = 2.5008e6
RHO_W = 1000.0
EARTH_RADIUS_M = 6.371e6
FILL = np.float32(1.0e20)
TIME_UNITS = "days since 1993-01-01 00:00:00"
CALENDAR = "proleptic_gregorian"
SAMPLES = 5            # per target-cell side
WET_THRESHOLD = 50.0   # percent
OCEAN_RES = 1.0        # degrees

# cmor name: (source variable, source units, scale, units, standard_name, long_name, comment)
ATMOS = {
    "rsdt": ("toa_incoming_shortwave", "W m-2", 1.0, "W m-2", "toa_incoming_shortwave_flux",
             "TOA Incident Shortwave Radiation", ""),
    "rsut": ("toa_outgoing_shortwave", "W m-2", 1.0, "W m-2", "toa_outgoing_shortwave_flux",
             "TOA Outgoing Shortwave Radiation", ""),
    "rlut": ("toa_outgoing_longwave", "W m-2", 1.0, "W m-2", "toa_outgoing_longwave_flux",
             "TOA Outgoing Longwave Radiation", "latest post-physics radiation solve"),
    "rsutcs": ("toa_clear_outgoing_shortwave", "W m-2", 1.0, "W m-2",
               "toa_outgoing_shortwave_flux_assuming_clear_sky", "TOA Outgoing Clear-Sky Shortwave Radiation", ""),
    "rlutcs": ("toa_clear_outgoing_longwave", "W m-2", 1.0, "W m-2",
               "toa_outgoing_longwave_flux_assuming_clear_sky", "TOA Outgoing Clear-Sky Longwave Radiation", ""),
    "rsds": ("surface_shortwave_down", "W m-2", 1.0, "W m-2", "surface_downwelling_shortwave_flux_in_air",
             "Surface Downwelling Shortwave Radiation", "grid-box blend of land, ocean and ice surfaces"),
    "rsus": ("surface_shortwave_up", "W m-2", 1.0, "W m-2", "surface_upwelling_shortwave_flux_in_air",
             "Surface Upwelling Shortwave Radiation", "grid-box blend of land, ocean and ice surfaces"),
    "rlds": ("surface_longwave_down", "W m-2", 1.0, "W m-2", "surface_downwelling_longwave_flux_in_air",
             "Surface Downwelling Longwave Radiation", "grid-box blend of land, ocean and ice surfaces"),
    "rlus": ("surface_longwave_up", "W m-2", 1.0, "W m-2", "surface_upwelling_longwave_flux_in_air",
             "Surface Upwelling Longwave Radiation", "grid-box blend of land, ocean and ice surfaces"),
    "hfss": ("surface_sensible_heat_upward", "W m-2", 1.0, "W m-2", "surface_upward_sensible_heat_flux",
             "Surface Upward Sensible Heat Flux", ""),
    "evspsbl": ("surface_water_vapor_upward", "kg m-2 s-1", 1.0, "kg m-2 s-1", "water_evapotranspiration_flux",
                "Evaporation Including Sublimation and Transpiration", ""),
    "hfls": ("surface_water_vapor_upward", "kg m-2 s-1", L_V, "W m-2", "surface_upward_latent_heat_flux",
             "Surface Upward Latent Heat Flux",
             "derived as 2.5008e6 J kg-1 times evspsbl, the ClimateBench convention E = hfls / L_v; "
             "sublimation is counted at the vaporisation constant"),
    "pr": ("total_precipitation_rate", "m s-1", RHO_W, "kg m-2 s-1", "precipitation_flux",
           "Precipitation", "rain plus snow, liquid water equivalent times 1000 kg m-3"),
    "prsn": ("snow_rate", "m s-1", RHO_W, "kg m-2 s-1", "snowfall_flux", "Snowfall Flux",
             "liquid water equivalent times 1000 kg m-3"),
    "ps": ("surface_pressure", "Pa", 1.0, "Pa", "surface_air_pressure", "Surface Air Pressure", ""),
    "prw": ("column_water", "kg m-2", 1.0, "kg m-2", "atmosphere_mass_content_of_water_vapor", "Water Vapor Path",
            "from specific humidity and sigma-layer mass"),
}
# cmor name: (source prefix, pressures in hPa, source units, units, standard_name, long_name)
PROFILES = {
    "ta": ("temperature", (850, 500, 200), "K", "K", "air_temperature", "Air Temperature"),
    "ua": ("zonal_wind", (850, 200), "m s-1", "m s-1", "eastward_wind", "Eastward Wind"),
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def sealed_files():
    seal = {}
    for line in (JOB / "result.sha256").read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        seal[name.removeprefix("./")] = digest
    files = []
    for year in YEARS:
        for month in range(1, 13):
            rel = f"run/monthly_{year}_{month:02d}.nc"
            digest = sha256(JOB / rel)
            if seal.get(rel) != digest:
                raise SystemExit(f"{rel} does not match the segment B result manifest")
            files.append((year, month, JOB / rel, digest))
    return files


def read_var(ds, name, units):
    v = ds[name]
    if v.units != units:
        raise SystemExit(f"{name} has units {v.units!r}, expected {units!r}")
    a = np.asarray(np.ma.filled(v[:].astype(float), np.nan), float)
    return a


# ---------------------------------------------------------------- atmosphere grid

def atmosphere_grid(ds):
    lon = np.asarray(ds["atmosphere_longitude"][:], float)
    lat = np.asarray(ds["atmosphere_latitude"][:], float)
    w = np.asarray(ds["atmosphere_area_weight"][:], float)
    lats = np.unique(lat)
    lons = np.unique(lon)
    if lats.size * lons.size != lat.size:
        raise SystemExit("atmosphere grid is not a full rectilinear Gaussian grid")
    j = np.searchsorted(lats, lat)
    i = np.searchsorted(lons, lon)
    if np.unique(j * lons.size + i).size != lat.size:
        raise SystemExit("atmosphere points do not map one to one onto the rectilinear grid")
    ring = np.zeros(lats.size)
    np.add.at(ring, j, w)
    ring /= ring.sum()
    if np.ptp(w.reshape(-1)[np.argsort(j * lons.size + i)].reshape(lats.size, lons.size), axis=1).max() > 1e-12:
        raise SystemExit("area weights vary along a ring")
    edges = np.concatenate([[-1.0], -1.0 + 2.0 * np.cumsum(ring)])
    edges[-1] = 1.0
    lat_bnds = np.rad2deg(np.arcsin(np.clip(edges, -1, 1)))
    lat_bnds = np.stack([lat_bnds[:-1], lat_bnds[1:]], axis=1)
    if not np.all((lat_bnds[:, 0] < lats) & (lats < lat_bnds[:, 1])):
        raise SystemExit("ring-weight latitude bounds do not enclose the ring latitudes")
    dlon = 360.0 / lons.size
    lon_bnds = np.stack([lons - dlon / 2, lons + dlon / 2], axis=1)
    return {"lat": lats, "lon": lons, "lat_bnds": lat_bnds, "lon_bnds": lon_bnds, "j": j, "i": i,
            "point_lon": lon, "point_lat": lat, "weight": w}


def to_grid(values, g):
    out = np.full((g["lat"].size, g["lon"].size), np.nan)
    out[g["j"], g["i"]] = values
    return out


def cell_area(lat_bnds, lon_bnds):
    s = np.sin(np.deg2rad(lat_bnds))
    band = (s[:, 1] - s[:, 0])[:, None]
    dl = np.deg2rad(lon_bnds[:, 1] - lon_bnds[:, 0])[None, :]
    return EARTH_RADIUS_M ** 2 * band * dl


# ---------------------------------------------------------------- ocean remap

def unit_vectors(lon, lat):
    lo, la = np.deg2rad(lon), np.deg2rad(lat)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def ocean_remap(ds):
    olon = np.asarray(ds["ocean_longitude"][:], float)
    olat = np.asarray(ds["ocean_latitude"][:], float)
    active = np.asarray(ds["ocean_active_mask"][:]) == 1
    area = np.asarray(ds["ocean_cell_area"][:], float)
    xyz = unit_vectors(olon, olat).reshape(-1, 3)
    flat_active = active.reshape(-1)
    # An inactive cell is a fold duplicate when an active cell sits at the same place, otherwise it is land.
    tree_active = cKDTree(xyz[flat_active])
    dist, _ = tree_active.query(xyz[~flat_active])
    duplicate = np.zeros(flat_active.size, bool)
    duplicate[np.where(~flat_active)[0][dist < 1e-7]] = True
    land = ~flat_active & ~duplicate
    candidates = np.where(flat_active | land)[0]
    tree = cKDTree(xyz[candidates])

    lat_c = np.arange(-90 + OCEAN_RES / 2, 90, OCEAN_RES)
    lon_c = np.arange(OCEAN_RES / 2, 360, OCEAN_RES)
    offs = (np.arange(SAMPLES) + 0.5) / SAMPLES - 0.5
    sub_lat = (lat_c[:, None] + offs[None, :] * OCEAN_RES)          # [ny, s]
    sub_lon = (lon_c[:, None] + offs[None, :] * OCEAN_RES)          # [nx, s]
    ny, nx = lat_c.size, lon_c.size
    LA = np.broadcast_to(sub_lat[:, None, :, None], (ny, nx, SAMPLES, SAMPLES))
    LO = np.broadcast_to(sub_lon[None, :, None, :], (ny, nx, SAMPLES, SAMPLES))
    wt = np.cos(np.deg2rad(LA))
    d, k = tree.query(unit_vectors(LO, LA).reshape(-1, 3))
    src = candidates[k]
    south_limit = olat.min() - 0.5 * OCEAN_RES
    wet = flat_active[src] & (LA.reshape(-1) >= south_limit)
    max_wet_chord = float(d[wet].max())
    target = np.repeat(np.arange(ny * nx), SAMPLES * SAMPLES)
    w_all = wt.reshape(-1)
    total = np.bincount(target, weights=w_all, minlength=ny * nx)
    wet_w = np.bincount(target, weights=w_all * wet, minlength=ny * nx)
    sftof = 100.0 * wet_w / total
    # Rows of the sparse operator: wet samples of each target cell, normalised to sum to one.
    active_index = -np.ones(flat_active.size, int)
    active_index[flat_active] = np.arange(flat_active.sum())
    rows, cols, vals = target[wet], active_index[src[wet]], w_all[wet] / wet_w[target[wet]]
    op = sp.csr_matrix((vals, (rows, cols)), shape=(ny * nx, int(flat_active.sum())))
    lat_bnds = np.stack([lat_c - OCEAN_RES / 2, lat_c + OCEAN_RES / 2], axis=1)
    lon_bnds = np.stack([lon_c - OCEAN_RES / 2, lon_c + OCEAN_RES / 2], axis=1)
    geom = {"lon": olon, "lat": olat, "active": active, "area": area}
    return {"op": op, "sftof": sftof.reshape(ny, nx), "lat": lat_c, "lon": lon_c, "lat_bnds": lat_bnds,
            "lon_bnds": lon_bnds, "active": active, "area": area, "geom": geom,
            "info": {"duplicates": int(duplicate.sum()), "land": int(land.sum()), "active": int(flat_active.sum()),
                     "samples_per_cell": SAMPLES * SAMPLES, "max_wet_sample_chord": max_wet_chord,
                     "max_wet_sample_distance_km": max_wet_chord * EARTH_RADIUS_M / 1e3}}


def remap(field, r):
    return (r["op"] @ field.reshape(-1)[r["active"].reshape(-1)]).reshape(r["sftof"].shape)


# ---------------------------------------------------------------- writing

def time_axis(months):
    starts = [datetime(y, m, 1) for y, m in months]
    ends = [datetime(y + m // 12, m % 12 + 1, 1) for y, m in months]
    b = np.stack([date2num(starts, TIME_UNITS, CALENDAR), date2num(ends, TIME_UNITS, CALENDAR)], axis=1)
    return 0.5 * (b[:, 0] + b[:, 1]), b


def write(path, name, data, attrs, grid, table, realm, grid_label, common, times=None, plev=None, cell_methods=""):
    with Dataset(path, "w", format="NETCDF4_CLASSIC") as nc:
        nc.createDimension("bnds", 2)
        dims = []
        if times is not None:
            nc.createDimension("time", times[0].size)
            t = nc.createVariable("time", "f8", ("time",))
            t.setncatts({"units": TIME_UNITS, "calendar": CALENDAR, "standard_name": "time", "long_name": "time",
                         "axis": "T", "bounds": "time_bnds"})
            t[:] = times[0]
            tb = nc.createVariable("time_bnds", "f8", ("time", "bnds"))
            tb[:] = times[1]
            dims.append("time")
        if plev is not None:
            nc.createDimension("plev", len(plev))
            p = nc.createVariable("plev", "f8", ("plev",))
            p.setncatts({"units": "Pa", "standard_name": "air_pressure", "long_name": "pressure", "axis": "Z",
                         "positive": "down"})
            p[:] = np.asarray(plev, float) * 100.0
            dims.append("plev")
        for axis, std, units, ax in (("lat", "latitude", "degrees_north", "Y"),
                                     ("lon", "longitude", "degrees_east", "X")):
            nc.createDimension(axis, grid[axis].size)
            c = nc.createVariable(axis, "f8", (axis,))
            c.setncatts({"units": units, "standard_name": std, "long_name": std, "axis": ax,
                         "bounds": f"{axis}_bnds"})
            c[:] = grid[axis]
            cb = nc.createVariable(f"{axis}_bnds", "f8", (axis, "bnds"))
            cb[:] = grid[f"{axis}_bnds"]
        dims += ["lat", "lon"]
        v = nc.createVariable(name, "f4", tuple(dims), zlib=True, complevel=1, fill_value=FILL)
        v.setncatts({k: val for k, val in attrs.items() if val != ""})
        v.missing_value = FILL
        if cell_methods:
            v.cell_methods = cell_methods
        out = np.where(np.isfinite(data), data, FILL).astype(np.float32)
        v[:] = np.ma.masked_equal(out, FILL)
        nc.setncatts({**common, "table_id": table, "realm": realm, "grid_label": grid_label,
                      "frequency": "mon" if times is not None else "fx", "variable_id": name,
                      "grid": grid["description"], "nominal_resolution": grid["nominal_resolution"]})


def filename(name, table, grid_label, timed=True):
    tail = "_199301-199612" if timed else ""
    return f"{name}_{table}_{SOURCE_ID}_{EXPERIMENT}_{VARIANT}_{grid_label}{tail}.nc"


# ---------------------------------------------------------------- main

def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    out_root = Path(sys.argv[1])
    if out_root.exists():
        raise SystemExit(f"{out_root} exists; refusing to overwrite")
    files = sealed_files()
    months = [(y, m) for y, m, _, _ in files]
    times = time_axis(months)

    atm = {n: [] for n in ATMOS}
    prof = {n: [] for n in PROFILES}
    tos, sic = [], []
    native = {"sst_mean": [], "nh_ice_area": [], "sh_ice_area": [], "nh_ice_extent": [], "sh_ice_extent": [],
              "toa_net_file": []}
    with Dataset(files[0][2]) as ds:
        g = atmosphere_grid(ds)
        r = ocean_remap(ds)
        source_attrs = {k: ds.getncattr(k) for k in ("source_commit", "configuration_sha256", "writer_sha256",
                                                     "writer_schema")}
    for year, month, path, _ in files:
        with Dataset(path) as ds:
            if int(ds.getncattr("complete")) != 1:
                raise SystemExit(f"{path} is not a complete month")
            if not (np.array_equal(np.asarray(ds["atmosphere_latitude"][:], float), g["point_lat"]) and
                    np.array_equal(np.asarray(ds["atmosphere_longitude"][:], float), g["point_lon"]) and
                    np.array_equal(np.asarray(ds["ocean_latitude"][:], float), r["geom"]["lat"]) and
                    np.array_equal(np.asarray(ds["ocean_active_mask"][:]) == 1, r["active"])):
                raise SystemExit(f"{path} geometry differs from the first month")
            for n, (src, units, scale, *_rest) in ATMOS.items():
                atm[n].append(to_grid(read_var(ds, src, units) * scale, g))
            for n, (prefix, levels, units, *_rest) in PROFILES.items():
                levs = []
                for p in levels:
                    val = read_var(ds, f"{prefix}_{p}hpa_contribution", units)
                    frac = read_var(ds, f"{prefix}_{p}hpa_valid_fraction", "1")
                    levs.append(to_grid(np.where(frac == 1.0, val, np.nan), g))
                prof[n].append(np.stack(levs))
            sst = read_var(ds, "sst", "degree_Celsius")
            conc = read_var(ds, "ice_concentration", "1")
            tos.append(remap(sst, r))
            sic.append(100.0 * remap(conc, r))
            a, act, olat = r["area"], r["active"], r["geom"]["lat"]
            native["sst_mean"].append(float((sst * a)[act].sum() / a[act].sum()))
            native["nh_ice_area"].append(float((conc * a)[act & (olat >= 0)].sum()))
            native["sh_ice_area"].append(float((conc * a)[act & (olat < 0)].sum()))
            native["nh_ice_extent"].append(float(a[act & (olat >= 0) & (conc > 0.15)].sum()))
            native["sh_ice_extent"].append(float(a[act & (olat < 0) & (conc > 0.15)].sum()))
            native["toa_net_file"].append(float(np.sum(read_var(ds, "toa_net_downward", "W m-2") * g["weight"])))

    out = out_root / EXPERIMENT / VARIANT
    out.mkdir(parents=True)
    created = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    common = {
        "Conventions": "CF-1.8",
        "title": f"{SOURCE_ID} {EXPERIMENT} {VARIANT}, CMOR-style monthly output for ClimateBench v2 and ClimateEval",
        "source_id": SOURCE_ID,
        "source": "ReadyESM v0.1.15 physics: SpeedyWeather T31 L27 atmosphere, Oceananigans 1 degree tripolar ocean "
                  "with 60 levels, "
                  "ClimaSeaIce sea ice, Terrarium land",
        "experiment_id": EXPERIMENT,
        "experiment": "free-running run with historical greenhouse gases and ozone from an observed start in "
                      "January 1993; not a CMIP6 historical run from 1850",
        "variant_label": VARIANT,
        "institution": "University of Reading",
        "product": "model-output",
        "creation_date": created,
        "comment": "Not written by the CMOR library and not a CMIP submission. Converted offline from sealed "
                   "ReadyESM monthly files by climatebench_cmor_v1/cmorize.py.",
        "source_run": "climate_run_v015_forced_v1, segment B, Slurm job 62507908",
        "source_commit": source_attrs["source_commit"],
        "source_configuration_sha256": source_attrs["configuration_sha256"],
        "source_writer_sha256": source_attrs["writer_sha256"],
        "history": f"{created} cmorize.py from 48 sealed monthly files (see manifest.json)",
    }
    g_desc = {**g, "description": "native T31 full Gaussian grid, 48 x 96, latitude bounds from model ring weights",
              "nominal_resolution": "500 km"}
    o_desc = {"lat": r["lat"], "lon": r["lon"], "lat_bnds": r["lat_bnds"], "lon_bnds": r["lon_bnds"],
              "description": f"regular {OCEAN_RES} degree grid; native tripolar cells sampled {SAMPLES} x {SAMPLES} "
                             f"per target cell by nearest native centre",
              "nominal_resolution": "100 km"}
    written = {}
    for n, (src, _u, scale, units, std, long_name, comment) in ATMOS.items():
        path = out / filename(n, "Amon", "gn")
        write(path, n, np.stack(atm[n]), {"units": units, "standard_name": std, "long_name": long_name,
              "comment": comment, "original_name": src + (f" x {scale}" if scale != 1.0 else "")},
              g_desc, "Amon", "atmos", "gn", common, times=times, cell_methods="area: time: mean")
        written[n] = path
    for n, (prefix, levels, _u, units, std, long_name) in PROFILES.items():
        path = out / filename(n, "Amon", "gn")
        write(path, n, np.stack(prof[n]), {"units": units, "standard_name": std, "long_name": long_name,
              "comment": "time mean on pressure levels, interpolated linearly in log pressure with no extrapolation; "
                         "masked where any time step of the month was below the surface",
              "original_name": f"{prefix}_<p>hpa_contribution where {prefix}_<p>hpa_valid_fraction == 1"},
              g_desc, "Amon", "atmos", "gn", common, times=times, plev=levels, cell_methods="time: mean")
        written[n] = path
    sea_mask = r["sftof"] < WET_THRESHOLD
    for n, data, units, std, long_name, src, table, realm in (
            ("tos", np.stack(tos), "degC", "sea_surface_temperature", "Sea Surface Temperature", "sst",
             "Omon", "ocean"),
            ("siconc", np.stack(sic), "%", "sea_ice_area_fraction", "Sea-Ice Area Percentage (Ocean Grid)",
             "ice_concentration x 100", "SImon", "seaIce")):
        data = np.where(sea_mask[None], np.nan, data)
        path = out / filename(n, table, "gr")
        write(path, n, data, {"units": units, "standard_name": std, "long_name": long_name, "original_name": src,
              "comment": f"remapped from the native tripolar grid; masked where sftof < {WET_THRESHOLD:g} percent"},
              o_desc, table, realm, "gr", common, times=times, cell_methods="area: mean where sea time: mean")
        written[n] = path
    path = out / filename("sftof", "Ofx", "gr", timed=False)
    write(path, "sftof", r["sftof"], {"units": "%", "standard_name": "sea_area_fraction",
          "long_name": "Sea Area Percentage", "original_name": "wet share of target cell from ocean_active_mask",
          "comment": "share of target-cell samples whose nearest native cell centre is an active ocean cell"},
          o_desc, "Ofx", "ocean", "gr", common)
    written["sftof"] = path

    checks = verify(written, native, g, r, months)
    manifest = {
        "schema": "readyesm_cmor_style_conversion_v1",
        "created_utc": created,
        "script_sha256": sha256(Path(__file__)),
        "source_job": str(JOB),
        "inputs": [{"file": str(p.relative_to(JOB)), "sha256": dg} for _, _, p, dg in files],
        "ocean_remap": r["info"],
        "outputs": {n: {"file": str(p.relative_to(out_root)), "sha256": sha256(p)} for n, p in written.items()},
        "checks": checks,
    }
    (out_root / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    print(json.dumps(checks, indent=1))


# ---------------------------------------------------------------- checks

def read_back(path, name):
    with Dataset(path) as nc:
        v = nc[name]
        coords_ok = all("standard_name" in nc[c].ncattrs() for c in v.dimensions)
        bounds_ok = all(nc[c].getncattr("bounds") in nc.variables for c in v.dimensions
                        if "bounds" in nc[c].ncattrs())
        return np.ma.filled(v[:].astype(float), np.nan), {a: v.getncattr(a) for a in v.ncattrs()}, \
            coords_ok and bounds_ok


def verify(written, native, g, r, months):
    """Read every file back and compare it with native-grid values."""
    area_a = cell_area(g["lat_bnds"], g["lon_bnds"])
    wa = area_a / area_a.sum()
    gm = lambda x: np.nansum(x * wa, axis=(-2, -1))
    fld = {n: read_back(p, n)[0] for n, p in written.items()}
    years = np.array([y for y, _ in months])
    n_toa = gm(fld["rsdt"] - fld["rsut"] - fld["rlut"])
    pe = gm(fld["pr"] - fld["evspsbl"]) * 86400.0
    month_len = np.array([(datetime(y + m // 12, m % 12 + 1, 1) - datetime(y, m, 1)).days for y, m in months])
    # Month lengths weight the annual and full-period means, as in the Tier I readout.
    annual_n = {int(y): float(np.sum((n_toa * month_len)[years == y]) / np.sum(month_len[years == y]))
                for y in sorted(set(years))}
    # The file's own TOA net (latest cache) against rsdt - rsut - rlut on the output grid.
    toa_consistency = float(np.max(np.abs(n_toa - np.array(native["toa_net_file"]))))

    area_o = cell_area(r["lat_bnds"], r["lon_bnds"])
    sftof = fld["sftof"] / 100.0
    lat_o = r["lat"][:, None] * np.ones_like(area_o)
    sea = np.isfinite(fld["tos"][0])
    tos_mean = np.nansum(fld["tos"] * area_o * sftof, axis=(1, 2)) / np.sum(area_o * sftof * sea)
    sic = np.nan_to_num(fld["siconc"]) / 100.0
    # Area: share of masked cells' ice is lost; extent-style area counts whole unmasked cells.
    nh_area = np.sum(sic * area_o * sftof * (lat_o >= 0), axis=(1, 2))
    sh_area = np.sum(sic * area_o * sftof * (lat_o < 0), axis=(1, 2))
    nh_ext = np.sum((sic > 0.15) * area_o * (lat_o >= 0), axis=(1, 2))
    sh_ext = np.sum((sic > 0.15) * area_o * (lat_o < 0), axis=(1, 2))
    meta = {}
    for n, p in written.items():
        _, attrs, coords_ok = read_back(p, n)
        meta[n] = all(k in attrs for k in ("units", "standard_name")) and coords_ok
    sel = lambda m: [i for i, (_, mm) in enumerate(months) if mm == m]
    return {
        "atmosphere_toa_net_annual_W_m2": annual_n,
        "atmosphere_toa_net_mean_W_m2": float(np.sum(n_toa * month_len) / month_len.sum()),
        "atmosphere_toa_net_vs_file_field_max_abs_W_m2": toa_consistency,
        "atmosphere_p_minus_e_mean_mm_day": float(np.sum(pe * month_len) / month_len.sum()),
        "ocean_wet_area_native_Mkm2": float(r["area"][r["active"]].sum() / 1e12),
        "ocean_wet_area_remapped_Mkm2": float(np.sum(area_o * sftof) / 1e12),
        "ocean_wet_area_unmasked_cells_Mkm2": float(np.sum(area_o * sea) / 1e12),
        "sst_global_mean_max_abs_diff_degC": float(np.max(np.abs(tos_mean - np.array(native["sst_mean"])))),
        "sst_global_mean_native_first_last_degC": [native["sst_mean"][0], native["sst_mean"][-1]],
        "ice_area_native_vs_remapped_Mkm2": {
            "nh_march": [[native["nh_ice_area"][i] / 1e12, float(nh_area[i] / 1e12)] for i in sel(3)],
            "nh_september": [[native["nh_ice_area"][i] / 1e12, float(nh_area[i] / 1e12)] for i in sel(9)],
            "sh_february": [[native["sh_ice_area"][i] / 1e12, float(sh_area[i] / 1e12)] for i in sel(2)],
            "sh_september": [[native["sh_ice_area"][i] / 1e12, float(sh_area[i] / 1e12)] for i in sel(9)],
        },
        "ice_extent_15pct_native_vs_remapped_Mkm2": {
            "nh_march": [[native["nh_ice_extent"][i] / 1e12, float(nh_ext[i] / 1e12)] for i in sel(3)],
            "nh_september": [[native["nh_ice_extent"][i] / 1e12, float(nh_ext[i] / 1e12)] for i in sel(9)],
            "sh_february": [[native["sh_ice_extent"][i] / 1e12, float(sh_ext[i] / 1e12)] for i in sel(2)],
            "sh_september": [[native["sh_ice_extent"][i] / 1e12, float(sh_ext[i] / 1e12)] for i in sel(9)],
        },
        "metadata_units_and_standard_names_present": meta,
    }


if __name__ == "__main__":
    raise SystemExit("This file provides compatibility helpers; use cmorize_v2.py.")
