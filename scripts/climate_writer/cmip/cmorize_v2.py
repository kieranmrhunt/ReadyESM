"""CMOR-style NetCDF from sealed ReadyESM monthly and annual files, version 2 (Claude, 2026-10-29).

Offline analysis code, not model code. Version 2 adds the fields that climate_writer_cmip_fields_v1 puts in the
writer. Version 1 (cmorize.py) is unchanged. Version 2 imports its grid, remap and time helpers after it checks the
file's SHA-256, so the existing variables keep the same data. A new variable is written only when every input file
has its source fields. Otherwise the manifest lists it under "not_available" with the reason.

Existing variables (as version 1): rsdt rsut rlut rsutcs rlutcs rsds rsus rlds rlus hfss evspsbl hfls pr prsn ps prw
(Amon), tos (Omon), siconc (SImon), sftof (Ofx). ta and ua come from the old 850, 500 and 200 hPa fields when the new
profile fields are missing.

New monthly atmosphere fields (Amon, native Gaussian grid gn):
  tas      cmip_tas, the lowest model level lowered dry-adiabatically to the surface. It is not a 2 m similarity
           value. The height coordinate is the nominal CMIP 2 m.
  ts       cmip_ts, radiative surface temperature
  clt      cmip_clt, maximum-random overlap of the RRTMGP cloud fraction
  clwvi    cmip_clwvi, liquid plus ice;  clivi  cmip_clivi
  tauu     cmip_tauu, tauv  cmip_tauv, the stress SpeedyWeather applies, downward positive
  ta ua va hus zg on the 19 CMIP plev19 levels. The writer samples 16 levels (1000 to 20 hPa). The 10, 5 and 1 hPa
           levels lie above the model top and are fully masked. A level is masked in a month where
           *_valid_fraction is below 1, as in version 1.
New monthly ocean fields:
  zos      cmip_zos minus its global area mean on the native grid each month (CMIP: dynamic sea level), Omon gr
  mlotst   cmip_mlotst, Omon gr
  phcint   four depth bands 0-300, 300-700, 700-2000 m and 2000 m to the bottom (cmip_ocean_heat_content_*), Omon gr.
           lev is the band midpoint with the band edges as bounds, so ClimateEval's extract_volume steps (0-300 m,
           0-2000 m, all) select whole bands.
  msftmz   cmip_msftmz, (time, basin, lev, lat), Omon gn. lev is the reference interface depth from the surface.
           Lines with no faces (cmip_overturning_face_count = 0) are masked.
  hfbasin  cmip_hfbasin, (time, basin, lat), Omon gn, masked as msftmz
  basin is a string coordinate "sector" with standard_name region, in the CMIP order atlantic_arctic_ocean,
  indian_pacific_ocean, global_ocean.
New annual ocean fields (Oyr gr, written only when the annual files have cmip_ocean_wet):
  thetao   ocean_temperature;  so  ocean_salinity;  phcint  cmip_ocean_layer_heat_content per model layer
  Each level is remapped with that level's wet mask, and masked where its wet share is below 50 percent.
  uo vo    east and north velocity from rotate_native_velocity.jl, which applies the O-015 observer to the annual
           native face velocities (--velocity-dir). Each rotated file must name the sealed annual file it came from by
           SHA-256. On the regular output grid x is east and y is north.

The reference interface depths come from the writer's cell-centre coordinate ocean_z: the top face is at 0 and each
centre is the midpoint of its two faces. In z-star runs the true depths move with the free surface. The coordinate is
the reference grid.

ClimateEval limits that apply to a run that starts in 1993: its ocean heat content diagnostic takes anomalies from
1979-1988, so it cannot run. Its ACC transport and ocean heat content items ask for monthly 3D fields, and the
writer saves 3D ocean fields once a year.

usage: python3 -I cmorize_v2.py RUN_DIR START END OUT_DIR --source-id ID --source TEXT --run-label TEXT
         --experiment-text TEXT [--experiment historical] [--variant r1i1p1f1] [--allow-incomplete]
         [--velocity-dir DIR] [--subdir run] [--allow-partial-year]
RUN_DIR holds result.sha256 and run/monthly_YYYY_MM.nc (and run/annual_YYYY.nc). START and END are YYYY-MM.
OUT_DIR receives EXPERIMENT/VARIANT/*.nc and manifest.json. It must not exist. --allow-incomplete accepts months
whose "complete" attribute is not 1 (for short test runs only). The manifest records it.
--subdir names the folder of the writer files (default run). --allow-partial-year also converts the
annual file of a year that is not complete (test runs only). The files and the manifest say so.
"""
import argparse
import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from netCDF4 import Dataset
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
V1_PATH = HERE / "cmorize.py"
V1_SHA256 = "e03eded1a0d9ad134d5d164d0e5dd0884148d1ab641882493aba6e89321e98d0"
PLEV19_HPA = (1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10, 5, 1)
WRITER_LEVELS_HPA = PLEV19_HPA[:16]
BASIN_NAMES = ("atlantic_arctic_ocean", "indian_pacific_ocean", "global_ocean")
BANDS = (("0_300m", 0.0, 300.0), ("300_700m", 300.0, 700.0), ("700_2000m", 700.0, 2000.0),
         ("below_2000m", 2000.0, None))


def load_v1():
    digest = hashlib.sha256(V1_PATH.read_bytes()).hexdigest()
    if digest != V1_SHA256:
        raise SystemExit(f"cmorize.py SHA-256 is {digest}, expected {V1_SHA256}")
    spec = importlib.util.spec_from_file_location("cmorize_v1", V1_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


V1 = load_v1()
FILL = V1.FILL

# cmor name: (source, source units, units, standard_name, long_name, comment)
NEW_SURFACE = {
    "tas": ("cmip_tas", "K", "K", "air_temperature", "Near-Surface Air Temperature",
            "lowest model level lowered to the surface with the dry adiabatic lapse rate (SpeedyWeather "
            "surface_air_temperature); not a 2 m similarity value; the height coordinate is the nominal CMIP value"),
    "ts": ("cmip_ts", "K", "K", "surface_temperature", "Surface Temperature",
           "radiative blend ((1-f) T_sea^4 + f T_land^4)^(1/4) passed to the radiation scheme"),
    "clt": ("cmip_clt", "%", "%", "cloud_area_fraction", "Total Cloud Cover Percentage",
            "maximum-random overlap of the latest RRTMGP cloud fraction"),
    "clwvi": ("cmip_clwvi", "kg m-2", "kg m-2", "atmosphere_mass_content_of_cloud_condensed_water",
              "Condensed Water Path", "grid-box mean cloud liquid plus ice, latest RRTMGP solve"),
    "clivi": ("cmip_clivi", "kg m-2", "kg m-2", "atmosphere_mass_content_of_cloud_ice", "Ice Water Path",
              "grid-box mean cloud ice, latest RRTMGP solve"),
    "tauu": ("cmip_tauu", "Pa", "Pa", "surface_downward_eastward_stress", "Surface Downward Eastward Wind Stress",
             "stress applied by SpeedyWeather SurfaceMomentumFlux"),
    "tauv": ("cmip_tauv", "Pa", "Pa", "surface_downward_northward_stress", "Surface Downward Northward Wind Stress",
             "stress applied by SpeedyWeather SurfaceMomentumFlux"),
}
# cmor name: (writer field, units, standard_name, long_name)
NEW_PROFILES = {
    "ta": ("ta", "K", "air_temperature", "Air Temperature"),
    "ua": ("ua", "m s-1", "eastward_wind", "Eastward Wind"),
    "va": ("va", "m s-1", "northward_wind", "Northward Wind"),
    "hus": ("hus", "1", "specific_humidity", "Specific Humidity"),
    "zg": ("zg", "m", "geopotential_height", "Geopotential Height"),
}
# Pairs that the old and new profile samplers both write: (new cmor name, level, old source prefix).
OVERLAP = (("ta", 850, "temperature"), ("ta", 500, "temperature"), ("ta", 200, "temperature"),
           ("ua", 850, "zonal_wind"), ("ua", 200, "zonal_wind"))


# ---------------------------------------------------------------- inputs

def months_between(start, end):
    (y0, m0), (y1, m1) = (tuple(int(x) for x in s.split("-")) for s in (start, end))
    out, y, m = [], y0, m0
    while (y, m) <= (y1, m1):
        out.append((y, m))
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    if not out:
        raise SystemExit("empty month range")
    return out


def sealed_inputs(run_dir, months, subdir="run", partial_year=False):
    seal = {}
    for line in (run_dir / "result.sha256").read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        seal[name.removeprefix("./")] = digest

    def checked(rel):
        digest = V1.sha256(run_dir / rel)
        if seal.get(rel) != digest:
            raise SystemExit(f"{rel} does not match the run's result manifest")
        return run_dir / rel, digest

    monthly = [(y, m, *checked(f"{subdir}/monthly_{y}_{m:02d}.nc")) for y, m in months]
    # Only complete years have a full annual mean. A test run may ask for a part-year file (--allow-partial-year).
    years = sorted({y for y, _ in months if partial_year or sum(1 for yy, _ in months if yy == y) == 12})
    annual = [(y, *checked(f"{subdir}/annual_{y}.nc")) for y in years if (run_dir / f"{subdir}/annual_{y}.nc").exists()]
    return monthly, annual


def has(ds, *names):
    return all(n in ds.variables for n in names)


def read_dims(ds, name, units, dims):
    """Read a variable, check its units, and return it with its dimensions in the order `dims`."""
    v = ds[name]
    if v.units != units:
        raise SystemExit(f"{name} has units {v.units!r}, expected {units!r}")
    if set(v.dimensions) != set(dims):
        raise SystemExit(f"{name} has dimensions {v.dimensions}, expected {dims}")
    a = np.asarray(np.ma.filled(v[:].astype(float), np.nan), float)
    return np.transpose(a, [v.dimensions.index(d) for d in dims])


def interface_depths(ocean_z):
    """Reference interface depths (m, positive down, surface first) from cell centres (bottom first, negative)."""
    centres = -np.asarray(ocean_z, float)[::-1]
    faces = np.zeros(centres.size + 1)
    for k, c in enumerate(centres):
        faces[k + 1] = 2.0 * c - faces[k]
    if not (np.all(np.diff(faces) > 0) and np.all((faces[:-1] < centres) & (centres < faces[1:]))):
        raise SystemExit("ocean_z centres do not give increasing interface depths")
    return centres, faces


# ---------------------------------------------------------------- per-level ocean remap

def level_remap(ds, wet3d):
    """Remap operators for each ocean level, built from the same samples as version 1's surface operator.

    wet3d is (k, j, i) with k = 0 at the bottom. Returns operators and wet shares with the surface level first.
    """
    olon = np.asarray(ds["ocean_longitude"][:], float)
    olat = np.asarray(ds["ocean_latitude"][:], float)
    flat_active = (np.asarray(ds["ocean_active_mask"][:]) == 1).reshape(-1)
    xyz = V1.unit_vectors(olon, olat).reshape(-1, 3)
    dist, _ = cKDTree(xyz[flat_active]).query(xyz[~flat_active])
    duplicate = np.zeros(flat_active.size, bool)
    duplicate[np.where(~flat_active)[0][dist < 1e-7]] = True
    candidates = np.where(flat_active | (~flat_active & ~duplicate))[0]
    tree = cKDTree(xyz[candidates])
    res, ns = V1.OCEAN_RES, V1.SAMPLES
    lat_c = np.arange(-90 + res / 2, 90, res)
    lon_c = np.arange(res / 2, 360, res)
    offs = (np.arange(ns) + 0.5) / ns - 0.5
    ny, nx = lat_c.size, lon_c.size
    LA = np.broadcast_to((lat_c[:, None] + offs[None, :] * res)[:, None, :, None], (ny, nx, ns, ns))
    LO = np.broadcast_to((lon_c[:, None] + offs[None, :] * res)[None, :, None, :], (ny, nx, ns, ns))
    w_all = np.cos(np.deg2rad(LA)).reshape(-1)
    _, k = tree.query(V1.unit_vectors(LO, LA).reshape(-1, 3))
    src = candidates[k]
    inside = LA.reshape(-1) >= olat.min() - 0.5 * res
    target = np.repeat(np.arange(ny * nx), ns * ns)
    total = np.bincount(target, weights=w_all, minlength=ny * nx)
    active_index = -np.ones(flat_active.size, int)
    active_index[flat_active] = np.arange(flat_active.sum())
    ops, shares = [], []
    for level in wet3d[::-1]:
        wet = flat_active[src] & inside & level.reshape(-1)[src]
        wet_w = np.bincount(target, weights=w_all * wet, minlength=ny * nx)
        with np.errstate(invalid="ignore", divide="ignore"):
            vals = w_all[wet] / wet_w[target[wet]]
        ops.append(sp.csr_matrix((vals, (target[wet], active_index[src[wet]])),
                                 shape=(ny * nx, int(flat_active.sum()))))
        shares.append((100.0 * wet_w / total).reshape(ny, nx))
    return ops, np.stack(shares), flat_active


def remap_levels(field, ops, shares, flat_active):
    """field is (k, j, i) bottom first; output is (lev, lat, lon) surface first, masked below 50 percent wet."""
    out = np.stack([(op @ level.reshape(-1)[flat_active]).reshape(shares.shape[1:])
                    for op, level in zip(ops, field[::-1])])
    return np.where(shares >= V1.WET_THRESHOLD, out, np.nan)


# ---------------------------------------------------------------- writing

def time_axis_years(years):
    starts = [(y, 1) for y in years]
    mid, b = V1.time_axis(starts)
    ends = V1.time_axis([(y, 12) for y in years])[1][:, 1]
    b = np.stack([b[:, 0], ends], axis=1)
    return 0.5 * (b[:, 0] + b[:, 1]), b


def write_nc(path, name, data, attrs, axes, common, table, realm, grid_label, frequency, cell_methods="",
             scalar_coords=()):
    """Write one variable. `axes` is a list of (dim, values, attributes, bounds or None). A dim named "basin"
    holds a list of strings and is written as the CMIP char coordinate "sector" with standard_name region."""
    with Dataset(path, "w", format="NETCDF4_CLASSIC") as nc:
        nc.createDimension("bnds", 2)
        coordinates = []
        for dim, values, a, bounds in axes:
            nc.createDimension(dim, len(values))
            if dim == "basin":
                width = max(len(s) for s in values)
                nc.createDimension("strlen", width)
                c = nc.createVariable("sector", "S1", ("basin", "strlen"))
                c.setncatts({"standard_name": "region", "long_name": "ocean basin", "units": ""})
                c[:] = np.array([list(s.ljust(width)) for s in values], "S1")
                coordinates.append("sector")
                continue
            c = nc.createVariable(dim, "f8", (dim,))
            c.setncatts({**a, **({"bounds": f"{dim}_bnds"} if bounds is not None else {})})
            c[:] = np.asarray(values, float)
            if bounds is not None:
                nc.createVariable(f"{dim}_bnds", "f8", (dim, "bnds"))[:] = np.asarray(bounds, float)
        for sname, value, a in scalar_coords:
            c = nc.createVariable(sname, "f8", ())
            c.setncatts(a)
            c.assignValue(value)
            coordinates.append(sname)
        v = nc.createVariable(name, "f4", tuple(d for d, *_ in axes), zlib=True, complevel=1, fill_value=FILL)
        v.setncatts({k: val for k, val in attrs.items() if val != ""})
        v.missing_value = FILL
        if coordinates:
            v.coordinates = " ".join(coordinates)
        if cell_methods:
            v.cell_methods = cell_methods
        out = np.where(np.isfinite(data), data, FILL).astype(np.float32)
        v[:] = np.ma.masked_equal(out, FILL)
        nc.setncatts({**common, "table_id": table, "realm": realm, "grid_label": grid_label, "frequency": frequency,
                      "variable_id": name})


TIME_ATTRS = {"units": V1.TIME_UNITS, "calendar": V1.CALENDAR, "standard_name": "time", "long_name": "time",
              "axis": "T"}
LAT_ATTRS = {"units": "degrees_north", "standard_name": "latitude", "long_name": "latitude", "axis": "Y"}
LON_ATTRS = {"units": "degrees_east", "standard_name": "longitude", "long_name": "longitude", "axis": "X"}
LEV_ATTRS = {"units": "m", "standard_name": "depth", "long_name": "ocean depth coordinate", "axis": "Z",
             "positive": "down"}
PLEV_ATTRS = {"units": "Pa", "standard_name": "air_pressure", "long_name": "pressure", "axis": "Z", "positive": "down"}


def horizontal(grid):
    return [("lat", grid["lat"], LAT_ATTRS, grid["lat_bnds"]), ("lon", grid["lon"], LON_ATTRS, grid["lon_bnds"])]


# ---------------------------------------------------------------- main

def main():
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("start")
    parser.add_argument("end")
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--source-id", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--run-label", required=True)
    parser.add_argument("--experiment-text", required=True, help="what the run is, for the 'experiment' attribute")
    parser.add_argument("--experiment", default="historical")
    parser.add_argument("--variant", default="r1i1p1f1")
    parser.add_argument("--allow-incomplete", action="store_true")
    parser.add_argument("--velocity-dir", type=Path, help="rotate_native_velocity.jl output for the annual files")
    parser.add_argument("--subdir", default="run", help="folder of the writer files inside RUN_DIR (default run)")
    parser.add_argument("--allow-partial-year", action="store_true",
                        help="also convert the annual file of a year that is not complete (test runs only)")
    args = parser.parse_args()
    if args.out_dir.exists():
        raise SystemExit(f"{args.out_dir} exists; refusing to overwrite")
    months = months_between(args.start, args.end)
    monthly, annual = sealed_inputs(args.run_dir, months, args.subdir, args.allow_partial_year)
    times = V1.time_axis(months)
    not_available = {}

    with Dataset(monthly[0][2]) as ds:
        g = V1.atmosphere_grid(ds)
        r = V1.ocean_remap(ds)
        centres, faces = interface_depths(ds["ocean_z"][:])
        source_attrs = {k: ds.getncattr(k) for k in ("source_commit", "configuration_sha256", "writer_sha256",
                                                     "writer_schema") if k in ds.ncattrs()}

    def present(name, *fields):
        ok = True
        for _, _, path, _ in monthly:
            with Dataset(path) as ds:
                ok = ok and has(ds, *fields)
        if not ok:
            not_available[name] = f"source fields missing in at least one monthly file: {', '.join(fields)}"
        return ok

    new_surface = {n: spec for n, spec in NEW_SURFACE.items() if present(n, spec[0])}
    new_profiles = {n: spec for n, spec in NEW_PROFILES.items() if present(
        n, *[f"cmip_{spec[0]}_{p}hpa_{s}" for p in WRITER_LEVELS_HPA for s in ("contribution", "valid_fraction")])}
    old_profiles = {n: spec for n, spec in V1.PROFILES.items() if n not in new_profiles}
    band_fields = [f"cmip_ocean_heat_content_{b}" for b, _, _ in BANDS]
    ocean2d = {n: present(n, *fields) for n, fields in (("zos", ("cmip_zos",)), ("mlotst", ("cmip_mlotst",)),
                                                        ("phcint (monthly bands)", band_fields))}
    overturning = present("msftmz", "cmip_msftmz", "cmip_hfbasin", "cmip_overturning_face_count")
    if not overturning:
        not_available["hfbasin"] = not_available["msftmz"]

    atm = {n: [] for n in V1.ATMOS} | {n: [] for n in new_surface}
    prof = {n: [] for n in old_profiles} | {n: [] for n in new_profiles}
    tos, sic, zos, mlotst, bands, psi, heat = [], [], [], [], [], [], []
    native = {"sst_mean": [], "nh_ice_area": [], "sh_ice_area": [], "nh_ice_extent": [], "sh_ice_extent": [],
              "toa_net_file": [], "zos_mean_removed_m": [], "column_heat_J": [], "overlap_max_abs": {}}
    face_count, psi_lat, complete = None, None, []
    for year, month, path, _ in monthly:
        with Dataset(path) as ds:
            complete.append(int(ds.getncattr("complete")) == 1)
            if not complete[-1] and not args.allow_incomplete:
                raise SystemExit(f"{path} is not a complete month")
            if not (np.array_equal(np.asarray(ds["atmosphere_latitude"][:], float), g["point_lat"]) and
                    np.array_equal(np.asarray(ds["atmosphere_longitude"][:], float), g["point_lon"]) and
                    np.array_equal(np.asarray(ds["ocean_latitude"][:], float), r["geom"]["lat"]) and
                    np.array_equal(np.asarray(ds["ocean_active_mask"][:]) == 1, r["active"])):
                raise SystemExit(f"{path} geometry differs from the first month")
            for n, (src, units, scale, *_rest) in V1.ATMOS.items():
                atm[n].append(V1.to_grid(V1.read_var(ds, src, units) * scale, g))
            for n, (src, units, *_rest) in new_surface.items():
                atm[n].append(V1.to_grid(V1.read_var(ds, src, units), g))
            for n, (prefix, levels, units, *_rest) in old_profiles.items():
                levs = []
                for p in levels:
                    val = V1.read_var(ds, f"{prefix}_{p}hpa_contribution", units)
                    frac = V1.read_var(ds, f"{prefix}_{p}hpa_valid_fraction", "1")
                    levs.append(V1.to_grid(np.where(frac == 1.0, val, np.nan), g))
                prof[n].append(np.stack(levs))
            for n, (field, units, *_rest) in new_profiles.items():
                levs = []
                for p in PLEV19_HPA:
                    if p not in WRITER_LEVELS_HPA:
                        levs.append(np.full((g["lat"].size, g["lon"].size), np.nan))
                        continue
                    val = V1.read_var(ds, f"cmip_{field}_{p}hpa_contribution", units)
                    frac = V1.read_var(ds, f"cmip_{field}_{p}hpa_valid_fraction", "1")
                    levs.append(V1.to_grid(np.where(frac == 1.0, val, np.nan), g))
                prof[n].append(np.stack(levs))
            for n, p, prefix in OVERLAP:
                if n in new_profiles and has(ds, f"{prefix}_{p}hpa_contribution"):
                    old = V1.read_var(ds, f"{prefix}_{p}hpa_contribution", NEW_PROFILES[n][1])
                    new = V1.read_var(ds, f"cmip_{n}_{p}hpa_contribution", NEW_PROFILES[n][1])
                    key = f"{n}_{p}hpa"
                    native["overlap_max_abs"][key] = max(native["overlap_max_abs"].get(key, 0.0),
                                                         float(np.max(np.abs(new - old))))
            sst = V1.read_var(ds, "sst", "degree_Celsius")
            conc = V1.read_var(ds, "ice_concentration", "1")
            tos.append(V1.remap(sst, r))
            sic.append(100.0 * V1.remap(conc, r))
            a, act, olat = r["area"], r["active"], r["geom"]["lat"]
            native["sst_mean"].append(float((sst * a)[act].sum() / a[act].sum()))
            native["nh_ice_area"].append(float((conc * a)[act & (olat >= 0)].sum()))
            native["sh_ice_area"].append(float((conc * a)[act & (olat < 0)].sum()))
            native["nh_ice_extent"].append(float(a[act & (olat >= 0) & (conc > 0.15)].sum()))
            native["sh_ice_extent"].append(float(a[act & (olat < 0) & (conc > 0.15)].sum()))
            native["toa_net_file"].append(float(np.sum(V1.read_var(ds, "toa_net_downward", "W m-2") * g["weight"])))
            if ocean2d["zos"]:
                eta = V1.read_var(ds, "cmip_zos", "m")
                mean = float((eta * a)[act].sum() / a[act].sum())
                native["zos_mean_removed_m"].append(mean)
                zos.append(V1.remap(eta - mean, r))
            if ocean2d["mlotst"]:
                mlotst.append(V1.remap(V1.read_var(ds, "cmip_mlotst", "m"), r))
            if ocean2d["phcint (monthly bands)"]:
                per_band = [V1.read_var(ds, f, "J m-2") for f in band_fields]
                native["column_heat_J"].append(float((sum(per_band) * a)[act].sum()))
                bands.append(np.stack([V1.remap(b, r) for b in per_band]))
            if overturning:
                dims = ("cmip_basin", "cmip_ocean_interface", "cmip_overturning_latitude")
                p = read_dims(ds, "cmip_msftmz", "kg s-1", dims)
                h = read_dims(ds, "cmip_hfbasin", "W", (dims[0], dims[2]))
                f = read_dims(ds, "cmip_overturning_face_count", "1", (dims[0], dims[2]))
                label = ds["cmip_msftmz"].long_name
                if "basins " + ", ".join(BASIN_NAMES) not in label or p.shape[0] != len(BASIN_NAMES):
                    raise SystemExit(f"cmip_msftmz basins are not {BASIN_NAMES}: {label!r}")
                if p.shape[1] != faces.size:
                    raise SystemExit("cmip_msftmz interfaces do not match ocean_z")
                lat = np.asarray(ds["cmip_overturning_latitude"][:], float)
                if face_count is None:
                    face_count, psi_lat = f, lat
                elif not (np.array_equal(f, face_count) and np.array_equal(lat, psi_lat)):
                    raise SystemExit(f"{path}: overturning lines differ from the first month")
                mask = (f > 0)
                psi.append(np.where(mask[:, None, :], p[:, ::-1, :], np.nan))  # surface interface first
                heat.append(np.where(mask, h, np.nan))

    out_root = args.out_dir
    out = out_root / args.experiment / args.variant
    out.mkdir(parents=True)
    created = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    common = {
        "Conventions": "CF-1.8",
        "title": f"{args.source_id} {args.experiment} {args.variant}, CMOR-style output for ClimateBench v2 and "
                 "ClimateEval",
        "source_id": args.source_id, "source": args.source, "experiment_id": args.experiment,
        "experiment": args.experiment_text,
        "variant_label": args.variant, "institution": "University of Reading", "product": "model-output",
        "creation_date": created,
        "comment": "Not written by the CMOR library and not a CMIP submission. Converted offline from sealed "
                   "ReadyESM writer files by climatebench_cmor_v1/cmorize_v2.py.",
        "source_run": args.run_label,
        **{("source_" + k if not k.startswith("source") else k): v for k, v in source_attrs.items()},
        "history": f"{created} cmorize_v2.py from {len(monthly)} monthly and {len(annual)} annual sealed files "
                   "(see manifest.json)",
    }
    if not all(complete):
        common["comment"] += " TEST OUTPUT: includes incomplete months."
    if args.allow_partial_year:
        common["comment"] += " TEST OUTPUT: annual files may cover part of a year."
    span = f"{months[0][0]}{months[0][1]:02d}-{months[-1][0]}{months[-1][1]:02d}"
    name_of = lambda n, table, grid, tail=span: (
        f"{n}_{table}_{args.source_id}_{args.experiment}_{args.variant}_{grid}" + (f"_{tail}" if tail else "") + ".nc")
    tax = [("time", times[0], TIME_ATTRS, times[1])]
    atm_grid = {"description": "native T31 full Gaussian grid, 48 x 96, latitude bounds from model ring weights",
                "nominal_resolution": "500 km"}
    ocn_grid = {"description": f"regular {V1.OCEAN_RES} degree grid; native tripolar cells sampled {V1.SAMPLES} x "
                               f"{V1.SAMPLES} per target cell by nearest native centre",
                "nominal_resolution": "100 km"}
    line_grid = {"description": "native tripolar grid, staircase lines of faces at the target latitudes",
                 "nominal_resolution": "100 km"}
    amon = {**common, "grid": atm_grid["description"], "nominal_resolution": atm_grid["nominal_resolution"]}
    ogr = {**common, "grid": ocn_grid["description"], "nominal_resolution": ocn_grid["nominal_resolution"]}
    ogn = {**common, "grid": line_grid["description"], "nominal_resolution": line_grid["nominal_resolution"]}
    written = {}

    def emit(key, n, table, grid_label, frequency, data, attrs, axes, glob, cell_methods="", scalar_coords=(),
             tail=span):
        path = out / name_of(n, table, grid_label, tail)
        write_nc(path, n, data, attrs, axes, glob, table, {"A": "atmos", "O": "ocean", "S": "seaIce"}[table[0]],
                 grid_label, frequency, cell_methods, scalar_coords)
        written[key] = path

    for n, (src, _u, scale, units, std, long_name, comment) in V1.ATMOS.items():
        emit(n, n, "Amon", "gn", "mon", np.stack(atm[n]),
             {"units": units, "standard_name": std, "long_name": long_name, "comment": comment,
              "original_name": src + (f" x {scale}" if scale != 1.0 else "")},
             tax + horizontal(g), amon, "area: time: mean")
    height = (("height", 2.0, {"units": "m", "standard_name": "height", "long_name": "height", "axis": "Z",
                               "positive": "up"}),)
    for n, (src, _u, units, std, long_name, comment) in new_surface.items():
        emit(n, n, "Amon", "gn", "mon", np.stack(atm[n]),
             {"units": units, "standard_name": std, "long_name": long_name, "comment": comment,
              "original_name": src}, tax + horizontal(g), amon, "area: time: mean",
             scalar_coords=height if n == "tas" else ())
    for n, (prefix, levels, _u, units, std, long_name) in old_profiles.items():
        emit(n, n, "Amon", "gn", "mon", np.stack(prof[n]),
             {"units": units, "standard_name": std, "long_name": long_name,
              "comment": "time mean on pressure levels, interpolated linearly in log pressure with no "
                         "extrapolation; masked where any time step of the month was below the surface",
              "original_name": f"{prefix}_<p>hpa_contribution where {prefix}_<p>hpa_valid_fraction == 1"},
             tax + [("plev", np.asarray(levels, float) * 100.0, PLEV_ATTRS, None)] + horizontal(g), amon,
             "time: mean")
    for n, (field, units, std, long_name) in new_profiles.items():
        emit(n, n, "Amon", "gn", "mon", np.stack(prof[n]),
             {"units": units, "standard_name": std, "long_name": long_name,
              "comment": "time mean on plev19, interpolated linearly in log pressure with no extrapolation; "
                         "masked where any time step of the month was outside the model's full-level range; "
                         "10, 5 and 1 hPa lie above the model top and are fully masked",
              "original_name": f"cmip_{field}_<p>hpa_contribution where cmip_{field}_<p>hpa_valid_fraction == 1"},
             tax + [("plev", np.asarray(PLEV19_HPA, float) * 100.0, PLEV_ATTRS, None)] + horizontal(g), amon,
             "time: mean")
    sea_mask = r["sftof"] < V1.WET_THRESHOLD
    o_axes = tax + horizontal({"lat": r["lat"], "lon": r["lon"], "lat_bnds": r["lat_bnds"],
                               "lon_bnds": r["lon_bnds"]})
    remap_comment = f"remapped from the native tripolar grid; masked where sftof < {V1.WET_THRESHOLD:g} percent"
    for n, data, units, std, long_name, src, table in (
            ("tos", tos, "degC", "sea_surface_temperature", "Sea Surface Temperature", "sst", "Omon"),
            ("siconc", sic, "%", "sea_ice_area_fraction", "Sea-Ice Area Percentage (Ocean Grid)",
             "ice_concentration x 100", "SImon"),
            ("zos", zos, "m", "sea_surface_height_above_geoid", "Sea Surface Height Above Geoid",
             "cmip_zos minus its native global area mean each month", "Omon"),
            ("mlotst", mlotst, "m", "ocean_mixed_layer_thickness_defined_by_sigma_t",
             "Ocean Mixed Layer Thickness Defined by Sigma T", "cmip_mlotst", "Omon")):
        if not data:
            continue
        emit(n, n, table, "gr", "mon", np.where(sea_mask[None], np.nan, np.stack(data)),
             {"units": units, "standard_name": std, "long_name": long_name, "original_name": src,
              "comment": remap_comment}, o_axes, ogr, "area: mean where sea time: mean")
    if bands:
        bottom = faces[-1]
        edges = np.array([[top, bottom if low is None else low] for _, top, low in BANDS])
        emit("phcint", "phcint", "Omon", "gr", "mon", np.where(sea_mask[None, None], np.nan, np.stack(bands)),
             {"units": "J m-2",
              "standard_name": "integral_wrt_depth_of_sea_water_potential_temperature_expressed_as_heat_content",
              "long_name": "Integrated Ocean Heat Content from Potential Temperature",
              "original_name": ", ".join(band_fields),
              "comment": "rho_0 c_p theta integrated over four depth bands below the free surface, theta in degC "
                         "(reference 0 degC), accumulated each step from the current cell thickness; lev is the "
                         "band midpoint and lev_bnds the band edges; " + remap_comment},
             tax[:1] + [("lev", edges.mean(axis=1), LEV_ATTRS, edges)] + o_axes[1:], ogr,
             "area: mean where sea time: mean")
    if psi:
        basin_axis = ("basin", list(BASIN_NAMES), {}, None)
        lat_axis = ("lat", psi_lat, LAT_ATTRS, None)
        emit("msftmz", "msftmz", "Omon", "gn", "mon", np.stack(psi),
             {"units": "kg s-1", "standard_name": "ocean_meridional_overturning_mass_streamfunction",
              "long_name": "Ocean Meridional Overturning Mass Streamfunction", "original_name": "cmip_msftmz",
              "comment": "staircase lines of native faces at each latitude, summed from the surface; lev is the "
                         "reference interface depth; masked where the basin has no line at that latitude"},
             tax + [basin_axis, ("lev", faces, LEV_ATTRS, None), lat_axis], ogn, "time: mean")
        emit("hfbasin", "hfbasin", "Omon", "gn", "mon", np.stack(heat),
             {"units": "W", "standard_name": "northward_ocean_heat_transport",
              "long_name": "Northward Ocean Heat Transport", "original_name": "cmip_hfbasin",
              "comment": "centred face temperature times the face volume flux on the overturning lines, reference "
                         "0 degC; not the model's own advective flux"},
             tax + [basin_axis, lat_axis], ogn, "time: mean")
    emit("sftof", "sftof", "Ofx", "gr", "fx", r["sftof"],
         {"units": "%", "standard_name": "sea_area_fraction", "long_name": "Sea Area Percentage",
          "original_name": "wet share of target cell from ocean_active_mask",
          "comment": "share of target-cell samples whose nearest native cell centre is an active ocean cell"},
         o_axes[1:], ogr, tail="")

    annual_checks = write_annual(annual, months, faces, centres, emit, ogr, not_available, native, args.velocity_dir)
    checks = verify(written, native, g, r, months, faces)
    checks["annual"] = annual_checks
    manifest = {
        "schema": "readyesm_cmor_style_conversion_v2",
        "created_utc": created,
        "script_sha256": V1.sha256(Path(__file__)),
        "version1_sha256": V1_SHA256,
        "run_dir": str(args.run_dir),
        "run_label": args.run_label,
        "complete_months": all(complete),
        "writer_subdir": args.subdir,
        "partial_year_annual_allowed": args.allow_partial_year,
        "inputs": [{"file": str(p.relative_to(args.run_dir)), "sha256": d} for _, _, p, d in monthly] +
                  [{"file": str(p.relative_to(args.run_dir)), "sha256": d} for _, p, d in annual],
        "ocean_remap": r["info"],
        "ocean_reference_interfaces_m": faces.tolist(),
        "outputs": {n: {"file": str(p.relative_to(out_root)), "sha256": V1.sha256(p)} for n, p in written.items()},
        "not_available": not_available,
        "checks": checks,
    }
    (out_root / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    print(json.dumps({"written": sorted(written), "not_available": not_available, "checks": checks}, indent=1))


def write_annual(annual, months, faces, centres, emit, ogr, not_available, native, velocity_dir=None):
    if not annual:
        not_available["thetao, so, phcint (annual)"] = "no complete year with an annual file"
        return {}
    with Dataset(annual[0][1]) as ds:
        if not has(ds, "cmip_ocean_wet", "ocean_temperature", "ocean_salinity", "cmip_ocean_layer_heat_content"):
            not_available["thetao, so, phcint (annual)"] = "annual files lack cmip_ocean_wet or the 3D fields"
            return {}
        dims = ("ocean_k", "ocean_j", "ocean_i")
        wet3d = read_dims(ds, "cmip_ocean_wet", "1", dims)
        if not np.all((wet3d == 0) | (wet3d == 1)):
            raise SystemExit("cmip_ocean_wet is not a 0/1 mask")
        wet3d = wet3d == 1
        ops, shares, flat_active = level_remap(ds, wet3d)
        reference = V1.ocean_remap(ds)
        # With every level wet, the per-level operator must equal version 1's surface operator.
        all_wet = level_remap(ds, np.ones_like(wet3d[:1], bool))[0][0]
    if (all_wet - reference["op"]).count_nonzero() != 0:
        raise SystemExit("per-level remap with all levels wet differs from version 1's surface operator")
    years = [y for y, _, _ in annual]
    fields = {"thetao": [], "so": [], "phcint": []}
    if velocity_dir is None:
        not_available["uo, vo (annual)"] = "no --velocity-dir from rotate_native_velocity.jl"
    else:
        fields |= {"uo": [], "vo": []}
    column_heat = []
    for year, path, digest in annual:
        with Dataset(path) as ds:
            if not np.array_equal(read_dims(ds, "cmip_ocean_wet", "1", dims) == 1, wet3d):
                raise SystemExit(f"{path}: wet mask differs from the first year")
            theta = read_dims(ds, "ocean_temperature", "degree_Celsius", dims)
            salt = read_dims(ds, "ocean_salinity", "psu", dims)
            layer = read_dims(ds, "cmip_ocean_layer_heat_content", "J m-2", dims)
            a = np.asarray(ds["ocean_cell_area"][:], float)
            act = np.asarray(ds["ocean_active_mask"][:]) == 1
            column_heat.append(float((layer.sum(axis=0) * a)[act].sum()))
            for n, x in (("thetao", theta), ("so", salt), ("phcint", layer)):
                fields[n].append(remap_levels(np.where(wet3d, x, np.nan), ops, shares, flat_active))
        if velocity_dir is not None:
            with Dataset(velocity_dir / f"geographic_velocity_{year}.nc") as vel:
                if vel.getncattr("source_sha256") != digest:
                    raise SystemExit(f"rotated velocity for {year} was not made from the sealed annual file")
                for n, src in (("uo", "eastward_velocity"), ("vo", "northward_velocity")):
                    x = read_dims(vel, src, "m s-1", dims)
                    fields[n].append(remap_levels(np.where(wet3d, x, np.nan), ops, shares, flat_active))
    times = time_axis_years(years)
    tail = f"{years[0]}-{years[-1]}"
    from_lat = lambda: {"lat": np.arange(-90 + V1.OCEAN_RES / 2, 90, V1.OCEAN_RES),
                        "lon": np.arange(V1.OCEAN_RES / 2, 360, V1.OCEAN_RES)}
    h = from_lat()
    hb = {"lat_bnds": np.stack([h["lat"] - V1.OCEAN_RES / 2, h["lat"] + V1.OCEAN_RES / 2], axis=1),
          "lon_bnds": np.stack([h["lon"] - V1.OCEAN_RES / 2, h["lon"] + V1.OCEAN_RES / 2], axis=1)}
    axes = [("time", times[0], TIME_ATTRS, times[1]),
            ("lev", centres, LEV_ATTRS, np.stack([faces[:-1], faces[1:]], axis=1))] + horizontal({**h, **hb})
    comment = (f"remapped level by level from the native tripolar grid with that level's wet mask; masked where the "
               f"level's wet share is below {V1.WET_THRESHOLD:g} percent; lev is the reference cell centre")
    for n, units, std, long_name, src in (
            ("thetao", "degC", "sea_water_potential_temperature", "Sea Water Potential Temperature",
             "ocean_temperature"),
            ("so", "0.001", "sea_water_salinity", "Sea Water Salinity", "ocean_salinity (psu)"),
            ("phcint", "J m-2", "integral_wrt_depth_of_sea_water_potential_temperature_expressed_as_heat_content",
             "Integrated Ocean Heat Content from Potential Temperature", "cmip_ocean_layer_heat_content"),
            ("uo", "m s-1", "sea_water_x_velocity", "Sea Water X Velocity",
             "cmip_native_u and cmip_native_v, centred and rotated east by the O-015 observer"),
            ("vo", "m s-1", "sea_water_y_velocity", "Sea Water Y Velocity",
             "cmip_native_u and cmip_native_v, centred and rotated north by the O-015 observer")):
        extra = ("; the rotation angle comes from Oceananigans rotation_angle, whose direction error for a solid-body "
                 "test flow is below 0.25 percent south of 85N, about 1.4 percent at 85-88N and about 5 percent at "
                 "88-89N") if n in ("uo", "vo") else ""
        if n not in fields:
            continue
        emit(f"{n}_annual", n, "Oyr", "gr", "yr", np.stack(fields[n]),
             {"units": units, "standard_name": std, "long_name": long_name, "original_name": src,
              "comment": comment + extra}, axes, ogr, "area: mean where sea time: mean", tail=tail)
    # The month-length mean of the monthly band totals equals the annual layer total when both see the same steps.
    checks = {"column_heat_annual_layers_J": column_heat}
    if native["column_heat_J"]:
        lengths = np.array([(datetime(y + m // 12, m % 12 + 1, 1) - datetime(y, m, 1)).days for y, m in months])
        yrs = np.array([y for y, _ in months])
        heat = np.array(native["column_heat_J"])
        from_bands = [float(np.sum((heat * lengths)[yrs == y]) / np.sum(lengths[yrs == y])) for y in years]
        checks["column_heat_from_monthly_bands_J"] = from_bands
        checks["column_heat_relative_difference"] = [abs(b - c) / abs(c) for b, c in zip(from_bands, column_heat)]
    checks["surface_level_wet_share_vs_sftof_max_abs_percent"] = float(np.max(np.abs(shares[0] - reference["sftof"])))
    checks["level_wet_share_never_increases_with_depth"] = bool(np.all(np.diff(shares, axis=0) <= 1e-9))
    return checks


# ---------------------------------------------------------------- checks

SV_KG_S = 1.025e9   # 1 Sv of sea water at 1025 kg m-3, for the readable checks only


def metadata_ok(path, name):
    """Units and standard_name on the variable, standard_name on every coordinate it uses, and bounds present."""
    with Dataset(path) as nc:
        v = nc[name]
        names = [d for d in v.dimensions if d in nc.variables] + getattr(v, "coordinates", "").split()
        dims_covered = all(d in nc.variables or (d == "basin" and "sector" in names) for d in v.dimensions)
        return (all(k in v.ncattrs() for k in ("units", "standard_name")) and dims_covered and
                all("standard_name" in nc[c].ncattrs() for c in names) and
                all(nc[c].getncattr("bounds") in nc.variables for c in names if "bounds" in nc[c].ncattrs()))

def verify(written, native, g, r, months, faces):
    """Read the files back and compare with native-grid values. Version 1 checks first, then the new fields."""
    v1_written = {n: p for n, p in written.items() if n in V1.ATMOS or n in ("ta", "ua", "tos", "siconc", "sftof")}
    v1_native = {k: native[k] for k in ("sst_mean", "nh_ice_area", "sh_ice_area", "nh_ice_extent", "sh_ice_extent",
                                        "toa_net_file")}
    checks = V1.verify(v1_written, v1_native, g, r, months)
    fld = lambda n: V1.read_back(written[n], n)[0]
    checks["new_metadata_units_and_standard_names_present"] = {
        n: metadata_ok(p, n.removesuffix("_annual")) for n, p in written.items() if n not in v1_written}
    area_a = V1.cell_area(g["lat_bnds"], g["lon_bnds"])
    wa = area_a / area_a.sum()
    gm = lambda x: float(np.mean(np.nansum(x * wa, axis=(-2, -1))))
    if "tas" in written and "ts" in written:
        checks["tas_minus_ts_global_mean_K"] = gm(fld("tas") - fld("ts"))
    if "clt" in written:
        c = fld("clt")
        checks["clt_range_percent"] = [float(np.nanmin(c)), float(np.nanmax(c))]
    if "clivi" in written and "clwvi" in written:
        checks["clivi_minus_clwvi_max_kg_m2"] = float(np.nanmax(fld("clivi") - fld("clwvi")))
    if "tauu" in written:
        checks["tauu_global_mean_Pa"] = gm(fld("tauu"))
    checks["new_vs_old_profile_sample_max_abs"] = native["overlap_max_abs"]
    if native["zos_mean_removed_m"]:
        checks["zos_native_global_mean_removed_m"] = [min(native["zos_mean_removed_m"]),
                                                      max(native["zos_mean_removed_m"])]
    if "msftmz" in written:
        with Dataset(written["msftmz"]) as nc:
            psi = np.ma.filled(nc["msftmz"][:].astype(float), np.nan)
            lat = nc["lat"][:]
            sector = ["".join(c.decode() for c in row).strip() for row in nc["sector"][:]]
        with Dataset(written["hfbasin"]) as nc:
            heat = np.ma.filled(nc["hfbasin"][:].astype(float), np.nan)
        j = int(np.argmin(np.abs(lat - 26.5)))
        deep = faces >= 500.0
        atl = psi[:, 0, :, j]
        checks["basin_names_read_back"] = sector
        amoc = np.nanmax(atl[:, deep], axis=1) / SV_KG_S
        checks["amoc_26N_max_below_500m_Sv_monthly_range"] = [float(np.nanmin(amoc)), float(np.nanmax(amoc))]
        checks["amoc_latitude_used"] = float(lat[j])
        checks["msftmz_surface_interface_max_abs_kg_s"] = float(np.nanmax(np.abs(psi[:, :, 0, :])))
        checks["msftmz_global_bottom_max_abs_Sv"] = float(np.nanmax(np.abs(psi[:, 2, -1, :])) / SV_KG_S)
        checks["hfbasin_26N_mean_PW"] = {s: float(np.nanmean(heat[:, b, j]) / 1e15) for b, s in enumerate(sector)}
    return checks


if __name__ == "__main__":
    main()
