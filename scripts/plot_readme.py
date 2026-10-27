#!/usr/bin/env python3
"""Plot final-state maps and daily timelines from ReadyESM diagnostics.

Requires numpy, scipy, matplotlib and netCDF4. Maps use spherical nearest
neighbours, including dry cells, so interpolation does not cross land masks.
"""

import argparse
import hashlib
import json
from pathlib import Path
import tomllib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
import netCDF4
import numpy as np
from scipy.spatial import cKDTree

BLUE, ORANGE, GREY = "#2365a4", "#c46a14", "#b7bbc0"


def read(ds, name, units=None):
    var = ds[name]
    if units is not None and var.units != units:
        raise ValueError(f"Unexpected units for {name}: {var.units}")
    return np.ma.filled(var[:], np.nan).astype(float)


def xyz(lon, lat):
    lon, lat = np.deg2rad(lon), np.deg2rad(lat)
    return np.column_stack((np.cos(lat) * np.cos(lon),
                            np.cos(lat) * np.sin(lon), np.sin(lat)))


def nearest(lon, lat, out_lon, out_lat):
    return cKDTree(xyz(lon.ravel(), lat.ravel())).query(
        xyz(out_lon.ravel(), out_lat.ravel()))[1].reshape(out_lon.shape)


def daily_mean(time, values, days):
    """Average equally spaced endpoint samples in (day-1, day]."""
    if time.shape != values.shape or not np.all(np.diff(time) > 0):
        raise ValueError("Invalid time axis")
    dt = np.diff(time)
    if not np.allclose(dt, dt[0], rtol=0, atol=1e-8):
        raise ValueError("Daily averaging requires regular samples")
    result = []
    for day in range(1, days + 1):
        selected = (time > day - 1 + 1e-8) & (time <= day + 1e-8)
        v = values[selected]
        if len(v) != round(1 / dt[0]) or not np.all(np.isfinite(v)):
            raise ValueError(f"Incomplete radiation samples on day {day}")
        result.append(float(v.mean()))
    return np.arange(1, days + 1), np.asarray(result)


def finite_range(values):
    values = values[np.isfinite(values)]
    if not len(values):
        raise ValueError("Empty map")
    return [float(values.min()), float(values.max())]


def render(path, prefix, version, days):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.titlesize": 13, "axes.titleweight": "bold",
                         "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.facecolor": "white"})
    prefix = Path(prefix)
    outputs = {name: Path(f"{prefix}-{name}.{ext}") for name, ext in
               [("maps", "png"), ("timelines", "png"), ("figures", "json")]}
    if any(p.exists() for p in outputs.values()):
        raise FileExistsError("Choose a new output prefix")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(path) as ds:
        clocks = {n: read(ds, n) for n in
                  ("time", "atmosphere_diagnostic_time", "coupled_budget_time")}
        for t in clocks.values():
            if not (np.all(np.isfinite(t)) and np.all(np.diff(t) > 0)
                    and np.isclose(t[0], 0) and np.isclose(t[-1], days)):
                raise ValueError("Diagnostics must cover the requested run")
        for name in ("atmosphere_diagnostic_time", "coupled_budget_time"):
            if not np.array_equal(clocks[name], np.arange(days + 1)):
                raise ValueError("Expected daily state diagnostics")
        end = netCDF4.num2date(days, ds["time"].units).strftime("%-d %B %Y")
        metadata = {"version_configuration": version, "model_day": days,
                    "start_date": ds.start_date, "final_date": end,
                    "diagnostics_sha256": hashlib.file_digest(open(path, "rb"), "sha256").hexdigest(),
                    "diagnostics_filename": Path(path).name,
                    "solar_distance_correction": bool(ds.applied_solar_distance_correction),
                    "maps": {}, "timelines": {}}
        report = Path(path).with_name("report.toml")
        if report.exists():
            with report.open("rb") as f:
                source = tomllib.load(f)
            metadata["model_source_commit"] = source.get("source_commit")
            metadata["production_diagnostics_pass"] = source.get("production_diagnostics_pass")
        metadata["map_method"] = "Spherical nearest native point, including dry cells. Equirectangular display."
        metadata["timeline_method"] = "Radiation: daily mean of regularly spaced samples. Other panels: saved daily instantaneous state/flux samples. Day-zero precipitation and evaporation omitted before flux evaluation. Missing values stay gaps. Extrema cover all active ocean depths and can move between cells."

        xx, yy = np.meshgrid(np.arange(-179.75, 180, .5), np.arange(-89.75, 90, .5))
        oi = nearest(read(ds, "ocean_longitude"), read(ds, "ocean_latitude"), xx, yy)
        ai = nearest(read(ds, "atmosphere_longitude"), read(ds, "atmosphere_latitude"), xx, yy)
        mask = read(ds, "ocean_surface_active_mask")
        if not np.all(np.isin(mask, [0, 1])):
            raise ValueError("Non-binary ocean mask")
        active = mask.astype(bool)
        temp = LinearSegmentedColormap.from_list("temperature", ["#174f7a", "#76a9cf", "#edf4f7", "#fff5e6", "#f4b86a", "#9c4709"])
        ice = LinearSegmentedColormap.from_list("ice", ["#bfdbef", "#4288be", "#08306b"])
        ice.set_under("white")
        blue = plt.get_cmap("Blues").copy()
        blue.set_under("white")
        concentration = read(ds, "sea_ice_concentration", "1")
        thickness = read(ds, "sea_ice_thickness", "m")
        rain = (read(ds, "atmosphere_rainfall_flux", "kg/m^2/s")
                + read(ds, "atmosphere_snowfall_flux", "kg/m^2/s")) * 86400
        if np.any(rain < 0):
            raise ValueError("Negative precipitation")
        panels = [
            ("Sea-surface temperature", read(ds, "ocean_surface_temperature", "degree_Celsius"), "ocean", "°C", temp, Normalize(-2, 32)),
            ("Sea-surface salinity", read(ds, "ocean_surface_salinity", "1e-3"), "ocean", "Practical salinity", "viridis", Normalize(20, 42)),
            ("Sea-ice concentration", concentration * 100, "ocean", "%", ice, Normalize(1e-12, 100)),
            ("Area-equivalent ice thickness", concentration * thickness, "ocean", "m (concentration × thickness)", ice, Normalize(1e-12, 5)),
            ("Surface-air temperature", read(ds, "atmosphere_surface_temperature", "K") - 273.15, "atmosphere", "°C", temp, Normalize(-60, 40)),
            ("Instantaneous precipitation", rain, "atmosphere", "mm day⁻¹ (log scale)", blue, LogNorm(.01, 100)),
            ("Cloud fraction", read(ds, "column_cloud_fraction", "1") * 100, "atmosphere", "%", blue, Normalize(0, 100)),
            ("Surface soil saturation", read(ds, "land_surface_moisture", "1") * 100, "atmosphere", "%", blue, Normalize(0, 100)),
            ("Outgoing longwave radiation", read(ds, "outgoing_longwave", "W/m^2"), "atmosphere", "W m⁻²", "magma", Normalize(100, 350)),
        ]
        fig, axes = plt.subplots(3, 3, figsize=(17, 10), layout="constrained")
        fig.suptitle(f"ReadyESM {version} · day {days} · {end}", fontsize=21, fontweight="bold")
        for k, (ax, panel) in enumerate(zip(axes.flat, panels)):
            title, data, grid, units, cmap, norm = panel
            if grid == "ocean":
                if data.shape != active.shape or not np.all(np.isfinite(data[active])):
                    raise ValueError(f"Invalid ocean field: {title}")
                data = np.where(active, data, np.nan)
            limits = finite_range(data)
            metadata["maps"][title] = {"units": units, "range": limits, "colour_range": [norm.vmin, norm.vmax]}
            raster = data.ravel()[oi if grid == "ocean" else ai]
            cmap = plt.get_cmap(cmap).copy() if isinstance(cmap, str) else cmap.copy()
            cmap.set_bad(GREY)
            # Exact zero ice and precipitation remain white, outside LogNorm.
            if title == "Instantaneous precipitation":
                raster = np.where(raster == 0, .001, raster)
            im = ax.imshow(raster, origin="lower", extent=(-180, 180, -90, 90),
                           cmap=cmap, norm=norm, interpolation="nearest", aspect="equal")
            ax.set_title(f"{chr(65+k)}  {title}", loc="left")
            ax.set_xticks([-180, -90, 0, 90, 180], ["180°", "90°W", "0°", "90°E", "180°"])
            ax.set_yticks([-60, 0, 60], ["60°S", "0°", "60°N"])
            extend = "both" if limits[0] < norm.vmin and limits[1] > norm.vmax else "max" if limits[1] > norm.vmax else "min" if limits[0] < norm.vmin and title not in ("Sea-ice concentration", "Area-equivalent ice thickness", "Instantaneous precipitation") else "neither"
            cb = fig.colorbar(im, ax=ax, orientation="horizontal", pad=.03, fraction=.07, aspect=35, extend=extend)
            cb.set_label(units, fontsize=10)
            cb.ax.tick_params(labelsize=10)
        fig.savefig(outputs["maps"], dpi=150)
        plt.close(fig)

        ta, tc = clocks["atmosphere_diagnostic_time"], clocks["coupled_budget_time"]
        fig, axes = plt.subplots(4, 2, figsize=(14, 12), layout="constrained")
        fig.suptitle(f"ReadyESM {version} · first {days} model days", fontsize=21, fontweight="bold")

        def line(ax, key, time, values, label=None, colour=BLUE, style="-"):
            ax.plot(time, values, color=colour, linestyle=style, linewidth=1.8, label=label)
            metadata["timelines"][key] = {"day": time.tolist(), "value": [float(x) if np.isfinite(x) else None for x in values]}

        titles = [("Net top-of-atmosphere radiation", "W m⁻² (positive downward)"),
                  ("Global surface-air temperature", "°C"),
                  ("Global precipitation and evaporation", "mm day⁻¹"),
                  ("Global column water vapour", "kg m⁻²"),
                  ("Maximum ocean temperature", "°C (all active depths)"),
                  ("Minimum ocean salinity", "Practical salinity (all active depths)"),
                  ("Maximum area-equivalent ice thickness", "m (concentration × thickness)"),
                  ("Land-mean surface soil saturation", "%")]
        a = axes.flat
        tr, nr = daily_mean(clocks["time"], read(ds, "toa_net_downward", "W/m^2"), days)
        line(a[0], "toa_net_daily_mean", tr, nr)
        a[0].axhline(0, color="#606060", linewidth=.8, linestyle=":")
        line(a[1], "surface_air_temperature_c", ta, read(ds, "global_surface_air_temperature", "K") - 273.15)
        p = (read(ds, "global_rainfall_flux", "kg/m^2/s") + read(ds, "global_snowfall_flux", "kg/m^2/s")) * 86400
        e = read(ds, "global_surface_water_vapor_flux", "kg/m^2/s") * 86400
        line(a[2], "precipitation_mm_day", ta[1:], p[1:], "Precipitation")
        line(a[2], "evaporation_mm_day", ta[1:], e[1:], "Evaporation", ORANGE, "--")
        a[2].legend(frameon=False, fontsize=10, loc="lower right")
        for ax, key, time, units, scale in [
            (a[3], "global_column_water_vapor", ta, "kg/m^2", 1),
            (a[4], "ocean_temperature_maximum_series", tc, "degree_Celsius", 1),
            (a[5], "ocean_salinity_minimum_series", tc, "1e-3", 1),
            (a[6], "sea_ice_max_area_equivalent_thickness", tc, "m", 1),
            (a[7], "global_land_surface_soil_moisture", ta, "1", 100)]:
            line(ax, key, time, read(ds, key, units) * scale)
        ticks = np.unique(np.r_[np.linspace(0, days, 7).round().astype(int), days])
        for k, (ax, (title, unit)) in enumerate(zip(a, titles)):
            ax.set_title(f"{chr(65+k)}  {title}", loc="left")
            ax.set_ylabel(unit)
            ax.set_xlabel("Model day")
            ax.set_xlim(0, days)
            ax.set_xticks(ticks)
            ax.grid(axis="y", color="#e3e5e8", linewidth=.6)
            ax.set_axisbelow(True)
        fig.savefig(outputs["timelines"], dpi=150)
        plt.close(fig)
    metadata["renderer_sha256"] = hashlib.file_digest(open(__file__, "rb"), "sha256").hexdigest()
    metadata["images"] = {k: {"filename": p.name, "sha256": hashlib.file_digest(open(p, "rb"), "sha256").hexdigest()} for k, p in outputs.items() if k != "figures"}
    outputs["figures"].write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"maps": len(panels), "timeline_panels": len(titles), "days": days, "outputs": {k: str(v) for k, v in outputs.items()}}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("diagnostics", type=Path)
    parser.add_argument("output_prefix", type=Path)
    parser.add_argument("--version", required=True)
    parser.add_argument("--expected-days", type=int, required=True)
    args = parser.parse_args()
    if args.expected_days <= 0:
        parser.error("--expected-days must be positive")
    render(args.diagnostics, args.output_prefix, args.version, args.expected_days)
