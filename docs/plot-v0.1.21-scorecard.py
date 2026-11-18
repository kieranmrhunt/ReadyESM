#!/usr/bin/env python3
"""Render the v0.1.21 release card from the sealed four-run Part C score."""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np


SCORE_SEAL = "42ad50448892d02bde93352d5dce80550a97701495a34475b35df57742734bb4"
SCORE_JSON = "69449672e24f2ae4fb4f23ce75d5baf0be6b147a2c28b08111f2ec5b32ba370a"
RUNS = (
    ("v0.1.19 published 66742450", "v0.1.19", "#64748b", "816fc17d1fa1955002393abdbe80da57ceeae139"),
    ("v0.1.20 reverted 66841309", "v0.1.20", "#087e8b", "fba77cd3ffc42f0f9227f9cd2ca24840bb8ad396"),
    ("v0.1.21 53279", "v0.1.21", "#7951a8", "396474675f78358401fcea445c4b805e1432a4d0"),
)
FIELDS = (
    ("t850", "Temperature 850 hPa"), ("t500", "Temperature 500 hPa"),
    ("t200", "Temperature 200 hPa"), ("u850", "Zonal wind 850 hPa"),
    ("u200", "Zonal wind 200 hPa"), ("cwv", "Column water vapour"),
    ("precip", "Precipitation"), ("olr", "Outgoing longwave"),
    ("osr", "Reflected shortwave"),
)
REGIONS = (
    ("nh_mid", "30–60°N"), ("subtropics", "20–30° N/S"),
    ("tropics", "20°S–20°N"), ("sh_mid", "30–60°S"),
    ("nh_polar", "60–90°N"), ("sh_polar", "60–90°S"),
)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def matrix(run):
    values = np.array([
        [run["per_field"][field][region]["readyesm_rmse"] /
         run["per_field"][field][region]["cmip6_median_of_model_medians"]
         for region, _ in REGIONS]
        for field, _ in FIELDS
    ])
    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("invalid regional error ratio")
    return values


def render(score_path, seal_path, png_path, pdf_path, record_path):
    if sha256(seal_path) != SCORE_SEAL or sha256(score_path) != SCORE_JSON:
        raise ValueError("frozen score receipt changed")
    subprocess.run(["sha256sum", "--strict", "--quiet", "-c", seal_path.name],
                   cwd=seal_path.parent, check=True)
    score = json.loads(score_path.read_text())
    if score["schema"] != "readyesm_partc_amended_full_v3" or score["cohort_size"] != 55:
        raise ValueError("unexpected score protocol or cohort")
    if set(score["runs"]) != {*(run[0] for run in RUNS), "v0.1.21 no cap 2495744"}:
        raise ValueError("expected the frozen four-run release comparison")
    for name, _, _, commit in RUNS:
        if score["runs"][name]["source_commit"] != commit:
            raise ValueError("release run source changed: " + name)
    if score["runs"][RUNS[2][0]]["result_sha256"] != "f0ba2b6a6257bfff746488371df2933162a0351ae43cb666e2c32d56fdcccfe0":
        raise ValueError("v0.1.21 run receipt changed")

    values = [score["runs"][name]["onr8"] for name, *_ in RUNS]
    typical = score["runs"][RUNS[2][0]]["cohort_median_of_model_medians_onr8"]
    maps = [matrix(score["runs"][RUNS[i][0]]) for i in (1, 2)]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig = plt.figure(figsize=(16, 11.8), facecolor="white")
    gs = fig.add_gridspec(2, 3, left=.13, right=.94, bottom=.17, top=.83,
                          width_ratios=[1, 1, .045], height_ratios=[.9, 2.6],
                          wspace=.18, hspace=.36)
    fig.text(.5, .965, "ReadyESM v0.1.21 · climate comparison", ha="center",
             va="top", fontsize=21, weight="bold", color="#24323c")
    fig.text(.5, .925, "One January–June 1993 run per version · March–June scored · 55 CMIP6 models",
             ha="center", va="top", fontsize=11, color="#4c5964")

    ax = fig.add_subplot(gs[0, :2])
    labels = ["Typical CMIP6 model", "v0.1.19", "v0.1.20", "v0.1.21"]
    bars = [typical, *values]
    ax.barh(range(4), bars, color=["#b4bdc4", *(run[2] for run in RUNS)], height=.62)
    ax.set_yticks(range(4), labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.55)
    ax.set_xlabel("Eight-field overall error ratio · lower is better")
    ax.set_title("Overall error against ERA5", loc="left", weight="bold", pad=12)
    ax.grid(axis="x", color="#e5e9ed", lw=.8)
    ax.set_axisbelow(True)
    for i, value in enumerate(bars):
        ax.text(value + .025, i, f"{value:.2f}", va="center", fontsize=11, weight="bold")
    ax.text(.98, 1.04, f"v0.1.21 − v0.1.20: {values[2]-values[1]:+.3f} (within single-start spread)",
            transform=ax.transAxes, ha="right", color=RUNS[2][2], weight="bold")

    cmap = plt.get_cmap("RdBu_r")
    norm = TwoSlopeNorm(vmin=.5, vcenter=1, vmax=2.2)
    for col, (title, values_map) in enumerate(zip(("v0.1.20", "v0.1.21"), maps)):
        hm = fig.add_subplot(gs[1, col])
        image = hm.imshow(values_map, cmap=cmap, norm=norm, aspect="auto")
        hm.set_title(title + " · regional error ratios", loc="left", weight="bold", pad=12)
        hm.set_xticks(np.arange(len(REGIONS)), [label for _, label in REGIONS], rotation=45, ha="right")
        hm.set_yticks(np.arange(len(FIELDS)),
                      [label for _, label in FIELDS] if col == 0 else [""] * len(FIELDS))
        hm.tick_params(length=0)
        hm.set_xticks(np.arange(-.5, len(REGIONS), 1), minor=True)
        hm.set_yticks(np.arange(-.5, len(FIELDS), 1), minor=True)
        hm.grid(which="minor", color="white", linewidth=1.4)
        for row in range(values_map.shape[0]):
            for column in range(values_map.shape[1]):
                v = values_map[row, column]
                rgb = cmap(norm(v))
                brightness = .2126 * rgb[0] + .7152 * rgb[1] + .0722 * rgb[2]
                hm.text(column, row, f"{v:.2f}", ha="center", va="center",
                        color="white" if brightness < .47 else "#24323c", fontsize=8)
    cb = fig.colorbar(image, cax=fig.add_subplot(gs[1, 2]))
    cb.set_label("Regional RMSE / CMIP6 model median")
    cb.set_ticks([.5, .75, 1, 1.5, 2.2])
    fig.text(.13, .091, "Blue below 1: smaller error than the CMIP6 median. Red above 1: larger error. Values above 2.2 use the top colour.",
             fontsize=9, color="#4c5964")
    fig.text(.13, .064, "850 hPa uses ERA5 ps > 87500 Pa and shared finite data; some model below-surface values may remain. Masking can change ranks.",
             fontsize=9, color="#4c5964")
    fig.text(.13, .037, "One observed start and four scored months do not establish long-term climate skill or equilibrium.",
             fontsize=9, color="#4c5964")

    for path in (png_path, pdf_path, record_path):
        path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=170, facecolor="white")
    fig.savefig(pdf_path, facecolor="white")
    plt.close(fig)
    record = {
        "schema": "readyesm_v021_partc_scorecard_v1",
        "score_schema": score["schema"],
        "score_json_sha256": SCORE_JSON,
        "score_result_manifest_sha256": SCORE_SEAL,
        "common_masks_sha256": score["common_masks_sha256"],
        "cohort_model_count": score["cohort_size"],
        "cohort_model_year_count": score["runs"][RUNS[2][0]]["pooled_model_years"],
        "source_commits": {name: commit for name, _, _, commit in RUNS},
        "source_run_jobs": [66742450, 66841309, 53279],
        "onr8": {name: value for (name, *_), value in zip(RUNS, values)},
        "onr9": {name: score["runs"][name]["onr9"] for name, *_ in RUNS},
        "typical_cmip6_onr8": typical,
        "regional_error_ratio": {
            RUNS[i+1][0]: {
                field: {region: float(maps[i][row, column])
                        for column, (region, _) in enumerate(REGIONS)}
                for row, (field, _) in enumerate(FIELDS)
            } for i in range(2)
        },
        "area_fraction_of_era5_mask": score["area_fraction_of_era5_mask"],
        "method_note": "One continuous January–June 1993 run per version from the same observed start. March–June monthly means are scored equally on a common-valid mask across four ReadyESM arms and 55 complete-case CMIP6 models. The 850-hPa rule uses ERA5 ps > 87500 Pa and shared finite data; some model values below their own surface may remain. The mask may change scores or rankings.",
        "png_sha256": sha256(png_path),
        "pdf_sha256": sha256(pdf_path),
    }
    record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print("V021_PARTC_SCORECARD_PASS", png_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("score_json", "score_result_manifest", "png", "pdf", "record"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    render(args.score_json, args.score_result_manifest, args.png, args.pdf, args.record)
