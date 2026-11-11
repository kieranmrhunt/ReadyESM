#!/usr/bin/env python3
"""Plot the sealed, common-mask six-month Part C comparison for v0.1.20."""
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

RESULT_SEAL = "802ec3617e06ee7b0e6b60aaee26c7aea33c0bcfbc76c49f0fcffa1bcb414c53"
FBA_RESULT_SEAL = "5be8e9032614f7a2785bac0fd04785ad93750402b2ee2e90a757f22dc40cd118"
RUNS = (
    ("v0.1.19 published 66742450", "v0.1.19", "#64748b"),
    ("v0.1.20 reverted 66841309", "v0.1.20", "#087e8b"),
    ("v0.1.20 all-on candidate 66632474", "all-on candidate", "#bd6a36"),
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
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ratio_matrix(run):
    matrix = np.array([
        [run["per_field"][field][region]["readyesm_rmse"] /
         run["per_field"][field][region]["cmip6_median_of_model_medians"]
         for region, _ in REGIONS]
        for field, _ in FIELDS
    ], dtype=float)
    if not np.all(np.isfinite(matrix)) or np.any(matrix <= 0):
        raise ValueError("non-finite or non-positive field error ratio")
    return matrix


def render(score_path, manifest_path, prefix):
    score_path, manifest_path, prefix = map(Path, (score_path, manifest_path, prefix))
    if sha256(manifest_path) != RESULT_SEAL:
        raise ValueError("joint score manifest seal changed")
    subprocess.run(["sha256sum", "--strict", "--quiet", "-c", manifest_path.name],
                   cwd=manifest_path.parent, check=True)
    score = json.loads(score_path.read_text())
    if score["schema"] != "readyesm_partc_amended_full_v3" or score["cohort_size"] != 55:
        raise ValueError("unexpected Part C score schema or cohort")
    if set(score["runs"]) != {name for name, _, _ in RUNS}:
        raise ValueError("unexpected scored runs")
    if score["runs"]["v0.1.20 reverted 66841309"]["result_sha256"] != FBA_RESULT_SEAL:
        raise ValueError("reverted run result seal changed")
    commits = {
        "v0.1.19 published 66742450": "816fc17d1fa1955002393abdbe80da57ceeae139",
        "v0.1.20 reverted 66841309": "fba77cd3ffc42f0f9227f9cd2ca24840bb8ad396",
        "v0.1.20 all-on candidate 66632474": "f08496a9a68eebf82f76242e321ee9a69386c4b9",
    }
    for name, commit in commits.items():
        if score["runs"][name]["source_commit"] != commit:
            raise ValueError("scored source commit changed: " + name)
    old = score["runs"][RUNS[0][0]]
    new = score["runs"][RUNS[1][0]]
    if not new["onr8"] < old["onr8"]:
        raise ValueError("declared release comparison no longer improves")
    matrices = (ratio_matrix(old), ratio_matrix(new))
    values = (score["runs"][RUNS[0][0]]["onr8"],
              score["runs"][RUNS[1][0]]["onr8"],
              score["runs"][RUNS[2][0]]["onr8"])
    typical = new["cohort_median_of_model_medians_onr8"]

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig = plt.figure(figsize=(16, 11.7), facecolor="white")
    gs = fig.add_gridspec(2, 3, left=.13, right=.94, bottom=.16, top=.84,
                          width_ratios=[1, 1, .045], height_ratios=[.92, 2.6],
                          wspace=.18, hspace=.34)
    fig.text(.5, .965, "ReadyESM v0.1.20 · six-month climate comparison",
             ha="center", va="top", fontsize=21, weight="bold", color="#24323c")
    fig.text(.5, .927,
             "January–June 1993 · same observed starting state · 55 CMIP6 models",
             ha="center", va="top", fontsize=11, color="#4c5964")
    ax = fig.add_subplot(gs[0, :2])
    labels = ["Typical CMIP6 model", "v0.1.19", "v0.1.20 released", "earlier v0.1.20 candidate"]
    bar_values = [typical, values[0], values[1], values[2]]
    colors = ["#b4bdc4", RUNS[0][2], RUNS[1][2], RUNS[2][2]]
    y = np.arange(len(labels))
    ax.barh(y, bar_values, color=colors, height=.62)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.9)
    ax.set_xlabel("Eight-field overall error ratio · lower is better")
    ax.set_title("Overall error against ERA5", loc="left", weight="bold", pad=12)
    ax.grid(axis="x", color="#e5e9ed", lw=.8)
    ax.set_axisbelow(True)
    for i, value in enumerate(bar_values):
        ax.text(value + .025, i, f"{value:.2f}", va="center", fontsize=11, weight="bold")
    improvement = 100 * (1 - values[1] / values[0])
    ax.text(.98, 1.04, f"v0.1.20 is {improvement:.1f}% lower than v0.1.19 on this test",
            transform=ax.transAxes, ha="right", color=RUNS[1][2], weight="bold")

    cmap = plt.get_cmap("RdBu_r")
    norm = TwoSlopeNorm(vmin=.5, vcenter=1, vmax=2.2)
    for col, (title, matrix) in enumerate(zip(("v0.1.19", "v0.1.20 released"), matrices)):
        hm = fig.add_subplot(gs[1, col])
        image = hm.imshow(matrix, cmap=cmap, norm=norm, aspect="auto")
        hm.set_title(title + " · regional error ratios", loc="left", weight="bold", pad=12)
        hm.set_xticks(np.arange(len(REGIONS)), [label for _, label in REGIONS], rotation=45, ha="right")
        hm.set_yticks(np.arange(len(FIELDS)),
                      [label for _, label in FIELDS] if col == 0 else [""] * len(FIELDS))
        hm.tick_params(length=0)
        hm.set_xticks(np.arange(-.5, len(REGIONS), 1), minor=True)
        hm.set_yticks(np.arange(-.5, len(FIELDS), 1), minor=True)
        hm.grid(which="minor", color="white", linewidth=1.4)
        for row in range(matrix.shape[0]):
            for column in range(matrix.shape[1]):
                value = matrix[row, column]
                red, green, blue, _ = cmap(norm(value))
                brightness = .2126 * red + .7152 * green + .0722 * blue
                hm.text(column, row, f"{value:.2f}", ha="center", va="center",
                        color="white" if brightness < .47 else "#24323c", fontsize=8)
    cax = fig.add_subplot(gs[1, 2])
    cb = fig.colorbar(image, cax=cax)
    cb.set_label("Regional RMSE / CMIP6 model median")
    cb.set_ticks([.5, .75, 1, 1.5, 2.2])
    fig.text(.13, .085,
             "Blue below 1: smaller error than the CMIP6 median. Red above 1: larger error. Values above 2.2 use the top colour.",
             fontsize=9, color="#4c5964")
    fig.text(.13, .060,
             "850 hPa is scored only above ground. Missing locations are dropped for every model, which can make all scores look better.",
             fontsize=9, color="#4c5964")
    fig.text(.13, .035,
             "One run from an observed starting state does not establish long-term climate skill. Earlier v0.1.19 release failures remain.",
             fontsize=9, color="#4c5964")

    png, pdf, metadata_path = (Path(str(prefix) + suffix) for suffix in (".png", ".pdf", ".json"))
    prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=170, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    plt.close(fig)
    metadata = {
        "schema": "readyesm_v020_reverted_partc_scorecard_v1",
        "score_schema": score["schema"], "cohort_model_count": score["cohort_size"],
        "cohort_model_year_count": score["runs"][RUNS[1][0]]["pooled_model_years"],
        "source_run_jobs": [66742450, 66841309, 66632474],
        "source_commits": commits,
        "score_json_sha256": sha256(score_path),
        "score_result_manifest_sha256": sha256(manifest_path),
        "common_masks_sha256": score["common_masks_sha256"],
        "fba_run_result_manifest_sha256": FBA_RESULT_SEAL,
        "onr8": {name: score["runs"][name]["onr8"] for name, _, _ in RUNS},
        "onr9": {name: score["runs"][name]["onr9"] for name, _, _ in RUNS},
        "typical_cmip6_onr8": typical,
        "fraction_model_years_worse_onr8": {
            name: score["runs"][name]["fraction_model_years_worse_onr8"] for name, _, _ in RUNS
        },
        "field_order": [field for field, _ in FIELDS],
        "region_order": [region for region, _ in REGIONS],
        "regional_error_ratio": {
            RUNS[i][0]: {
                field: {region: float(matrices[i][row, column])
                        for column, (region, _) in enumerate(REGIONS)}
                for row, (field, _) in enumerate(FIELDS)
            } for i in range(2)
        },
        "area_fraction_of_era5_mask": score["area_fraction_of_era5_mask"],
        "comparison_scope": "One continuous January-June 1993 run per version from the same observed starting state. The comparison includes 55 CMIP6 models and 1650 model-years. It uses locations valid for every model.",
        "850_hpa_rule": "retained only where defined, with no below-ground extrapolation",
        "mask_caveat": "Locations with any missing data are dropped for every model, especially 850 hPa near high ground. This can make every model look slightly better.",
        "png_sha256": sha256(png), "pdf_sha256": sha256(pdf),
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print("V020_REVERTED_PARTC_SCORECARD_PASS", png)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("score_json")
    parser.add_argument("score_result_manifest")
    parser.add_argument("output_prefix")
    args = parser.parse_args()
    render(args.score_json, args.score_result_manifest, args.output_prefix)
