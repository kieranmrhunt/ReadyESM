#!/usr/bin/env python3
"""Add a separate v0.1.20 six-month panel above the historical version card."""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


OLD_IMAGE_SHA256 = "b31036fe64be758d0a401a91a76dc1e8400211ca2cf041800805ba745b8205db"
OLD_RECORD_SHA256 = "f7a3c4f2f1a6e4385315c0462ef402ba57d9c9f1ac4fab7a8207ae2eb35f98d1"
SCORE_MANIFEST_SHA256 = "802ec3617e06ee7b0e6b60aaee26c7aea33c0bcfbc76c49f0fcffa1bcb414c53"
RUNS = (
    ("Typical CMIP6 model", "typical_cmip6_onr8", "#adb8c2"),
    ("v0.1.19, six months", "v0.1.19 published 66742450", "#64748b"),
    ("v0.1.20, six months", "v0.1.20 reverted 66841309", "#087e8b"),
    ("Earlier candidate, not released", "v0.1.20 all-on candidate 66632474", "#bd6a36"),
)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def render(old_image, old_record, score_record, prefix):
    for path, expected in ((old_image, OLD_IMAGE_SHA256), (old_record, OLD_RECORD_SHA256)):
        if sha256(path) != expected:
            raise ValueError(f"historical source changed: {path}")
    old = json.loads(old_record.read_text())
    if old["figure_sha256"] != OLD_IMAGE_SHA256 or old["publication_status"] != "qualified_release_with_failed_score":
        raise ValueError("historical scorecard status changed")
    score = json.loads(score_record.read_text())
    if (score["schema"] != "readyesm_v020_reverted_partc_scorecard_v1"
            or score["score_result_manifest_sha256"] != SCORE_MANIFEST_SHA256
            or score["cohort_model_count"] != 55
            or score["cohort_model_year_count"] != 1650):
        raise ValueError("v0.1.20 scorecard source changed")
    if score["source_commits"]["v0.1.20 reverted 66841309"] != "fba77cd3ffc42f0f9227f9cd2ca24840bb8ad396":
        raise ValueError("released run source changed")
    values = [score["typical_cmip6_onr8"]] + [score["onr8"][name] for _, name, _ in RUNS[1:]]
    if not (values[0] < values[2] < values[1] < values[3]):
        raise ValueError("unexpected six-month score ordering")

    historical = plt.imread(old_image)
    height, width = historical.shape[:2]
    if (width, height) != (3695, 2785):
        raise ValueError("historical figure dimensions changed")
    total_height = 3810
    old_fraction = height / total_height
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig = plt.figure(figsize=(width / 200, total_height / 200), dpi=200, facecolor="white")
    ax_old = fig.add_axes([0, 0, 1, old_fraction])
    ax_old.imshow(historical, interpolation="nearest", aspect="auto")
    ax_old.set_axis_off()

    fig.text(.5, .985, "ReadyESM v0.1.20 · six-month climate score",
             ha="center", va="top", fontsize=23, weight="bold", color="#24323c")
    fig.text(.5, .956,
             "January–June 1993 from the same observed starting state · 55 CMIP6 models",
             ha="center", va="top", fontsize=12, color="#4c5964")
    ax = fig.add_axes([.19, .785, .69, .15])
    y = list(range(len(RUNS)))
    ax.barh(y, values, color=[row[2] for row in RUNS], height=.60)
    ax.set_yticks(y, [row[0] for row in RUNS], fontsize=11)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.86)
    ax.set_xlabel("Eight-field overall error ratio · lower is better", fontsize=11)
    ax.grid(axis="x", color="#e5e9ed", linewidth=.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    for index, value in enumerate(values):
        ax.text(value + .025, index, f"{value:.2f}", va="center", fontsize=12, weight="bold")
    fig.text(.5, .755,
             "Six-month scores use locations valid for every model, and 850 hPa only above ground. This can make all scores look better.",
             ha="center", va="top", fontsize=10, color="#4c5964")
    fig.text(.5, .738,
             "The January/July version history below uses a different 31-day method. The v0.1.19 release notes list the checks it did not pass.",
             ha="center", va="top", fontsize=10, color="#4c5964")

    png = Path(str(prefix) + ".png")
    pdf = Path(str(prefix) + ".pdf")
    meta_path = Path(str(prefix) + ".json")
    if any(path.exists() for path in (png, pdf, meta_path)):
        raise FileExistsError("choose an unused output prefix")
    png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=200, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    plt.close(fig)
    metadata = {
        "schema": "readyesm_v020_version_history_separate_methods_v1",
        "historical_v019_image_sha256": sha256(old_image),
        "historical_v019_record_sha256": sha256(old_record),
        "v020_scorecard_record_sha256": sha256(score_record),
        "v020_score_result_manifest_sha256": SCORE_MANIFEST_SHA256,
        "v020_run_source_commit": score["source_commits"]["v0.1.20 reverted 66841309"],
        "cmip6_model_count": 55,
        "cmip6_model_year_count": 1650,
        "six_month_onr8": {label: value for (label, _, _), value in zip(RUNS, values)},
        "method_note": "The new January-June 1993 panel uses one continuous run per ReadyESM version and locations valid for every model. The older January/July history remains an unchanged, separate 31-day comparison. Neither establishes long-term climate skill.",
        "png_sha256": sha256(png),
        "pdf_sha256": sha256(pdf),
    }
    meta_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print("V020_VERSION_HISTORY_FIGURE_PASS", png)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("historical_png", type=Path)
    parser.add_argument("historical_record", type=Path)
    parser.add_argument("v020_score_record", type=Path)
    parser.add_argument("output_prefix", type=Path)
    args = parser.parse_args()
    render(args.historical_png, args.historical_record, args.v020_score_record, args.output_prefix)
