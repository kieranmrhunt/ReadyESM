#!/usr/bin/env python3
"""Render the v0.1.17 README comparison from its pinned numeric record.

Usage: python docs/plot-v0.1.17-comparison.py \
    docs/readiesm-v0.1.17-comparison.json \
    docs/readiesm-v0.1.17-comparison.png
"""

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


BLUE = "#2365a4"
ORANGE = "#b66016"
INK = "#26313b"


def load_record(path):
    record = json.loads(Path(path).read_text())
    if record.get("schema") != "readiesm_v017_comparison_figure_v1":
        raise ValueError("Unexpected figure record")
    members = record["members"]
    expected = {(season, member) for season in ("JAN", "JUL") for member in (1, 2, 3)}
    if len(members) != 6 or {(r["season"], r["member"]) for r in members} != expected:
        raise ValueError("Expected exactly three paired members in each season")
    for row in members:
        if not all(math.isfinite(row[key]) for key in ("control", "candidate")):
            raise ValueError("Non-finite score")
    for season in ("JAN", "JUL"):
        subset = [r for r in members if r["season"] == season]
        for arm in ("control", "candidate"):
            mean = sum(r[arm] for r in subset) / 3
            if not math.isclose(mean, record["means"][season][arm], abs_tol=1e-12):
                raise ValueError(f"{season} {arm} mean does not match the frozen gate")
    eddy = record["eddy_screen"]
    if eddy["members_per_season"] != 1 or not math.isclose(eddy["target_fall"], 0.10):
        raise ValueError("Unexpected eddy-screen scope or target")
    if any(not 0 <= eddy["fall"][season] < eddy["target_fall"] for season in ("JAN", "JUL")):
        raise ValueError("Eddy-screen finding differs from the saved decision")
    return record


def render(record, output):
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.facecolor": "white",
    })
    fig, (score_ax, eddy_ax) = plt.subplots(
        1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [1.55, 1]},
        layout="constrained",
    )
    fig.suptitle("ReadyESM v0.1.17 compared with v0.1.16", fontsize=17, fontweight="bold")

    positions = {("JAN", 1): 6, ("JAN", 2): 5, ("JAN", 3): 4,
                 ("JUL", 1): 2.5, ("JUL", 2): 1.5, ("JUL", 3): 0.5}
    score_ax.axvline(0, color=INK, linewidth=1.1)
    score_ax.axhline(3.25, color="#d7dce0", linewidth=0.9)
    for row in record["members"]:
        y = positions[(row["season"], row["member"])]
        delta = row["candidate"] - row["control"]
        colour = BLUE if delta < 0 else ORANGE
        score_ax.barh(y, delta, height=0.55, color=colour, zorder=2)
        label = f"{delta:+.4f}" if abs(delta) < 0.001 else f"{delta:+.3f}"
        score_ax.text(delta + (0.0013 if delta >= 0 else -0.0013), y,
                      label, ha="left" if delta >= 0 else "right",
                      va="center", fontsize=9, color=INK)
    score_ax.set_yticks([6, 5, 4, 2.5, 1.5, 0.5],
                        ["Jan · member 1", "Jan · member 2", "Jan · member 3",
                         "Jul · member 1", "Jul · member 2", "Jul · member 3"])
    score_ax.set_xlim(-0.063, 0.047)
    score_ax.set_xticks([-0.06, -0.04, -0.02, 0, 0.02, 0.04])
    score_ax.set_xlabel("Change in normalised error (candidate − control)")
    score_ax.set_title("A  Original 7-field score · three pairs per season", loc="left")
    score_ax.grid(axis="x", color="#e3e7ea", linewidth=0.7)
    score_ax.set_axisbelow(True)
    score_ax.legend(handles=[
        Line2D([0], [0], color=BLUE, linewidth=7, label="Improved"),
        Line2D([0], [0], color=ORANGE, linewidth=7, label="Worsened"),
    ], loc="upper right", ncol=2, frameon=False, fontsize=9)

    eddy = record["eddy_screen"]
    for y, season in [(1, "JAN"), (0, "JUL")]:
        fall = eddy["fall"][season]
        eddy_ax.barh(y, fall, height=0.5, color=BLUE, zorder=2)
        eddy_ax.text(fall + 0.003, y, f"{fall:.3f}", va="center", color=INK)
    eddy_ax.axvline(eddy["target_fall"], color=ORANGE, linewidth=1.7, linestyle="--")
    eddy_ax.text(eddy["target_fall"] + 0.002, 1.52, "Target 0.10", color=ORANGE, fontsize=9)
    eddy_ax.set_xlim(0, 0.125)
    eddy_ax.set_ylim(-0.65, 1.8)
    eddy_ax.set_yticks([1, 0], ["January", "July"])
    eddy_ax.set_xlabel("Fall in mean absolute log EKE250 error")
    eddy_ax.set_title("B  Eddy screen · one pair per season", loc="left")
    eddy_ax.grid(axis="x", color="#e3e7ea", linewidth=0.7)
    eddy_ax.set_axisbelow(True)
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("record", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    render(load_record(args.record), args.output)
