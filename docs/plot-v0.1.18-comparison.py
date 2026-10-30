#!/usr/bin/env python3
"""Build and render the v0.1.18 README comparison from a sealed neutral readout.

Barry wrote the first version (29 Oct 2026). Claude made this plain-label successor (30 Oct 2026) after Kieran asked
that no internal jargon appear on public figures. The scored numbers, sources and checks are unchanged. This version
also records the same-host run cost from the paired reports. Version 3 of the record (30 Oct) removes the maintainer
statement from the figure and its record, as Kieran asked, and instead requires the record to state the rule result.

Two steps, as for v0.1.17, so the figure is drawn from a pinned numeric record:

  python3 -I plot-v0.1.18-comparison.py build READOUT_DIR RECORD.json \
      --array CANDIDATE_ARRAY --score-job SCORE_JOB --score-seal RESULT_SHA256 --guard-evidence GUARD_RESULTS_DIR
  python3 -I plot-v0.1.18-comparison.py render RECORD.json OUTPUT.png [--stamp-record]

READOUT_DIR is the `readout` folder of the sealed neutral score job
(v018_neutral_readout_v1). The frozen v3 and v4 scorers are only read.
Their verdicts are copied as scored. The neutral release decision is a
separate field and is drawn in its own box.

Controls are six matched months with the v0.1.17 settings on an earlier commit, a2f2ac4.
They are not the published v0.1.17 tag.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

SCHEMA = "readiesm_v018_comparison_figure_v3"
SEASONS = ("JAN", "JUL")
MEMBERS = (1, 2, 3)
CONTROL_LABEL = "matched v0.1.17 run (v0.1.17 settings on earlier commit a2f2ac4)"
CANDIDATE_LABEL = "v0.1.18"
CONTROL_SOURCE = "a2f2ac4d4d6cf8fadd17d43808acf54d87c7f91d"
CANDIDATE_SOURCE = "48c232cc6d313370c9d5c8fd2a52745af4e6c2ac"  # b3: guard repair, Field diagnostic fix (tree eb1ea01c088cbe39e5e22f29e2d377d3e27f1fe4)
GUARDS = ("olr_abs_bias", "albedo_abs_bias", "net_abs_bias")
NEUTRAL_DECISIONS = {"ACCEPT_NEUTRAL_OR_IMPROVED", "REJECT"}
RANK_CAVEAT = (
    "These runs start from observations. In a separate test of six scored fields, "
    "a climatology and persistence beat every CMIP6 model. The rank shows "
    "agreement with the observed month, not free-running climate skill."
)
RANK_BASELINE_SOURCE = "development_wave_20261025/initialisation_advantage_v1/baselines.json"
RANK_BASELINE_SHA256 = "417fe8872b678b397b1ae253d91960d313fb0295c813fdaceaff8f6836163aa0"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def close(a, b, message):
    require(math.isclose(a, b, rel_tol=0, abs_tol=1e-12), message)


def load_guard_evidence(directory):
    """Verify the real-GPU guard wrapper result for the candidate source and return a small record of it."""
    import subprocess
    import tomllib
    directory = Path(directory)
    require(directory.is_dir() and (directory / "result.sha256").is_file(), "guard evidence directory or manifest missing")
    subprocess.run(["sha256sum", "--strict", "--quiet", "-c", "result.sha256"], cwd=directory, check=True)
    info = dict(line.split("=", 1) for line in (directory / "run_info").read_text().splitlines() if "=" in line)
    require(info.get("commit") == CANDIDATE_SOURCE, "guard evidence is for another source commit")
    require(info.get("exit_code") == "0", "guard wrapper run did not exit 0")
    result = tomllib.loads((directory / "result.toml").read_text())
    require(result.get("source_commit") == CANDIDATE_SOURCE, "guard result names another source commit")
    require(result.get("verdict") == "GPU_GUARD_WRAPPER_PASS", "guard wrapper verdict is not PASS")
    require(result.get("failures") == [] and result.get("setup_errors") == [], "guard wrapper has failures or setup errors")
    # cases are the entries that state an expected outcome. Other dict entries (for example recorded field types) are information.
    checked = {name: c for name, c in result["cases"].items() if isinstance(c, dict) and "expected" in c}
    require(checked and all(c.get("ok") is True and c.get("got") == c.get("expected") for c in checked.values()),
            "a guard wrapper case did not match its expected outcome")
    timing = result["cases"]
    return {
        "source_commit": CANDIDATE_SOURCE, "device": result.get("device"), "verdict": result["verdict"],
        "cases_checked": len(checked), "run_start_utc": info.get("start"), "run_end_utc": info.get("end"),
        "result_toml_sha256": sha(directory / "result.toml"), "result_manifest_sha256": sha(directory / "result.sha256"),
        "single_wrapper_sample_ms_per_full_pass": float(timing["timing/milliseconds_per_full_pass"]),
        "reference_step_ms": float(timing["timing/control_step_milliseconds"]),
        "arrays_per_pass": int(timing["timing/arrays_per_pass"]),
        "timing_note": "one timing sample of the state check on a shared GPU; it is not the run-cost estimate (see run_cost)",
    }


def member_values(readout, season, member, version, field):
    path = Path(readout) / f"{season}_m{member}" / version / "scorecard.json"
    runs = load(path)["aggregate"]["runs"]
    return {"control": float(runs[f"c{member}"][field]), "candidate": float(runs[f"x{member}"][field])}, path


def run_cost(verification):
    """Per-model-day cost of each pair from the sealed report.toml day clocks (day 1 to day 31). Same host only."""
    import tomllib
    pairs = []
    for pair in verification["pairs"]:
        rate, host = {}, {}
        for arm in ("control", "candidate"):
            job = Path(pair[arm])
            elapsed = tomllib.loads((job / "output" / "report.toml").read_text())["elapsed_day_seconds"]
            first, last = min(int(k) for k in elapsed), max(int(k) for k in elapsed)
            require((first, last) == (1, 31), "expected day clocks for days 1 to 31")
            rate[arm] = (elapsed[str(last)] - elapsed[str(first)]) / (last - first)
            host[arm] = (job / "host").read_text().strip()
        pairs.append({"season": pair["season"], "member": pair["member"], "control_host": host["control"],
                      "candidate_host": host["candidate"], "same_host": host["control"] == host["candidate"],
                      "control_seconds_per_day": rate["control"], "candidate_seconds_per_day": rate["candidate"],
                      "percent_longer": 100.0 * (rate["candidate"] / rate["control"] - 1.0)})
    same = [p["percent_longer"] for p in pairs if p["same_host"]]
    require(len(same) >= 2, "too few same-host pairs for a cost range")
    return {"method": "report.toml elapsed_day_seconds, day 1 to day 31, same GPU host only", "pairs": pairs,
            "same_host_count": len(same), "same_host_min_percent": min(same), "same_host_max_percent": max(same),
            "same_host_mean_percent": sum(same) / len(same)}


def build(readout, array, score_job, score_seal, guard_evidence):
    readout = Path(readout)
    gate3 = load(readout / "version_gate_v3.json")
    gate4 = load(readout / "version_gate.json")
    physical = load(readout / "breakage.json")
    neutral = load(readout / "neutral_adjudication.json")
    verification = load(readout / "verification.json")
    require(gate3.get("schema") == "readiesm_version_gate_v3", "wrong v3 gate schema")
    require(gate4.get("schema") == "readiesm_version_gate_v4_draft", "wrong v4 gate schema")
    require(neutral.get("schema") == "v018_neutral_skill_adjudication_v1", "wrong neutral schema")
    require(neutral["skill_adjudication"] in NEUTRAL_DECISIONS, "invalid neutral decision")
    require(verification["sources"] == {"control": CONTROL_SOURCE, "candidate": CANDIDATE_SOURCE},
            "unexpected control or candidate source")
    require(str(verification["candidate_array"]) == str(array), "readout belongs to another candidate array")
    require(len(verification["pairs"]) == 6, "six pairs required")
    members = []
    for season in SEASONS:
        for member in MEMBERS:
            onr7, p3 = member_values(readout, season, member, "v3", "onr_all_fields")
            onr9, p4 = member_values(readout, season, member, "v4", "onr9")
            pos3, _ = member_values(readout, season, member, "v3", "cmip6_models_worse")
            pos4, _ = member_values(readout, season, member, "v4", "cmip6_models_worse")
            members.append({"season": season, "member": member, "onr7": onr7, "onr9": onr9,
                            "cmip6_position_v3": pos3, "cmip6_position_v4": pos4,
                            "v3_source_file": f"{season}_m{member}/v3/scorecard.json", "v3_source_sha256": sha(p3),
                            "v4_source_file": f"{season}_m{member}/v4/scorecard.json", "v4_source_sha256": sha(p4)})
    seasons = {}
    for season in SEASONS:
        r3, r4 = gate3["seasons"][season], gate4["seasons"][season]
        n = neutral["seasons"][season]
        seasons[season] = {
            "onr7": {k: float(r3["onr7"][k]) for k in ("control_mean", "candidate_mean", "delta_mean", "noise")},
            "onr9": {k: float(r4["skill"]["onr9"][k]) for k in ("control_mean", "candidate_mean", "delta_mean", "noise")},
            "cmip6_position_v3": {k: float(r3["cmip6_position"][k]) for k in ("control_mean", "candidate_mean", "delta_mean", "noise")},
            "cmip6_position_v4": {k: float(r4["skill"]["cmip6_position"][k]) for k in ("control_mean", "candidate_mean", "delta_mean", "noise")},
            "original_v3_verdict": r3.get("verdict"),
            "original_v4_verdict": r4.get("verdict"),
            "noise_outcome": n["season_outcome"],
        }
    guards = {name: {k: gate4["two_season_global_guards"][name][k]
                     for k in ("control_mean", "candidate_mean", "allowance", "passing")} for name in GUARDS}
    guards["two_season_net"] = {"passing": gate4["two_season_net_guard"]["passing"]}
    guards["physical_breakage"] = {"passing": physical["passing"]}
    record = {
        "schema": SCHEMA,
        "scope": "Three matched 31-day runs from observed starts in each of January and July 1993. The original "
                 "score uses seven fields and the wider score nine. Both are error relative to the CMIP6 median.",
        "cmip6_rank_caveat": RANK_CAVEAT,
        "cmip6_rank_caveat_source": RANK_BASELINE_SOURCE,
        "cmip6_rank_caveat_source_sha256": RANK_BASELINE_SHA256,
        "control": CONTROL_LABEL, "candidate": CANDIDATE_LABEL,
        "control_source": CONTROL_SOURCE, "candidate_source": CANDIDATE_SOURCE,
        "candidate_array": str(array), "score_job": str(score_job), "score_result_manifest_sha256": score_seal,
        "members": members, "seasons": seasons, "guards": guards,
        "release_decision": {
            "original_v3_metric_verdict": gate3["metric_verdict"],
            "original_v4_metric_verdict": gate4["metric_verdict"],
            "neutral_rule_decision": neutral["skill_adjudication"],
            "no_clear_worse": neutral["no_clear_worse"],
            "all_required_guards_pass": neutral["all_required_guards_pass"],
            "neutral_adjudication_sha256": sha(readout / "neutral_adjudication.json"),
        },
        "runtime_guard_evidence": load_guard_evidence(guard_evidence),
        "run_cost": run_cost(verification),
        "rule_result_shown": "pass" if neutral["skill_adjudication"] == "ACCEPT_NEUTRAL_OR_IMPROVED" else "reject",
        "release_status": "pending: restart and release documentation gates",
        "renderer": "docs/plot-v0.1.18-comparison.py",
        "renderer_sha256": sha(__file__),
        "image_sha256": None,
    }
    check_record(record)
    return record


def check_record(record):
    require(record.get("schema") == SCHEMA, "unexpected figure record")
    require(record.get("cmip6_rank_caveat") == RANK_CAVEAT, "missing or changed CMIP6 rank caveat")
    require(record.get("cmip6_rank_caveat_source") == RANK_BASELINE_SOURCE and
            record.get("cmip6_rank_caveat_source_sha256") == RANK_BASELINE_SHA256,
            "missing or changed six-field baseline provenance")
    members = record["members"]
    require(len(members) == 6 and {(r["season"], r["member"]) for r in members} ==
            {(s, m) for s in SEASONS for m in MEMBERS}, "expected exactly three paired members in each season")
    for row in members:
        for metric in ("onr7", "onr9"):
            require(all(math.isfinite(row[metric][arm]) for arm in ("control", "candidate")), "non-finite score")
    for season in SEASONS:
        subset = [r for r in members if r["season"] == season]
        for metric in ("onr7", "onr9", "cmip6_position_v3", "cmip6_position_v4"):
            gate = record["seasons"][season][metric]
            for arm in ("control", "candidate"):
                close(sum(r[metric][arm] for r in subset) / 3, gate[arm + "_mean"],
                      f"{season} {metric} {arm} mean does not match the frozen gate")
            close(gate["candidate_mean"] - gate["control_mean"], gate["delta_mean"],
                  f"{season} {metric} delta does not match the frozen gate")
            require(math.isfinite(gate["noise"]) and gate["noise"] >= 0, f"{season} {metric} invalid noise bound")
    evidence = record.get("runtime_guard_evidence") or {}
    require(evidence.get("verdict") == "GPU_GUARD_WRAPPER_PASS" and evidence.get("source_commit") == record["candidate_source"],
            "the record needs passing guard wrapper evidence for the candidate source")
    require(record.get("release_status", "").startswith("pending"), "release status must stay pending in this record")
    decision = record["release_decision"]
    require(decision["neutral_rule_decision"] in NEUTRAL_DECISIONS, "invalid neutral decision")
    require(set(record["guards"]) == set(GUARDS) | {"two_season_net", "physical_breakage"}, "incomplete guard set")
    accept = decision["no_clear_worse"] is True and decision["all_required_guards_pass"] is True
    require((decision["neutral_rule_decision"] == "ACCEPT_NEUTRAL_OR_IMPROVED") == accept,
            "neutral decision contradicts its own flags")
    expected_shown = "pass" if decision["neutral_rule_decision"] == "ACCEPT_NEUTRAL_OR_IMPROVED" else "reject"
    require(record.get("rule_result_shown") == expected_shown,
            "the figure must state the release rule result exactly as scored")
    require("maintainer_acceptance" not in record, "the public figure record carries no maintainer statement")
    cost = record.get("run_cost") or {}
    require(cost.get("same_host_count", 0) >= 2 and
            cost["same_host_min_percent"] <= cost["same_host_mean_percent"] <= cost["same_host_max_percent"],
            "missing or inconsistent run cost")
    guard_ok = all(record["guards"][g]["passing"] is True for g in record["guards"])
    require(guard_ok == decision["all_required_guards_pass"], "guard rows disagree with the neutral decision")
    return record


def render(record, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyBboxPatch

    blue, orange, ink, grey = "#2365a4", "#b66016", "#26313b", "#6b7782"
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 12,
                         "axes.titleweight": "bold", "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(14.2, 8.2), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], width_ratios=[1, 1])
    fig.suptitle("ReadyESM v0.1.18 compared with matched v0.1.17 runs",
                 fontsize=17, fontweight="bold")
    positions = {("JAN", 1): 6, ("JAN", 2): 5, ("JAN", 3): 4, ("JUL", 1): 2.5, ("JUL", 2): 1.5, ("JUL", 3): 0.5}
    spans = {"JAN": (3.6, 6.4), "JUL": (0.1, 2.9)}

    def bars(ax, metric, title, xlabel):
        values = [r[metric]["candidate"] - r[metric]["control"] for r in record["members"]]
        limit = max(max(abs(v) for v in values), 0.005)  # keep the axis open when every change and bound is zero
        for season, (low, high) in spans.items():
            noise = record["seasons"][season][metric]["noise"]
            limit = max(limit, noise)
            ax.fill_betweenx([low, high], -noise, noise, color="#e9edf1", zorder=0)
        ax.axvline(0, color=ink, linewidth=1.1)
        ax.axhline(3.25, color="#d7dce0", linewidth=0.9)
        for row in record["members"]:
            y = positions[(row["season"], row["member"])]
            delta = row[metric]["candidate"] - row[metric]["control"]
            ax.barh(y, delta, height=0.55, color=blue if delta < 0 else orange, zorder=2)
            ax.text(delta + (0.04 if delta >= 0 else -0.04) * limit, y,
                    f"{delta:+.4f}" if abs(delta) < 0.001 else f"{delta:+.3f}",
                    ha="left" if delta >= 0 else "right", va="center", fontsize=9, color=ink)
        for season, y in (("JAN", 6.55), ("JUL", 3.05)):
            mean = record["seasons"][season][metric]["delta_mean"]
            ax.plot([mean], [y], marker="D", color=ink, markersize=6, zorder=3)
        ax.set_yticks(list(positions.values()),
                      [f"{s.title()} · member {m}" for (s, m) in positions])
        ax.set_xlim(-1.55 * limit, 1.55 * limit)
        ax.set_ylim(-0.1, 7.0)
        ax.set_xlabel(xlabel)
        ax.set_title(title, loc="left")
        ax.grid(axis="x", color="#e3e7ea", linewidth=0.7)
        ax.set_axisbelow(True)

    ax7 = fig.add_subplot(grid[0, 0])
    bars(ax7, "onr7", "A  Original score (seven fields)",
         "Change in error, scaled by the CMIP6 median (v0.1.18 minus v0.1.17)")
    ax9 = fig.add_subplot(grid[0, 1])
    bars(ax9, "onr9", "B  Wider score (nine fields)",
         "Change in error, scaled by the CMIP6 median (v0.1.18 minus v0.1.17)")
    fig.legend(handles=[Line2D([0], [0], color=blue, linewidth=7, label="Smaller error"),
                        Line2D([0], [0], color=orange, linewidth=7, label="Larger error"),
                        Line2D([0], [0], marker="D", color=ink, linestyle="", label="Season mean"),
                        Line2D([0], [0], color="#e9edf1", linewidth=9, label="Spread of the matched v0.1.17 runs")],
               loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=4, frameon=False, fontsize=9)

    axg = fig.add_subplot(grid[1, 0])
    axg.axis("off")
    axg.set_title("C  Energy and safety checks", loc="left")
    names = {"olr_abs_bias": "Outgoing longwave, W m⁻²", "albedo_abs_bias": "Reflectivity",
             "net_abs_bias": "Top-of-atmosphere net, W m⁻²"}
    columns = ((0.56, "v0.1.17"), (0.69, "v0.1.18"),
               (0.84, "Allowed\nincrease"), (0.99, "Result"))
    y = 0.86
    axg.text(0.0, y, "Mean absolute error", fontweight="bold", transform=axg.transAxes, va="bottom")
    for x, head in columns:
        axg.text(x, y, head, fontweight="bold", transform=axg.transAxes,
                 ha="right" if head == "Result" else "center", fontsize=9, va="bottom")
    for name in GUARDS:
        y -= 0.15
        g = record["guards"][name]
        axg.text(0.0, y, names[name], transform=axg.transAxes)
        for x, value in ((0.56, g["control_mean"]), (0.69, g["candidate_mean"]), (0.84, g["allowance"])):
            axg.text(x, y, f"{value:.3f}", transform=axg.transAxes, ha="center", fontsize=9)
        axg.text(0.99, y, "pass" if g["passing"] else "fail", transform=axg.transAxes,
                 color=blue if g["passing"] else orange, fontweight="bold", ha="right", fontsize=9)
    y -= 0.15
    ok = record["guards"]["physical_breakage"]["passing"]
    axg.text(0.0, y, "Model state checks", transform=axg.transAxes)
    axg.text(0.99, y, "pass" if ok else "fail", transform=axg.transAxes,
             color=blue if ok else orange, fontweight="bold", ha="right", fontsize=9)

    axd = fig.add_subplot(grid[1, 1])
    axd.axis("off")
    axd.add_patch(FancyBboxPatch((0.0, 0.30), 1.0, 0.62, boxstyle="round,pad=0.0,rounding_size=0.02",
                                 transform=axd.transAxes, facecolor="#f4f6f8", edgecolor=grey, linewidth=1.4))
    axd.set_title("D  Release decision", loc="left")
    d = record["release_decision"]
    plain = {"NEUTRAL": "within spread", "WORSE": "worse", "BETTER": "better", "UNDECIDED": "no decision",
             "IMPROVED": "better", "REJECT": "reject", "ACCEPT": "accept", "NON_WORSE": "not worse"}
    word = lambda v: plain.get(str(v), str(v).lower())
    s3, s4 = record["seasons"], record["seasons"]
    jan4, jul4 = word(s4['JAN']['original_v4_verdict']), word(s4['JUL']['original_v4_verdict'])
    wider = (f"Wider score: {jan4} in either season" if jan4 == jul4 == "no decision" else
             f"Wider score: {jan4} in both seasons" if jan4 == jul4 else f"Wider score: January {jan4}, July {jul4}")
    lines = [("Scores", True),
             (f"Original score: January {word(s3['JAN']['original_v3_verdict'])}, "
              f"July {word(s3['JUL']['original_v3_verdict'])}, so {word(d['original_v3_metric_verdict'])}", False),
             (wider, False),
             ("Release rule", True),
             ("Pass if no season is clearly worse and every check passes.", False),
             ("Rule result: " + record["rule_result_shown"], True)]
    y = 0.86
    for text, bold in lines:
        axd.text(0.04, y, text, transform=axd.transAxes, fontweight="bold" if bold else "normal",
                 color=ink, fontsize=10.0 if bold else 9.2, va="top")
        y -= 0.093
    cost = record["run_cost"]
    # Figure-level text, so the layout engine does not shrink the panels to make room for it.
    fig.text(0.01, -0.045,
             "Controls: three matched runs per season with the v0.1.17 settings on an earlier commit "
             f"({record['control_source'][:7]}), not the published v0.1.17 tag. v0.1.18 commit {record['candidate_source'][:7]}.\n"
             "All runs start from observations. In a separate test, two simple guesses scored better than every CMIP6 model,\n"
             "so these scores show agreement with the observed month, not long-term climate skill. "
             f"The new per-step state check passed its GPU test ({record['runtime_guard_evidence']['cases_checked']} cases).\n"
             f"GPU runs took {cost['same_host_min_percent']:.0f} to {cost['same_host_max_percent']:.0f} percent longer "
             f"per model day than the matched v0.1.17 runs ({cost['same_host_count']} pairs on the same GPU host).",
             fontsize=8.8, color=grey, va="top", ha="left")
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build")
    b.add_argument("readout", type=Path)
    b.add_argument("record", type=Path)
    b.add_argument("--array", required=True)
    b.add_argument("--score-job", required=True)
    b.add_argument("--score-seal", required=True)
    b.add_argument("--guard-evidence", type=Path, required=True,
                   help="results folder of the real-GPU guard wrapper test for the candidate source")
    r = sub.add_parser("render")
    r.add_argument("record", type=Path)
    r.add_argument("output", type=Path)
    r.add_argument("--stamp-record", action="store_true",
                   help="write the image sha256 into the record after rendering")
    args = parser.parse_args()
    if args.command == "build":
        require(not args.record.exists(), "record already exists")
        require(len(args.score_seal) == 64, "score seal must be a sha256")
        record = build(args.readout, args.array, args.score_job, args.score_seal, args.guard_evidence)
        args.record.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    else:
        record = check_record(load(args.record))
        require(not args.output.exists(), "output exists")
        render(record, args.output)
        if args.stamp_record:
            record["image_sha256"] = sha(args.output)
            args.record.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
        print(sha(args.output))


if __name__ == "__main__":
    main()
