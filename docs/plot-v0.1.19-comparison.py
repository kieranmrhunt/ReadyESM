#!/usr/bin/env python3
"""Build and render the v0.1.19 comparison from a verified six-pair score.

The graphic keeps the published v0.1.18 layout and score definitions. Its
record binds to the terminal exact-source v0.1.19 restart and paired scorer.
The original v3/v4 and separate neutral verdicts remain as scored.

  python3 -I -B plot-v0.1.19-comparison.py build RECORD.json \
      --array CANDIDATE_ARRAY --score-job SCORE_JOB --score-seal RESULT_SHA256 \
      --restart-job RESTART_JOB
  python3 -I -B plot-v0.1.19-comparison.py render RECORD.json OUTPUT.png --stamp-record

The controls are matched v0.1.18 scored months from source 48c232cc, not a
free-running climate reference. All runs begin from observations.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import tomllib

ROOT = Path('/gws/ssde/j25a/ncas_climate/vol1/users/kieran/readiesm')
SCORE_BUNDLE = ROOT / 'artifacts/development_wave_20261101/v019_paired_score_gate_stdout_fix_codex_v5'
SCORE_SEAL = '0a691e8e6674e83d216378ba5da0969a6caec342c6ba2a738e7f4f598a92161f'
SCORE_WORK = Path('/work/scratch-nopw2/kieran/readiesm_recovery/v019_paired_score_codex_v5/jobs')
RESTART_BUNDLE = ROOT / 'artifacts/development_wave_20261031/v019_exact_gpu_restart_codex_v3'
RESTART_SEAL = '2df0ff60e4ba267b0b30aac2ff5c36efc5020218ef6ff07859a02c6c940f4772'
RESTART_WORK = Path('/work/scratch-nopw2/kieran/readiesm_recovery/v019_exact_gpu_restart_codex_v3/jobs')
PYTHON = '/home/users/kieran/miniconda3/envs/py311/bin/python3'
SCHEMA = "readiesm_v019_comparison_figure_v1"
SEASONS = ("JAN", "JUL")
MEMBERS = (1, 2, 3)
CONTROL_LABEL = "matched v0.1.18 run"
CANDIDATE_LABEL = "v0.1.19"
CONTROL_SOURCE = "48c232cc6d313370c9d5c8fd2a52745af4e6c2ac"
CANDIDATE_SOURCE = "4caab6896ee4780f9fe283d03d955407bdfa3748"
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


def verified_score(array, score_job, expected_result):
    require(re.fullmatch(r'[1-9][0-9]*', str(score_job)) and
            re.fullmatch(r'[1-9][0-9]*', str(array)) and
            re.fullmatch(r'[0-9a-f]{64}', expected_result), 'score IDs or result SHA-256 invalid')
    require(sha(SCORE_BUNDLE / 'frozen_inputs.sha256') == SCORE_SEAL,
            'frozen scorer input seal differs')
    checked = subprocess.run([PYTHON, '-I', '-B', str(SCORE_BUNDLE / 'verify_score.py'),
                              str(score_job), str(array), SCORE_SEAL],
                             check=True, capture_output=True, text=True)
    marker = re.fullmatch(r'V019_PAIRED_SCORE_VERIFIED (ACCEPT_NEUTRAL_OR_IMPROVED|REJECT) '
                          r'result_manifest_sha256=([0-9a-f]{64})\n', checked.stdout)
    require(marker is not None and marker.group(2) == expected_result,
            'score verifier or result seal differs')
    work = SCORE_WORK / f'score_{score_job}_{array}'
    return work / 'readout', marker.group(1)


def load_restart_evidence(job):
    """Bind the figure's restart statement to the exact-source A100 gate."""
    require(re.fullmatch(r'[1-9][0-9]*', str(job)), 'restart job ID invalid')
    require(sha(RESTART_BUNDLE / 'frozen_inputs.sha256') == RESTART_SEAL,
            'A100 restart input seal differs')
    directory = RESTART_WORK / str(job)
    require(directory.is_dir() and not directory.is_symlink(), 'restart work directory missing')
    require((directory / 'frozen_manifest_sha256').read_text().strip() == RESTART_SEAL,
            'executed A100 restart input seal differs')
    accounting = subprocess.check_output(
        ['sacct', '-M', 'sakura', '-X', '-n', '-P', '-j', str(job),
         '--format=JobID,State,ExitCode'], text=True).splitlines()
    require(accounting == [f'{job}|COMPLETED|0:0'], 'A100 restart job is not terminal PASS')
    subprocess.run(['sha256sum', '--strict', '--quiet', '-c', 'result.sha256'],
                   cwd=directory, check=True)
    status = dict(line.split('=', 1) for line in (directory / 'status').read_text().splitlines())
    require(status == {'phase': 'complete', 'exit_code': '0'} and
            (directory / 'verdict').read_text().strip() == 'V019_RESTART_GPU_PASS',
            'A100 restart status or verdict differs')
    verified = subprocess.run([PYTHON, '-I', '-B', str(RESTART_BUNDLE / 'verify.py'), str(directory)],
                              check=True, capture_output=True, text=True)
    require(verified.stdout.splitlines()[-1:] == ['V019_RESTART_GPU_VERIFIED'],
            'A100 restart verifier did not PASS')
    reference = tomllib.loads((directory / 'reference/report.toml').read_text())
    resume = tomllib.loads((directory / 'resume/report.toml').read_text())
    require(reference['source_commit'] == resume['source_commit'] == CANDIDATE_SOURCE and
            reference['production_diagnostics_pass'] is True and
            resume['production_diagnostics_pass'] is True,
            'A100 restart source or production diagnostics differ')
    return {'job': str(job), 'source_commit': CANDIDATE_SOURCE,
            'verdict': 'V019_RESTART_GPU_PASS',
            'result_manifest_sha256': sha(directory / 'result.sha256'),
            'input_manifest_sha256': RESTART_SEAL,
            'reference_report_sha256': sha(directory / 'reference/report.toml'),
            'resume_report_sha256': sha(directory / 'resume/report.toml'),
            'production_diagnostics_pass': True,
            'scope': 'fresh-process exact-state restart and identity-refusal controls on one A100'}


def member_values(readout, season, member, version, field):
    path = Path(readout) / f"{season}_m{member}" / version / "scorecard.json"
    runs = load(path)["aggregate"]["runs"]
    return {"control": float(runs[f"c{member}"][field]), "candidate": float(runs[f"x{member}"][field])}, path


def run_cost(verification):
    """Report paired run rates, comparing only members run on one GPU host."""
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
    return {"method": "report.toml elapsed_day_seconds, day 1 to day 31; percent comparison only on the same GPU host",
            "pairs": pairs, "same_host_count": len(same),
            "same_host_percent_range": [min(same), max(same)] if len(same) >= 2 else None,
            "same_host_mean_percent": sum(same) / len(same) if len(same) >= 2 else None,
            "comparison_status": "comparable" if len(same) >= 2 else "insufficient same-host pairs"}


def build(array, score_job, score_seal, restart_job):
    require(sha(ROOT / 'artifacts' / RANK_BASELINE_SOURCE) == RANK_BASELINE_SHA256,
            'six-field baseline provenance changed')
    readout, verified_verdict = verified_score(array, score_job, score_seal)
    gate3 = load(readout / "version_gate_v3.json")
    gate4 = load(readout / "version_gate.json")
    physical = load(readout / "breakage.json")
    neutral = load(readout / "neutral_adjudication.json")
    treatment = load(readout / "treatment_gate.json")
    verification = load(readout / "verification.json")
    require(gate3.get("schema") == "readiesm_version_gate_v3", "wrong v3 gate schema")
    require(gate4.get("schema") == "readiesm_version_gate_v4_draft", "wrong v4 gate schema")
    require(neutral.get("schema") == "v018_neutral_skill_adjudication_v1", "wrong neutral schema")
    require(neutral["skill_adjudication"] in NEUTRAL_DECISIONS, "invalid neutral decision")
    require(neutral["skill_adjudication"] == verified_verdict,
            "readout decision differs from terminal scorer verifier")
    require(verification["sources"] == {"control": CONTROL_SOURCE, "candidate": CANDIDATE_SOURCE},
            "unexpected control or candidate source")
    require(str(verification["candidate_array"]) == str(array), "readout belongs to another candidate array")
    require(len(verification["pairs"]) == 6, "six pairs required")
    require(treatment.get('schema') == 'v019_paired_albedo_treatment_gate_v1' and
            treatment.get('passed') is True and len(treatment.get('rows', [])) == 6,
            'the six-month treatment-applied gate did not pass')
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
        "treatment_applied": {"passed": True, "members_checked": 6,
                              "receipt_sha256": sha(readout / "treatment_gate.json")},
        "release_decision": {
            "original_v3_metric_verdict": gate3["metric_verdict"],
            "original_v4_metric_verdict": gate4["metric_verdict"],
            "neutral_rule_decision": neutral["skill_adjudication"],
            "no_clear_worse": neutral["no_clear_worse"],
            "all_required_guards_pass": neutral["all_required_guards_pass"],
            "neutral_adjudication_sha256": sha(readout / "neutral_adjudication.json"),
        },
        "gpu_restart_evidence": load_restart_evidence(restart_job),
        "run_cost": run_cost(verification),
        "rule_result_shown": "pass" if neutral["skill_adjudication"] == "ACCEPT_NEUTRAL_OR_IMPROVED" else "reject",
        "release_status": "score result only; publication and final-tree checks are recorded separately",
        "renderer": "docs/plot-v0.1.19-comparison.py",
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
    evidence = record.get("gpu_restart_evidence") or {}
    require(evidence.get("verdict") == "V019_RESTART_GPU_PASS" and
            evidence.get("source_commit") == record["candidate_source"] and
            evidence.get("production_diagnostics_pass") is True and
            evidence.get("input_manifest_sha256") == RESTART_SEAL,
            "the record needs a passing exact-source A100 restart receipt")
    require(record.get("release_status") in (
            "score result only; publication and final-tree checks are recorded separately",
            "qualified release with failed frozen score"),
            "figure record must not imply publication")
    decision = record["release_decision"]
    require(decision["neutral_rule_decision"] in NEUTRAL_DECISIONS, "invalid neutral decision")
    require(set(record["guards"]) == set(GUARDS) | {"two_season_net", "physical_breakage"}, "incomplete guard set")
    require(record.get('treatment_applied', {}).get('passed') is True and
            record['treatment_applied'].get('members_checked') == 6 and
            re.fullmatch(r'[0-9a-f]{64}', record['treatment_applied'].get('receipt_sha256', '')),
            'missing six-month treatment-applied receipt')
    accept = decision["no_clear_worse"] is True and decision["all_required_guards_pass"] is True
    require((decision["neutral_rule_decision"] == "ACCEPT_NEUTRAL_OR_IMPROVED") == accept,
            "neutral decision contradicts its own flags")
    expected_shown = "pass" if decision["neutral_rule_decision"] == "ACCEPT_NEUTRAL_OR_IMPROVED" else "reject"
    require(record.get("rule_result_shown") == expected_shown,
            "the figure must state the release rule result exactly as scored")
    require("maintainer_acceptance" not in record, "the public figure record carries no maintainer statement")
    cost = record.get("run_cost") or {}
    require(len(cost.get("pairs", [])) == 6 and
            cost.get("same_host_count") == sum(p["same_host"] for p in cost["pairs"]),
            "missing or inconsistent run cost")
    if cost["same_host_count"] >= 2:
        bounds = cost["same_host_percent_range"]
        require(bounds[0] <= cost["same_host_mean_percent"] <= bounds[1] and
                cost["comparison_status"] == "comparable", "bad comparable run cost")
    else:
        require(cost["same_host_percent_range"] is None and
                cost["same_host_mean_percent"] is None and
                cost["comparison_status"] == "insufficient same-host pairs",
                "bad uncomparable run cost")
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
    fig.suptitle("ReadyESM v0.1.19 compared with matched v0.1.18 runs",
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
         "Change in error relative to 59 climate models (v0.1.19 minus v0.1.18)")
    ax9 = fig.add_subplot(grid[0, 1])
    bars(ax9, "onr9", "B  Wider score (nine fields)",
         "Change in error relative to 59 climate models (v0.1.19 minus v0.1.18)")
    fig.legend(handles=[Line2D([0], [0], color=blue, linewidth=7, label="Smaller error"),
                        Line2D([0], [0], color=orange, linewidth=7, label="Larger error"),
                        Line2D([0], [0], marker="D", color=ink, linestyle="", label="Season mean"),
                        Line2D([0], [0], color="#e9edf1", linewidth=9, label="Spread of the matched v0.1.18 runs")],
               loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=4, frameon=False, fontsize=9)

    axg = fig.add_subplot(grid[1, 0])
    axg.axis("off")
    axg.set_title("C  Energy and safety checks", loc="left")
    names = {"olr_abs_bias": "Outgoing longwave, W m⁻²", "albedo_abs_bias": "Reflectivity",
             "net_abs_bias": "Net energy at top of atmosphere, W m⁻²"}
    columns = ((0.56, "v0.1.18"), (0.69, "v0.1.19"),
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
    y -= 0.15
    axg.text(0.0, y, "Treatment applied, six months", transform=axg.transAxes)
    axg.text(0.99, y, "pass", transform=axg.transAxes,
             color=blue, fontweight="bold", ha="right", fontsize=9)

    axd = fig.add_subplot(grid[1, 1])
    axd.axis("off")
    axd.add_patch(FancyBboxPatch((0.0, 0.30), 1.0, 0.62, boxstyle="round,pad=0.0,rounding_size=0.02",
                                 transform=axd.transAxes, facecolor="#f4f6f8", edgecolor=grey, linewidth=1.4))
    axd.set_title("D  Release scoring rule", loc="left")
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
    # Figure-level text leaves the panel geometry unchanged.
    cost_note = (f"Time change on the same graphics card: {cost['same_host_percent_range'][0]:+.0f}% to "
                 f"{cost['same_host_percent_range'][1]:+.0f}% per model day "
                 f"({cost['same_host_count']} same-host pairs)."
                 if cost['comparison_status'] == 'comparable' else
                 "Graphics-card time comparison omitted: fewer than two pairs ran on the same host.")
    fig.text(0.01, -0.045,
             "Controls: three matched v0.1.18 scored runs per season (source "
             f"{record['control_source'][:7]}). v0.1.19 scored source {record['candidate_source'][:7]}.\n"
             "All runs start from observations. Simple forecasts beat all 59 comparison climate models\n"
             "in a separate six-field test. These scores show agreement over one month, not long-term climate skill. "
             "The save-and-restore test and model checks passed.\n" + cost_note,
             fontsize=8.8, color=grey, va="top", ha="left")
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build")
    b.add_argument("record", type=Path)
    b.add_argument("--array", required=True)
    b.add_argument("--score-job", required=True)
    b.add_argument("--score-seal", required=True)
    b.add_argument("--restart-job", required=True,
                   help="terminal exact-source A100 restart gate job")
    r = sub.add_parser("render")
    r.add_argument("record", type=Path)
    r.add_argument("output", type=Path)
    r.add_argument("--stamp-record", action="store_true",
                   help="write the image sha256 into the record after rendering")
    args = parser.parse_args()
    if args.command == "build":
        require(not args.record.exists(), "record already exists")
        require(len(args.score_seal) == 64, "score seal must be a sha256")
        record = build(args.array, args.score_job, args.score_seal, args.restart_job)
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
