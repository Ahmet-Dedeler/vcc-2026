#!/usr/bin/env python3
"""Recompute a dated board review and core-input audit without fitting a new model.

Reads a caller-supplied public snapshot; writes only a compact evidence JSON.
No network, submissions, reserve-data access, training, or production code imports.
"""
from __future__ import annotations

import argparse
import ast
import datetime as dt
import json
from pathlib import Path

import numpy as np

from check_submission_novelty import sha256

ROOT = Path(__file__).resolve().parents[1]
MEMBERS = ("pds", "mse", "nmae", "fid", "reach", "jac")
IDS = ("Xv0ayvuoHpsKYWv2iu42", "Q519I2pLWRCW7PAEx5ul",
       "1n5LFk6PSEHq15UHENIE", "2K5ppxkowMSj2ugoZqYb", "xuDelEik7zGWMXdNu0R4")


def decompose(before, after):
    parts = {m: (after[f"score_{m}"] - before[f"score_{m}"]) / 6 for m in MEMBERS}
    delta = after["score_avg"] - before["score_avg"]
    assert abs(sum(parts.values()) - delta) < 1e-12
    return {"overall_delta": delta, "overall_contributions": parts}


def core_inputs():
    tree = ast.parse((ROOT / "scripts/make_submission_pm_cd4.py").read_text())
    frozen = next(ast.literal_eval(node.value) for node in tree.body
                  if isinstance(node, ast.Assign)
                  and any(isinstance(t, ast.Name) and t.id == "FROZEN_INPUT_SHA256"
                          for t in node.targets))
    rows = []
    # These aliases are the production builder's frozen_input_paths() mapping.
    roots = {"bundle": ROOT / "data/challenge", "donors": ROOT / "data/donors",
             "perturbmatch": ROOT / "data/external/perturbmatch",
             "cd4": ROOT / "data/external/cd4t/cloud-compact/cd4_compact"}
    bundle = json.loads((ROOT / "scripts/frozen_inputs.json").read_text())
    table = {**bundle["files"], **frozen}
    registry = next(ast.literal_eval(node.value) for node in tree.body
                    if isinstance(node, ast.AnnAssign)
                    and isinstance(node.target, ast.Name)
                    and node.target.id == "WHITENER_REGISTRY")
    # Match whitener_expected_hashes(), not the legacy full-union defaults.
    for suffix, digest in zip((".npz", ".json"), registry["union_pruned"]):
        table["whitener/whitener_union" + suffix] = digest
    for key, expected in table.items():
        kind, rest = key.split("/", 1)
        if kind == "whitener":
            suffix = ".npz" if rest.endswith(".npz") else ".json"
            path = ROOT / "data/donors" / ("whitener_union_pruned" + suffix)
        else:
            path = roots[kind] / rest
        actual = sha256(path) if path.is_file() else None
        rows.append({"key": key, "path": str(path.relative_to(ROOT)), "expected": expected,
                     "actual": actual, "match": actual == expected})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--board", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    board = json.loads(args.board.read_text())
    rows = sorted(board["entries"], key=lambda r: -r["score_avg"])
    ours = next(r for r in rows if r["team_name"] == "Ahmet Dedeler")
    report = {"checked_at": dt.datetime.now(dt.timezone.utc).isoformat(),
              "board_path": str(args.board), "board_sha256": sha256(args.board),
              "teams": len(rows), "rank": rows.index(ours) + 1, "entry": ours,
              "top_six": rows[:6], "member_ranks": {
                  m: 1 + sum(r[f"score_{m}"] > ours[f"score_{m}"] for r in rows)
                  for m in MEMBERS},
              "leader_gap": decompose(ours, rows[0]), "round7": []}
    captures = [json.loads((ROOT / f"data/board-entries/{i}.json").read_text()) for i in IDS]
    for cap in captures:
        entry = cap["entry"]
        report["round7"].append({"entry_id": entry["entry_id"],
                                "package": cap["package"], "scores": entry,
                                "versus_cd4svd5": decompose(captures[0]["entry"], entry),
                                "versus_plain15": decompose(captures[1]["entry"], entry)})
    team_id = ours["team_id"]
    history = {}
    for path in (ROOT / "data/leaderboard").glob("*.json"):
        snapshot = json.loads(path.read_text())
        for r in snapshot.get("entries", []):
            if r.get("team_id") == team_id and r.get("score_avg") is not None:
                history[r["id"]] = r
    history[ours["id"]] = ours
    report["best_recorded_entry"] = max(history.values(), key=lambda r: r["score_avg"])
    report["recorded_distinct_entries"] = len(history)
    previous_leader = next(r for r in json.loads(
        (ROOT / "data/leaderboard/2026-09-08T1253Z.json").read_text())["entries"]
                          if r["team_id"] == rows[0]["team_id"])
    report["leader_change_since_20260908T1253Z"] = decompose(previous_leader, rows[0])
    # A fixed-gallery counterexample to "higher PDS implies higher own-target cosine".
    truths = np.array([[1., 0.], [.8, .6]])
    high_cos = np.array([.9, np.sqrt(1 - .9**2)])
    lower_cos = np.array([.5, -np.sqrt(1 - .5**2)])
    a, b = truths @ high_cos, truths @ lower_cos
    assert a[0] > b[0] and a[1] > a[0] and b[0] > b[1]
    report["pds_is_not_cosine_counterexample"] = {
        "fixed_truths": truths.tolist(),
        "prediction_A": {"own_cosine": a[0], "decoy_cosine": a[1], "pds": 0},
        "prediction_B": {"own_cosine": b[0], "decoy_cosine": b[1], "pds": 1}}
    # Likewise a pre-operator row-norm identity does not survive an anisotropic map.
    operator = np.diag([2., 1.])
    report["pre_operator_norm_counterexample"] = {
        "input_norms": [1., 1.], "output_norms": np.linalg.norm(operator, axis=1).tolist()}
    report["core_inputs"] = core_inputs()
    report["all_core_inputs_match"] = all(r["match"] for r in report["core_inputs"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("rank", "teams", "member_ranks", "leader_gap",
                                           "all_core_inputs_match")}, indent=2))
    return int(not report["all_core_inputs_match"])


if __name__ == "__main__":
    raise SystemExit(main())
