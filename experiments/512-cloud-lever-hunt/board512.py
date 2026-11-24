#!/usr/bin/env python3
"""Board forensics for experiment 512.

Load every tracked leaderboard snapshot, dedupe entries by id, and reconstruct
per-team trajectories. Stdlib only. Writes work/forensics512.json.
"""
from __future__ import annotations

import datetime as dt
import glob
import json
import os
import re
from collections import defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(ROOT, "..", ".."))
OUT = os.path.join(ROOT, "work", "forensics512.json")

RAW = [
    ("pds", "pds_cosine"),
    ("mse", "expr_mse_unbiased_capped_norm"),
    ("nmae", "de_wilcoxon_lfc_nmae"),
    ("fid", "de_wilcoxon_direction_fidelity_yield_raw"),
    ("reach", "de_wilcoxon_direction_reach_raw"),
    ("jac", "de_wilcoxon_sig_jaccard"),
]
SCALED = [
    ("pds", "score_pds"),
    ("mse", "score_mse"),
    ("nmae", "score_nmae"),
    ("fid", "score_fid"),
    ("reach", "score_reach"),
    ("jac", "score_jac"),
]
# A raw member "moves" above this; below BIT it is treated as published-identical.
MOVE = {"pds": 5e-4, "mse": 1e-3, "nmae": 1e-3, "fid": 5e-4, "reach": 1e-3, "jac": 3e-4}
BIT = 1e-9


def snap_time(path: str) -> dt.datetime:
    name = os.path.basename(path)
    m = re.search(r"(20\d{2})-(\d{2})-(\d{2})T(\d{4})Z?", name)
    if m:
        stamp = "".join(m.groups())
    else:
        m = re.search(r"(20\d{6})T(\d{4})", name)
        if not m:
            raise ValueError(name)
        stamp = m.group(1) + m.group(2)
    return dt.datetime.strptime(stamp, "%Y%m%d%H%M").replace(tzinfo=dt.timezone.utc)


def fnum(x):
    if x is None:
        return None
    return float(x)


def load():
    paths = sorted(glob.glob(os.path.join(REPO, "data", "leaderboard", "*.json")), key=snap_time)
    by_id = {}
    n_rows = 0
    for p in paths:
        with open(p) as f:
            doc = json.load(f)
        entries = doc["entries"] if isinstance(doc, dict) else doc
        n_rows += len(entries)
        ts = snap_time(p).strftime("%Y-%m-%dT%H:%MZ")
        for e in entries:
            eid = e["id"]
            rec = by_id.get(eid)
            if rec is None:
                slim = {
                    "id": eid,
                    "team_id": e.get("team_id"),
                    "team_name": e.get("team_name"),
                    "organization": e.get("organization"),
                    "members": e.get("team_members") or [],
                    "submission_date": e.get("submission_date"),
                    "submission_count": e.get("submission_count"),
                    "model_name": e.get("model_name"),
                    "description": e.get("description") or "",
                    "score_avg": fnum(e.get("score_avg")),
                    "first_snap": ts,
                    "last_snap": ts,
                    "n_snaps": 1,
                }
                for short, key in RAW + SCALED:
                    # scaled keys collide on short name; store separately
                    pass
                for short, key in RAW:
                    slim["raw_" + short] = fnum(e.get(key))
                for short, key in SCALED:
                    slim["sc_" + short] = fnum(e.get(key))
                by_id[eid] = slim
            else:
                rec["last_snap"] = ts
                rec["n_snaps"] += 1
                # submission_count can rise while the same id is displayed; keep the max
                sc = e.get("submission_count")
                if sc is not None and (rec["submission_count"] is None or sc > rec["submission_count"]):
                    rec["submission_count"] = sc
    return paths, n_rows, by_id


def parse_sub(s: str) -> dt.datetime:
    # 2026-09-27T01:44:07.865265Z
    if s.endswith("Z"):
        s = s[:-1] + "+00:00"
    return dt.datetime.fromisoformat(s)


def delta_pattern(a, b):
    moved, flat, bits = [], [], []
    draw, dsc = {}, {}
    for short, _ in RAW:
        d = b["raw_" + short] - a["raw_" + short]
        draw[short] = d
        if abs(d) <= BIT:
            bits.append(short)
            flat.append(short)
        elif abs(d) < MOVE[short]:
            flat.append(short)
        else:
            moved.append(short)
    for short, _ in SCALED:
        dsc[short] = b["sc_" + short] - a["sc_" + short]
    dsc["overall"] = b["score_avg"] - a["score_avg"]
    # identity of the six-mean
    mean6 = sum(b["sc_" + s] for s, _ in SCALED) / 6.0
    return {
        "moved": moved,
        "flat_or_bit": flat,
        "bit": bits,
        "draw": draw,
        "dsc": dsc,
        "score_avg_minus_mean6": b["score_avg"] - mean6,
    }


def pattern_code(moved):
    order = [s for s, _ in RAW]
    return "".join(s[0] if s in moved else "." for s in order)


def main():
    paths, n_rows, by_id = load()
    entries = list(by_id.values())
    # group by team_id
    teams = defaultdict(list)
    for e in entries:
        teams[e["team_id"]].append(e)
    for tid, rows in teams.items():
        rows.sort(key=lambda r: (parse_sub(r["submission_date"]), r["id"]))

    # latest entry per team
    latest = []
    for tid, rows in teams.items():
        last = max(rows, key=lambda r: parse_sub(r["submission_date"]))
        latest.append(last)
    latest.sort(key=lambda r: -r["score_avg"])

    # member-set identity: same people => same human team even if team_id differs
    def people_key(r):
        return tuple(sorted(r["members"])) if r["members"] else (r["team_name"],)

    # top 10 distinct team_ids
    top10 = latest[:10]

    # also a distinct-people ranking
    seen_people = set()
    distinct_people = []
    for r in latest:
        k = people_key(r)
        if k in seen_people:
            continue
        seen_people.add(k)
        distinct_people.append(r)
        if len(distinct_people) >= 12:
            break

    # find our team
    ours = [r for r in latest if "dedeler" in r["team_name"].lower() or any("dedeler" in (m or "").lower() for m in r["members"])]

    def pack(rows_latest, ref, tol=0.008):
        """Teams whose latest raw vector is close to ref."""
        out = []
        for r in rows_latest:
            if r["team_id"] == ref["team_id"]:
                continue
            dists = [abs(r["raw_" + s] - ref["raw_" + s]) for s, _ in RAW]
            if max(dists) <= tol and abs(r["score_avg"] - ref["score_avg"]) <= 0.02:
                out.append({"team": r["team_name"], "overall": r["score_avg"], "max_abs_raw": max(dists)})
        return out

    trajectories = []
    interesting = []
    for r0 in top10 + ours:
        rows = teams[r0["team_id"]]
        steps = []
        # first entry
        first = rows[0]
        steps.append({
            "kind": "first",
            "date": first["submission_date"],
            "id": first["id"],
            "model": first["model_name"],
            "desc": (first["description"] or "")[:400],
            "overall": first["score_avg"],
            "raw": {s: first["raw_" + s] for s, _ in RAW},
            "scaled": {s: first["sc_" + s] for s, _ in SCALED},
            "submission_count": first["submission_count"],
            "captured_n": len(rows),
        })
        for a, b in zip(rows, rows[1:]):
            pat = delta_pattern(a, b)
            rec = {
                "kind": "step",
                "date": b["submission_date"],
                "from_model": a["model_name"],
                "to_model": b["model_name"],
                "from_id": a["id"],
                "to_id": b["id"],
                "desc": (b["description"] or "")[:400],
                "overall": b["score_avg"],
                "d_overall": pat["dsc"]["overall"],
                "code": pattern_code(pat["moved"]),
                "moved": pat["moved"],
                "bit": pat["bit"],
                "draw": {k: round(v, 6) for k, v in pat["draw"].items()},
                "dsc": {k: round(v, 6) for k, v in pat["dsc"].items()},
                "mixed": bool(pat["bit"]) and bool(pat["moved"]),
            }
            steps.append(rec)
            if rec["mixed"] or (pat["bit"] and any(abs(pat["draw"][s]) >= 0.0005 for s in pat["draw"] if s not in pat["bit"])):
                interesting.append({
                    "team": r0["team_name"],
                    "date": b["submission_date"][:16],
                    "models": f"{a['model_name']} -> {b['model_name']}",
                    "code": rec["code"],
                    "bit": pat["bit"],
                    "moved": pat["moved"],
                    "draw": rec["draw"],
                    "dsc": rec["dsc"],
                    "d_overall": round(pat["dsc"]["overall"], 5),
                    "desc": (b["description"] or "")[:240],
                })
        # collapse consecutive identical-code steps for the summary trajectory
        trajectories.append({
            "team": r0["team_name"],
            "team_id": r0["team_id"],
            "n_captured": len(rows),
            "submission_count_latest": rows[-1]["submission_count"],
            "latest_overall": rows[-1]["score_avg"],
            "latest_model": rows[-1]["model_name"],
            "first_date": rows[0]["submission_date"],
            "last_date": rows[-1]["submission_date"],
            "steps": steps,
        })

    # gap vs our best captured (max score) and vs latest
    our_best = max(ours, key=lambda r: r["score_avg"]) if ours else None
    # actually best ENTRY not latest team row
    our_entries = [e for e in entries if e["team_id"] in {r["team_id"] for r in ours}] if ours else []
    our_best_entry = max(our_entries, key=lambda r: r["score_avg"]) if our_entries else None

    def gap_row(other, base):
        if other is None or base is None:
            return None
        d_sc = {s: other["sc_" + s] - base["sc_" + s] for s, _ in SCALED}
        return {
            "team": other["team_name"],
            "model": other["model_name"],
            "date": other["submission_date"],
            "overall": other["score_avg"],
            "d_overall": other["score_avg"] - base["score_avg"],
            "raw": {s: other["raw_" + s] for s, _ in RAW},
            "d_raw": {s: other["raw_" + s] - base["raw_" + s] for s, _ in RAW},
            "d_scaled": d_sc,
            "d_overall_points": {s: d_sc[s] / 6.0 for s in d_sc},
            "de_points": sum(d_sc[s] / 6.0 for s in ("nmae", "fid", "reach", "jac")),
            "bulk_points": sum(d_sc[s] / 6.0 for s in ("pds", "mse")),
        }

    # score identity check
    max_mean6_err = 0.0
    n_scored = 0
    for e in entries:
        if e["score_avg"] is None or e["sc_pds"] is None:
            continue
        mean6 = sum(e["sc_" + s] for s, _ in SCALED) / 6.0
        max_mean6_err = max(max_mean6_err, abs(e["score_avg"] - mean6))
        n_scored += 1

    # global mixed-step census at BIT and MOVE, all teams, not just top 10
    census = defaultdict(int)
    n_steps_all = 0
    n_mixed_all = 0
    examples_nmae_only = []
    for tid, rows in teams.items():
        for a, b in zip(rows, rows[1:]):
            n_steps_all += 1
            pat = delta_pattern(a, b)
            if pat["bit"] and pat["moved"]:
                n_mixed_all += 1
                key = tuple(sorted(pat["moved"]))
                census["+".join(pat["moved"]) if pat["moved"] else "none"] += 1
            # nmae moved >= 0.001, pds mse fid jac all bit
            if (
                abs(pat["draw"]["nmae"]) >= 0.001
                and all(s in pat["bit"] for s in ("pds", "mse", "fid", "jac"))
            ):
                examples_nmae_only.append({
                    "team": b["team_name"],
                    "date": b["submission_date"][:16],
                    "models": f"{a['model_name']} -> {b['model_name']}",
                    "dnmae": round(pat["draw"]["nmae"], 5),
                    "dreach": round(pat["draw"]["reach"], 5),
                    "d_overall": round(pat["dsc"]["overall"], 5),
                    "reach_bit": "reach" in pat["bit"],
                })

    # chill descriptions: every captured entry with a description
    def team_rows_by_name(substr):
        hits = []
        for tid, rows in teams.items():
            if rows and substr.lower() in rows[-1]["team_name"].lower():
                hits.extend(rows)
        return hits

    def desc_dump(substr):
        rows = team_rows_by_name(substr)
        out = []
        for r in rows:
            if r["description"]:
                out.append({
                    "date": r["submission_date"][:16],
                    "model": r["model_name"],
                    "overall": r["score_avg"],
                    "raw_nmae": r["raw_nmae"],
                    "raw_pds": r["raw_pds"],
                    "raw_mse": r["raw_mse"],
                    "desc": r["description"][:500],
                })
        return out

    # latest raw table for top 10 + ours
    table = []
    for r in top10:
        table.append({
            "team": r["team_name"],
            "model": r["model_name"],
            "date": r["submission_date"],
            "overall": r["score_avg"],
            "raw": {s: r["raw_" + s] for s, _ in RAW},
            "scaled": {s: r["sc_" + s] for s, _ in SCALED},
            "n_captured": len(teams[r["team_id"]]),
            "submission_count": r["submission_count"],
            "org": r["organization"],
        })

    doc = {
        "n_snapshots": len(paths),
        "snapshot_first": snap_time(paths[0]).strftime("%Y-%m-%dT%H:%MZ"),
        "snapshot_last": snap_time(paths[-1]).strftime("%Y-%m-%dT%H:%MZ"),
        "n_rows_with_dupes": n_rows,
        "n_unique_entries": len(entries),
        "n_teams": len(teams),
        "n_scored_entries": n_scored,
        "max_abs_score_avg_minus_mean6": max_mean6_err,
        "move_thresholds": MOVE,
        "bit_threshold": BIT,
        "top10_latest": table,
        "distinct_people_top": [
            {"team": r["team_name"], "overall": r["score_avg"], "people": r["members"][:8]}
            for r in distinct_people
        ],
        # Filled after this dict: latest captured row for our team. Gaps below are
        # against that displayed row, which is not 446bj2 (see RESULTS.md).
        "ours_latest": None,
        "our_best_entry": None if not our_best_entry else {
            "id": our_best_entry["id"],
            "model": our_best_entry["model_name"],
            "date": our_best_entry["submission_date"],
            "overall": our_best_entry["score_avg"],
            "raw": {s: our_best_entry["raw_" + s] for s, _ in RAW},
            "scaled": {s: our_best_entry["sc_" + s] for s, _ in SCALED},
        },
        "gaps_vs_our_best": [
            gap_row(r, our_best_entry) for r in top10
        ],
        "n_steps_all_teams": n_steps_all,
        "n_mixed_bit_and_move": n_mixed_all,
        "mixed_moved_set_counts": dict(sorted(census.items(), key=lambda kv: -kv[1])[:40]),
        "nmae_only_pds_mse_fid_jac_bit": examples_nmae_only,
        "interesting_top10": interesting,
        "trajectories": trajectories,
        "descriptions": {
            "chill": desc_dump("chill"),
            "testh": desc_dump("testh"),
            "cqawesome": desc_dump("cqawesome"),
        },
    }
    # fix ours_latest properly
    if ours:
        cur = max(ours, key=lambda r: parse_sub(r["submission_date"]))
        doc["ours_latest"] = {
            "team": cur["team_name"],
            "team_id": cur["team_id"],
            "model": cur["model_name"],
            "date": cur["submission_date"],
            "overall": cur["score_avg"],
            "raw": {s: cur["raw_" + s] for s, _ in RAW},
            "scaled": {s: cur["sc_" + s] for s, _ in SCALED},
            "n_captured": len(teams[cur["team_id"]]),
            "submission_count": cur["submission_count"],
        }

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(doc, f)
    print(json.dumps({
        "snapshots": doc["n_snapshots"],
        "range": [doc["snapshot_first"], doc["snapshot_last"]],
        "unique": doc["n_unique_entries"],
        "teams": doc["n_teams"],
        "mean6_err": doc["max_abs_score_avg_minus_mean6"],
        "steps": n_steps_all,
        "mixed": n_mixed_all,
        "nmae_only": len(examples_nmae_only),
        "top1": table[0]["team"] if table else None,
        "our_best": None if not our_best_entry else our_best_entry["score_avg"],
    }, indent=2))


if __name__ == "__main__":
    main()
