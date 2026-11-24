"""480 resume-2: what the board's raw members say about where rivals' nmae comes from (context, not a test).

Reads every data/leaderboard/*.json snapshot, deduplicates entries by id, and reports:
  1. the latest entry of the top DE teams with their raw DE members (nmae, fid, reach, jac) next to ours;
  2. "nmae-only" within-team steps: consecutive entries of one team whose pds, mse, fid and jac are all bit-identical
     (bulk, call set and the signs of called genes unchanged) while nmae moves by >= 0.001 raw. Such a step can only
     change called genes' magnitudes or uncalled genes' arithmetic means (sign or magnitude); the board cannot say which.
  3. "DE-only" steps: pds and mse bit-identical, |Δnmae| >= 0.005 (the DE channel moved, anything in it).

    uv run python board480.py     # -> work/board480.json
"""
from __future__ import annotations

import glob
import json
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
M = ["pds_cosine", "expr_mse_unbiased_capped_norm", "de_wilcoxon_lfc_nmae", "de_wilcoxon_direction_fidelity_yield_raw",
     "de_wilcoxon_direction_reach_raw", "de_wilcoxon_sig_jaccard"]
SHORT = dict(zip(M, ["pds", "mse", "nmae", "fid", "reach", "jac"]))
TEAMS = ["Illumina AI", "cqawesome", "CellWM", "chill", "Testh", "nailong king", "Ava Brown", "Anna Xie", "Vivian Xu",
         "SQ W", "HY H", "Ahmet Dedeler"]

ents = {}
for f in sorted(glob.glob(str(ROOT / "data/leaderboard/*.json"))):
    try:
        d = json.load(open(f))
    except Exception:
        continue
    for e in (d.get("entries", []) if isinstance(d, dict) else d):
        if isinstance(e, dict) and "id" in e and all(isinstance(e.get(m), float) for m in M):
            ents[e["id"]] = e
by_team = defaultdict(list)
for e in ents.values():
    by_team[e["team_name"]].append(e)
for t in by_team:
    by_team[t].sort(key=lambda e: e["submission_date"])
print(f"{len(ents)} unique scored entries, {len(by_team)} teams")

latest = {}
print("\n1. latest entry, raw DE members")
print(f"{'team':16s} {'date':16s} {'overall':>7s} {'nmae':>6s} {'fid':>6s} {'reach':>6s} {'jac':>6s} {'pds':>6s} {'mse':>6s}")
for t in TEAMS:
    if t not in by_team:
        continue
    e = by_team[t][-1]
    latest[t] = {SHORT[m]: e[m] for m in M} | dict(date=e["submission_date"], overall=e["score_avg"],
                                                     model=e.get("model_name"))
    print(f"{t:16s} {e['submission_date'][:16]:16s} {e['score_avg']:7.4f} " +
          " ".join(f"{e[m]:6.4f}" for m in M[2:] + M[:2]))

steps_only, steps_de = [], []
for t, L in by_team.items():
    for a, b in zip(L, L[1:]):
        same_bulk = a[M[0]] == b[M[0]] and a[M[1]] == b[M[1]]
        dn = b[M[2]] - a[M[2]]
        if not same_bulk:
            continue
        rec = dict(team=t, date=b["submission_date"][:16], model_from=a.get("model_name"), model_to=b.get("model_name"),
                   nmae_from=a[M[2]], d_nmae=dn, d_fid=b[M[3]] - a[M[3]], d_reach=b[M[4]] - a[M[4]],
                   d_jac=b[M[5]] - a[M[5]])
        if a[M[3]] == b[M[3]] and a[M[5]] == b[M[5]] and abs(dn) >= 0.001:
            steps_only.append(rec)
        if abs(dn) >= 0.005:
            steps_de.append(rec)

print("\n2. nmae-only steps (pds, mse, fid, jac bit-identical; |Δnmae| >= 0.001), all teams")
for r in sorted(steps_only, key=lambda r: r["d_nmae"]):
    print(f"  {r['team']:18s} {r['date']} {str(r['model_from'])[:22]:22s} -> {str(r['model_to'])[:22]:22s} "
          f"nmae {r['nmae_from']:.4f} {r['d_nmae']:+.4f}  reach {r['d_reach']:+.4f}")
print("\n3. DE-only steps (pds, mse bit-identical; |Δnmae| >= 0.005) for the listed teams")
for r in sorted([r for r in steps_de if r["team"] in TEAMS], key=lambda r: (r["team"], r["date"])):
    print(f"  {r['team']:18s} {r['date']} nmae {r['nmae_from']:.4f} {r['d_nmae']:+.4f}  fid {r['d_fid']:+.4f} "
          f"reach {r['d_reach']:+.4f} jac {r['d_jac']:+.4f}  ({str(r['model_to'])[:24]})")
json.dump(dict(n_entries=len(ents), latest=latest, nmae_only_steps=steps_only,
               de_only_steps=[r for r in steps_de if r["team"] in TEAMS]),
          open(HERE / "work/board480.json", "w"), indent=1)
