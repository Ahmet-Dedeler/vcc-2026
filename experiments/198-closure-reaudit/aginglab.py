#!/usr/bin/env python3
"""198 — adjudicate the AgingLab GeroAI v9 -> v11 reach move against 172's null.

Claim in the mandate: v9->v11 gained +0.1002 raw board `reach` with `pds` going DOWN
0.008 -- the largest ordering-shaped move on the board -- while 172 tried their
published feature family and got nothing. Reconcile.
"""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "work"
OUT.mkdir(exist_ok=True)

rows = json.loads((ROOT / "experiments/186-intel-0903/board/leaderboard_2026-09-03.json"
                   ).read_text())["entries"]

M = [("pds_cosine", "pds"), ("expr_mse_unbiased_capped_norm", "mse"),
     ("de_wilcoxon_lfc_nmae", "nmae"),
     ("de_wilcoxon_direction_fidelity_yield_raw", "fid"),
     ("de_wilcoxon_direction_reach_raw", "reach"),
     ("de_wilcoxon_sig_jaccard", "jac")]
S = [("score_pds", "pds"), ("score_mse", "mse"), ("score_nmae", "nmae"),
     ("score_fid", "fid"), ("score_reach", "reach"), ("score_jac", "jac")]


def show(sel, label):
    sel = sorted(sel, key=lambda r: r["submission_date"])
    print(f"\n=== {label}  ({len(sel)} entries) ===")
    print(f"{'date':11s} {'model':26s} {'overall':>9s} " +
          " ".join(f"{n:>8s}" for _, n in M))
    for r in sel:
        print(f"{r['submission_date'][:10]:11s} {str(r['model_name'])[:26]:26s} "
              f"{r['score_avg']:9.5f} " +
              " ".join(f"{r[k]:8.4f}" if r.get(k) is not None else "     ---"
                       for k, _ in M))
    print(f"{'':11s} {'--- scaled ---':26s} {'':9s} " +
          " ".join(f"{n:>8s}" for _, n in S))
    for r in sel:
        print(f"{r['submission_date'][:10]:11s} {str(r['model_name'])[:26]:26s} "
              f"{r['score_avg']:9.5f} " +
              " ".join(f"{r[k]:8.4f}" if r.get(k) is not None else "     ---"
                       for k, _ in S))
    return sel


ag = [r for r in rows if "aging" in str(r.get("team_name", "")).lower()
      or "gero" in str(r.get("model_name", "")).lower()]
sel = show(ag, "AgingLab / GeroAI")

# the specific pair named in the mandate
by = {str(r["model_name"]): r for r in sel}
print("\n=== the v9 -> v11 pair ===")
cands = {k: v for k, v in by.items() if "v9" in k or "v10" in k or "v11" in k}
for k in sorted(cands):
    r = cands[k]
    print(f"  {k}: reach {r['de_wilcoxon_direction_reach_raw']:.4f}  "
          f"pds {r['pds_cosine']:.4f}  mse {r['expr_mse_unbiased_capped_norm']:.4f}  "
          f"fid {r['de_wilcoxon_direction_fidelity_yield_raw']:.4f}  "
          f"jac {r['de_wilcoxon_sig_jaccard']:.4f}  "
          f"nmae {r['de_wilcoxon_lfc_nmae']:.4f}  overall {r['score_avg']:.5f}")

json.dump(sel, open(OUT / "aginglab_trajectory.json", "w"), indent=1)
