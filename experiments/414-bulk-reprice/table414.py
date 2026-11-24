"""Render the ranked markdown tables for RESULTS.md from work/analyse414.json (no numbers typed by hand)."""
import json
from pathlib import Path

import numpy as np

W = Path(__file__).resolve().parent / "work"
an = json.loads((W / "analyse414.json").read_text())
K = an["plolo_kolf2"]
ORDER = ["Kinf+sh@1.5", "Kinf+sh@1.35", "K0.5_inf+sh@1.25", "Kinf+sh@1.25", "K0.5_inf+sh@1.15", "Kinf+sh@1.15",
         "K3+sh@1.15", "Kinf+shL@1.15", "Kinf+sh@1.1", "Kinf+sh@1.05", "K0.7_inf", "K0.5_inf", "K0.4_inf", "K0.5_3",
         "Kinf", "K3", "Kinf+gf@0.1", "Kinf+gf@0.25", "gi@0.1+Kinf", "gi@0.25+Kinf", "Kinf+gf@0.5", "gi@0.5+Kinf",
         "Kinf+gf@1", "gi@1+Kinf", "Kinf+gf@0.25+sh@1.15", "Kinf@x0.97", "Kinf@x1.03",
         "Kinf+shperm@1.15", "Kinf+sh@0.87", "Kinf_shuf", "K3_shuf", "Kinf+gf_shuf@0.25", "sh@1.15", "gf@0.25", "base"]
print("| arm (bulk only, on 410c) | rel cos vs base | Δrel vs 410c | Δpds literal | Δpds board-emul. | Δraw mse 197 (k_key) | "
      "Δraw mse cliff | Δoverall mandate | **Δoverall primary** | Δoverall if levelled | render A ratio | folds + (rel/pds) | kolf2 Δrel / Δpds |")
print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|")
for a in ORDER:
    if a not in an["table"]:
        continue
    r = an["table"][a]["vs410c"]
    ref = an["plolo_cancer"]["Kinf"]
    nrel = sum(x - y > 0 for x, y in zip(an["plolo_cancer"][a]["rel_f"], ref["rel_f"]))
    npd = sum(x > 0 for x in r["dpds_board_f"] if np.isfinite(x))
    kf = ""
    if a in K:
        kf = f"{(1 + K[a]['rel']) / (1 + K['Kinf']['rel']) - 1:+.3f} / {K[a]['dpds'] - K['Kinf']['dpds']:+.4f}"
    pb = f"{r['dpds_board']:+.4f}" if np.isfinite(r["dpds_board"]) else "–"
    print(f"| `{a}` | {r['rel']:+.3f} | {r['d_rel']:+.3f} | {r['dpds']:+.4f} | {pb} | {r['d197']:+.4f} | {r['dcliff']:+.4f} | "
          f"{r['ov_mandate']:+.4f} | **{r['ov_primary']:+.4f}** | {r['ov_lev']:+.4f} | {np.mean(list(r['Ar'].values())):.3f} | "
          f"{nrel}/3, {npd if np.isfinite(r['dpds_board']) else '–'}/3 | {kf} |")
