#!/usr/bin/env python3
"""638b - CLEAN4-ext: CLEAN4's own definition (unconditional P_down, equal-weight mean of k562_gwps, k562_essential,
nadig_jurkat, nadig_hepg2) recomputed on each bank's FULL gene axis instead of 371's 7-bank intersection (6,346 genes).
Question: on gate genes that 371's CLEAN4 does NOT cover, how accurate are CLEAN4-ext's signs (at the hook's tau 0.05)?
LOLO per 636 (K562 banks dropped on the k562 surface). Analysis only."""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
from inst636 import BANK_FORBIDDEN, Inst, clean4_lolo, bank_uncond
I = Inst()
U = {b: bank_uncond(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")}
U = {b: v for b, v in U.items() if v[0] is not None}
print("banks with delta on disk:", list(U))
rows = []
for s in I.S:
    C4, used = clean4_lolo(I, s)
    use = [b for b in U if s not in BANK_FORBIDDEN[b]]
    EXT = np.nanmean(np.stack([I.gene_vec(s, *U[b]) for b in use]), 0)
    D = I.S[s]; gate = D["gate"]; truth_dn = D["lfc"] < 0
    for tau in (0.0, 0.05):
        old = np.isfinite(C4) & (np.abs(C4 - 0.5) >= max(tau, 1e-12))
        new = ~np.isfinite(C4) & np.isfinite(EXT) & (np.abs(EXT - 0.5) >= max(tau, 1e-12))
        for lab, P, msk in (("CLEAN4 (371) on its coverage", C4, old), ("CLEAN4-ext on NEWLY covered genes", EXT, new),
                            ("all-down on newly covered genes", np.full(len(D["G"]), 0.9), new)):
            m = gate & msk[None, :]
            ok = ((P > 0.5)[None, :] == truth_dn) & m
            w = np.abs(np.nan_to_num(D["lfc"])) * m
            rows.append(dict(surface=s, tau=tau, source=lab, banks=use, gate_cov=float(m.sum() / gate.sum()), n=int(m.sum()),
                             acc=float(ok.sum() / max(m.sum(), 1)), acc_w=float((w * ok).sum() / max(w.sum(), 1e-12)),
                             gain_units=float(((2 * ok - 1) * m * np.abs(np.nan_to_num(D["lfc"]))).sum() / max(np.abs(np.nan_to_num(D["lfc"]) * gate).sum(), 1e-12))))
print(f"{'surface':8s} {'tau':>4s} {'source':36s} {'gate_cov':>8s} {'n':>7s} {'acc':>6s} {'acc_w':>6s} {'(2a-1)|l| share':>15s}")
for r in rows:
    print(f"{r['surface']:8s} {r['tau']:4.2f} {r['source']:36s} {r['gate_cov']:8.3f} {r['n']:7d} {r['acc']:6.4f} {r['acc_w']:6.4f} {r['gain_units']:15.4f}")
(HERE / "work/ext638.json").write_text(json.dumps(rows, indent=1, default=float))
