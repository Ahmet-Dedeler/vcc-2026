"""470 side analysis: across public CRISPRi banks, do WEAK perturbations move their changed genes mostly DOWN?

H1 (467): the reference DE genes of the weakest third of h1-test targets are down 63 % of the time per target (strong:
55 %), and the four reach-loss targets' are ~all down. If weak CRISPRi responses are down-dominated in every line, a
down-leaning prior for design-weak targets should transfer to A/B/C and October.

Proxy DE set per perturbation (banks hold pseudobulk log-ratios, not p-values): genes among the bank's 2,000 most
expressed controls (b_ctrl) with |delta| > 0.3 (~4 sd of the pseudobulk noise at >= 1 count/cell and >= 75 cells), the
target itself excluded. Weak = 1..20 such genes, mid = 21..100, strong = > 100.
usage: uv run python experiments/470-weak-target-prior/downfrac470.py
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BANKS = ["k562_gwps", "rpe1", "nadig_jurkat", "orion_hct116", "orion_hek293t", "h1_2025_full"]
out = {}
for bank in BANKS:
    d = ROOT / "data/donors" / bank
    delta = np.load(d / "delta.npy", mmap_mode="r")
    perts = [l.strip() for l in open(d / "perts.txt")]
    genes = np.array([l.strip() for l in open(d / "genes.txt")])
    bpath = d / "b_ctrl.npy" if (d / "b_ctrl.npy").exists() else ROOT / "data/donors/h1_2025/b_ctrl.npy"
    b = np.load(bpath)
    assert b.shape[0] == genes.shape[0], (bank, b.shape, genes.shape)
    top = np.argsort(-b)[:2000]
    gtop = genes[top]
    rows = []
    for s in range(0, len(perts), 2000):
        X = np.asarray(delta[s:s + 2000][:, top], dtype=np.float32)
        for i in range(X.shape[0]):
            p = perts[s + i]
            m = np.abs(X[i]) > 0.3
            m &= gtop != p
            n = int(m.sum())
            if n == 0:
                continue
            rows.append((n, int((X[i][m] < 0).sum())))
    r = np.array(rows)
    res = {}
    for name, lo, hi in (("weak", 1, 20), ("mid", 21, 100), ("strong", 101, 10**9)):
        sel = (r[:, 0] >= lo) & (r[:, 0] <= hi)
        if sel.sum() == 0:
            continue
        per = r[sel, 1] / r[sel, 0]
        res[name] = {"perturbations": int(sel.sum()), "mean_per_pert_down": round(float(per.mean()), 3),
                     "pooled_down": round(float(r[sel, 1].sum() / r[sel, 0].sum()), 3),
                     "share_majority_down": round(float((per > 0.5).mean()), 3)}
    out[bank] = res
    print(bank, json.dumps(res), flush=True)
(Path(__file__).parent / "work/downfrac470.json").write_text(json.dumps(out, indent=1))
