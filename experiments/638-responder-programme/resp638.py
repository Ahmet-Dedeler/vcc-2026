#!/usr/bin/env python3
"""638 — responder-weighted programmes as gate-sign sources, scored leave-one-LINE-out on 636's instrument.

Idea: gate genes are responders, but CLEAN4's P_down (371) counts every perturbation of a bank equally, most of them
non-responders. A programme taken over STRONG perturbations (the bank's own response strength, no truth used) should be
closer to "which way does g go when the cell responds". It weights perturbations, not gene-level DE events (634 tested
the latter and lost).

Per bank b (only banks allowed on the surface; 636's BANK_FORBIDDEN), with s_p = RMS of perturbation p's delta over the
bank's genes:
  U   : P_down = mean_p 1[delta_pg < 0]                       (CLEAN4's per-bank definition)
  W1  : P_down weighted by s_p                                (responder-weighted)
  W2  : P_down weighted by s_p^2
  TOPq: P_down over the top q (10 %, 25 %) of perturbations by s_p
  MEAN: 0.5 - 0.4 * sign(mean_p delta_pg)                     (|delta|-weighted direction)
Combined across a surface's allowed banks with equal weights (nanmean), and also with CLEAN4's four banks only. Scored on
gate signs against CLEAN4-LOLO (paired, target bootstrap). Analysis only; streaming mmap; < 4 GB.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
from inst636 import BANK_FORBIDDEN, Inst, clean4_lolo  # noqa: E402

WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
CLEAN4_BANKS = ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")


def bank_programmes(bank: str):
    d = ROOT / "data/donors" / bank
    if not (d / "delta.npy").exists():
        return None
    D = np.load(d / "delta.npy", mmap_mode="r")
    genes = [ln.strip() for ln in open(d / "genes.txt")]
    n_t, n_g = D.shape
    s = np.empty(n_t)
    for a in range(0, n_t, 2000):
        C = np.asarray(D[a:a + 2000], dtype=np.float32)
        s[a:a + 2000] = np.sqrt((C.astype(np.float64) ** 2).mean(1))
    w1, w2 = s / s.sum(), s ** 2 / (s ** 2).sum()
    q10, q25 = np.quantile(s, 0.90), np.quantile(s, 0.75)
    acc = {k: np.zeros(n_g) for k in ("U", "W1", "W2", "TOP10", "TOP25", "SUM")}
    n10, n25 = int((s >= q10).sum()), int((s >= q25).sum())
    for a in range(0, n_t, 2000):
        C = np.asarray(D[a:a + 2000], dtype=np.float32)
        neg = (C < 0).astype(np.float64)
        sl = slice(a, a + C.shape[0])
        acc["U"] += neg.sum(0)
        acc["W1"] += w1[sl] @ neg
        acc["W2"] += w2[sl] @ neg
        acc["TOP10"] += neg[s[sl] >= q10].sum(0)
        acc["TOP25"] += neg[s[sl] >= q25].sum(0)
        acc["SUM"] += C.astype(np.float64).sum(0)
    out = {"U": acc["U"] / n_t, "W1": acc["W1"], "W2": acc["W2"], "TOP10": acc["TOP10"] / n10, "TOP25": acc["TOP25"] / n25,
           "MEAN": 0.5 - 0.4 * np.sign(acc["SUM"])}
    print(f"bank {bank}: {n_t} perts, s median {np.median(s):.4f}, top10 {n10}", flush=True)
    return genes, out


def main():
    I = Inst()
    B = {}
    for b in BANK_FORBIDDEN:
        r = bank_programmes(b)
        if r is not None:
            B[b] = r
    C4 = {s: clean4_lolo(I, s)[0] for s in I.S}
    rows = []
    rows += I.score("CLEAN4 LOLO (reference)", C4)
    for form in ("U", "W1", "W2", "TOP10", "TOP25", "MEAN"):
        for pool, label in ((None, "allowed banks"), (CLEAN4_BANKS, "CLEAN4's banks")):
            P = {}
            for s in I.S:
                use = [b for b in B if s not in BANK_FORBIDDEN[b] and (pool is None or b in pool)]
                P[s] = np.nanmean(np.stack([I.gene_vec(s, B[b][0], B[b][1][form]) for b in use]), 0)
            rows += I.score(f"{form} over {label}", P, ref=C4)
    print(f"\n{'surface':8s} {'source':34s} {'tgt':6s} {'cov':>5s} {'acc':>6s} {'acc_w':>6s} {'up':>6s} {'dn':>6s} {'d_vs_C4 [95% CI]':>26s}")
    for r in rows:
        if r["targets"] != "all" and r["surface"] != "h1":
            continue
        d = (f"{r['d_vs_ref']:+.4f} [{r['d_ci'][0]:+.4f},{r['d_ci'][1]:+.4f}]" if "d_vs_ref" in r else "")
        print(f"{r['surface']:8s} {r['source'][:34]:34s} {r['targets']:6s} {r['coverage']:5.2f} {r['acc']:6.4f} "
              f"{r['acc_w']:6.4f} {r['acc_up']:6.4f} {r['acc_dn']:6.4f} {d:>26s}")
    (WORK / "resp638.json").write_text(json.dumps(rows, indent=1, default=float))


if __name__ == "__main__":
    sys.exit(main())
