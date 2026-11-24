#!/usr/bin/env python3
"""397d -- H1 known-answer row for the SHIPPED 399 prior construction: six-bank family (CLEAN4 + raw Orion), K=25, shrink 400,
family-marginal default; signs on the called block only (q80 protection not modelled here, as in 397)."""
import numpy as np, json
exec(open("~/Code/virtual-cell-challenge-2026/experiments/397-knn-tilt/knn397.py").read().split("# ---- kNN tilts")[0].replace('print("KNOWN ANSWER', 'pass  # print("KNOWN ANSWER'))
for b in ("orion_hct116", "orion_hek293t", "nadig_hepg2"):
    d = (REPO / "data/donors" / b) if b != "nadig_hepg2" else (REPO / "experiments/166-donor-acquisition/donors/nadig_hepg2")
    D = np.asarray(np.load(d / "delta.npy", mmap_mode="r"), dtype=np.float32); P = np.array([l.strip() for l in open(d / "perts.txt")]); G = np.array([l.strip() for l in open(d / "genes.txt")])
    pi = {p: i for i, p in enumerate(P)}; gi_b = {g: i for i, g in enumerate(G)}; extra[b] = (D, pi, np.array([gi_b.get(g, -1) for g in genes]))
FAM = ["k562_essential", "nadig_jurkat", "nadig_hepg2", "orion_hct116", "orion_hek293t"]
P_marg6, _ = votes(np.arange(len(GWP)), FAM); sg6 = signs_from(P_marg6, sg4)
has = np.array([t in gwp for t in targets]); K = 25
for m in (50.0, 400.0):
    SG = SG0.copy()
    for i in range(nT):
        if not has[i]: SG[i] = np.where(cb[i], sg6, SG0[i]); continue
        v = GWn[gwp[t := targets[i]]]; s_ = GWn @ v; s_[gwp[t]] = -np.inf; nb = np.argsort(-s_)[:K]; P, nn = votes(nb, FAM)
        Ps = np.where(np.isfinite(P) & np.isfinite(P_marg6), (nn * np.nan_to_num(P) + m * np.nan_to_num(P_marg6)) / (nn + m), np.where(np.isfinite(P), P, P_marg6))
        SG[i] = np.where(cb[i], signs_from(Ps, sg6), SG0[i])
    eval_arm(f"SHIPPED-like kNN25 six-bank shrink{int(m)}", SG, {})
eval_arm("six-bank marginal", np.where(cb, sg6[None, :], SG0), {})
