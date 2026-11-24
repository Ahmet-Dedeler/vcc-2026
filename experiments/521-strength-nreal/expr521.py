"""521 — the recipient line's own expression of the target as a null / negligible predictor (CRISPRi cannot knock down a
gene the line barely expresses). Build-time observable for any context (its own control cells).

H1: target CPM in H1's controls (323 reference b_ctrl) vs n_real; folds: target CPM in the held line's controls (the gate
files' ctrl_cpm; a target absent from the tested axis is < 5 CPM there). A/B/C: the distribution of the 300 panel targets'
control CPM per context (360's challenge_cpm.npz), to see whether A/B/C even HAS low-expression targets.
Also the strength + expression combination (mean rank).
    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/521-strength-nreal/expr521.py
"""
from __future__ import annotations

import json
import sys

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

import lib521 as L  # noqa: E402
from pred521 import BANKS, GATES, PANEL, REF, auc, cd4_feats, crank, evaluate, fold_surface, meanrank, strength  # noqa: E402


def main():
    out = {}
    ref = np.load(REF)
    g_full = ref["genes"].astype(str)
    b = np.asarray(ref["b_ctrl"], float)
    cpm = 1e6 * np.expm1(b) / np.expm1(b).sum()
    pos = {g: i for i, g in enumerate(g_full)}
    P = np.load(L.WORK / "pred521_h1.npz")
    T = [str(t) for t in P["targets"]]
    n = P["nreal"]
    e = np.array([cpm[pos[t]] if t in pos else 0.0 for t in T])
    le = np.log1p(e)
    print("H1: target control CPM vs n_real")
    for i in np.argsort(e)[:10]:
        print(f"   {T[i]:9s} cpm {e[i]:8.1f}  n_real {n[i]}")
    print(f"   null targets' CPM: " + ", ".join(f"{T[i]} {e[i]:.1f}" for i in np.flatnonzero(n == 0)))
    rows = [evaluate(le, n, "H1 log target CPM"), evaluate(P["S4__S3_p_vip_thr4"], n, "H1 S4 (reference)"),
            evaluate(meanrank([le, P["S4__S3_p_vip_thr4"]]), n, "H1 S4 + target CPM (mean rank)")]
    out["h1"] = rows
    for fold in ("k562", "hct116", "hek293t"):
        z = np.load(GATES[fold])
        Tf, nf = fold_surface(fold)
        gg = z["genes"].astype(str)
        cc = np.asarray(z["ctrl_cpm"], float)
        gp = {g: i for i, g in enumerate(gg)}
        ef = np.array([cc[gp[t]] if t in gp else 0.0 for t in Tf])
        print(f"\nfold {fold}: targets absent from the tested axis (< 5 CPM): {int(np.sum(ef == 0))}; their n_real "
              f"{sorted(nf[ef == 0].tolist())[:12]}")
        G = {bk: strength(BANKS[bk], Tf, 4.0) for bk in ("k562", "hct116", "hek293t", "vip", "kolf2")}
        G["cd4rob"], G["cd4nde"] = cd4_feats(Tf)
        prod = [bk for bk in ("k562", "hct116", "hek293t", "vip") if bk != fold and not (fold == "k562" and bk == "vip")]
        sall = meanrank([G[bk] for bk in prod + ["kolf2", "cd4rob", "cd4nde"]])
        fr = [evaluate(np.log1p(ef), nf, f"{fold} log target CPM"), evaluate(sall, nf, f"{fold} strength (all sources)"),
              evaluate(meanrank([np.log1p(ef), sall]), nf, f"{fold} strength + target CPM")]
        out[fold] = fr
    c = np.load(L.ROOT / "experiments/360-depmap-three-lines/work/challenge_cpm.npz", allow_pickle=True)
    cg = {g: i for i, g in enumerate(c["genes"].astype(str))}
    print("\nA/B/C: the 300 panel targets' own control CPM per context (quantiles 0 / 5 / 10 / 25 / 50 %) and counts < 5 / < 10")
    for ctx in ("A", "B", "C"):
        v = np.array([c[f"cpm_{ctx}"][cg[t]] if t in cg else np.nan for t in PANEL])
        q = np.nanquantile(v, [0, 0.05, 0.1, 0.25, 0.5])
        out[f"abc_{ctx}"] = dict(q=q.tolist(), lt5=int(np.sum(v < 5)), lt10=int(np.sum(v < 10)))
        print(f"   {ctx}: " + " ".join(f"{x:.1f}" for x in q) + f"   < 5 CPM: {int(np.sum(v < 5))}, < 10: {int(np.sum(v < 10))}")
    print("   H1's 85 targets for comparison: " + " ".join(f"{x:.1f}" for x in np.quantile(e, [0, 0.05, 0.1, 0.25, 0.5])))
    (L.WORK / "expr521.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
