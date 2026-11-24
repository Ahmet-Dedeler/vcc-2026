#!/usr/bin/env python3
"""517 — 462's control-halves rule for the RECIPIENT-SIDE sources of §2: the identity axis (-I: H1 controls' log mean CPM
minus A/B/C's) is computed from h1-test's OWN controls, the same cells the truth is taken against, so its H1 accuracy could
be control noise. Split h1-test's controls in halves (halves517's md5 rule): I_X from half X; truth T_Y against half Y.
  shared = -I_A vs T_A (and B/B);  disjoint = -I_A vs T_B (and B/A).   Also z(CLEAN4) + 2 z(-I_X), 517's c4id2.
KA: truth against all controls == the reference lfc.  usage: uv run python experiments/517-flex-programme/halves517b.py
"""
from __future__ import annotations

import csv
import glob
import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from stream517 import CSRStream  # noqa: E402


def main():
    z = np.load(glob.glob(str(ROOT / "experiments/323-*/work/reference.npz"))[0], allow_pickle=False)
    genes, keep, targets = z["genes"].astype(str), np.asarray(z["keep"], bool), list(z["targets"].astype(str))
    q_ref, lfc_ref, bctrl = np.asarray(z["q"], float), np.asarray(z["lfc"], float), np.asarray(z["b_ctrl"], float)
    kg = genes[keep]
    cols = np.full(genes.size, -1, np.int64); cols[np.flatnonzero(keep)] = np.arange(int(keep.sum()))
    it = CSRStream(ROOT / "data/harness/h1-test.h5ad")
    f = it.f
    tg = f["obs/target_gene"]; tv = tg["categories"][:].astype(str)[tg["codes"][:]]
    ix = f["obs/_index"]; ids = (ix[:] if not isinstance(ix, h5py.Group) else ix["values"][:]).astype(str)
    is_c = tv == "non-targeting"
    half = np.array([int(hashlib.md5(s.encode()).hexdigest()[:8], 16) & 1 for s in ids])
    tid = {t: i for i, t in enumerate(targets)}
    T = len(targets)
    lab = np.array([tid[t] if t in tid else (T + half[k] if is_c[k] else -1) for k, t in enumerate(tv)])
    lib = it.libraries()
    S = it.group_sums(lab, T + 2, cols=cols, cpm=True, lib=lib)
    it.close()
    n = np.bincount(lab[lab >= 0], minlength=T + 2).astype(float)
    Mt = S[:T] / n[:T, None]; mA, mB = S[T] / n[T], S[T + 1] / n[T + 1]
    mF = (S[T] + S[T + 1]) / (n[T] + n[T + 1])
    lg = lambda m, c: np.log2((m + 1e-9) / (c + 1e-9))   # noqa: E731
    TF, TA, TB = lg(Mt, mF), lg(Mt, mA), lg(Mt, mB)
    fin = np.isfinite(lfc_ref) & np.isfinite(TF)
    ka = float(np.abs(TF[fin] - lfc_ref[fin]).max())
    print(f"KA truth(all controls) vs reference lfc: max |d| {ka:.2e}; halves {int(n[T])} / {int(n[T + 1])}", flush=True)
    assert ka < 1e-6
    offi = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
    opos = {g: i for i, g in enumerate(offi)}
    abc = []
    for c in "ABC":
        mc = np.load(ROOT / f"experiments/491-own-control-programme/work/ctrl_meancpm_{c}.npy")
        abc.append(np.array([mc[opos[g]] if g in opos else np.nan for g in kg]))
    with np.errstate(invalid="ignore"):
        abc = np.nanmean(np.vstack([np.log(np.maximum(v, 1e-3)) for v in abc]), 0)
    negI = {X: -(np.log(np.maximum(m, 1e-3)) - abc) for X, m in (("A", mA), ("B", mB), ("F", mF))}
    src = np.load(HERE / "work/sources517.npz", allow_pickle=False)
    assert list(src["genes_keep"]) == list(kg)
    zc4 = np.where(np.isfinite(src["P_c4"]), (0.5 - src["P_c4"]) / np.nanstd(src["P_c4"]), 0.0)
    e_ = np.expm1(bctrl[keep]); cpm = 1e6 * e_ / e_.sum()
    order = np.argsort(-cpm)
    block = np.zeros(kg.size, bool); block[[i for i in order if kg[i] not in set(targets)][:2200]] = True
    own = np.array([[g == t for g in kg] for t in targets])
    gate = (q_ref < 0.05) & np.isfinite(lfc_ref) & (lfc_ref != 0) & ~own

    def acc(score, Tm, region, gated):
        s = np.where(np.isfinite(score) & (score != 0), np.sign(score), np.nan)
        hit = nn = 0
        for t in range(T):
            m = region & np.isfinite(s) & np.isfinite(Tm[t]) & (Tm[t] != 0) & ~own[t]
            if gated:
                m &= gate[t]
            hit += int((np.sign(Tm[t][m]) == s[m]).sum()); nn += int(m.sum())
        return hit / max(nn, 1)
    res = {"ka": ka}
    for sname, mk in (("-I", lambda X: negI[X]), ("CLEAN4 + 2(-I)", lambda X: zc4 + 2 * negI[X] / np.nanstd(negI[X]))):
        for rname, (region, gated) in {"fid-like block": (block, False), "gate block": (block, True),
                                       "gate nonblock": (~block, True)}.items():
            sh = np.mean([acc(mk("A"), TA, region, gated), acc(mk("B"), TB, region, gated)])
            dj = np.mean([acc(mk("A"), TB, region, gated), acc(mk("B"), TA, region, gated)])
            full = acc(mk("F"), TF, region, gated)
            c4 = np.mean([acc(zc4, TA, region, gated), acc(zc4, TB, region, gated)])
            res[f"{sname} | {rname}"] = {"full": full, "shared": sh, "disjoint": dj, "inflation": sh - dj, "clean4_halves": c4}
            print(f"{sname:16s} {rname:15s} full {full:.4f} | shared {sh:.4f} disjoint {dj:.4f} (inflation {sh - dj:+.4f}) "
                  f"| CLEAN4 on halves {c4:.4f}", flush=True)
    (HERE / "work/halves517b.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
