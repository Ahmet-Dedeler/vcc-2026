#!/usr/bin/env python3
"""364 oracle endpoint: flip ALL wrong-sign gate pairs (p=1.0 pooled, deterministic)
and recompute nmae/fid/reach at the amplitude grid; plus mean-over-targets achieved
p for the 5 grid points (same seeds as breakeven364.py). Exact arithmetic, saved arrays only."""
import json, os, time
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
BASE = "~/Code/virtual-cell-challenge-2026/experiments/323-current-family-validation/work"
ALPHA = 0.05
t0 = time.monotonic()

s = np.load(os.path.join(BASE, "stats-cap15.npz"))
targets = s["targets"].astype(str)
pp, qp, lp0 = s["p"].astype(float), s["q"].astype(float), s["lfc"].astype(float)
s.close()
r = np.load(os.path.join(BASE, "reference.npz"), allow_pickle=False)
genes_full = r["genes"].astype(str)
keep = r["keep"].astype(bool)
r_q, r_lfc = r["q"].astype(float), r["lfc"].astype(float)
r.close()
ugenes = genes_full[keep]
n, g = len(targets), len(ugenes)
own = np.equal.outer(targets, ugenes)
adjud = np.isfinite(r_lfc) & (r_lfc != 0)
ref_sig = (r_q < ALPHA) & adjud
gate = ref_sig & ~own
elig = gate.sum(1) >= 10
w = np.where(gate & elig[:, None], np.abs(r_lfc), 0.0)
Dtot = w.sum()
correct0 = (np.sign(lp0) == np.sign(r_lfc)) & (lp0 != 0)
abs_lp0 = np.abs(lp0)

def fid_per_target(qp_, lp_):
    ps = np.isfinite(qp_) & (qp_ < ALPHA)
    scored = ps & adjud & ~own
    n_pred = scored.sum(1).astype(float)
    n_real = ((ref_sig & ~own)).sum(1).astype(float)
    match = scored & np.isfinite(lp_) & (lp_ != 0) & (np.sign(lp_) == np.sign(r_lfc))
    k = match.sum(1).astype(float)
    den = np.maximum(n_pred, n_real)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, k / den, np.nan)

def nmae_per_target(lp_):
    gt = ref_sig & ~own
    cnt = gt.sum(1)
    num = np.where(gt, np.abs(lp_ - r_lfc), 0.0).sum(1)
    den = np.where(gt, np.abs(r_lfc), 0.0).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where((cnt >= 10) & np.isfinite(den) & (den != 0), num / den, np.nan)

def reach_per_target(lp_, qp_, pp_):
    pool = ref_sig & ~own
    n_conf = pool.sum(1).astype(float)
    R = np.full(n, np.nan)
    qm = np.where(np.isfinite(qp_), qp_, np.inf)
    pm = np.where(np.isfinite(pp_), pp_, np.inf)
    ala = np.abs(lp_)
    sig = np.isfinite(qp_) & (qp_ < ALPHA)
    feat = ugenes
    for i in range(n):
        if n_conf[i] == 0:
            continue
        idx = np.lexsort((feat, -np.where(np.isfinite(ala[i]), ala[i], -np.inf),
                          pm[i], qm[i], ~sig[i]))
        m = ((np.isfinite(lp_[i]) & (lp_[i] != 0)
              & (np.sign(lp_[i]) == np.sign(r_lfc[i]))) & pool[i])[idx]
        a = (adjud[i] & pool[i])[idx]
        ma = m[a]
        if ma.size == 0:
            R[i] = 0.0
            continue
        pur = np.cumsum(ma) / np.arange(1, ma.size + 1)
        ok = np.flatnonzero(pur >= 0.9)
        R[i] = (int(ok.max() + 1) / n_conf[i]) if ok.size else 0.0
    return R

# oracle: every gate pair correct sign (eligible or not — fid/reach read all)
s_or = np.sign(lp0).copy()
s_or[gate] = np.sign(r_lfc[gate])
# note: delivered lh==0 pairs (if any) now commit; count them
n_zero = int(((lp0 == 0) & gate).sum())
out = {"n_delivered_zero_on_gate": n_zero}
for a in (1.0, 1.2, 1.45, 2.0):
    L = lp0.copy(); L[gate] = a * s_or[gate] * abs_lp0[gate]
    N = nmae_per_target(L); F = fid_per_target(qp, L); R = reach_per_target(L, qp, pp)
    out[str(a)] = {"nmae": float(np.nanmean(N[elig])), "fid": float(np.nanmean(F)),
                   "reach": float(np.nanmean(R))}
    print(f"oracle a={a}: nmae={out[str(a)]['nmae']:.5f} fid={out[str(a)]['fid']:.6f} reach={out[str(a)]['reach']:.6f}", flush=True)

# mean-over-targets achieved p per grid point (same seeds 1000+d as breakeven364.py)
wrong_pairs = [(i, j) for i in np.flatnonzero(elig) for j in np.flatnonzero(gate[i] & ~correct0[i])]
wrong_mass = np.array([abs(r_lfc[i, j]) for (i, j) in wrong_pairs])
need_base = (w * correct0).sum()
for p in (0.60, 0.62, 0.65, 0.70, 0.75):
    need = p * Dtot - need_base
    mot = []
    for d in range(20):
        rng = np.random.default_rng(1000 + d)
        perm = rng.permutation(len(wrong_pairs))
        s_new = np.sign(lp0).copy()
        if need > 0:
            cum = np.cumsum(wrong_mass[perm])
            k = min(int(np.searchsorted(cum, need, side="left")) + 1, len(perm))
            for t in perm[:k]:
                i, j = wrong_pairs[t]
                s_new[i, j] = np.sign(r_lfc[i, j])
        pw = np.array([float((w[i] * (s_new[i] == np.sign(r_lfc[i]))).sum() / w[i].sum())
                       for i in range(n) if elig[i]])
        mot.append(float(pw.mean()))
    print(f"p_pooled={p}: mean-over-targets achieved {np.mean(mot):.5f}±{np.std(mot, ddof=1):.5f}", flush=True)
    out[f"mot_{p}"] = {"mean": float(np.mean(mot)), "sd": float(np.std(mot, ddof=1))}

json.dump(out, open(os.path.join(WORK, "oracle364.json"), "w"), indent=1)
print(f"done in {time.monotonic()-t0:.1f}s", flush=True)
