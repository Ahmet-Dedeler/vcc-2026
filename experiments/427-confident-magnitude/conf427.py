"""427 analytic: large magnitudes where the donors agree. For each h1-test target p and gene g, the non-H1 donor tables
(K562-GWPS, K562-ess, RPE1, Nadig Jurkat, Orion HCT116, Orion HEK293T; pseudobulk log1p-CP50k deltas) give a vote.
Where >= n donors cover (p, g), all agree in sign and min |delta| >= tau, predict lfc = sign * kappa * mean|delta| / ln 2;
elsewhere keep our delivered lfc (426's base = the shipped 410b construction). nmae per the scorer on h1-test's gates.
Also reports the sign accuracy of those confident cells on the gate (the quantity that decides whether magnitude pays)."""
from pathlib import Path

import numpy as np
import polars as pl

ROOT = Path(__file__).resolve().parents[2]
W323 = ROOT / "experiments/323-current-family-validation/work"
r = np.load(W323 / "reference.npz", allow_pickle=True)
genes = r["genes"].astype(str)[r["keep"].astype(bool)]; gi = {g: i for i, g in enumerate(genes)}
s = np.load(ROOT / "experiments/426-mde-magnitude/work/h1/stats-base_s0_bps.npz"); targets = s["targets"].astype(str)
L = np.nan_to_num(s["lfc"]); Q = s["q"]
DON = ["k562_gwps", "k562_essential", "rpe1", "nadig_jurkat", "orion_hct116", "orion_hek293t"]
D = np.full((len(DON), len(targets), len(genes)), np.nan, dtype=np.float32)
for di, d in enumerate(DON):
    base = ROOT / "data/donors" / d
    P = [l.strip() for l in open(base / "perts.txt")]; G = [l.strip() for l in open(base / "genes.txt")]
    pi = {p: i for i, p in enumerate(P)}; gcol = {g: i for i, g in enumerate(G)}
    M = np.load(base / "delta.npy", mmap_mode="r")
    cols = np.array([gcol.get(g, -1) for g in genes]); okc = cols >= 0
    for k, t in enumerate(targets):
        if t in pi:
            row = np.asarray(M[pi[t]], dtype=np.float32)
            D[di, k, okc] = row[cols[okc]]
    print(f"{d}: targets covered {sum(t in pi for t in targets)}/{len(targets)}, genes {okc.sum()}/{len(genes)}", flush=True)
ref = pl.read_parquet(W323 / "reference_de.parquet")
gates = []
for k, t in enumerate(targets):
    m = ref.filter(pl.col("target") == t)
    g = m["feature"].to_numpy().astype(str); q = m["p_adj"].to_numpy().astype(float); l = m["log2_fold_change"].to_numpy().astype(float)
    ok = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
    if ok.sum() >= 10:
        gates.append((k, np.array([gi[x] for x in g[ok]]), l[ok]))
def nmae(Lh):
    return float(np.mean([np.abs(Lh[k, idx] - l).sum() / np.abs(l).sum() for k, idx, l in gates]))
base = nmae(L)
ncov = np.isfinite(D).sum(0)
pos = (np.nan_to_num(D) > 0).sum(0); neg = (np.nan_to_num(D) < 0).sum(0)
agree_sign = np.where(pos >= neg, 1.0, -1.0)
allagree = (np.maximum(pos, neg) == ncov) & (ncov > 0)
absD = np.abs(D)
minabs = np.where(np.isfinite(D), absD, np.inf).min(0)
meanabs = np.nanmean(np.where(np.isfinite(D), absD, np.nan), 0)
K = np.concatenate([np.full(i.size, k) for k, i, _ in gates]); I = np.concatenate([i for _, i, _ in gates]); LL = np.concatenate([l for _, _, l in gates])
print(f"\nbase nmae {base:.5f}")
print(" n>=  tau   cells/target  gate-cells  sign-acc(gate)  med|l| gate  med mean|d|/ln2  | kappa 0.5   1.0")
for nmin in (2, 3, 4):
    for tau in (0.1, 0.2, 0.3, 0.5, 0.8):
        conf = (ncov >= nmin) & allagree & (minabs >= tau)
        cg = conf[K, I]
        if cg.sum() == 0:
            continue
        acc = np.mean(agree_sign[K, I][cg] == np.sign(LL[cg]))
        out = []
        for kap in (0.5, 1.0):
            Lh = np.where(conf, agree_sign * kap * np.nan_to_num(meanabs) / np.log(2), L)
            out.append(nmae(Lh) - base)
        print(f" {nmin:3d} {tau:5.2f}  {conf.sum() / len(targets):11.1f}  {cg.sum():10d}  {acc:14.3f}  {np.median(np.abs(LL[cg])):11.3f}"
              f"  {np.median(meanabs[K, I][cg] / np.log(2)):15.3f}  | {out[0]:+.5f} {out[1]:+.5f}")

print("\nonly RAISE magnitude where confident: lfc = s * max(|ours| if sign(ours)==s else 0, kappa*mean|d|/ln2)")
for nmin in (2, 3):
    for tau in (0.05, 0.1, 0.2):
        conf = (ncov >= nmin) & allagree & (minabs >= tau)
        cg = conf[K, I]
        ours_g = L[K, I][cg]
        print(f" n>={nmin} tau {tau}: gate cells {cg.sum()}, our |lfc| median {np.median(np.abs(ours_g)):.3f}, our sign acc {np.mean(np.sign(ours_g) == np.sign(LL[cg])):.3f}", end="")
        res = []
        for kap in (1.0, 1.5, 2.0, 3.0):
            tgt = kap * np.nan_to_num(meanabs) / np.log(2)
            keep_ours = np.where(np.sign(L) == agree_sign, np.abs(L), 0.0)
            Lh = np.where(conf, agree_sign * np.maximum(keep_ours, tgt), L)
            res.append(nmae(Lh) - base)
        print("  | k1 %+.5f k1.5 %+.5f k2 %+.5f k3 %+.5f" % tuple(res))

# ---- reach: the walk over the reference's DE pool in OUR order (called first by q, then p, then -|lfc|, then name)
P_ = s["p"]
def reach(order_first=None, signs=None):
    vals = []
    for k, t in enumerate(targets):
        m = ref.filter(pl.col("target") == t)
        g = m["feature"].to_numpy().astype(str); q = m["p_adj"].to_numpy().astype(float); l = m["log2_fold_change"].to_numpy().astype(float)
        ok = (q < .05) & (g != t)
        if ok.sum() == 0:
            continue
        idx = np.array([gi[x] for x in g[ok]]); lr = l[ok]
        adj = np.isfinite(lr) & (lr != 0)
        sg_ = np.sign(L[k, idx]) if signs is None else np.sign(signs[k, idx])
        called = Q[k, idx] < .05
        key = [(0 if called[j] else 1, Q[k, idx[j]], P_[k, idx[j]], -abs(L[k, idx[j]]), genes[idx[j]]) for j in range(idx.size)]
        if order_first is not None:
            key = [((-1, -order_first[k, idx[j]]) if order_first[k, idx[j]] > 0 else (0,)) + kk for j, kk in enumerate(key)]
        order = sorted(range(idx.size), key=lambda j: key[j])
        c = n = best = 0
        for j in order:
            if not adj[j]:
                continue
            n += 1; c += int(sg_[j] == np.sign(lr[j]))
            if c / n >= 0.9:
                best = n
        vals.append(best / ok.sum())
    return float(np.mean(vals)), len(vals)
r0, nt = reach()
print(f"\nreach (delivered order) {r0:.5f} over {nt} targets")
for nmin, tau in ((3, 0.1), (3, 0.2), (2, 0.2), (4, 0.1)):
    conf = (ncov >= nmin) & allagree & (minabs >= tau)
    strength = np.where(conf, np.nan_to_num(minabs) + 1e-6, 0.0)
    sig = np.where(conf, agree_sign, L)   # confident cells carry the donors' sign
    r1, _ = reach(order_first=strength, signs=sig)
    print(f"  confident first (n>={nmin}, tau {tau}): reach {r1:.5f}  delta {r1 - r0:+.5f}")
