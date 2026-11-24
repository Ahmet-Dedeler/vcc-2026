"""157 — a power-matching instrument that `reach` can actually be read through.

`debench.Surface(n_real_target=N)` replaces every perturbation's BH set with the reference's
top-N by p-value.  That is the right emulation for `fid` and `jac`, whose artefact is the
SIZE of `n_real` against `n_pred`.  It is the wrong one for `reach`, and h1-test shows why:
h1-test's native median `n_real` is 222, so `N=240` should be nearly a no-op, and instead
the emitted order's reach falls 0.183 -> 0.008.

The reason is that `reach = mean_p(k*_p / n_real_p)` and `n_real_p` is strongly heterogeneous.
153 §11 measured reach by `n_real` tercile at 0.150 / 0.070 / 0.052 — most of our reach comes
from the perturbations with SMALL reference DE sets.  Forcing every perturbation to the same
large N deletes that stratum by construction, so the collapse measures the flattening, not
the power.

This scales each perturbation's pool by a common FACTOR instead:
    pool_p = the reference's top-round(f * n_real_p) genes by its own p-value,
which raises power while preserving the per-perturbation heterogeneity the metric averages
over.  `f = 1` reproduces the native BH set to within p-value ties (a known-answer row, per
AGENTS.md 11.1 rule 10's corollary).

Deterministic.  No emitter here, no cell-eval2, no build, no submission.
"""
from __future__ import annotations
import glob
import os
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
from debench import Surface, overall_delta  # noqa: E402

DONORS = os.environ.get("DONORS", "k562_gwps,orion_hek293t,h1_2025").split(",")
FACTORS = [0.5, 1.0, 1.5, 2.0, 3.5, 5.0]
ALPHA = 0.05


class PowerSurface(Surface):
    """Surface whose reference pool is scaled by a factor of its own BH set size."""

    def __init__(self, path, factor=1.0):
        super().__init__(path)
        z = np.load(path, allow_pickle=True)
        rp = np.where(np.isfinite(z["real_p"]), z["real_p"], np.inf)
        base = (z["real_padj"] < ALPHA) & np.isfinite(z["real_padj"])
        sig = np.zeros(rp.shape, bool)
        for i in range(rp.shape[0]):
            k = int(round(factor * int(base[i].sum())))
            k = max(0, min(k, int(np.isfinite(z["real_p"][i]).sum())))
            if k:
                sig[i, np.argsort(rp[i], kind="stable")[:k]] = True
        self.real_sig = sig & np.isfinite(z["real_p"])
        for p, j in enumerate(self.self_idx):
            if j >= 0:
                self.real_sig[p, j] = False


def best_key(genes, perts, pred_lfc):
    P, G = pred_lfc.shape
    mats, haves = [], []
    for name in DONORS:
        root = f"{ROOT}/data/donors/{name}"
        dg = np.array(open(f"{root}/genes.txt").read().split())
        dp = np.array(open(f"{root}/perts.txt").read().split())
        delta = np.load(f"{root}/delta.npy", mmap_mode="r")
        gi = {g: i for i, g in enumerate(dg)}
        cols = np.array([gi.get(g, -1) for g in genes])
        ok = cols >= 0
        pi = {p: i for i, p in enumerate(dp)}
        M = np.zeros((P, G)); have = np.zeros(P, bool)
        for j, p in enumerate(perts):
            r = pi.get(p)
            if r is None:
                continue
            M[j, ok] = np.asarray(delta[r], float)[cols[ok]]
            have[j] = True
        mats.append(M); haves.append(have)
    stack = np.stack(mats)
    mask = np.stack([np.broadcast_to(h[:, None], (P, G)) for h in haves])
    na = mask.sum(0).astype(float)
    z = np.zeros((P, G))
    con = np.divide(np.abs((np.sign(stack) * mask).sum(0)), na, out=z.copy(), where=na > 0)
    am = np.divide((np.abs(stack) * mask).sum(0), na, out=z.copy(), where=na > 0)
    return np.sqrt(np.maximum(am, 0.0)) * con


def tercile_reach(s, key=None):
    nr = np.array([int(s.real_sig[p].sum()) for p in range(len(s.perts))])
    vals, sizes = [], []
    for p in range(len(s.perts)):
        pool = np.flatnonzero(s.real_sig[p])
        if pool.size == 0:
            continue
        o = pool[s._order(p, pool, key, None)]
        o = o[s.adj[p, o]]
        if o.size == 0:
            vals.append(0.0); sizes.append(nr[p]); continue
        ok = np.sign(s.pred_lfc[p, o]) == np.sign(s.real_lfc[p, o])
        pur = np.cumsum(ok) / np.arange(1, ok.size + 1)
        d = np.flatnonzero(pur >= 0.9)
        vals.append((int(d[-1]) + 1) / nr[p] if d.size else 0.0)
        sizes.append(nr[p])
    v, sz = np.array(vals), np.array(sizes)
    q = np.quantile(sz, [1 / 3, 2 / 3])
    return [float(v[sz <= q[0]].mean()), float(v[(sz > q[0]) & (sz <= q[1])].mean()),
            float(v[sz > q[1]].mean())]


def main() -> None:
    pat = sys.argv[1] if len(sys.argv) > 1 else f"{ROOT}/experiments/157-reach-order/de_shipped.*.npz"
    paths = sorted(glob.glob(pat))
    s0 = Surface(paths[0])
    K = best_key(s0.genes, s0.perts, s0.pred_lfc)
    print(f"donors used for the key: {DONORS}")
    print(f"{'arm / order':<36}" + "".join(f"{'f=' + str(f):>10}" for f in FACTORS))
    print("-" * (36 + 10 * len(FACTORS)))
    cap = {}
    for path in paths:
        arm = os.path.basename(path).replace("de_shipped.", "").replace("de_h1test.", "")[:-4]
        rows = {"EMITTED": [], "|pred_lfc|": [], "154 key": [], "ORACLE": []}
        med = []
        for f in FACTORS:
            s = PowerSurface(path, f)
            med.append(float(np.median([int(s.real_sig[p].sum()) for p in range(len(s.perts))])))
            orc = (np.sign(s.pred_lfc) == np.sign(s.real_lfc)).astype(float)
            rows["EMITTED"].append(s.reach())
            rows["|pred_lfc|"].append(s.reach(key=np.abs(s.pred_lfc)))
            rows["154 key"].append(s.reach(key=K))
            rows["ORACLE"].append(s.reach(key=orc))
        for nm, v in rows.items():
            print(f"{arm + '  ' + nm:<36}" + "".join(f"{x:>10.5f}" for x in v))
        print(f"{'    median n_real':<36}" + "".join(f"{x:>10.0f}" for x in med))
        cap[arm] = rows

    print(f"\n--- headroom capture, best buildable vs oracle, over the EMITTED order ---")
    print(f"{'arm':<30}" + "".join(f"{'f=' + str(f):>10}" for f in FACTORS))
    for arm, r in cap.items():
        c = [max(r["|pred_lfc|"][i], r["154 key"][i]) for i in range(len(FACTORS))]
        pct = [(c[i] - r["EMITTED"][i]) / max(r["ORACLE"][i] - r["EMITTED"][i], 1e-12) * 100
               for i in range(len(FACTORS))]
        print(f"{arm + ' capture %':<30}" + "".join(f"{x:>10.1f}" for x in pct))
        d = [overall_delta(d_reach=c[i] - r["EMITTED"][i]) for i in range(len(FACTORS))]
        print(f"{arm + ' d_overall':<30}" + "".join(f"{x:>+10.5f}" for x in d))

    print("\n--- reach by n_real tercile, EMITTED order, f=1 (small / mid / large) ---")
    for path in paths:
        arm = os.path.basename(path).replace("de_shipped.", "").replace("de_h1test.", "")[:-4]
        s = PowerSurface(path, 1.0)
        t = tercile_reach(s)
        tk = tercile_reach(s, K)
        print(f"{arm:<30} emitted {t[0]:.4f} / {t[1]:.4f} / {t[2]:.4f}"
              f"    154 key {tk[0]:.4f} / {tk[1]:.4f} / {tk[2]:.4f}")


if __name__ == "__main__":
    main()
