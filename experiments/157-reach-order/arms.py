"""157 — score every realised emission arm's ORDER, natively and power-matched.

Consumes the tables written by `realize_shipped.py` / `realize_h1test.py`.  For each arm it
reports the emitted-order reach (what the scorer reads), the reach of the best counterfactual
key from 154, the oracle, and the anatomy (mean/median k*, k*=0 count, head purity) at
several `n_real_target` values.

The comparison that matters: 154's keys were scored against the **phi=0** emitted table and
against `|pred_lfc|` as their baseline.  The shipped emitter is phi=4.  If the shipped
emitted order already beats 154's best key, the buildable-key channel is behind the
incumbent, not ahead of it.

Deterministic given the tables.  No cell-eval2, no build, no submission.
"""
from __future__ import annotations
import glob
import os
import sys

import numpy as np
from scipy import stats

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
from debench import Surface, overall_delta  # noqa: E402

# The donor that IS the evaluation surface must be excluded. `h1_2025` covers 85/85
# h1-test targets, so leaving it in makes every donor key an oracle (it read 0.404
# reach before this was caught). `orion_hct116` is likewise excluded on de_thinned.
DONORS = os.environ.get("DONORS", "k562_gwps,orion_hek293t,h1_2025").split(",")
NS = [None, 68, 132, 240, 400]


def best_key(genes, perts, pred_lfc):
    """154's winner: geometric donor magnitude x cross-donor sign concordance."""
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
    n_avail = mask.sum(0).astype(float)
    z = np.zeros((P, G))
    concord = np.divide(np.abs((np.sign(stack) * mask).sum(0)), n_avail, out=z.copy(),
                        where=n_avail > 0)
    absmean = np.divide((np.abs(stack) * mask).sum(0), n_avail, out=z.copy(), where=n_avail > 0)
    return np.sqrt(np.maximum(absmean, 0.0)) * concord


def anatomy(s: Surface, key=None):
    """mean/median k*, k*=0 count, and purity of the first 10 pool genes in our order."""
    ks, ns, head = [], [], []
    for p in range(len(s.perts)):
        pool = np.flatnonzero(s.real_sig[p])
        if pool.size == 0:
            continue
        pool = pool[s._order(p, pool, key, None)]
        pool = pool[s.adj[p, pool]]
        if pool.size == 0:
            ks.append(0); ns.append(int(s.real_sig[p].sum())); continue
        ok = np.sign(s.pred_lfc[p, pool]) == np.sign(s.real_lfc[p, pool])
        pur = np.cumsum(ok) / np.arange(1, ok.size + 1)
        d = np.flatnonzero(pur >= 0.9)
        ks.append(int(d[-1] + 1) if d.size else 0)
        ns.append(int(s.real_sig[p].sum()))
        head.append(ok[:10].mean())
    ks = np.array(ks)
    return ks.mean(), float(np.median(ks)), int((ks == 0).sum()), len(ks), float(np.mean(head))


def main() -> None:
    pat = sys.argv[1] if len(sys.argv) > 1 else f"{ROOT}/experiments/157-reach-order/de_shipped.*.npz"
    paths = sorted(glob.glob(pat))
    if not paths:
        sys.exit(f"no tables match {pat}")
    lab = {None: "native BH"}
    s0 = Surface(paths[0])
    K = best_key(s0.genes, s0.perts, s0.pred_lfc)

    print(f"{'arm / order':<34}{'n_pred':>8}" + "".join(f"{lab.get(N, f'N={N}'):>11}" for N in NS))
    print("-" * (42 + 11 * len(NS)))
    store = {}
    for path in paths:
        arm = os.path.basename(path).split(".")[-2]
        row_emit, row_key, row_lfc, row_orc = [], [], [], []
        for N in NS:
            s = Surface(path, n_real_target=N)
            orc = (np.sign(s.pred_lfc) == np.sign(s.real_lfc)).astype(float)
            row_emit.append(s.reach())
            row_key.append(s.reach(key=K))
            row_lfc.append(s.reach(key=np.abs(s.pred_lfc)))
            row_orc.append(s.reach(key=orc))
        npd = float(np.median(Surface(path).counts()[1]))
        store[arm] = (row_emit, row_key, row_lfc, row_orc, npd)
        print(f"{arm + '  [EMITTED order]':<34}{npd:>8.0f}"
              + "".join(f"{v:>11.5f}" for v in row_emit))
        print(f"{'   + |pred_lfc| key':<34}{'':>8}" + "".join(f"{v:>11.5f}" for v in row_lfc))
        print(f"{'   + 154 best key (geo x concord)':<34}{'':>8}"
              + "".join(f"{v:>11.5f}" for v in row_key))
        print(f"{'   + ORACLE':<34}{'':>8}" + "".join(f"{v:>11.5f}" for v in row_orc))

    print("\n--- anatomy of the EMITTED order (native BH) ---")
    print(f"{'arm':<26}{'mean k*':>9}{'med k*':>8}{'k*=0':>7}{'perts':>7}{'top-10 purity':>15}")
    for path in paths:
        arm = os.path.basename(path).split(".")[-2]
        s = Surface(path)
        m, md, z, n, h = anatomy(s)
        print(f"{arm:<26}{m:>9.2f}{md:>8.0f}{z:>7d}{n:>7d}{h:>15.4f}")

    print("\n--- inversion of the EMITTED order (native BH) ---")
    print(f"{'arm':<26}{'acc in call':>13}{'acc out':>10}{'pooled rho':>12}{'pool acc':>10}")
    for path in paths:
        arm = os.path.basename(path).split(".")[-2]
        s = Surface(path)
        m = s.real_sig & s.adj
        ok = (np.sign(s.pred_lfc) == np.sign(s.real_lfc))
        sig = s.pred_sig & m
        out = m & ~s.pred_sig
        q = np.where(np.isfinite(s.pred_padj), s.pred_padj, 1.0)[m]
        rho = stats.spearmanr(q, ok[m].astype(float)).statistic
        print(f"{arm:<26}{ok[sig].mean():>13.4f}{ok[out].mean():>10.4f}{rho:>+12.4f}"
              f"{ok[m].mean():>10.4f}")

    print("\n--- deployable delta at N=240: (key order) - (shipped emitted order) ---")
    ship = [a for a in store if "a15" in a and not a.endswith("_b")]
    if ship:
        e, k, l, o, _ = store[ship[0]]
        i = NS.index(240)
        for nm, v in (("|pred_lfc|", l[i]), ("154 best key", k[i]), ("ORACLE", o[i])):
            print(f"  {nm:<20}{v - e[i]:>+10.5f} raw reach   "
                  f"{overall_delta(d_reach=v - e[i]):>+10.5f} overall")


if __name__ == "__main__":
    main()
