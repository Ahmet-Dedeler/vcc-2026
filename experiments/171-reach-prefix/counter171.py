"""171 — counterfactual reach: what is the CEILING of moving the pole position?

Every arm here is an ORDERING or CALL-SET change applied to the SAME cached DE table, so the
emitted pseudobulk is literally the same array for every row -- `pds`, `mse` and `nmae` are
bit-identical by construction, not by argument.  This is the ceiling calculation: it says what
reach WOULD be if an emission law could realise the order, before any realisation loss.  A
counterfactual that does not pay here cannot pay after the emitter takes its cut (095: a
+0.00583 counterfactual key realised at -0.00021).

Also reported at proportional power `f` (157 section 3b), because h1-test's native n_real is
1-5,484 with 94% of reach coming from perturbations with n_real <= 13, and the challenge sits
near n_real ~ 240.  Uniform `n_real_target` is the wrong instrument (rule 9 / 157 section 3a).

KNOWN-ANSWER GATES:
  * `incumbent` must read 0.18283 at f = 1 (157 `power_h1test.log`, re-derived by 164).
  * `incumbent` on `armC_a8_flat` must read 0.20291.
  * the f = 1 pool must equal the native BH set (157 section 3b's own gate).
  * ORACLE must sit near `sign_accuracy / 0.9` -- 157 section 3b's structural identity; if it
    does not, the walk is wrong.
"""
from __future__ import annotations

import os
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
sys.path.insert(0, f"{ROOT}/experiments/171-reach-prefix")
from debench import Surface, overall_delta  # noqa: E402
from diag171 import Walk, kstar, spec_order  # noqa: E402

D164 = f"{ROOT}/experiments/164-callset-artefact"
ALPHA = 0.05
FACTORS = [1.0, 2.0, 3.5]


def pool_at_factor(W, f):
    """The reference's top round(f * n_real_p) genes by its own p-value (157 section 3b)."""
    if f == 1.0:
        return W.real_sig.copy()
    z = np.load(W._path, allow_pickle=True)
    rp = np.where(np.isfinite(z["real_p"]), z["real_p"], np.inf)
    sig = np.zeros(rp.shape, bool)
    for i in range(rp.shape[0]):
        k = int(round(f * int(W.real_sig[i].sum())))
        k = max(0, min(k, int(np.isfinite(z["real_p"][i]).sum())))
        if k:
            sig[i, np.argsort(rp[i], kind="stable")[:k]] = True
    sig &= np.isfinite(z["real_p"])
    for p, j in enumerate(W.self_idx):
        if j >= 0:
            sig[p, j] = False
    return sig


# --------------------------------------------------------------------------------------
# orderings.  Each returns the permutation of `pool`; higher priority first.
# --------------------------------------------------------------------------------------
def make_orders(W):
    genes = W.genes

    def incumbent(p, pool, sig):
        return spec_order(pool, genes, sig, W.pred_padj[p], W.pred_p[p], W.pred_lfc[p])

    def keepfirst(p, pool, sig):
        """Our significant calls first, but KEEP-block calls ahead of TAIL-block calls."""
        tail = W.spread[p, pool] > 0.5
        return np.lexsort((genes[pool], -np.abs(W.pred_lfc[p, pool]), W.pred_p[p, pool],
                           W.pred_padj[p, pool], tail, ~sig[pool]))

    def keeponly(p, pool, sig):
        """Only KEEP-block calls take pole position; TAIL-block calls join the non-called."""
        s = sig[pool] & ~(W.spread[p, pool] > 0.5)
        return np.lexsort((genes[pool], -np.abs(W.pred_lfc[p, pool]), W.pred_p[p, pool],
                           W.pred_padj[p, pool], ~s))

    def nosigpole(p, pool, sig):
        """Drop the significant-first rule entirely; order by p_adj alone."""
        return np.lexsort((genes[pool], -np.abs(W.pred_lfc[p, pool]), W.pred_p[p, pool],
                           W.pred_padj[p, pool]))

    def lfcfirst(p, pool, sig):
        """Significant calls first, then |predicted lfc| descending (ignore p_adj)."""
        return np.lexsort((genes[pool], W.pred_padj[p, pool],
                           -np.abs(W.pred_lfc[p, pool]), ~sig[pool]))

    def lfconly(p, pool, sig):
        """|predicted lfc| descending over the whole pool -- no significance pole at all."""
        return np.lexsort((genes[pool], -np.abs(W.pred_lfc[p, pool])))

    def keep_lfc(p, pool, sig):
        """KEEP-block calls first, ordered by |lfc|; then tail calls; then the rest."""
        tail = W.spread[p, pool] > 0.5
        return np.lexsort((genes[pool], W.pred_padj[p, pool],
                           -np.abs(W.pred_lfc[p, pool]), tail, ~sig[pool]))

    def oracle(p, pool, sig):
        ok = np.sign(W.pred_lfc[p, pool]) == np.sign(W.real_lfc[p, pool])
        return np.lexsort((genes[pool], ~ok))

    def antioracle(p, pool, sig):
        ok = np.sign(W.pred_lfc[p, pool]) == np.sign(W.real_lfc[p, pool])
        return np.lexsort((genes[pool], ok))

    return dict(incumbent=incumbent, keepfirst=keepfirst, keeponly=keeponly,
                nosigpole=nosigpole, lfcfirst=lfcfirst, lfconly=lfconly,
                keep_lfc=keep_lfc, ORACLE=oracle, ANTIORACLE=antioracle)


def reach_of(W, real_sig, order_fn, call=None):
    out, ks, heads = [], [], []
    for p in range(W.perts.size):
        pool = np.flatnonzero(real_sig[p])
        if pool.size == 0:
            continue
        sig = W.pred_sig[p] if call is None else call[p]
        o = order_fn(p, pool, sig)
        walked = pool[o]
        walked = walked[W.adj[p, walked]]
        if walked.size == 0:
            out.append(0.0)
            ks.append(0)
            continue
        ok = np.sign(W.pred_lfc[p, walked]) == np.sign(W.real_lfc[p, walked])
        k, _ = kstar(ok)
        out.append(k / pool.size)
        ks.append(k)
        heads.append(ok[:1].mean())
    return float(np.mean(out)), np.array(ks), float(np.mean(heads)), np.array(out)


def main():
    arm = os.environ.get("ARM", "a15")
    path = f"{D164}/de164.{arm}.npz"
    W = Walk(path)
    W._path = path
    S = Surface(path)
    orders = make_orders(W)

    print(f"=== 171 counterfactual reach ceiling — {arm}, h1-test")
    print()
    base = reach_of(W, W.real_sig, orders["incumbent"])[0]
    expect = {"a15": 0.18283, "armC_a8_flat": 0.20291}.get(arm)
    print(f"GATE incumbent f=1  {base:.5f}   debench {S.reach():.5f}   "
          f"expected {expect}   {'PASS' if expect and abs(base-expect) < 5e-5 else '?'}")
    acc = S.sign_accuracy()
    print(f"GATE oracle ~ accuracy/0.9 : accuracy {acc:.4f} -> {acc/0.9:.4f} ; "
          f"ORACLE reads {reach_of(W, W.real_sig, orders['ORACLE'])[0]:.4f}")
    print()

    # ---------------------------------------------------------------- ordering ceilings
    print("ORDERING COUNTERFACTUALS (call set FIXED at the shipped one; pseudobulk identical)")
    print(f"  {'ordering':<14} " + " ".join(f"{'f='+str(f):>9}" for f in FACTORS)
          + f" {'k* med':>7} {'P(1)':>7} {'d overall @f=1':>15}")
    pools = {f: pool_at_factor(W, f) for f in FACTORS}
    res = {}
    for name, fn in orders.items():
        vals, ks, p1 = [], None, None
        for f in FACTORS:
            r, k, h, _ = reach_of(W, pools[f], fn)
            vals.append(r)
            if f == 1.0:
                ks, p1 = k, h
        res[name] = vals
        d = overall_delta(d_reach=vals[0] - base)
        print(f"  {name:<14} " + " ".join(f"{v:9.5f}" for v in vals)
              + f" {np.median(ks):7.0f} {p1:7.4f} {d:+15.5f}")
    print()

    # ---------------------------------------------------------------- call-set arms
    print("CALL-SET COUNTERFACTUALS (which genes get the significant pole position)")
    print("  These also move fid/jac, so all four DE members are shown.  pds/nmae are")
    print("  functions of pred_lfc alone and are IDENTICAL across every row by construction.")
    calls = {
        "shipped": W.pred_sig,
        "keep only": W.pred_sig & (W.spread <= 0.5),
        "tail only": W.pred_sig & (W.spread > 0.5),
        "call nothing": np.zeros_like(W.pred_sig),
        "call everything": np.isfinite(W.pred_padj),
    }
    print(f"  {'call set':<16} {'n_pred med':>11} {'reach':>9} {'fid':>8} {'jac':>8} "
          f"{'nmae':>8} {'pds':>8} {'d overall':>10}")
    for name, c in calls.items():
        r = reach_of(W, W.real_sig, orders["incumbent"], call=c)[0]
        f_ = S.fid(call=c)
        j_ = S.jac(call=c)
        npd = np.median([int((c[p] & W.adj[p]).sum()) for p in range(W.perts.size)])
        d = overall_delta(d_reach=r - base, d_fid=f_ - S.fid(), d_jac=j_ - S.jac())
        print(f"  {name:<16} {npd:11.0f} {r:9.5f} {f_:8.5f} {j_:8.5f} "
              f"{S.nmae():8.5f} {S.pds():8.5f} {d:+10.5f}")
    print()

    # ---------------------------------------------------------------- tercile detail
    print("WHERE THE ORDERING GAIN LANDS, BY n_real TERCILE (f = 1)")
    nr = np.array([int(W.real_sig[p].sum()) for p in range(W.perts.size)])
    nr = nr[nr > 0]
    q = np.quantile(nr, [1 / 3, 2 / 3])
    print(f"  {'ordering':<14} {'small':>9} {'mid':>9} {'large':>9}   "
          f"(n_real <= {q[0]:.0f} / <= {q[1]:.0f} / >)")
    for name in ("incumbent", "keepfirst", "keeponly", "nosigpole", "lfconly", "ORACLE"):
        _, _, _, per = reach_of(W, W.real_sig, orders[name])
        m = [per[nr <= q[0]], per[(nr > q[0]) & (nr <= q[1])], per[nr > q[1]]]
        print(f"  {name:<14} " + " ".join(f"{v.mean():9.5f}" for v in m))
    print()

    # ---------------------------------------------------------------- head purity
    print("PURITY AT THE HEAD, BY ORDERING (mean over perturbations of P(k))")
    print(f"  {'ordering':<14} " + " ".join(f"{'P('+str(k)+')':>8}" for k in (1, 2, 3, 5, 10)))
    for name in ("incumbent", "keepfirst", "keeponly", "nosigpole", "lfcfirst", "lfconly",
                 "ORACLE", "ANTIORACLE"):
        fn = orders[name]
        cols = []
        for k in (1, 2, 3, 5, 10):
            v = []
            for p in range(W.perts.size):
                pool = np.flatnonzero(W.real_sig[p])
                if pool.size == 0:
                    continue
                walked = pool[fn(p, pool, W.pred_sig[p])]
                walked = walked[W.adj[p, walked]]
                if walked.size < k:
                    continue
                ok = np.sign(W.pred_lfc[p, walked[:k]]) == np.sign(W.real_lfc[p, walked[:k]])
                v.append(ok.mean())
            cols.append(np.mean(v))
        print(f"  {name:<14} " + " ".join(f"{c:8.4f}" for c in cols))


if __name__ == "__main__":
    main()
