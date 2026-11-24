"""171 — decompose the reach walk on h1-test: where does P(k) die, and what kills it?

Reads 164's cached h1-test DE tables (`de164.*.npz`), which carry the reference DE table,
our emitted p / p_adj / lfc, AND the per-target emission spread vector -- the last is the
exact block label for "did our model give this gene a real predicted effect".

Nothing is emitted here.  Deterministic, seed-free, seconds.

KNOWN-ANSWER GATES (AGENTS.md 11.1 rule 10's corollary).  Two independent ones:
  * `debench.Surface(de164.a15.npz).reach()` must read 0.18283 (157 `power_h1test.log`,
    re-derived by 164 from a fresh emission).  armC must read 0.20291.
  * my own walk, written from BRIEF section 5 rather than from debench, must agree with
    debench's to 1e-12 on the production order.  If the two implementations agree only
    because they are the same code, the check is vacuous -- so `reach_walk` below is written
    against the spec text and reproduces the ordering chain independently.

`spread == tail_spread` is the block label.  a15 ranks the gated axis by
`|log2(pred_cpm / ctrl_cpm)|` and gives the top 8,000 `keep_spread = 0.25`; everything else
gets `tail_spread = 1.0`.  So "spread == 1.0" is exactly "our model put this gene outside its
top-8,000 by predicted effect magnitude" -- the operational form of zero/near-zero predicted
effect.  On h1-test's 10,779-gene scored axis that is 2,779 genes per target.
"""
from __future__ import annotations

import os
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
from debench import Surface  # noqa: E402

D164 = f"{ROOT}/experiments/164-callset-artefact"
ALPHA = 0.05
P0 = 0.9


# --------------------------------------------------------------------------------------
# the walk, written from docs/reference/vcc2026-metrics-BRIEF.md section 5
# --------------------------------------------------------------------------------------
def spec_order(pool, genes, sig, padj, p, lfc):
    """The submission's confidence order over `pool`.

    BRIEF section 5: "its significant calls first, then ascending predicted adjusted
    p-value, then ascending predicted p-value, then descending |lfc|, then gene name."
    np.lexsort takes the LAST key as primary.
    """
    return np.lexsort((genes[pool], -np.abs(lfc[pool]), p[pool], padj[pool], ~sig[pool]))


def kstar(ok):
    """k* = max{k : P(k) >= 0.9}, 0 if no prefix clears it.  "Deepest", not "first"."""
    if ok.size == 0:
        return 0, np.zeros(0)
    purity = np.cumsum(ok) / np.arange(1, ok.size + 1)
    deep = np.flatnonzero(purity >= P0)
    return (int(deep[-1] + 1) if deep.size else 0), purity


class Walk:
    """One surface, one arm: everything the reach walk touches, per perturbation."""

    def __init__(self, path):
        z = np.load(path, allow_pickle=True)
        self.genes = np.asarray([str(g) for g in z["genes"]])
        self.perts = np.asarray([str(p) for p in z["perts"]])
        self.real_padj, self.real_lfc = z["real_padj"], z["real_lfc"]
        self.pred_padj, self.pred_p, self.pred_lfc = z["pred_padj"], z["pred_p"], z["pred_lfc"]
        self.spread = z["spread"] if "spread" in z.files else None
        gpos = {g: i for i, g in enumerate(self.genes)}
        self.self_idx = np.array([gpos.get(p, -1) for p in self.perts])
        self.adj = np.isfinite(self.real_lfc) & (self.real_lfc != 0.0)
        self.real_sig = (self.real_padj < ALPHA) & np.isfinite(self.real_padj)
        self.pred_sig = (self.pred_padj < ALPHA) & np.isfinite(self.pred_padj)
        for i, j in enumerate(self.self_idx):
            if j >= 0:
                self.adj[i, j] = self.real_sig[i, j] = self.pred_sig[i, j] = False

    def rows(self, call=None, order_fn=None):
        """Per-perturbation anatomy of the walk.  `order_fn(p, pool)` overrides the order."""
        out = []
        for p in range(self.perts.size):
            pool = np.flatnonzero(self.real_sig[p])
            n_real = pool.size
            if n_real == 0:
                continue
            sig = self.pred_sig[p] if call is None else call[p]
            if order_fn is None:
                o = spec_order(pool, self.genes, sig, self.pred_padj[p],
                               self.pred_p[p], self.pred_lfc[p])
            else:
                o = order_fn(p, pool, sig)
            walked = pool[o]
            walked = walked[self.adj[p, walked]]
            ok = (np.sign(self.pred_lfc[p, walked]) == np.sign(self.real_lfc[p, walked]))
            k, purity = kstar(ok)
            out.append(dict(
                p=p, n_real=n_real, n_adj=walked.size, walked=walked, ok=ok,
                purity=purity, kstar=k, R=k / n_real,
                sig_walked=sig[walked],
                tail_walked=(self.spread[p, walked] > 0.5) if self.spread is not None else None,
                n_sig_in_pool=int(sig[pool].sum()),
                n_sig_in_walk=int(sig[walked].sum()),
            ))
        return out

    def reach(self, rows=None):
        rows = self.rows() if rows is None else rows
        return float(np.mean([r["R"] for r in rows]))


def pct(x):
    return "  n/a" if x is None else f"{100*x:5.1f}"


def main():
    arm = os.environ.get("ARM", "a15")
    W = Walk(f"{D164}/de164.{arm}.npz")
    S = Surface(f"{D164}/de164.{arm}.npz")

    rows = W.rows()
    mine, theirs = W.reach(rows), S.reach()
    print(f"=== 171 reach-walk decomposition — arm {arm}, surface h1-test "
          f"({W.perts.size} perturbations x {W.genes.size} gated genes)")
    print()
    print("KNOWN-ANSWER GATES")
    expect = {"a15": 0.18283, "armC_a8_flat": 0.20291}.get(arm)
    print(f"  debench Surface.reach()          {theirs:.5f}"
          + (f"   expected {expect:.5f}  "
             f"{'PASS' if expect is None or abs(theirs-expect) < 5e-5 else 'FAIL'}"
             if expect is not None else ""))
    print(f"  independent spec walk            {mine:.5f}   "
          f"|diff| {abs(mine-theirs):.3e}  "
          f"{'PASS' if abs(mine-theirs) < 1e-12 else 'FAIL'}")
    print()

    nr = np.array([r["n_real"] for r in rows])
    na = np.array([r["n_adj"] for r in rows])
    ks = np.array([r["kstar"] for r in rows])
    R = np.array([r["R"] for r in rows])
    nsig = np.array([r["n_sig_in_walk"] for r in rows])

    print("PER-PERTURBATION ANATOMY")
    print(f"  n_real                 median {np.median(nr):7.0f}   mean {nr.mean():8.1f}   "
          f"min {nr.min()}  max {nr.max()}")
    print(f"  adjudicable in pool    median {np.median(na):7.0f}   mean {na.mean():8.1f}")
    print(f"  of those, WE CALL SIG  median {np.median(nsig):7.0f}   mean {nsig.mean():8.1f}   "
          f"= {100*nsig.sum()/max(na.sum(),1):.1f}% of the walked pool")
    print(f"  k*                     median {np.median(ks):7.0f}   mean {ks.mean():8.1f}   "
          f"max {ks.max()}")
    print(f"  k* == 0                {100*np.mean(ks == 0):5.1f}% of perturbations")
    print(f"  R_p = k*/n_real        median {np.median(R):7.4f}   mean {R.mean():8.4f}")
    print()

    # ---- the head of the walk: who takes the pole positions -----------------------------
    print("HEAD COMPOSITION — what sits at positions 1..h of the walk, pooled over perts")
    print(f"  {'h':>4} {'% our sig calls':>16} {'% TAIL block':>14} {'% sig AND tail':>16} "
          f"{'mean P(h)':>10} {'% perts P(h)>=0.9':>18}")
    for h in (1, 2, 3, 5, 10, 20, 50, 100):
        s = t = st = 0
        n = 0
        pk, clear = [], []
        for r in rows:
            m = min(h, r["n_adj"])
            if m == 0:
                continue
            n += m
            s += int(r["sig_walked"][:m].sum())
            if r["tail_walked"] is not None:
                t += int(r["tail_walked"][:m].sum())
                st += int((r["sig_walked"][:m] & r["tail_walked"][:m]).sum())
            pk.append(r["purity"][m - 1])
            clear.append(r["purity"][m - 1] >= P0)
        print(f"  {h:>4} {100*s/n:15.1f}% {100*t/n:13.1f}% {100*st/n:15.1f}% "
              f"{np.mean(pk):10.4f} {100*np.mean(clear):17.1f}%")
    print()

    # ---- the prefix that sets k* --------------------------------------------------------
    print("THE PREFIX AHEAD OF k*  (positions 1..k*, i.e. the genes that earned the score)")
    tot_k = sum(r["kstar"] for r in rows)
    tail_k = sum(int(r["tail_walked"][:r["kstar"]].sum()) for r in rows if r["kstar"])
    sig_k = sum(int(r["sig_walked"][:r["kstar"]].sum()) for r in rows if r["kstar"])
    print(f"  genes ahead of k*, pooled            {tot_k}")
    print(f"  of those, ZERO-EFFECT (tail block)   {tail_k}  = {100*tail_k/max(tot_k,1):.1f}%")
    print(f"  of those, our significant calls      {sig_k}  = {100*sig_k/max(tot_k,1):.1f}%")
    print()
    # per-perturbation version (rule 8: the metric averages per perturbation)
    fr = [r["tail_walked"][:r["kstar"]].mean() for r in rows if r["kstar"]]
    print(f"  per-perturbation zero-effect fraction ahead of k*: "
          f"mean {np.mean(fr):.3f}  median {np.median(fr):.3f}  "
          f"({len(fr)}/{len(rows)} perturbations with k* > 0)")
    print()

    # ---- where the walk first goes wrong ------------------------------------------------
    print("FIRST SIGN ERROR IN THE WALK")
    first = []
    for r in rows:
        w = np.flatnonzero(~r["ok"])
        first.append(int(w[0] + 1) if w.size else r["n_adj"] + 1)
    first = np.array(first)
    print(f"  position of the first wrong sign: median {np.median(first):.0f}  "
          f"mean {first.mean():.1f}")
    for h in (1, 2, 3, 5, 10):
        print(f"    first error at position <= {h:<3d}  {100*np.mean(first <= h):5.1f}% "
              f"of perturbations")
    print()

    # ---- accuracy by block, restricted to the reach pool --------------------------------
    print("SIGN ACCURACY ON THE REACH POOL, BY BLOCK  (pooled and per-perturbation)")
    for lab, sel in (("our sig calls, KEEP block", lambda r: r["sig_walked"] & ~r["tail_walked"]),
                     ("our sig calls, TAIL block", lambda r: r["sig_walked"] & r["tail_walked"]),
                     ("not called, KEEP block", lambda r: ~r["sig_walked"] & ~r["tail_walked"]),
                     ("not called, TAIL block", lambda r: ~r["sig_walked"] & r["tail_walked"])):
        num = den = 0
        per = []
        for r in rows:
            m = sel(r)
            num += int(r["ok"][m].sum())
            den += int(m.sum())
            if m.sum():
                per.append(r["ok"][m].mean())
        print(f"  {lab:<26} n={den:8d}  pooled {num/max(den,1):.4f}   "
              f"per-pert {np.mean(per):.4f} (sd of mean {np.std(per)/np.sqrt(len(per)):.4f})")
    print()

    # ---- the n_real terciles (157: 98% of reach comes from the small tercile) -----------
    print("REACH BY n_real TERCILE")
    q = np.quantile(nr, [1 / 3, 2 / 3])
    for lab, m in (("small", nr <= q[0]), ("mid", (nr > q[0]) & (nr <= q[1])),
                   ("large", nr > q[1])):
        idx = np.flatnonzero(m)
        tf = [rows[i]["tail_walked"][:rows[i]["kstar"]].mean()
              for i in idx if rows[i]["kstar"]]
        print(f"  {lab:<6} n={idx.size:3d}  n_real {nr[m].min():5d}-{nr[m].max():5d}  "
              f"reach {R[m].mean():.4f}  contributes {R[m].sum()/R.sum()*100:5.1f}% of the mean"
              f"   k* median {np.median(ks[m]):5.0f}   zero-eff ahead of k* "
              f"{np.mean(tf) if tf else float('nan'):.3f}")
    print()

    # ---- the P(k) curve itself ----------------------------------------------------------
    print("THE P(k) CURVE — mean over perturbations of cumulative purity at depth k")
    print(f"  {'k':>5} {'mean P(k)':>10} {'median P(k)':>12} {'% perts >= 0.9':>15} "
          f"{'n perts with n_adj>=k':>22}")
    for k in (1, 2, 3, 4, 5, 7, 10, 15, 20, 30, 50, 75, 100, 150, 200, 300, 500):
        v = [r["purity"][k - 1] for r in rows if r["n_adj"] >= k]
        if not v:
            continue
        print(f"  {k:>5} {np.mean(v):10.4f} {np.median(v):12.4f} "
              f"{100*np.mean(np.array(v) >= P0):14.1f}% {len(v):22d}")
    print()

    np.savez_compressed(
        f"{ROOT}/experiments/171-reach-prefix/work/walk_{arm}.npz",
        n_real=nr, n_adj=na, kstar=ks, R=R, n_sig=nsig, first_err=first,
    )


if __name__ == "__main__":
    main()
