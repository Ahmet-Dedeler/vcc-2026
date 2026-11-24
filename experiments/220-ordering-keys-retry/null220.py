"""220 Phase 3 F1 — the permutation control ON THE DELIVERED INSTRUMENT, at the new precision.

The mandate pre-committed: *"190's shuffled regulon scored higher than the real network; 189's
shuffled control tied its call-set arm.  At 77x tighter noise these controls become sharp — run
them, and if a control still ties or beats the arm, the family is dead for good and that is the
result."*

Twelve theta shuffles and six regulon-sign shuffles, each DELIVERED through the shipped
`rc_callset` emitter (`emit220.py`, two runs), against the three-seed real arms.  Both pools.

Known-answer row: `de220x.absH.npz` must be byte-equal to `de220.absH.npz` -- the second run's
baseline is the same deterministic object as the first's, and both are `rc_callset`.

    uv run python null220.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/153-de-bench"))
from debench import ALPHA, Surface                                        # noqa: E402
from six220 import WaldSurface                                            # noqa: E402

E185 = ROOT / "experiments/185-rank-controlled-emission"
SH189 = [(f"shuf189_{i}", "de220" if i <= 4 else "de220x") for i in range(1, 13)]
SH190 = [(f"shuf190_{i}", "de220" if i <= 2 else "de220x") for i in range(1, 7)]
REAL189 = [("k189_l050", "de220"), ("k189_l050_s2", "de220"), ("k189_l050_s3", "de220")]
REAL190 = [("k190_km25", "de220"), ("k190_km25_s2", "de220"), ("k190_km25_s3", "de220")]
BASE = [("absH", "de220"), ("absH_s2", "de220"), ("absH_s3", "de220")]


def main() -> None:
    meta = np.load(E185 / "work/ctrl_meta.npz", allow_pickle=True)
    cpm = meta["cpm"][meta["keep"]]

    a = np.load(WORK / "de220.absH.npz", allow_pickle=True)
    b = np.load(WORK / "de220x.absH.npz", allow_pickle=True)
    same = all(np.array_equal(a[k], b[k]) for k in
               ("pred_p", "pred_padj", "pred_lfc", "pred_z", "digest"))
    print(f"KNOWN-ANSWER ROW (rule 10): the two runs' `absH` are the same object: "
          f"{'PASS' if same else 'FAIL'}")
    if not same and not os.environ.get("NOABORT"):
        raise SystemExit("ABORT: the two emission runs disagree on the baseline")

    out = {}
    for pool in ("native", "wald89"):
        def reach(name, tag):
            p = str(WORK / f"{tag}.{name}.npz")
            s = WaldSurface(p, cpm) if pool == "wald89" else Surface(p)
            return s.reach()

        base = reach(*BASE[0])
        print(f"\n=== POOL: {pool} — DELIVERED, `absH` = {base:.6f} ===")
        seedsp = np.asarray([reach(n, t) - base for n, t in BASE])
        print(f"  the incumbent's own three-seed spread: " +
              " ".join(f"{x:+.5f}" for x in seedsp) + f"   sd {seedsp.std(ddof=1):.5f}")
        res = {}
        for nm, real, shuf in (("189 |H| x exp(0.50 theta)", REAL189, SH189),
                               ("190 ct-prop-k25 lam -0.25", REAL190, SH190)):
            r = np.asarray([reach(n, t) - base for n, t in real])
            s = np.asarray([reach(n, t) - base for n, t in shuf])
            z = (r.mean() - s.mean()) / s.std(ddof=1)
            p = float((s >= r.mean()).mean())
            print(f"\n  {nm}")
            print(f"    REAL    3 seeds  " + " ".join(f"{x:+.5f}" for x in r) +
                  f"   mean {r.mean():+.5f}   seed sem {r.std(ddof=1)/np.sqrt(3):.5f}")
            print(f"    SHUFFLE {len(s)} draws " +
                  " ".join(f"{x:+.5f}" for x in np.sort(s)[::-1]))
            print(f"            mean {s.mean():+.5f}   sd {s.std(ddof=1):.5f}   "
                  f"beats the real arm in {int((s >= r.mean()).sum())} of {len(s)}")
            print(f"    ==>  z {z:+.2f}   P(shuffle >= real) {p:.3f}   "
                  f"{'THE CONTROL TIES OR BEATS THE ARM' if p >= 0.05 else 'separable'}")
            res[nm] = dict(real=list(map(float, r)), real_mean=float(r.mean()),
                           seed_sem=float(r.std(ddof=1) / np.sqrt(3)),
                           shuffle=list(map(float, s)), shuf_mean=float(s.mean()),
                           shuf_sd=float(s.std(ddof=1)), z=float(z), p=p,
                           n_beat=int((s >= r.mean()).sum()))
        out[pool] = dict(base=float(base), seed_spread=list(map(float, seedsp)), arms=res)

    # ---- the two axes, on the DELIVERED arm itself ----------------------------------
    print("\n" + "=" * 92)
    print("THE TWO AXES ON THE DELIVERED ARM — the experiment in one table")
    print("=" * 92)
    sys.path.insert(0, str(HERE))
    import common220 as C220

    def rvec(path):
        S = Surface(str(path))
        v = []
        for i in range(len(S.perts)):
            pool = np.flatnonzero(S.real_sig[i])
            if pool.size == 0:
                continue
            q = pool[S._order(i, pool, None, None)]
            q = q[S.adj[i, q]]
            if q.size == 0:
                v.append(0.0)
                continue
            ok = np.sign(S.pred_lfc[i, q]) == np.sign(S.real_lfc[i, q])
            pur = np.cumsum(ok) / np.arange(1, q.size + 1)
            d = np.flatnonzero(pur >= 0.9)
            v.append((int(d[-1] + 1) if d.size else 0) / int(S.real_sig[i].sum()))
        return np.asarray(v)

    b = rvec(WORK / "de220.absH.npz")
    print(f"{'arm':<14}{'delta':>11}{'paired sem':>12}{'target se':>12}{'ratio':>8}"
          f"{'95% CI over targets':>26}{'P(better)':>11}")
    axes = {}
    for nm in ("k189_l050", "k190_km25"):
        bb = C220.paired_bootstrap(rvec(WORK / f"de220.{nm}.npz"), b)
        # the PAIRED seed sem: each seed's arm against its OWN `absH` seed (214's convention)
        pr = np.asarray([Surface(str(WORK / f"de220.{nm}{t}.npz")).reach()
                         - Surface(str(WORK / f"de220.absH{t}.npz")).reach()
                         for t in ("", "_s2", "_s3")])
        ss = float(pr.std(ddof=1) / np.sqrt(3))
        print(f"{nm:<14}{bb['mean']:>+11.6f}{ss:>11.6f}{bb['se']:>12.5f}"
              f"{bb['se']/ss:>8.0f}x"
              f"{f'[{bb[chr(108)+chr(111)]:+.5f}, {bb[chr(104)+chr(105)]:+.5f}]':>25}"
              f"{bb['p_better']:>11.3f}")
        axes[nm] = dict(seed_sem=ss, **bb)
    out["axes"] = axes

    json.dump(out, open(WORK / "null220.json", "w"), indent=1, default=float)
    print(f"\nwrote {WORK/'null220.json'}")


if __name__ == "__main__":
    main()
