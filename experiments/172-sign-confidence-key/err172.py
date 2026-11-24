"""172 — how big is a `reach` difference on h1-test before it means anything?

`anatomy172.py` measures that **94.4% of h1-test's `reach` under `|H|` is contributed by the
32 targets whose `n_real` is <= 20**, at a median `k*` of 3.  So the whole endpoint rests on a
few hundred (target, gene) sign decisions, and a table of keys differing by +0.01 is not
readable without an error bar.  Nobody in this project has put one on `reach` from the target
axis (157 §7 quotes SEED spreads of 0.008-0.011, which is a different and smaller source of
variation: it holds the target set fixed and re-draws the emitter).

Two estimators, both over the TARGET axis, which is the axis `reach` averages over:

  * jackknife-over-targets SE of `reach(key)` for one key;
  * a PAIRED bootstrap over targets of `reach(A) - reach(B)`, which is the quantity every row
    of the results table actually claims.  Paired matters: the same 32 small-pool targets
    dominate both keys, so the unpaired SE (~0.036) is 4-5x the SE of the difference.

`reach` is a mean of per-target `R_p`, so resampling targets is exact -- no re-walk is needed
beyond computing `R_p` once per key.

Two comparators, and they answer different questions:
  vs `|H|`      -- is the learned key better than the best simple key we already have?
  vs `|cons|`   -- did LEARNING add anything on top of the exact same information the model
                   was given?  `|cons|` is the un-whitened consensus magnitude, i.e. the model
                   with all its coefficients collapsed onto one feature.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from anatomy172 import kstar          # noqa: E402
from harness172 import H1Test, check  # noqa: E402

B = 20000


def main() -> None:
    T = H1Test()
    print("known-answer rows:")
    check(T, verbose=False)
    print("  OK (harness172.check)")
    Z = np.load(os.path.join(HERE, "work", "allkeys.npz"), allow_pickle=True)
    order = ["production p-values", "|cons| (un-whitened)", "geo x concord", "|H| (whitened)"]
    order += [k for k in Z.files if k not in order + ["genes", "targets", "ORACLE"]
              and not k.startswith("null")]
    nulls = sorted(k for k in Z.files if k.startswith("null"))
    order += nulls + ["ORACLE"]
    R = {k: kstar(T, np.asarray(Z[k], dtype=np.float64))[:, 2] for k in order}
    # debench's production order is a 5-way lexsort, not a single scalar; `allkeys`' one-scalar
    # approximation of it reads 0.17941 against the true 0.18283, so the exact row is added here.
    R = {"production (debench order)": kstar(T, None)[:, 2], **R}
    order = ["production (debench order)"] + order
    n = len(next(iter(R.values())))
    idx = np.random.default_rng(7).integers(0, n, size=(B, n))

    out = {}
    for comp in ("|H| (whitened)", "|cons| (un-whitened)", "M1 logistic @ C=0.1 w=top25"):
        base = R[comp]
        bb = base[idx].mean(1)
        print(f"\n{n} scored targets.  jackknife SE over targets, and the PAIRED bootstrap "
              f"({B} resamples) of the difference against {comp}:\n")
        print(f"{'key':<40}{'reach':>10}{'jackSE':>9}{'d':>10}{'SE(d)':>9}"
              f"{'95% CI of d':>24}{'P(d>0)':>9}")
        for k in order:
            r = R[k]
            m = r.mean()
            jk = np.array([np.delete(r, i).mean() for i in range(n)])
            se = np.sqrt((n - 1) / n * ((jk - jk.mean()) ** 2).sum())
            bd = r[idx].mean(1) - bb
            lo, hi = np.quantile(bd, [0.025, 0.975])
            print(f"{k:<40}{m:>10.5f}{se:>9.5f}{m-base.mean():>+10.5f}{bd.std():>9.5f}"
                  f"{f'[{lo:+.5f}, {hi:+.5f}]':>24}{float((bd > 0).mean()):>9.3f}")
            out.setdefault(k, {})["reach"] = float(m)
            out[k]["jack_se"] = float(se)
            out[k][f"vs {comp}"] = dict(d=float(m - base.mean()), se=float(bd.std()),
                                        ci=[float(lo), float(hi)],
                                        p_gt0=float((bd > 0).mean()))

    nr = np.asarray([int(T.S.real_sig[p].sum()) for p in range(len(T.perts))
                     if T.S.real_sig[p].sum() > 0])
    small = nr <= 20
    print(f"\nrestricted to the {int(small.sum())} targets that carry 94.4% of the mass "
          f"(n_real <= 20), and to the {int((~small).sum())} that carry 5.6%:")
    print(f"{'key':<40}{'n_real<=20':>12}{'n_real>20':>12}")
    for k in order:
        print(f"{k:<40}{R[k][small].mean():>12.5f}{R[k][~small].mean():>12.5f}")
        out[k]["small_pool"] = float(R[k][small].mean())
        out[k]["large_pool"] = float(R[k][~small].mean())
    with open(os.path.join(HERE, "work", "err172.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
