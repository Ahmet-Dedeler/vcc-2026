"""193 — known-answer row for the nonlinear-shrinkage implementations themselves.

`analytical_nls` and `rie_lp` are new code.  Before any statement like "nonlinear shrinkage
loses on PLOLO" is allowed to stand, the estimators have to be shown to WORK on data where the
right answer is known -- otherwise a broken port and a real negative result are
indistinguishable (rule 10's corollary).

The test: draw `X` (n x p) with a known population covariance `Sigma`, form the sample
covariance, and compare the Frobenius (percentage relative improvement in average loss, PRIAL)
of each estimator against the ORACLE rotationally-invariant estimator
    d_oracle_i = u_i^T Sigma u_i
which is the best any `V f(ev) V^T` estimator can do.  Ledoit & Wolf 2020 report PRIAL well
above 90 % for the analytical estimator in this regime; anything near it is a working port.

Two population spectra are used:
  * `spike`   -- a few large eigenvalues over a flat bulk (the model MP-clipping assumes)
  * `powerlaw`-- eigenvalues decaying like a power law over 6 decades, i.e. the shape the DONOR
                 covariance actually has.  Reporting both is how the PLOLO result gets its
                 mechanism rather than just its sign.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import maps193  # noqa: E402


def prial(d_est, d_oracle, ev):
    """Percentage relative improvement in average loss vs the sample eigenvalues."""
    l_sample = float(((ev - d_oracle) ** 2).sum())
    l_est = float(((d_est - d_oracle) ** 2).sum())
    return 100.0 * (l_sample - l_est) / l_sample if l_sample > 0 else float("nan")


def trial(pop, p, n, seed):
    rng = np.random.default_rng(seed)
    tau = pop(p)
    A = np.sqrt(tau)[None, :] * rng.standard_normal((n, p))
    C = (A.T @ A) / n
    ev, V = np.linalg.eigh(C)
    ev = np.maximum(ev, ev.max() * 1e-12)
    # oracle RIE: d_i = u_i^T Sigma u_i, and Sigma = diag(tau) here
    d_or = (V * V * tau[:, None]).sum(0)
    res = {}
    res["sample"] = prial(ev, d_or, ev)
    res["nls-LW2020"] = prial(maps193.analytical_nls(ev, n), d_or, ev)
    res["rie-LP-eta0.5"] = prial(maps193.rie_lp(ev, n, 0.5), d_or, ev)
    res["clipMP"] = prial(maps193.clip_mp(ev, p / n), d_or, ev)
    best = max(("ridge%g" % lam, prial(maps193.ridge(ev, lam), d_or, ev))
               for lam in (1e-4, 1e-3, 1e-2, 3e-2, 1e-1))
    res[f"best-ridge({best[0]})"] = best[1]
    return res, ev, d_or


def spike(p):
    tau = np.ones(p)
    tau[: max(1, p // 50)] = 20.0
    return tau


def powerlaw(p):
    return np.logspace(0, -6, p)


def main() -> None:
    out = {}
    for name, pop in (("spike", spike), ("powerlaw", powerlaw)):
        for p, n in ((500, 1000), (1000, 1600)):
            res, ev, d_or = trial(pop, p, n, seed=0)
            key = f"{name} p={p} n={n} (c={p/n:.2f})"
            out[key] = res
            print(f"\n{key}   PRIAL vs oracle RIE (higher is better, sample = 0)")
            for k, v in res.items():
                print(f"   {k:24s} {v:8.2f} %")
            print(f"   spectrum span {ev.max()/ev.min():.3e}")
    (HERE / "work/nlscheck.json").write_text(json.dumps(out, indent=1))

    ok = out["spike p=500 n=1000 (c=0.50)"]["nls-LW2020"] > 60.0
    print("\nPORT GATE (spike, analytical NLS PRIAL > 60%): " + ("PASS" if ok else "FAIL"))
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
