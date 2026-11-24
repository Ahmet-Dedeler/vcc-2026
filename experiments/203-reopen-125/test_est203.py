"""203 — the re-implementation's known-answer suite (rule 10).

`est203` re-implements 125's estimator on the production gene axis.  A re-implementation
that is not bit-equal to the original measures a different object, and the numbers would
look exactly as reasonable as the right ones.  So every function is checked against
125's own `run.py` on 125's own frozen inputs:

  T1  row scales and normalised observations == `transform_training`
  T2  raw moments                            == `fit_raw_moments`
  T3  EB-shrunk parameters                   == `posterior_parameters`
  T4  posterior effect / confidence          == `posterior_prediction`
  T5  the 122 comparator                     == `equal_mean_agreement`
  T6  the F-general moment code reduces to the F=2 code exactly
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import est203  # noqa: E402

sys.path.insert(0, str(ROOT / "experiments/125-hierarchical-panel-posterior"))
spec = importlib.util.spec_from_file_location(
    "e125run", ROOT / "experiments/125-hierarchical-panel-posterior/run.py")
E125 = importlib.util.module_from_spec(spec)
sys.modules["e125run"] = E125
spec.loader.exec_module(E125)

W122 = ROOT / "experiments/122-panel-response-atlas/work/source-only-phase1-v1"
F125 = ROOT / "experiments/125-hierarchical-panel-posterior/work/preregistered-freeze-v1"

effects = np.asarray(np.load(W122 / "family-effects.npy"), dtype=np.float64)
valid = np.asarray(np.load(W122 / "family-valid.npy"), dtype=bool)
folds = np.load(F125 / "target-folds.npy")

FAILS = []


def check(name, a, b, tol=0.0):
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = float(np.max(np.abs(a - b))) if a.size else 0.0
    ok = d <= tol
    print(f"  {'PASS' if ok else 'FAIL'}  {name:<46} max abs diff {d:.3e}")
    if not ok:
        FAILS.append(name)


for held in range(3):
    tr = tuple(i for i in range(3) if i != held)
    ref = E125.transform_training(effects[list(tr)], valid[list(tr)],
                                  target_permutation=None, response_permutation=None)
    mine = est203.prepare(ref.effects, valid_mask=ref.valid)
    print(f"\nheld family {held}")
    check("T1 row scales", mine.row_scale, ref.row_scales)
    check("T1 normalized", mine.normalized, ref.normalized)
    check("T1 valid mask", mine.valid, ref.valid)

    for k in (0, 3):
        test = folds == k
        rm = E125.fit_raw_moments(ref, ~test, test)
        mm = est203.moments(mine, ~test)
        check(f"T2 pair count (fold {k})", mm["n"], rm.pair_count)
        check(f"T2 prior mean raw (fold {k})", mm["prior_mean_raw"], rm.prior_mean_raw)
        check(f"T2 latent raw (fold {k})", mm["latent"], rm.latent_raw)
        check(f"T2 noise0 raw (fold {k})", mm["noise"][0], rm.noise0_raw)
        check(f"T2 noise1 raw (fold {k})", mm["noise"][1], rm.noise1_raw)
        check(f"T2 global latent (fold {k})", mm["g_latent"], rm.global_latent)
        check(f"T2 global total (fold {k})", mm["g_total"], rm.global_total)

        for name, s in est203.SETTINGS.items():
            ps = E125.PosteriorSetting(name=name, pseudo_count=s["pseudo"],
                                       variance_floor_fraction=s["floor_frac"],
                                       prior_scale=s["prior_scale"])
            rp = E125.posterior_parameters(rm, ps)
            mp = est203.parameters(mm, s)
            check(f"T3 prior mean [{name}] (fold {k})", mp["prior_mean"], rp.prior_mean)
            check(f"T3 latent var [{name}] (fold {k})", mp["latent"], rp.latent_variance)
            check(f"T3 noise0 [{name}] (fold {k})", mp["noise"][0], rp.noise0_variance)
            check(f"T3 noise1 [{name}] (fold {k})", mp["noise"][1], rp.noise1_variance)
            re_, rc, rs = E125.posterior_prediction(ref, (0, 1), rp)
            me_, mc, ms = est203.posterior(mine, mp)
            # 125 forces the target-self diagonal to zero after the posterior; est203
            # relies on production having already zeroed the panel columns, so compare
            # off-diagonal only.
            off = ~np.eye(300, dtype=bool)
            check(f"T4 effect [{name}] (fold {k})", me_[off], re_[off])
            check(f"T4 confidence [{name}] (fold {k})", mc[off], rc[off])
            check(f"T4 support [{name}] (fold {k})", ms[off], rs[off])

    be, bc, bs = E125.equal_mean_agreement(ref, (0, 1))
    me, mc, msup = est203.equal_mean_agreement(mine)
    off = ~np.eye(300, dtype=bool)
    check("T5 comparator effect", me[off], be[off])
    check("T5 comparator confidence", mc[off], bc[off])
    check("T5 comparator support", msup[off], bs[off])

# T6 -- F-general reduction: three identical copies of a two-source panel must give the
# same latent/noise as the two-source fit (mean of pairwise covs == the single cov).
print("\nT6 F-general reduction")
rng = np.random.default_rng(7)
e2 = rng.standard_normal((2, 40, 60))
v2 = np.ones_like(e2, dtype=bool)
p2 = est203.prepare(e2, valid_mask=v2)
m2 = est203.moments(p2, np.ones(40, dtype=bool))
e3 = np.concatenate([e2, e2[:1]], axis=0)
p3 = est203.prepare(e3, valid_mask=np.ones_like(e3, dtype=bool))
m3 = est203.moments(p3, np.ones(40, dtype=bool))
# with a duplicated source 0, mean pairwise cov = (c01 + v0 + c01)/3
c01 = (np.where(v2[0] & v2[1], p2.normalized[0] * p2.normalized[1], 0.0).sum(0)
       - p2.normalized[0].sum(0) * p2.normalized[1].sum(0) / 40) / 39
v0 = (np.square(p2.normalized[0]).sum(0) - p2.normalized[0].sum(0) ** 2 / 40) / 39
check("T6 3-source cov == (c01+v0+c01)/3",
      m3["latent"], np.maximum((2 * c01 + v0) / 3.0, 0.0), tol=1e-12)
check("T6 2-source cov == c01", m2["latent"], np.maximum(c01, 0.0), tol=1e-12)

print("\n" + ("ALL KNOWN-ANSWER ROWS PASS" if not FAILS else f"FAILURES: {FAILS}"))
sys.exit(1 if FAILS else 0)
