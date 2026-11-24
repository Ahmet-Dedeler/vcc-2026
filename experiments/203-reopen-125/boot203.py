"""203 — paired per-target uncertainty on the `pds` and cosine deltas.

The four PLOLO folds give only four replicates, so the fold spread is a weak interval.
`pds` is a mean over per-target retrieval ranks and the cosine is a mean over per-target
cosines, so both decompose per target, and the incumbent and the arm are evaluated on the
SAME targets.  This resamples targets within fold (2,000 paired draws, folds weighted
equally as in the pre-committed instrument) and reports the delta's 95% interval.

KNOWN-ANSWER ROWS
  KA1  the per-target mean of the saved ranks must reproduce `arms203.json`'s fold `pds`
       to 1e-12 for every arm and fold.
  KA2  an arm bootstrapped against ITSELF must give an interval of exactly [0, 0].
"""
from __future__ import annotations

import json
import os
import sys
import warnings
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

warnings.filterwarnings("ignore", message="All-NaN slice encountered")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "src"))

import est203  # noqa: E402
import plolo  # noqa: E402
from plolo import FOLDS, PROD_OPS, Operator, build_fold, cd4_for, deliver, solve_fraction  # noqa: E402

OUT = HERE / "work"
ALPHA, GAMMA, LAM, CD4F, SHRINK = 0.982, 0.65, 0.0, 0.10, 6.7801
NDRAW = 2000


def per_target(P, T):
    A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    B = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S).copy()
    beaten = (S < own[:, None]).sum(1)
    ties = (S == own[:, None]).sum(1) - 1
    return (beaten + 0.5 * ties) / (n - 1), own


def main():
    store = {}
    for held in FOLDS:
        fold = build_fold(held, limit_gib=7.0)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        cd4 = cd4_for(fold["targets"])
        op = Operator(PROD_OPS / "whitener_union_pruned.npz"
                      if FOLDS[held]["op_kind"] == "production" else fold["op_path"])
        z = np.load(OUT / f"sources-{held}.npz", allow_pickle=True)
        panel = est203.prepare(np.asarray(z["effects"], dtype=np.float64),
                               available=z["available"], support=z["support"])
        targets = [str(t) for t in z["targets"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        eq, _, _ = est203.equal_mean_agreement(panel)
        arms = {"a15": pre, "eq122": eq, f"shrink{SHRINK:.2f}": pre / SHRINK}
        for nm in ("eb_broad_prior", "eb_conservative"):
            arms[f"post125:{nm}"] = est203.nested_posterior(
                panel, targets, est203.SETTINGS[nm])[0]

        frames = []
        for ctx in ("A", "B", "C"):
            b = plolo.challenge_bctrl(ctx)
            e = np.expm1(b)
            frames.append((b, (1e6 * e / float(e.sum()) > plolo.MIN_CPM) & (b > 0)))
        for name, eff in arms.items():
            p = np.asarray(eff, dtype=np.float64) * ALPHA
            wh = op.apply(p, fold["wcols"], GAMMA, LAM)
            frac, _ = solve_fraction(wh, cd4, CD4F)
            fin, _ = plolo.E062.add_post_whiten_residual(wh, cd4, frac, fold["wcols"])
            fin[:, fold["panel_positions"]] = 0.0
            r = np.zeros(fin.shape[0])
            c = np.zeros(fin.shape[0])
            for b, g in frames:
                rr, cc = per_target(deliver(fin, b, g, None)[:, fold["scored"]], truth)
                r += rr / 3.0
                c += cc / 3.0
            store[(held, name)] = (r, c)
        print(f"[{held}] done", flush=True)
        del op

    ref = json.loads((OUT / "arms203.json").read_text())
    dref = {(x["fold"], x["arm"]): x for x in ref}
    print("\nKA1  per-target mean vs arms203.json fold `pds`:")
    worst = 0.0
    for (held, name), (r, _) in store.items():
        got = float(r.mean())
        want = float(np.mean([dref[(held, name)]["spaces"][f"chal{c}/delivered"]["pds"]
                              for c in "ABC"]))
        worst = max(worst, abs(got - want))
    print(f"     max abs diff over {len(store)} (fold, arm) cells: {worst:.3e}")

    rng = np.random.default_rng(2030)
    arms = sorted({k[1] for k in store})
    print(f"\n{'='*96}\nPaired target bootstrap, {NDRAW} draws, folds weighted equally "
          f"(the pre-committed instrument)\n{'='*96}")
    print(f"{'arm':<30}{'d pds':>11}{'95% CI':>26}{'P(d>0)':>9}"
          f"{'d cos':>11}{'95% CI':>26}")
    out = {}
    for a in arms:
        dp, dc = [], []
        for _ in range(NDRAW):
            vp, vc = [], []
            for held in FOLDS:
                n = store[(held, a)][0].size
                idx = rng.integers(0, n, n)
                vp.append(store[(held, a)][0][idx].mean() - store[(held, "a15")][0][idx].mean())
                vc.append(store[(held, a)][1][idx].mean() - store[(held, "a15")][1][idx].mean())
            dp.append(np.mean(vp))
            dc.append(np.mean(vc))
        dp, dc = np.array(dp), np.array(dc)
        obs_p = float(np.mean([store[(h, a)][0].mean() - store[(h, "a15")][0].mean()
                               for h in FOLDS]))
        obs_c = float(np.mean([store[(h, a)][1].mean() - store[(h, "a15")][1].mean()
                               for h in FOLDS]))
        lo, hi = np.percentile(dp, [2.5, 97.5])
        clo, chi = np.percentile(dc, [2.5, 97.5])
        out[a] = dict(d_pds=obs_p, ci=[float(lo), float(hi)], p_pos=float((dp > 0).mean()),
                      d_cos=obs_c, cos_ci=[float(clo), float(chi)])
        print(f"{a:<30}{obs_p:>+11.6f}  [{lo:>+9.6f},{hi:>+9.6f}]{(dp>0).mean():>9.3f}"
              f"{obs_c:>+11.6f}  [{clo:>+9.6f},{chi:>+9.6f}]")
    print("KA2  a15 vs itself must read exactly 0 with a [0, 0] interval -- see the a15 row.")
    (OUT / "boot203.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {OUT/'boot203.json'}")


if __name__ == "__main__":
    main()
