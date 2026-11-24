"""208 — `pds` AND the transfer cosine for the ROUNDING RULE, on PLOLO (162's instrument).

Reproduces 196's `pdsham196.py` row for row (that is the known-answer gate) and adds the
DIRECTION endpoint the mandate asks for: the macro per-target cosine between the DELIVERED
prediction and its own truth, which 204 measured as the sharper instrument on this surface
(resolvable to ~0.2% relative, i.e. ~±0.00006 overall through 197's +0.00365-per-0.01
exchange rate).

Three rules at production's `n L`:

  floor   `T = floor(n L p)`                    -- what `emit_counts` ships
  ham     largest-remainder over the same `n L p`
  inf     no rounding at all                    -- the ceiling

The SPREAD POLICY does not appear here and cannot: `emit_counts` computes `totals` from
`(n, L, profile)` BEFORE `spread` is read, so flat-1.00 leaves the delivered per-gene totals
bit-identical and `pds`/cos exactly invariant.  `emit208.py` asserts that identity per block
(digest, 85/85) rather than this file assuming it.

    uv run python pds208.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
E188 = ROOT / "experiments/188-magnitude-calibration"
E196 = ROOT / "experiments/196-emitted-depth"
sys.path.insert(0, str(E196))
sys.path.insert(0, str(E188))
from budget188 import FOLDS, PDS_SPAN, TS, deliver, load, mass_consistent_effect, pds_of  # noqa: E402
from ham196 import hamilton  # noqa: E402

NL_1X = 400 * 20_109
# 169's cached PLOLO chalA rows at m = inf -- 196's own known-answer gate.
KA_INF = {"k562_gwps": 0.790712, "orion_hct116": 0.790896,
          "orion_hek293t": 0.771621, "kolf2": 0.615317}
# 196's measured rows, which this file must reproduce (pdsham196.json).
KA_196 = {"floor": 0.7390163577417317, "ham": 0.74222366757526}
COS_RATE = 0.3653          # 197/204: overall per unit transfer cosine (+0.00365 per +0.01)


def macro_cos(P, T):
    """204/162's endpoint: the mean over targets of cos(prediction_i, truth_i)."""
    a = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    b = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    return float((a * b).sum(1).mean())


def deliver_r(effects, b, gate, nL, rule):
    out = np.empty_like(effects)
    err = np.zeros(effects.shape[0])
    for i in range(effects.shape[0]):
        g = mass_consistent_effect(np.where(gate, effects[i], 0.0), b, gate, "minnorm")
        p = np.clip(np.expm1(b + g), 0.0, None) / TS
        p = p / float(p.sum())
        x = nL * p
        T = np.floor(x) if rule == "floor" else hamilton(x).astype(np.float64)
        nz = x > 0
        err[i] = float(np.abs((T[nz] - x[nz]) / x[nz]).mean())
        out[i] = np.log1p(TS * (T / T.sum())) - b
    return out, float(err.mean())


def main() -> None:
    z = np.load(E188 / "work/chalframes.npz")
    bA, gA = z["bA"], z["gA"]
    F = {f: load(f) for f in FOLDS}
    res, cos, quant, relerr, boot = {}, {}, {}, {}, {}
    print("PLOLO, challenge frame A, delivered — `pds` and macro cosine vs the ROUNDING RULE\n")
    for f in FOLDS:
        P0 = np.where(gA[None, :], F[f]["pred"], 0.0)
        Tsc = F[f]["truth"][:, F[f]["scored"]]
        inf = deliver(P0, bA, gA)
        res[(f, "inf")] = pds_of(inf[:, F[f]["scored"]], Tsc)
        cos[(f, "inf")] = macro_cos(inf[:, F[f]["scored"]], Tsc)
        percos = {}
        for rule in ("floor", "ham"):
            D, e = deliver_r(P0, bA, gA, NL_1X, rule)
            Ds = D[:, F[f]["scored"]]
            res[(f, rule)] = pds_of(Ds, Tsc)
            cos[(f, rule)] = macro_cos(Ds, Tsc)
            relerr[(f, rule)] = e
            a = Ds / np.maximum(np.linalg.norm(Ds, axis=1, keepdims=True), 1e-300)
            b = Tsc / np.maximum(np.linalg.norm(Tsc, axis=1, keepdims=True), 1e-300)
            percos[rule] = (a * b).sum(1)
            quant[(f, rule)] = (float(np.linalg.norm(Ds - inf[:, F[f]["scored"]]))
                                / float(np.linalg.norm(inf[:, F[f]["scored"]])))
        # 204's PAIRED target bootstrap on the cosine difference -- the same targets are
        # resampled for both rules, so the donor-to-donor spread cancels.
        d = percos["ham"] - percos["floor"]
        rng = np.random.default_rng(208)
        bs = np.array([d[rng.integers(0, d.size, d.size)].mean() for _ in range(4000)])
        boot[f] = dict(mean=float(d.mean()), lo=float(np.quantile(bs, 0.025)),
                       hi=float(np.quantile(bs, 0.975)),
                       p_gt0=float((bs > 0).mean()), n=int(d.size))
        print(f"  [{f}] {F[f]['pred'].shape[0]} perturbations, "
              f"{F[f]['scored'].size} scored genes   dcos {d.mean():+.6f} "
              f"[{boot[f]['lo']:+.6f}, {boot[f]['hi']:+.6f}]  p>0 {boot[f]['p_gt0']:.3f}",
              flush=True)

    print("\n--- KNOWN-ANSWER GATE 1: m = inf must reproduce 169's cached PLOLO chalA rows ---")
    worst = max(abs(res[(f, "inf")] - KA_INF[f]) for f in FOLDS)
    for f in FOLDS:
        print(f"  {f:<16} 169 {KA_INF[f]:.6f}   here {res[(f,'inf')]:.6f}   "
              f"|d| {abs(res[(f,'inf')]-KA_INF[f]):.2e}")
    print(f"  worst |d pds| = {worst:.2e}")
    if worst > 2e-6:
        raise SystemExit("ABORT (rule 10): PLOLO identity row does not reproduce")

    print("\n--- KNOWN-ANSWER GATE 2: the two rules must reproduce 196's own pdsham196 rows ---")
    bad = []
    for rule in ("floor", "ham"):
        mu = float(np.mean([res[(f, rule)] for f in FOLDS]))
        d = abs(mu - KA_196[rule])
        ok = d <= 2e-6
        print(f"  {rule:<6} 196 {KA_196[rule]:.6f}   here {mu:.6f}   |d| {d:.2e}  "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            bad.append(rule)
    if bad:
        raise SystemExit(f"ABORT (rule 10): {bad}")

    print(f"\n{'rule':>7}" + "".join(f"{f[:11]:>13}" for f in FOLDS) +
          f"{'MEAN pds':>11}{'d pds':>10}{'pds ovr':>10}"
          f"{'MEAN cos':>11}{'d cos':>11}{'rel %':>8}{'cos ovr':>10}")
    rows = []
    base_p = float(np.mean([res[(f, "floor")] for f in FOLDS]))
    base_c = float(np.mean([cos[(f, "floor")] for f in FOLDS]))
    for rule in ("floor", "ham", "inf"):
        vp = [res[(f, rule)] for f in FOLDS]
        vc = [cos[(f, rule)] for f in FOLDS]
        mp, mc = float(np.mean(vp)), float(np.mean(vc))
        dp, dc = mp - base_p, mc - base_c
        rows.append(dict(rule=rule, pds=mp, d_pds=dp, pds_overall=dp / PDS_SPAN / 6,
                         cos=mc, d_cos=dc, cos_rel=dc / base_c,
                         cos_overall=dc * COS_RATE,
                         per_fold_pds={f: res[(f, rule)] for f in FOLDS},
                         per_fold_cos={f: cos[(f, rule)] for f in FOLDS},
                         quant=float(np.mean([quant[(f, rule)] for f in FOLDS]))
                         if rule != "inf" else 0.0))
        print(f"{rule:>7}" + "".join(f"{v:>13.6f}" for v in vp) +
              f"{mp:>11.6f}{dp:>+10.6f}{dp/PDS_SPAN/6:>+10.6f}"
              f"{mc:>11.6f}{dc:>+11.6f}{100*dc/base_c:>+8.2f}{dc*COS_RATE:>+10.6f}")

    print(f"\n  PDS_SPAN = {PDS_SPAN:.5f}; 162's PLOLO LOO sigma 0.00149 raw pds = "
          f"{0.00149/PDS_SPAN/6:.6f} overall.")
    print(f"  cosine: 204's resolution on this surface is ~0.20% relative "
          f"(~{0.002*base_c*COS_RATE:.6f} overall); ham moves it "
          f"{100*(rows[1]['d_cos'])/base_c:+.2f}% relative.")
    print("\n  PAIRED TARGET BOOTSTRAP on d cos (ham - floor), per fold:")
    for f in FOLDS:
        b = boot[f]
        print(f"    {f:<16} {b['mean']:+.6f}  95% [{b['lo']:+.6f}, {b['hi']:+.6f}]  "
              f"p>0 {b['p_gt0']:.4f}  n {b['n']}")
    json.dump(dict(rows=rows, bootstrap=boot, cos_rate=COS_RATE, pds_span=PDS_SPAN),
              open(WORK / "pds208.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
