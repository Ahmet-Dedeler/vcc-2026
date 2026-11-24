"""169 — `reach` (and the other three DE members) with an error bar, from MATCHED seeds.

`reach` is the member rule 13 charges a blanket -0.0146 scaled toll for whenever the emitted
profile moves, and it is also the noisiest member this project measures: 157 recorded a
same-policy re-seed spread of 0.008-0.011 raw on h1-test, which is the same size as every arm
delta ever measured here.  A one-seed reading of `reach` is not a measurement.

So the incumbent and the arm are emitted at FIVE MATCHED SEEDS and compared pairwise -- the
reference DE table, the control block, the gene axis and the seed are all identical between
the two members of a pair, so the only thing that differs is `q`.  Two estimators are
reported:

  * the seed-paired difference of the aggregate `reach`, `n = 5` pairs;
  * the perturbation-level paired difference of the seed-averaged per-perturbation reach,
    `n = 85` perturbations, which is the higher-powered of the two.

Both are then converted to board units through the a8<->a15 transfer factor (157 §1), and the
resulting confidence interval is what replaces -- or fails to replace -- rule 13's blanket toll.

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
sys.path.insert(0, HERE)
from debench import ANCHOR, Surface  # noqa: E402
from six169 import PowerSurface, overlap_ratio, precision  # noqa: E402

GROUPS = {
    1.00: ["q1.00_mf", "q1.00_mf_s2", "q1.00_mf_s3", "q1.00_mf_s4", "q1.00_mf_s5"],
    1.15: ["q1.15_mf", "q1.15_mf_b", "q1.15_mf_s3", "q1.15_mf_s4", "q1.15_mf_s5"],
    1.25: ["q1.25_mf", "q1.25_mf_s2", "q1.25_mf_s3"],
    # the lfc-space arm (`sharpen_lfc`, h-space row norm restored), matched to the same seeds
    11.15: ["q1.15L_mf", "q1.15L_mf_s2", "q1.15L_mf_s3", "q1.15L_mf_s4", "q1.15L_mf_s5"],
    11.25: ["q1.25L_mf", "q1.25L_mf_s2"],
    # and the lfc-NORM-preserving variant, which carries a 17% drop in `A` (two variables)
    21.15: ["q1.15Lnorm_mf", "q1.15Lnorm_mf_s2", "q1.15Lnorm_mf_s3", "q1.15Lnorm_mf_s4",
            "q1.15Lnorm_mf_s5"],
}
K_REACH = (0.1890 - 0.1811) / (0.2029 - 0.1828)     # board / h1-test, a8 <-> a15 (157 §1)
SPAN_R = ANCHOR["reach"][1] - ANCHOR["reach"][0]
SPAN_F = ANCHOR["fid"][1] - ANCHOR["fid"][0]
FACTORS = (1.0, 2.0, 3.5)


def per_pert_reach(path, factor=1.0):
    s = PowerSurface(path, factor)
    out = np.zeros(len(s.perts))
    for p in range(len(s.perts)):
        pool = np.flatnonzero(s.real_sig[p])
        if pool.size == 0:
            continue
        o = pool[s._order(p, pool, None, None)]
        o = o[s.adj[p, o]]
        if o.size == 0:
            continue
        ok = np.sign(s.pred_lfc[p, o]) == np.sign(s.real_lfc[p, o])
        pur = np.cumsum(ok) / np.arange(1, ok.size + 1)
        d = np.flatnonzero(pur >= 0.9)
        out[p] = (int(d[-1]) + 1) / int(s.real_sig[p].sum()) if d.size else 0.0
    return out


def main() -> None:
    have = {}
    for q, names in GROUPS.items():
        got = [n for n in names if os.path.exists(os.path.join(HERE, f"de169.{n}.npz"))]
        if len(got) >= 2:
            have[q] = got
        elif q in (1.00, 1.15):
            sys.exit(f"q={q}: only {got} present -- run realize169.py with the seed arms")
        else:
            print(f"  (skipping q={q}: {len(got)} arms present)")
    LABEL = lambda q: ("h-space q=%.2f" % q if q < 10 else
                   ("lfc-space q=%.2f" % (q - 10) if q < 20 else
                    "lfc-NORM q=%.2f" % (q - 20)))  # noqa: E731
    print("matched seeds: " + ", ".join(f"{LABEL(q)}:{len(v)}" for q, v in have.items()))

    R = {}          # (q, seed index) -> per-perturbation reach vector, per factor
    agg = {}
    for q, names in have.items():
        for i, n in enumerate(names):
            path = os.path.join(HERE, f"de169.{n}.npz")
            for f in FACTORS:
                R[(q, i, f)] = per_pert_reach(path, f)
            s = Surface(path)
            agg[(q, i)] = dict(reach=float(R[(q, i, 1.0)].mean()),
                               precision=precision(s), recall=overlap_ratio(s),
                               nmae=s.nmae(), pds=s.pds(),
                               n_pred=float(np.median(s.counts()[1])))

    print(f"\nper-seed aggregate `reach` (f = 1, the board-calibrated setting)")
    print(f"{'arm':>16} " + " ".join(f"{'seed' + str(i + 1):>10}" for i in range(5))
          + f" {'mean':>10} {'sd':>9}")
    for q in sorted(have):
        v = [R[(q, i, 1.0)].mean() for i in range(len(have[q]))]
        print(f"{LABEL(q):>16} " + " ".join(f"{x:10.5f}" for x in v)
              + " " * 11 * (5 - len(v)) + f" {np.mean(v):10.5f} {np.std(v, ddof=1):9.5f}")

    NM = {q: min(len(have[q]), len(have[1.00])) for q in have}
    print(f"\nSEED-PAIRED d reach vs q = 1.00, matched pairs per q: "
          + ", ".join(f"{LABEL(q)}:{NM[q]}" for q in sorted(NM) if q != 1.00))
    print(f"{'arm':>16}{'f':>6}{'mean d':>10}{'sd of d':>10}{'sem':>10}{'t':>8}"
          f"{'95% CI raw':>22}{'95% CI overall':>24}")
    out = {}
    for q in sorted(have):
        if q == 1.00:
            continue
        nmin = NM[q]
        tcrit = {2: 12.706, 3: 4.303, 4: 3.182, 5: 2.776}[nmin]
        for f in FACTORS:
            d = np.array([R[(q, i, f)].mean() - R[(1.00, i, f)].mean() for i in range(nmin)])
            sem = d.std(ddof=1) / np.sqrt(d.size)
            lo, hi = d.mean() - tcrit * sem, d.mean() + tcrit * sem
            g = K_REACH / SPAN_R / 6
            print(f"{LABEL(q):>16}{f:>6.1f}{d.mean():>+10.5f}{d.std(ddof=1):>10.5f}{sem:>10.5f}"
                  f"{d.mean() / max(sem, 1e-12):>8.2f}"
                  f"  [{lo:+.5f}, {hi:+.5f}]  [{lo * g:+.6f}, {hi * g:+.6f}]")
            if f == 1.0:
                out[q] = dict(mean_raw=float(d.mean()), sem_raw=float(sem),
                              mean_overall=float(d.mean() * g),
                              lo=float(lo * g), hi=float(hi * g))

    print(f"\nPERTURBATION-PAIRED d reach (f = 1): the seed-averaged per-perturbation reach, "
          f"85 pairs")
    print(f"{'arm':>16}{'mean d':>10}{'sd':>9}{'sem':>10}{'t':>8}{'n>0':>7}{'n<0':>7}"
          f"{'95% CI overall':>26}")
    for q in sorted(have):
        if q == 1.00:
            continue
        nmin = NM[q]
        a = np.mean([R[(q, i, 1.0)] for i in range(nmin)], axis=0)
        b = np.mean([R[(1.00, i, 1.0)] for i in range(nmin)], axis=0)
        d = a - b
        sem = d.std(ddof=1) / np.sqrt(d.size)
        g = K_REACH / SPAN_R / 6
        lo, hi = d.mean() - 1.989 * sem, d.mean() + 1.989 * sem       # t_{.975, 84}
        print(f"{LABEL(q):>16}{d.mean():>+10.5f}{d.std(ddof=1):>9.5f}{sem:>10.5f}"
              f"{d.mean() / max(sem, 1e-12):>8.2f}{int((d > 0).sum()):>7}"
              f"{int((d < 0).sum()):>7}   [{lo * g:+.6f}, {hi * g:+.6f}]")

    print(f"\nRULE 13's BLANKET TOLL for an arm that changes the emitted profile: "
          f"-0.014569 scaled = {-0.014569 / 6:+.6f} overall")
    print(f"h1-test -> board transfer factor {K_REACH:.4f}; reach span {SPAN_R:.5f}")

    print(f"\nthe other three members, seed-paired")
    print(f"{'arm':>16}{'d precision':>13}{'sem':>9}{'d recall':>10}{'sem':>9}"
          f"{'d nmae':>10}{'sem':>9}{'d pds':>10}{'sem':>9}{'d n_pred':>10}")
    for q in sorted(have):
        if q == 1.00:
            continue
        nmin = NM[q]
        cells = []
        for k in ("precision", "recall", "nmae", "pds"):
            d = np.array([agg[(q, i)][k] - agg[(1.00, i)][k] for i in range(nmin)])
            cells.append((d.mean(), d.std(ddof=1) / np.sqrt(d.size)))
        dn = np.mean([agg[(q, i)]["n_pred"] - agg[(1.00, i)]["n_pred"] for i in range(nmin)])
        print(f"{LABEL(q):>16}" + "".join(f"{m:>+13.5f}{s:>9.5f}" if i == 0 else
                                     f"{m:>+10.5f}{s:>9.5f}" for i, (m, s) in enumerate(cells))
              + f"{dn:>+10.0f}")
    for q in sorted(have):
        if q == 1.00:
            continue
        n = NM[q]
        dp = np.mean([agg[(q, i)]['precision'] - agg[(1.00, i)]['precision']
                      for i in range(n)])
        dr = np.mean([agg[(q, i)]['recall'] - agg[(1.00, i)]['recall'] for i in range(n)])
        # `jac` rebuilt on the board's own counts (rule 9): n_real 240, n_pred 4424.
        NR, NP = 240.0, 4424.0
        r0 = np.mean([agg[(1.00, i)]['recall'] for i in range(n)])
        j0 = r0 * NR / (NP + NR - r0 * NR)
        j1 = (r0 + dr) * NR / (NP + NR - (r0 + dr) * NR)
        print(f"  {LABEL(q)}: d fid(precision) overall {dp / SPAN_F / 6:+.6f} | "
              f"d jac overall {(j1 - j0) / (ANCHOR['jac'][1] - ANCHOR['jac'][0]) / 6:+.6f}")

    with open(os.path.join(WORK, "reach169.json"), "w") as fh:
        json.dump({str(k): v for k, v in out.items()}, fh, indent=1, default=float)


if __name__ == "__main__":
    sys.exit(main())
