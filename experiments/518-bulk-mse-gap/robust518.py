"""518 — robustness of the strength key (exact PLOLO rows on the dumps): permutation control, target bootstrap, and the
strength statistic's threshold.

  CTRL   the meta strength ranks permuted over targets (seeded), theta 1.25, on 446bj2's key -> must lose.
  BOOT   4000 target resamples of the per-target (x, a) arrays of base and arm (the aggregate cosine's ratio), per fold.
  THR    the robust count at |z| > 3 / 4 / 5 noise units (analytic route, LOFO theta), to show the statistic is not tuned.

    uv run python experiments/518-bulk-mse-gap/robust518.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import json  # noqa: E402
import warnings  # noqa: E402

import numpy as np  # noqa: E402

from lofo518 import crank, fit, gain  # noqa: E402
from score518 import Fold  # noqa: E402
from skey518 import zfeat  # noqa: E402

warnings.filterwarnings("ignore")
WORK = HERE / "work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
CANCER = FOLDS[:3]


def cos_from(arr, d, idx=None):
    x, a = arr["x"], arr["a"]
    if idx is not None:
        x, a, d = x[:, idx], a[:, idx], d[idx]
    return float(np.mean([x[j].sum() / np.sqrt(a[j].sum() * d.sum()) for j in range(3)]))


def main():
    out = {}
    for f in FOLDS:
        F = Fold(f)
        d = (F.truth * F.truth).sum(1)
        b, ab = F.score()
        z = zfeat(F, "meta")
        rng = np.random.default_rng(518)
        s_arm, aa = F.score(w=F.w_base * np.exp(1.25 * z))
        s_ctl, _ = F.score(w=F.w_base * np.exp(1.25 * z[rng.permutation(z.size)]))
        n = z.size
        boots = []
        rb = np.random.default_rng(5180)
        for _ in range(4000):
            idx = rb.integers(0, n, n)
            boots.append(cos_from(aa, d, idx) / cos_from(ab, d, idx) - 1)
        lo, hi = np.quantile(boots, [0.025, 0.975])
        out[f] = dict(r_arm=s_arm["cos_agg"] / b["cos_agg"] - 1, r_ctrl=s_ctl["cos_agg"] / b["cos_agg"] - 1,
                      boot95=[float(lo), float(hi)], p_le0=float(np.mean(np.asarray(boots) <= 0)),
                      dpds_ctrl=s_ctl["pds"] - b["pds"])
        print(f"{f:14s} meta theta 1.25: r {100 * out[f]['r_arm']:+.2f} % [95 % target CI {100 * lo:+.2f}, {100 * hi:+.2f}], "
              f"P(<=0) {out[f]['p_le0']:.3f} | CONTROL permuted: r {100 * out[f]['r_ctrl']:+.2f} %, dpds {out[f]['dpds_ctrl']:+.5f}", flush=True)
    # threshold sensitivity (analytic, LOFO over the cancer folds; the sources' strength only)
    import strength518 as S5
    base_fn = S5.robust_strength
    thr_out = {}
    for thr in (3.0, 4.0, 5.0):
        def rs(path, targets=S5.PANEL, _t=thr):
            perts = [l.strip() for l in (path / "perts.txt").read_text().splitlines() if l.strip()]
            genes = [l.strip() for l in (path / "genes.txt").read_text().splitlines() if l.strip()]
            bb = np.load(path / "b_ctrl.npy").astype(np.float64)
            cpm = 1e6 * np.expm1(bb) / np.expm1(bb).sum()
            keep = (cpm > 10.0) & np.asarray([g not in S5.PSET for g in genes])
            pi = {p: i for i, p in enumerate(perts)}
            rows = [pi[t] for t in targets if t in pi]
            cov = [t for t in targets if t in pi]
            delta = np.load(path / "delta.npy", mmap_mode="r")
            R = np.asarray(delta[np.asarray(sorted(rows))], np.float64)[np.argsort(np.argsort(rows))][:, keep]
            sc = 1.4826 * np.median(np.abs(R - np.median(R, 0)), 0)
            ok = sc > 0
            Z = R[:, ok] / sc[ok]
            nz = np.median(np.abs(Z), 1) / 0.6745
            cnt = (np.abs(Z) > _t * nz[:, None]).sum(1)
            o = dict(zip(cov, np.log1p(cnt).tolist()))
            return np.asarray([o.get(t, np.nan) for t in targets])
        strength = {k: rs(p) for k, p in S5.BANKS.items() if k in ("k562_gwps", "orion_hct116", "orion_hek293t")}
        data = {}
        for f in FOLDS:
            zz = np.load(S5.W540 / f"arrays540-{f}.npz")
            x, a, dd = zz[f"{S5.BASE}__x"], zz[f"{S5.BASE}__a"], zz["d"]
            targets = [str(t) for t in np.load(WORK / f"dump518-{f}.npz")["targets"]]
            idx = np.asarray([S5.PANEL.index(t) for t in targets])
            cols = []
            for k in [s for s in strength if s != f]:
                v = strength[k][idx]
                r = crank(v)
                r[~np.isfinite(v)] = np.nan
                cols.append(r)
            m = np.nanmean(np.vstack(cols), 0)
            data[f] = (x, a, dd, np.where(np.isfinite(m), m, 0.0)[:, None],
                       float(np.mean([x[j].sum() / np.sqrt(a[j].sum() * dd.sum()) for j in range(3)])))
        held = {h: gain(fit(data, [g for g in CANCER if g != h], [0]), data[h], [0]) for h in CANCER}
        th = fit(data, list(CANCER), [0])
        thr_out[thr] = dict(held=held, kolf2=gain(th, data["kolf2"], [0]), theta=float(th[0]))
        print(f"threshold {thr:g}: LOFO cancer {100 * np.mean(list(held.values())):+.2f} % (" + " / ".join(f"{100 * v:+.2f}" for v in held.values())
              + f"), kolf2 {100 * thr_out[thr]['kolf2']:+.2f} %, theta {th[0]:+.2f}")
    (WORK / "robust518.json").write_text(json.dumps(dict(exact=out, threshold=thr_out), indent=1, default=float))


if __name__ == "__main__":
    main()
