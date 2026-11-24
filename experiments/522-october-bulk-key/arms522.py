"""522 — the October key forms, exact, on every surface: the A/B/C PLOLO folds (518's dumps; k562, hct116, hek293t,
kolf2) and 438's October-like panels' cancer folds (522's dumps; efp, efr, real x k562, hct116, hek293t).

Every arm is a per-target weight w on the surface's unkeyed composition (key_scale restores the non-panel Frobenius norm),
then q 1.75, render, ballast KMAX 16 (fold522). The base is the LODO key at clip (0.5, 3) (446bj2 on A/B/C; the October
default `lodo_hi_pct` off). Families:
  cap H          clip(key/median, 0.5, H), H in 4 ... 16 and uncapped
  pct P          clip(key/median, 0.5, hi(P)), hi(P) = the P-th percentile of carried key/median (543; P_PROD =
                 96.64697522769474 is the candidate October form, == cap 6 on the A/B/C artifact)
  S th           518's strength key, src form, on the cap-3 key: w_cap3 x exp(th z), z = zsrc (bounded: z in [-.5, .5])
  S+cd4 th       the same with CD4's robust strength added to the rank mean (October-portable: CD4 is extracted per panel)
  Sonly th       exp(th z) alone, no LODO key (518's substitute test)
  cap6 x S th    518's S2 construction on the src form;   pct x S th: the same on the percentile cap
  gate           cap 6, and the percentile cap only where z > 0 (the tail allowed where the donors' own strength agrees)
  CTRL           S th 1 with z permuted over targets (seed 522): must not help
KNOWN ANSWERS (abort): each surface's load KA (fold522); Oct: fixed cap 6 == 540o's `V1_0.5_6` row (L.cmp < 1e-12);
Abc: zsrc == skey518.zfeat(F, 'src') exactly and `src|theta1|on_cap3` == skey518.json (r and dpds, |d| < 1e-12).

    uv run python experiments/522-october-bulk-key/arms522.py <panel> [folds...]      # ~1-2 GB, a few min per fold
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import fold522 as FD  # noqa: E402
import lib438 as L  # noqa: E402

WORK = HERE / "work"
W540 = ROOT / "experiments/540-lodo-clip-scan/work"
W518 = ROOT / "experiments/518-bulk-mse-gap/work"
TOL = 1e-12
CAPS = (4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 16.0, 1e9)
PCTS = (90.0, 92.5, 95.0, 96.0, FD.P_546, FD.P_PROD, 97.5, 98.0, 99.0)
THS = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
COMBO_THS = (0.25, 0.5, 0.75, 1.0)
KEEP_ARRAYS = ("cap3", "cap6", "pctPROD", "cap1e+09", "S1", "S0.5", "cap6xS0.5", "pctPRODxS0.5", "pctPRODxS0.25", "gate")


def pname(P):
    """arm name of a percentile cap: 'pctPROD' for 497's unrounded P, 'pct96.647' for 546's rounded one."""
    return "pctPROD" if P == FD.P_PROD else f"pct{P:g}"


def arms_for(S, z, zc):
    kv, car = S.kv, S.carried
    w3 = FD.wclip(kv, car, 0.5, 3.0)
    w6 = FD.wclip(kv, car, 0.5, 6.0)
    hp = FD.pct_hi(kv, car, FD.P_PROD)
    wp = FD.wclip(kv, car, 0.5, hp)
    A = {"cap3": w3}
    for H in CAPS:
        A[f"cap{H:g}"] = FD.wclip(kv, car, 0.5, H)
    for P in PCTS:
        A[pname(P)] = FD.wclip(kv, car, 0.5, FD.pct_hi(kv, car, P))
    for th in THS:
        A[f"S{th:g}"] = w3 * np.exp(th * z)
    for th in (0.5, 1.0, 1.5):
        A[f"Scd4_{th:g}"] = w3 * np.exp(th * zc)
    for th in (1.0, 2.0):
        A[f"Sonly{th:g}"] = np.exp(th * z)
    for th in COMBO_THS:
        A[f"cap6xS{th:g}"] = w6 * np.exp(th * z)
        A[f"pctPRODxS{th:g}"] = wp * np.exp(th * z)
    A["gate"] = np.where(z > 0, wp, w6)
    zp = np.random.default_rng(522).permutation(z)
    A["CTRL S1 permuted"] = w3 * np.exp(1.0 * zp)
    return A, dict(hi_pct_prod=hp, hi_pct_546=FD.pct_hi(kv, car, FD.P_546),
                   his={pname(P): FD.pct_hi(kv, car, P) for P in PCTS})


def run(panel: str, held: str) -> None:
    t0 = time.time()
    log = lambda *a: print(f"[{panel}/{held}]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    S = FD.load(panel, held)
    log(f"loaded; KA on load {S.ka}")
    z = FD.zsrc(S)
    zc = FD.zsrc(S, extra=("cd4_robust",))
    ka = dict(load=S.ka)
    if panel == "abc":
        import skey518
        z518 = skey518.zfeat(S.F, "src")
        assert np.array_equal(z, z518), "KA zsrc != skey518.zfeat"
        ka["zsrc == skey518.zfeat"] = "PASS (bit-identical)"
    A, his = arms_for(S, z, zc)
    res = dict(panel=panel, fold=held, n=S.n, n_carried=int(S.carried.sum()), sources=S.sources, his=his, ka=ka, arms={},
               z_quantiles=np.quantile(z, [0, .1, .5, .9, 1]).tolist(), frac_z_zero=float(np.mean(z == 0)))
    arrays = {}
    base = None
    for nm, w in A.items():
        s, arr = S.score(w)
        if base is None:
            base = s
        rec = dict(cos_agg=s["cos_agg"], pds=s["pds"], pds_lam05=s["pds_lam0.5"], K_mean=s["K_mean"],
                   r=s["cos_agg"] / base["cos_agg"] - 1, dpds=s["pds"] - base["pds"],
                   dpds_bd=s["pds_lam0.5"] - base["pds_lam0.5"], w_max=float(np.max(w)),
                   w_q=np.quantile(w, [0, .5, .9, .99, 1]).tolist(),
                   AB={c: s["AB"][c] for c in ("A", "B", "C")})
        res["arms"][nm] = rec
        arrays[f"w__{nm}"] = np.asarray(w, float)
        if nm in KEEP_ARRAYS:
            for k in ("x", "a", "pds", "d"):
                arrays[f"{nm}__{k}"] = arr[k]
        log(f"{nm:22s} cos {s['cos_agg']:.6f} r {100 * rec['r']:+7.3f} % dpds {rec['dpds']:+.5f} K {s['K_mean']:.2f} wmax {rec['w_max']:.2f}")
    # ---- known answers that need the finished rows
    if panel == "abc":
        ref = json.loads((W518 / "skey518.json").read_text())[held]["arms"]["src|theta1|on_cap3"]
        d1 = abs(res["arms"]["S1"]["r"] - ref["r_vs_own_base"])
        d2 = abs(res["arms"]["S1"]["dpds"] - ref["dpds"])
        log(f"KA S1 (src th1 on cap3) vs skey518: r {100 * res['arms']['S1']['r']:+.6f} % vs {100 * ref['r_vs_own_base']:+.6f} %; |d| {d1:.1e} / {d2:.1e}")
        assert d1 < TOL and d2 < TOL, "KA FAILED skey518 src th1"
        ka["S1 == skey518 src|theta1|on_cap3"] = max(d1, d2)
        r540 = json.loads((W540 / f"arms540-{held}.json").read_text())["arms"]
        d6 = abs(res["arms"]["cap6"]["cos_agg"] - r540["V1_0.5_6"]["realised"]["16"]["cos_agg"])
        assert d6 < TOL, "KA cap6 vs arms540"
        ka["cap6 == arms540 V1_0.5_6 cos_agg"] = d6
    else:
        ref = json.loads((W540 / f"comp540o-{panel}-{held}.json").read_text())["arms"]["V1_0.5_6"]["realised"]["16"]
        s6 = {k: res["arms"]["cap6"][k] for k in ("cos_agg", "pds")}
        d6 = max(abs(s6["cos_agg"] - ref["cos_agg"]), abs(s6["pds"] - ref["pds"]))
        assert d6 < TOL, "KA cap6 vs 540o"
        ka["cap6 == 540o V1_0.5_6 (cos_agg, pds)"] = d6
        ref546 = json.loads((FD.W546 / f"comp546o-{panel}-{held}.json").read_text())["arms"]["p96.6"]["realised"]["16"]
        dp = max(abs(res["arms"][f"pct{FD.P_546:g}"]["cos_agg"] - ref546["cos_agg"]), abs(res["arms"][f"pct{FD.P_546:g}"]["pds"] - ref546["pds"]))
        assert dp < TOL, "KA pct96.647 vs 546o"
        ka["pct96.647 == 546o p96.6 (cos_agg, pds)"] = dp
    log(f"KA {ka}")
    res["peak_rss_gib"] = L.E062.peak_rss_gib()
    res["seconds"] = time.time() - t0
    (WORK / f"arms522-{panel}-{held}.json").write_text(json.dumps(res, indent=1, default=float))
    np.savez_compressed(WORK / f"arms522-{panel}-{held}.npz", z=z, zc=zc, **arrays)
    log(f"done; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    panel = sys.argv[1]
    folds = sys.argv[2:] or [f for p, f in FD.surfaces((panel,))]
    for f in folds:
        run(panel, f)
