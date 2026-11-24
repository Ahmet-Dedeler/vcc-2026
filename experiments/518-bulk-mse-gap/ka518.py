"""518 — known answers, run before any decomposition number is read (abort on failure).

KA1  187's cliff identity (analyse418.Cliff) at 446bj2's measured energies reproduces the board raw mse 0.6264081539046438
     at its solved X (|d| < 1e-12).
KA2  540's recorded cap-6 price: re-deriving `V1_0.5_6` from the saved inputs (abc540.json, arms540-<fold>.json,
     analyse492g.json, ballast492g.json) through 540's arithmetic reproduces analyse540.json's cancer r, dmse at k_X 1.31 and
     1.0, and the net overall (|d| < 1e-12). Nothing is written to 540's directory.
KA3  420's saved per-target arrays (x = <pred, truth>, a = |pred|^2 per frame, d = |truth|^2) rebuild comp420's realised-16
     cos_agg for `stack|V1_5_3+sh1.75` (446bj2's bulk, unlevelled) and comp420L's levelled row
     `stack|V1_5_3+sh1.75:A=1.101,B=1.155,C=1.183` on all four folds (|d| < 1e-12).
KA4  H1: 096's scorer-exact msekit on 435's harness stats (410b-construction bulk, coarse MDE, seed 0) and 323's reference
     reproduces the harness's recorded raw mse 0.9253584049461503 (|d| < 1e-12).

    uv run python experiments/518-bulk-mse-gap/ka518.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/418-bulk-composition"))
sys.path.insert(0, str(ROOT / "experiments/492-per-target-sharpening"))
W420 = ROOT / "experiments/420-bulk-stack/work"
W540 = ROOT / "experiments/540-lodo-clip-scan/work"
W492 = ROOT / "experiments/492-per-target-sharpening/work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
CANCER = FOLDS[:3]
TOL = 1e-12


def cos_agg_from_arrays(z, name):
    x, a, d = z[f"{name}__x"], z[f"{name}__a"], z["d"]
    return float(np.mean([x[j].sum() / np.sqrt(a[j].sum() * d.sum()) for j in range(3)]))


def ka1_cliff():
    import analyse418 as A8
    from analyse492g import RAW_446BJ2
    cl = A8.Cliff()
    g = json.loads((W492 / "analyse492g.json").read_text())
    X0 = cl.solve_X(g["AB_446bj2"], RAW_446BJ2)
    d = abs(cl.raw(g["AB_446bj2"], X0) - RAW_446BJ2)
    assert d < TOL and abs(X0 - g["X_446bj2"]) < TOL, f"KA1 FAILED {d}"
    return dict(raw=cl.raw(g["AB_446bj2"], X0), board=RAW_446BJ2, absdiff=d, X=X0)


def ka2_540():
    """analyse540.main()'s arithmetic for the cap-6 arm only, reading the same inputs, writing nothing there."""
    import analyse418 as A8
    from analyse492g import CTX, RAW_446BJ2, level_opt
    cl = A8.Cliff()
    g = json.loads((W492 / "analyse492g.json").read_text())
    bal = json.loads((W492 / "ballast492g.json").read_text())
    slope_K, slope_D = bal["fits"]["dK_mean"]["slope"], bal["fits"]["D_from_design"]["slope"]
    AB0 = g["AB_446bj2"]
    X0 = cl.solve_X(AB0, RAW_446BJ2)
    abc = json.loads((W540 / "abc540.json").read_text())
    fs = [0.92, 0.95, 0.98, 1.0, 1.02, 1.05, 1.08]

    def ratio_map(tag, key):
        mp = abc["maps"][tag]
        return {c: [mp[f"{f:g}"][c][key] / mp["1"][c][key] for f in fs] for c in CTX}

    base, tg = "V1_0.5_3", "V1_0.5_6"
    f0, raw0 = level_opt(cl, AB0, X0, fs, ratio_map(base, "A_own"), ratio_map(base, "B_own"))
    R = {f: json.loads((W540 / f"arms540-{f}.json").read_text()) for f in FOLDS}
    a = abc["arms"][tg]
    ABa = {c: dict(A=AB0[c]["A"] * a[c]["A_ratio"], B=AB0[c]["B"] * a[c]["B_ratio"]) for c in CTX}
    amp = np.sqrt(np.mean([ABa[c]["A"] for c in CTX]) / np.mean([AB0[c]["A"] for c in CTX]))
    ram, rbm = ratio_map(tg, "A_own"), ratio_map(tg, "B_own")
    rows = {}
    for f in CANCER:
        b, x = R[f]["arms"][base]["realised"]["16"], R[f]["arms"][tg]["realised"]["16"]
        r = x["cos_agg"] / b["cos_agg"] - 1
        row = dict(r=r, dpds_lit=x["pds"] - b["pds"])
        for k in (1.0, 1.31):
            row[k] = level_opt(cl, ABa, X0 * (1 + k * r) * amp, fs, ram, rbm)[1] - raw0
        rows[f] = row
    r = float(np.mean([rows[f]["r"] for f in CANCER]))
    dm131 = float(np.mean([rows[f][1.31] for f in CANCER]))
    dm1 = float(np.mean([rows[f][1.0] for f in CANCER]))
    dp = 0.34 * float(np.mean([rows[f]["dpds_lit"] for f in CANCER]))
    net_K = (A8.PDS_SC * dp + A8.MSE_SC * dm131) / 6 + (-0.2808) * slope_K * a["dK_mean"]
    rec = json.loads((W540 / "analyse540.json").read_text())["arms"][tg]
    want = dict(r=rec["cancer"]["r"], dm131=rec["cancer"]["dmse_k1.310"], dm1=rec["cancer"]["dmse_k1.000"],
                net_K=rec["calibrated"]["k_X 1.31"]["net_K"])
    got = dict(r=r, dm131=dm131, dm1=dm1, net_K=net_K)
    d = max(abs(got[k] - want[k]) for k in want)
    assert d < TOL, f"KA2 FAILED {got} vs {want}"
    return dict(got=got, recorded=want, absdiff=d, base_level=f0, base_raw=raw0)


def ka3_arrays():
    out = {}
    for f in FOLDS:
        z = np.load(W420 / f"arrays420-{f}.npz")
        c = json.loads((W420 / f"comp420-{f}.json").read_text())["arms"]["stack|V1_5_3+sh1.75"]["realised"]["16"]["cos_agg"]
        zl = np.load(W420 / f"arrays420L-{f}.npz")
        lv = "stack|V1_5_3+sh1.75:A=1.101,B=1.155,C=1.183"
        cL = json.loads((W420 / f"comp420L-{f}.json").read_text())["arms"][lv]["realised"]["16"]["cos_agg"]
        g1, g2 = cos_agg_from_arrays(z, "stack|V1_5_3+sh1.75@R16"), cos_agg_from_arrays(zl, lv + "@R16")
        d = max(abs(g1 - c), abs(g2 - cL))
        assert d < TOL, f"KA3 FAILED {f}: {g1} vs {c}; {g2} vs {cL}"
        out[f] = dict(unlevelled=g1, levelled=g2, absdiff=d)
    return out


def ka4_h1():
    sys.path.insert(0, str(ROOT / "experiments/096-mse-floor"))
    import msekit
    z = np.load(ROOT / "experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz")
    r = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz")
    pred_sums = np.vstack([z["sums"], r["sums"][-1]])
    pred_jk = np.append(z["jk"], r["jk"][-1])
    assert list(z["targets"]) == list(r["targets"]) and str(r["perts"][-1]) == "non-targeting"
    mse, num, den = msekit.score(r["perts"], pred_sums, r["counts"], pred_jk, r["perts"], r["sums"], r["counts"], r["jk"], r["genes"])
    want = json.loads((ROOT / "experiments/435-zs-on-coarse/work/h1/metrics-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.json").read_text())["raw"]["mse"]
    d = abs(mse - want)
    assert d < TOL, f"KA4 FAILED {mse} vs {want}"
    return dict(mse=mse, recorded=want, absdiff=d)


def main():
    res = dict(KA1=ka1_cliff(), KA2=ka2_540(), KA3=ka3_arrays(), KA4=ka4_h1())
    for k, v in res.items():
        print(k, json.dumps(v, default=float)[:400])
    print("ALL KNOWN ANSWERS PASS")
    (WORK / "ka518.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
