#!/usr/bin/env python3
"""363 — paired six-member contrasts + tilt extrapolation + sign-tracking.

Reads only this experiment's work/ + 323's reference/anchors (read-only).
Writes results.json + per-target-363.csv into this experiment's work/.
"""
import json
import os
import sys

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np
import polars as pl

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
ROOT = "~/Code/virtual-cell-challenge-2026"
BASE = os.path.join(ROOT, "experiments/323-current-family-validation")
W323 = os.path.join(BASE, "work")
for _p in (HERE, BASE, os.path.join(ROOT, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import evaluate323 as E  # noqa: E402
import report323 as B  # noqa: E402
from cell_eval2.prep import bulk_lognorm_means  # noqa: E402

ARMS = ("cap15", "block_down", "block_random", "block_up")
DRAWS, SEED = 10000, 20260909
MEMBERS = list(E.MEMBERS)
assert MEMBERS == ["pds", "mse", "nmae", "fid", "reach", "jac"]


def resample_counts(n):
    rng = np.random.default_rng(SEED)
    draws = rng.integers(0, n, (DRAWS, n))
    counts = np.zeros((DRAWS, n), dtype=np.int64)
    np.add.at(counts, (np.repeat(np.arange(DRAWS), n), draws.ravel()), 1)
    return np.vstack([np.ones(n, dtype=np.int64), counts])


def main():
    anchors = json.load(open(os.path.join(W323, "h1-anchors.json")))
    base = np.array([anchors["baseline"][m] for m in MEMBERS])
    repl = np.array([anchors["replicate"][m] for m in MEMBERS])
    r = np.load(os.path.join(W323, "reference.npz"), allow_pickle=False)
    targets = r["targets"].astype(str)
    genes = r["genes"].astype(str)
    truth_all = bulk_lognorm_means(r["sums"].astype(float), E.TS)
    truth, control = truth_all[:-1], truth_all[-1]
    real_jk, control_jk = r["jk"][:-1], r["jk"][-1]
    lookup = {g: i for i, g in enumerate(genes)}
    own = np.array([lookup[t] for t in targets])
    n = len(targets)
    counts = resample_counts(n)

    raw, overall = {}, {}
    for name in ARMS:
        metrics = json.load(open(os.path.join(WORK, f"metrics-{name}.json")))
        z = np.load(os.path.join(WORK, f"stats-{name}.npz"), allow_pickle=False)
        assert np.array_equal(z["targets"].astype(str), targets)
        pm = bulk_lognorm_means(z["sums"].astype(float), E.TS)
        mse, budget, rho, den = B.mse_bootstrap(
            counts, pm, truth, control, z["jk"], real_jk, control_jk, own)
        assert abs(mse[0] - metrics["raw"]["mse"]) <= 1e-10, f"mse identity {name}"
        assert np.all(den > 0)
        table = pl.read_csv(os.path.join(WORK, f"per-target-{name}.csv"))
        mat = np.empty((len(counts), len(MEMBERS)))
        for j, m in enumerate(MEMBERS):
            if m == "mse":
                mat[:, j] = mse
                continue
            rows = table.filter(pl.col("metric") == E.MEMBERS[m])
            mp = dict(zip(rows["perturbation"], rows["value"]))
            vals = np.array([mp.get(str(t), np.nan) for t in targets], float)
            valid = np.isfinite(vals)
            mat[:, j] = (counts[:, valid] @ vals[valid]) / counts[:, valid].sum(1)
            assert abs(mat[0, j] - metrics["raw"][m]) <= 1e-12, f"mean identity {name}/{m}"
        raw[name] = mat
        overall[name] = B.scaled(mat, base, repl).mean(1)
        z.close()

    out_arms = {}
    for name in ARMS:
        metrics = json.load(open(os.path.join(WORK, f"metrics-{name}.json")))
        sc = (np.array([metrics["raw"][m] for m in MEMBERS]) - base) / (repl - base)
        sc[1] = min(max(sc[1], 0.0), 1.0)
        sc[2] = max(sc[2], -6.0)
        out_arms[name] = {"raw": metrics["raw"], "local_scaled": dict(zip(MEMBERS, map(float, sc))),
                          "local_h1_overall": float(sc.mean()),
                          "seconds": metrics["seconds"]}

    contrasts = {}
    for changed in ("block_down", "block_random", "block_up"):
        d_all = overall[changed] - overall["cap15"]
        entry = {"overall_point": float(d_all[0]), "overall_ci": B.ci(d_all[1:]), "raw": {}}
        for j, m in enumerate(MEMBERS):
            d = raw[changed][:, j] - raw["cap15"][:, j]
            e = {"point": float(d[0]), **B.ci(d[1:])}
            if m in ("mse", "nmae"):
                e["beneficial"] = bool(e["hi"] < 0)
            else:
                e["beneficial"] = bool(e["lo"] > 0)
            entry["raw"][m] = e
        contrasts[f"{changed}-minus-cap15"] = entry

    # ---- tilt extrapolation: per-target fid delta vs own block down-tilt ----
    s_cap = np.load(os.path.join(WORK, "stats-cap15.npz"))
    q_cap, l_cap = s_cap["q"], s_cap["lfc"]
    s_cap.close()
    rr = np.load(os.path.join(W323, "reference.npz"), allow_pickle=False)
    keep = rr["keep"]
    ugenes = rr["genes"].astype(str)[keep]
    r_q, r_lfc = rr["q"], rr["lfc"]
    ref_lin = rr["ref_lin"]
    rtargets = rr["targets"].astype(str)
    rr.close()
    bl = np.load(os.path.join(WORK, "blocks.npz"))
    block = bl["block"]
    bl.close()
    assert np.array_equal(rtargets, targets)
    ze = np.load(os.path.join(W323, "effects.npz"), allow_pickle=False)
    egenes = ze["genes"].astype(str)
    et = ze["targets"].astype(str)
    cap15 = ze["cap15"].astype(float)
    ze.close()
    pos = {g: i for i, g in enumerate(egenes)}
    epos = {t: i for i, t in enumerate(et)}
    D = np.zeros((len(targets), len(ugenes)))
    for j, g in enumerate(ugenes):
        if g in pos:
            D[:, j] = cap15[[epos[t] for t in targets], pos[g]]
    own_u = np.equal.outer(targets, ugenes)
    adjud = np.isfinite(r_lfc) & (r_lfc != 0)
    ref_sig = r_q < 0.05

    def fid_vec(qp, lp):
        ps = np.isfinite(qp) & (qp < 0.05)
        scored = ps & adjud & ~own_u
        n_pred = scored.sum(1).astype(float)
        n_real = (ref_sig & ~own_u).sum(1).astype(float)
        k = (scored & np.isfinite(lp) & (lp != 0)
             & (np.sign(lp) == np.sign(r_lfc))).sum(1).astype(float)
        den = np.maximum(n_pred, n_real)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(den > 0, k / den, np.nan), n_pred, n_real, k

    # tilt = frac of cap15 CALLED block genes with reference lfc < 0
    calls_cap = (np.isfinite(q_cap) & (q_cap < 0.05)) & ~own_u
    tilt = np.where((calls_cap & block).sum(1) > 0,
                    ((r_lfc < 0) & calls_cap & block).sum(1) / (calls_cap & block).sum(1), np.nan)
    F_cap, n_pred_cap, n_real_cap, _ = fid_vec(q_cap, l_cap)
    rows = {"target": [], "tilt_block": [], "n_pred_cap15": [], "n_real": []}
    for a in ("block_down", "block_random", "block_up"):
        rows[f"fid_{a}"] = []
        rows[f"dfid_{a}"] = []
    for a in ("block_down", "block_random", "block_up"):
        z = np.load(os.path.join(WORK, f"stats-{a}.npz"))
        F, n_pred, _, _ = fid_vec(z["q"], z["lfc"])
        z.close()
        rows[f"fid_{a}"].extend(map(float, F))
        rows[f"dfid_{a}"].extend(map(float, F - F_cap))
    rows["target"].extend(map(str, targets))
    rows["tilt_block"].extend(map(float, tilt))
    rows["n_pred_cap15"].extend(map(float, n_pred_cap))
    rows["n_real"].extend(map(float, n_real_cap))
    pl.DataFrame(rows).write_csv(os.path.join(WORK, "per-target-363.csv"))

    tilt_corr = {}
    for a in ("block_down", "block_random", "block_up"):
        d = np.array(rows[f"dfid_{a}"])
        m = np.isfinite(d) & np.isfinite(tilt)
        tilt_corr[a] = {"pearson_r": float(np.corrcoef(tilt[m], d[m])[0, 1]),
                        "slope_per_0.01_tilt": float(
                            np.polyfit(tilt[m], d[m], 1)[0] * 0.01),
                        "n": int(m.sum()),
                        "tilt_mean": float(np.nanmean(tilt)),
                        "tilt_median": float(np.nanmedian(tilt)),
                        "tilt_min": float(np.nanmin(tilt)),
                        "tilt_max": float(np.nanmax(tilt))}

    # ---- sign tracking: does delivered lfc follow the design sign on block genes? ----
    track = {}
    for a in ARMS:
        z = np.load(os.path.join(WORK, f"stats-{a}.npz"))
        lp = z["lfc"]
        z.close()
        # design matrix for this arm on universe
        if a == "cap15":
            DD = D
        else:
            DD = D.copy()
            sgn = {"block_down": -1.0, "block_up": 1.0}[a] if a != "block_random" else None
            if sgn is not None:
                DD[block] = sgn * np.abs(D[block])
            else:
                rng = np.random.default_rng(20260910)
                S = rng.choice(np.array([-1.0, 1.0]), size=D.shape)
                DD[block] = S[block] * np.abs(D[block])
        agree = ((np.sign(lp[block]) == np.sign(DD[block]))
                 & (lp[block] != 0) & (DD[block] != 0))
        denom = ((lp[block] != 0) & (DD[block] != 0)).sum()
        track[a] = {"delivered_eq_design_sign_frac": float(agree.sum() / denom),
                    "denom": int(denom)}

    # ---- call-set stability ----
    call = {}
    for a in ARMS:
        z = np.load(os.path.join(WORK, f"stats-{a}.npz"))
        ps = (np.isfinite(z["q"]) & (z["q"] < 0.05) & ~own_u)
        z.close()
        call[a] = {"n_pred_mean": float(ps.sum(1).mean()),
                   "n_pred_median": float(np.median(ps.sum(1).astype(float)))}
    zc = np.load(os.path.join(WORK, "stats-cap15.npz"))
    pc = (np.isfinite(zc["q"]) & (zc["q"] < 0.05) & ~own_u)
    zc.close()
    for a in ("block_down", "block_random", "block_up"):
        z = np.load(os.path.join(WORK, f"stats-{a}.npz"))
        ps = (np.isfinite(z["q"]) & (z["q"] < 0.05) & ~own_u)
        z.close()
        inter = (ps & pc).sum(1).astype(float)
        union = (ps | pc).sum(1).astype(float)
        call[a]["jaccard_vs_cap15_mean"] = float(np.mean(inter / union))

    result = {"arms": out_arms, "contrasts": contrasts, "tilt": tilt_corr,
              "sign_tracking": track, "call_set": call,
              "anchors": "experiments/323-current-family-validation/work/h1-anchors.json",
              "bootstrap": {"draws": DRAWS, "seed": SEED, "unit": "target", "paired": True},
              "scope": "85-target reused H1 development panel; not an A/B/C forecast"}
    with open(os.path.join(WORK, "results.json"), "w") as f:
        json.dump(result, f, indent=2, sort_keys=True)
    print(json.dumps({"arms": {k: v["local_h1_overall"] for k, v in out_arms.items()},
                      "contrasts": {k: {"overall_point": v["overall_point"],
                                        "overall_ci": v["overall_ci"]} for k, v in contrasts.items()},
                      "tilt": tilt_corr, "sign_tracking": track, "call_set": call}, indent=2))


if __name__ == "__main__":
    main()
