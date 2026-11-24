"""473c — price every 473 arm vs 446b's bulk in overall units (preregistered in PREREGISTRATION.md).

Energy (mse's A_c/B_c): 446b's MEASURED render energies (ntclevel-446b.json, the artifact) x the arm/base realised-route
ratio on the real A/B/C bulk (abc473.json, `A_own`/`B_own`: every gene but the row's own target, ntclevel's convention).
X (the cross term): solved from 446b's board raw mse (0.6408807601) with 187's cliff model (analyse418.Cliff); the arm's
X = X_base (1 + k r) sqrt(mean A_arm / mean A_base), r = the arm's relative cos_agg change on PLOLO (realised, KMAX 16),
k in {1, 1.407}. pds: the board-emulated pds (414's beta, lam 0.5) change on PLOLO (realised).
Re-level: base and arm each at their own EV-optimal per-context factor (420's level_opt: 187's cliff sd, the realised
energy map from abc473's `level_maps`), so the price is what the production re-level would deliver.
overall = (2.2138 dpds_bd - 1.051 dmse_raw) / 6.  Target bootstrap: 4,000 draws per fold on (r, dpds_bd), energies and
factors held at the point estimate.

    uv run python experiments/473-410d2-reprice/analyse473.py
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
sys.path.insert(0, str(ROOT / "experiments/418-bulk-composition"))
import analyse418 as A8  # noqa: E402

CTX = ("A", "B", "C")
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
FOLDS = (*CANCER, "kolf2")
BASE = "stack|K3+sh1.75"
ARMS_PLOLO = {"A1": "A1 clip2 pre rg", "A2": "A2 unclipped pre rg", "A3": "A3 clip3 pre rg",
              "A1post": "A1post clip2 post-sharpen", "C1": "C1 clip2 permuted"}
ARMS_ABC = {"A1": "A1 m clip[1/2,2] pre-sharpen rg", "A2": "A2 m unclipped pre-sharpen rg",
            "A3": "A3 m clip[1/3,3] pre-sharpen rg", "A1post": "A1post m clip[1/2,2] post-sharpen row-renorm",
            "C1": "C1 A1 with m permuted over genes"}
MAP_OF = {"A1": "A1", "A2": "A1", "A3": "A1", "A1post": "A1post", "C1": "A1"}   # energy-map shape used for the re-level
RAW_446B = 0.6408807601228349          # data/board-entries/UDgjOKL9gI1udAexq6ue.json
SD187 = {"A": 0.0025, "B": 0.0004, "C": 0.0013}
KS = (1.0, A8.K_KEY)
NBOOT = 4000


def level_opt(cl, AB, X, fs, ra_map, rb_map, seed=473):
    """EV-optimal per-context factor f on the realised energy map (420's level_opt arithmetic)."""
    grid = np.round(np.arange(max(0.92, min(fs)), min(1.08, max(fs)) + 1e-9, 0.001), 3)
    rng = np.random.default_rng(seed)
    best, raw_parts = {}, {}
    for c in CTX:
        ra = np.interp(grid, fs, ra_map[c])
        rb = np.interp(grid, fs, rb_map[c])
        A, B, Kc, claim = AB[c]["A"], AB[c]["B"], cl.Kc[c], cl.claim[c]
        draws = claim * np.exp(2 * rng.normal(0, SD187[c], 400))
        ev = np.array([np.mean(A * ra[i] + Kc - np.minimum(draws, B * rb[i]) - 2 * X * np.sqrt(ra[i])) for i in range(grid.size)])
        i = int(np.argmin(ev))
        best[c] = float(grid[i])
        raw_parts[c] = float(1.0 + ev[i] / cl.sumD)                 # EV mse_c at the optimum: mse_c = 1 + [..]/sumD
    return best, float(np.mean(list(raw_parts.values())))


def main():
    cl = A8.Cliff()
    ntc = json.loads((ROOT / "experiments/446-q175-level/work/ntclevel-446b.json").read_text())["per_context"]
    AB0 = {c: dict(A=ntc[c]["A"], B=ntc[c]["B"]) for c in CTX}
    X0 = cl.solve_X(AB0, RAW_446B)
    chk = cl.raw(AB0, X0)
    assert abs(chk - RAW_446B) < 1e-12, "cliff model does not reproduce 446b's raw at its solved X"
    abc = json.loads((WORK / "abc473.json").read_text())
    maps = abc["level_maps"]
    fs = [float(f) for f in maps["factors"]]

    def ratio_map(tag, key):
        mp = maps["maps"][tag]
        return {c: [mp[f"{f:g}"][c]["realised"][key] / mp["1"][c]["realised"][key] for f in fs] for c in CTX}

    ra0, rb0 = ratio_map("base", "A_own"), ratio_map("base", "B_own")
    f0, raw0 = level_opt(cl, AB0, X0, fs, ra0, rb0)
    out = dict(X_446b=X0, base_level=f0, base_raw_levelled_ev=raw0, arms={})
    print(f"446b: X {X0:.6f}; re-level factors {f0} (1.0 = as shipped); EV raw at optimum {raw0:.5f} (board {RAW_446B:.5f})")

    R = {}
    for f in FOLDS:
        p = WORK / f"arms473-{f}.json"
        if p.exists():
            R[f] = json.loads(p.read_text())
    print(f"folds present: {list(R)}")
    Z = {f: np.load(WORK / f"arrays473-{f}.npz") for f in R}

    for arm, pname in ARMS_PLOLO.items():
        e_arm = abc["arms"][ARMS_ABC[arm]]["E"]
        e_b = abc["arms"][ARMS_ABC[arm]]["E_base"]
        er = {c: dict(A=e_arm[c]["realised"]["A_own"] / e_b[c]["realised"]["A_own"],
                      B=e_arm[c]["realised"]["B_own"] / e_b[c]["realised"]["B_own"]) for c in CTX}
        ABa = {c: dict(A=AB0[c]["A"] * er[c]["A"], B=AB0[c]["B"] * er[c]["B"]) for c in CTX}
        amp = np.sqrt(np.mean([ABa[c]["A"] for c in CTX]) / np.mean([AB0[c]["A"] for c in CTX]))
        ram, rbm = ratio_map(MAP_OF[arm], "A_own"), ratio_map(MAP_OF[arm], "B_own")
        res = dict(energy_ratio={c: er[c]["A"] for c in CTX}, folds={})
        for f, Rf in R.items():
            b, a = Rf["arms"][BASE]["realised"]["16"], Rf["arms"][pname]["realised"]["16"]
            r = a["cos_agg"] / b["cos_agg"] - 1
            dpb, dpl = a["pds_lam0.5"] - b["pds_lam0.5"], a["pds"] - b["pds"]
            row = dict(r=r, dpds_bd=dpb, dpds_lit=dpl)
            for k in KS:
                Xa = X0 * (1 + k * r) * amp
                fa, rawa = level_opt(cl, ABa, Xa, fs, ram, rbm)
                dm = rawa - raw0
                row[f"k{k:.3f}"] = dict(dmse=dm, overall=(A8.PDS_SC * dpb + A8.MSE_SC * dm) / 6, level=fa)
            # target bootstrap on (r, dpds_bd) with energies/factors at the point estimate, k = 1
            z = Z[f]
            d = z["d"]
            xb, ab, pb = z[f"{BASE}@R16__x"], z[f"{BASE}@R16__a"], z[f"{BASE}@R16__pdsl"]
            xa, aa, pa = z[f"{pname}@R16__x"], z[f"{pname}@R16__a"], z[f"{pname}@R16__pdsl"]
            rng = np.random.default_rng(4730 + FOLDS.index(f))
            n = d.size
            bo_r, bo_p = np.empty(NBOOT), np.empty(NBOOT)
            for i in range(NBOOT):
                ix = rng.integers(0, n, n)
                dd = d[ix].sum()
                cb = np.mean([xb[j, ix].sum() / np.sqrt(ab[j, ix].sum() * dd) for j in range(3)])
                ca = np.mean([xa[j, ix].sum() / np.sqrt(aa[j, ix].sum() * dd) for j in range(3)])
                bo_r[i] = ca / cb - 1
                bo_p[i] = float(np.mean(pa[:, ix] - pb[:, ix]))
            k1 = row["k1.000"]
            # linearise overall in (r, dpds) around the point estimate for the bootstrap (mse is linear in X at fixed f)
            dmse_dr = -2 * X0 * amp * np.mean([np.sqrt(np.interp(k1["level"][c], fs, ram[c])) for c in CTX]) / cl.sumD
            ov = (A8.PDS_SC * bo_p + A8.MSE_SC * (k1["dmse"] + dmse_dr * (bo_r - r))) / 6
            row["boot_k1"] = dict(r_ci=[float(np.quantile(bo_r, q)) for q in (0.025, 0.5, 0.975)],
                                  dpds_ci=[float(np.quantile(bo_p, q)) for q in (0.025, 0.5, 0.975)],
                                  overall_ci=[float(np.quantile(ov, q)) for q in (0.025, 0.5, 0.975)],
                                  P_overall_le0=float(np.mean(ov <= 0)), r_draws=bo_r, ov_draws=ov)
            res["folds"][f] = row
        canc = [f for f in CANCER if f in res["folds"]]
        if canc:
            for k in KS:
                kk = f"k{k:.3f}"
                res[f"cancer_mean_{kk}"] = float(np.mean([res["folds"][f][kk]["overall"] for f in canc]))
                res[f"cancer_pos_{kk}"] = int(sum(res["folds"][f][kk]["overall"] > 0 for f in canc))
            pooled_r = np.mean([res["folds"][f]["boot_k1"]["r_draws"] for f in canc], axis=0)
            pooled_ov = np.mean([res["folds"][f]["boot_k1"]["ov_draws"] for f in canc], axis=0)
            res["cancer_pooled_boot_k1"] = dict(r_ci=[float(np.quantile(pooled_r, q)) for q in (0.025, 0.5, 0.975)],
                                                overall_ci=[float(np.quantile(pooled_ov, q)) for q in (0.025, 0.5, 0.975)])
        for f in res["folds"]:
            for key in ("r_draws", "ov_draws"):
                res["folds"][f]["boot_k1"].pop(key)
        out["arms"][arm] = res
        line = f"{arm:7s} E {[round(er[c]['A'], 4) for c in CTX]}"
        for f in res["folds"]:
            rr = res["folds"][f]
            line += (f" | {f[:8]} r {rr['r']:+.4f} dpds {rr['dpds_bd']:+.5f} ov1 {rr['k1.000']['overall']:+.5f} "
                     f"CI {np.round(rr['boot_k1']['overall_ci'], 5).tolist()}")
        if canc:
            line += (f" || cancer mean k1 {res['cancer_mean_k1.000']:+.5f} ({res['cancer_pos_k1.000']}/{len(canc)} +), "
                     f"k1.407 {res['cancer_mean_k1.407']:+.5f}; pooled r CI {np.round(res['cancer_pooled_boot_k1']['r_ci'], 4).tolist()}")
        print(line, flush=True)

    # ---------------- the preregistered GO rule
    a1 = out["arms"].get("A1", {})
    go = {}
    abcA1 = abc["arms"][ARMS_ABC["A1"]]["go_item1"]
    go["item1_render_safety"] = abcA1
    if a1.get("folds"):
        canc = [f for f in CANCER if f in a1["folds"]]
        go["item2_folds_pos_k1"] = a1.get("cancer_pos_k1.000")
        go["item2_pooled_r_lower"] = a1["cancer_pooled_boot_k1"]["r_ci"][0] if canc else None
        go["item2"] = bool(canc and a1["cancer_pos_k1.000"] >= 2 and a1["cancer_pooled_boot_k1"]["r_ci"][0] > 0)
        go["item3_cancer_k1"] = a1.get("cancer_mean_k1.000")
        go["item3_kolf2_k1"] = a1["folds"].get("kolf2", {}).get("k1.000", {}).get("overall")
        go["item3"] = bool(canc and a1["cancer_mean_k1.000"] >= 0.0015 and (go["item3_kolf2_k1"] is not None and go["item3_kolf2_k1"] >= 0))
        c1 = out["arms"].get("C1", {})
        go["item4_c1_neg"] = int(sum(c1["folds"][f]["k1.000"]["overall"] < 0 for f in canc)) if c1.get("folds") else None
        go["item4"] = bool(go["item4_c1_neg"] is not None and go["item4_c1_neg"] >= 2)
        go["folds_present"] = list(a1["folds"])
    out["go"] = go
    print("GO rule:", json.dumps(go, default=float))
    (WORK / "analyse473.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
