"""492g pricer — 479's arms (P2g, S_P2, …) priced vs 446bj2 (the LODO-key stack), board-calibrated, for 479's preregistration.

analyse492.py's arithmetic with the anchor moved from 446b to 446bj2:
  energy  A_c/B_c of 446bj2 = ntclevel-446bj.json's MEASURED values (the LODO dump at 446b's factors) x this route's realised
          ratio 446bj2/446bj (abc492g `bj2_over_bj`); an arm's = 446bj2's x the arm/base realised ratio (abc492g).
  X       solved from 446bj2's board raw mse (0.6264081539) with 187's cliff model (analyse418.Cliff);
          X_arm = X_0 (1 + k r) sqrt(mean A_arm / mean A_0), r = the arm's relative cos_agg change on PLOLO (arms492g,
          realised KMAX 16) vs `stack|V1_5_3+sh1.75`.
  level   base and arm each at their EV-optimal per-context factor on abc492g's realised energy maps (420's level_opt).
  pds     the board-emulated pds change (414's beta, lam 0.5) on PLOLO; the literal pds change is reported beside it.
  price   overall = (2.2138 dpds - 1.051 dmse_raw) / 6 at k 1 / 1.407 (the instrument), then BOARD-CALIBRATED:
          k_pds 0.361 (both levelled q steps 433h->433k->446b; calib492) on dpds_bd, k_X 0.61 (the step adjacent to 446b)
          and 0.12 (pessimistic) on X, minus the ballast's nmae cost (ballast492g: the displacement model fitted on the
          three MDE-era bulk-only board steps; the K-count model beside it), -0.2808 overall per unit raw nmae.
Also: the same instrument retrodicts the board's key step 446b -> 446bj2 from comp420's rows (no fit), so the anchor's own
transfer is on the page.

    uv run python experiments/492-per-target-sharpening/analyse492g.py
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
BASE = "stack|V1_5_3+sh1.75"
RAW_446BJ2 = 0.6264081539046438        # data/board-entries/YvVru3DAF89Jiyr34olR.json
RAW_446B = 0.6408807601228349          # data/board-entries/UDgjOKL9gI1udAexq6ue.json
BOARD_KEY = dict(dpds=0.8674284652545521 - 0.8672203641768859, dmse=RAW_446BJ2 - RAW_446B, dnmae=0.000060)
SD187 = {"A": 0.0025, "B": 0.0004, "C": 0.0013}
K_PDS, K_PDS_LIT, K_X, K_X_PESS = 0.3612872636253093, 0.34, 0.61, 0.12     # calib492.json (lit: mean of 0.346 / 0.334)
KS = (K_X_PESS, K_X, 1.0, A8.K_KEY)
NMAE_OV = -0.2808                       # overall per unit raw nmae (446b's board record)
NBOOT = 4000
# PLOLO arm -> (abc492g arm for the energy ratio, abc492g map for the re-level)
ARMS = {"P2g": ("P2g", "P2g"), "S_P2": ("S_P2", "S_P2"), "P2 (no guard)": ("P2 (no guard)", "P2g"),
        "L30g": ("L30g", "L30g"), "P2g_K3key": ("P2g_K3key", "P2g"), "CTRL P2g permuted": ("CTRL P2g permuted", "G2"),
        "stack|V1_5_3+sh2": ("G2", "G2")}


def level_opt(cl, AB, X, fs, ra_map, rb_map, seed=473):
    """EV-optimal per-context factor f on the realised energy map (420's level_opt arithmetic, as analyse473/492)."""
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
        raw_parts[c] = float(1.0 + ev[i] / cl.sumD)
    return best, float(np.mean(list(raw_parts.values())))


def main():
    cl = A8.Cliff()
    abc = json.loads((WORK / "abc492g.json").read_text())
    bal = json.loads((WORK / "ballast492g.json").read_text())
    ntj = json.loads((ROOT / "experiments/446-q175-level/work/ntclevel-446bj.json").read_text())["per_context"]
    ntb = json.loads((ROOT / "experiments/446-q175-level/work/ntclevel-446b.json").read_text())["per_context"]
    AB0 = {c: dict(A=ntj[c]["A"] * abc["bj2_over_bj"][c]["A_own"], B=ntj[c]["B"] * abc["bj2_over_bj"][c]["B_own"]) for c in CTX}
    X0 = cl.solve_X(AB0, RAW_446BJ2)
    assert abs(cl.raw(AB0, X0) - RAW_446BJ2) < 1e-12, "cliff model does not reproduce 446bj2's raw at its solved X"
    ABb = {c: dict(A=ntb[c]["A"], B=ntb[c]["B"]) for c in CTX}
    Xb = cl.solve_X(ABb, RAW_446B)
    maps = abc["level_maps"]
    fs = [float(f) for f in maps["factors"]]

    def ratio_map(tag, key):
        mp = maps["maps"][tag]
        return {c: [mp[f"{f:g}"][c]["realised"][key] / mp["1"][c]["realised"][key] for f in fs] for c in CTX}

    ra0, rb0 = ratio_map("base", "A_own"), ratio_map("base", "B_own")
    f0, raw0 = level_opt(cl, AB0, X0, fs, ra0, rb0)
    out = dict(X_446bj2=X0, X_446b=Xb, AB_446bj2=AB0, base_level=f0, base_raw_levelled_ev=raw0, arms={})
    print(f"446bj2: X {X0:.6f} (446b {Xb:.6f}); energies A {[round(AB0[c]['A'], 5) for c in CTX]}; re-level {f0}; "
          f"EV raw at optimum {raw0:.5f} (board {RAW_446BJ2:.5f})")

    # ---------------- the anchor's own step, retrodicted by this instrument (comp420 rows, measured energies; no fit)
    kr = {}
    for f in FOLDS:
        R = json.loads((ROOT / f"experiments/420-bulk-stack/work/comp420-{f}.json").read_text())["arms"]
        b, a = R["stack|K3+sh1.75"]["realised"]["16"], R[BASE]["realised"]["16"]
        kr[f] = dict(r=a["cos_agg"] / b["cos_agg"] - 1, dpds_bd=a["pds_lam0.5"] - b["pds_lam0.5"], dpds_lit=a["pds"] - b["pds"])
    rk = float(np.mean([kr[f]["r"] for f in CANCER]))
    amp_k = np.sqrt(np.mean([AB0[c]["A"] for c in CTX]) / np.mean([ABb[c]["A"] for c in CTX]))
    pred = {f"k{k:.3f}": cl.raw(AB0, Xb * (1 + k * rk) * amp_k) - RAW_446B for k in KS}
    kgrid = np.linspace(0, 6, 6001)
    kfit = float(kgrid[int(np.argmin([abs(cl.raw(AB0, Xb * (1 + k * rk) * amp_k) - RAW_446B - BOARD_KEY["dmse"]) for k in kgrid]))])
    dpb = float(np.mean([kr[f]["dpds_bd"] for f in CANCER]))
    dpl = float(np.mean([kr[f]["dpds_lit"] for f in CANCER]))
    out["key_step_retrodiction"] = dict(plolo=kr, r_cancer=rk, dpds_bd=dpb, dpds_lit=dpl, energy_amp=float(amp_k),
                                        pred_dmse=pred, k_fit_to_board=kfit, board=BOARD_KEY,
                                        pred_dpds_bd_cal=K_PDS * dpb, pred_dpds_lit_cal=K_PDS_LIT * dpl)
    print(f"KEY STEP 446b -> 446bj2 (retrodiction): PLOLO r {rk:+.4f}, dpds bd {dpb:+.5f} / lit {dpl:+.5f}, energy amp {amp_k:.4f}")
    print(f"  predicted dmse raw " + " ".join(f"{k}: {v:+.5f}" for k, v in pred.items()) + f"  | board {BOARD_KEY['dmse']:+.5f} "
          f"-> k that reproduces it {kfit:.2f}")
    print(f"  predicted dpds raw: k_pds x bd {K_PDS * dpb:+.5f}, k_lit x lit {K_PDS_LIT * dpl:+.5f} | board {BOARD_KEY['dpds']:+.5f}")

    # ---------------- the arms
    R, Z = {}, {}
    for f in FOLDS:
        p = WORK / f"arms492g-{f}.json"
        if p.exists():
            R[f] = json.loads(p.read_text())
            Z[f] = np.load(WORK / f"arrays492g-{f}.npz")
    print(f"folds present: {list(R)}")
    for arm, (eab, mtag) in ARMS.items():
        if not R or arm not in next(iter(R.values()))["arms"]:
            continue
        e_arm, e_b = abc["arms"][eab]["E"], abc["arms"][eab]["E_base"]
        er = {c: dict(A=e_arm[c]["realised"]["A_own"] / e_b[c]["realised"]["A_own"],
                      B=e_arm[c]["realised"]["B_own"] / e_b[c]["realised"]["B_own"]) for c in CTX}
        ABa = {c: dict(A=AB0[c]["A"] * er[c]["A"], B=AB0[c]["B"] * er[c]["B"]) for c in CTX}
        amp = np.sqrt(np.mean([ABa[c]["A"] for c in CTX]) / np.mean([AB0[c]["A"] for c in CTX]))
        ram, rbm = ratio_map(mtag, "A_own"), ratio_map(mtag, "B_own")
        res = dict(energy_ratio={c: er[c]["A"] for c in CTX}, folds={})
        for f, Rf in R.items():
            b, a = Rf["arms"][BASE]["realised"]["16"], Rf["arms"][arm]["realised"]["16"]
            r = a["cos_agg"] / b["cos_agg"] - 1
            row = dict(r=r, dpds_bd=a["pds_lam0.5"] - b["pds_lam0.5"], dpds_lit=a["pds"] - b["pds"])
            for k in KS:
                fa, rawa = level_opt(cl, ABa, X0 * (1 + k * r) * amp, fs, ram, rbm)
                dm = rawa - raw0
                row[f"k{k:.3f}"] = dict(dmse=dm, overall=(A8.PDS_SC * row["dpds_bd"] + A8.MSE_SC * dm) / 6, level=fa)
            z = Z[f]
            d = z["d"]
            xb, ab_, pb = z[f"{BASE}@R16__x"], z[f"{BASE}@R16__a"], z[f"{BASE}@R16__pdsl"]
            xa, aa, pa = z[f"{arm}@R16__x"], z[f"{arm}@R16__a"], z[f"{arm}@R16__pdsl"]
            rng = np.random.default_rng(4920 + FOLDS.index(f))
            nn = d.size
            bo_r, bo_p = np.empty(NBOOT), np.empty(NBOOT)
            for i in range(NBOOT):
                ix = rng.integers(0, nn, nn)
                dd = d[ix].sum()
                cb = np.mean([xb[j, ix].sum() / np.sqrt(ab_[j, ix].sum() * dd) for j in range(3)])
                ca = np.mean([xa[j, ix].sum() / np.sqrt(aa[j, ix].sum() * dd) for j in range(3)])
                bo_r[i] = ca / cb - 1
                bo_p[i] = float(np.mean(pa[:, ix] - pb[:, ix]))
            k1 = row["k1.000"]
            dmse_dr = -2 * X0 * amp * np.mean([np.sqrt(np.interp(k1["level"][c], fs, ram[c])) for c in CTX]) / cl.sumD
            ov = (A8.PDS_SC * bo_p + A8.MSE_SC * (k1["dmse"] + dmse_dr * (bo_r - r))) / 6
            row["boot_k1"] = dict(r_ci=[float(np.quantile(bo_r, q)) for q in (0.025, 0.5, 0.975)],
                                  dpds_ci=[float(np.quantile(bo_p, q)) for q in (0.025, 0.5, 0.975)],
                                  overall_ci=[float(np.quantile(ov, q)) for q in (0.025, 0.5, 0.975)],
                                  r_draws=bo_r, ov_draws=ov)
            res["folds"][f] = row
        canc = [f for f in CANCER if f in res["folds"]]
        if canc:
            for k in KS:
                kk = f"k{k:.3f}"
                res[f"cancer_mean_{kk}"] = float(np.mean([res["folds"][f][kk]["overall"] for f in canc]))
                res[f"cancer_dmse_{kk}"] = float(np.mean([res["folds"][f][kk]["dmse"] for f in canc]))
                res[f"cancer_pos_{kk}"] = int(sum(res["folds"][f][kk]["overall"] > 0 for f in canc))
            res["cancer_r"] = float(np.mean([res["folds"][f]["r"] for f in canc]))
            res["cancer_dpds_bd"] = float(np.mean([res["folds"][f]["dpds_bd"] for f in canc]))
            res["cancer_dpds_lit"] = float(np.mean([res["folds"][f]["dpds_lit"] for f in canc]))
            pooled_r = np.mean([res["folds"][f]["boot_k1"]["r_draws"] for f in canc], axis=0)
            pooled_ov = np.mean([res["folds"][f]["boot_k1"]["ov_draws"] for f in canc], axis=0)
            res["cancer_pooled_boot_k1"] = dict(r_ci=[float(np.quantile(pooled_r, q)) for q in (0.025, 0.5, 0.975)],
                                                overall_ci=[float(np.quantile(pooled_ov, q)) for q in (0.025, 0.5, 0.975)])
            # ---- board-calibrated prediction vs 446bj2 (raw member deltas and overall)
            bl = bal["arms"].get(eab if eab != "P2 (no guard)" else "P2g", {})
            nm_D, nm_K = bl.get("nmae_raw_D", float("nan")), bl.get("nmae_raw_K", float("nan"))
            cal = {}
            for tag, kx in (("k_X 0.61", K_X), ("k_X 0.12", K_X_PESS)):
                dp = K_PDS * res["cancer_dpds_bd"]
                dm = res[f"cancer_dmse_k{kx:.3f}"]
                gross = (A8.PDS_SC * dp + A8.MSE_SC * dm) / 6
                dpl = K_PDS_LIT * res["cancer_dpds_lit"]
                gl = (A8.PDS_SC * dpl + A8.MSE_SC * dm) / 6          # PRIMARY: literal pds x 0.34 (the key step's check)
                cal[tag] = dict(dpds_raw=dp, dmse_raw=dm, gross=gross, net_D=gross + NMAE_OV * nm_D, net_K=gross + NMAE_OV * nm_K,
                                dpds_raw_lit=dpl, gross_lit=gl, net_D_lit=gl + NMAE_OV * nm_D, net_K_lit=gl + NMAE_OV * nm_K)
            res["calibrated"] = dict(ballast_nmae_raw_D=nm_D, ballast_nmae_raw_K=nm_K, **cal)
        for f in res["folds"]:
            for key in ("r_draws", "ov_draws"):
                res["folds"][f]["boot_k1"].pop(key)
        out["arms"][arm] = res
        line = f"{arm:18s} E {[round(er[c]['A'], 4) for c in CTX]}"
        for f in res["folds"]:
            rr = res["folds"][f]
            line += f" | {f[:8]} r {rr['r']:+.4f} dpds {rr['dpds_bd']:+.5f}/{rr['dpds_lit']:+.5f} ov1 {rr['k1.000']['overall']:+.5f}"
        if canc:
            c61 = res["calibrated"]["k_X 0.61"]
            c12 = res["calibrated"]["k_X 0.12"]
            line += (f" || cancer k1 {res['cancer_mean_k1.000']:+.5f} ({res['cancer_pos_k1.000']}/{len(canc)} +); "
                     f"CAL (lit pds x0.34): pds {c61['dpds_raw_lit']:+.5f} raw, mse {c61['dmse_raw']:+.5f} (k_X .61) / "
                     f"{c12['dmse_raw']:+.5f} (.12), net(D) {c61['net_D_lit']:+.5f} / {c12['net_D_lit']:+.5f}, net(K) "
                     f"{c61['net_K_lit']:+.5f} / {c12['net_K_lit']:+.5f}; [bd pds: {c61['dpds_raw']:+.5f}, net(D) {c61['net_D']:+.5f}]")
        print(line, flush=True)
    (WORK / "analyse492g.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
