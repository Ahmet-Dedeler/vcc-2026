"""492 calibration: the instrument's transfer to the board, fitted on the board's own levelled sharpening steps
(433h -> 433k: q 1.25 -> 1.5; 433k -> 446b: q 1.5 -> 1.75; same stack, key and DE channel). PLOLO rows are 420's
(`comp420-<fold>.json`, stack|K3+sh{q}, realised KMAX 16), A/B/C energies 420's (`abc420.json`, realised 16), the cliff
model 187's anchored at 446b (as analyse473/492). For each step: k_pds = board dpds / PLOLO dpds (cancer mean, literal and
board-emulated), and k_X such that the levelled-cliff price of the step reproduces the board's dmse. Then 492/492d/493's
arms are priced with those k.

    uv run python experiments/492-per-target-sharpening/calib492.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/418-bulk-composition"))
sys.path.insert(0, str(ROOT / "experiments/473-410d2-reprice"))
import analyse418 as A8  # noqa: E402
from analyse473 import level_opt  # noqa: E402

CTX = ("A", "B", "C")
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
RAW_446B = 0.6408807601228349
BOARD = {("1.25", "1.5"): dict(dpds=0.00462281, dmse=-0.01789753), ("1.5", "1.75"): dict(dpds=0.00316983, dmse=-0.00764239)}


def main():
    cl = A8.Cliff()
    ntc = json.loads((ROOT / "experiments/446-q175-level/work/ntclevel-446b.json").read_text())["per_context"]
    AB446 = {c: dict(A=ntc[c]["A"], B=ntc[c]["B"]) for c in CTX}
    X446 = cl.solve_X(AB446, RAW_446B)
    a420 = json.loads((ROOT / "experiments/420-bulk-stack/work/abc420.json").read_text())["designs"]
    E = {q: {c: a420[f"stack|K3+sh{q}"]["realised"]["16"][c] for c in CTX} for q in ("1.25", "1.5", "1.75")}
    R = {f: json.loads((ROOT / f"experiments/420-bulk-stack/work/comp420-{f}.json").read_text())["arms"] for f in CANCER}
    cos = {q: np.mean([R[f][f"stack|K3+sh{q}"]["realised"]["16"]["cos_agg"] for f in CANCER]) for q in E}
    abc = json.loads((HERE / "work/abc492.json").read_text())
    mp = abc["level_maps"]["maps"]["base"]
    fs = [float(x) for x in abc["level_maps"]["factors"]]
    ra = {c: [mp[f"{f:g}"][c]["realised"]["A_own"] / mp["1"][c]["realised"]["A_own"] for f in fs] for c in CTX}
    rb = {c: [mp[f"{f:g}"][c]["realised"]["B_own"] / mp["1"][c]["realised"]["B_own"] for f in fs] for c in CTX}

    def raw_lev(q, k):
        """levelled EV raw of design q, anchored at 446b (q 1.75): energies by 420's realised ratios, X by the PLOLO cos."""
        AB = {c: dict(A=AB446[c]["A"] * E[q][c]["A"] / E["1.75"][c]["A"], B=AB446[c]["B"] * E[q][c]["B"] / E["1.75"][c]["B"]) for c in CTX}
        r = cos[q] / cos["1.75"] - 1
        X = X446 * (1 + k * r) * np.sqrt(np.mean([AB[c]["A"] for c in CTX]) / np.mean([AB446[c]["A"] for c in CTX]))
        return level_opt(cl, AB, X, fs, ra, rb)[1]

    out = {"steps": {}}
    for (qb, qa), bd in BOARD.items():
        dl = np.mean([R[f][f"stack|K3+sh{qa}"]["realised"]["16"]["pds"] - R[f][f"stack|K3+sh{qb}"]["realised"]["16"]["pds"] for f in CANCER])
        db = np.mean([R[f][f"stack|K3+sh{qa}"]["realised"]["16"]["pds_lam0.5"] - R[f][f"stack|K3+sh{qb}"]["realised"]["16"]["pds_lam0.5"] for f in CANCER])
        grid = np.linspace(0, 2.0, 201)
        dm = np.array([raw_lev(qa, k) - raw_lev(qb, k) for k in grid])
        kx = float(grid[int(np.argmin(np.abs(dm - bd["dmse"])))])
        out["steps"][f"{qb}->{qa}"] = dict(k_pds_lit=bd["dpds"] / dl, k_pds_bd=bd["dpds"] / db, k_X=kx,
                                           dmse_k1=float(dm[100 // 2 * 1 if False else int(np.argmin(np.abs(grid - 1.0)))]),
                                           board=bd, plolo_dpds_lit=float(dl), plolo_dpds_bd=float(db))
        print(f"step q {qb}->{qa}: board dpds {bd['dpds']:+.5f} dmse {bd['dmse']:+.5f} | PLOLO dpds lit {dl:+.5f} bd {db:+.5f} "
              f"-> k_pds lit {bd['dpds'] / dl:.3f} bd {bd['dpds'] / db:.3f}; levelled dmse at k 1 {out['steps'][f'{qb}->{qa}']['dmse_k1']:+.5f} "
              f"-> k_X fitted {kx:.2f}", flush=True)
    kp = float(np.mean([s["k_pds_bd"] for s in out["steps"].values()]))
    # k_X: the step adjacent to 446b (q 1.5 -> 1.75) is the calibration that applies at 446b's base; the first step's
    # energy route is anchored far from its own base and fits poorly (reported as the pessimistic bound)
    kx = out["steps"]["1.5->1.75"]["k_X"]
    kx_lo = out["steps"]["1.25->1.5"]["k_X"]
    out["k_pds_bd"], out["k_X"], out["k_X_pessimistic"] = kp, kx, kx_lo
    # apply to 492 / 492d arms (their analyse JSONs hold per-fold r, dpds_bd; re-price dmse at k_X through the same pricer)
    sys.argv = [sys.argv[0]]
    for tag, fn in (("492", "analyse492.json"), ("492b", "analyse492b.json"), ("492c", "analyse492c.json"), ("492d", "analyse492d.json")):
        d = json.loads((HERE / "work" / fn).read_text())
        for arm, x in d["arms"].items():
            f3 = [f for f in CANCER if f in x["folds"]]
            dp = float(np.mean([x["folds"][f]["dpds_bd"] for f in f3]))
            # dmse scales ~linearly in k around k 1 (X enters linearly at a fixed level): interpolate between k 1 and 1.407
            d1 = float(np.mean([x["folds"][f]["k1.000"]["dmse"] for f in f3]))
            d2 = float(np.mean([x["folds"][f]["k1.407"]["dmse"] for f in f3]))
            dmk = d1 + (kx - 1.0) * (d2 - d1) / 0.407
            dml = d1 + (kx_lo - 1.0) * (d2 - d1) / 0.407
            ov = (A8.PDS_SC * kp * dp + A8.MSE_SC * dmk) / 6
            ovl = (A8.PDS_SC * kp * dp + A8.MSE_SC * dml) / 6
            out.setdefault("calibrated", {})[f"{tag}:{arm}"] = dict(overall=ov, overall_pessimistic=ovl, dpds_board=kp * dp,
                                                                   dmse_board=dmk, dmse_board_pessimistic=dml)
            print(f"{tag:5s} {arm:7s} calibrated board (k_pds {kp:.3f}, k_X {kx:.2f} / pess. {kx_lo:.2f}): pds {kp * dp:+.5f} raw, "
                  f"mse {dmk:+.5f} / {dml:+.5f} raw -> overall {ov:+.5f} / {ovl:+.5f}", flush=True)
    (HERE / "work/calib492.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
