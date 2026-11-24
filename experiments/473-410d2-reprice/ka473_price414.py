"""473 KA-4: 414's bounded-arm price, recomputed from 414's saved inputs with 414's own functions (price410d.main's
arithmetic, output redirected to 473/work). Must reproduce price410d.json's `matched_overall`:
cancer +0.0045474, kolf2 +0.0025158 (|d| < 1e-12).

    uv run python experiments/473-410d2-reprice/ka473_price414.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True          # importing 414's modules must not drop __pycache__ into their dirs
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W414 = ROOT / "experiments/414-bulk-reprice"
sys.path.insert(0, str(W414))
from analyse414 import CTX, MSE_SC, PDS_SC, joint187, fold_rows, mse_c  # noqa: E402

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")


def main():
    ref = json.loads((W414 / "work/price410d.json").read_text())
    an = json.loads((W414 / "work/analyse414.json").read_text())
    dec = json.loads((W414 / "work/decomp414.json").read_text())
    mdl = joint187.Model(joint187.F_FIT)
    claim = {c: mdl.t[c]["claim"] for c in CTX}
    Kc = {c: mdl.t[c]["Kc"] for c in CTX}
    sumD = joint187.sumD0
    Xb, k = an["Xb"], an["k_key"]
    T = an["plolo_cancer"]
    A410b = {c: dec["410b"][c] for c in CTX}
    AB_c = {c: dict(A=A410b[c]["A"] * T["Kinf"]["Ar"][c] / T["K3"]["Ar"][c],
                    B=A410b[c]["B"] * T["Kinf"]["Br"][c] / T["K3"]["Br"][c]) for c in CTX}
    X_c = an["amplitude"]["X410c"]
    out = {}
    worst = 0.0
    arm = "Kinf+gf@0.25"
    for label, folds in (("cancer folds", CANCER), ("kolf2 (= the shipped 410d multiplier for the gain arm)", ("kolf2",))):
        rows = [fold_rows(f) for f in folds]
        ar = {c: float(np.mean([r[arm]["AB"][c]["A"] / r["Kinf"]["AB"][c]["A"] for r in rows])) for c in CTX}
        br = {c: float(np.mean([(r[arm]["AB"][c]["B"] if np.isfinite(r[arm]["AB"][c]["B"]) else r[arm]["AB"][c]["A"])
                                / r["Kinf"]["AB"][c]["B"] for r in rows])) for c in CTX}
        cr = float(np.mean([r[arm]["cos_agg"] / r["Kinf"]["cos_agg"] for r in rows]))
        dpb = float(np.mean([r[arm]["pds_l05"] - r["Kinf"]["pds_l05"] for r in rows]))
        AB_d = {c: dict(A=AB_c[c]["A"] * ar[c], B=AB_c[c]["B"] * br[c]) for c in CTX}
        rel_c = T["Kinf"]["rel"]
        s = {c: 1 / np.sqrt(ar[c]) for c in CTX}
        X_m = X_c * (1 + k * ((1 + rel_c) * cr - 1)) / (1 + k * rel_c)
        raw_c = float(np.mean([mse_c(AB_c[c]["A"], AB_c[c]["B"], Kc[c], claim[c], X_c, sumD) for c in CTX]))
        raw_m = float(np.mean([mse_c(AB_d[c]["A"] * s[c] ** 2, AB_d[c]["B"] * s[c] ** 2, Kc[c], claim[c], X_m, sumD)
                               for c in CTX]))
        mo = (PDS_SC * dpb + MSE_SC * (raw_m - raw_c)) / 6
        want = ref[f"{arm} | {label}"]["matched_overall"]
        worst = max(worst, abs(mo - want))
        out[label] = dict(matched_overall=mo, ref=want, cos_ratio=cr, dpds_board=dpb, A_ratio=ar)
        print(f"KA-4 {label:55s} matched_overall {mo:+.7f} (414: {want:+.7f}, |d| {abs(mo - want):.1e})", flush=True)
    assert worst < 1e-12, f"KA-4 FAILED: {worst}"
    out["worst"] = worst
    (HERE / "work/ka473_price414.json").write_text(json.dumps(out, indent=1, default=float))
    print("KA-4 PASS")


if __name__ == "__main__":
    main()
