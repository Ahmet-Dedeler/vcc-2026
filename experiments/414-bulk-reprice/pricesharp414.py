"""414: the recommended sharpening arm vs 410c with the render-space energy ratio read on the REAL 410c bulk (sharp410c.py)
instead of PLOLO's; direction (cos ratio) and pds from PLOLO (cancer folds and kolf2)."""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from analyse414 import CTX, MSE_SC, PDS_SC, joint187, fold_rows, mse_c
W = HERE / "work"
an = json.loads((W / "analyse414.json").read_text()); dec = json.loads((W / "decomp414.json").read_text())
sh = json.loads((W / "sharp410c.json").read_text())
m = joint187.Model(joint187.F_FIT); claim = {c: m.t[c]["claim"] for c in CTX}; Kc = {c: m.t[c]["Kc"] for c in CTX}; sumD = joint187.sumD0
T = an["plolo_cancer"]; k = an["k_key"]
AB_c = {c: dict(A=dec["410b"][c]["A"] * T["Kinf"]["Ar"][c] / T["K3"]["Ar"][c], B=dec["410b"][c]["B"] * T["Kinf"]["Br"][c] / T["K3"]["Br"][c]) for c in CTX}
X_c = an["amplitude"]["X410c"]; rel_c = T["Kinf"]["rel"]
raw_c = np.mean([mse_c(AB_c[c]["A"], AB_c[c]["B"], Kc[c], claim[c], X_c, sumD) for c in CTX])
out = {}
for q in ("1.15", "1.25", "1.5"):
    arm = f"Kinf+sh@{q}"
    ar = {c: sh[c][f"q{float(q):g}"]["energy_vs_q1"] for c in CTX}
    for label, folds in (("cancer", ("k562_gwps", "orion_hct116", "orion_hek293t")), ("kolf2", ("kolf2",))):
        rows = [fold_rows(f) for f in folds]
        cr = float(np.mean([r[arm]["cos_agg"] / r["Kinf"]["cos_agg"] for r in rows]))
        dpb = float(np.mean([r[arm]["pds_l05"] - r["Kinf"]["pds_l05"] for r in rows]))
        X = X_c * (1 + k * ((1 + rel_c) * cr - 1)) / (1 + k * rel_c) * np.sqrt(np.mean(list(ar.values())))
        raw = np.mean([mse_c(AB_c[c]["A"] * ar[c], AB_c[c]["B"] * ar[c], Kc[c], claim[c], X, sumD) for c in CTX])
        ov = (PDS_SC * dpb + MSE_SC * (raw - raw_c)) / 6
        out[f"{arm}|{label}"] = dict(energy_ratio=ar, cos_ratio=cr, dpds_board=dpb, d_raw=float(raw - raw_c), overall=float(ov))
        print(f"{arm:14s} {label:6s} real-bulk energy {[round(ar[c], 4) for c in CTX]} cos x{cr:.4f} dpds {dpb:+.4f} "
              f"d raw mse {raw - raw_c:+.4f} overall {ov:+.4f}")
(W / "pricesharp414.json").write_text(json.dumps(out, indent=1))
