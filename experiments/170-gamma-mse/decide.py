"""170 — apply AGENTS.md 11.1 rule 13 (the composed slot gate) to the measured deltas.

The gate, as written in `docs/CURRENT_STATE.md` (2026-09-02):

    predicted = PLOLO d_pds/6 * 1.011
              + 096 d_mse/6                (if amplitude / spread / phi / dispersion moves)
              - 0.0146/6                   (if the EMITTED PROFILE or emission law changes)
    SHIP only if predicted > +0.0019 (shape-only) or > +0.0089 (amplitude-changing)

Two adjustments, both stated rather than assumed:

* Standing consequence 1 of the same section: "any arm touching the emitted profile is charged
  -0.0146 on `reach` unless it is priced on `h1-test` DIRECTLY". This experiment prices reach
  on h1-test directly, so the measured h1-test `reach` delta replaces the prior -- but the
  prior's own error bar does not vanish, and h1-test reproduces board `reach` only to within
  0.014 raw, i.e. +-0.0026 overall. That is carried as the residual risk.
* A gamma arm is an OPERATOR arm. 162 classes the four operator/amplitude arms (a11, a12,
  a14, a19) together at an unpriced RMS of 0.00890, against 0.00188 for the two shape-only
  arms. Both thresholds are reported; the operator threshold is the one that governs.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
from report import ANCHOR, mse_scaled, scaled  # noqa: E402

BASE = 0.65
SHAPE_ONLY_FLOOR, AMPLITUDE_FLOOR = 0.0019, 0.0089
H1_REACH_TRANSFER_SIGMA = 0.014 / (ANCHOR["reach"][1] - ANCHOR["reach"][0]) / 6


def main() -> None:
    rows = {r["gamma"]: r for r in json.loads((WORK / "h1sweep.json").read_text())}
    base = rows[BASE]
    twin = None
    if (WORK / "h1sweep_s1.json").exists():
        twin = {r["gamma"]: r for r in json.loads((WORK / "h1sweep_s1.json").read_text())}

    chal = None
    cs = WORK / "chal_minnorm.json"
    cv = WORK / "chal_minnorm_solve.json"
    if cs.exists() and cv.exists():
        crows = {r["gamma"]: r for r in json.loads(cs.read_text())}
        sol = json.loads(cv.read_text())
        chal = {}
        for gam, r in crows.items():
            Xc = sol["X_at_065"] * (rows[gam]["X"] / base["X"]) if gam in rows else None
            if Xc is None:
                continue
            raw = float(np.mean([
                (r["per"][c]["A"] + r["per"][c]["Kc"]
                 - min(r["per"][c]["claim"], r["per"][c]["B"]) - 2 * Xc) / sol["sumD"] + 1.0
                for c in ("A", "B", "C")]))
            chal[gam] = raw
        base_score = mse_scaled(chal[BASE])

    print(f"{'gamma':>6}{'d pds':>10}{'d mse(h1)':>11}{'d mse(board)':>13}{'d nmae':>10}"
          f"{'d reach':>10}{'d fid':>10}{'d jac':>10}{'PREDICTED':>11}{'verdict':>26}")
    for gam in sorted(rows):
        r = rows[gam]
        d_pds = (scaled("pds", r["pds"]) - scaled("pds", base["pds"])) / 6
        d_mse_h1 = (mse_scaled(r["mse"]) - mse_scaled(base["mse"])) / 6
        d_mse_bd = ((mse_scaled(chal[gam]) - base_score) / 6
                    if chal and gam in chal else float("nan"))
        d_nmae = (scaled("nmae", r["nmae"]) - scaled("nmae", base["nmae"])) / 6
        d = {m: (scaled(m, r[m]) - scaled(m, base[m])) / 6 for m in ("reach", "fid", "jac")}
        pred = (d_pds * 1.011 + (d_mse_bd if d_mse_bd == d_mse_bd else 0.0)
                + d_nmae + d["reach"] + d["fid"] + d["jac"])
        verdict = ("SHIP (operator floor)" if pred > AMPLITUDE_FLOOR else
                   "below operator floor" if pred > SHAPE_ONLY_FLOOR else
                   "VETO")
        print(f"{gam:>6.2f}{d_pds:>+10.5f}{d_mse_h1:>+11.5f}{d_mse_bd:>+13.5f}"
              f"{d_nmae:>+10.5f}{d['reach']:>+10.5f}{d['fid']:>+10.5f}{d['jac']:>+10.5f}"
              f"{pred:>+11.5f}{verdict:>26}")

    if twin:
        sig = []
        for gam, tw in twin.items():
            b = rows[gam]
            tot = ((mse_scaled(tw["mse"]) - mse_scaled(b["mse"]))
                   + (scaled("pds", tw["pds"]) - scaled("pds", b["pds"]))
                   + (scaled("nmae", tw["nmae"]) - scaled("nmae", b["nmae"]))
                   + sum(scaled(m, tw[m]) - scaled(m, b[m])
                         for m in ("reach", "fid", "jac"))) / 6
            sig.append(abs(tot))
        print(f"\nre-seeded twin |delta overall| at fixed gamma: "
              f"{['%+.5f' % s for s in sig]}  -- the emission noise floor on this surface")
    print(f"h1-test -> board `reach` transfer sigma: +-{H1_REACH_TRANSFER_SIGMA:.5f} overall")
    print(f"rule 13 floors: shape-only +{SHAPE_ONLY_FLOOR}, operator/amplitude "
          f"+{AMPLITUDE_FLOOR}")


if __name__ == "__main__":
    main()
