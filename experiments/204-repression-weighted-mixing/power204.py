"""204 — how big a gain does the minimum-variance model itself predict, and does the
measurement have the power to see it?

A negative is only worth something if the instrument could have seen the positive.  The
model under test says the consensus of `e_i = A_i s + n_i` should weight by `A_i / N_i`
rather than uniformly.  Everything needed to price that is already measured: the per-donor
signal and noise energies from the Gram decomposition (`atten2_204.json`).

  variance of the equal-weight consensus of s:     (1/k^2) * sum_i v_i
  variance of the optimal-weight consensus:        1 / sum_i (1/v_i)      (Cauchy-Schwarz)
  with v_i = N_i / A_i^2  (each donor rescaled to unit signal, which is what the ratio of
  weights does anyway)

and the endpoint the board reads is the consensus's own reliability
`c = A / sqrt(A^2 + N)`, so a fractional drop `f` in the consensus noise raises the macro
cosine by `sqrt(1 / (1 - f * (1 - c_pred^2)))`.  `c_pred` comes from the measured PLOLO
macro cosine divided by the held-out donor's own measured reliability.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")


def main() -> None:
    at = json.loads((WORK / "atten2_204.json").read_text())
    fa = json.loads((WORK / "factor2_204.json").read_text())
    six = json.loads((WORK / "six204.json").read_text())

    out = {}
    for cols in ("common_genes", "all_official"):
        for scope in ("shared", "panel"):
            r = at["routeB_triple_gene_controlled"][cols][scope]
            A = np.array([np.sqrt(r["signal_energy"][d]) for d in DONORS])
            N = np.array([r["noise_energy"][d] for d in DONORS])
            v = N / A ** 2
            k = len(DONORS)
            var_eq = v.sum() / k ** 2
            var_opt = 1.0 / (1.0 / v).sum()
            f = 1.0 - var_opt / var_eq
            out[f"{cols}/{scope}"] = dict(
                snr=r["snr"], w_minvar_rel=r["w_minvar_rel"],
                noise_variance_drop=float(f))
            print(f"[{cols}/{scope}] optimal weighting removes {100*f:.1f}% of the "
                  f"consensus noise variance; w_minvar_rel="
                  f"{ {d: round(r['w_minvar_rel'][d], 3) for d in DONORS} }")

    # consensus reliability implied by the measured PLOLO macro cosine
    cos_plolo = six["cos_base_plolo"]
    c_truth = float(np.mean([fa["panel_block"]["c_rel"][d] for d in DONORS]))  # relative
    # absolute reliabilities from the panel block's own fit
    cs = json.loads((WORK / "factor204.json").read_text())["blocks"]["panel_block"]["c"]
    c_donor = float(np.mean([cs[d] for d in DONORS]))
    c_pred = cos_plolo / c_donor
    print(f"\nPLOLO macro cosine (equal weights, 4-fold mean) {cos_plolo:.6f}; "
          f"mean single-donor reliability {c_donor:.4f} -> implied consensus "
          f"reliability c_pred ~ {c_pred:.4f}")

    for key, rec in out.items():
        f = rec["noise_variance_drop"]
        gain = np.sqrt(1.0 / (1.0 - f * (1.0 - c_pred ** 2))) - 1.0
        rec["predicted_rel_cosine_gain"] = float(gain)
        rec["predicted_abs_cosine_gain_plolo"] = float(gain * cos_plolo)
        rec["predicted_overall_from_mse"] = float(0.3653 * 0.079249 * gain)
        print(f"  {key:26s} predicted cosine gain {100*gain:+.1f}% relative "
              f"(= {gain*cos_plolo:+.6f} PLOLO absolute, "
              f"{0.3653*0.079249*gain:+.5f} overall through `mse` alone)")

    meas = {a: six["arms"][a]["rel_cos"] for a in
            ("minvar", "minvar_self", "a_self", "deatten", "deatten_slope", "k13")}
    print("\nMEASURED relative cosine move:")
    for a, v in meas.items():
        print(f"  {a:16s} {100*v:+.2f}%")
    out["measured_rel_cos"] = meas
    out["cos_base_plolo"] = cos_plolo
    out["c_donor_panel_block"] = c_donor
    out["c_pred_implied"] = c_pred
    (WORK / "power204.json").write_text(json.dumps(out, indent=2))
    print("wrote work/power204.json")


if __name__ == "__main__":
    main()
