"""246 — the six-member price of the surviving arms, every member explicit, never at zero.

Route (218 §7.4 / AGENTS.md §8.4):
  * the five cosine-driven members (`mse`, `reach`, `nmae`, `fid`, `jac`) from `cd4f034`'s
    realised board rate card per +0.01 ABSOLUTE cosine, on the MACRO PLOLO cosine delta
    (uncalibrated, k_cos = 1; the factor is not identified, 1.000 on `h1union`, 1.384 on
    `cd4f034`);
  * `pds` from PLOLO directly (raw Δ / 0.452 / 6 × 1.011) at THREE transfer factors: 0.7×
    (`h1union` realised 70 % of its predicted `pds` gain, NEXT-BUILD-SPEC §4), 1.0×, and 2.0×
    (§8.4's realised COST factor on `cd4f034`, k_pds = 1.987) -- the last is what
    `price218.price` charges and is reported as the optimistic end for a GAIN;
  * the exact `nmae` identity (six246) as a MAGNITUDE BOUND alongside the rate-card `nmae`;
  * the sign channel: Δq_w × 232's +0.00399 per +0.01 (one-sided 234: +0.00257), reported
    for information, not added (it is already inside the rate card's `nmae`/`reach`/`fid`).

    uv run python price246.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))
from price218 import BOARD_COS, GATE_AMPL, GATE_SHAPE, RATE, SIGMA_PDS, SPAN_PDS  # noqa: E402

ARMS = ["B2_cd4_svd2", "B2_cd4_svd3", "B2_cd4_svd5", "B2_cd4_svd7", "B2_cd4_svd10",
        "B2_cd4_svd3_resolved", "B2_cd4_svd5_resolved", "B2_cd4_svd10_resolved",
        "B3_cd4_svd5_splithalf", "B3_cd4_svd5_splithalfb", "B3_cd4_svd10_splithalf",
        "C2_partial@0.34", "E_cd4_svd3+partial@0.34", "E_cd4_svd5+partial@0.34",
        "A_svd2@0.1", "A_svd2@0.2", "DOSE_cd4@0.28", "B2_cd4_randE5", "B2_cd4_randE5b"]
SIX = {"B2_cd4_svd2": "B2_cd4_svd2", "B2_cd4_svd3": "B2_cd4_svd3", "B2_cd4_svd5": "B2_cd4_svd5",
       "B2_cd4_svd7": "B2_cd4_svd7", "B2_cd4_svd10": "B2_cd4_svd10",
       "B2_cd4_svd3_resolved": "B2_cd4_svd3_resolved",
       "B2_cd4_svd5_resolved": "B2_cd4_svd5_resolved", "C2_partial@0.34": "C_partial@0.34",
       "E_cd4_svd3+partial@0.34": "E_cd4_svd3+partial@0.34",
       "E_cd4_svd5+partial@0.34": "E_cd4_svd5+partial@0.34"}


def main() -> None:
    a = json.loads((WORK / "analyse246.json").read_text())
    six = json.loads((WORK / "six246.json").read_text())["arms"]
    summ = {}
    for title, fam in a["report"].items():
        if "kolf2" in title:
            continue
        summ.update(fam)
    out = {}
    print(f"{'arm':26s}{'macro%':>8s}{'d pds':>10s}{'sig':>5s}{'mse':>9s}{'reach':>9s}"
          f"{'nmae(rc)':>9s}{'nmae(ex)':>9s}{'fid':>9s}{'jac':>9s}{'pds x0.7':>9s}"
          f"{'pds x1':>8s}{'pds x2':>8s}{'TOT x0.7':>9s}{'TOT x1':>8s}{'TOT x2':>8s}"
          f"{'dq_w':>9s}{'q_w ovr':>8s}{'e':>6s}")
    for arm in ARMS:
        s = summ.get(arm)
        if s is None:
            continue
        rm, dp, dq = s["rel_macro"], s["d_pds"], s["d_qw"]
        dcos = BOARD_COS * rm
        legs = {m: RATE[m] * dcos / 0.01 for m in RATE if m != "pds"}
        pds1 = dp / SPAN_PDS / 6.0 * 1.011
        tot = {f: sum(legs.values()) + f * pds1 for f in (0.7, 1.0, 2.0)}
        ex = six.get(SIX.get(arm, ""), {})
        nm_ex = ex.get("overall_nmae_vs_inc", float("nan"))
        e = ex.get("e_vs_inc", float("nan"))
        qw_ovr = dq / 0.01 * 0.00399
        out[arm] = dict(rel_macro=rm, d_pds=dp, pds_sigma=abs(dp) / SIGMA_PDS,
                        legs=legs, pds_x1=pds1, totals={str(f): v for f, v in tot.items()},
                        nmae_exact_overall=nm_ex, e_exact=e, d_qw=dq, qw_overall_232=qw_ovr,
                        delivered_ratio_vs_a15=ex.get("delivered_norm_ratio"),
                        d_precision89=ex.get("d_precision_vs_inc"))
        print(f"{arm:26s}{100*rm:+8.3f}{dp:+10.6f}{abs(dp)/SIGMA_PDS:5.1f}"
              f"{legs['mse']:+9.5f}{legs['reach']:+9.5f}{legs['nmae']:+9.5f}{nm_ex:+9.5f}"
              f"{legs['fid']:+9.5f}{legs['jac']:+9.5f}{0.7*pds1:+9.5f}{pds1:+8.5f}"
              f"{2*pds1:+8.5f}{tot[0.7]:+9.5f}{tot[1.0]:+8.5f}{tot[2.0]:+8.5f}"
              f"{dq:+9.5f}{qw_ovr:+8.5f}{e:6.3f}")
    print(f"\ngates: shape-only {GATE_SHAPE:+.4f}   amplitude-changing {GATE_AMPL:+.4f}   "
          f"PLOLO pds LOO sigma {SIGMA_PDS} raw")
    (WORK / "price246.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {WORK/'price246.json'}")


if __name__ == "__main__":
    main()
