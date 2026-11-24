"""234 — the six-member price of every realised flip arm.  Routes are 228/232's, unchanged.

  `pds`   d pds / 0.452 / 6 * 1.011                          (PLOLO 3 cancer folds, a23)
  `mse`   197's inversion, `k_cos` carried as the RANGE {1.000, 1.384} -- not identified
  `nmae`  188's identity, / 0.5963 / 6
  `fid`   call-set sign precision / 0.30 / 6 (yield pinned at 1 by the shipped emitter)
  `jac`   relative, rescaled to the board's own raw `jac` -- EXACTLY zero for every arm
          here, which is 232 §11(a) re-confirmed on a different arm family

It also reports what the SIGN-ACCURACY EXCHANGE RATE would have predicted from `d q_w`
alone, so the gap between the two is visible: 232's rate was fitted on an oracle confined
to the reference's DE pool, and a real key is pool-blind.

    uv run python price234.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common234 import CANCER, GATE, QW_NEEDED, RATE_PER_QW, WORK  # noqa: E402
from phase0_234 import six_member  # noqa: E402


def main() -> None:
    rows = json.loads((WORK / "arm234.json").read_text())
    R = {(r["fold"], r["arm"]): r for r in rows}
    arms = sorted({r["arm"] for r in rows if r["arm"] != "BASE"})
    base = [R[(f, "BASE")] for f in CANCER]
    out = {}
    for arm in arms:
        if not all((f, arm) in R for f in CANCER):
            continue
        a = [R[(f, arm)] for f in CANCER]
        p = six_member(base, a, 1.384)
        p1 = six_member(base, a, 1.000)
        p["composite_k1000"] = p1["composite"]
        p["o_mse_k1000"] = p1["o_mse"]
        p["predicted_from_rate"] = {k: p["d_qw"] * v for k, v in RATE_PER_QW.items()}
        p["n_flipped"] = int(sum(R[(f, arm)]["n_flipped"] for f in CANCER))
        p["max_abs_d_jac"] = max(R[(f, arm)].get("max_abs_d_jac", 0.0) for f in CANCER)
        p["pds_positive"] = sum(1 for x, b in zip(a, base) if x["pds"] > b["pds"])
        if ("kolf2", arm) in R:
            k, kb = R[("kolf2", arm)], R[("kolf2", "BASE")]
            p["kolf2"] = dict(
                d_qw=k["q_w_pool"]["wald89"] - kb["q_w_pool"]["wald89"],
                d_pds=k["pds"] - kb["pds"],
                d_nmae=k["nmae_pools"]["wald89"] - kb["nmae_pools"]["wald89"],
                rel_cos_pct=100.0 * (k["cos_agg"] - kb["cos_agg"]) / abs(kb["cos_agg"]))
        out[arm] = p

    hdr = (f"{'arm':<22}{'flips':>9}{'d_qw':>10}{'pred(rate)':>12}{'ACTUAL':>10}"
           f"{'d_nmae':>10}{'rel cos%':>10}{'d_pds':>10}{'o_pds':>10}{'|djac|':>9}{'k+':>5}")
    print(f"\n{'='*len(hdr)}\nSIX-MEMBER PRICE, 3 cancer folds, shipped h1union stack, "
          f"wald89.  Gate +{GATE}\n{'='*len(hdr)}")
    print(hdr)
    print("-" * len(hdr))
    for arm in sorted(out, key=lambda a: -out[a]["composite"]):
        p = out[arm]
        print(f"{arm:<22}{p['n_flipped']:>9d}{p['d_qw']:>+10.5f}"
              f"{p['predicted_from_rate']['k1.384']:>+12.5f}{p['composite']:>+10.5f}"
              f"{p['d_nmae']:>+10.5f}{100*p['rel_cos']:>+10.3f}"
              f"{p['o_pds']*0.452*6/1.011:>+10.5f}{p['o_pds']:>+10.5f}"
              f"{p['max_abs_d_jac']:>9.1e}{p['pds_positive']:>4d}/3")

    print(f"\n{'arm':<22}{'o_nmae':>10}{'o_pds':>10}{'o_mse k1':>11}{'o_mse k1.4':>12}"
          f"{'o_fid':>10}{'o_jac':>9}{'TOT k1':>10}{'TOT k1.4':>10}{'kolf2 dqw':>11}")
    print("-" * 115)
    for arm in sorted(out, key=lambda a: -out[a]["composite"]):
        p = out[arm]
        kk = f"{p['kolf2']['d_qw']:+.5f}" if "kolf2" in p else "-"
        print(f"{arm:<22}{p['o_nmae']:>+10.5f}{p['o_pds']:>+10.5f}"
              f"{p['o_mse_k1000']:>+11.5f}{p['o_mse']:>+12.5f}{p['o_fid']:>+10.5f}"
              f"{p['o_jac']:>+9.5f}{p['composite_k1000']:>+10.5f}{p['composite']:>+10.5f}"
              f"{kk:>11}")

    best = max(out, key=lambda a: out[a]["composite"])
    print(f"\n  BEST composite anywhere: {best}  {out[best]['composite']:+.5f} "
          f"(k1.384) / {out[best]['composite_k1000']:+.5f} (k1.000)  against the "
          f"+{GATE} gate")
    print(f"  Gate needs d q_w >= +{QW_NEEDED['k1.384']:.5f}; the exchange rate is valid "
          "only for POOL-CONFINED sign changes.")
    (WORK / "price234.json").write_text(json.dumps(dict(arms=out, gate=GATE), indent=1))
    print(f"\n-> {WORK/'price234.json'}")


if __name__ == "__main__":
    main()
