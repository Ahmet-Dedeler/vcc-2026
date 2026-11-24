"""232 Phase 3 — the six-member price, three cancer folds, `kolf2` apart.  Led by `nmae`.

ROUTES ARE 228's, UNCHANGED (rule 12).  Reusing them is the point: the result has to be
comparable with 228's line for line, which is why `M`, `U` and `e = dM/dU` are carried
alongside every row.

  `pds`   PLOLO 3 cancer folds, a23 render.  `d pds / 0.452 / 6 * 1.011`.
  `mse`   197's inversion `d raw mse = -2.074495 * d cos_board`, `/0.950162566 / 6`,
          `d cos_board = 0.099469 * rel * k_cos`, `k_cos` carried as the RANGE
          {1.000, 1.384} because §8.1c found a 7.84x spread and it is NOT identified.
  `nmae`  188's exact identity, `/0.5963 / 6`.
  `fid`   call-set sign precision, yield pinned at 1 by the shipped emitter.  `/0.30 / 6`.
  `jac`   `/0.365 / 6`, absolute and rescaled to the board's own raw `jac`.
  `reach` ordering half only, indicative (its yield half carries `n_real`, rule 9).

ALSO HERE: the SIGN-ACCURACY EXCHANGE RATE.  The `oracle`/`antioracle` rows are a
controlled two-sided dose response in |t|-weighted sign accuracy `q_w` at BIT-IDENTICAL
magnitude.  Regressing each member on `d q_w` across those rows gives, for the first time,
what a unit of sign accuracy is worth in overall points -- which is what decides whether
ANY sign-combination rule can pay, not just the ones measured here.

    uv run python price232.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
SPAN_PDS, CALIB, SIGMA_PDS = 0.452, 1.011, 0.001494
SPAN_MSE, DMSE_DCOS = 0.950162566, 2.074495
SPAN_NMAE = 0.5963
SPAN_FID, SPAN_JAC, SPAN_REACH = 0.30, 0.365, 0.895
OURS_COS = 0.09946902740184964
OURS_JAC_RAW, OURS_REACH_RAW = 0.04095656108311888, 0.27654544083789157
K_COS = (1.000, 1.384)
GATE_SHAPE = 0.0019
ORDERING_NOISE = 0.004502          # §8.1b


def price(rows_base, rows_arm) -> dict:
    def mean(fn):
        return float(np.mean([fn(a, b) for a, b in zip(rows_arm, rows_base)]))

    rel = mean(lambda a, b: 100.0 * (a["cos_agg"] - b["cos_agg"]) / abs(b["cos_agg"]))
    dpds = mean(lambda a, b: a["pds"] - b["pds"])
    dnm = {k: mean(lambda a, b, k=k: a["nmae_pools"][k] - b["nmae_pools"][k])
           for k in ("wald89", "lfc89", "h89")}
    dU = mean(lambda a, b: a["nmae_U"]["wald89"] - b["nmae_U"]["wald89"])
    dM = mean(lambda a, b: a["nmae_M"]["wald89"] - b["nmae_M"]["wald89"])
    dqw = mean(lambda a, b: a["q_w_pool"]["wald89"] - b["q_w_pool"]["wald89"])
    dpr = mean(lambda a, b: a["reach_prec89"] - b["reach_prec89"])
    dfi = mean(lambda a, b: a["fid_precision"] - b["fid_precision"])
    dja = mean(lambda a, b: a["jac"]["89"] - b["jac"]["89"])
    rja = mean(lambda a, b: a["jac"]["89"] / b["jac"]["89"] - 1.0)
    dnorm = mean(lambda a, b: a["delivered_median_norm"] / b["delivered_median_norm"] - 1.0)
    npos = sum(1 for a, b in zip(rows_arm, rows_base) if a["pds"] > b["pds"])
    nneg_nmae = sum(1 for a, b in zip(rows_arm, rows_base)
                    if a["nmae_pools"]["wald89"] < b["nmae_pools"]["wald89"])

    o_pds = dpds / SPAN_PDS / 6.0 * CALIB
    o_mse = {f"k{k:g}": (DMSE_DCOS * OURS_COS * (rel / 100.0) * k) / SPAN_MSE / 6.0
             for k in K_COS}
    o_nmae = {k: -v / SPAN_NMAE / 6.0 for k, v in dnm.items()}
    o_fid = dfi / SPAN_FID / 6.0
    o_jac_abs = dja / SPAN_JAC / 6.0
    o_jac_rel = rja * OURS_JAC_RAW / SPAN_JAC / 6.0
    o_reach = dpr * (OURS_REACH_RAW / 0.62) / SPAN_REACH / 6.0

    comp = {f"k{k:g}|{pool}|{jm}": (o_pds + o_mse[f"k{k:g}"] + o_nmae[pool] + o_fid
                                    + (o_jac_abs if jm == "abs" else o_jac_rel))
            for k in K_COS for pool in ("wald89", "lfc89", "h89") for jm in ("abs", "rel")}
    return dict(rel_cos_agg_pct=rel, d_pds=dpds, pds_sigma=dpds / SIGMA_PDS,
                d_nmae=dnm, dU=dU, dM=dM, e=(dM / dU if dU != 0 else float("nan")),
                d_qw=dqw, d_prec89=dpr, d_fid=dfi, d_jac89=dja, rel_jac89=rja,
                d_delivered_norm_pct=100.0 * dnorm,
                pds_positive_cancer=npos, nmae_improved_cancer=nneg_nmae,
                o_pds=o_pds, o_mse=o_mse, o_nmae=o_nmae, o_fid=o_fid,
                o_jac_abs=o_jac_abs, o_jac_rel=o_jac_rel, o_reach_indicative=o_reach,
                composite=comp,
                composite_with_reach={k: v + o_reach for k, v in comp.items()})


def main() -> None:
    rows = json.loads((WORK / "arm232.json").read_text())
    R = {(r["fold"], r["arm"]): r for r in rows}
    arms = sorted({r["arm"] for r in rows if r["arm"] != "BASE"})
    out, empty = {}, []

    for arm in arms:
        if not all((f, arm) in R for f in CANCER):
            continue
        if any(R[(f, arm)].get("empty") for f in CANCER):
            empty.append(arm)
            continue
        rb = [R[(f, "BASE")] for f in CANCER]
        ra = [R[(f, arm)] for f in CANCER]
        p = price(rb, ra)
        p["n_flipped_cancer"] = int(sum(R[(f, arm)]["n_flipped"] for f in CANCER))
        if ("kolf2", arm) in R and not R[("kolf2", arm)].get("empty"):
            k, kb = R[("kolf2", arm)], R[("kolf2", "BASE")]
            p["kolf2"] = dict(
                rel_cos_agg_pct=100.0 * (k["cos_agg"] - kb["cos_agg"]) / abs(kb["cos_agg"]),
                d_pds=k["pds"] - kb["pds"],
                d_nmae=k["nmae_pools"]["wald89"] - kb["nmae_pools"]["wald89"],
                d_qw=k["q_w_pool"]["wald89"] - kb["q_w_pool"]["wald89"])
        out[arm] = p

    hdr = (f"{'arm':<34}{'flips':>8}{'d_qw':>10}{'d_nmae':>10}{'e':>8}{'rel cos%':>10}"
           f"{'d_pds':>10}{'pds+':>6}")
    print("\nPHASE 3 — three cancer folds, shipped h1union stack, wald89")
    print(hdr)
    print("-" * len(hdr))
    for arm in sorted(out, key=lambda a: -out[a]["composite"]["k1.384|wald89|rel"]):
        p = out[arm]
        print(f"{arm:<34}{p['n_flipped_cancer']:>8d}{p['d_qw']:>+10.5f}"
              f"{p['d_nmae']['wald89']:>+10.5f}{p['e']:>8.3f}"
              f"{p['rel_cos_agg_pct']:>+10.3f}{p['d_pds']:>+10.5f}"
              f"{p['pds_positive_cancer']:>4d}/3")

    hdr2 = (f"{'arm':<34}{'o_nmae':>10}{'o_pds':>10}{'o_mse k1':>10}{'o_mse k1.4':>11}"
            f"{'o_fid':>9}{'o_jac':>9}{'TOT k1':>10}{'TOT k1.4':>10}")
    print("\nSIX-MEMBER PRICE (overall points).  Gate +0.0019.")
    print(hdr2)
    print("-" * len(hdr2))
    for arm in sorted(out, key=lambda a: -out[a]["composite"]["k1.384|wald89|rel"]):
        p = out[arm]
        print(f"{arm:<34}{p['o_nmae']['wald89']:>+10.5f}{p['o_pds']:>+10.5f}"
              f"{p['o_mse']['k1']:>+10.5f}{p['o_mse']['k1.384']:>+11.5f}"
              f"{p['o_fid']:>+9.5f}{p['o_jac_rel']:>+9.5f}"
              f"{p['composite']['k1|wald89|rel']:>+10.5f}"
              f"{p['composite']['k1.384|wald89|rel']:>+10.5f}")

    # ---------------- the sign-accuracy exchange rate, from the oracle dose response
    ex = {}
    orc = [a for a in out if "oracle@" in a or "antioracle@" in a]
    if orc:
        x = np.array([out[a]["d_qw"] for a in orc])
        for key, fn in (("o_nmae", lambda p: p["o_nmae"]["wald89"]),
                        ("o_pds", lambda p: p["o_pds"]),
                        ("o_mse_k1.384", lambda p: p["o_mse"]["k1.384"]),
                        ("o_fid", lambda p: p["o_fid"]),
                        ("o_jac", lambda p: p["o_jac_rel"]),
                        ("composite_k1.384", lambda p: p["composite"]["k1.384|wald89|rel"]),
                        ("d_nmae", lambda p: p["d_nmae"]["wald89"]),
                        ("rel_cos_agg_pct", lambda p: p["rel_cos_agg_pct"])):
            y = np.array([fn(out[a]) for a in orc])
            slope = float(np.linalg.lstsq(x[:, None], y, rcond=None)[0][0])
            r = float(np.corrcoef(x, y)[0, 1]) if len(x) > 2 else float("nan")
            ex[key] = dict(per_unit_qw=slope, per_0p01_qw=slope * 0.01, r=r)
        print("\nSIGN-ACCURACY EXCHANGE RATE (oracle/antioracle dose response, "
              f"{len(orc)} rows, magnitude bit-identical)")
        print(f"{'member':<22}{'per +0.01 q_w':>16}{'r':>8}")
        print("-" * 46)
        for k, v in ex.items():
            print(f"{k:<22}{v['per_0p01_qw']:>+16.6f}{v['r']:>8.3f}")

    if empty:
        print(f"\nEMPTY-MASK arms (identical to BASE by construction): {', '.join(empty)}")
    (WORK / "price232.json").write_text(json.dumps(
        dict(arms=out, exchange_rate=ex, empty_mask=empty,
             gate=GATE_SHAPE, ordering_noise=ORDERING_NOISE), indent=1))
    print(f"\n-> {WORK / 'price232.json'}")


if __name__ == "__main__":
    main()
