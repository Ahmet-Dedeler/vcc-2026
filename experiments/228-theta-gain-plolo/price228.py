"""228 Phase 5 — the six-member price, at fixed alpha and at amplitude-matched alpha.

ROUTES, and why each member takes the one it takes (rule 12):

  `pds`   PLOLO, 3 cancer folds, a23 render.  `d pds / 0.452 / 6 * 1.011` (162; LOO sigma
          0.001494 raw, sign 6/6).  `kolf2` reported apart.
  `mse`   197's inversion, which is the SAME object as 096's `claim` route and is the one
          validated to 4.9e-7 out of sample:  `d raw mse = -2.074495 * d cos_board`, then
          `/0.950162566 / 6`.  `d cos_board = ours_cos * rel_cos_agg * k_cos`, and `k_cos` is
          NOT identified -- 1.000 on `h1union`, 1.384 on `cd4f034`, both on this functional --
          so it is carried as a RANGE and never as a point.
  `nmae`  188's EXACT identity on the PLOLO folds (`|lfc|*sqrt(cpm)` pool at the reference's
          BH depth, 89 genes/row), all three pools; cross-checked against 187's `q_w` route
          `-2 W_hi d q_w`.  `/0.5963 / 6`.
  `fid`   the SHIPPED call set is prescribed (400 + 2,200 = 2,600), so `min(1, n_pred/n_real)`
          is pinned at 1 and `fid` moves ONLY through the sign precision of the call set,
          which this instrument measures.  `/0.30 / 6` (span 0.28-0.33, conservative end).
  `jac`   `|pred|` pinned at 2,600 and 2,200 of its members are a function of the control
          cells alone, so only the 400-gene head can move the intersection.  Reported both
          on the absolute proxy and rescaled to the board's own raw `jac`; `/0.365 / 6`.
  `reach` ORDERING half only (|truth|-weighted sign precision of the top-89 under key |H|).
          The yield half carries `n_real` and is unreadable on every harness we hold (rule 9),
          so the composite is quoted WITH and WITHOUT an indicative reach term.

    uv run python price228.py
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
SPAN_MSE, DMSE_DCOS = 0.950162566, 2.074495          # 197
SPAN_NMAE, W_HI = 0.5963, 0.0153                     # 187/188
SPAN_FID, SPAN_JAC, SPAN_REACH = 0.30, 0.365, 0.895  # AGENTS.md §4
OURS_COS = 0.09946902740184964                       # 197's inversion of `h1union`
OURS_JAC_RAW, OURS_REACH_RAW = 0.04095656108311888, 0.27654544083789157
K_COS = (1.000, 1.384)                               # NOT identified: h1union / cd4f034
GATE_SHAPE, GATE_AMP = 0.0019, 0.0089


def main() -> None:
    adel = json.loads((WORK / "adel228.json").read_text())
    six = json.loads((WORK / "six228.json").read_text())
    arms = json.loads((WORK / "analyse228.json").read_text())
    A = {(r["fold"], r["arm"]): r for r in adel}
    S = {(r["fold"], r["arm"]): r for r in six}

    def mean(rows, fn):
        return float(np.mean([fn(r) for r in rows]))

    names = [a for a in sorted({r["arm"] for r in adel}) if a != "BASE"]
    out = {}
    print(f"{'arm':<20}{'alpha*':>8}{'rel cos%':>10}{'d_pds':>10}{'d_nmae':>10}"
          f"{'d_prec89':>10}{'d_fid':>9}{'d_jac89':>10}{'pds+':>6}")
    print("-" * 93)
    for arm in names:
        for tag, src in (("adel", A), ("a1", S)):
            if not all((f, arm) in src for f in CANCER):
                continue
            rows = [(src[(f, "BASE")], src[(f, arm)]) for f in CANCER]
            rel = mean(rows, lambda p: 100.0 * (p[1]["cos_agg"] - p[0]["cos_agg"])
                       / abs(p[0]["cos_agg"])) if tag == "adel" else \
                arms["per_variant"][arm]["cancer_rel_cos_agg"]
            dpds = mean(rows, lambda p: p[1]["pds"] - p[0]["pds"]) if tag == "adel" else \
                arms["per_variant"][arm]["cancer_d_pds"]
            dnm = {k: mean(rows, lambda p, k=k: p[1]["nmae_pools"][k] - p[0]["nmae_pools"][k])
                   for k in ("wald89", "lfc89", "h89")}
            dpr = mean(rows, lambda p: p[1]["reach_prec89"] - p[0]["reach_prec89"])
            dfi = mean(rows, lambda p: p[1]["fid_precision"] - p[0]["fid_precision"])
            dja = {k: mean(rows, lambda p, k=k: p[1]["jac"][k] - p[0]["jac"][k])
                   for k in ("89", "200", "400")}
            rja = {k: mean(rows, lambda p, k=k: (p[1]["jac"][k] / p[0]["jac"][k]) - 1.0)
                   for k in ("89", "200", "400")}
            npos = (sum(1 for p in rows if p[1]["pds"] - p[0]["pds"] > 0) if tag == "adel"
                    else arms["per_variant"][arm]["pds_positive_cancer"])
            alpha = float(np.mean([src[(f, arm)].get("alpha_effective", 0.982)
                                   for f in CANCER]))
            dqw = arms["per_variant"][arm]["cancer_d_qw"] if arm in arms["per_variant"] else 0.0

            o_pds = dpds / SPAN_PDS / 6.0 * CALIB
            o_mse = {f"k{k:g}": (DMSE_DCOS * OURS_COS * (rel / 100.0) * k) / SPAN_MSE / 6.0
                     for k in K_COS}
            o_nmae = {k: -v / SPAN_NMAE / 6.0 for k, v in dnm.items()}
            o_nmae_qw = (2.0 * W_HI * dqw) / SPAN_NMAE / 6.0
            o_fid = dfi / SPAN_FID / 6.0
            o_jac_abs = dja["89"] / SPAN_JAC / 6.0
            o_jac_rel = rja["89"] * OURS_JAC_RAW / SPAN_JAC / 6.0
            o_reach = dpr * (OURS_REACH_RAW / 0.62) / SPAN_REACH / 6.0   # indicative only

            key = f"{arm}|{tag}"
            out[key] = dict(
                alpha=alpha, rel_cos_agg_pct=rel, d_pds=dpds, pds_sigma=dpds / SIGMA_PDS,
                d_nmae=dnm, d_prec89=dpr, d_fid_precision=dfi, d_jac=dja, rel_jac=rja,
                pds_positive_cancer=npos,
                kolf2_rel_cos_agg=(100.0 * (A[("kolf2", arm)]["cos_agg"]
                                            - A[("kolf2", "BASE")]["cos_agg"])
                                   / abs(A[("kolf2", "BASE")]["cos_agg"]))
                if ("kolf2", arm) in A else None,
                kolf2_d_pds=(A[("kolf2", arm)]["pds"] - A[("kolf2", "BASE")]["pds"])
                if ("kolf2", arm) in A else None,
                o_pds=o_pds, o_mse=o_mse, o_nmae=o_nmae, o_nmae_qw=o_nmae_qw,
                o_fid=o_fid, o_jac_abs=o_jac_abs, o_jac_rel=o_jac_rel,
                o_reach_indicative=o_reach,
                composite={
                    f"k{k:g}|{pool}|{jm}": (
                        o_pds + o_mse[f"k{k:g}"] + o_nmae[pool] + o_fid
                        + (o_jac_abs if jm == "abs" else o_jac_rel))
                    for k in K_COS for pool in ("wald89", "lfc89", "h89")
                    for jm in ("abs", "rel")},
            )
            out[key]["composite_with_reach"] = {
                k: v + o_reach for k, v in out[key]["composite"].items()}
            if tag == "adel":
                print(f"{arm:<20}{alpha:>8.4f}{rel:>+10.2f}{dpds:>+10.5f}"
                      f"{dnm['wald89']:>+10.5f}{dpr:>+10.5f}{dfi:>+9.5f}"
                      f"{dja['89']:>+10.6f}{npos:>4d}/3")

    print(f"\n{'SIX-MEMBER PRICE (overall points), three cancer folds, @adel':^100}")
    print(f"{'arm':<20}{'o_pds':>10}{'o_mse k1':>10}{'o_mse k1.4':>11}{'o_nmae':>10}"
          f"{'o_fid':>9}{'o_jac':>9}{'TOTAL k1':>11}{'TOTAL k1.4':>11}{'+reach':>10}")
    print("-" * 111)
    for arm in names:
        k = f"{arm}|adel"
        if k not in out:
            continue
        r = out[k]
        c1 = r["composite"]["k1|wald89|rel"]
        c2 = r["composite"]["k1.384|wald89|rel"]
        cr = r["composite_with_reach"]["k1.384|wald89|rel"]
        print(f"{arm:<20}{r['o_pds']:>+10.5f}{r['o_mse']['k1']:>+10.5f}"
              f"{r['o_mse']['k1.384']:>+11.5f}{r['o_nmae']['wald89']:>+10.5f}"
              f"{r['o_fid']:>+9.5f}{r['o_jac_rel']:>+9.5f}{c1:>+11.5f}{c2:>+11.5f}"
              f"{cr:>+10.5f}")

    print(f"\nGate (rule 13, shape-only class, profile-class arm -> reach toll 0 per 169): "
          f"+{GATE_SHAPE}")
    (WORK / "price228.json").write_text(json.dumps(out, indent=2))
    print(f"-> {WORK / 'price228.json'}")


if __name__ == "__main__":
    main()
