"""302 — descriptive tails, POST-HOC and LABELLED.  Reads only the sums access302.py wrote.

Every number here is descriptive: it is not part of C6a/C6b/C6c and was not preregistered as a
rule.  Its job is to say WHY the three verdict cells read what they read.

    PYTHONDONTWRITEBYTECODE=1 uv run python desc302.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib302 as L  # noqa: E402
from analyse302 import cos_pair  # noqa: E402

LOG: list[str] = []


def P(m: str = "") -> None:
    print(m, flush=True)
    LOG.append(m)


def main() -> None:
    rows = list(csv.DictReader(open(L.WORK / "symmetry302.csv")))
    pri = np.asarray([r["cohort"] == "primary" for r in rows])
    sym = np.asarray([r["symbol"] for r in rows])
    S = np.asarray([float(r["S_g"]) for r in rows])
    rel = np.asarray([float(r["orf_reliability"]) for r in rows])
    cells = np.asarray([float(r["orf_cells"]) for r in rows])
    umi = np.asarray([float(r["orf_umi_per_cell"]) for r in rows])
    C_held = np.asarray([float(r["C_held"]) if r["C_held"] else np.nan for r in rows])
    Sp, relp, cellsp = S[pri], rel[pri], cells[pri]
    out: dict = {}

    q = [0, 5, 10, 25, 50, 75, 90, 95, 100]
    P("S_g on the 148-TF primary cohort, percentiles "
      + ", ".join(f"p{a}={np.percentile(Sp, a):+.4f}" for a in q))
    out["S_percentiles"] = {f"p{a}": float(np.percentile(Sp, a)) for a in q}

    # --- the measurement ceiling ------------------------------------------------------------
    r_half = float(np.median(relp))
    sb = 2 * r_half / (1 + r_half)
    P(f"\nORF split-half reliability (half-size samples): median {r_half:+.4f} "
      f"→ Spearman–Brown full-sample reliability {sb:+.4f}")
    P("attenuation ceiling  max|S| = sqrt(rho_ORF · rho_CRISPRi):")
    ceil = {}
    for rk in (0.10, 0.243, 0.30, 0.50, 1.00):
        ceil[str(rk)] = float(np.sqrt(sb * rk))
        P(f"   rho_CRISPRi = {rk:.3f}  ->  max|S| = {np.sqrt(sb * rk):.4f}  "
          f"(observed median {np.median(Sp):+.4f} implies a true mirror cosine "
          f"{np.median(Sp) / np.sqrt(sb * rk):+.4f})")
    out["reliability"] = dict(median_half=r_half, spearman_brown=sb, ceiling=ceil,
                              q10=float(np.percentile(relp, 10)),
                              q90=float(np.percentile(relp, 90)),
                              share_le0=float((relp <= 0).mean()))
    P(f"   {(relp <= 0).mean():.1%} of the 148 TFs have a split-half reliability <= 0")

    # --- reliability / depth strata (post-hoc) ----------------------------------------------
    P("\nstrata (post-hoc, labelled) — median S with a 20,000-draw paired TF bootstrap:")
    strata = {}
    ord_rel = np.argsort(relp)
    t = len(Sp) // 3
    groups = [("reliability T1 (lowest)", ord_rel[:t]), ("reliability T2", ord_rel[t:2 * t]),
              ("reliability T3 (highest)", ord_rel[2 * t:]),
              ("rel >= 0.10", np.flatnonzero(relp >= 0.10)),
              ("rel >= 0.20", np.flatnonzero(relp >= 0.20)),
              ("cells >= 1000", np.flatnonzero(cellsp >= 1000)),
              ("cells >= 2000", np.flatnonzero(cellsp >= 2000))]
    for name, ii in groups:
        if ii.size < 5:
            continue
        v = Sp[ii]
        lo, hi, _ = L.boot_ci(lambda k: float(np.median(v[k])), len(v), 20000, L.SEED)
        strata[name] = dict(n=int(v.size), median=float(np.median(v)), ci=[lo, hi],
                            rel_median=float(np.median(relp[ii])),
                            cells_median=float(np.median(cellsp[ii])),
                            share_neg=float((v < 0).mean()))
        P(f"   {name:26s} n {v.size:3d}  rel {np.median(relp[ii]):+.3f}  cells "
          f"{np.median(cellsp[ii]):6.0f}  median S {np.median(v):+.4f} [{lo:+.4f}, {hi:+.4f}]  "
          f"{(v < 0).mean():.0%} negative")
    out["strata"] = strata

    # --- the null's own tail shares ----------------------------------------------------------
    z = np.load(L.WORK / "sums302.npz", allow_pickle=True)
    st = json.load(open(L.WORK / "stats302.json"))
    P(f"\nnull (200 TF-shuffled pairings): median {st['C6a']['null_median']:+.4f} ± "
      f"{st['C6a']['null_sd']:.4f}; observed {st['C6a']['median']:+.4f} "
      f"(z {st['C6a']['z']:+.2f}, {st['C6a']['null_frac_le_obs']:.0%} of draws <= observed)")

    # --- correlations of |S| and S with the covariates ---------------------------------------
    P("\n|S| against the covariates (post-hoc):")
    for name, v in [("reliability", relp), ("cells", cellsp), ("umi/cell", umi[pri]),
                    ("C_held", C_held[pri])]:
        P(f"   rho(|S|, {name:12s}) {L.spearman(v, np.abs(Sp)):+.4f}   "
          f"rho(S, {name:12s}) {L.spearman(v, Sp):+.4f}")

    # --- the tails ---------------------------------------------------------------------------
    o = np.argsort(Sp)
    P("\nmost MIRRORED (S most negative; post-hoc, labelled):")
    for i in o[:20]:
        P(f"   {sym[pri][i]:10s} S {Sp[i]:+.4f}  cells {cellsp[i]:6.0f}  rel {relp[i]:+.3f}"
          f"  C_held {C_held[pri][i]:+.4f}" if np.isfinite(C_held[pri][i]) else
          f"   {sym[pri][i]:10s} S {Sp[i]:+.4f}  cells {cellsp[i]:6.0f}  rel {relp[i]:+.3f}")
    P("most ALIGNED (S most positive; post-hoc, labelled):")
    for i in o[::-1][:20]:
        P(f"   {sym[pri][i]:10s} S {Sp[i]:+.4f}  cells {cellsp[i]:6.0f}  rel {relp[i]:+.3f}"
          f"  C_held {C_held[pri][i]:+.4f}" if np.isfinite(C_held[pri][i]) else
          f"   {sym[pri][i]:10s} S {Sp[i]:+.4f}  cells {cellsp[i]:6.0f}  rel {relp[i]:+.3f}")
    out["top_mirrored"] = [(sym[pri][i], float(Sp[i]), float(cellsp[i]), float(relp[i]))
                           for i in o[:20]]
    out["top_aligned"] = [(sym[pri][i], float(Sp[i]), float(cellsp[i]), float(relp[i]))
                          for i in o[::-1][:20]]

    # --- isoform agreement vs the mirror statistic -------------------------------------------
    iso = list(csv.DictReader(open(L.WORK / "isoform302.csv")))
    P(f"\nper-isoform rows {len(iso)} over {len({r['symbol'] for r in iso})} TFs "
      f"(secondary, declared): median S {np.median([float(r['S']) for r in iso]):+.4f}; "
      f"isoform agreement median cos {st['secondaries']['per_isoform']['isoform_agreement_median']:+.4f} "
      f"over {st['secondaries']['per_isoform']['n_pairs']} pairs — the two-construct ceiling "
      f"for this library, measured in the same cells")
    out["isoform"] = st["secondaries"]["per_isoform"]

    (L.WORK / "desc302.json").write_text(json.dumps(out, indent=1, default=float))
    P("\n-> work/desc302.json")


if __name__ == "__main__":
    sys.exit(main())
