"""302 — two more POST-HOC, LABELLED readings from the same sums: the isoform (two-construct)
ceiling as a function of depth, and the disattenuated mirror cosine per reliability stratum.

    PYTHONDONTWRITEBYTECODE=1 uv run python strat302.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib302 as L  # noqa: E402

RHO_KD = 0.243   # the H1 Flex reference's split-half residual reliability (273 §0.6) — an
                 # ASSUMPTION: the gse295214_h1 bank's own per-target reliability is not measured.
LOG: list[str] = []


def P(m: str = "") -> None:
    print(m, flush=True)
    LOG.append(m)


def main() -> None:
    rows = list(csv.DictReader(open(L.WORK / "symmetry302.csv")))
    pri = np.asarray([r["cohort"] == "primary" for r in rows])
    sym = np.asarray([r["symbol"] for r in rows])
    S = np.asarray([float(r["S_g"]) for r in rows])[pri]
    rel = np.asarray([float(r["orf_reliability"]) for r in rows])[pri]
    cells = np.asarray([float(r["orf_cells"]) for r in rows])[pri]
    sym = sym[pri]
    out: dict = {}

    def sb(r):        # Spearman–Brown: half-sample reliability -> full-sample reliability
        return 2 * r / (1 + r)

    P(f"disattenuation assumes rho_CRISPRi = {RHO_KD} (the H1 Flex reference, 273 §0.6); the "
      f"gse295214_h1 bank's own per-target reliability is NOT measured here.")
    P("\nstratum                     n   rel_half  rho_ORF  ceiling  median S    implied true cos  "
      "[bootstrap 95 %]")
    strata = {}
    ord_r = np.argsort(rel)
    t = len(S) // 3
    groups = [("all 148", np.arange(len(S))),
              ("reliability T1", ord_r[:t]), ("reliability T2", ord_r[t:2 * t]),
              ("reliability T3", ord_r[2 * t:]),
              ("rel >= 0.10", np.flatnonzero(rel >= 0.10)),
              ("rel >= 0.20", np.flatnonzero(rel >= 0.20)),
              ("rel >= 0.30", np.flatnonzero(rel >= 0.30))]
    for name, ii in groups:
        if ii.size < 5:
            continue
        v, rr = S[ii], float(np.median(rel[ii]))
        ceil = float(np.sqrt(sb(max(rr, 1e-6)) * RHO_KD))
        lo, hi, _ = L.boot_ci(lambda k: float(np.median(v[k])), len(v), 20000, L.SEED)
        strata[name] = dict(n=int(v.size), rel_half=rr, rho_orf=float(sb(rr)), ceiling=ceil,
                            median=float(np.median(v)), ci=[lo, hi],
                            implied=float(np.median(v) / ceil),
                            implied_ci=[lo / ceil, hi / ceil])
        P(f"{name:24s} {v.size:4d}   {rr:+.4f}  {sb(rr):+.4f}  {ceil:.4f}  {np.median(v):+.4f}   "
          f"{np.median(v) / ceil:+.4f}  [{lo / ceil:+.4f}, {hi / ceil:+.4f}]")
    out["strata_disattenuated"] = strata
    out["rho_kd_assumed"] = RHO_KD

    # ---- the two-construct (isoform) ceiling by depth ---------------------------------------
    iso = list(csv.DictReader(open(L.WORK / "isoform302.csv")))
    st = json.load(open(L.WORK / "stats302.json"))
    P(f"\nper-isoform mirror rows: {len(iso)} over {len({r['symbol'] for r in iso})} TFs; "
      f"median S {np.median([float(r['S']) for r in iso]):+.4f} "
      f"(pooled-symbol median {np.median(S):+.4f})")
    P(f"isoform agreement (two ORF constructs of one TF, same cells, same batch): median cos "
      f"{st['secondaries']['per_isoform']['isoform_agreement_median']:+.4f} over "
      f"{st['secondaries']['per_isoform']['n_pairs']} pairs")

    # per-TF: does the mirror statistic agree between a TF's own isoforms?
    by_sym: dict[str, list[float]] = {}
    for r in iso:
        by_sym.setdefault(r["symbol"], []).append(float(r["S"]))
    agree = [(np.std(v), np.mean(v)) for v in by_sym.values() if len(v) >= 2]
    P(f"per-TF spread of the isoform-level S: median sd {np.median([a for a, _ in agree]):.4f} "
      f"over {len(agree)} TFs — for comparison the whole cohort's sd is {np.std(S):.4f}")
    out["isoform_S_spread"] = dict(median_sd=float(np.median([a for a, _ in agree])),
                                   cohort_sd=float(np.std(S)), n_tfs=len(agree))

    # ---- what a mirror WOULD look like: the best-measured single TFs -------------------------
    best = np.argsort(-rel)[:10]
    P("\nthe 10 best-measured TFs (post-hoc, labelled):")
    for i in best:
        P(f"   {sym[i]:10s} rel {rel[i]:+.3f}  cells {cells[i]:6.0f}  S {S[i]:+.4f}  "
          f"implied true cos {S[i] / np.sqrt(sb(max(rel[i], 1e-6)) * RHO_KD):+.4f}")
    out["best_measured"] = [(sym[i], float(rel[i]), float(cells[i]), float(S[i])) for i in best]
    out["log"] = LOG
    (L.WORK / "strat302.json").write_text(json.dumps(out, indent=1, default=float))
    P("\n-> work/strat302.json")


if __name__ == "__main__":
    sys.exit(main())
