"""246 — two small derived tables from the sweep JSONs (no fold load):

  1. MECHANISM: the D family (deflated self-residual on the BARE whitened base) against `f0`,
     i.e. does deflation keep 099/218's self-residual cosine gain and remove its `pds` cost
     when the relaxation is NOT already in the stack?
  2. EXCHANGE RATE: `pds` bought per unit of aggregate cosine given up, on the DOSE axis
     (216's lever, `DOSE_cd4@f`) versus the DEFLATION axis (`B2_cd4_svd<k>`), same folds,
     same base, same instrument; and the same rate for the energy-matched random controls.

    uv run python mech246.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

WORK = Path(__file__).resolve().parent / "work"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")


def load():
    rows = {}
    for held in CANCER + ("kolf2",):
        for tag in ("sweep246", "sweep246b", "sweep246c"):
            p = WORK / f"{tag}-{held}.json"
            if p.exists():
                for r in json.loads(p.read_text())["rows"]:
                    rows[(held, r["arm"])] = r
    return rows


def table(rows, arms, base, folds, title):
    print("\n" + "=" * 100 + f"\n{title}\n" + "=" * 100)
    print(f"{'arm':26s} {'rel agg%':>9s} {'rel macro%':>10s} {'d pds':>10s} {'d q_w':>10s} "
          f"{'pds/cos':>9s}")
    out = {}
    for a in arms:
        sel = [(rows[(f, a)], rows[(f, base)]) for f in folds if (f, a) in rows]
        if not sel:
            continue
        ra = np.mean([r["cos_agg"] / b["cos_agg"] - 1 for r, b in sel])
        rm = np.mean([r["a23_cos"] / b["a23_cos"] - 1 for r, b in sel])
        dp = np.mean([r["pds"] - b["pds"] for r, b in sel])
        dq = np.mean([r["a15_sgn_w"] - b["a15_sgn_w"] for r, b in sel])
        rate = dp / ra if abs(ra) > 1e-9 else float("nan")
        out[a] = dict(rel_agg=float(ra), rel_macro=float(rm), d_pds=float(dp), d_qw=float(dq),
                      pds_per_unit_rel_cos=float(rate), n=len(sel))
        print(f"{a:26s} {100*ra:+9.3f} {100*rm:+10.3f} {dp:+10.6f} {dq:+10.6f} {rate:+9.4f}")
    return out


def main():
    rows = load()
    res = {}
    res["D_vs_f0_cancer"] = table(
        rows, ["self0.1_bare"] + [f"D_{m}@0.1" for m in
                                  ("none", "centre", "svd1", "svd2", "svd5", "svd10", "svd20",
                                   "randE5", "randE5b", "svd5_SHUF")],
        "f0", CANCER, "MECHANISM: deflated SELF-residual on the BARE base, vs f0 (3 cancer folds)")
    res["D_vs_f0_kolf2"] = table(
        rows, [f"D_{m}@0.1" for m in ("none", "svd2", "svd5", "svd10", "randE5")],
        "f0", ("kolf2",), "same, kolf2 apart")
    res["exchange_cancer"] = table(
        rows, ["DOSE_cd4@0.28", "DOSE_cd4@0.2"] + [f"B2_cd4_svd{k}" for k in (1, 2, 3, 5, 7, 10)]
        + [f"B2_cd4_svd{k}_resolved" for k in (3, 5)]
        + [f"B2_cd4_randE{k}{s}" for k in (3, 5, 10) for s in ("", "b")],
        "STACK", CANCER,
        "EXCHANGE RATE vs STACK (3 cancer folds): d pds per unit relative aggregate cosine given up")
    (WORK / "mech246.json").write_text(json.dumps(res, indent=1))
    print(f"\nwrote {WORK/'mech246.json'}")


if __name__ == "__main__":
    main()
