"""223 Phase 2/3 — reduce `sweep223.json` to the reported table, all four falsifiers, and
the price under AGENTS.md 8.4's CORRECTED rule.

FRAME.  Everything is a delta against `base222` = `cd4@0.34 + h1pair@0.28` on the same fold,
so every number is the MARGINAL value of the source on top of the arm 222 is building.
The `kolf2` fold is reported but EXCLUDED from every headline, exactly as 218 excluded it:
its recipient is pluripotent, which flatters any pluripotent source, and its truth is a
z-score matrix with a placeholder `b_ctrl`.

PRIMARY endpoint: macro cosine (`cos`, per-target unweighted), instrument resolution
+-0.20% relative (204).  SECONDARY: `pds`.  THIRD: `q_w` (`a15_sgn_w`).

PRICE.  218 section 7.4's rate card, from `cd4f034`'s own board row divided by its measured
cosine move -- +0.007259 overall per +0.01 ABSOLUTE cosine, all six members, with `pds`
also priced independently through 207's route and doubled per 8.4.  `k = 1` (uncalibrated)
is the headline because the PLOLO->board transfer factor was measured on arms that move
~93% of panel rows and these move 5-42%.

    uv run python analyse223.py
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
G = (0.17, 0.28, 0.45)
# 218 section 7.4, from cd4f034's board row / its measured cosine move (+0.017990)
RATE = {"mse": 0.003634, "reach": 0.002446, "nmae": 0.000980, "fid": 0.000847,
        "jac": 0.000100, "pds_board": -0.000748}
RATE_TOTAL_EX_PDS = RATE["mse"] + RATE["reach"] + RATE["nmae"] + RATE["fid"] + RATE["jac"]
PDS_SPAN = 0.452          # 207's route: raw pds delta / span / 6 members * 1.011, x2.0 (8.4)


def price(d_cos_abs: float, d_pds_raw: float, k: float = 1.0) -> dict:
    """All six members. Never carry four at zero (8.4 error 1)."""
    cos_gain = d_cos_abs * k
    five = RATE_TOTAL_EX_PDS * (cos_gain / 0.01)
    pds_local = (d_pds_raw / PDS_SPAN) / 6.0 * 1.011 * 2.0
    return dict(five_members=five, pds=pds_local, total=five + pds_local)


def main() -> None:
    rows = json.loads((WORK / "sweep223.json").read_text())
    base = {(r["fold"], r["arm"]): r for r in rows if r["kind"] == "baseline"}
    arms = [r for r in rows if r["kind"] == "arm"]
    folds = sorted({r["fold"] for r in rows})

    print("### baselines (macro cos / cos_agg / pds), per fold")
    print(f"{'fold':16s}" + "".join(f"{a:>34s}" for a in ("f0", "cd4", "base222")))
    for f in folds:
        line = f"{f:16s}"
        for a in ("f0", "cd4", "base222"):
            b = base.get((f, a))
            line += ("".rjust(34) if b is None else
                     f"{b['cos']:>12.6f}{b['cos_agg']:>10.6f}{b['pds']:>12.6f}")
        print(line)

    # ------------------------------------------------------------------ index
    idx = defaultdict(dict)
    for r in arms:
        idx[(r["source"], r["variant"], r["rule"], r["g"])][r["fold"]] = r

    sources = sorted({r["source"] for r in arms})

    def agg(src, var, rule, g, folds_use):
        got = [idx[(src, var, rule, g)].get(f) for f in folds_use]
        got = [x for x in got if x is not None]
        if not got:
            return None
        return dict(
            n=len(got),
            rel_cos=float(np.mean([x["rel_cos_pct"] for x in got])),
            d_cos=float(np.mean([x["d_cos"] for x in got])),
            rel_cos_agg=float(np.mean([x["rel_cos_agg_pct"] for x in got])),
            d_pds=float(np.mean([x["d_pds"] for x in got])),
            d_qw=float(np.mean([x["d_qw"] for x in got])),
            pos_cos=int(sum(x["d_cos"] > 0 for x in got)),
            pos_pds=int(sum(x["d_pds"] > 0 for x in got)),
            frac=float(np.mean([x["audit"]["fraction"] for x in got])),
            carry=float(np.mean([x["audit"]["carrying_rows"] for x in got])),
            cov=float(np.mean([x["audit"]["coverage_frac"] for x in got])),
            noncarry_pp=float(np.max([abs(x["rowsplit"].get("non_carrying_rel_pct", 0.0))
                                      for x in got])),
            carry_rel=float(np.mean([x["rowsplit"].get("carrying_rel_pct", 0.0)
                                     for x in got])),
        )

    print("\n### PHASE 2 — marginal value on top of `cd4@0.34 + h1pair@0.28`, "
          "CARRYING-ROW sizing, THREE CANCER FOLDS")
    hdr = (f"{'source':18s}{'g':>6s}{'cov':>7s}{'frac':>8s}"
           f"{'rel dcos%':>11s}{'dcos':>11s}{'dpds':>10s}{'dq_w':>10s}"
           f"{'+cos':>5s}{'+pds':>5s}{'price k=1':>11s}")
    print(hdr)
    print("-" * len(hdr))
    table = []
    for src in sources:
        for g in G:
            a = agg(src, "", "CARRY", g, CANCER)
            if a is None:
                continue
            p = price(a["d_cos"], a["d_pds"])
            table.append(dict(source=src, g=g, **a, price_k1=p["total"],
                              price_five=p["five_members"], price_pds=p["pds"]))
            print(f"{src:18s}{g:>6.2f}{a['cov']:>7.1%}{a['frac']:>8.4f}"
                  f"{a['rel_cos']:>11.3f}{a['d_cos']:>11.2e}{a['d_pds']:>10.5f}"
                  f"{a['d_qw']:>10.5f}{a['pos_cos']:>4d}/3{a['pos_pds']:>4d}/3"
                  f"{p['total']:>11.5f}")

    print("\n### PHASE 3 FALSIFIERS (three cancer folds)")
    fh = (f"{'source':18s}{'g':>6s}{'real rel%':>11s}{'SHUF rel%':>11s}"
          f"{'CD4rows rel%':>14s}{'ratio real/CD4':>16s}"
          f"{'non-carry |pp|':>16s}{'carry rel%':>12s}")
    print(fh)
    print("-" * len(fh))
    fal = []
    for src in sources:
        for g in G:
            a = agg(src, "", "CARRY", g, CANCER)
            s = agg(src, "SHUF", "CARRY", g, CANCER)
            c = agg(src, "CD4ROWS", "CARRY", g, CANCER)
            if a is None:
                continue
            ratio = (a["rel_cos"] / c["rel_cos"]) if c and abs(c["rel_cos"]) > 1e-9 else None
            fal.append(dict(source=src, g=g, real=a["rel_cos"],
                            shuf=(s["rel_cos"] if s else None),
                            cd4rows=(c["rel_cos"] if c else None), ratio=ratio,
                            noncarry_pp=a["noncarry_pp"], carry_rel=a["carry_rel"],
                            F1_beats_shuffle=bool(s and a["rel_cos"] > s["rel_cos"]),
                            F2_noncarry_quiet=bool(a["noncarry_pp"] <= 0.05),
                            F3_beats_cd4rows=bool(c and a["rel_cos"] > c["rel_cos"])))
            print(f"{src:18s}{g:>6.2f}{a['rel_cos']:>11.3f}"
                  f"{(s['rel_cos'] if s else float('nan')):>11.3f}"
                  f"{(c['rel_cos'] if c else float('nan')):>14.3f}"
                  f"{(f'{ratio:.2f}' if ratio else 'n/a'):>16s}"
                  f"{a['noncarry_pp']:>16.4f}{a['carry_rel']:>12.3f}")

    print("\n### THE SILENT ZERO, in this experiment's own numbers "
          "(OLD median rule, all four folds)")
    oh = (f"{'source':18s}{'coverage':>10s}{'OLD frac':>10s}"
          f"{'bit-identical to base222':>26s}{'CARRY frac':>12s}")
    print(oh)
    print("-" * len(oh))
    void = []
    for src in sources:
        got = [idx[(src, "", "OLD", 0.28)].get(f) for f in folds]
        got = [x for x in got if x is not None]
        gc = [idx[(src, "", "CARRY", 0.28)].get(f) for f in folds]
        gc = [x for x in gc if x is not None]
        if not got:
            continue
        allbit = all(x["bit_identical_to_base222"] for x in got)
        void.append(dict(source=src, coverage=float(np.mean([x["audit"]["coverage_frac"]
                                                            for x in got])),
                         old_frac=float(np.mean([x["audit"]["fraction"] for x in got])),
                         bit_identical=allbit,
                         carry_frac=float(np.mean([x["audit"]["fraction"] for x in gc]))))
        print(f"{src:18s}{np.mean([x['audit']['coverage_frac'] for x in got]):>10.1%}"
              f"{np.mean([x['audit']['fraction'] for x in got]):>10.5f}"
              f"{str(allbit):>26s}{np.mean([x['audit']['fraction'] for x in gc]):>12.5f}")

    print("\n### FOUR-FOLD view (incl. the excluded kolf2 fold), CARRY, g = 0.28")
    print(f"{'source':18s}" + "".join(f"{f[:12]:>14s}" for f in folds) + f"{'4-fold':>10s}")
    for src in sources:
        line = f"{src:18s}"
        vals = []
        for f in folds:
            r = idx[(src, "", "CARRY", 0.28)].get(f)
            line += (f"{r['rel_cos_pct']:>14.3f}" if r else "".rjust(14))
            if r:
                vals.append(r["rel_cos_pct"])
        print(line + f"{np.mean(vals):>10.3f}")

    out = dict(table=table, falsifiers=fal, void=void,
               rate_card=RATE, gate_shape_only=0.0019, gate_amplitude=0.0089)
    (WORK / "analyse223.json").write_text(json.dumps(out, indent=2))
    print(f"\n-> {WORK / 'analyse223.json'}")


if __name__ == "__main__":
    main()
