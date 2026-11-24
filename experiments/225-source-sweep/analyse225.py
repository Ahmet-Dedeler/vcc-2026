"""225 — reduce `sweep225.json` / `stack225.json` to the reported tables.

PRIMARY endpoint macro cosine (`cos`, per-target unweighted; 204's instrument resolution is
+-0.20% relative).  SECONDARY `pds`.  PRICE: 218 §7.4's rate card through §8.4's CORRECTED
rule, all six members, `k = 1` headline — 223's `price()` verbatim, imported not re-derived.

Headline folds are the three CANCER folds; `kolf2` is excluded exactly as 218 and 223
excluded it (pluripotent recipient, z-score truth with a placeholder `b_ctrl`).

    uv run python analyse225.py
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from harness225 import CANCER, GATE_AMPLITUDE, GATE_SHAPE_ONLY, WORK, price

RES = 0.20      # macro-cosine instrument resolution, % relative (204)


def agg(rows, key):
    """Mean over the cancer folds of every per-fold delta, plus fold-level sign counts."""
    out = {}
    for r in rows:
        if r["fold"] not in CANCER:
            continue
        out.setdefault(key(r), []).append(r)
    res = {}
    for k, v in out.items():
        if len({r["fold"] for r in v}) != len(CANCER):
            continue                        # never average an incomplete fold set
        dcos = float(np.mean([r["d_cos"] for r in v]))
        rel = float(np.mean([r["rel_cos_pct"] for r in v]))
        dpds = float(np.mean([r["d_pds"] for r in v]))
        dqw = float(np.mean([r["d_qw"] for r in v]))
        p = price(dcos, dpds)
        res[k] = dict(n_folds=len(v), rel_cos_pct=rel, d_cos=dcos, d_pds=dpds, d_qw=dqw,
                      cos_pos=sum(r["d_cos"] > 0 for r in v),
                      pds_pos=sum(r["d_pds"] > 0 for r in v),
                      price_k1=p["total"], price_five=p["five_members"], price_pds=p["pds"],
                      carrying=int(np.mean([r["audit"]["carrying_rows"] for r in v])),
                      max_noncarry_pp=float(max(
                          abs(r["rowsplit"].get("non_carrying_rel_pct", 0.0)) for r in v)),
                      per_fold={r["fold"]: dict(rel_cos_pct=r["rel_cos_pct"],
                                                d_pds=r["d_pds"]) for r in v})
    return res


def main() -> None:
    rows = json.loads((WORK / "sweep225.json").read_text())
    if (WORK / "sweep225b.json").exists():
        # the `jiang_5line_noK562` leak-repair run, measured in a separate process; its
        # `dld1_lowmoi` rows are a third independent reproduction of the known answer
        rows += [r for r in json.loads((WORK / "sweep225b.json").read_text())
                 if r.get("source") == "jiang_5line_noK562"]
    arms = [r for r in rows if r["kind"] == "arm"]
    A = agg(arms, lambda r: (r["source"], r["variant"], r["arm"].split("_")[-2], r["g"]))

    print("### Phase 2 — every census source, marginal over base222, three cancer folds")
    print(f"{'source':30s}{'carry':>6s}{'g':>6s}{'relΔcos%':>10s}{'Δpds':>10s}"
          f"{'Δq_w':>10s}{'+cos':>5s}{'+pds':>5s}{'price k=1':>11s}{'×gate':>7s}")
    best = {}
    for (src, var, base, g), v in sorted(A.items()):
        if var or base != "base222":
            continue
        print(f"{src:30s}{v['carrying']:>6d}{g:>6.2f}{v['rel_cos_pct']:>+10.3f}"
              f"{v['d_pds']:>+10.5f}{v['d_qw']:>+10.5f}{v['cos_pos']:>4d}/3"
              f"{v['pds_pos']:>4d}/3{v['price_k1']:>+11.5f}"
              f"{v['price_k1'] / GATE_SHAPE_ONLY:>7.2f}")
        cur = best.get(src)
        if cur is None or v["price_k1"] > cur[1]["price_k1"]:
            best[src] = (g, v)

    print("\n### Phase 3 — falsifiers at each source's own price-argmax dose")
    print(f"{'source':30s}{'g*':>6s}{'real':>9s}{'SHUF':>9s}{'SHUF2':>9s}"
          f"{'CD4ROWS':>9s}{'ratio':>8s}{'F1':>5s}{'F2':>5s}{'F3':>5s}")
    fals = {}
    for src, (g, v) in sorted(best.items()):
        row = dict(g=g, real=v["rel_cos_pct"], price=v["price_k1"], d_pds=v["d_pds"],
                   pds_pos=v["pds_pos"], cos_pos=v["cos_pos"],
                   max_noncarry_pp=v["max_noncarry_pp"])
        for var, tag in (("SHUF", "shuf"), ("SHUF2", "shuf2"), ("CD4ROWS", "cd4rows")):
            w = A.get((src, var, "base222", g))
            row[tag] = w["rel_cos_pct"] if w else None
            if var == "CD4ROWS" and w:
                row["cd4rows_price"] = w["price_k1"]
        row["ratio_vs_cd4rows"] = (row["real"] / row["cd4rows"]
                                   if row.get("cd4rows") not in (None, 0.0) else None)
        row["F1"] = (row["real"] > max(x for x in (row["shuf"], row["shuf2"])
                                       if x is not None)) and row["real"] > 0
        row["F2"] = row["max_noncarry_pp"] <= 0.018
        # F3 is "beat the FREE alternative", not "beat it by a positive ratio": when the
        # row-matched CD4 control is itself negative the ratio flips sign and is meaningless,
        # so the substantive comparison is the difference.  (223 reported ratios because
        # every one of its controls happened to be positive.)
        row["F3"] = bool(row["cd4rows"] is not None and row["real"] > row["cd4rows"]
                         and row["real"] > 0)
        fals[src] = row
        rr = row["ratio_vs_cd4rows"]
        print(f"{src:30s}{g:>6.2f}{row['real']:>+9.3f}"
              f"{(row['shuf'] or 0):>+9.3f}{(row['shuf2'] or 0):>+9.3f}"
              f"{(row['cd4rows'] or 0):>+9.3f}"
              f"{(f'{rr:.2f}x' if rr else '   -'):>8s}"
              f"{'  ok' if row['F1'] else ' FAIL':>5s}"
              f"{'  ok' if row['F2'] else ' FAIL':>5s}"
              f"{'  ok' if row['F3'] else ' FAIL':>5s}")

    print("\n### the same sources on `cd4@0.34` alone (the pre-222 frame)")
    for src, (g, _) in sorted(best.items()):
        v = A.get((src, "", "cd4", g))
        if v:
            print(f"{src:30s} g={g:.2f}  rel {v['rel_cos_pct']:+.3f}%  "
                  f"dpds {v['d_pds']:+.5f}  price {v['price_k1']:+.5f}")

    out = dict(per_source=[dict(source=s, variant=v, base=b, g=g, **d)
                           for (s, v, b, g), d in A.items()],
               falsifiers=fals, gate_shape_only=GATE_SHAPE_ONLY,
               gate_amplitude=GATE_AMPLITUDE, resolution_rel_pct=RES)

    sp = WORK / "stack225.json"
    if sp.exists():
        srows = json.loads(sp.read_text())
        S = agg(srows, lambda r: (r["tag"], r["variant"], r["g"]))
        print("\n### Phase 4 — the stack curve (marginal over base222, cancer folds)")
        print(f"{'stack':26s}{'n':>3s}{'g':>6s}{'carry':>6s}{'relΔcos%':>10s}"
              f"{'Δpds':>10s}{'price':>10s}{'×dld1':>7s}{'×gate':>7s}")
        dld1 = max((v["price_k1"] for (t, vr, g), v in S.items()
                    if t == "CUM1" and not vr), default=None)
        bests = {}
        for (tag, var, g), v in S.items():
            if var:
                continue
            cur = bests.get(tag)
            if cur is None or v["price_k1"] > cur[1]["price_k1"]:
                bests[tag] = (g, v)
        order = (sorted([t for t in bests if t.startswith("CUM")],
                        key=lambda t: int(t[3:]))
                 + sorted(t for t in bests if t.startswith("LOO"))
                 + [t for t in ("ALLSMALL", "ALLSMALL-nodld1") if t in bests])
        for tag in order:
            g, v = bests[tag]
            n = next((r["n_members"] for r in srows if r["tag"] == tag), 0)
            print(f"{tag:26s}{n:>3d}{g:>6.2f}{v['carrying']:>6d}"
                  f"{v['rel_cos_pct']:>+10.3f}{v['d_pds']:>+10.5f}{v['price_k1']:>+10.5f}"
                  f"{(v['price_k1'] / dld1 if dld1 else 0):>7.2f}"
                  f"{v['price_k1'] / GATE_SHAPE_ONLY:>7.2f}")
        out["stack"] = {t: dict(g=g, **v) for t, (g, v) in bests.items()}
        out["stack_all"] = [dict(tag=t, variant=vr, g=g, **d) for (t, vr, g), d in S.items()]

    (WORK / "analyse225.json").write_text(json.dumps(out, indent=2))
    print(f"\n-> {WORK / 'analyse225.json'}")


if __name__ == "__main__":
    main()
