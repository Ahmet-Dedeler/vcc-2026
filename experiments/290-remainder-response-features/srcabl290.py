"""290 step 6 — POST-HOC, labelled: is the half-1 signal just the two Orion banks?

Half 1's 14 sources are wildly unequal in coverage (orion_hct116 16,594 rows and orion_hek293t
16,721 vs 23 for tian2019_ipsc), so the aggregated feature block could be a property of Orion
rather than of responses in general.  Re-aggregates the frozen 18 columns over source SUBSETS from
`work/persrc290.npz` (no bank is re-read) and refits, with the same folds and target.

Known-answer row: the full half-1 subset must reproduce the frozen +0.3174 at 0.0e+00.

    PYTHONDONTWRITEBYTECODE=1 uv run python srcabl290.py

Writes work/srcabl290.json, work/srcabl290.log.
"""
from __future__ import annotations

import json
import resource
import time

import numpy as np

import lib290 as L
from rep290 import aggregate

T0 = time.time()
LOG: list[str] = []


def P(msg: str = "") -> None:
    line = f"[{time.time() - T0:7.1f}s] {msg}" if msg else ""
    print(line, flush=True)
    LOG.append(line)


def main() -> None:
    h = L.check_prereg()
    P(f"prereg {h} MATCH")

    z = np.load(L.WORK / "persrc290.npz", allow_pickle=True)
    keys = np.asarray(z["keys"]).astype(str)
    on_axis = np.asarray(z["on_axis"])
    Q = {k: np.asarray(z[k]) for k in
         (["mag", "breadth", "gen_cos", "pc5", "pc20"] + L.PE_COLS
          + ["proj_prolif", "proj_stress", "ontarget"])}
    h1 = np.asarray(z["half1"])

    comp = np.load(L.SPLITS284, allow_pickle=True)["comp"]
    fold = L.component_folds(comp, 5)
    ref = np.load(L.REM284, allow_pickle=True)
    u = np.flatnonzero(np.asarray(ref["ok"]))
    y = np.asarray(ref["r2"])[u]
    r1 = np.asarray(ref["r1"])[u]
    foldu, compu = fold[u], comp[u]
    P(f"unit {len(u)} genes; 284's remainder rho {L.spearman(r1, y):+.4f} (KA: +0.3277)")

    n_rows = {k: int(np.isfinite(Q["mag"][i]).sum()) for i, k in enumerate(keys)}
    orion = np.array([i for i in h1 if keys[i].startswith("orion_")])
    big = np.array([i for i in h1 if n_rows[keys[i]] >= 1000])
    small = np.array([i for i in h1 if n_rows[keys[i]] < 1000])

    SUBSETS = {"half 1, all 14 sources (frozen)": h1,
               "the 2 Orion banks only": orion,
               "half 1 WITHOUT the 2 Orion banks": np.array([i for i in h1 if i not in set(orion)]),
               "half 1, the 6 banks with >= 1000 rows": big,
               "half 1, the 8 banks with < 1000 rows": small}
    for i in h1:
        SUBSETS[f"single source: {keys[i]}"] = np.array([i])

    out = {"prereg_sha256": h, "n_rows": n_rows, "cells": {}}
    for name, rows in SUBSETS.items():
        X = aggregate(Q, rows, on_axis)[u]
        have = np.isfinite(X[:, [L.FEATS290.index(c) for c in
                                 ("logmag", "breadth", "gen_cos")]]).all(axis=1)
        if have.sum() < 200:
            P(f"  {name:44s} only {int(have.sum())} of the unit covered — skipped")
            continue
        p, _ = L.oof_ridge(X[have], y[have], foldu[have])
        b = L.boot_rho(p, y[have], compu[have], np.random.default_rng(L.SEED + 161))
        j = L.FEATS290.index
        cell = dict(n_sources=len(rows), n_covered=int(have.sum()), rho=b,
                    gen_pc5=L.spearman(X[have, j("gen_pc5")], y[have]),
                    gen_pc20=L.spearman(X[have, j("gen_pc20")], y[have]),
                    breadth=L.spearman(X[have, j("breadth")], y[have]),
                    logmag=L.spearman(X[have, j("logmag")], y[have]),
                    pe_prolif_abs=L.spearman(X[have, j("pe_prolif_abs")], y[have]),
                    remainder_half1_same_genes=L.spearman(r1[have], y[have]))
        out["cells"][name] = cell
        P(f"  {name:44s} n_src {len(rows):2d}  covered {int(have.sum()):5d}  model "
          f"{b['point']:+.4f} [{b['lo']:+.4f}, {b['hi']:+.4f}]  pc5 {cell['gen_pc5']:+.4f}  "
          f"pc20 {cell['gen_pc20']:+.4f}  breadth {cell['breadth']:+.4f}  "
          f"logmag {cell['logmag']:+.4f}  (r1 on the same genes {cell['remainder_half1_same_genes']:+.4f})")

    frozen = out["cells"]["half 1, all 14 sources (frozen)"]["rho"]["point"]
    out["ka_frozen_reproduced"] = dict(value=frozen, want=0.3174,
                                       pass_=bool(abs(frozen - 0.31738512466156) < 1e-9))
    P(f"\nKA (internal): the full half-1 subset reproduces the frozen cell {frozen:+.6f} "
      f"-> {'PASS' if out['ka_frozen_reproduced']['pass_'] else 'FAIL'}")

    out["seconds"] = time.time() - T0
    out["peak_rss_gib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30
    (L.WORK / "srcabl290.json").write_text(json.dumps(out, indent=1, default=float))
    (L.WORK / "srcabl290.log").write_text("\n".join(LOG))
    P(f"peak RSS {out['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main()
