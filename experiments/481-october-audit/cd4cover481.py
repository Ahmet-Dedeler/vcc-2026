"""481 F1: which panel targets the CD4 atlas covers, per condition, split by K562-essential status (obs/var only, < 1 GB).

Uses extract_cd4.py's own `resolve_target_rows` (exact symbol, casefold, Ensembl id), so the counts equal its dry runs.
Then asks: on an ESSENTIAL-FREE Arc-style panel (A/B/C is essential-free; Arc's 2025 H1 panel carries 78 K562-essential
targets), what is the chance that a fresh 300-target draw trips the extraction gate (any condition < 275 or union < 290)?
The per-target coverage rates are taken from the non-essential targets of the two Arc panels; the draw is binomial with
the observed per-condition correlation handled by sampling targets (bootstrap over the pooled non-essential targets).

usage: uv run python cd4cover481.py -> work/cd4cover481.json
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SRC = ROOT / "experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad"
spec = importlib.util.spec_from_file_location("x047", ROOT / "experiments/047-cd4-direct/extract_cd4.py")
x047 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(x047)
PANELS = ["abc", "h1_2025", "efr", "efp", "real199", "pess199"]
ESS = {ln.strip() for ln in (ROOT / "data/donors/k562_essential/perts.txt").read_text().split("\n") if ln.strip()}


def main() -> None:
    with h5py.File(SRC, "r") as h:
        obs = read_elem(h["obs"])
        var = read_elem(h["var"])
    out = {}
    pooled = []   # (panel, target, essential, covered per condition)
    for name in PANELS:
        panel = [str(t) for t in pd.read_csv(HERE / "work/panels" / f"{name}.csv")["target_gene"]]
        sel, _ = x047.resolve_target_rows(obs, var, panel)
        cov = {c: set(sel[c]) for c in x047.CONDITIONS}
        union = set().union(*cov.values())
        ess = [t for t in panel if t in ESS]
        miss = [t for t in panel if t not in union]
        rec = {"per_condition": {c: len(v) for c, v in cov.items()}, "union": len(union),
               "complete": len(set.intersection(*cov.values())), "n_essential": len(ess),
               "union_missing": len(miss), "union_missing_essential": sum(t in ESS for t in miss),
               "nonessential_union_rate": (sum(t in union for t in panel if t not in ESS)
                                           / max(1, sum(t not in ESS for t in panel))),
               "nonessential_min_condition_rate": min(sum(t in v for t in panel if t not in ESS)
                                                      for v in cov.values()) / max(1, sum(t not in ESS for t in panel))}
        out[name] = rec
        print(f"{name:8s} {rec}", flush=True)
        if name in ("abc", "h1_2025"):
            for t in panel:
                if t not in ESS:
                    pooled.append((name, [t in cov[c] for c in x047.CONDITIONS]))
    # a fresh essential-free Arc-style draw: bootstrap 300 targets from the Arc panels' non-essential targets
    rng = np.random.default_rng(481)
    for tag, keep in (("abc_only", ("abc",)), ("h1_2025_nonessential_only", ("h1_2025",)), ("both", ("abc", "h1_2025"))):
        P = np.asarray([row for (pn, row) in pooled if pn in keep], dtype=bool)
        fire, unions, B = 0, [], 20000
        for _ in range(B):
            d = P[rng.integers(0, len(P), 300)]
            u = int(d.any(1).sum())
            unions.append(u)
            fire += (u < 290) or (d.sum(0).min() < 275)
        out[f"bootstrap_{tag}"] = {"pool_targets": int(len(P)), "draws": B, "p_gate_fires": fire / B,
                                   "union_q05_25_50_75_95": np.quantile(unions, [0.05, 0.25, 0.5, 0.75, 0.95]).tolist()}
        print(tag, out[f"bootstrap_{tag}"], flush=True)
    (HERE / "work/cd4cover481.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
