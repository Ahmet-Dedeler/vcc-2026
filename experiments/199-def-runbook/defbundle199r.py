"""A SECOND synthetic D/E/F bundle with a REALISTIC panel draw.

`defbundle199.py` draws its 300 from the union of the three production donors' coverage,
which is dominated by the two ~18k-target Orion lines, so K562 covers only 150/300 and the
CD4 atlas 174/300 -- a deliberately pessimistic stress test of the coverage-drop path.

AGENTS.md §2.1 says the real panel is selected for "adequate representation in existing
public datasets", which on A/B/C gave K562 269/300, Orion 290-294 and CD4 291/300. This
bundle draws from the INTERSECTION of K562-gwps, both Orion lines and the CD4 atlas's
target set (read from the pinned source), again disjoint from the current 300, so its
census numbers are what Oct 22 should look like if §2.1 holds. It is used for the cheap
stages only (preflight / identity / census / CD4 dry-run); the heavy rehearsal runs on the
pessimistic bundle, where more can go wrong.

    uv run python experiments/199-def-runbook/defbundle199r.py
"""

from __future__ import annotations

import importlib.util
import json
import shutil
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SRC = ROOT / "data/challenge"
OUT = HERE / "work/def-bundle-r"
CD4_SOURCE = HERE / "work/GWCD4i.DE_stats.h5ad"
SEED = 1990

spec = importlib.util.spec_from_file_location("defbundle199", HERE / "defbundle199.py")
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)


def cd4_targets() -> set[str]:
    with h5py.File(CD4_SOURCE, "r") as h:
        node = h["obs"]["target_contrast_gene_name"]
        if isinstance(node, h5py.Group):
            cats = [c.decode() if isinstance(c, bytes) else str(c) for c in node["categories"][:]]
            return set(cats)
        return {v.decode() if isinstance(v, bytes) else str(v) for v in node[:]}


def main() -> None:
    t0 = time.time()
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    genes = [str(g) for g in pd.read_csv(SRC / "gene_names.csv")["gene_name"]]
    current = {str(t) for t in pd.read_csv(SRC / "pert_counts.csv")["target_gene"]}
    cover = None
    for donor in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        s = {ln.strip() for ln in (ROOT / f"data/donors/{donor}/perts.txt").read_text().splitlines() if ln.strip()}
        cover = s if cover is None else cover & s
    cover &= cd4_targets()
    eligible = sorted(set(genes) & cover - current)
    rng = np.random.default_rng(SEED)
    panel = sorted(rng.choice(np.asarray(eligible), size=300, replace=False).tolist())
    assert not (set(panel) & current)
    pd.DataFrame({"target_gene": panel}).to_csv(OUT / "pert_counts.csv", index=False)
    shutil.copyfile(SRC / "gene_names.csv", OUT / "gene_names.csv")
    manifest = json.loads((SRC / "manifest.json").read_text())
    manifest["partition"] = "test"
    manifest["panel_id"] = "vcc2026-synthetic-def-199r"
    manifest["contexts"] = [base.RENAME[c] for c in manifest["contexts"]]
    manifest["per_context"] = {base.RENAME[k]: v for k, v in manifest["per_context"].items()}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    files = {}
    for old, new in base.RENAME.items():
        rec = base.relabel(SRC / f"context_{old}.h5ad", OUT / f"context_{new}.h5ad", new)
        files[f"context_{new}.h5ad"] = rec["sha256"]
    out = {"seed": SEED, "panel_id": manifest["panel_id"], "eligible_pool": len(eligible),
           "panel_first5": panel[:5], "panel_last5": panel[-5:],
           "files": {**files, **{f: base.sha256(OUT / f) for f in ("manifest.json", "gene_names.csv", "pert_counts.csv")}},
           "seconds": round(time.time() - t0, 1)}
    (HERE / "work/def-bundle199r.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: out[k] for k in ("panel_id", "eligible_pool", "panel_first5", "seconds")}, indent=2))


if __name__ == "__main__":
    main()
