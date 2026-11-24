"""238 Phase 1 — THE HONEST PANEL: 199's draw with the common-essential list removed.

199's two synthetic bundles draw their 300 from DONOR COVERAGE.  Donor libraries screen
common-essential genes, so both draws carry them: 34 (`real`) and 30 (`pess`).  **Arc's one
observed panel carries 0** (`DEF-RUNBOOK` §5, reproduced by `census238.py`), and there is an
assay reason to expect that again: an essential knockdown kills the cell, so there is no
transcriptional response to score.

So 233's and 236's October-source prices were measured on panels that contain a stratum Arc's
does not.  This builds the third panel that removes exactly that difference.

WHAT CHANGES, AND WHAT DOES NOT.  `defbundle199.py` / `defbundle199r.py` are reused verbatim
for everything except one line: the eligible pool has `data/donors/k562_essential/perts.txt`
subtracted before the draw.  Same seeds (199 / 1990), same gene axis, same control `.h5ad`
files (A/B/C relabelled D/E/F, byte-identical), same manifest surgery.  So the contrast
`pess -> efp` and `real -> efr` isolates ONE variable: whether the panel contains essential
genes.

  efp   pool = union of the three production donors' coverage,  minus abc, minus essential
  efr   pool = K562 n Orion-HCT116 n Orion-HEK293T n CD4 atlas, minus abc, minus essential

ASSERTED AT BUILD TIME (not checked later and hoped for):
  * exactly 300 unique targets, all on the official gene axis
  * |panel n abc| == 0                          -- disjoint, like both 199 draws
  * |panel n essential| == 0                    -- like Arc's real panel, unlike both 199 draws
  * the three control files are byte-identical to 199's, so nothing but the panel differs

    uv run python bundle238.py
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
SRC = ROOT / "data/challenge"
E199 = ROOT / "experiments/199-def-runbook"
CD4_SOURCE = E199 / "work/GWCD4i.DE_stats.h5ad"
ESSENTIAL_LIST = ROOT / "data/donors/k562_essential/perts.txt"
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
sys.path.insert(0, str(HERE))

# 199's module, imported for `relabel` and `sha256` so the control files are made the same way
spec = importlib.util.spec_from_file_location("defbundle199", E199 / "defbundle199.py")
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

SPECS = {
    # key      out dir          seed  pool          panel_id
    "efp": dict(out=WORK / "ef-bundle", seed=199, pool="union",
                panel_id="vcc2026-synthetic-ef-238"),
    "efr": dict(out=WORK / "ef-bundle-r", seed=1990, pool="intersect",
                panel_id="vcc2026-synthetic-ef-238r"),
}


def essential() -> set[str]:
    return {ln.strip() for ln in ESSENTIAL_LIST.read_text().splitlines() if ln.strip()}


def cd4_targets() -> set[str]:
    """199r's own reader, verbatim."""
    with h5py.File(CD4_SOURCE, "r") as h:
        node = h["obs"]["target_contrast_gene_name"]
        if isinstance(node, h5py.Group):
            return {c.decode() if isinstance(c, bytes) else str(c)
                    for c in node["categories"][:]}
        return {v.decode() if isinstance(v, bytes) else str(v) for v in node[:]}


def donor_set(name: str) -> set[str]:
    return {ln.strip() for ln in
            (ROOT / f"data/donors/{name}/perts.txt").read_text().splitlines() if ln.strip()}


def build(key: str) -> dict:
    spec_ = SPECS[key]
    out: Path = spec_["out"]
    t0 = time.time()
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    genes = [str(g) for g in pd.read_csv(SRC / "gene_names.csv")["gene_name"]]
    current = {str(t) for t in pd.read_csv(SRC / "pert_counts.csv")["target_gene"]}
    ess = essential()

    if spec_["pool"] == "union":                      # defbundle199.py's pool
        cover: set[str] = set()
        for d in DONORS:
            cover |= donor_set(d)
    else:                                             # defbundle199r.py's pool
        cover = None
        for d in DONORS:
            s = donor_set(d)
            cover = s if cover is None else (cover & s)
        cover &= cd4_targets()

    eligible_with = sorted(set(genes) & cover - current)          # 199's own pool
    eligible = sorted(set(eligible_with) - ess)                   # THE ONE CHANGE
    rng = np.random.default_rng(spec_["seed"])
    panel = sorted(rng.choice(np.asarray(eligible), size=300, replace=False).tolist())

    # --------------------------------------------------------------- assertions, at build time
    assert len(panel) == 300 and len(set(panel)) == 300
    assert not (set(panel) & current), "panel is not disjoint from A/B/C"
    assert not (set(panel) & ess), "panel still carries essential genes"
    assert set(panel) <= set(genes), "a target is not on the official gene axis"

    pd.DataFrame({"target_gene": panel}).to_csv(out / "pert_counts.csv", index=False)
    shutil.copyfile(SRC / "gene_names.csv", out / "gene_names.csv")

    manifest = json.loads((SRC / "manifest.json").read_text())
    manifest["partition"] = "test"
    manifest["panel_id"] = spec_["panel_id"]
    manifest["contexts"] = [base.RENAME[c] for c in manifest["contexts"]]
    manifest["per_context"] = {base.RENAME[k]: v for k, v in manifest["per_context"].items()}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    files = {}
    for old, new in base.RENAME.items():
        rec = base.relabel(SRC / f"context_{old}.h5ad", out / f"context_{new}.h5ad", new)
        files[f"context_{new}.h5ad"] = rec["sha256"]

    # the control files must be byte-identical to 199's, or something other than the panel moved
    ref = (E199 / "work/def-bundle") if key == "efp" else (E199 / "work/def-bundle-r")
    same = {}
    for f in ("context_D.h5ad", "context_E.h5ad", "context_F.h5ad", "gene_names.csv"):
        a, b = out / f, ref / f
        same[f] = (base.sha256(a) == base.sha256(b)) if b.exists() else None
    assert all(v is not False for v in same.values()), f"control files differ from 199's: {same}"

    ref_panel = ({t.strip() for t in (ref / "pert_counts.csv").read_text().splitlines()[1:]
                  if t.strip()} if (ref / "pert_counts.csv").exists() else set())
    rec = dict(
        key=key, panel_id=spec_["panel_id"], seed=spec_["seed"], pool=spec_["pool"],
        eligible_with_essential=len(eligible_with), eligible_essential_free=len(eligible),
        essential_removed_from_pool=len(eligible_with) - len(eligible),
        n_panel=300, panel_first5=panel[:5], panel_last5=panel[-5:],
        overlap_abc=0, overlap_essential=0,
        overlap_with_199_sibling=len(set(panel) & ref_panel),
        control_files_identical_to_199=same,
        files={**files, **{f: base.sha256(out / f) for f in
                           ("manifest.json", "gene_names.csv", "pert_counts.csv")}},
        essential_list_sha256=hashlib.sha256(ESSENTIAL_LIST.read_bytes()).hexdigest(),
        seconds=round(time.time() - t0, 1))
    print(f"[{key}] pool {spec_['pool']}: {len(eligible_with)} eligible -> "
          f"{len(eligible)} essential-free (removed {len(eligible_with) - len(eligible)});  "
          f"panel 300, n_abc 0, n_essential 0, overlap with 199's sibling "
          f"{rec['overlap_with_199_sibling']};  controls identical: "
          f"{sum(1 for v in same.values() if v)}/{len(same)}   ({rec['seconds']}s)")
    return rec


def main() -> None:
    out = {k: build(k) for k in SPECS}
    (WORK / "bundle238.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"-> {WORK / 'bundle238.json'}")


if __name__ == "__main__":
    main()
