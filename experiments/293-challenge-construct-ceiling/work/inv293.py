"""293 -- INVENTORY (metadata only): every target with >= 2 constructs in each 2025 split.

Reads only obs (guide_id, target_gene, batch) from the three 2025 h5ad files. No expression value.

    PYTHONDONTWRITEBYTECODE=1 uv run python inv293.py
"""
from __future__ import annotations

import json
import re
import time
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
WORK = Path(__file__).resolve().parent

SPLITS = {
    "train": ROOT / "data/external/vcc2025/train/adata_Training.h5ad",
    "validation": ROOT / "data/external/vcc2025/validation/adata_Validation.h5ad",
    "test": ROOT / "data/external/vcc2025/test/adata_Test.h5ad",
}


def labels_of(h, path):
    o = h[path]
    if isinstance(o, h5py.Group) and "categories" in o:
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in o["categories"][:]])
        return cats[o["codes"][:]]
    v = o["values"][:] if isinstance(o, h5py.Group) else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def construct_of(gid: str) -> str:
    """A dual-guide construct id: 'GENE_P1_A|GENE_P1_B' -> 'GENE_P1'. Falls back to the raw id."""
    parts = gid.split("|")
    stems = []
    for p in parts:
        # strip a trailing _A / _B guide-position suffix
        m = re.match(r"^(.*)_([AB])$", p)
        stems.append(m.group(1) if m else p)
    if len(set(stems)) == 1:
        return stems[0]
    return "|".join(parts)  # unexpected shape: keep verbatim


def main() -> None:
    t0 = time.time()
    out = {}
    for split, path in SPLITS.items():
        with h5py.File(path, "r") as h:
            gid = labels_of(h, "obs/guide_id")
            tgt = labels_of(h, "obs/target_gene")
            bat = labels_of(h, "obs/batch")
            shape = [int(x) for x in h["X"].attrs["shape"]]
            genes = labels_of(h, "var/_index")
        n = len(gid)
        ctrl_mask = np.asarray([t in ("non-targeting", "NTC", "control") for t in tgt])
        # map target -> constructs
        cons = np.asarray([construct_of(g) for g in gid])
        by_t = defaultdict(set)
        for t, c in zip(tgt, cons):
            by_t[t].add(c)
        multi = {t: sorted(cs) for t, cs in by_t.items() if len(cs) > 1 and not (t in ("non-targeting",))}
        rows = []
        for t, cs in sorted(multi.items()):
            for c in cs:
                m = (tgt == t) & (cons == c)
                bb = Counter(bat[m].tolist())
                rows.append(dict(split=split, target=t, construct=c, n_cells=int(m.sum()),
                                 n_batches=len(bb), batches=dict(sorted(bb.items())),
                                 guide_ids=sorted(set(gid[m].tolist()))))
        gset = set(genes.tolist())
        out[split] = dict(
            path=str(path), shape=shape, n_cells=n,
            n_targets=int(len(set(tgt.tolist()) - {"non-targeting"})),
            n_guide_categories=int(len(set(gid.tolist()))),
            n_constructs=int(len(set(cons.tolist()))),
            n_control_cells=int(ctrl_mask.sum()),
            control_guide_ids=sorted(set(gid[ctrl_mask].tolist()))[:20],
            n_control_guide_ids=int(len(set(gid[ctrl_mask].tolist()))),
            control_batches=dict(sorted(Counter(bat[ctrl_mask].tolist()).items())),
            n_batches=int(len(set(bat.tolist()))),
            batches=dict(sorted(Counter(bat.tolist()).items())),
            multi_targets=sorted(multi),
            n_multi=len(multi),
            multi_rows=rows,
            targets_in_var=int(sum(1 for t in set(tgt.tolist()) - {"non-targeting"} if t in gset)),
        )
        print(f"{split}: {n} cells, {out[split]['n_targets']} targets, "
              f"{out[split]['n_guide_categories']} guide cats, {out[split]['n_constructs']} constructs, "
              f"{len(multi)} multi-construct targets: {sorted(multi)}", flush=True)
        for r in rows:
            print(f"    {r['target']:10s} {r['construct']:28s} n={r['n_cells']:5d} batches={r['n_batches']} "
                  f"{r['guide_ids']}", flush=True)
    out["seconds"] = round(time.time() - t0, 1)
    (WORK / "inventory293.json").write_text(json.dumps(out, indent=1))
    print("done", out["seconds"], "s")


if __name__ == "__main__":
    main()
