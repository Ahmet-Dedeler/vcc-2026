"""305 step 0 -- input manifest and the split/construct mapping, from METADATA ONLY.

Opens the three 2025 h5ads for `obs/{guide_id,target_gene,batch}`, `var/_index` and `X.attrs['shape']`.
**No element of `X/data` is read here.**  Emits `work/input_manifest.json` (paths, bytes, sha256, shapes,
gene-list digests) and `work/mapping305.json` (per split: constructs, cells per construct, batch/lane
composition, the NTC parity partition -- 293's `ntc_par[c] = i % 2` over the sorted NTC constructs).

    PYTHONDONTWRITEBYTECODE=1 uv run python manifest305.py
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib305 as X  # noqa: E402


def main() -> None:
    t0 = time.time()
    man: dict = {"generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "files": {}}
    mapping: dict = {}
    gene_digests = {}
    for s, path in X.SPLITS.items():
        gid, tgt, bat, genes, shape = X.read_obs(path)
        cons = np.asarray([X.construct_of(g) for g in gid])
        is_ctrl = tgt == X.CTRL
        ntc_cons = sorted(set(cons[is_ctrl].tolist()))
        ntc_par = {c: i % 2 for i, c in enumerate(ntc_cons)}          # 293 run293.py:92, verbatim
        allc = sorted(set(cons.tolist()))
        cells = {c: int((cons == c).sum()) for c in allc}
        tgt_of = {}
        for c, t in zip(cons, tgt):
            tgt_of[c] = t
        lanes = sorted(set(bat.tolist()))
        batch_of_c = {c: dict(sorted(
            {b: int(((cons == c) & (bat == b)).sum()) for b in lanes}.items()) ) for c in allc}
        batch_of_c = {c: {b: n for b, n in v.items() if n} for c, v in batch_of_c.items()}
        gd = hashlib.sha256("\n".join(genes.tolist()).encode()).hexdigest()
        gene_digests[s] = gd
        man["files"][str(path.relative_to(X.ROOT))] = dict(
            bytes=path.stat().st_size, sha256=X.sha256_file(path),
            n_cells=int(shape[0]), n_genes=int(shape[1]), n_lanes=len(lanes),
            n_constructs=len(allc), n_ntc_constructs=len(ntc_cons),
            n_targeting_constructs=len(allc) - len(ntc_cons),
            n_ctrl_cells=int(is_ctrl.sum()), gene_list_sha256=gd)
        mapping[s] = dict(
            n_cells=int(shape[0]), n_genes=int(shape[1]), lanes=lanes,
            constructs=allc, cells_per_construct=cells, target_of_construct=tgt_of,
            ntc_constructs=ntc_cons, ntc_parity=ntc_par,
            ntc_cells_by_parity={p: int(sum(cells[c] for c in ntc_cons if ntc_par[c] == p)) for p in (0, 1)},
            batch_of_construct=batch_of_c)
        print(f"{s:11s} {shape[0]:>8,} cells x {shape[1]} genes | {len(allc)} constructs "
              f"({len(allc) - len(ntc_cons)} targeting + {len(ntc_cons)} NTC) | {len(lanes)} lanes | "
              f"NTC cells {mapping[s]['ntc_cells_by_parity'][0]:,}/{mapping[s]['ntc_cells_by_parity'][1]:,} "
              f"by parity | genes sha {gd[:12]}")

    man["gene_list_identical_across_splits"] = len(set(gene_digests.values())) == 1
    man["derived_inputs"] = {
        str(p.relative_to(X.ROOT)): dict(bytes=p.stat().st_size, sha256=X.sha256_file(p))
        for p in (X.HERE / "candidates305.csv", X.HERE / "candidates305_meta.json",
                  X.ROOT / "experiments/298-reviewer-checks/work/seed298.csv",
                  X.ROOT / "experiments/301-metadata-acquisitions/work/seed_hits301_cage.csv",
                  X.ROOT / "experiments/301-metadata-acquisitions/work/vcc2025_spacers.csv",
                  X.ROOT / "experiments/301-metadata-acquisitions/work/promoters_cage_tss_hg19.fa.gz",
                  X.ROOT / "experiments/293-challenge-construct-ceiling/work/run293_train.json",
                  X.ROOT / "experiments/293-challenge-construct-ceiling/work/run293_validation.json",
                  X.ROOT / "experiments/293-challenge-construct-ceiling/work/run293_test.json",
                  X.ROOT / "experiments/293-challenge-construct-ceiling/work/lib293.py",
                  X.ROOT / "experiments/285-reliability-benchmark/work/lib285.py")}
    (X.HERE / "input_manifest.json").write_text(json.dumps(man, indent=1))
    (X.HERE / "mapping305.json").write_text(json.dumps(mapping, indent=1))
    print(f"\ngene list identical across the three splits: {man['gene_list_identical_across_splits']}")
    print(f"done {time.time() - t0:.0f}s, peak {X.peak_gib():.2f} GiB")


if __name__ == "__main__":
    main()
