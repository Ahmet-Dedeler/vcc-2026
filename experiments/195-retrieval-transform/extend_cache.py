"""195 — add two never-scored donor lines to the cache, on the same 7,760-gene feat axis.

Stage B's negative has one structural escape route.  After a LODO fold holds out its scored
line, exactly **one** pair of lines remains, so the training set contains no line boundary at
all: nothing in the fit ever asks the correction to survive a change of cell line.  A
correction that is pair-specific would look perfect on donor dev (a perturbation holdout) and
still die on the panel, which is exactly the observed pattern.

RPE1 and Nadig-Jurkat are never a LODO fold, so adding them gives four training lines per fold
and six line pairs, and an inner leave-one-LINE-out becomes possible without ever touching the
scored line.  K562-essential is deliberately NOT added: it is the same cell line as
`k562_gwps`, which IS a scored fold, so it would leak that fold's biology.

Genes on the feat axis that a line does not measure are stored as exact zeros.  A zero query
coordinate contributes nothing to the inner product, so the retrieval problem is effectively
posed on the pair's shared sub-axis, which is the honest thing to do with a ragged gene panel.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
HERE = Path(__file__).resolve().parent
CACHE = HERE / "work" / "cache"
EXTRA = ["rpe1", "nadig_jurkat"]


def main() -> None:
    feat = (CACHE / "feat.txt").read_text().split("\n")
    fpos = {g: i for i, g in enumerate(feat)}
    man = json.loads((CACHE / "manifest.json").read_text())
    for n in EXTRA:
        d = ROOT / "data" / "donors" / n
        genes = (d / "genes.txt").read_text().split("\n")
        delta = np.load(d / "delta.npy")
        m = delta.mean(axis=0)
        s = delta - m
        sc = np.abs(s).mean(axis=0)
        live = sc[sc > 0]
        full_scale = np.maximum(sc, np.percentile(live, 25.0))

        have = [(fpos[g], i) for i, g in enumerate(genes) if g in fpos]
        col_f = np.fromiter((a for a, _ in have), np.int64, len(have))
        col_d = np.fromiter((b for _, b in have), np.int64, len(have))
        s_raw = np.zeros((s.shape[0], len(feat)), np.float32)
        s_raw[:, col_f] = s[:, col_d]
        scale = np.ones(len(feat))
        scale[col_f] = full_scale[col_d]

        np.save(CACHE / f"{n}.s_raw.npy", s_raw)
        np.save(CACHE / f"{n}.scale.npy", scale)
        man["donors"][n] = {
            "delta_path": str((d / "delta.npy").relative_to(ROOT)),
            "delta_sha256": hashlib.sha256((d / "delta.npy").read_bytes()).hexdigest(),
            "n_perts": int(s.shape[0]), "n_genes_full": int(s.shape[1]),
            "feat_genes_present": len(have),
            "s_raw_sha256": hashlib.sha256(s_raw.tobytes()).hexdigest(),
            "scale_sha256": hashlib.sha256(scale.tobytes()).hexdigest(),
            "note": "never a LODO fold; missing feat genes stored as exact zeros",
        }
        print(f"{n:16} {s_raw.shape}  feat genes present {len(have)}/{len(feat)}  "
              f"sha {man['donors'][n]['s_raw_sha256'][:16]}")
        del delta, s, s_raw
    (CACHE / "manifest.json").write_text(json.dumps(man, indent=2))


if __name__ == "__main__":
    main()
