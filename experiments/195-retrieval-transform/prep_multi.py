"""195 — build the multi-line training pool for one fold: every ordered pair of donor lines.

Stage B fits on the single line pair that survives a LODO fold, so its training signal has no
line boundary in it.  Here each fold gets four training lines -- the two non-held LODO lines
plus RPE1 and Nadig-Jurkat, neither of which is ever scored -- giving twelve ordered
(gallery line <- query line) directions and six line pairs.  That makes an inner
leave-one-LINE-out possible, which is the only way to ask, without touching the scored line,
whether a fitted correction survives a change of cell line at all.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from harness import LODO, Bank, build_fold, common_perts, fold_basis_rows, load_banks, sha256_arr  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OPS = HERE / "work" / "operators"
PREP = HERE / "work" / "prep_multi"
EXTRA = ["rpe1", "nadig_jurkat"]
BLOCK = 2048
MIN_PERTS = 512


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fold", required=True, choices=LODO)
    args = ap.parse_args()
    held = args.fold
    PREP.mkdir(parents=True, exist_ok=True)

    banks = load_banks()
    f = build_fold(held, banks, common_perts(), fold_basis_rows(banks)[held])
    del f["basis"]
    z = np.load(OPS / f"op.{held}.npz")
    V = z["evecs"]
    op_sha = sha256_arr(z["evals"], V)
    assert op_sha == json.loads((HERE / "work/phase0.json").read_text())["operators"][held]["sha256"]

    lines = f["others"] + EXTRA
    for n in EXTRA:
        banks[n] = Bank(n)
    vcc = set(pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"])
    ref = f["ref"]
    sc = {n: banks[n].gn_median_norm() for n in lines}

    meta = {"held": held, "lines": lines, "operator": {"path": str((OPS / f'op.{held}.npz').relative_to(ROOT)),
                                                       "sha256": op_sha},
            "ref": ref, "sc": sc, "directions": {}}
    for tgt, pred in itertools.permutations(lines, 2):
        shared = [p for p in banks[tgt].perts if p in banks[pred].idx and p not in vcc]
        if len(shared) < MIN_PERTS:
            continue
        rt = np.fromiter((banks[tgt].idx[p] for p in shared), np.int64, len(shared))
        rp = np.fromiter((banks[pred].idx[p] for p in shared), np.int64, len(shared))
        t0 = time.time()
        tag = f"{held}.{tgt}_from_{pred}"
        Yq = np.empty((len(shared), V.shape[0]), np.float32)
        Yt = np.empty((len(shared), V.shape[0]), np.float32)
        Q = (ref / sc[pred]) * banks[pred].gn(rp)
        for i in range(0, len(shared), BLOCK):
            Yq[i:i + BLOCK] = (np.asarray(Q[i:i + BLOCK], np.float64) @ V).astype(np.float32)
        del Q
        R = banks[tgt].raw
        for i in range(0, len(shared), BLOCK):
            Yt[i:i + BLOCK] = (np.asarray(R[rt[i:i + BLOCK]], np.float64) @ V).astype(np.float32)
        np.save(PREP / f"{tag}.Yq.npy", Yq)
        np.save(PREP / f"{tag}.Yt.npy", Yt)
        meta["directions"][f"{tgt}_from_{pred}"] = {
            "tag": tag, "n_perts": len(shared), "target_line": tgt, "query_line": pred,
            "perts_sha256": hashlib.sha256("\n".join(shared).encode()).hexdigest()}
        print(f"  {tag:52} {len(shared):6d} queries ({time.time()-t0:.1f}s)", flush=True)
        del Yq, Yt

    (PREP / f"{held}.meta.json").write_text(json.dumps(meta, indent=2))
    print("wrote", PREP / f"{held}.meta.json")


if __name__ == "__main__":
    main()
