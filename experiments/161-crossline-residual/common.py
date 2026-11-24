"""161 — shared loaders for the cross-line residual work.

Every donor table on disk is the SAME object in the SAME space (built by
`scripts/build_donor_table.py`):

    b       = log1p(TS * P / sum P)          TS = 5e4, so  sum(expm1(b)) == TS exactly
    delta_p = b_p - b_ctrl

so `b_ctrl + delta_p` is itself a valid `b`, i.e. **the reference effect lives exactly on the
mass manifold experiment 156 identified**.  That is what makes a donor table a legitimate
delivered-space evaluation surface: we can render a prediction against the recipient's own
`b_ctrl` on the recipient's own full gene axis and the truth satisfies the renderer's
constraint by construction.

Line GROUPS matter, not tables.  `k562_gwps` and `k562_essential` are two libraries of ONE
cell line and `cd4_rest`/`cd4_stim48hr`/`cd4_stim8hr` are one donor in three states; a
leave-one-line-out fold that splits either pair is leaking the answer (159 measured the
same-line pair at +0.2712 against +0.03..+0.16 for genuine cross-line pairs).
"""
from __future__ import annotations

import json
import os

import numpy as np

SEED = 20260902
TS = 5.0e4
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")

CD4 = "data/external/cd4t/cloud-compact/cd4_compact"

# table id -> (path, cell LINE group, human label)
TABLES = {
    "k562_gwps":     ("data/donors/k562_gwps",      "K562",    "K562 (Replogle gwps)"),
    "k562_essential": ("data/donors/k562_essential", "K562",    "K562 (Replogle essential)"),
    "rpe1":          ("data/donors/rpe1",           "RPE1",    "RPE1"),
    "nadig_jurkat":  ("data/donors/nadig_jurkat",   "JURKAT",  "Jurkat"),
    "orion_hct116":  ("data/donors/orion_hct116",   "HCT116",  "HCT116 (X-Atlas/Orion)"),
    "orion_hek293t": ("data/donors/orion_hek293t",  "HEK293T", "HEK293T (X-Atlas/Orion)"),
    "h1_2025":       ("data/donors/h1_2025",        "H1",      "H1 hESC (VCC2025)"),
    "kolf2":         ("data/external/kolf2",        "KOLF2",   "KOLF2.1J iPSC"),
    "cd4_rest":      (f"{CD4}/cd4_rest",            "CD4T",    "CD4 T, resting"),
    "cd4_stim8hr":   (f"{CD4}/cd4_stim8hr",         "CD4T",    "CD4 T, stim 8h"),
    "cd4_stim48hr":  (f"{CD4}/cd4_stim48hr",        "CD4T",    "CD4 T, stim 48h"),
}

# 159 §4 measured a large target-independent component ("blob") inside these two tables:
# within-line off-diagonal cosine +0.17..+0.28 against +0.002..+0.009 for the clean ones.
BLOBBY = {"rpe1", "kolf2"}


def panel_targets() -> set[str]:
    p = os.path.join(REPO, "data/challenge/pert_counts.csv")
    import csv
    with open(p) as fh:
        return {r["target_gene"] for r in csv.DictReader(fh)}


def load_table(tid: str, mmap: bool = True) -> dict | None:
    """Load a donor table.  `delta` is memory-mapped by default: orion is 2.7 GB each."""
    rel = TABLES[tid][0]
    p = os.path.join(REPO, rel)
    if not os.path.exists(os.path.join(p, "delta.npy")):
        return None
    perts = [l.strip() for l in open(os.path.join(p, "perts.txt")) if l.strip()]
    genes = [l.strip() for l in open(os.path.join(p, "genes.txt")) if l.strip()]
    delta = np.load(os.path.join(p, "delta.npy"), mmap_mode="r" if mmap else None)
    b = np.load(os.path.join(p, "b_ctrl.npy"))
    assert delta.shape == (len(perts), len(genes)), f"{tid}: shape mismatch"
    assert b.shape == (len(genes),), f"{tid}: b_ctrl shape mismatch"
    meta = {}
    mp = os.path.join(p, "meta.json")
    if os.path.exists(mp):
        meta = json.load(open(mp))
    return dict(id=tid, path=rel, line=TABLES[tid][1], label=TABLES[tid][2],
                perts=perts, genes=genes, delta=delta, b=b, meta=meta)


def submatrix(tab: dict, targets: list[str], genes: list[str]) -> np.ndarray:
    """Materialise a (targets x genes) float64 block from a memory-mapped delta."""
    pi = {p: i for i, p in enumerate(tab["perts"])}
    gi = {g: i for i, g in enumerate(tab["genes"])}
    ridx = np.asarray([pi[t] for t in targets])
    cidx = np.asarray([gi[g] for g in genes])
    order = np.argsort(ridx)                      # mmap row reads want ascending order
    out = np.empty((len(targets), len(genes)), dtype=np.float64)
    for k in order:
        out[k] = np.nan_to_num(np.asarray(tab["delta"][ridx[k]])[cidx], copy=False)
    return out


def bvec(tab: dict, genes: list[str]) -> np.ndarray:
    gi = {g: i for i, g in enumerate(tab["genes"])}
    return tab["b"][np.asarray([gi[g] for g in genes])].astype(np.float64)


def unit(A: np.ndarray) -> np.ndarray:
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)
