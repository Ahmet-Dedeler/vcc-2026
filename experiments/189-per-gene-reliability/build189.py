"""189 Phase 1b — per-gene sign-agreement counts (k_g, n_g) for every cross-LINE pair.

One pass over the donor tables.  For each unordered pair of distinct cell LINES, on that pair's
own shared target set and on the working gene axis:

    centre across targets -> gate each line at its own per-gene median |effect| ->
    count supported pairs `n_g` and agreeing supported pairs `k_g`

exactly as `pilot189.py` reproduced.  Everything is computed TWICE: once on all shared targets
and once with the 300 challenge-panel targets removed, so the panel-disjointness falsifier is a
column rather than a re-run.

Genes absent from either table of a pair get n_g = 0 and are handled downstream by the
empirical-Bayes posterior (which returns the prior mean at n = 0).

Output: `work/pairs189.npz`.
"""
from __future__ import annotations

import itertools
import json
import os

import numpy as np

import argparse
import sys

from common189 import WORK, bvec, centre, gate_mask, load_table, pair_agreement, panel_targets

# one table per distinct cell LINE (161: k562_gwps and k562_essential are the same line, and
# gwps carries 4.8x the targets).  k562_essential is used only in the sensitivity check.
LINE_TABLE = {
    "HCT116": "orion_hct116",
    "HEK293T": "orion_hek293t",
    "K562": "k562_gwps",
    "RPE1": "rpe1",
    "JURKAT": "nadig_jurkat",
}
LINES = list(LINE_TABLE)
CHUNK = 1200          # genes per pass; see the note in the pair loop


def h1test_targets() -> set[str]:
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                                    "172-sign-confidence-key"))
    from harness172 import H1Test
    return set(map(str, H1Test().perts))


def main(drop: set[str] | None = None, out_name: str = "pairs189") -> None:
    os.makedirs(WORK, exist_ok=True)
    drop = drop or set()
    work = list(np.load(os.path.join(WORK, "work_genes.npy")))
    wi = {g: i for i, g in enumerate(work)}
    G = len(work)
    panel = panel_targets()
    tabs = {L: load_table(t) for L, t in LINE_TABLE.items()}

    out: dict[str, np.ndarray] = {}
    meta = []
    for A, B in itertools.combinations(LINES, 2):
        ta, tb = tabs[A], tabs[B]
        targets = sorted((set(ta["perts"]) & set(tb["perts"])) - drop)
        genes = [g for g in work if g in set(ta["genes"]) and g in set(tb["genes"])]
        gcol = np.array([wi[g] for g in genes])
        keep = np.array([t not in panel for t in targets])
        K = np.zeros(G)
        N = np.zeros(G)
        Kn = np.zeros(G)
        Nn = np.zeros(G)
        # Centring, the median gate and the agreement count are all COLUMN-separable, so a
        # gene-chunked pass is bit-identical to the whole-block one and keeps the anonymous
        # footprint at ~200 MB per pair instead of ~4 GB (RSS cap, four sibling agents).
        for c0 in range(0, len(genes), CHUNK):
            gs = genes[c0:c0 + CHUNK]
            col = gcol[c0:c0 + CHUNK]
            MA = centre(np.nan_to_num(_sub(ta, targets, gs)))
            MB = centre(np.nan_to_num(_sub(tb, targets, gs)))
            k, n = pair_agreement(MA, MB)
            K[col], N[col] = k, n
            kn, nn = pair_agreement(centre(MA[keep]), centre(MB[keep]))
            Kn[col], Nn[col] = kn, nn
            del MA, MB
        k, n, kn, nn = K[gcol], N[gcol], Kn[gcol], Nn[gcol]

        tag = f"{A}|{B}"
        out[f"k::{tag}"], out[f"n::{tag}"] = K, N
        out[f"kNP::{tag}"], out[f"nNP::{tag}"] = Kn, Nn
        rate = float(k.sum() / max(n.sum(), 1))
        meta.append(dict(pair=tag, A=A, B=B, targets=len(targets),
                         panel=int((~keep).sum()), genes=len(genes),
                         rate=rate, rate_np=float(kn.sum() / max(nn.sum(), 1)),
                         med_n=float(np.median(n)), med_n_np=float(np.median(nn))))
        print(f"{tag:<20} targets {len(targets):6d} (panel {int((~keep).sum()):4d})  "
              f"genes {len(genes):6d}  pooled agreement {rate:.4f}  "
              f"median n_g {np.median(n):7.0f}")

    # per-gene covariates on the working axis, per line (for the confounder controls)
    for L, t in LINE_TABLE.items():
        b = np.full(G, np.nan)
        gs = [g for g in work if g in set(tabs[L]["genes"])]
        b[[wi[g] for g in gs]] = bvec(tabs[L], gs)
        out[f"b::{L}"] = b

    out["genes"] = np.array(work)
    np.savez_compressed(os.path.join(WORK, f"{out_name}.npz"), **out)
    json.dump(meta, open(os.path.join(WORK, f"{out_name}.json"), "w"), indent=2)
    print(f"\nwrote work/{out_name}.npz ({len(meta)} cross-line pairs, {G} genes)")


def _sub(tab, targets, genes) -> np.ndarray:
    pi = {p: i for i, p in enumerate(tab["perts"])}
    gi = {g: i for i, g in enumerate(tab["genes"])}
    ridx = np.asarray([pi[t] for t in targets])
    cidx = np.asarray([gi[g] for g in genes])
    out = np.empty((len(targets), len(genes)), dtype=np.float32)
    for k in np.argsort(ridx):
        out[k] = np.asarray(tab["delta"][ridx[k]])[cidx]
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--drop-h1test-targets", action="store_true",
                    help="strict leak control: remove h1-test's 85 targets from every donor "
                         "pair before the counts are taken")
    ap.add_argument("--out", default="pairs189")
    a = ap.parse_args()
    d = h1test_targets() if a.drop_h1test_targets else set()
    if d:
        print(f"dropping {len(d)} h1-test targets from every donor pair")
    main(drop=d, out_name=a.out)
