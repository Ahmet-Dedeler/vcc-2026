"""309 -- shared helpers.  Everything is read from 305's saved integer sufficient sums
(`sums305_<split>.npz`) and 301's / 307's saved sequence products.  No h5ad is opened.

The response formula is 305's, re-typed here exactly as 307's `lib307.py` has it, so that nothing
in 305's or 307's directory is invoked for writing:

    prof(S)[g] = log1p( 50000 * S[g] / sum_j S[j] )        # sum over ALL 18,080 genes
    e_c(g)     = prof(S_c)[g] - prof(S_ctrl)[g]
"""
from __future__ import annotations

import csv
import json
import resource
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
E301 = ROOT / "experiments/301-metadata-acquisitions/work"
E298 = ROOT / "experiments/298-reviewer-checks/work"
E305 = ROOT / "experiments/305-t7-expression-test/work"
E307 = ROOT / "experiments/307-offtarget-enrichment/work"
HERE = Path(__file__).resolve().parent
SPLITS = ("train", "validation", "test")
CTRL = "non-targeting"
SEED309 = 20260907
LS = (12, 11, 10, 9)


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def prof(S: np.ndarray) -> np.ndarray:
    S = np.asarray(S, dtype=np.float64)
    return np.log1p(50000.0 * S / S.sum())


class Split:
    """One split's construct-level integer sums, straight out of 305's npz."""

    def __init__(self, split: str):
        z = np.load(E305 / f"sums305_{split}.npz", allow_pickle=True)
        self.split = split
        self.genes = np.asarray([str(g) for g in z["genes"]])
        names = np.asarray([str(n) for n in z["names"]])
        halves = np.asarray(z["halves"])
        sums = np.asarray(z["sums"], dtype=np.int64)
        cells = np.asarray(z["cells"], dtype=np.int64)
        self.constructs = list(dict.fromkeys(names.tolist()))
        self.S, self.cells = {}, {}
        for c in self.constructs:
            m = names == c
            order = np.argsort(halves[m])
            self.S[c] = sums[m][order]
            self.cells[c] = cells[m][order]
        self.meta = json.load(open(E305 / f"sums305_{split}.json"))
        self.target_of = self.meta["target_of_construct"]
        self.ntc = sorted(self.meta["ntc_constructs"])
        self.ntc_parity = self.meta["ntc_parity"]
        self.targeting = [c for c in self.constructs if self.target_of[c] != CTRL]
        self.ctrlA = np.sum([self.S[c].sum(0) for c in self.ntc if self.ntc_parity[c] == 0], axis=0)
        self.ctrlB = np.sum([self.S[c].sum(0) for c in self.ntc if self.ntc_parity[c] == 1], axis=0)
        self.ctrl_pooled = self.ctrlA + self.ctrlB
        self.basal = prof(self.ctrl_pooled)
        self.raw_cp50k = 50000.0 * self.ctrl_pooled / self.ctrl_pooled.sum()
        self.gix = {g: i for i, g in enumerate(self.genes)}

    def effect(self, construct: str, loo_control: bool = False) -> np.ndarray:
        Sc = self.S[construct].sum(0)
        ref = self.ctrl_pooled - (Sc if loo_control else 0)
        return prof(Sc) - prof(ref)

    def half_effects(self, construct: str, loo_control: bool = False):
        S0, S1 = self.S[construct]
        a = self.ctrlA - (S0 if loo_control and self.ntc_parity.get(construct) == 0 else 0)
        b = self.ctrlB - (S1 if loo_control and self.ntc_parity.get(construct) == 1 else 0)
        if loo_control and construct in self.ntc_parity:
            whole = self.S[construct].sum(0)
            if self.ntc_parity[construct] == 0:
                a, b = self.ctrlA - whole, self.ctrlB
            else:
                a, b = self.ctrlA, self.ctrlB - whole
        return prof(S0) - prof(a), prof(S1) - prof(b)


def load_constructs() -> list[dict]:
    """307's frozen construct table (301's rule, L*, candidate sets) -- read only."""
    return list(csv.DictReader(open(E307 / "constructs307.csv")))


def load_genes() -> dict[str, dict]:
    out = {}
    for r in csv.DictReader(open(E307 / "genes307.csv")):
        out[r["gene"]] = r
    return out
