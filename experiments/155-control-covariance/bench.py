"""155 — shared scoring frame.

Every arm in this experiment is scored on ONE fixed masking of 153's `Surface`, so the
incumbent and the candidate are never compared across different gene sets.  The masking is
the one 153 §8 used to produce the target numbers (sign acc 0.5366 / head purity 0.6646 /
reach 0.09040 / pds 0.7953): 104's frozen HCT116 axis intersected with the whitener's 17,942
canonical columns intersected with the harness's >5-CPM gate.

Metric definitions come from `../153-de-bench/debench.py` verbatim; nothing is redefined here.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "153-de-bench"))
from debench import Surface, overall_delta  # noqa: E402,F401

REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
SURFACE = os.path.join(REPO, "experiments/095-clean-harness/de_thinned.npz")
FROZEN104 = os.path.join(
    REPO, "experiments/104-context-heldout-credit-adapter/work/hct116/frozen.npz"
)
WCOLS = ("/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/"
         "aa534728-8641-48d7-a125-17e75d9c9399/scratchpad/wcols.npy")


def load_frame():
    """Return (surface, incumbent_on_surface_axis, gene_mask, pert_mask).

    `incumbent` is 104's `baseline` mapped onto the surface's (perts x genes) grid, zero
    outside the shared axis -- exactly what 153 q9 scores at effective gamma 0.65.
    """
    import csv

    s = Surface(SURFACE)
    canon = [r["gene_name"] for r in csv.DictReader(
        open(os.path.join(REPO, "data/challenge/gene_names.csv")))]
    wgenes = [canon[i] for i in np.load(WCOLS)]

    z = np.load(FROZEN104, allow_pickle=True)
    g104 = {str(g): i for i, g in enumerate(z["genes"])}
    t104 = {str(t): i for i, t in enumerate(z["targets"])}

    gpos = {g: i for i, g in enumerate(s.genes)}
    gsrc = np.array([g104[g] for g in wgenes if g in gpos])       # index into 104's axis
    gdst = np.array([gpos[g] for g in wgenes if g in gpos])       # index into surface axis
    prow = np.array([t104.get(p, -1) for p in s.perts])
    pok = prow >= 0

    keep = np.zeros_like(s.real_sig)
    keep[np.ix_(pok, gdst)] = True
    s.real_sig &= keep
    s.adj &= keep

    inc = np.zeros_like(s.pred_lfc)
    inc[np.ix_(pok, gdst)] = z["baseline"][np.ix_(prow[pok], gsrc)]
    return s, inc.astype(np.float64), gdst, pok


def head_purity(s, P, q: float = 90.0) -> tuple[float, int]:
    """Top-decile (by |predicted effect|) sign accuracy over the reference DE pool.

    Identical to 153 q8/q9: per perturbation, keep the reference-significant adjudicable
    genes, take those at or above the 90th percentile of |P|, pool across perturbations.
    Perturbations with fewer than 10 pool genes are skipped.
    """
    num = den = 0
    for p in range(len(s.perts)):
        i = np.flatnonzero(s.real_sig[p] & s.adj[p])
        if i.size < 10:
            continue
        m = i[np.abs(P[p, i]) >= np.percentile(np.abs(P[p, i]), q)]
        num += int((np.sign(P[p, m]) == np.sign(s.real_lfc[p, m])).sum())
        den += m.size
    return num / max(den, 1), den


def score(s, P, key=None):
    """The four endpoints this experiment reports, on one prediction matrix."""
    s.pred_lfc = P
    k = key if key is not None else [np.abs(P[p]) for p in range(len(s.perts))]
    hp, n = head_purity(s, P)
    return {
        "sign_acc": s.sign_accuracy(),
        "head_purity": hp,
        "head_n": n,
        "reach": s.reach(key=k),
        "pds": s.pds(),
        "nmae": s.nmae(),
    }


ROW = "{name:<34} {sign_acc:>9.4f} {head_purity:>12.4f} {reach:>9.5f} {pds:>8.4f} {nmae:>8.4f}"
HDR = ("{:<34} {:>9} {:>12} {:>9} {:>8} {:>8}"
       .format("arm", "sign acc", "head purity", "reach", "pds", "nmae"))
