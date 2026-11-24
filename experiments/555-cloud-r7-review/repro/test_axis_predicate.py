"""The axis digest's collision, and what r7a's child actually compares.

A miniature of axes534.axis_sha: sha256 of the donor-gene union minus the panel, on a fixed
gene order. Two panels that differ only by genes outside that union share an axis. The child
in r7a551.axis_check does not compare T panels to each other, and ABC_AXES is a fixed pair of
digests (axes534.py:28-29), not a live recompute of data/challenge.
"""
from __future__ import annotations

import hashlib


def axis_sha(gene_order, panel, donor_genes):
    tindex = {g: i for i, g in enumerate(gene_order)}
    union = {g for g in donor_genes if g in tindex} - set(panel)
    wcols = sorted(tindex[g] for g in union)
    syms = [gene_order[i] for i in wcols]
    return hashlib.sha256("\n".join(syms).encode()).hexdigest(), len(wcols)


genes = [f"G{i}" for i in range(6)]
donors = ["G0", "G1", "G2", "G3"]
# G4, G5 are off the donor union. Panels that differ only there share an axis.
a, na = axis_sha(genes, {"G4"}, donors)
b, nb = axis_sha(genes, {"G5"}, donors)
c, nc = axis_sha(genes, {"G1"}, donors)
assert a == b and na == nb == 4
assert a != c
print("panels {{G4}} and {{G5}} share an axis", a[:12], "n", na)
print("panel {{G1}} differs", c[:12], "n", nc)
print("r7a child compares each operator to ABC_AXES (2 hardcoded digests), tuned_on_axes() "
      "(9 fixed bundles x 4 donor sets) and a recomputed R6 axis; it does not compare T1..T4 to each other")
print('{"shared_axis_outside_union": true}')
