"""413: panel-overlap helper.

exact(genes)    = |panel ∩ genes| on literal symbols (the number the mandate asks for).
resolved(genes) = same after mapping HGNC previous/alias symbols to the current approved symbol
                  (only unambiguous maps; reported separately so the literal number stays primary).
"""
import csv, pathlib, collections
HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PANEL = [r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))]
P = set(PANEL)

def _hgnc():
    prev, alias = collections.defaultdict(set), collections.defaultdict(set)
    approved = set()
    with open(HERE / "dl/hgnc_complete_set.txt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            s = r["symbol"]; approved.add(s)
            for p in filter(None, r["prev_symbol"].split("|")): prev[p].add(s)
            for a in filter(None, r["alias_symbol"].split("|")): alias[a].add(s)
    return approved, prev, alias
APPROVED, PREV, ALIAS = _hgnc()

def to_current(g):
    if g in APPROVED: return g
    for m in (PREV, ALIAS):
        if g in m and len(m[g]) == 1: return next(iter(m[g]))
    return g

def exact(genes):
    return len(P & set(genes))

PC = {to_current(p): p for p in P}

def resolved(genes):
    return len(set(PC) & {to_current(g) for g in genes})

def report(name, genes):
    genes = [g for g in genes if g]
    return {"source": name, "n_unique": len(set(genes)), "panel_exact": exact(genes), "panel_resolved": resolved(genes),
            "panel_hits": sorted(PC[c] for c in set(PC) & {to_current(g) for g in genes})}
