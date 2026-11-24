"""279 item 1 — the HELD-set inventory: which NON-panel perturbed genes exist in >= 2 cell lines
on disk, how many of them are HGNC zinc-finger genes, and what the common gene axis would be.

METADATA ONLY: reads perts.txt / genes.txt and the HGNC gene-group table.  No delta.npy is opened,
so no response is looked at before the preregistration is hashed.

    PYTHONDONTWRITEBYTECODE=1 uv run python inv279.py
"""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/279-znf-state-dependence"
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/276-conservation-map"))
import lib276 as L  # noqa: E402

HGNC = ROOT / "data/external/xcell/hgnc_protein_coding.txt"

# Every single-line CRISPRi / CRISPR-KO bank in 276's scan (work/inventory276_sources.csv),
# INCLUDING the five 276 excluded only because they carry 0-1 PANEL targets.  Aggregates
# (jiang_6line, jiang_5line_noK562, the three CMap tables), mouse, CRISPRa and the two
# single-replicate DLD-1 duplicates of dld1_lowmoi are out.
CANDIDATES: dict[str, tuple[str, str, str]] = dict(L.SOURCES)
CANDIDATES.update({
    "k562_essential":  ("data/donors/k562_essential", "K562", "Replogle22-essential"),
    "rpe1":            ("data/donors/rpe1", "RPE1", "Replogle22-RPE1"),
    "nadig_jurkat":    ("data/donors/nadig_jurkat", "Jurkat", "Nadig25"),
    "nadig_hepg2":     ("experiments/166-donor-acquisition/donors/nadig_hepg2", "HepG2", "Nadig25"),
    "tian2019_ipsc":   ("experiments/225-source-sweep/work/tables/tian2019_ipsc", "iPSC-Tian19", "Tian19"),
    "gasperini_k562":  ("experiments/225-source-sweep/work/tables/gasperini_k562", "K562", "Gasperini19"),
    "wessels_thp1":    ("experiments/225-source-sweep/work/tables/wessels_thp1", "THP-1", "Wessels23"),
})

ZNF_GROUPS = ("Zinc fingers C2H2-type", "KRAB domain containing")


def hgnc_groups() -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    with open(HGNC) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = {x.strip() for x in (r["gene_group"] or "").split("|") if x.strip()}
            out[r["symbol"]] = g
    return out


def main() -> None:
    groups = hgnc_groups()
    panel = set(L.PANEL)

    per_src, lines_of = {}, {}
    for k, (rel, line, ds) in CANDIDATES.items():
        p = ROOT / rel
        if not (p / "perts.txt").exists():
            print(f"  MISSING {k} {rel}")
            continue
        perts, genes, _ = L.read_bank(rel)
        per_src[k] = dict(line=line, dataset=ds, rel=rel,
                          perts=set(perts), n_perts=len(perts),
                          genes=set(genes), n_genes=len(genes),
                          n_official=len(set(genes) & set(L.OFFICIAL)),
                          n_panel=len(set(perts) & panel),
                          n_nonpanel=len(set(perts) - panel))
        lines_of[k] = line

    # a gene's LINE count = number of distinct cell lines that perturbed it
    lines_per_gene: dict[str, set[str]] = defaultdict(set)
    srcs_per_gene: dict[str, set[str]] = defaultdict(set)
    for k, d in per_src.items():
        for g in d["perts"]:
            if g in panel:
                continue
            lines_per_gene[g].add(d["line"])
            srcs_per_gene[g].add(k)

    held2 = {g for g, ls in lines_per_gene.items() if len(ls) >= 2}
    held2_official = {g for g in held2 if g in set(L.OFFICIAL)}

    def cls(g: str) -> str:
        gg = groups.get(g, set())
        if any(x in gg for x in ZNF_GROUPS):
            return "ZNF"
        if any("inc finger" in x or x.startswith("Zinc") for x in gg):
            return "OTHER_ZF"
        return "NONZNF"

    counts = Counter(cls(g) for g in held2_official)
    print(f"non-panel genes in >= 2 lines: {len(held2)} ({len(held2_official)} on the challenge axis)")
    print(f"  classes: {dict(counts)}")
    for n in (2, 3, 4, 5, 6):
        sub = {g for g in held2_official if len(lines_per_gene[g]) >= n}
        c = Counter(cls(g) for g in sub)
        print(f"  >= {n} lines: {len(sub):6d}  ZNF {c['ZNF']:4d}  OTHER_ZF {c['OTHER_ZF']:4d}  nonZNF {c['NONZNF']:5d}")

    # per-source held-gene counts + ZNF counts
    print("\nper source (held = non-panel, on the challenge axis, in >= 2 lines):")
    rows = []
    for k, d in sorted(per_src.items(), key=lambda x: -len(x[1]["perts"] & held2_official)):
        h = d["perts"] & held2_official
        z = {g for g in h if cls(g) == "ZNF"}
        rows.append(dict(key=k, line=d["line"], dataset=d["dataset"], n_perts=d["n_perts"],
                         n_genes_official=d["n_official"], n_panel=d["n_panel"],
                         n_held=len(h), n_held_znf=len(z)))
        print(f"  {k:22s} {d['line']:16s} perts {d['n_perts']:6d}  axis {d['n_official']:6d}"
              f"  panel {d['n_panel']:3d}  held {len(h):6d}  heldZNF {len(z):4d}")

    # candidate common axes (metadata only)
    def axis_of(keys):
        s = set(L.OFFICIAL)
        for k in keys:
            s &= per_src[k]["genes"]
        return len(s)

    sets = {
        "276-26-sources": list(L.SOURCES),
        "all-candidates": [k for k in per_src],
        "held>=5-only": None,   # filled below
        "wide-6": ["k562_gwps", "orion_hct116", "orion_hek293t", "dld1_lowmoi",
                   "igvf_telohaec", "igvf_cardiomyocyte"],
        "wide-8": ["k562_gwps", "orion_hct116", "orion_hek293t", "dld1_lowmoi",
                   "igvf_telohaec", "igvf_cardiomyocyte", "gse295214_h1", "altair_a549"],
        "wide-8+replogle": ["k562_gwps", "orion_hct116", "orion_hek293t", "dld1_lowmoi",
                            "igvf_telohaec", "igvf_cardiomyocyte", "gse295214_h1", "altair_a549",
                            "rpe1", "nadig_jurkat", "nadig_hepg2", "k562_essential"],
    }
    keep5 = [k for k, d in per_src.items() if len(d["perts"] & held2_official) >= 5]
    sets["held>=5-only"] = keep5
    print("\ncommon-axis size for candidate source sets:")
    axes = {}
    for name, keys in sets.items():
        keys = [k for k in keys if k in per_src]
        axes[name] = dict(n_sources=len(keys), n_axis=axis_of(keys), keys=keys)
        print(f"  {name:18s} {len(keys):3d} sources -> {axis_of(keys):6d} genes")

    # how many held genes / ZNFs are covered in >=2 lines *within* each candidate set
    print("\nheld coverage inside each candidate set:")
    for name, d in axes.items():
        lp: dict[str, set[str]] = defaultdict(set)
        for k in d["keys"]:
            for g in per_src[k]["perts"]:
                if g not in panel:
                    lp[g].add(per_src[k]["line"])
        h = {g for g, ls in lp.items() if len(ls) >= 2 and g in set(L.OFFICIAL)}
        c = Counter(cls(g) for g in h)
        d.update(n_held=len(h), n_held_znf=c["ZNF"], n_held_nonznf=c["NONZNF"],
                 n_held_otherzf=c["OTHER_ZF"])
        print(f"  {name:18s} held {len(h):6d}  ZNF {c['ZNF']:4d}  OTHER_ZF {c['OTHER_ZF']:4d}"
              f"  nonZNF {c['NONZNF']:6d}")

    # the other large TF families available as controls
    fam_counts = Counter()
    for g in held2_official:
        for x in groups.get(g, set()):
            fam_counts[x] += 1
    print("\nlargest HGNC gene groups inside the held set (control-family candidates):")
    for f, c in fam_counts.most_common(30):
        print(f"  {c:5d}  {f}")

    WORK.mkdir(exist_ok=True)
    (WORK / "inventory279.json").write_text(json.dumps(dict(
        n_held_2lines=len(held2), n_held_2lines_official=len(held2_official),
        classes=dict(counts), per_source=rows, axes=axes,
        top_families=fam_counts.most_common(60)), indent=1))
    with open(WORK / "inventory279_sources.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
