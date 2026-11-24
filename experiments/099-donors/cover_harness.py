"""Coverage of each harness's own scored perturbation panel, per candidate donor."""
import sys
from pathlib import Path
import numpy as np, anndata as ad, scipy.sparse as sp
ROOT = Path("~/Code/virtual-cell-challenge-2026")

CAND = {
    "k562_gwps":     ROOT/"data/donors/k562_gwps",
    "orion_hct116":  ROOT/"data/donors/orion_hct116",
    "orion_hek293t": ROOT/"data/donors/orion_hek293t",
    "h1_2025":       ROOT/"data/donors/h1_2025",
    "k562_essential":ROOT/"data/donors/k562_essential",
    "nadig_jurkat":  ROOT/"data/donors/nadig_jurkat",
    "rpe1":          ROOT/"data/donors/rpe1",
    "kolf2":         ROOT/"data/external/kolf2",
    "perturbmatch":  ROOT/"data/external/perturbmatch",
}
H = {"h1-test": ROOT/"data/harness/h1-test.h5ad",
     "orion":   ROOT/"data/harness/orion-hct116-thinned.h5ad"}

chal = [l.strip() for l in open(ROOT/"data/challenge/pert_counts.csv")][1:]
chal = set(t for t in chal if t)
axis = set(g.strip() for g in open(ROOT/"data/challenge/gene_names.csv") if g.strip())

panels = {}
for hn, hp in H.items():
    a = ad.read_h5ad(hp, backed="r")
    lab = set(a.obs["target_gene"].astype(str)) - {"non-targeting"}
    genes = set(str(g) for g in a.var_names)
    panels[hn] = (lab, genes)
    print(f"{hn}: {len(lab)} perturbations, {len(genes)} genes")
    del a

feng = ROOT/"data/external/feng-ipsc/genomewide_vcc_panel"
extra = {"feng_ipsc": (set(l.strip() for l in open(feng/"targets.txt") if l.strip()),
                       set(l.strip() for l in open(feng/"genes.txt") if l.strip()))}

print(f"\n{'donor':16s}{'perts':>8}{'genes':>8}{'chal300':>9}{'axis18k':>9}"
      + "".join(f"{h+'_p':>12s}{h+'_g':>10s}" for h in H))
rows = {}
for n, p in list(CAND.items()):
    perts = set(l.strip() for l in open(p/"perts.txt") if l.strip())
    genes = set(l.strip() for l in open(p/"genes.txt") if l.strip())
    rows[n] = (perts, genes)
for n, (perts, genes) in list(rows.items()) + list(extra.items()):
    s = f"{n:16s}{len(perts):>8}{len(genes):>8}{len(perts&chal):>9}{len(genes&axis):>9}"
    for hn in H:
        lab, hg = panels[hn]
        s += f"{len(perts&lab):>8}/{len(lab):<3}{len(genes&hg):>10}"
    print(s)
