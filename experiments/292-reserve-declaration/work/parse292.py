"""292 item 2b -- resolve each perturbation LABEL to a gene symbol, per dataset, from the label
vocabulary alone.

Input: 278's obs-only label counts (`experiments/278-depth-inventory/work/h5ad_counts/*.json`,
metadata written by an obs-only scanner) and 292's own `scan292.json` (obs/var/uns only).
No h5ad file is opened here, so no expression value can be touched.

The per-dataset label grammars below are FROZEN in this file and were written by reading the label
vocabularies only (the top labels printed in RESULTS.md).  A resolved symbol counts as a genetic
perturbation of a gene iff it is in the union of 279's 17,106 held map genes and the 300-target
challenge panel.

Output: work/perts292.json, work/inventory292.csv
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("~/Code/virtual-cell-challenge-2026")
WORK = ROOT / "experiments/292-reserve-declaration/work"
C278 = ROOT / "experiments/278-depth-inventory/work/h5ad_counts"

HELD = set(np.asarray(np.load(ROOT / "experiments/279-znf-state-dependence/work/held279.npz",
                              allow_pickle=True)["genes"]).astype(str).tolist())
PANEL = set(pd.read_csv(ROOT / "experiments/276-conservation-map/work/conservation276.csv")
            ["target"].astype(str).tolist())
UNIVERSE = HELD | PANEL

CTRL_PAT = re.compile(r"^(control|ctrl|non-?targeting.*|nontargeting.*|ntc|nt|none|nan|multiplet|"
                      r"p-INTERGENIC\d+|sgLACZ|sgOR2J2|GFP.*|Gal4.*|63\(mod\).*|safe-?targeting.*|"
                      r"unassigned|\*|)$", re.I)


def _ctrl(s: str) -> bool:
    return bool(CTRL_PAT.match(s.strip()))


# --------------------------------------------------------------- frozen per-dataset grammars
def g_identity(s):
    return s if re.fullmatch(r"[A-Za-z0-9\-\.]+", s) else None


def g_suffix_underscore(s):          # "<GENE>_pDS033", "<GENE>_2"
    return s.split("_")[0] if "_" in s else None


def g_tcrlibrary(s):                 # "Tcrlibrary_<GENE>_<n>"
    p = s.split("_")
    return p[1] if len(p) >= 3 and p[0] == "Tcrlibrary" else None


def g_dixit(s):                      # "p-sg<GENE>-<n>"
    m = re.fullmatch(r"p-sg(.+?)-\d+", s)
    return m.group(1) if m else None


def g_papalexi(s):                   # "<GENE>g<n>"
    m = re.fullmatch(r"(.+?)g\d+", s)
    return m.group(1) if m else None


def g_shifrut(s):                    # "ES.sg<n>.<GENE>"
    m = re.fullmatch(r"ES\.sg\d+\.(.+)", s)
    return m.group(1) if m else None


def g_tapseq(s):                     # "<GENE>_<+|->_<coord>.23-P1P2"; enhancers are "chr..:.."
    m = re.fullmatch(r"([A-Za-z0-9\-\.]+)_[+-]_[\d\.]+-P\S*", s)
    return m.group(1) if m else None


def g_single_only(s):                # plain symbols; drop double perturbations "A_B"
    return s if re.fullmatch(r"[A-Za-z0-9\-\.]+", s) else None


def g_mcfarland(s):                  # "sg<GENE>" / "sg<GENE>-<n>"
    m = re.fullmatch(r"sg([A-Za-z0-9]+?)(?:-\d+)?", s)
    return m.group(1) if m else None


def g_none(s):
    return None


GRAMMAR = {
    "AdamsonWeissman2016_GSM2406675_10X001": ("<GENE>_p<plasmid>", g_suffix_underscore),
    "AdamsonWeissman2016_GSM2406677_10X005": ("<GENE>_p<plasmid>", g_suffix_underscore),
    "AdamsonWeissman2016_GSM2406681_10X010": ("<GENE>_p<plasmid>", g_suffix_underscore),
    "DatlingerBock2017": ("Tcrlibrary_<GENE>_<guide>", g_tcrlibrary),
    "DatlingerBock2021": ("<GENE>_<guide>", g_suffix_underscore),
    "DixitRegev2016": ("p-sg<GENE>-<guide>", g_dixit),
    "PapalexiSatija2021_eccite_RNA": ("<GENE>g<guide>", g_papalexi),
    "PapalexiSatija2021_eccite_protein": ("<GENE>g<guide>", g_papalexi),
    "PapalexiSatija2021_eccite_arrayed_RNA": ("<GENE>", g_identity),
    "PapalexiSatija2021_eccite_arrayed_protein": ("<GENE>", g_identity),
    "ShifrutMarson2018": ("ES.sg<n>.<GENE>", g_shifrut),
    "SchraivogelSteinmetz2020_TAP_SCREEN__chromosome_11_screen": ("<GENE>_<strand>_<TSS>-P1P2", g_tapseq),
    "SchraivogelSteinmetz2020_TAP_SCREEN__chromosome_8_screen": ("<GENE>_<strand>_<TSS>-P1P2", g_tapseq),
    "TianKampmann2019_day7neuron": ("<GENE> (doubles '<A>_<B>' dropped)", g_single_only),
    "TianKampmann2019_iPSC": ("<GENE> (doubles '<A>_<B>' dropped)", g_single_only),
    "McFarlandTsherniak2020": ("sg<GENE>[-<n>] (drugs dropped)", g_mcfarland),
    "JoungZhang2023_atlas": ("numeric ORF id -- NOT RESOLVABLE from this file", g_none),
    "JoungZhang2023_combinatorial": ("<TF_A>+<TF_B> -- combinatorial only", g_none),
    "GasperiniShendure2019_atscale": ("multi-guide enhancer combination -- NOT RESOLVABLE", g_none),
    "GasperiniShendure2019_highMOI": ("multi-guide enhancer combination -- NOT RESOLVABLE", g_none),
    "GasperiniShendure2019_lowMOI": ("multi-guide enhancer combination + single TSS", g_identity),
    "XieHon2017": ("chr;start;end;strand enhancer coordinate -- NOT a gene", g_none),
}
DEFAULT = ("<GENE>", g_identity)


def resolve(stem: str, counts: dict[str, int]) -> dict:
    desc, fn = GRAMMAR.get(stem, DEFAULT)
    per_gene: dict[str, list[int]] = {}
    ctrl = 0
    unresolved = 0
    for lab, n in counts.items():
        s = str(lab).strip()
        if _ctrl(s):
            ctrl += n
            continue
        g = fn(s)
        if g is None or g not in UNIVERSE:
            unresolved += n
            continue
        per_gene.setdefault(g, []).append(n)
    cells = {g: int(sum(v)) for g, v in per_gene.items()}
    guides = {g: len(v) for g, v in per_gene.items()}
    arr = np.asarray(sorted(cells.values())) if cells else np.asarray([])
    out = {
        "grammar": desc,
        "n_gene_perts": len(cells),
        "n_in_held": len(set(cells) & HELD),
        "n_in_panel": len(set(cells) & PANEL),
        "control_cells": int(ctrl),
        "unresolved_cells": int(unresolved),
        "cells_on_gene_perts": int(arr.sum()) if len(arr) else 0,
        "median_cells": float(np.median(arr)) if len(arr) else None,
        "q10_cells": float(np.quantile(arr, 0.10)) if len(arr) else None,
        "q90_cells": float(np.quantile(arr, 0.90)) if len(arr) else None,
        "n_ge50": int((arr >= 50).sum()) if len(arr) else 0,
        "n_ge100": int((arr >= 100).sum()) if len(arr) else 0,
        "n_ge50_in_held": len([g for g, n in cells.items() if n >= 50 and g in HELD]),
        "n_ge50_in_panel": len([g for g, n in cells.items() if n >= 50 and g in PANEL]),
        "median_guides_per_gene": float(np.median(list(guides.values()))) if guides else None,
        "n_genes_multi_guide": int(sum(1 for v in guides.values() if v >= 2)),
        "genes": sorted(cells),
    }
    return out


def main() -> None:
    scan = {r["file"][:-5]: r for r in json.loads((WORK / "scan292.json").read_text())}
    out = {}
    for stem in sorted(scan):
        j = C278 / f"scperturb__{stem}.json"
        counts = json.loads(j.read_text()).get("counts", {})
        out[stem] = resolve(stem, counts)
        r = out[stem]
        print(f"{stem:<58} genes={r['n_gene_perts']:>5} held={r['n_in_held']:>5} "
              f"panel={r['n_in_panel']:>4} ge50={r['n_ge50']:>5} med={r['median_cells']} "
              f"ctrl={r['control_cells']:>7} guides/gene={r['median_guides_per_gene']}")
    (WORK / "perts292.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
