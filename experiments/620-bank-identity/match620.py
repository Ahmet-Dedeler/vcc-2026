"""620 — can a context's controls be matched to a BANK line's own control profile (not CCLE)?

Question: def_identity matches D/E/F controls against CCLE bulk RNA-seq only, and the 590/597 gate (registered AND
unanimous AND min_top1 >= 0.65) deliberately refuses lines absent from CCLE (HEK293, iPSC/ESC, RPE1-like, primary T).
Several of those cell types have large perturbation banks on disk (Orion HEK293T genome-wide, kolf2 10k, H1, RPE1,
CD4 T). This script asks whether the bank lines' own control pseudobulks (b_ctrl.npy, log1p CP50k) can identify a
probe's line or cell type across studies and chemistries, with known-answer probes:
  - cross-study same line: xucao_hek293 <-> orion_hek293t; viperturb_k562 (Flex) -> 3' K562 banks;
    h1_2025 (Arc Flex) -> gse295214 (Pan H1, 5'); challenge A (Arc Flex Jurkat) -> nadig_jurkat; A549 trio.
  - negative controls: challenge B (HeLa) and C (CAL33) have no bank line.
Analysis only; small memory (vectors only). Writes work/match620.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
T225 = ROOT / "experiments/225-source-sweep/work/tables"
BANKS = {  # name -> (dir, line group, chemistry)
    "orion_hek293t": (ROOT / "data/donors/orion_hek293t", "HEK293", "3p-FiCS"),
    "xucao_hek293": (T225 / "xucao_hek293", "HEK293", "3p"),
    "orion_hct116": (ROOT / "data/donors/orion_hct116", "HCT116", "3p-FiCS"),
    "k562_gwps": (ROOT / "data/donors/k562_gwps", "K562", "3p"),
    "viperturb_k562": (ROOT / "data/donors/viperturb_k562", "K562", "Flex"),
    "jiang_K562": (T225 / "jiang_K562", "K562", "3p"),
    "gasperini_k562": (T225 / "gasperini_k562", "K562", "3p"),
    "rpe1": (ROOT / "data/donors/rpe1", "RPE1", "3p"),
    "nadig_jurkat": (ROOT / "data/donors/nadig_jurkat", "Jurkat", "3p"),
    "nadig_hepg2": (ROOT / "experiments/166-donor-acquisition/donors/nadig_hepg2", "HepG2", "3p"),
    "frangieh_melanoma": (ROOT / "experiments/166-donor-acquisition/donors/frangieh_melanoma", "melanoma", "3p"),
    "h1_2025": (ROOT / "data/donors/h1_2025", "PSC", "Flex"),
    "pan_h1": (ROOT / "experiments/115-gse295214-h1/work/compact-donor", "PSC", "5p"),
    "kolf2": (ROOT / "data/external/kolf2", "PSC", "3p"),
    "feng_ipsc": (ROOT / "experiments/099-donors/work/feng_ipsc", "PSC", "3p"),
    "tian_ipsc_neuron": (T225 / "tian_ipsc_neuron", "neuron", "3p"),
    "igvf_cardiomyocyte": (T225 / "igvf_cardiomyocyte", "cardiomyocyte", "3p"),
    "telohaec": (ROOT / "data/external/igvf-telohaec", "endothelial", "3p"),
    "myllia_a549": (ROOT / "experiments/120-myllia-a549/work/compact-donor", "A549", "?"),
    "altair_a549": (ROOT / "experiments/223-silent-zero-audit/work/tables/altair_a549", "A549", "?"),
    "jiang_A549": (T225 / "jiang_A549", "A549", "3p"),
    "jiang_HAP1": (T225 / "jiang_HAP1", "HAP1", "3p"),
    "jiang_BXPC3": (T225 / "jiang_BXPC3", "BXPC3", "3p"),
    "jiang_HT29": (T225 / "jiang_HT29", "HT29", "3p"),
    "jiang_MCF7": (T225 / "jiang_MCF7", "MCF7", "3p"),
    "sunshine_calu3": (T225 / "sunshine_calu3", "Calu3", "3p"),
    "wessels_thp1": (T225 / "wessels_thp1", "THP1", "3p"),
    "dld1_lowmoi": (ROOT / "experiments/223-silent-zero-audit/work/tables/dld1_lowmoi", "DLD1", "3p"),
    "cd4_rest": (ROOT / "experiments/216-cd4-hard/work/cd4_compact_ka/cd4_rest", "CD4T", "3p"),
}
CHALLENGE = {c: ROOT / f"experiments/540-lodo-clip-scan/work/plolo/bctrl-challenge-{c}.npy" for c in "ABC"}
CHALLENGE_TRUTH = {"A": "Jurkat", "B": "HeLa(no bank)", "C": "CAL33(no bank)"}


def load_bank(d: Path) -> dict[str, float]:
    b = np.load(d / "b_ctrl.npy").astype(np.float64)
    genes = (d / "genes.txt").read_text().split("\n")
    genes = [g for g in genes if g != ""]
    assert len(genes) == b.shape[0], (d, len(genes), b.shape)
    out: dict[str, float] = {}
    for g, v in zip(genes, b):
        if g not in out:  # first occurrence wins
            out[g] = float(v)
    return out


def main() -> None:
    profiles: dict[str, dict[str, float]] = {}
    for n, (d, _, _) in BANKS.items():
        try:
            profiles[n] = load_bank(d)
        except Exception as e:  # noqa: BLE001
            print(f"skip {n}: {e}")
    chall_genes = [l.strip() for l in (ROOT / "data/challenge/gene_names.csv").read_text().splitlines()[1:] if l.strip()]
    for c, p in CHALLENGE.items():
        v = np.load(p)
        profiles[f"challenge_{c}"] = dict(zip(chall_genes, map(float, v)))
    names = list(profiles)
    # common axis: genes present in >= 85 % of profiles
    from collections import Counter
    cnt = Counter(g for p in profiles.values() for g in p)
    axis = sorted(g for g, k in cnt.items() if k >= 0.85 * len(names))
    M = np.full((len(names), len(axis)), np.nan)
    for i, n in enumerate(names):
        p = profiles[n]
        M[i] = [p.get(g, np.nan) for g in axis]
    ref_rows = [i for i, n in enumerate(names) if not n.startswith("challenge_")]
    # expressed in the reference set on average (log1p CP50k > 1 ~ 35 CP50k ~ 34 CPM... keep the informative range)
    mu = np.nanmean(M[ref_rows], axis=0)
    sd = np.nanstd(M[ref_rows], axis=0)
    keep = (mu > 1.0) & (sd > 0.1)
    Z = (M[:, keep] - mu[keep]) / sd[keep]
    print(f"profiles {len(names)}; axis {len(axis)} genes; scored {int(keep.sum())}")

    def corr(i: int, j: int) -> float:
        a, b = Z[i], Z[j]
        ok = np.isfinite(a) & np.isfinite(b)
        a, b = a[ok], b[ok]
        ra = np.argsort(np.argsort(a)); rb = np.argsort(np.argsort(b))
        return float(np.corrcoef(ra, rb)[0, 1])

    group = {n: BANKS[n][1] for n in BANKS if n in profiles}
    for c in "ABC":
        group[f"challenge_{c}"] = CHALLENGE_TRUTH[c]
    S = np.array([[corr(i, j) if i != j else np.nan for j in range(len(names))] for i in range(len(names))])
    out = {"axis_genes": len(axis), "scored_genes": int(keep.sum()), "probes": {}}
    for i, n in enumerate(names):
        order = [j for j in np.argsort(-np.nan_to_num(S[i], nan=-9)) if j != i and not names[j].startswith("challenge_")]
        top = [(names[j], group.get(names[j]), round(S[i, j], 4)) for j in order[:4]]
        same = [names[j] for j in order if group.get(names[j]) == group.get(n)]
        rank_same = (order.index(names.index(same[0])) + 1) if same else None
        out["probes"][n] = {"group": group.get(n), "top": top, "first_same_group_rank": rank_same}
        print(f"{n:20s} [{group.get(n)}] top: " + "; ".join(f"{a}({b}) {s:+.3f}" for a, b, s in top)
              + f"   first same-group rank: {rank_same}")
    (Path(__file__).parent / "work/match620.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
