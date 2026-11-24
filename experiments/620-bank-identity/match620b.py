"""620b — cross-chemistry bank identity with a leave-one-pair-out Flex-vs-3' gene bias.

Real control profiles only (banks whose b_ctrl is a placeholder are dropped: jiang_*, feng_ipsc, altair, dld1,
igvf_cardiomyocyte are all-zero; kolf2 is all-one). Flex probes: challenge A/B/C (Arc), h1_2025 (Arc), viperturb_k562
(Satija), cd4_rest (Marson). The chemistry bias per gene = mean over KNOWN same-line Flex-vs-non-Flex pairs of
(Flex - other): K562 (viperturb_k562 - k562_gwps), PSC (h1_2025 - pan_h1), Jurkat (challenge_A - nadig_jurkat).
Each Flex probe is corrected with a bias fitted WITHOUT its own pair (leave-one-pair-out), then ranked against the
non-Flex references by Spearman on reference-standardised genes. Negative controls: challenge B (HeLa), C (CAL33).
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np

import importlib.util
spec = importlib.util.spec_from_file_location("m620", Path(__file__).with_name("match620.py"))
m620 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m620)  # type: ignore[union-attr]

ROOT = m620.ROOT
DROP = {"jiang_K562", "feng_ipsc", "altair_a549", "jiang_A549", "jiang_HAP1", "jiang_BXPC3", "jiang_HT29", "jiang_MCF7",
        "dld1_lowmoi", "igvf_cardiomyocyte", "kolf2"}
FLEX = {"challenge_A", "challenge_B", "challenge_C", "h1_2025", "viperturb_k562", "cd4_rest"}
PAIRS = {"K562": ("viperturb_k562", "k562_gwps"), "PSC": ("h1_2025", "pan_h1"), "Jurkat": ("challenge_A", "nadig_jurkat")}


def main() -> None:
    prof: dict[str, dict[str, float]] = {}
    for n, (d, _, _) in m620.BANKS.items():
        if n in DROP:
            continue
        prof[n] = m620.load_bank(d)
    genes = [l.strip() for l in (ROOT / "data/challenge/gene_names.csv").read_text().splitlines()[1:] if l.strip()]
    for c, p in m620.CHALLENGE.items():
        prof[f"challenge_{c}"] = dict(zip(genes, map(float, np.load(p))))
    group = {n: m620.BANKS[n][1] for n in prof if n in m620.BANKS}
    group.update({f"challenge_{c}": m620.CHALLENGE_TRUTH[c] for c in "ABC"})
    names = list(prof)
    cnt = Counter(g for p in prof.values() for g in p)
    axis = sorted(g for g, k in cnt.items() if k >= 0.85 * len(names))
    M = np.array([[prof[n].get(g, np.nan) for g in axis] for n in names])
    idx = {n: i for i, n in enumerate(names)}
    refs = [n for n in names if n not in FLEX]
    out = {"axis": len(axis), "refs": refs, "probes": {}}

    def rank_corr(a: np.ndarray, b: np.ndarray) -> float:
        ok = np.isfinite(a) & np.isfinite(b)
        ra = np.argsort(np.argsort(a[ok])); rb = np.argsort(np.argsort(b[ok]))
        return float(np.corrcoef(ra, rb)[0, 1])

    for probe in sorted(FLEX):
        held = [k for k, (f, _) in PAIRS.items() if f == probe]
        use = [k for k in PAIRS if k not in held]
        bias = np.nanmean([M[idx[PAIRS[k][0]]] - M[idx[PAIRS[k][1]]] for k in use], axis=0)
        x = M[idx[probe]] - bias
        R = M[[idx[r] for r in refs]]
        mu = np.nanmean(R, axis=0); sd = np.nanstd(R, axis=0)
        keep = (mu > 1.0) & (sd > 0.1) & np.isfinite(bias)
        z_ref = (R[:, keep] - mu[keep]) / sd[keep]
        z_x = (x[keep] - mu[keep]) / sd[keep]
        scores = sorted(((rank_corr(z_x, z_ref[i]), r) for i, r in enumerate(refs)), reverse=True)
        same = [i for i, (_, r) in enumerate(scores) if group.get(r) == group.get(probe)]
        out["probes"][probe] = {"group": group.get(probe), "bias_from": use, "genes": int(keep.sum()),
                                 "top": [(r, group.get(r), round(s, 4)) for s, r in scores[:4]],
                                 "first_same_group_rank": (same[0] + 1) if same else None,
                                 "margin_top1_top2": round(scores[0][0] - scores[1][0], 4)}
        print(f"{probe:16s} [{group.get(probe)}] bias<-{use} genes {int(keep.sum())}: "
              + "; ".join(f"{r}({group.get(r)}) {s:+.3f}" for s, r in scores[:4])
              + f" | same-group rank {out['probes'][probe]['first_same_group_rank']}")
    (Path(__file__).parent / "work/match620b.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
