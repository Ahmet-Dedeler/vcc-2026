"""620e — within-lab centring instead of a cross-lab chemistry bias.

Each lab/chemistry group with >= 2 lines is centred on its own mean, which removes the shared lab + chemistry offset
(on Oct 22 the Arc group is A-F, six lines). Probes: Arc Flex {A, B, C} and Satija Flex {K562, HAP1, HEK293}.
References: Weissman 3' {k562_gwps, rpe1, nadig_jurkat, nadig_hepg2} and Orion FiCS {orion_hct116, orion_hek293t}.
Score = Spearman over genes of the centred profiles (renormalised over the shared axis as in 620d).
Known answers: A -> nadig_jurkat, satija_K562 -> k562_gwps, satija_HEK293 -> orion_hek293t (near-same-line).
Negatives: B (HeLa), C (CAL33), satija_HAP1.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
import importlib.util

here = Path(__file__).parent
spec = importlib.util.spec_from_file_location("m620", here / "match620.py")
m620 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m620)  # type: ignore[union-attr]
ROOT, TS = m620.ROOT, 50000.0
LABS = {
    "arc": ["challenge_A", "challenge_B", "challenge_C"],
    "satija": ["satija_K562", "satija_HAP1", "satija_HEK293"],
    "weissman": ["k562_gwps", "rpe1", "nadig_jurkat", "nadig_hepg2"],
    "orion": ["orion_hct116", "orion_hek293t"],
}
TRUTH = {"challenge_A": "nadig_jurkat", "satija_K562": "k562_gwps", "satija_HEK293": "orion_hek293t",
         "challenge_B": None, "challenge_C": None, "satija_HAP1": None}


def main() -> None:
    prof = {n: m620.load_bank(m620.BANKS[n][0]) for n in LABS["weissman"] + LABS["orion"]}
    genes = [l.strip() for l in (ROOT / "data/challenge/gene_names.csv").read_text().splitlines()[1:] if l.strip()]
    for c, p in m620.CHALLENGE.items():
        prof[f"challenge_{c}"] = dict(zip(genes, map(float, np.load(p))))
    z = np.load(ROOT / "experiments/597-flex-hek293-floor/work/pseudobulk597.npz", allow_pickle=True)
    for i, ln in enumerate([str(x) for x in z["lines"]]):
        cnt = z["counts"][i].astype(float)
        prof[f"satija_{ln}"] = dict(zip([str(x) for x in z["genes"]], np.log1p(TS * cnt / cnt.sum())))
    names = list(prof)
    cnt = Counter(g for p in prof.values() for g in p)
    axis = sorted(g for g, k in cnt.items() if k == len(names))  # genes present in every profile
    raw = np.array([[prof[n][g] for g in axis] for n in names])
    P = np.expm1(raw) / TS
    P = P / P.sum(axis=1, keepdims=True)
    M = np.log1p(TS * P)
    idx = {n: i for i, n in enumerate(names)}
    expressed = (M.mean(axis=0) > 1.0)
    C = {}
    for lab, members in LABS.items():
        mu = M[[idx[n] for n in members]].mean(axis=0)
        for n in members:
            C[n] = (M[idx[n]] - mu)[expressed]

    def rc(a, b):
        return float(np.corrcoef(np.argsort(np.argsort(a)), np.argsort(np.argsort(b)))[0, 1])

    refs = LABS["weissman"] + LABS["orion"]
    out = {"genes": int(expressed.sum()), "rows": []}
    for probe, truth in TRUTH.items():
        sc = sorted(((rc(C[probe], C[r]), r) for r in refs), reverse=True)
        top, r1 = sc[0]
        margin = top - sc[1][0]
        out["rows"].append({"probe": probe, "truth": truth, "top": r1, "score": round(top, 4), "margin": round(margin, 4),
                            "ranking": [(r, round(s, 4)) for s, r in sc]})
        print(f"{probe:14s} truth {str(truth):14s} | " + "; ".join(f"{r} {s:+.3f}" for s, r in sc))
    (here / "work/lab620e.json").write_text(json.dumps(out, indent=1))
    print("genes", int(expressed.sum()))


if __name__ == "__main__":
    main()
