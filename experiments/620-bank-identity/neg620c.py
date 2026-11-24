"""620c — negative controls by removing the truth: for each Flex probe whose line HAS a bank, drop every reference of
its own line group (and its own pair from the bias fit) and record the best wrong-line score. This is the false-positive
regime a D/E/F line with no bank would face. Also reports the 3'-vs-3' within-chemistry cross-line maxima (lab
confounding) for context.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
import importlib.util

here = Path(__file__).parent
spec = importlib.util.spec_from_file_location("m620b", here / "match620b.py")
m620b = importlib.util.module_from_spec(spec); spec.loader.exec_module(m620b)  # type: ignore[union-attr]
m620 = m620b.m620


def main() -> None:
    prof = {n: m620.load_bank(d) for n, (d, _, _) in m620.BANKS.items() if n not in m620b.DROP}
    genes = [l.strip() for l in (m620.ROOT / "data/challenge/gene_names.csv").read_text().splitlines()[1:] if l.strip()]
    for c, p in m620.CHALLENGE.items():
        prof[f"challenge_{c}"] = dict(zip(genes, map(float, np.load(p))))
    group = {n: m620.BANKS[n][1] for n in prof if n in m620.BANKS}
    group.update({f"challenge_{c}": m620.CHALLENGE_TRUTH[c] for c in "ABC"})
    names = list(prof)
    cnt = Counter(g for p in prof.values() for g in p)
    axis = sorted(g for g, k in cnt.items() if k >= 0.85 * len(names))
    M = np.array([[prof[n].get(g, np.nan) for g in axis] for n in names])
    idx = {n: i for i, n in enumerate(names)}

    def rc(a, b):
        ok = np.isfinite(a) & np.isfinite(b)
        return float(np.corrcoef(np.argsort(np.argsort(a[ok])), np.argsort(np.argsort(b[ok])))[0, 1])

    out = {}
    for probe in ["challenge_A", "h1_2025", "viperturb_k562", "challenge_B", "challenge_C", "cd4_rest"]:
        g0 = group[probe]
        use = [k for k, (f, o) in m620b.PAIRS.items() if f != probe and group.get(o) != g0]
        bias = np.nanmean([M[idx[m620b.PAIRS[k][0]]] - M[idx[m620b.PAIRS[k][1]]] for k in use], axis=0)
        refs = [n for n in names if n not in m620b.FLEX and group.get(n) != g0]
        R = M[[idx[r] for r in refs]]
        mu = np.nanmean(R, axis=0); sd = np.nanstd(R, axis=0)
        keep = (mu > 1.0) & (sd > 0.1) & np.isfinite(bias)
        zx = (M[idx[probe]] - bias)[keep]
        zx = (zx - mu[keep]) / sd[keep]
        zr = (R[:, keep] - mu[keep]) / sd[keep]
        sc = sorted(((rc(zx, zr[i]), r) for i, r in enumerate(refs)), reverse=True)
        out[probe] = {"truth_removed": g0, "bias_pairs": use, "top": [(r, group[r], round(s, 4)) for s, r in sc[:3]]}
        print(f"{probe:15s} truth group {g0!r} removed; bias<-{use}: " + "; ".join(f"{r}({group[r]}) {s:+.3f}" for s, r in sc[:3]))
    (here / "work/neg620c.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
