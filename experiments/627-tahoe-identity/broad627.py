#!/usr/bin/env python3
"""627c — the CRISPR-KO analogue of 627's Tahoe test on the Broad-16 (12 working lines; reserved 4 never read).

Each line is screened separately (its own OR-control guides), so there is no cross-line ambient soup. To avoid the
shared-control trap, the 10 OR control KOs of each line are split into two halves (5 + 5, fixed by name):
  response: lfc[t, L, g] = log2((pb_t + 1) / (ctrl_half_A + 1)) on CPM-normalised pseudobulks (authors' sums by KO);
  identity: I_L = log1p CPM(ctrl_half_B of L) - mean over the other 11 lines (their half B), genes > 5 CPM in L.
Line-specific deviation r[t, L] = lfc[t, L] - mean over lines of lfc[t, .] (removes the target's own shared response).
Stat per (t, L): Spearman(r, -I) over expressed genes; null: a random other line's identity. By response strength
(|lfc| norm) quintile. Analysis only, < 4 GB.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/622-identity-track"))
import broad622 as B  # noqa: E402

RNG = np.random.default_rng(6271)


def sp(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(rankdata(a[ok]), rankdata(b[ok]))[0, 1]) if ok.sum() > 100 else np.nan


def main() -> None:
    gp, PB, cp = B.load_matrix(B.P / "pseudobulk_sum_by_ko.csv", 0)
    gl, LFC, cl = B.load_matrix(B.P / "lfc_matrix.csv", 1)
    lines = sorted({ln for ln, _ in cp})
    assert len(lines) == 12 and not (set(lines) & B.RESERVED)
    all_t = sorted(set(k for _, k in cl))
    ors = sorted({k for _, k in cp if k not in set(all_t)})
    halfA, halfB = ors[0::2], ors[1::2]
    print(f"OR controls {len(ors)}: half A {halfA}, half B {halfB}")
    col = {(l2, k): j for j, (l2, k) in enumerate(cp)}

    def cpm(ln, kos):
        cols = [col[(ln, k)] for k in kos if (ln, k) in col]
        s = np.nansum(PB[:, cols], 1)
        return 1e6 * s / np.nansum(s)
    ctrlA = {ln: cpm(ln, halfA) for ln in lines}
    ctrlB = {ln: cpm(ln, halfB) for ln in lines}
    I = {}
    for ln in lines:
        lg = np.log1p(ctrlB[ln]); ref = np.mean([np.log1p(ctrlB[o]) for o in lines if o != ln], 0)
        I[ln] = np.where(ctrlB[ln] > 5, lg - ref, np.nan)
    targets = [t for t in all_t if all((ln, t) in col for ln in lines)]
    lfc = {}
    for ln in lines:
        lfc[ln] = np.stack([np.log2((cpm(ln, [t]) + 1) / (ctrlA[ln] + 1)) for t in targets])   # T x G
    mean_lfc = np.mean([lfc[ln] for ln in lines], 0)
    rows = []
    for ln in lines:
        r = lfc[ln] - mean_lfc
        for ti in range(len(targets)):
            other = RNG.choice([o for o in lines if o != ln])
            rows.append((float(np.linalg.norm(np.nan_to_num(lfc[ln][ti]))), sp(r[ti], -I[ln]), sp(r[ti], -I[other])))
    a = np.array(rows)
    out = {"n_pairs": len(a), "mean_rho": float(np.nanmean(a[:, 1])), "frac_pos": float(np.nanmean(a[:, 1] > 0)),
           "null": float(np.nanmean(a[:, 2])), "null_frac_pos": float(np.nanmean(a[:, 2] > 0)), "by_strength": []}
    print(f"pairs {len(a)} (12 lines x {len(targets)} targets): mean Spearman(r, -I) {out['mean_rho']:+.4f} "
          f"(positive {out['frac_pos']:.3f}); null {out['null']:+.4f} ({out['null_frac_pos']:.3f})")
    q = np.nanquantile(a[:, 0], [0.2, 0.4, 0.6, 0.8])
    for k in range(5):
        s = np.digitize(a[:, 0], q) == k
        row = dict(q=k + 1, n=int(s.sum()), rho=float(np.nanmean(a[s, 1])), null=float(np.nanmean(a[s, 2])))
        out["by_strength"].append(row)
        print(f"  strength Q{k + 1}: rho {row['rho']:+.4f} (null {row['null']:+.4f}), n {row['n']}")
    (HERE / "work/broad627.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
