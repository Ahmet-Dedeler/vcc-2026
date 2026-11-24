"""414: would sharpening blow up any REAL 410c block the way 228's gain blew up SYVN1-C in 410d?

The 410c effect matrix is not on disk, so its rendered bulk stands in: d = log1p(5e4 tot/sum) - b_ctrl per (context,
target) block of the 410c artifact, panel coordinates zeroed. Each row is re-rendered through production's
profile_from_effect (minnorm mass fix, the context's 10-CPM gate), with and without 169's row-renormalised sharpening.
Reports the render-space energy ratio per context (sharpened / q = 1 re-render), the max per-row energy ratio, and the
largest single-gene share of a block's mass. Known answer: the q = 1 re-render reproduces the artifact's own row
energies (the render is idempotent up to its mass fix); its ratio to the artifact is printed.

    uv run python sharp410c.py
"""
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
os.environ.setdefault("NUMBA_CACHE_DIR", str(HERE / "work" / "numba-cache"))
import numpy as np  # noqa: E402

sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(HERE.parents[1] / "experiments/169-sharpening-six-member"))
from shared414 import rows_of  # noqa: E402
from plolo import TS, MIN_CPM, profile_from_effect  # noqa: E402
from sharp import sharpen  # noqa: E402

ROOT = HERE.parents[1]
ART = ROOT / "experiments/410-pertscale-revisit/work/build_c/the-kettle-applied-for-tenure.h5ad"
QS = (1.0, 1.15, 1.25, 1.5)


def main():
    cache = {}
    genes, R = rows_of(ART, cache)
    targets = sorted(R["A"])
    gpos = {g: i for i, g in enumerate(genes)}
    panel = np.zeros(genes.size, bool)
    panel[[gpos[t] for t in targets]] = True
    out = {}
    for ctx in "ABC":
        b = cache["b"][ctx]
        expr = np.expm1(b)
        cpm = 1e6 * expr / expr.sum()
        gate = (cpm > MIN_CPM) & (b > 0)
        T = np.vstack([R[ctx][t] for t in targets])
        d = np.log1p(TS * T / T.sum(1, keepdims=True)) - b[None, :]
        e_art = np.einsum("ij,ij->i", np.where(panel, 0, d), np.where(panel, 0, d))
        d[:, panel] = 0.0
        res = {}
        E = {}
        for q in QS:
            F = sharpen(d, q, gate, mode="row")
            dl, share = np.empty_like(F), np.empty(F.shape[0])
            for i in range(F.shape[0]):
                p = profile_from_effect(F[i], b, gate, "minnorm")
                dl[i] = np.log1p(TS * p) - b
                share[i] = p.max()
            dl[:, panel] = 0.0
            E[q] = np.einsum("ij,ij->i", dl, dl)
            res[f"q{q:g}"] = dict(energy_vs_q1=float(E[q].sum() / E[1.0].sum()) if q != 1.0 else 1.0,
                                  max_row_ratio_vs_q1=float((E[q] / E[1.0]).max()) if q != 1.0 else 1.0,
                                  q99_row_ratio_vs_q1=float(np.quantile(E[q] / E[1.0], 0.99)) if q != 1.0 else 1.0,
                                  max_gene_mass_share=float(share.max()),
                                  target_of_max_share=targets[int(share.argmax())])
        res["ka_q1_rerender_over_artifact_energy"] = float(E[1.0].sum() / e_art.sum())
        out[ctx] = res
        print(f"{ctx}: KA q=1 re-render / artifact energy {res['ka_q1_rerender_over_artifact_energy']:.4f}", flush=True)
        for q in QS[1:]:
            r = res[f"q{q:g}"]
            print(f"   q {q:<4}: energy vs q1 {r['energy_vs_q1']:.4f}; row ratio max {r['max_row_ratio_vs_q1']:.3f} q99 "
                  f"{r['q99_row_ratio_vs_q1']:.3f}; max single-gene mass share {r['max_gene_mass_share']:.3f} "
                  f"({r['target_of_max_share']}) [q1: {res['q1']['max_gene_mass_share']:.3f}]", flush=True)
    (HERE / "work" / "sharp410c.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
