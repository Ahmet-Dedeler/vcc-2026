"""414: where does 410d's render-space energy come from? Per-gene / per-target energy of 410d vs 410c, context by context,
against the shipped multiplier m (read-only stream of both artifacts).
    uv run python diag410d.py [<ref.h5ad> <artifact.h5ad>]      (default: 410c vs 410d)
m is the 410d multiplier (a label only; meaningless for other arms). Reports max per-target energy ratio per context.
"""
import json, os, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
os.environ.setdefault("NUMBA_CACHE_DIR", str(HERE / "work" / "numba-cache"))
import numpy as np
sys.path.insert(0, str(HERE))
from shared414 import rows_of
ROOT = HERE.parents[1]
TS = 5e4
C = Path(sys.argv[1]) if len(sys.argv) > 2 else ROOT / "experiments/410-pertscale-revisit/work/build_c/the-kettle-applied-for-tenure.h5ad"
D = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "experiments/410-pertscale-revisit/work/build_d/the-lighthouse-keeper-went-remote.h5ad"
cache = {}
genes, Rc = rows_of(C, cache)
_, Rd = rows_of(D, cache)
m = np.asarray(np.load(ROOT / "experiments/410-pertscale-revisit/work/bulk_gene_gain_g228reg_t025.npz")["m"], float)
targets = sorted(Rc["A"])
gpos = {g: i for i, g in enumerate(genes)}
out = {}
for ctx in "ABC":
    b = cache["b"][ctx]
    def delta(R):
        T = np.vstack([R[ctx][t] for t in targets]); d = np.log1p(TS * T / T.sum(1, keepdims=True)) - b[None, :]
        for i, t in enumerate(targets): d[i, gpos[t]] = 0.0
        return d
    dc, dd = delta(Rc), delta(Rd)
    ge_c, ge_d = (dc ** 2).sum(0), (dd ** 2).sum(0)
    te_c, te_d = (dc ** 2).sum(1), (dd ** 2).sum(1)
    top = np.argsort(-(ge_d - ge_c))[:10]
    cpm = 1e6 * np.expm1(b) / np.expm1(b).sum()
    res = dict(energy_ratio=float(ge_d.sum() / ge_c.sum()),
               top_genes=[dict(gene=str(genes[g]), m=float(m[g]), cpm=float(cpm[g]), share_of_increase=float((ge_d[g] - ge_c[g]) / (ge_d.sum() - ge_c.sum())),
                               max_abs_delta_410d=float(np.abs(dd[:, g]).max()), max_abs_delta_410c=float(np.abs(dc[:, g]).max())) for g in top],
               top_targets=[dict(target=targets[i], e410c=float(te_c[i]), e410d=float(te_d[i])) for i in np.argsort(-(te_d - te_c))[:5]],
               share_increase_top10_genes=float((ge_d[top] - ge_c[top]).sum() / (ge_d.sum() - ge_c.sum())),
               median_target_ratio=float(np.median(te_d / np.maximum(te_c, 1e-12))),
               max_target_ratio=float(np.max(te_d / np.maximum(te_c, 1e-12))))
    out[ctx] = res
    print(f"\n{ctx}: energy ratio {res['energy_ratio']:.3f}; median per-target ratio {res['median_target_ratio']:.3f} (max "
          f"{res['max_target_ratio']:.2f}); top-10 genes carry "
          f"{100 * res['share_increase_top10_genes']:.1f}% of the increase")
    for r in res["top_genes"][:6]:
        print(f"   {r['gene']:12s} m {r['m']:+7.2f} cpm {r['cpm']:9.2f} share {100 * r['share_of_increase']:5.1f}%  max|d| 410d {r['max_abs_delta_410d']:.3f} (410c {r['max_abs_delta_410c']:.3f})")
    print("   top targets:", [(r["target"], round(r["e410c"], 2), round(r["e410d"], 2)) for r in res["top_targets"]])
(HERE / "work" / f"diag-{Path(C).stem}-vs-{Path(D).stem}.json").write_text(json.dumps(out, indent=1))
