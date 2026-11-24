"""414: the 410d gain's effect-space energy constant on PLOLO folds vs the real build's logged 0.532508, and where the
energy goes (share of bulk energy on the top-1% m genes, and on genes gated OUT of every challenge context).
    uv run python gainconst414.py kolf2
"""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from arms414 import (ALPHA, CD4_F, E062, FOLDS, GAMMA, LAM, PROD_OPS, UNION, UNION_G, Operator, build_fold, cd4_for,
                     frames, key_scale, merged, multiplier, shared_energy_key, solve_fraction, solve_fraction_carrying,
                     source_residual, union_rows, weights)
ROOT = HERE.parents[1]
out = {}
for held in sys.argv[1:]:
    fold = build_fold(held); wcols = np.asarray(fold["wcols"]); pre = np.asarray(fold["pre"], float)
    targets = list(fold["targets"])
    op = Operator(PROD_OPS / "whitener_union_pruned.npz" if FOLDS[held]["op_kind"] == "production" else Path(FOLDS[held]["op"]))
    cd4 = cd4_for(targets); h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    wh = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    f1, _ = solve_fraction(wh, cd4, CD4_F); f2, _ = solve_fraction_carrying(wh, h1, UNION_G)
    base, _ = E062.add_post_whiten_residual(wh, merged(cd4, f1, h1, f2), 1.0, wcols); base[:, fold["panel_positions"]] = 0.0
    key, carried, _ = shared_energy_key(targets, fold["sources"], op, wcols, cols=None)
    Ki = key_scale(base, weights(key, carried, 0.3, 1e9)[0])[0]
    z = np.load(ROOT / f"experiments/228-theta-gain-plolo/work/gain-{held}.npz", allow_pickle=True)
    m = multiplier(np.asarray(z["gain_reg"], float), np.asarray(z["live"], bool), 0.25)
    G = Ki * m[None, :]
    c = float(np.linalg.norm(Ki) / np.linalg.norm(G))
    gated_any = np.zeros(Ki.shape[1], bool)
    for ctx, (b, gate) in frames().items():
        gated_any |= gate
    e0, e1 = (Ki ** 2).sum(0), ((G * c) ** 2).sum(0)
    top = m >= np.quantile(m, 0.99)
    res = dict(energy_constant=c, share_top1pct_m_before=float(e0[top].sum() / e0.sum()),
               share_top1pct_m_after=float(e1[top].sum() / e1.sum()),
               share_ungated_before=float(e0[~gated_any].sum() / e0.sum()), share_ungated_after=float(e1[~gated_any].sum() / e1.sum()))
    out[held] = res
    print(held, json.dumps({k: round(v, 4) for k, v in res.items()}), "(real 410d build logged energy constant 0.532508)", flush=True)
(HERE / "work" / f"gainconst414-{'_'.join(sys.argv[1:])}.json").write_text(json.dumps(out, indent=1))
