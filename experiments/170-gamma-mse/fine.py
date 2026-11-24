"""170 — the SHAPE of the gamma response on h1-test, at 0.025 resolution, without emission.

`A`, `X` and `B` are functionals of the predicted PSEUDOBULK, and the quasi emitter's
exact-total repair reproduces that pseudobulk to rel L2 ~0.0014 (a21, board-verified at
`Delta pds = 0.000000`). So the whole `mse` numerator except `claim` -- and `pds`, and the
aggregate cosine -- can be evaluated straight off the rendered profile, with no cells and no
seed. `claim` is held at the value the emitted sweep measures, which is the same at every
gamma to five decimals (the emission policy is unchanged and `min(C_hat, C_p)` saturates).

That buys 23 gamma points instead of 9, which is what decides "monotone or interior optimum".
The emitted sweep is the ground truth; this is the shape, and `h1sweep.json`'s points are
overlaid on it as a check.

    VCC_MASS_FIX=minnorm uv run python fine.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
from h1kit import MIN_CPM, TS, Surface170, deliver, to_h1_axis  # noqa: E402

WHITENER = REPO / "data/donors/whitener_union_pruned.npz"
GRID = [round(0.300 + 0.025 * i, 3) for i in range(23)]     # 0.300 .. 0.850


def main() -> None:
    from cell_eval2.metrics.delta import _across_pert_budget, _row_weights
    from cell_eval2.metrics.discrimination import discrimination_score

    path = REPO / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("prod_builder_170f", path)
    B = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_170f"] = B
    spec.loader.exec_module(B)
    tag = "minnorm" if B.MASS_FIX else "none"
    print(f"[fine] MASS_FIX={tag}", flush=True)

    S = Surface170()
    pre = np.load(WORK / "pre.npz", allow_pickle=True)
    cd4z = np.load(WORK / "cd4.npz", allow_pickle=True)
    cols = pre["operator_columns"].astype(np.int64)
    h1_pos = pre["h1_positions"].astype(np.int64)
    cand = pre["cand_h1"]
    cd4 = B.EXP062.CD4Prepared(mode=str(cd4z["mode"][0]), mean_effect=cd4z["h1_mean"],
                               gene_votes=cd4z["h1_votes"],
                               target_condition_count=cd4z["h1_counts"],
                               condition_details=[])
    cd4c = B.EXP062.CD4Prepared(mode=str(cd4z["mode"][0]), mean_effect=cd4z["chal_mean"],
                                gene_votes=cd4z["chal_votes"],
                                target_condition_count=cd4z["chal_counts"],
                                condition_details=[])
    cand_c = pre["cand_chal"]
    panel_pos = pre["panel_positions"].astype(np.int64)
    whit = B.Whitener.load(WHITENER)

    gj = json.loads((WORK / "gate.json").read_text())
    fg = np.asarray(sorted(float(k) for k in gj))
    ff = np.asarray([gj[f"{v:.2f}"]["chal"]["cd4_fraction"] for v in fg])

    order = {str(t): i for i, t in enumerate(pre["h1_targets"])}
    row_order = np.asarray([order[str(t)] for t in S.targets])

    # the reference side, on h1-test's own axis
    rpos = {str(p): i for i, p in enumerate(S.rp)}
    d_true = S.rm[[rpos[str(t)] for t in S.targets]] - S.b_ctrl[None, :]
    gpos = {str(g): i for i, g in enumerate(S.genes)}
    cols_own = {i: gpos[str(t)] for i, t in enumerate(S.targets) if str(t) in gpos}
    w = _row_weights(S.targets.size, S.genes.size, cols_own)
    r = np.einsum("ij,ij->i", d_true, d_true)
    for k, c in cols_own.items():
        r[k] -= d_true[k, c] ** 2
    R = float((r * w).sum())

    out = []
    for gamma in GRID:
        t0 = time.time()
        wh, waudit = B.apply_v9_whitening(cand, whit, cols, gamma=gamma, lam=0.0)
        # The CD4 fraction production solves on the PANEL, interpolated from `gate.json`'s
        # ten exactly-measured gammas. It varies by 0.15% over the whole range and reaches
        # only 4 of h1-test's 85 rows, so interpolating it costs nothing and saves a
        # 300 x 17,942 @ 17,942 x 17,942 GEMM per grid point.
        frac = float(np.interp(gamma, fg, ff))
        final, _ = B.EXP062.add_post_whiten_residual(wh, cd4, frac, cols)
        eff, _ = to_h1_axis(final, pre["genes"], pre["h1_genes"])
        eff = eff[row_order]
        P, dh = deliver(B, eff, S.b_ctrl)
        pm = S.b_ctrl[None, :] + dh

        a = np.einsum("ij,ij->i", dh, dh)
        x = np.einsum("ij,ij->i", dh, d_true)
        for k, c in cols_own.items():
            a[k] -= dh[k, c] ** 2
            x[k] -= dh[k, c] * d_true[k, c]
        A = float((a * w).sum())
        X = float((x * w).sum())
        bud = float(_across_pert_budget(pm, cols_own, w))
        pds = float(np.mean(list(discrimination_score(
            pred_bulk=(S.targets, pm), real_bulk=(S.rp, S.rm), pert_col="target_gene",
            control="non-targeting", distance="cosine", rank_denominator="n-1",
            tie_policy="midrank", exclude_target_gene=True, exclusion_scope="panel",
            control_source="real", genes=S.genes, target_gene_map=None).values())))
        out.append(dict(gamma=gamma, A=A, X=X, B=bud, R=R, pds=pds,
                        cos=float(X / np.sqrt(A * R)), AmB=float(A - bud),
                        row_norm=float(np.median(np.linalg.norm(dh, axis=1))),
                        after_op=float(waudit["after_operator_median_norm"]),
                        secs=round(time.time() - t0, 1)))
        print(f"[fine] gamma {gamma:.3f}  A {A:.6f}  X {X:.6f}  B {bud:.6f}  "
              f"A-B {A - bud:.6f}  pds {pds:.5f}  cos {X / np.sqrt(A * R):.5f} "
              f"[{out[-1]['secs']:.0f}s]", flush=True)
        with open(WORK / f"fine_{tag}.json", "w") as fh:
            json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
