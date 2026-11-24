"""415pp — the LODO-V1 per-target scale artifact built on the SHIPPED composition ("PP"), 410's exact format.

Same key as build_weights415.py (LODO-V1: mean_j <P^{-j}_p, d_{j,p}> / |P_p|^2, raw bank rows, carried = >= 2 donors,
linear LOO with P's constants), but P is the bulk row the builder actually scales (sens415.py `compose_pp`): normrel cap
-1.5 donor weighting with PM50 at the uniform K562 share on every PM-observed gene, v9 whitening with the panel columns
kept, CD4 SVD-5 deflated at f 0.34 (all rows), H1 union g 0.28 and CD4 partial rows g 0.5 (carrying), one preserve.

Known answers (abort on failure):
  KA-410  410's key/weights recomputed == 410's shipped artifacts (exact)
  KA-PP   the PP constants reproduce the 410b build meta (whitening medians, operator preserve, three fractions,
          post-add preserve) to <= 1e-5 relative, and the CD4 deflation's U_k sha256 == meta (bit-exact)
PLOLO support for the PP-built key: work/sens415.json (beats K410 on all three cancer folds at both clips).

usage: uv run python experiments/415-lodo-key/build_weights415pp.py
writes work/pert_scale_415pp_lo0.3_hi3.npz, work/pert_scale_415pp_lo0.3_hi1e+09.npz (+ .json sidecars)
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from sens415 import (DEFLATE_K, MAIN, META, ROOT, WORK, Donors, cd4_frozen_and_partial, compose_pp,  # noqa: E402
                     deflate, lodo_key, loo_rows, plolo, source_residual, union_rows, UNION, OFFICIAL, E062)
from build_weights415 import A410, OP_PATH, OP_SHA, RULES  # noqa: E402
from key410 import shared_energy_key, weights  # noqa: E402
from lodo415 import spearman  # noqa: E402
from size410 import donor_rows  # noqa: E402


def main() -> None:
    targets = [str(t) for t in pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"]]
    assert len(targets) == 300
    wcols = np.asarray(np.load(ROOT / "experiments/162-retrodicting-validator/work/fold-kolf2.npz")["wcols"])
    op = plolo.Operator(OP_PATH, name="union_pruned")
    assert op.sha256 == OP_SHA
    # OCT-22 GUARD (DEF-RUNBOOK 8.3), verbatim from build_weights410.py.
    in_axis = sorted(set(targets) & {plolo.OFFICIAL[i] for i in wcols})
    if in_axis:
        raise SystemExit(f"{len(in_axis)} panel genes sit on the operator axis (e.g. {in_axis[:5]}): this operator was fitted "
                         "for a different panel. Refit it for this panel (DEF-RUNBOOK) and point OP_PATH/wcols at the refit.")
    pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([pos[g] for g in targets if g in pos], dtype=np.int64)

    # KA-410
    k410, car410, nd410 = shared_energy_key(targets, MAIN, op, wcols, cols=None)
    for (lo, hi), (fname, sha) in A410.items():
        path = ROOT / "experiments/410-pertscale-revisit/work" / fname
        assert hashlib.sha256(path.read_bytes()).hexdigest() == sha
        z = np.load(path)
        w410, _ = weights(k410, car410, lo, hi)
        assert np.abs(z["key"] - np.nan_to_num(k410, nan=0.0)).max() < 1e-15 and np.abs(z["w"] - w410).max() < 1e-12
    print("KA-410 PASS", flush=True)

    # the shipped composition + KA-PP
    dn = Donors(targets, MAIN)
    frozen, partial = cd4_frozen_and_partial(targets)
    cd4_def, daud = deflate(frozen, DEFLATE_K)
    assert daud["basis_u_k_sha256"] == META["cd4"]["deflation"]["basis_u_k_sha256"], "KA-PP deflation basis"
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    pre = dn.pre(MAIN, True)
    final, _, const, radd, c = compose_pp(op, pre, cd4_def, partial, h1, wcols, panel_positions)
    prod = dict(before_operator_median_norm=META["whitening"]["before_operator_median_norm"],
                operator_preserve_scale=META["whitening"]["operator_preserve_scale"],
                before_all_median_norm=META["whitening"]["before_all_median_norm"],
                after_all_median_norm=META["whitening"]["after_all_median_norm"],
                frac_cd4=META["cd4"]["fraction"],
                frac_h1=META["h1_union"]["residual_frame"]["applied_fraction"],
                frac_partial=META["cd4_partial"]["residual_frame"]["applied_fraction"],
                s=META["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"])
    rel = {k: const[k] / prod[k] - 1 for k in prod}
    worst = max(abs(v) for v in rel.values())
    assert worst < 1e-5, f"KA-PP failed: {rel}"
    print(f"KA-PP PASS worst relative constant error {worst:.2e}; deflation U_k sha matches", flush=True)

    carried = sum(dn.A[s].astype(int) for s in MAIN) >= 2
    assert np.array_equal(carried, car410)
    nd = sum(dn.A[s].astype(int) for s in MAIN).astype(np.int64)
    assert np.array_equal(nd, nd410.astype(np.int64))
    pp = (final * final).sum(1)
    num = loo_rows(op, dn, MAIN, True, c, const["s"], radd, wcols, panel_positions, None, None)
    d = {s: donor_rows(s, targets)[0] for s in MAIN}
    key = lodo_key(num, d, pp, dn.A, carried, MAIN)
    art228 = np.load(WORK / "pert_scale_415_lo0.3_hi3.npz")["key"]

    for lo, hi in RULES:
        w, med = weights(key, carried, lo, hi)
        w410, _ = weights(k410, car410, lo, hi)
        w228, _ = weights(np.where(carried, art228, np.nan), carried, lo, hi)
        tag = f"lo{lo:g}_hi{hi:g}"
        out = WORK / f"pert_scale_415pp_{tag}.npz"
        np.savez(out, targets=np.asarray(targets), w=w.astype(np.float64), key=np.nan_to_num(key, nan=0.0),
                 carried=carried, n_donors=nd)
        sha = hashlib.sha256(out.read_bytes()).hexdigest()
        q = lambda x: np.round(np.quantile(x, [0, .1, .25, .5, .75, .9, 1]), 4).tolist()
        side = dict(sha256=sha, lo=lo, hi=hi, donors=MAIN, operator_sha256=OP_SHA, gamma=0.65,
                    key="LODO-V1 on the SHIPPED composition (PP): mean_j <P^{-j}_p, d_j,p> / |P_p|^2, P = normrel cap "
                        "-1.5 + PM50 (061 broad genes) base, v9 whitening, CD4 deflated k5 f0.34 + H1 union g0.28 + CD4 "
                        "partial g0.5, one preserve; linear LOO with P's constants; d_j raw bank rows; all columns",
                    composition=const, composition_rel_to_410b_meta=rel, deflation=daud,
                    key_median_carried=med, carried=int(carried.sum()),
                    n_donor_hist={int(k): int((nd == k).sum()) for k in np.unique(nd)},
                    key_quantiles_carried=q(key[carried]), frac_key_le0=float((key[carried] <= 0).mean()),
                    w_quantiles=q(w), frac_at_floor=float((w[carried] <= lo + 1e-12).mean()),
                    frac_at_cap=float((w[carried] >= hi - 1e-12).mean()),
                    vs_410=dict(spearman_key=spearman(key[carried], k410[carried]),
                                frac_w_changed_gt_25pct=float((np.abs(np.log(w / w410)) > np.log(1.25)).mean()),
                                median_abs_log_ratio=float(np.median(np.abs(np.log(w / w410)))),
                                artifact_410=A410[(lo, hi)][0], artifact_410_sha256=A410[(lo, hi)][1]),
                    vs_415_228composition=dict(spearman_key=spearman(key[carried], art228[carried]),
                                               frac_w_changed_gt_25pct=float((np.abs(np.log(w / w228)) > np.log(1.25)).mean())),
                    top10=[(targets[i], round(float(key[i] / med), 2)) for i in np.argsort(-np.nan_to_num(key, nan=-9))[:10]],
                    bottom10=[(targets[i], round(float(key[i] / med), 2)) for i in np.argsort(np.nan_to_num(key, nan=9))[:10]],
                    validation="work/sens415.json (PLOLO on the PP prediction, three cancer folds + kolf2)")
        out.with_suffix(".json").write_text(json.dumps(side, indent=1))
        print(json.dumps({k: side[k] for k in ("sha256", "lo", "hi", "carried", "w_quantiles", "frac_at_floor",
                                               "frac_at_cap", "vs_410", "vs_415_228composition", "top10")}), flush=True)
    print(f"peak RSS {E062.peak_rss_gib():.2f} GiB")


if __name__ == "__main__":
    main()
