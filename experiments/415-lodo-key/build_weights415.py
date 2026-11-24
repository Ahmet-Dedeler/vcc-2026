"""415 — the LODO per-target scale artifact for the REAL A/B/C panel, in 410's exact format (targets, w, key, carried,
n_donors npz + json sidecar with sha256). Written only because LODO-V1 beat 410's key on all three cancer PLOLO folds
(work/lodo415.json; PREREGISTRATION.md decision rule).

key_p = mean_j <P^{-j}_p, d_{j,p}> / |P_p|^2   over the production donors j carrying p (carried = >= 2 donors)
  P      = 228's production composition on the A/B/C panel: c*W(ALPHA*pre) on the operator axis (ALPHA*pre off it),
           + frac_cd4*CD4 (f 0.34, all-rows sizing) + frac_u*H1union (g 0.28, carrying sizing), post-add preserve s,
           panel zeroed; pre = build_baseline of the three prepared donor rows (k562 PM50-blended), W = union_pruned
           (sha e2f65ecf..., the file the A/B/C whitener stage reproduces)
  P^{-j} = the same with pre -> mean of the OTHER carrying donors' prepared rows, constants c, fracs, s fixed at P's
  d_j    = donor j's raw bank row (410's donor_rows), inner products over the whole official axis (410's cols=None)
w = clip(key / median(key over carried), lo, hi); targets carried by < 2 donors get 1.0.

Known answer first (aborts on failure): 410's key recomputed here reproduces 410's shipped artifacts' key and w arrays.
The builder is NOT touched; wiring this file into VCC_PERT_SCALE is a separate, reviewed step.

usage: uv run python experiments/415-lodo-key/build_weights415.py
writes work/pert_scale_415_lo0.3_hi3.npz, work/pert_scale_415_lo0.3_hi1e+09.npz (+ .json sidecars)
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
for p in ("experiments/410-pertscale-revisit", "experiments/162-retrodicting-validator",
          "experiments/218-independent-donors"):
    sys.path.insert(0, str(ROOT / p))
import plolo  # noqa: E402
from resid218 import source_residual  # noqa: E402
from compose218 import union_rows  # noqa: E402

plolo.WORK = WORK
from plolo import E062, OFFICIAL, cd4_for  # noqa: E402
from key410 import raw_whiten, shared_energy_key, weights  # noqa: E402
from size410 import donor_rows  # noqa: E402
from lodo415 import ALPHA, PREP_NAME, UNION, compose, whiten_split, spearman  # noqa: E402

DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
OP_PATH = ROOT / "data/donors/whitener_union_pruned.npz"
OP_SHA = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"
RULES = ((0.3, 3.0), (0.3, 1e9))
A410 = {(0.3, 3.0): ("pert_scale_410_lo0.3_hi3.npz", "d7f76a22f62ad6991d7220ce8a65218db38b2d71dcfe3c3781ec6a8f70409e0e"),
        (0.3, 1e9): ("pert_scale_410_lo0.3_hi1e+09.npz", "94c048e99aecd40004e70f559fe26bdd3b8b9aef9cb8c8cd0efbf8a72d2fe69e")}


def main() -> None:
    targets = [str(t) for t in pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"]]
    assert len(targets) == 300
    wcols = np.asarray(np.load(ROOT / "experiments/162-retrodicting-validator/work/fold-kolf2.npz")["wcols"])
    op = plolo.Operator(OP_PATH, name="union_pruned")
    assert op.sha256 == OP_SHA
    # OCT-22 GUARD (DEF-RUNBOOK 8.3), verbatim from build_weights410.py: the operator axis is union(donor coverage) -
    # PANEL. A new panel on the old axis puts panel genes inside wcols and every weight silently comes out on the wrong
    # geometry. Refuse instead.
    in_axis = sorted(set(targets) & {plolo.OFFICIAL[i] for i in wcols})
    if in_axis:
        raise SystemExit(f"{len(in_axis)} panel genes sit on the operator axis (e.g. {in_axis[:5]}): this operator was fitted "
                         "for a different panel. Refit it for this panel (DEF-RUNBOOK) and point OP_PATH/wcols at the refit.")
    pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([pos[g] for g in targets if g in pos], dtype=np.int64)

    # ---------------- KA: 410's key and weights, recomputed, equal 410's shipped artifacts
    k410, car410, nd410 = shared_energy_key(targets, DONORS, op, wcols, cols=None)
    for (lo, hi), (fname, sha) in A410.items():
        path = ROOT / "experiments/410-pertscale-revisit/work" / fname
        assert hashlib.sha256(path.read_bytes()).hexdigest() == sha, f"410 artifact {fname} moved"
        z = np.load(path)
        w410, _ = weights(k410, car410, lo, hi)
        assert [str(t) for t in z["targets"]] == targets
        assert np.array_equal(z["carried"], car410) and np.array_equal(z["n_donors"], nd410.astype(np.int64))
        dk = float(np.abs(z["key"] - np.nan_to_num(k410, nan=0.0)).max())
        dw = float(np.abs(z["w"] - w410).max())
        assert dk < 1e-15 and dw < 1e-12, f"KA 410 {fname}: key {dk} w {dw}"
        print(f"KA PASS 410 artifact {fname}: max|dkey| {dk:.1e} max|dw| {dw:.1e}", flush=True)

    # ---------------- the production composition on the A/B/C panel
    E, A = {}, {}
    mean_on_axis = {}
    for s in DONORS:
        z = np.load(WORK / f"prep-{PREP_NAME[s]}.npz", allow_pickle=True)
        assert [str(t) for t in z["targets"]] == targets
        E[s], A[s], mean_on_axis[s] = z["effects"], z["available"], z["mean_on_axis"]
    votes = sum(A[s].astype(int) for s in DONORS)
    assert (votes > 0).all(), "a panel target has no donor row (build_baseline fallback not expected on A/B/C)"
    pre = np.zeros_like(E[DONORS[0]])
    for s in DONORS:
        pre[A[s]] += E[s][A[s]]
    pre /= votes[:, None]
    cd4 = cd4_for(targets)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    final, const, radd, whitened = compose(op, pre, cd4, h1, wcols, panel_positions)
    c, WM = whiten_split(op, pre * ALPHA, wcols)
    chk = pre * ALPHA
    chk[:, wcols] = c * WM[:, wcols]
    assert np.allclose(chk, whitened, rtol=1e-10, atol=1e-12)
    s_add = const["s"]
    const["c"] = c
    const["whitened_median_row_norm"] = float(np.median(np.linalg.norm(whitened, axis=1)))
    Rr = s_add * radd
    Rr[:, panel_positions] = 0.0
    const["residual_norm_share_median"] = float(np.median(np.linalg.norm(Rr, axis=1) / np.linalg.norm(final, axis=1)))
    const["h1_carrying_targets"] = int((np.asarray(h1.target_condition_count) > 0).sum())
    print(f"composition: {json.dumps({k: round(v, 5) for k, v in const.items()})}", flush=True)

    carried = votes >= 2
    assert np.array_equal(carried, car410)
    raw = {s: donor_rows(s, targets) for s in DONORS}
    for s in DONORS:
        assert np.array_equal(raw[s][1], A[s]), f"availability mismatch {s}"
    pp = (final * final).sum(1)
    num = {}
    for j in DONORS:
        prej = np.zeros_like(pre)
        cnt = np.zeros(len(targets), int)
        for i in DONORS:
            if i != j:
                prej[A[i]] += E[i][A[i]]
                cnt += A[i].astype(int)
        ok = cnt > 0
        prej[ok] /= cnt[ok, None]
        bj = prej * ALPHA
        bj[:, wcols] = c * raw_whiten(op, bj, wcols)[:, wcols]
        Pj = s_add * (bj + radd)
        Pj[:, panel_positions] = 0.0
        num[j] = (Pj * raw[j][0]).sum(1)
    key = np.full(len(targets), np.nan)
    for p in np.flatnonzero(carried):
        key[p] = float(np.mean([num[j][p] for j in DONORS if A[j][p]])) / pp[p]
    nd = votes.astype(np.int64)

    for lo, hi in RULES:
        w, med = weights(key, carried, lo, hi)
        w410, med410 = weights(k410, car410, lo, hi)
        tag = f"lo{lo:g}_hi{hi:g}"
        out = WORK / f"pert_scale_415_{tag}.npz"
        np.savez(out, targets=np.asarray(targets), w=w.astype(np.float64), key=np.nan_to_num(key, nan=0.0),
                 carried=carried, n_donors=nd)
        sha = hashlib.sha256(out.read_bytes()).hexdigest()
        q = lambda x: np.round(np.quantile(x, [0, .1, .25, .5, .75, .9, 1]), 4).tolist()
        side = dict(sha256=sha, lo=lo, hi=hi, donors=DONORS, operator_sha256=OP_SHA, gamma=0.65,
                    key="LODO-V1: mean_j <P^{-j}_p, d_j,p> / |P_p|^2, P = 228 production composition (whitened base + "
                        "CD4 f0.34 + H1 union g0.28), linear LOO with P's constants, d_j raw bank rows, all columns",
                    composition=const,
                    key_median_carried=med, carried=int(carried.sum()),
                    n_donor_hist={int(k): int((nd == k).sum()) for k in np.unique(nd)},
                    key_quantiles_carried=q(key[carried]), frac_key_le0=float((key[carried] <= 0).mean()),
                    w_quantiles=q(w), frac_at_floor=float((w[carried] <= lo + 1e-12).mean()),
                    frac_at_cap=float((w[carried] >= hi - 1e-12).mean()),
                    vs_410=dict(spearman_key=spearman(key[carried], k410[carried]),
                                pearson_logw=float(np.corrcoef(np.log(w), np.log(w410))[0, 1]),
                                frac_w_changed_gt_25pct=float((np.abs(np.log(w / w410)) > np.log(1.25)).mean()),
                                median_abs_log_ratio=float(np.median(np.abs(np.log(w / w410)))),
                                artifact_410=A410[(lo, hi)][0], artifact_410_sha256=A410[(lo, hi)][1]),
                    top10=[(targets[i], round(float(key[i] / med), 2)) for i in np.argsort(-np.nan_to_num(key, nan=-9))[:10]],
                    bottom10=[(targets[i], round(float(key[i] / med), 2)) for i in np.argsort(np.nan_to_num(key, nan=9))[:10]],
                    validation="work/lodo415.json (PLOLO, three cancer folds + kolf2)")
        out.with_suffix(".json").write_text(json.dumps(side, indent=1))
        print(json.dumps(side, indent=1), flush=True)
    print(f"peak RSS {E062.peak_rss_gib():.2f} GiB")


if __name__ == "__main__":
    main()
