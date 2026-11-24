#!/usr/bin/env python3
"""416 prep b (follow-ups written after seed 0 showed the carrier-gamma dose response; same composition code as prep416):
  w065nocd4  whitened at gamma 0.65, NO CD4 residual (with `pre` = gamma 0 without CD4 this completes the 2x2
             {whitened, unwhitened} x {CD4, no CD4} on the carriers)
  gm030      gamma -0.30 (anti-whitening: C^{+0.15}), + CD4 -- does the dose response continue past gamma 0?
Rows norm-matched to the shipped cap15 rows. Writes work/carriers416b.npz (+ .json); carriers416.npz is not touched.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT / "experiments/323-current-family-validation", ROOT / "src"):
    sys.path.insert(0, str(p))

import numpy as np  # noqa: E402
import run323 as R  # noqa: E402  (read-only)


def main():
    b = R.builder("cap15")
    comp = json.loads((R.WORK / "compose.json").read_text())
    assert R.sha(R.OPERATOR) == comp["new_operator_sha256"]
    z = np.load(R.WORK / "stage.npz", allow_pickle=False)
    genes, targets, cols = z["genes"], z["targets"], z["columns"]
    with np.load(R.WORK / "effects.npz", allow_pickle=False) as ze:
        eff = ze["cap15"].copy()
    assert b.array_sha256(eff) == comp["arms"]["cap15"]["final_sha256"]
    model = SimpleNamespace(
        names=list(R.DONORS), target_genes=genes,
        tcols={n: np.flatnonzero(np.isin(genes, (ROOT / f"data/donors/{n}/genes.txt").read_text().splitlines()))
               for n in R.DONORS})
    op = b.load_production_operator(model, targets, R.OPERATOR, R.OPERATOR.with_suffix(".json"),
                                    expected_sha256=R.sha(R.OPERATOR))
    cd4 = b.EXP062.CD4Prepared(mode="gene_norm", mean_effect=z["cd4_mean"], gene_votes=z["cd4_votes"],
                              target_condition_count=z["cd4_counts"], condition_details=[])
    panel_positions = np.flatnonzero(np.isin(genes, targets))
    pre = np.asarray(z["cap15"], dtype=np.float64)
    ref_norm = np.linalg.norm(eff, axis=1)

    def norm_match(x):
        return x * (ref_norm / np.maximum(np.linalg.norm(x, axis=1), 1e-12))[:, None]

    # known answer again: 0.65 + CD4 == cap15
    w65, _ = b.apply_v9_whitening(pre, op.whitener, cols, gamma=0.65, lam=0.0)
    f65, _, _ = b.apply_frozen_cd4_residual(w65, cd4, cols, panel_positions)
    assert b.array_sha256(f65) == comp["arms"]["cap15"]["final_sha256"], "composition known answer failed"
    out, rec = {}, {"known_answer_g065_reproduces_cap15": True, "designs": {}}
    out["w065nocd4"] = norm_match(w65)
    wm, wa = b.apply_v9_whitening(pre, op.whitener, cols, gamma=-0.30, lam=0.0)
    fm, ca, _ = b.apply_frozen_cd4_residual(wm, cd4, cols, panel_positions)
    out["gm030"] = norm_match(fm)
    for k, v in out.items():
        cos = (v * eff).sum(1) / (np.linalg.norm(v, axis=1) * ref_norm)
        rec["designs"][k] = dict(sha256=b.array_sha256(v), row_cos_vs_cap15_q=np.quantile(cos, [0, .1, .5, .9, 1]).tolist())
        print(k, json.dumps(rec["designs"][k]), flush=True)
    np.savez(HERE / "work/carriers416b.npz", genes=genes, targets=targets, **out)
    rec["npz_sha256"] = hashlib.sha256((HERE / "work/carriers416b.npz").read_bytes()).hexdigest()
    (HERE / "work/carriers416b.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
