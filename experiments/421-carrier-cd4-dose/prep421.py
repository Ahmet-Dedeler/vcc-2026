#!/usr/bin/env python3
"""421 prep (copied from 416/419): gamma-0 carrier designs at a HIGHER CD4 dose (416's 2x2: CD4 is worth ~+0.006 of DE),
f = 0.34 x {1, 1.5, 2, 2.5}; f 0.34 must reproduce 416's g000 sha (known answer). 416 prep: the CARRIER effect matrices for the H1 menu, composed exactly like 323's cap15 (stage -> whiten -> CD4),
one variable each, every row norm-matched to the shipped cap15 row (so only the direction changes).

  gXXX     323's pre-whitening cap15 design whitened at gamma XXX/100 with 323's H1 operator, then 323's frozen CD4
           residual (its dose re-solved against that base, as production does). g065 must reproduce effects.npz cap15
           bit-for-bit (known answer for the composition path; compose.json final_sha256).
  gain     cap15 x 228's per-gene gain (gain-kolf2.npz: fitted without kolf2 and without h1_2025; gain_reg, t 0.25,
           m = 1 + t (g / median_live - 1), 1 off the live set) -- 410d's recipe applied to the carriers.
  signrep  387's VCC_SIGN_REPAIR=cons on the carriers: whitened (0.65) entries whose sign disagrees with the pre-whitening
           consensus take the consensus sign (magnitude kept), before the CD4 residual (the builder's order).
Writes work/carriers421.npz + work/carriers421.json. Reads only; no builder edits.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT / "experiments/323-current-family-validation", ROOT / "src"):
    sys.path.insert(0, str(p))

import numpy as np  # noqa: E402
import run323 as R  # noqa: E402  (read-only)

GAMMAS = (0.65,)  # kept only for the cap15 known answer
CD4_FS = (0.34, 0.51, 0.68, 0.85)
GAIN = ROOT / "experiments/228-theta-gain-plolo/work/gain-kolf2.npz"
GAIN_T = 0.25


def main():
    b = R.builder("cap15")
    comp = json.loads((R.WORK / "compose.json").read_text())
    assert R.sha(R.OPERATOR) == comp["new_operator_sha256"], "H1 operator drift"
    z = np.load(R.WORK / "stage.npz", allow_pickle=False)
    genes, targets, cols = z["genes"], z["targets"], z["columns"]
    with np.load(R.WORK / "effects.npz", allow_pickle=False) as ze:
        eff = ze["cap15"].copy()
        assert (ze["genes"] == genes).all() and (ze["targets"] == targets).all()
    assert b.array_sha256(eff) == comp["arms"]["cap15"]["final_sha256"]
    model = SimpleNamespace(
        names=list(R.DONORS), target_genes=genes,
        tcols={n: np.flatnonzero(np.isin(genes, (ROOT / f"data/donors/{n}/genes.txt").read_text().splitlines()))
               for n in R.DONORS})
    op = b.load_production_operator(model, targets, R.OPERATOR, R.OPERATOR.with_suffix(".json"),
                                    expected_sha256=R.sha(R.OPERATOR))
    assert np.array_equal(cols, op.columns)
    cd4 = b.EXP062.CD4Prepared(mode="gene_norm", mean_effect=z["cd4_mean"], gene_votes=z["cd4_votes"],
                              target_condition_count=z["cd4_counts"], condition_details=[])
    panel_positions = np.flatnonzero(np.isin(genes, targets))
    pre = np.asarray(z["cap15"], dtype=np.float64)
    ref_norm = np.linalg.norm(eff, axis=1)

    def norm_match(x):
        s = ref_norm / np.maximum(np.linalg.norm(x, axis=1), 1e-12)
        return x * s[:, None], s

    out, rec = {}, {"eff_sha256": b.array_sha256(eff), "operator_sha256": R.sha(R.OPERATOR), "designs": {}}
    for g in GAMMAS:
        w, wa = b.apply_v9_whitening(pre, op.whitener, cols, gamma=g, lam=0.0)
        final, ca, ua = b.apply_frozen_cd4_residual(w, cd4, cols, panel_positions)
        assert ua is None
        name = f"g{int(round(g * 100)):03d}" if g >= 0 else f"gm{int(round(-g * 100)):03d}"
        if name == "g065":
            ok = b.array_sha256(final) == comp["arms"]["cap15"]["final_sha256"]
            rec["known_answer_g065_reproduces_cap15"] = bool(ok)
            print("KNOWN ANSWER g065 == cap15:", ok, flush=True)
            if not ok:
                raise SystemExit("g065 composition does not reproduce 323's cap15")
            continue
        nm, s = norm_match(final)
        cosrow = (nm * eff).sum(1) / (np.linalg.norm(nm, axis=1) * ref_norm)
        out[name] = nm
        rec["designs"][name] = dict(gamma=g, cd4_fraction=float(ca["fraction"]), norm_scale_q=np.quantile(s, [0, .5, 1]).tolist(),
                                    row_cos_vs_cap15_q=np.quantile(cosrow, [0, .1, .5, .9, 1]).tolist(),
                                    sha256=b.array_sha256(nm))
        print(name, json.dumps(rec["designs"][name]), flush=True)
    # 421: gamma-0 carriers at several CD4 doses. The CD4 target f is a MODULE global read at import (VCC_CD4_TARGET_F), so each
    # dose loads its own builder instance under 323's env with only VCC_CD4_TARGET_F changed.
    import importlib.util as _u
    for f in CD4_FS:
        with R.model_env("cap15"):
            os.environ["VCC_CD4_TARGET_F"] = str(f)
            spec = _u.spec_from_file_location(f"builder421_{f}", ROOT / "scripts/make_submission_pm_cd4.py")
            bf = _u.module_from_spec(spec); sys.modules[spec.name] = bf; spec.loader.exec_module(bf)
        w0, _ = bf.apply_v9_whitening(pre, op.whitener, cols, gamma=0.0, lam=0.0)
        final, ca, ua = bf.apply_frozen_cd4_residual(w0, cd4, cols, panel_positions)
        assert ua is None
        nm, s = norm_match(final)
        name = f"c{int(round(f * 100)):03d}"
        cosrow = (nm * eff).sum(1) / (np.linalg.norm(nm, axis=1) * ref_norm)
        out[name] = nm
        rec["designs"][name] = dict(gamma=0.0, cd4_target_f=f, cd4_fraction=float(ca["fraction"]), norm_scale_q=np.quantile(s, [0, .5, 1]).tolist(),
                                    row_cos_vs_cap15_q=np.quantile(cosrow, [0, .1, .5, .9, 1]).tolist(), sha256=b.array_sha256(nm))
        print(name, json.dumps(rec["designs"][name]), flush=True)
    g416 = json.loads((ROOT / "experiments/416-de-carrier-menu/work/carriers416.json").read_text())["designs"]["g000"]["sha256"]
    ok = rec["designs"]["c034"]["sha256"] == g416
    rec["known_answer_c034_eq_416_g000"] = bool(ok)
    print("KNOWN ANSWER c034 == 416 g000:", ok, flush=True)
    if not ok:
        raise SystemExit("c034 does not reproduce 416's g000")
    np.savez(HERE / "work/carriers421.npz", genes=genes, targets=targets, **out)
    rec["npz_sha256"] = hashlib.sha256((HERE / "work/carriers421.npz").read_bytes()).hexdigest()
    (HERE / "work/carriers421.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
