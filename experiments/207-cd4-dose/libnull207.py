"""207 Phase 4 — the FORM's falsifiers: shuffled-CD4 library, and an unmatched-lineage library.

Same construction as `library207.py`, with the fourth member replaced by
  * `libshuf` -- CD4 with whole rows permuted across targets (falsifier 1, for the library form);
  * `libkolf` -- KOLF2, an iPSC line covering 273/300 of the panel, as the fourth member
                 (falsifier 2: is the gain LINEAGE, or is it just a fourth independent line?).
KOLF2 is never used on the fold where it is the truth.

Both are priced at their own re-solved alpha, exactly as the real library arms are.

    uv run python libnull207.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"

from plolo import FOLDS, PROD_OPS, Operator, build_fold, cd4_for  # noqa: E402
from sweep207 import (ALPHA, GAMMA, LAM, alt_line_residual, sha8,  # noqa: E402
                      score_effects, shuffled_cd4)
from library207 import votes_for  # noqa: E402

W = [0.25, 0.5, 0.75, 1.0]


def main() -> None:
    out_path = WORK / "libnull207.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}

    for held in FOLDS:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        votes = votes_for(held, fold["targets"])
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        n0 = float(np.median(np.linalg.norm(pre[:, fold["wcols"]], axis=1)))

        members = {"shuf": shuffled_cd4(cd4)}
        if held != "kolf2":
            members["kolf"] = alt_line_residual(cd4, "kolf2", fold["targets"])
        for tag, mem in members.items():
            eff = np.asarray(mem.mean_effect, dtype=np.float64)
            carry = np.asarray(mem.target_condition_count) > 0
            for w in W:
                name = f"lib{tag}{w:g}@a*"
                if (held, name) in done:
                    continue
                t0 = time.time()
                sel = carry & (votes > 0)
                pre4 = pre.copy()
                pre4[sel] = ((pre[sel] * votes[sel, None] + w * eff[sel])
                             / (votes[sel, None] + w))
                n1 = float(np.median(np.linalg.norm(pre4[:, fold["wcols"]], axis=1)))
                a = ALPHA * n0 / n1 if n1 else ALPHA
                final = op.apply(pre4 * a, fold["wcols"], GAMMA, LAM)
                final[:, fold["panel_positions"]] = 0.0
                e = score_effects(final, fold, truth)
                rows.append(dict(fold=held, arm=name, n_targets=len(fold["targets"]),
                                 note=f"falsifier library member {tag} w={w:g}",
                                 **e,
                                 audit=dict(lib_w=w, residual_f=0.0, alpha=a,
                                            pre_norm_ratio=n1 / n0, member=tag,
                                            operator=dict(path=str(op_path.relative_to(ROOT)),
                                                          sha8=sha8(op_path))),
                                 seconds=round(time.time() - t0, 1)))
                out_path.write_text(json.dumps(rows, indent=2))
                print(f"[{held}/{name}] a {a:.4f} pds {e['pds']:.6f} "
                      f"cos_agg {e['cos_agg']:.6f} q_w {e['a15_sgn_w']:.6f} "
                      f"prenorm {n1/n0:.4f}", flush=True)
        del op


if __name__ == "__main__":
    main()
