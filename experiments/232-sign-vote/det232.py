"""232 — the determinism check, in a SECOND, INDEPENDENT process launch (228 §9's protocol).

Rebuilds one arm from scratch and compares every reported quantity, plus the per-gene
panel-total digest, against the stored row.  Nothing is cached between the two runs except
the fold and source `.npz` files, which are inputs, not results.

    uv run python det232.py 'production|vote_z@0' orion_hct116
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from arm232 import _mask_target_counts, score_arm  # noqa: E402
from common232 import (WORK, FoldBundle, control_block, ctx_gate, de_pools,  # noqa: E402
                       flip_mask, impose_sign, load_sources, stack_sources, vote_bundle)

KEYS = ("pds", "cos", "cos_agg", "sgn_w", "reach_prec89", "fid_precision",
        "delivered_median_norm", "digest")


def main() -> None:
    arm = sys.argv[1] if len(sys.argv) > 1 else "production|vote_z@0"
    held = sys.argv[2] if len(sys.argv) > 2 else "orion_hct116"
    stored = {(r["fold"], r["arm"]): r
              for r in json.loads((WORK / "arm232.json").read_text())}[(held, arm)]

    ctrl = {c: control_block(c) for c in ("A", "B", "C")}
    bctx, gate = {}, {}
    for c in ("A", "B", "C"):
        bctx[c], gate[c] = ctx_gate(c)

    b = FoldBundle(held)
    src = load_sources(held)
    pools = de_pools(b)
    base = b.final()
    tag, _, spec = arm.partition("|")
    name, _, tt = spec.partition("@")
    _, E, avail, supp = stack_sources(src, use_extra=(tag != "production"))
    vb = vote_bundle(E, avail, supp)
    del E
    m, s = flip_mask(vb, "z" if name.endswith("_z") else "maj", float(tt))
    changed, _ = _mask_target_counts(m, s, base)
    fin = impose_sign(base, m, s)
    e = score_arm(b, fin, pools, ctrl, gate, bctx)
    e.pop("_block_digests")

    ok = True
    print(f"determinism re-run: {held} / {arm}")
    print(f"  {'n_flipped':<26}{stored['n_flipped']:>22d}{int(changed.sum()):>22d}  "
          f"{'OK' if int(changed.sum()) == stored['n_flipped'] else 'MISMATCH'}")
    ok &= int(changed.sum()) == stored["n_flipped"]
    for k in KEYS:
        a, c = stored[k], e[k]
        same = (a == c) if isinstance(a, str) else (abs(a - c) == 0.0)
        ok &= same
        fmt = (f"  {k:<26}{a:>22}{c:>22}" if isinstance(a, str)
               else f"  {k:<26}{a:>22.15f}{c:>22.15f}")
        print(fmt + f"  {'OK' if same else 'MISMATCH'}")
    for k in ("nmae_pools", "nmae_U", "nmae_M", "q_w_pool", "jac"):
        for kk in stored[k]:
            same = abs(stored[k][kk] - e[k][kk]) == 0.0
            ok &= same
            print(f"  {k + '.' + kk:<26}{stored[k][kk]:>22.15f}{e[k][kk]:>22.15f}  "
                  f"{'OK' if same else 'MISMATCH'}")
    print(f"\n{'DETERMINISM PASS -- bit-identical' if ok else 'DETERMINISM FAIL'}")
    (WORK / "det232.json").write_text(json.dumps(
        dict(fold=held, arm=arm, PASS=bool(ok),
             rerun={k: e[k] for k in KEYS}), indent=1))
    assert ok


if __name__ == "__main__":
    main()
