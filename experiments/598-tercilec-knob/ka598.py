#!/usr/bin/env python3
"""598 KA-c checker: the tercilec build of 479p2g's configuration (RESULTS.md, KA plan, step 3).

    uv run python experiments/598-tercilec-knob/ka598.py experiments/527-arc-hook-arm/work/build_KA598c/ka527-KA598c.h5ad.meta.json

PASS requires all of:
  1. config["bulk_sharpen_pt"] == work/ka598c_expected.json (test598 T6, computed on the same dump with the applied
     builder's own functions): every non-float field exactly (mode, q, counts 98/104/98, floor_to_low, tie_to_mid,
     nonfinite, n_carried 294, noncarried_to_mid 6, moved_vs_tercile = exactly E2F4 2->1, EPHB2 0->1, SLC25A13 2->1,
     TMEM209 0->1, carried_key lodo415_lo0.5_hi3, carried_key_sha256 8c82bef3...); edges and key_ratio_c within 1e-12
     relative; key_ratio_maxdev <= 1e-12;
  2. the config differs from 479p2g's meta in bulk_sharpen_pt only (rule 10), plus de_sign_emitter_sha256 when the
     emitter moved since 479p2g (edb5b8dc then; 535's 0f56209d since 09-30);
  3. the output sha256 differs from 479p2g's dd82740d (four targets' q moved; an identical artifact would mean the knob
     never reached the render).
run_ka_527.sh's own line for this build reads "KA FAIL" (it compares against 479p2g's sha): that is expected here.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXPECTED = ROOT / "experiments/598-tercilec-knob/work/ka598c_expected.json"
REF479 = ROOT / "submissions/479p2g-cap15-partial50-clean4q80-pertscale3-ballast-sharpenp2g-stacklodo05load-ctxlevelr-mde010coarse-zr5.h5ad.meta.json"
FLOATS = {"edges", "key_ratio_c"}


def close(a, b) -> bool:
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(close(x, y) for x, y in zip(a, b))
    try:
        return abs(float(a) - float(b)) <= 1e-12 * max(1.0, abs(float(b)))
    except (TypeError, ValueError):
        return False


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    meta = json.loads(Path(sys.argv[1]).read_text())
    ref = json.loads(REF479.read_text())
    want = json.loads(EXPECTED.read_text())
    got = meta["config"].get("bulk_sharpen_pt") or {}
    bad = []
    for k, v in want.items():
        if k == "key_ratio_maxdev":
            ok = isinstance(got.get(k), (int, float)) and got[k] <= 1e-12
        elif k in FLOATS:
            ok = close(got.get(k), v)
        else:
            ok = got.get(k) == v
        if not ok:
            bad.append(f"{k}: got {json.dumps(got.get(k))} want {json.dumps(v)}")
    extra = sorted(set(got) - set(want))
    if extra:
        bad.append(f"unexpected audit keys {extra}")
    print(f"[1] bulk_sharpen_pt {'== expected' if not bad else 'DIFFERS'}: counts {got.get('counts')} moved {got.get('moved_vs_tercile')}")
    for b in bad:
        print(f"    {b}")
    c, r = meta["config"], ref["config"]
    diff = sorted(k for k in set(c) | set(r) if c.get(k) != r.get(k))
    allowed = {"bulk_sharpen_pt"} | ({"de_sign_emitter_sha256"} if c.get("de_sign_emitter_sha256") != r.get("de_sign_emitter_sha256") else set())
    ok2 = set(diff) == allowed
    print(f"[2] config diff vs 479p2g: {diff} ({'exactly the intended' if ok2 else 'NOT the intended ' + str(sorted(allowed))})")
    sa, sb = meta["output"]["sha256"], ref["output"]["sha256"]
    ok3 = sa != sb
    print(f"[3] output {sa[:16]} vs 479p2g {sb[:16]}: {'differs (expected)' if ok3 else 'IDENTICAL: the knob did not reach the render'}")
    ok = not bad and ok2 and ok3
    print("KA598c PASS" if ok else "KA598c FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
