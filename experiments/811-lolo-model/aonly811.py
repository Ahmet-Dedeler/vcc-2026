"""811 LABC instrument read (lab lead, 2026-10-16 ~10:50Z): c0_l25, c1_l25 and src_l25 vs inc on 801's LABC/<held> folds
(A/B/C's 300 targets built with the fresh-fold L recipe, D/E/F rules, inc = pct+P2g).

A WRAPPER, so the frozen m811.py (sha256 e00192e6...) stays byte-identical. As frozen, m811.main treats any panel starting
with "L" as an L fold and requires 804 / base811 reference arrays that do not exist for LABC (it would void every fold).
This file adds NO arithmetic: it calls m811.fold_U (c0, c1), base811.consensus over the fold's sources (src, exactly as
base811 built src_lX), base811.blend, engine801.named("pct+P2g") and engine801.price, and saves with ka804.save.
Known answers: the engine's own fold KAs (E.load check=True) on every fold; and on the known-answer fold L1/k562_gwps the
wrapper must reproduce m811's stage-1 c0_l25 / c1_l25 arrays and base811's inc / src_l25 arrays at |d| <= 1e-12.

    uv run python experiments/811-lolo-model/labc811.py <fold-id>          # writes $AONLY811_OUT/{bulk,meta}
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import m811 as M  # noqa: E402

B, E, K = M.B, M.E, M.K
OUT = Path(os.environ.get("AONLY811_OUT", "/work/aonly811"))
KA_FOLD = "L1/k562_gwps"
KA_REF = {"c0_l25": HERE / "work/run-m811s1/bulk", "c1_l25": HERE / "work/run-m811s1/bulk",
          "inc": HERE / "work/run-base811/bulk", "src_l25": HERE / "work/run-base811/bulk"}


def main(fold: str):
    t0 = time.time()
    log = lambda *x: print(f"[{fold}]", *x, f"({time.time() - t0:.0f}s)", flush=True)
    assert fold[:1] == "A" and fold[1:2].isdigit() or fold == KA_FOLD, "aonly811: A<k>/ folds (+ the KA fold) only"; B.BANKS = HERE / "work" / "banks_ext"
    F = E.load(fold, check=True)
    K.OUT = OUT
    inc0 = E.named(F, "pct+P2g")
    inc = E.Arm(w=inc0.w, q=inc0.q, name="inc")
    fin = np.asarray(F.fin, np.float64)
    npn = F.nonpanel
    U0, info0 = M.fold_U(F, False)
    U1, info1 = M.fold_U(F, True)
    Us, _ = B.consensus(F, {f"src:{s}": (s,) for s in F.sources})
    arms = {"inc": inc,
            "c0_l25": E.Arm(w=inc.w, q=inc.q, fin=B.blend(fin, U0, 0.25, npn), name="c0_l25"),
            "c1_l25": E.Arm(w=inc.w, q=inc.q, fin=B.blend(fin, U1, 0.25, npn), name="c1_l25"),
            "src_l25": E.Arm(w=inc.w, q=inc.q, fin=B.blend(fin, Us, 0.25, npn), name="src_l25")}
    rec = dict(fold=fold, n=F.n, held=F.held, sources=F.sources, ka_engine=F.ka, arms={}, ka811={},
               U=dict(c0=M.shared_share(U0), c1=M.shared_share(U1), src=M.shared_share(Us)),
               script_sha256=B.sha(Path(__file__)), m811_sha256=B.sha(HERE / "m811.py"), base811_sha256=B.sha(HERE / "base811.py"))
    for nm, arm in arms.items():
        b = None if nm == "inc" else "inc"
        s, arr = E.price(F, arm, base=arms[b] if b else None)
        arr["targets"] = np.asarray(F.targets)
        rec["arms"][nm] = K.save(fold, nm, s, arr, b)
        if fold == KA_FOLD:
            keep = {k: v for k, v in arr.items() if not k.startswith("base_")}
            try:
                rec["ka811"][nm] = B.compare(fold, nm, keep, KA_REF[nm] / fold.replace("/", "__"))
                log(f"KA {nm} == {KA_REF[nm].parent.name}/{nm} (worst |d| {rec['ka811'][nm]['worst']:.1e})")
            except Exception as e:
                rec["ka811"][nm] = dict(FAILED=f"{type(e).__name__}: {e}")
                log(f"KA {nm} FAILED: {e}")
        if b:
            v = s["vs_base"]
            log(f"{nm:7s} r {100 * v['r']:+.3f}% dpds {v['dpds']:+.5f} dpdsl {v['dpds_lam05']:+.5f} dK {100 * v['dK_rel']:+.2f}% "
                f"off-axis {v['share_unscored']:.2f}")
    (OUT / "meta").mkdir(parents=True, exist_ok=True)
    (OUT / "meta" / f"{fold.replace('/', '__')}.json").write_text(json.dumps(rec, indent=1, default=float))
    bad = [k for k, v in rec["ka811"].items() if "FAILED" in v]
    assert not bad, f"KA FAILED on {bad}"
    log("done")


if __name__ == "__main__":
    main(sys.argv[1])
