#!/usr/bin/env python3
"""630 -- `VCC_DE_UPRN=<s|b>_k<K*100>`, the rank-neutral UP magnitude pass (after zs), as an UNAPPLIED builder patch.

The coordinator applies it AFTER 612: 538 -> 554 -> 590 (10-05 bundle) -> 598 -> 612 -> 630. It refuses a builder without
612's mark, and a second apply.

Why (council 600, [d4cd5f]; experiments/630-rank-neutral-up/RESULTS.md): on h1-test, with the real production emitter and
446bj2z's construction, b_k100 moves nmae -0.0066 raw (f 1; 3-seed f 0.35 overall-equiv +0.0030 [+0.0014, +0.0044]) with
fid, jac, reach, pds and mse EXACTLY unchanged and 0 call flips.

What it changes:
  1. scripts/make_submission_pm_cd4.py:
     K. after VCC_DE_MAG's parse: the knob (needs VCC_DE_MAG; `<s|b>_k<K*100>`, K in [0.5, 3]) and the module's sha pin;
     C. per context, after zs's scorer-order control block: the sha check, the scorer's tie sums, the median control
        library and the kept-axis position map (only when the knob is set);
     B. per block, right after the zs pass: `rank_up.apply(block, n_car, rc_scorer_control, DE_UPRN, ...)`;
     M. the meta gains `de_uprn` ONLY when the knob is set.
  2. src/vcc2026/rank_up.py = experiments/630-rank-neutral-up/src_rank_up.py (installed; refused if a different file is there).
Knob unset: no code runs and nothing is written, so the config and every block are unchanged (KA0: unset rebuilds the base
byte for byte).

    uv run python experiments/630-rank-neutral-up/patch_builder630.py --root <repo> --check   # anchors + order; writes nothing
    uv run python experiments/630-rank-neutral-up/patch_builder630.py --root <repo>           # apply (refuses a second time)
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
BUILDER = "scripts/make_submission_pm_cd4.py"
MODULE_DST = "src/vcc2026/rank_up.py"
MODULE_SRC = HERE / "src_rank_up.py"
MODULE_SHA = hashlib.sha256(MODULE_SRC.read_bytes()).hexdigest()
MARK = "630: rank-neutral UP"
MARK612 = "612: key gate"

EDITS = [
    # K. the knob, right after VCC_DE_MAG's
    ("    DE_MAG = DE_MAG_ENV\n",
     "    DE_MAG = DE_MAG_ENV\n"
     "# `VCC_DE_UPRN=<s|b>_k<K*100>[_q<n>r<n>p<n>]` -- 630: rank-neutral UP magnitude (experiments/630-rank-neutral-up), after the zs pass.\n"
     "# For uncalled genes the design signs UP: one carrier per gene is raised to just above the control max and receives the\n"
     "# gene's counts from the ballast rows; each count is paid by a sub-gate count of that carrier moved to the same ballast\n"
     "# row, so every per-gene total (pds, mse's profile half) and every library is exact; each gene's net rank-sum change\n"
     "# stays within |dAUC| <= 1e-4 (s: ballast donors inside their control intervals; b: budgeted crossings) and a call-flip\n"
     "# guard re-runs without any gene whose BH call flipped (0 flips). Target t = max(l, K * 0.24 / sqrt(control count per\n"
     "# cell)). Module src/vcc2026/rank_up.py, sha pinned. Needs VCC_DE_MAG. Unset is a no-op.\n"
     "DE_UPRN_ENV = os.environ.get(\"VCC_DE_UPRN\", \"\").strip()\n"
     "DE_UPRN = None\n"
     "DE_UPRN_MODULE_PATH = ROOT / \"src/vcc2026/rank_up.py\"\n"
     f"DE_UPRN_MODULE_SHA256 = \"{MODULE_SHA}\"\n"
     "if DE_UPRN_ENV:\n"
     "    import re as _re_ur  # noqa: PLC0415\n"
     "    _ur = _re_ur.fullmatch(r\"(s|b)_k(\\d{2,3})(?:_q([1-9]\\d*)r([1-9]\\d*)p([1-9]\\d*))?\", DE_UPRN_ENV)\n"
     "    if _ur is None or not 50 <= int(_ur.group(2)) <= 300:\n"
     "        raise SystemExit(f\"VCC_DE_UPRN must look like '<s|b>_k<K*100>[_q<quantum>r<rounds>p<pool>]' with K in [0.5, 3], e.g. b_k100_q20r3p10: {DE_UPRN_ENV!r}\")\n"
     "    if DE_MAG is None:\n"
     "        raise SystemExit(\"VCC_DE_UPRN needs VCC_DE_MAG (it runs after the zs pass, on its evaluator-order control block)\")\n"
     "    DE_UPRN = DE_UPRN_ENV\n"),
    # C. per context: sha pin, the scorer's tie sums, the median control library, the kept-axis position map
    ("                    memory_checkpoint(f\"de_mag scorer control block {context}\", max_rss_gib)\n",
     "                    memory_checkpoint(f\"de_mag scorer control block {context}\", max_rss_gib)\n"
     "                    if DE_UPRN is not None:\n"
     "                        _ur_sha = hashlib.sha256(DE_UPRN_MODULE_PATH.read_bytes()).hexdigest()\n"
     "                        if _ur_sha != DE_UPRN_MODULE_SHA256:\n"
     "                            raise SystemExit(f\"rank_up.py sha256 {_ur_sha} != pinned {DE_UPRN_MODULE_SHA256}\")\n"
     "                        from vcc2026 import rank_up as _rank_up  # noqa: PLC0415\n"
     "                        rc_uprn_tie = _rank_up.scorer_tie(rc_scorer_control)\n"
     "                        rc_uprn_lib_med = float(np.median(libraries))\n"
     "                        rc_uprn_kpos = np.cumsum(np.asarray(rc_scorer_control.keep, dtype=bool)) - 1\n"
     "                        rc_uprn_gpos = {str(g): i for i, g in enumerate(spec.genes)}\n"
     "                        print(f\"[de_uprn] {context}: {DE_UPRN}; median control library {rc_uprn_lib_med:.0f}\", flush=True)\n"),
    # B. per block, right after zs
    ("                                rc_diag = {**rc_diag, **_dm_diag}\n",
     "                                rc_diag = {**rc_diag, **_dm_diag}\n"
     "                            if DE_UPRN is not None:\n"
     "                                # 630: rank-neutral UP magnitude after zs (see the knob's comment); asserts totals,\n"
     "                                # libraries, per-gene U within the cap, and 0 call flips\n"
     "                                from vcc2026 import rank_up as _rank_up  # noqa: PLC0415\n"
     "                                _ur_full = rc_uprn_gpos.get(str(target), -1)\n"
     "                                _ur_own = (int(rc_uprn_kpos[_ur_full]) if _ur_full >= 0 and bool(np.asarray(\n"
     "                                    rc_scorer_control.keep, dtype=bool)[_ur_full]) else -1)\n"
     "                                block, _ur_diag = _rank_up.apply(np.asarray(block), int(spec.cells_per_pert) - int(ballast_k),\n"
     "                                                                 rc_scorer_control, DE_UPRN, tie=rc_uprn_tie,\n"
     "                                                                 lib_med=rc_uprn_lib_med, own_col=_ur_own)\n"
     "                                rc_diag = {**rc_diag, **_ur_diag}\n"),
    # M. the meta record
    ("                                                             \"(evaluator-order control block)\"}}),\n",
     "                                                             \"(evaluator-order control block)\"}}),\n"
     "            **({} if DE_UPRN is None else {\"de_uprn\": {\"variant\": DE_UPRN, \"module_sha256\": DE_UPRN_MODULE_SHA256,\n"
     "                                                      \"where\": \"after the zs pass; one promoted carrier per uncalled UP \"\n"
     "                                                               \"gene fed from the ballast rows; totals, libraries exact; \"\n"
     "                                                               \"|dAUC| <= 1e-4 per gene; 0 call flips\"}}),\n"),
]


def patched_text(text: str) -> str:
    if MARK in text:
        raise SystemExit(f"{BUILDER}: already carries {MARK!r} (a second apply is refused)")
    if MARK612 not in text:
        raise SystemExit(f"{BUILDER}: lacks 612's mark {MARK612!r}: apply 598 and 612 first "
                         f"(order 538 -> 554 -> 590 -> 598 -> 612 -> 630)")
    for old, new in EDITS:
        if text.count(old) != 1:
            raise SystemExit(f"{BUILDER}: anchor found {text.count(old)} times (need exactly 1): {old[:80]!r}")
        text = text.replace(old, new, 1)
    if MARK not in text:
        raise SystemExit("internal: the mark is not in the patched text")
    return text


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--check", action="store_true", help="verify anchors and order; write nothing")
    a = ap.parse_args()
    p = a.root / BUILDER
    new = patched_text(p.read_text())
    compile(new, str(p), "exec")
    dst = a.root / MODULE_DST
    if dst.exists() and hashlib.sha256(dst.read_bytes()).hexdigest() != MODULE_SHA:
        raise SystemExit(f"{MODULE_DST} exists with a different sha: refusing to overwrite")
    if a.check:
        print(f"630 check OK: {BUILDER} carries 612's mark, all {len(EDITS)} anchors found once, the patched text compiles; "
              f"module sha {MODULE_SHA}")
        return 0
    fd, tmp = tempfile.mkstemp(dir=p.parent, prefix=".patch630-")
    with os.fdopen(fd, "w") as f:
        f.write(new)
    os.chmod(tmp, os.stat(p).st_mode)
    os.replace(tmp, p)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_bytes(MODULE_SRC.read_bytes())
    print(f"630 applied to {p}; installed {dst} (sha {MODULE_SHA})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
