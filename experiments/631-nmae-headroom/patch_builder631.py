#!/usr/bin/env python3
"""631 -- `VCC_DE_DNMAG=d_m<m*100>[_q<n>]`, the DOWN magnitude pass (after zs and after 630's UP pass), as an UNAPPLIED
builder patch.

Order: ... -> 612 -> 630 -> 631 (it refuses a builder without 630's mark, and a second apply). When both knobs are set the
UP pass runs first and the DOWN pass reads its output (631: the two passes share one currency, and 630's flow doubles the
DOWN pass's delivery). Either knob works alone.

Why (council 600, [d4cd5f]; experiments/631-nmae-headroom/RESULTS.md): called DOWN genes are delivered at about half the
right magnitude. On h1-test's real emitter, ×2 on every called DOWN gene is worth −0.023 / −0.025 raw nmae at full
delivery. The pass delivers 14 % of it alone (+0.0008 overall-eq at f 0.35) and 28 % after 630 (630 + DN: +0.0037 vs
630's +0.0029 at f 0.35; +0.0049 vs +0.0034 cmp), with fid / jac unchanged, reach ≤ 0.0014 raw, and 0 call flips.

What it changes:
  1. scripts/make_submission_pm_cd4.py:
     K. after 630's knob: `VCC_DE_DNMAG` (needs VCC_DE_MAG; d_m<m*100>[_q<n>], m in [1.25, 4]) and the module sha pin;
     C. per context, after zs's scorer-order control block: the sha check, the scorer's tie sums, the kept-axis position map;
     B. per block, right after 630's UP pass (or after zs when 630 is unset): `dn_mag.apply(block, n_car, ...)`;
     M. the meta gains `de_dnmag` ONLY when the knob is set.
  2. src/vcc2026/dn_mag.py = experiments/631-nmae-headroom/src_dn_mag631.py (installed; refused if a different file is there).
Knob unset: no code runs and nothing is written (KA0: unset rebuilds the base byte for byte).

    uv run python experiments/631-nmae-headroom/patch_builder631.py --root <repo> --check
    uv run python experiments/631-nmae-headroom/patch_builder631.py --root <repo>
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
MODULE_DST = "src/vcc2026/dn_mag.py"
MODULE_SRC = HERE / "src_dn_mag631.py"
MODULE_SHA = hashlib.sha256(MODULE_SRC.read_bytes()).hexdigest()
MARK = "631: DOWN magnitude"
MARK630 = "630: rank-neutral UP"

EDITS = [
    # K. the knob, right after 630's
    ("    DE_UPRN = DE_UPRN_ENV\n",
     "    DE_UPRN = DE_UPRN_ENV\n"
     "# `VCC_DE_DNMAG=d_m<m*100>[_q<n>]` -- 631: DOWN magnitude (experiments/631-nmae-headroom), after zs and 630's UP pass.\n"
     "# Called DOWN genes are delivered at about half the right magnitude. Each called DOWN gene's counts move from its top\n"
     "# carrier to the ballast row whose rise crosses the fewest controls, toward lfc = m x today's; the gene's net rank-sum\n"
     "# change stays within |dAUC| <= 1e-4; the carrier's library is paid back with sub-gate counts from that ballast row, so\n"
     "# every per-gene total and every library is exact; a call-flip guard re-runs without any flipped gene; asserts fail\n"
     "# SAFE (the block comes back unchanged, dn_invariant_failed = 1). Module src/vcc2026/dn_mag.py, sha pinned. Needs\n"
     "# VCC_DE_MAG. Unset is a no-op.\n"
     "DE_DNMAG_ENV = os.environ.get(\"VCC_DE_DNMAG\", \"\").strip()\n"
     "DE_DNMAG = None\n"
     "DE_DNMAG_MODULE_PATH = ROOT / \"src/vcc2026/dn_mag.py\"\n"
     f"DE_DNMAG_MODULE_SHA256 = \"{MODULE_SHA}\"\n"
     "if DE_DNMAG_ENV:\n"
     "    import re as _re_dn  # noqa: PLC0415\n"
     "    _dn = _re_dn.fullmatch(r\"d_m(\\d{3})(?:_q([1-9]\\d*))?\", DE_DNMAG_ENV)\n"
     "    if _dn is None or not 125 <= int(_dn.group(1)) <= 400:\n"
     "        raise SystemExit(f\"VCC_DE_DNMAG must look like 'd_m<m*100>[_q<quantum>]' with m in [1.25, 4], e.g. d_m200: {DE_DNMAG_ENV!r}\")\n"
     "    if DE_MAG is None:\n"
     "        raise SystemExit(\"VCC_DE_DNMAG needs VCC_DE_MAG (it runs after the zs pass, on its evaluator-order control block)\")\n"
     "    DE_DNMAG = DE_DNMAG_ENV\n"),
    # C. per context: sha pin, tie sums, kept-axis map (independent of 630's setup)
    ("                    memory_checkpoint(f\"de_mag scorer control block {context}\", max_rss_gib)\n",
     "                    memory_checkpoint(f\"de_mag scorer control block {context}\", max_rss_gib)\n"
     "                    if DE_DNMAG is not None:\n"
     "                        _dn_sha = hashlib.sha256(DE_DNMAG_MODULE_PATH.read_bytes()).hexdigest()\n"
     "                        if _dn_sha != DE_DNMAG_MODULE_SHA256:\n"
     "                            raise SystemExit(f\"dn_mag.py sha256 {_dn_sha} != pinned {DE_DNMAG_MODULE_SHA256}\")\n"
     "                        from vcc2026 import dn_mag as _dn_mag  # noqa: PLC0415\n"
     "                        rc_dnmag_tie = _dn_mag.scorer_tie(rc_scorer_control)\n"
     "                        rc_dnmag_kpos = np.cumsum(np.asarray(rc_scorer_control.keep, dtype=bool)) - 1\n"
     "                        rc_dnmag_gpos = {str(g): i for i, g in enumerate(spec.genes)}\n"
     "                        print(f\"[de_dnmag] {context}: {DE_DNMAG}\", flush=True)\n"),
    # B. per block, right after 630's UP pass (its block ends with this line; the DOWN pass is at the zs level of indent)
    ("                                rc_diag = {**rc_diag, **_ur_diag}\n",
     "                                rc_diag = {**rc_diag, **_ur_diag}\n"
     "                            if DE_DNMAG is not None:\n"
     "                                # 631: DOWN magnitude after zs (and after 630's UP pass when set); see the knob's comment\n"
     "                                from vcc2026 import dn_mag as _dn_mag  # noqa: PLC0415\n"
     "                                _dn_full = rc_dnmag_gpos.get(str(target), -1)\n"
     "                                _dn_own = (int(rc_dnmag_kpos[_dn_full]) if _dn_full >= 0 and bool(np.asarray(\n"
     "                                    rc_scorer_control.keep, dtype=bool)[_dn_full]) else -1)\n"
     "                                block, _dn_diag = _dn_mag.apply(np.asarray(block), int(spec.cells_per_pert) - int(ballast_k),\n"
     "                                                                rc_scorer_control, DE_DNMAG, tie=rc_dnmag_tie,\n"
     "                                                                own_col=_dn_own)\n"
     "                                rc_diag = {**rc_diag, **_dn_diag}\n"),
    # M. the meta record, after 630's
    ("                                                               \"|dAUC| <= 1e-4 per gene; 0 call flips\"}}),\n",
     "                                                               \"|dAUC| <= 1e-4 per gene; 0 call flips\"}}),\n"
     "            **({} if DE_DNMAG is None else {\"de_dnmag\": {\"variant\": DE_DNMAG, \"module_sha256\": DE_DNMAG_MODULE_SHA256,\n"
     "                                                        \"where\": \"after zs and 630's UP pass; called DOWN genes' \"\n"
     "                                                                 \"counts top carrier -> ballast; totals, libraries \"\n"
     "                                                                 \"exact; |dAUC| <= 1e-4 per gene; 0 call flips\"}}),\n"),
]


def patched_text(text: str) -> str:
    if MARK in text:
        raise SystemExit(f"{BUILDER}: already carries {MARK!r} (a second apply is refused)")
    if MARK630 not in text:
        raise SystemExit(f"{BUILDER}: lacks 630's mark {MARK630!r}: apply 630 first (order ... -> 612 -> 630 -> 631)")
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
        print(f"631 check OK: {BUILDER} carries 630's mark, all {len(EDITS)} anchors found once, the patched text compiles; "
              f"module sha {MODULE_SHA}")
        return 0
    fd, tmp = tempfile.mkstemp(dir=p.parent, prefix=".patch631-")
    with os.fdopen(fd, "w") as f:
        f.write(new)
    os.chmod(tmp, os.stat(p).st_mode)
    os.replace(tmp, p)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_bytes(MODULE_SRC.read_bytes())
    print(f"631 applied to {p}; installed {dst} (sha {MODULE_SHA})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
