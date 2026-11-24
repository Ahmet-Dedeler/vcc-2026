#!/usr/bin/env python3
"""612 -- `VCC_DE_MDE_SIGN_KEYGATE=<q>`, the KEY GATE on the uncalled-gene sign hook, as an UNAPPLIED builder patch.

The coordinator applies it AFTER 598 (10-06): 538 -> 554 -> 590 (10-05 bundle) -> 598 -> 612. It refuses a builder without
598's mark, and a second apply.

Why (council 600, [21dd4c]; experiments/612-key-gated-hook/RESULTS.md):
- The hook (`VCC_DE_MDE_SIGN`, 496 / 517) hands the emitter a sign for every complement gene its source covers.
- On h1-test's uncalled gate genes our design key beats CLEAN4 where |key| is high and loses where it is low
  (610/gatesign610).
- On 480's held-out cancer folds the key-GATED rule (the source's sign on each target's weak-key genes, the delivered key
  sign elsewhere) is the min-max-regret sign rule (612/kgh612, qscan612). With identity-augmented signs in the weak half
  (622's form) it beats the CLEAN4 version on 3/4 folds (622/hyb622).
- So the gate is ORTHOGONAL to the source. One knob masks whatever per-target sign the hook delivers: CLEAN4 (`clean4:`),
  a same-line/identity programme (`programme:`), or a file (`file:`, e.g. 604's identity artifact's pdown_<c>).

What it changes (ONE file, scripts/make_submission_pm_cd4.py; the emitter stays at 0f56209d; the VCC_DE_MDE_SIGN parse is
NOT touched):
  D. after the hook's parse: `VCC_DE_MDE_SIGN_KEYGATE=<q>` (0 < q < 1; it needs VCC_DE_MDE_SIGN). There are two helpers.
     `_kg_mask_sign(sign_gated, key_full, keep, q)` = where(|key[keep]| == 0 or |key[keep]| < the q-quantile of the
     target's nonzero |key[keep]|, sign_gated, 0); this is exactly run612.py's gate, which the H1 real-emitter prices used.
     `_kg_apply(sign, key_full, keep)` returns `sign` itself (the same object) when the knob is unset;
  F. the per-target call passes `_kg_apply(<the sign it passed before>, key, rc_control.keep)`. `key` is the design key
     the same call passes as `key=` (np.where(gate, context_rows[t], 0) on the full axis; the emitter reads
     key[control.keep]), and 0 = the emitter keeps its delivered sign;
  M. the meta gains `de_mde_sign_keygate` ONLY when the knob is set.
Knob unset: the same arrays reach the emitter and the config is unchanged (KA-1 / KA-2 byte for byte).

    uv run python experiments/612-key-gated-hook/patch612.py --root <repo> --check   # anchors + order only, writes nothing
    uv run python experiments/612-key-gated-hook/patch612.py --root <repo>           # apply (refuses a second time)
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path

MARK = "612: key gate"
MARK598 = "598: tercilec"
BUILDER = Path("scripts/make_submission_pm_cd4.py")

D_OLD = '''DE_MDE_SIGN_AUDIT: dict[str, dict[str, object]] = {}   # 517: per-context record of the programme form only
'''
D_NEW = '''DE_MDE_SIGN_AUDIT: dict[str, dict[str, object]] = {}   # 517: per-context record of the programme form only
# `VCC_DE_MDE_SIGN_KEYGATE=<q>` -- 612: key gate (council 600, 2026-10-02). It masks the hook's per-target sign (any
# VCC_DE_MDE_SIGN source) to the target's WEAK-key genes: |design key| on the gated axis is zero or below the q-quantile
# of the target's nonzero |key|. The strong half gets 0 (the emitter keeps its delivered sign, i.e. the key's). This is
# the min-max-regret sign rule on 480's held-out cancer folds (612). Unset = the hook's sign passes unchanged.
DE_MDE_SIGN_KEYGATE_ENV = os.environ.get("VCC_DE_MDE_SIGN_KEYGATE", "").strip()
DE_MDE_SIGN_KEYGATE = None
if DE_MDE_SIGN_KEYGATE_ENV:   # 612: key gate
    try:
        DE_MDE_SIGN_KEYGATE = float(DE_MDE_SIGN_KEYGATE_ENV)
    except ValueError:
        raise SystemExit(f"VCC_DE_MDE_SIGN_KEYGATE must be a quantile 0 < q < 1: {DE_MDE_SIGN_KEYGATE_ENV!r}")
    if not (0.0 < DE_MDE_SIGN_KEYGATE < 1.0):
        raise SystemExit(f"VCC_DE_MDE_SIGN_KEYGATE must be a quantile 0 < q < 1: {DE_MDE_SIGN_KEYGATE_ENV!r}")
    if DE_MDE_SIGN is None:
        raise SystemExit("VCC_DE_MDE_SIGN_KEYGATE gates the hook's sign: it needs VCC_DE_MDE_SIGN")


def _kg_mask_sign(sign_gated, key_full, keep, q):
    """612: key gate. The hook's sign on the target's WEAK-key genes, 0 elsewhere (rc_emit_k `mde_add_sign`: 0 = keep the
    delivered sign). Weak = |key[keep]| == 0 or below the q-quantile of the target's nonzero |key[keep]|; genes TIED at the
    cut (capped or rounded keys) fill the weak count up to floor(q * n_nonzero) in stable gene order, so a tied target
    cannot lose its whole weak half (605's cloud review). `key_full` is the design key this target's emitter call
    receives (full gene axis); `keep` is the control's gated axis, on which `sign_gated` and the emitter's key_gated
    live. Without a tie at the cut this is run612.py's gate exactly (test612 T4a/T4b)."""
    k = np.abs(np.asarray(key_full, dtype=np.float64)[np.asarray(keep, dtype=bool)])
    s = np.asarray(sign_gated, dtype=np.float64)
    if k.shape != s.shape:
        raise ValueError(f"612: gated key {k.shape} and hook sign {s.shape} differ")
    nz = k > 0
    if not nz.any():
        return np.ascontiguousarray(s)
    cut = float(np.quantile(k[nz], q))
    weak = (~nz) | (k < cut)
    need = int(np.floor(q * int(nz.sum()) + 1e-9)) - int((nz & (k < cut)).sum())
    if need > 0:
        weak[np.flatnonzero(nz & (k == cut))[:need]] = True
    return np.ascontiguousarray(np.where(weak, s, 0.0))


def _kg_apply(sign_gated, key_full, keep):
    """612: key gate. The hook's per-target sign, masked when VCC_DE_MDE_SIGN_KEYGATE is set; the SAME object otherwise."""
    if DE_MDE_SIGN_KEYGATE is None:
        return sign_gated
    return _kg_mask_sign(sign_gated, key_full, keep, DE_MDE_SIGN_KEYGATE)
'''
F_OLD = '''                                **({} if de_mde_sign is None else {"mde_add_sign": (
                                    de_mde_sign_c4 if (de_mde_sign_c4 is not None and str(target_value) in de_mde_sign_weak)
                                    else de_mde_sign)}),
'''
F_NEW = '''                                **({} if de_mde_sign is None else {"mde_add_sign": _kg_apply(   # 612: key gate (no-op if unset)
                                    de_mde_sign_c4 if (de_mde_sign_c4 is not None and str(target_value) in de_mde_sign_weak)
                                    else de_mde_sign, key, rc_control.keep)}),
'''
M_OLD = '''                    **({} if DE_MDE_SIGN is None else {"de_mde_sign": {
                        **dict(DE_MDE_SIGN), **({"contexts": dict(DE_MDE_SIGN_AUDIT)} if DE_MDE_SIGN_AUDIT else {})}}),
'''
M_NEW = '''                    **({} if DE_MDE_SIGN is None else {"de_mde_sign": {
                        **dict(DE_MDE_SIGN), **({"contexts": dict(DE_MDE_SIGN_AUDIT)} if DE_MDE_SIGN_AUDIT else {})}}),
                    **({} if DE_MDE_SIGN_KEYGATE is None else {"de_mde_sign_keygate": {   # 612: key gate (absent when unset)
                        "q": DE_MDE_SIGN_KEYGATE, "rule": "hook sign kept where |design key| == 0 or < the target's q-quantile "
                                                          "of nonzero |key| on the gated axis; 0 (delivered sign) elsewhere"}}),
'''
EDITS = [("D knob + helpers", D_OLD, D_NEW), ("F per-target call", F_OLD, F_NEW), ("M meta", M_OLD, M_NEW)]


def patched_text(text: str) -> str:
    if MARK in text:
        raise SystemExit(f"{BUILDER}: already carries {MARK!r} (a second apply is refused)")
    if MARK598 not in text:
        raise SystemExit(f"{BUILDER}: lacks 598's mark {MARK598!r}: apply 598 first (order 538 -> 554 -> 590 -> 598 -> 612)")
    for name, old, _ in EDITS:
        n = text.count(old)
        if n != 1:
            raise SystemExit(f"{BUILDER}: anchor {name} found {n} times (expected exactly 1)")
    for name, old, new in EDITS:
        text = text.replace(old, new, 1)
    if MARK not in text:
        raise SystemExit("612: the patched text lacks the mark (patch bug)")
    return text


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--check", action="store_true", help="verify anchors and order; write nothing")
    a = ap.parse_args()
    p = a.root / BUILDER
    text = p.read_text()
    new = patched_text(text)
    compile(new, str(p), "exec")
    if a.check:
        print(f"612 check OK: {BUILDER} carries 598's mark, all {len(EDITS)} anchors found once, the patched text compiles")
        return 0
    fd, tmp = tempfile.mkstemp(dir=p.parent, prefix=".patch612-")
    with os.fdopen(fd, "w") as f:
        f.write(new)
    os.chmod(tmp, os.stat(p).st_mode)
    os.replace(tmp, p)
    print(f"612 applied: {BUILDER} ({len(EDITS)} edits)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
