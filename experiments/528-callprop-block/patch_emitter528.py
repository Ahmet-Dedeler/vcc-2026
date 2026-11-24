#!/usr/bin/env python3
"""528 — the emitter patch: `rc_hybrid_block_k(..., block_key=None, block_keep=0)`, the CallSet's CPM block chosen by a
call-propensity key instead of control CPM. NOT applied to production; it writes a patched COPY.

What it changes, and only when `block_key` is given (unset = the verbatim code path, bit-identical by construction):
  * the prescribed call set: head = top n_head by |key| (unchanged); block = the top `block_keep` genes by control CPM
    (ref_lin) among the non-head genes, then the top n_cpm - block_keep genes by `block_key` among the rest (stable
    argsort: ties toward the lower gene index, like the CPM key). block_keep = n_cpm reproduces the shipped block for any
    key; block_keep = 0 is the full swap.
  * the prior branch's steered block `blk` (the CLEAN4-shrunk block signs, q80 protection, the MDE's block/complement split
    and 496's complement hook) is the same block by default (`block_steer="key"`), so the block that is called and the
    block that is steered stay one object, as they are today; `block_steer="cpm"` keeps the steered block on control CPM
    and moves only the call set; `block_steer="kept"` steers only the called block's CPM part.
Nothing else moves: pass 1, the |z| multiset, the complement cap, the order law, the class split, placement.

    uv run python experiments/528-callprop-block/patch_emitter528.py [--src src/vcc2026/rc_emit_k.py] [--out <path>]
The default output is this directory's k_emit528.py (the harness's emitter). The source sha is asserted.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SRC_SHA = "1658ae2c400f75ee1659a8ffcaebca569bf3c8afd55ffd1c3006414344ef2791"   # production rc_emit_k.py (496's hook in)

HELPER = '''

def block_from_key(head, ref_lin, n_cpm, block_key, block_keep=0):
    """528: the CallSet's CPM block chosen by a call-propensity key (gated gene axis).

    The top `block_keep` genes by control CPM (`ref_lin`) among the non-head genes are kept, and the remaining
    `n_cpm - block_keep` slots go to the top genes by `block_key` among the rest. Both rankings use a stable argsort on
    the negated key, so ties go to the lower gene index exactly as in the shipped CPM ranking. `block_keep = n_cpm`
    reproduces the shipped block for any key; `block_keep = 0` is the full swap.
    """
    head = np.asarray(head, dtype=bool)
    bk = np.asarray(block_key, dtype=np.float64)
    if bk.shape != head.shape or not np.all(np.isfinite(bk)):
        raise ValueError("block_key must be a finite array over the gated gene axis")
    n_cpm = int(n_cpm)
    nk = int(block_keep)
    if not 0 <= nk <= n_cpm:
        raise ValueError("block_keep must lie in [0, n_cpm]")
    blk = np.zeros(head.size, dtype=bool)
    if nk:
        blk[np.argsort(-np.where(head, -np.inf, np.asarray(ref_lin, dtype=np.float64)), kind="stable")[:nk]] = True
    if n_cpm - nk:
        blk[np.argsort(-np.where(head | blk, -np.inf, bk), kind="stable")[: n_cpm - nk]] = True
    if int(blk.sum()) != n_cpm or bool((blk & head).any()):
        raise AssertionError("the propensity block did not come out at its size")
    return blk
'''

REPLACEMENTS = [
    # 1. the helper, right after prescribed_mask's definition block (before build_tables)
    ("\n\ndef build_tables(cs, m, Mk, Lk):", HELPER + "\n\ndef build_tables(cs, m, Mk, Lk):"),
    # 2. prescribed_mask: optional block key (only 394's harness ever called it)
    ('''def prescribed_mask(key, control, call_set):
    """The CallSet membership recipe of `rc_hybrid_block` as a gated bool mask."""
    key_gated = np.asarray(key, dtype=np.float64)[control.keep]
    key_abs = np.abs(key_gated)
    n_keep = key_abs.size
    head = np.zeros(n_keep, dtype=bool)
    head[np.argsort(-key_abs, kind="stable")[: int(call_set.n_head)]] = True
    cpm_rank = np.where(head, -np.inf,
                        np.asarray(control.ref_lin, dtype=np.float64))
    called = head.copy()
    called[np.argsort(-cpm_rank, kind="stable")[: int(call_set.n_cpm)]] = True
''',
     '''def prescribed_mask(key, control, call_set, block_key=None, block_keep=0):
    """The CallSet membership recipe of `rc_hybrid_block` as a gated bool mask (528: optional propensity block)."""
    key_gated = np.asarray(key, dtype=np.float64)[control.keep]
    key_abs = np.abs(key_gated)
    n_keep = key_abs.size
    head = np.zeros(n_keep, dtype=bool)
    head[np.argsort(-key_abs, kind="stable")[: int(call_set.n_head)]] = True
    if block_key is not None:
        return head | block_from_key(head, control.ref_lin, int(call_set.n_cpm), block_key, block_keep)
    cpm_rank = np.where(head, -np.inf,
                        np.asarray(control.ref_lin, dtype=np.float64))
    called = head.copy()
    called[np.argsort(-cpm_rank, kind="stable")[: int(call_set.n_cpm)]] = True
'''),
    # 3. the signature
    ('''    emit_counts=None,
    order_key=None,
):
    """K-class-library variant of `rc_emit.rc_hybrid_block`.''',
     '''    emit_counts=None,
    order_key=None,
    block_key=None,
    block_keep=0,
    block_steer="key",
):
    """K-class-library variant of `rc_emit.rc_hybrid_block`.'''),
    # 4. docstring note
    ('''    Invariants, asserted per block: (1) per-gene totals bit-identical to''',
     '''    `block_key` (528, optional, gated gene axis): the CallSet's CPM block is chosen by this call-propensity key (see
    `block_from_key`; `block_keep` CPM genes kept first). None (the default) is the shipped CPM block, verbatim.
    `block_steer` (528): "key" (default) = the steered block (`lfc_target["prior"]`: CLEAN4 signs, q80, kf, the MDE's
    block/complement split and 496's complement hook) is the called block, as today; "cpm" = only the CALL set follows the
    key and the steered block stays the CPM block (genes swapped in are called with their unsteered design lfc + MDE,
    genes swapped out keep their CLEAN4 steering but are not called: H1 shows ~30 % of them cannot be held below BH);
    "kept" = the steered block is the called block's CPM part only (genes swapped in are called unsteered, genes swapped
    out are neither called nor steered).

    Invariants, asserted per block: (1) per-gene totals bit-identical to'''),
    # 5. the prescribed call set in rc_hybrid_block_k
    ('''        head = np.zeros(n_keep, dtype=bool)
        head[np.argsort(-key_abs, kind="stable")[: int(call_set.n_head)]] = True
        cpm_rank = np.where(head, -np.inf, np.asarray(control.ref_lin, dtype=np.float64))
        called = head.copy()
        called[np.argsort(-cpm_rank, kind="stable")[: int(call_set.n_cpm)]] = True
        if int(called.sum()) != call_set.size:
            raise AssertionError("the prescribed call set did not come out at its size")
    block_bonus = (''',
     '''        head = np.zeros(n_keep, dtype=bool)
        head[np.argsort(-key_abs, kind="stable")[: int(call_set.n_head)]] = True
        if block_key is None:
            cpm_rank = np.where(head, -np.inf, np.asarray(control.ref_lin, dtype=np.float64))
            called = head.copy()
            called[np.argsort(-cpm_rank, kind="stable")[: int(call_set.n_cpm)]] = True
        else:
            if block_steer not in ("key", "cpm", "kept"):
                raise ValueError("block_steer must be 'key', 'cpm' or 'kept'")
            called = head | block_from_key(head, control.ref_lin, int(call_set.n_cpm), block_key, block_keep)
        if int(called.sum()) != call_set.size:
            raise AssertionError("the prescribed call set did not come out at its size")
    elif block_key is not None:
        raise ValueError("block_key needs a prescribed call set")
    block_bonus = ('''),
    # 6. the steered block in the prior branch
    ('''            blk = np.zeros(n_keep, dtype=bool)
            blk[np.argsort(-np.where(head, -np.inf, ref_lin), kind="stable")[: int(call_set.n_cpm)]] = True
            blk &= (Ms > 0) & np.isfinite(E0v)''',
     '''            blk = np.zeros(n_keep, dtype=bool)
            if block_key is None or block_steer in ("cpm", "kept"):
                blk[np.argsort(-np.where(head, -np.inf, ref_lin), kind="stable")[: int(call_set.n_cpm)]] = True
            if block_key is not None and block_steer == "key":
                blk |= block_from_key(head, ref_lin, int(call_set.n_cpm), block_key, block_keep)
            elif block_key is not None and block_steer == "kept":
                blk &= block_from_key(head, ref_lin, int(call_set.n_cpm), block_key, block_keep)
            blk &= (Ms > 0) & np.isfinite(E0v)'''),
]


def patch_text(src: str) -> str:
    out = src
    for old, new in REPLACEMENTS:
        n = out.count(old)
        if n != 1:
            raise SystemExit(f"patch anchor found {n} times (want 1): {old[:80]!r}")
        out = out.replace(old, new)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=str(ROOT / "src/vcc2026/rc_emit_k.py"))
    ap.add_argument("--out", default=str(HERE / "k_emit528.py"))
    ap.add_argument("--any-sha", action="store_true", help="skip the source sha assert (tests on copies)")
    a = ap.parse_args()
    src = Path(a.src).read_bytes()
    sha = hashlib.sha256(src).hexdigest()
    if not a.any_sha and sha != SRC_SHA:
        raise SystemExit(f"source sha {sha} != the pinned production emitter {SRC_SHA}")
    text = patch_text(src.decode())
    Path(a.out).write_text(text)
    print(f"patched {a.src} ({sha[:8]}) -> {a.out} ({hashlib.sha256(text.encode()).hexdigest()[:8]})")


if __name__ == "__main__":
    main()
