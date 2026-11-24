#!/usr/bin/env python3
"""672 -- `VCC_BULK_LOAD_CTX=<path.npz>:<sha256>`: a bulk matrix PER CONTEXT on the ballast render, as an UNAPPLIED builder patch.

Why (coordinator [400b6e], 2026-10-14): `VCC_BULK_LOAD` (417) swaps the ballast bulk for ONE shared (300 x 18533) effect matrix,
and contexts differ only by `VCC_BULK_CTX_SCALE`'s scalar. A source that predicts each cell line separately (Arc's PIE, exp 700:
per (context, perturbation, gene) lfc from the context's own control mean) cannot be delivered that way. The writer already
renders the ballast per context from `bulk_rows` (the same-line blend, 451, overrides it per context), so a per-context bulk is a
one-line override there plus a loader.

What it changes in scripts/make_submission_pm_cd4.py:
  K. the knob, after VCC_BULK_LOAD's parse: `VCC_BULK_LOAD_CTX="<path.npz>:<sha256>"`; needs VCC_BULK_LOAD; refused with
     VCC_BULK_DUMP. The npz holds `perts`, `genes` (bundle order), `contexts`, and one float64 (300 x 18533) `bulk_<ctx>` per
     bundle context, in the same effect space as the 417 dump: log1p(TS * p_pert) - log1p(TS * p_ctrl).
  L. load + checks, right after `effects = loaded` in VCC_BULK_LOAD's block: sha, perts/genes == the bundle's, contexts == the
     bundle's, shapes, finiteness. The SHARED loaded matrix stays `effects` (carriers' design inputs, tercile keys, every audit),
     so only the ballast render changes.
  W. in the writer's context loop, right after `bulk_rows = bulk_effects`: `bulk_rows = BULK_LOAD_CTX_ROWS[context]`
     (with VCC_SAMELINE, a matched context's same-line blend of the shared matrix replaces it: same-line wins). VCC_BULK_CTX_SCALE and the sharpen q still apply on top.
  M. the meta config gains `bulk_load_ctx` + `bulk_load_ctx_sha256` ONLY when the knob is set.
Knob unset: no code runs and nothing is written (KA0: the all-unset rebuild == 446bj2z ec3a5ae6, and any VCC_BULK_LOAD arm
rebuilds byte for byte).

KA1 (with the knob): an npz whose bulk_A/B/C all equal the shared dump must rebuild the shared-dump arm byte for byte except the
two new config keys (the override is then the identity).

    uv run python experiments/672-bulk-load-ctx/patch_builder672.py --root <repo> --check
    uv run python experiments/672-bulk-load-ctx/patch_builder672.py --root <repo>
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path

BUILDER = "scripts/make_submission_pm_cd4.py"
MARK = "672: per-context bulk load"

EDITS = [
    # K. the knob
    ("if BULK_DUMP is not None and BULK_LOAD is not None:\n"
     "    raise SystemExit(\"VCC_BULK_DUMP and VCC_BULK_LOAD are exclusive\")\n",
     "if BULK_DUMP is not None and BULK_LOAD is not None:\n"
     "    raise SystemExit(\"VCC_BULK_DUMP and VCC_BULK_LOAD are exclusive\")\n"
     "# `VCC_BULK_LOAD_CTX=\"<path.npz>:<sha256>\"` -- 672: per-context bulk load (experiments/672-bulk-load-ctx). One bulk matrix\n"
     "# PER CONTEXT (npz keys bulk_<ctx>, perts, genes, contexts; the 417 dump's effect space) replaces the shared loaded bulk on\n"
     "# that context's ballast render only; the carriers, the tercile keys and every audit keep VCC_BULK_LOAD's shared matrix.\n"
     "# Needs VCC_BULK_LOAD. With VCC_SAMELINE, same-line wins on matched contexts. Unset is a no-op.\n"
     "BULK_LOAD_CTX = os.environ.get(\"VCC_BULK_LOAD_CTX\", \"\").strip() or None\n"
     "BULK_LOAD_CTX_ROWS: dict[str, np.ndarray] = {}\n"
     "if BULK_LOAD_CTX is not None:\n"
     "    if \":\" not in BULK_LOAD_CTX:\n"
     "        raise SystemExit(f\"VCC_BULK_LOAD_CTX must be '<path.npz>:<sha256>': {BULK_LOAD_CTX!r}\")\n"
     "    _blc_path, _blc_sha = BULK_LOAD_CTX.rsplit(\":\", 1)\n"
     "    BULK_LOAD_CTX = (Path(_blc_path), _blc_sha.strip().lower())\n"
     "    if len(BULK_LOAD_CTX[1]) != 64:\n"
     "        raise SystemExit(f\"VCC_BULK_LOAD_CTX sha must be a full sha256: {BULK_LOAD_CTX[1]!r}\")\n"
     "    if BULK_LOAD is None:\n"
     "        raise SystemExit(\"VCC_BULK_LOAD_CTX needs VCC_BULK_LOAD (the shared matrix stays the design and audit input)\")\n"),
    # L. load + checks
    ("        effects = loaded\n",
     "        effects = loaded\n"
     "        if BULK_LOAD_CTX is not None:\n"
     "            _blc_got = file_sha256(BULK_LOAD_CTX[0])\n"
     "            if _blc_got != BULK_LOAD_CTX[1]:\n"
     "                raise SystemExit(f\"VCC_BULK_LOAD_CTX sha mismatch: {_blc_got} != {BULK_LOAD_CTX[1]}\")\n"
     "            with np.load(BULK_LOAD_CTX[0], allow_pickle=False) as _zc:\n"
     "                if list(map(str, _zc[\"perts\"])) != list(map(str, spec.perts)) or list(map(str, _zc[\"genes\"])) != list(map(str, spec.genes)):\n"
     "                    raise SystemExit(\"VCC_BULK_LOAD_CTX perts/genes differ from the bundle's\")\n"
     "                if sorted(map(str, _zc[\"contexts\"])) != sorted(map(str, spec.contexts)):\n"
     "                    raise SystemExit(f\"VCC_BULK_LOAD_CTX contexts {sorted(map(str, _zc['contexts']))} != bundle {sorted(map(str, spec.contexts))}\")\n"
     "                for _c in map(str, spec.contexts):\n"
     "                    _m = np.asarray(_zc[f\"bulk_{_c}\"], dtype=np.float64)\n"
     "                    if _m.shape != loaded.shape or not np.isfinite(_m).all():\n"
     "                        raise SystemExit(f\"VCC_BULK_LOAD_CTX bulk_{_c} shape/finiteness: {_m.shape}\")\n"
     "                    BULK_LOAD_CTX_ROWS[_c] = _m\n"
     "            effect_audit[\"bulk_load_ctx\"] = {\"path\": str(BULK_LOAD_CTX[0]), \"sha256\": _blc_got,\n"
     "                                             \"max_abs_diff_vs_shared\": {c: float(np.max(np.abs(m - loaded))) for c, m in BULK_LOAD_CTX_ROWS.items()}}\n"
     "            print(f\"[bulk-load-ctx] 672: per-context bulk load {BULK_LOAD_CTX[0]} sha={_blc_got[:16]} \"\n"
     "                  f\"max|ctx-shared|={effect_audit['bulk_load_ctx']['max_abs_diff_vs_shared']}\", flush=True)\n"),
    # W. the writer's per-context override
    ("            bulk_rows = bulk_effects\n"
     "            sl_prior: dict[str, np.ndarray] = {}\n",
     "            bulk_rows = bulk_effects\n"
     "            if BULK_LOAD_CTX_ROWS:\n"
     "                # 672: with VCC_SAMELINE on, a context the identity call matched to a same-line bank gets the same-line\n"
     "                # blend of the SHARED matrix below (same-line wins there); unmatched contexts keep this per-context row.\n"
     "                bulk_rows = BULK_LOAD_CTX_ROWS[str(context)]\n"
     "                print(f\"[bulk-load-ctx] {context}: per-context bulk row (a same-line match below replaces it)\", flush=True)\n"
     "            sl_prior: dict[str, np.ndarray] = {}\n"),
    # M. the meta key, only when set
    ("            **({} if BULK_LOAD is None else {\"bulk_load\": str(BULK_LOAD[0].name), \"bulk_load_sha256\": BULK_LOAD[1]}),\n",
     "            **({} if BULK_LOAD is None else {\"bulk_load\": str(BULK_LOAD[0].name), \"bulk_load_sha256\": BULK_LOAD[1]}),\n"
     "            **({} if BULK_LOAD_CTX is None else {\"bulk_load_ctx\": str(BULK_LOAD_CTX[0].name),\n"
     "                                                 \"bulk_load_ctx_sha256\": BULK_LOAD_CTX[1]}),\n"),
]


def patched_text(text: str) -> str:
    if MARK in text:
        raise SystemExit(f"{BUILDER}: already carries {MARK!r} (a second apply is refused)")
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
    ap.add_argument("--check", action="store_true", help="verify anchors; write nothing")
    ap.add_argument("--out", type=Path, help="write the patched builder here instead of in place (for tests)")
    a = ap.parse_args()
    p = a.root / BUILDER
    new = patched_text(p.read_text())
    compile(new, str(p), "exec")
    if a.check:
        print(f"672 check OK: all {len(EDITS)} anchors found once in {BUILDER}; the patched text compiles")
        return 0
    if a.out is not None:
        a.out.write_text(new)
        print(f"672 patched text written to {a.out} (live builder untouched)")
        return 0
    fd, tmp = tempfile.mkstemp(dir=p.parent, prefix=".patch672-")
    with os.fdopen(fd, "w") as f:
        f.write(new)
    os.chmod(tmp, os.stat(p).st_mode)
    os.replace(tmp, p)
    print(f"672 applied to {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
