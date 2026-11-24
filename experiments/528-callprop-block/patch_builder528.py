#!/usr/bin/env python3
"""528: the call-propensity CPM block, `VCC_RC_CALLSET_BLOCK`. NOT APPLIED: a patch for the coordinator, to go on the
builder after 527 (517's hook in) together with patch_emitter528.py on the emitter (production 1658ae2c -> the patched
emitter, whose sha this patch pins in DE_SIGN_EMITTER_SHA256).

`VCC_RC_CALLSET_BLOCK` unset or `cpm`: the shipped block, bit-identical (no new keyword reaches the emitter; the emitter's
unset path is verbatim). `prop:<npz>:<sha256>[:<swap>[:<steer>]]`: per context, the CallSet's CPM block is chosen by the
artifact's `key_<context>` (official gene axis; callblock528.py writes it from the context's OWN controls: log10 control
CPM + beta x excess dropout, median-rescaled to H1's dropout scale) instead of control CPM:
  - `swap` in (0, 1] (default 1) = the share of the CPM block the key may replace: the top n_cpm (1 - swap) genes by control
    CPM are kept, the rest filled by the key (the emitter's block_keep); 1 = the full propensity block;
  - `steer` = kept (default) | key | cpm, the emitter's block_steer: which genes take the block's CLEAN4 steering. `kept`
    steers only the called block's CPM part (genes swapped in are called unsteered, genes swapped out are neither called
    nor steered: H1's best form); `key` steers the whole called block (the swapped-in genes take the kf shrink: H1 nmae
    +0.004); `cpm` keeps steering the CPM block (~30 % of the steered-not-called genes leak through BH on H1);
  - the head (top 400 by |H|), the block size, z_cap, the order law, the MDE and the complement hook are unchanged;
  - needs the K-class path (VCC_DE_SIGN_PRIOR) and the patched emitter (its `block_key` / `block_steer` keywords);
  - a context absent from the artifact is refused (the artifact is per panel: October computes D/E/F's on the day);
  - config `rc_callset_block` = {"source": "prop", "path", "sha256", "swap", "block_keep", "steer", "contexts": audit}.
usage: uv run python experiments/528-callprop-block/patch_builder528.py <builder.py> [--emitter-sha <sha256>]
"""
import argparse
import hashlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument("builder")
ap.add_argument("--emitter-sha", default=None, help="the patched emitter's sha256 (default: this directory's k_emit528.py)")
a = ap.parse_args()
path = Path(a.builder)
src = path.read_text()
if "VCC_RC_CALLSET_BLOCK" in src:
    raise SystemExit("already patched")
em_sha = a.emitter_sha or hashlib.sha256((HERE / "k_emit528.py").read_bytes()).hexdigest()
OLD_EM = 'DE_SIGN_EMITTER_SHA256 = "1658ae2c400f75ee1659a8ffcaebca569bf3c8afd55ffd1c3006414344ef2791"'
edits = [
    # 1. the emitter pin: 528's patched emitter (block_key); the unset path is verbatim
    (OLD_EM + "   # 496: + optional mde_add_sign (edb5b8dc... before; 426h: + mde_coarse)\n",
     f'DE_SIGN_EMITTER_SHA256 = "{em_sha}"   # 528: + optional block_key / block_keep / block_steer (1658ae2c... before; 496: + mde_add_sign; 426h: + mde_coarse)\n'),
    # 2. the parse, after the call set's three numbers (and 523's size override)
    ("RC_CALLSET_Z_CAP = 1.90\n",
     "RC_CALLSET_Z_CAP = 1.90\n"
     "# 528 (2026-09-29): VCC_RC_CALLSET_BLOCK = cpm (default, the shipped CPM block, bit-identical) |\n"
     "# prop:<npz>:<sha256>[:<swap>[:<steer>]]: the CallSet's CPM block chosen per context by a call-propensity key from the\n"
     "# context's own controls (the artifact's key_<context>, official axis; experiments/528-callprop-block/callblock528.py).\n"
     "# swap in (0, 1] = the share of the block the key may replace (the top n_cpm (1 - swap) genes by control CPM are kept);\n"
     "# steer = kept (default: only the called block's CPM part takes the CLEAN4 steering) | key | cpm (the emitter's\n"
     "# block_steer). Needs VCC_DE_SIGN_PRIOR (K-class path).\n"
     'RC_CALLSET_BLOCK_ENV = os.environ.get("VCC_RC_CALLSET_BLOCK", "").strip()\n'
     "RC_CALLSET_BLOCK = None\n"
     "RC_CALLSET_BLOCK_AUDIT: dict[str, object] = {}\n"
     "_RC_CALLSET_BLOCK_CACHE: dict[str, object] = {}\n"
     'if RC_CALLSET_BLOCK_ENV and RC_CALLSET_BLOCK_ENV.lower() != "cpm":\n'
     '    _cb = [s.strip() for s in RC_CALLSET_BLOCK_ENV.split(":")]\n'
     '    if _cb[0].lower() != "prop" or len(_cb) not in (3, 4, 5) or len(_cb[2]) != 64 or not _cb[1]:\n'
     "        raise SystemExit(f\"VCC_RC_CALLSET_BLOCK must be 'cpm' or 'prop:<npz>:<sha256>[:<swap>[:<steer>]]': {RC_CALLSET_BLOCK_ENV!r}\")\n"
     "    try:\n"
     "        _cb_swap = float(_cb[3]) if len(_cb) >= 4 else 1.0\n"
     "    except ValueError:\n"
     "        raise SystemExit(f\"VCC_RC_CALLSET_BLOCK swap must be a number in (0, 1]: {RC_CALLSET_BLOCK_ENV!r}\")\n"
     "    if not (0.0 < _cb_swap <= 1.0):\n"
     "        raise SystemExit(f\"VCC_RC_CALLSET_BLOCK swap must be in (0, 1]: {RC_CALLSET_BLOCK_ENV!r}\")\n"
     '    _cb_steer = _cb[4].lower() if len(_cb) == 5 else "kept"\n'
     '    if _cb_steer not in ("kept", "key", "cpm"):\n'
     "        raise SystemExit(f\"VCC_RC_CALLSET_BLOCK steer must be kept, key or cpm: {RC_CALLSET_BLOCK_ENV!r}\")\n"
     '    RC_CALLSET_BLOCK = {"source": "prop", "path": _cb[1], "sha256": _cb[2].lower(), "swap": _cb_swap,\n'
     '                        "block_keep": int(RC_CALLSET_CPM - int(round(_cb_swap * RC_CALLSET_CPM))), "steer": _cb_steer}\n'),
    # 3. a per-context slot beside rc_control
    ("            rc_control = None\n"
     "            rc_scorer_control = None   # 429: VCC_DE_MAG's evaluator-order control block\n",
     "            rc_control = None\n"
     "            rc_scorer_control = None   # 429: VCC_DE_MAG's evaluator-order control block\n"
     "            rc_block_key = None        # 528: VCC_RC_CALLSET_BLOCK's per-context key on the gated axis\n"),
    # 4. the per-context key, right after the control block is built
    ('                memory_checkpoint(f"rc control block {context}", max_rss_gib)\n',
     '                memory_checkpoint(f"rc control block {context}", max_rss_gib)\n'
     "                if RC_CALLSET_BLOCK is not None:\n"
     "                    # 528: the CPM block chosen by the artifact's call-propensity key for THIS context (its own controls)\n"
     "                    if DE_SIGN_PRIOR is None:\n"
     '                        raise SystemExit("VCC_RC_CALLSET_BLOCK needs the K-class emitter path (VCC_DE_SIGN_PRIOR)")\n'
     "                    import hashlib as _hashlib  # noqa: PLC0415\n"
     "                    import inspect as _inspect  # noqa: PLC0415\n"
     "                    from vcc2026 import rc_emit_k as _rck  # noqa: PLC0415\n"
     '                    if "block_steer" not in _inspect.signature(_rck.rc_hybrid_block_k).parameters:\n'
     '                        raise SystemExit("VCC_RC_CALLSET_BLOCK needs the 528 emitter (rc_hybrid_block_k(block_key=..., block_steer=...))")\n'
     '                    _cbp = Path(RC_CALLSET_BLOCK["path"])\n'
     "                    _cbp = _cbp if _cbp.is_absolute() else ROOT / _cbp\n"
     "                    if not _RC_CALLSET_BLOCK_CACHE:\n"
     '                        if _hashlib.sha256(_cbp.read_bytes()).hexdigest() != RC_CALLSET_BLOCK["sha256"]:\n'
     '                            raise SystemExit(f"VCC_RC_CALLSET_BLOCK: {_cbp} does not match its sha256")\n'
     "                        with np.load(_cbp, allow_pickle=False) as _cbz:\n"
     '                            if list(map(str, _cbz["genes"])) != list(map(str, spec.genes)):\n'
     '                                raise SystemExit("VCC_RC_CALLSET_BLOCK: the artifact gene axis differs from the official order")\n'
     '                            _RC_CALLSET_BLOCK_CACHE.update({k[4:]: np.asarray(_cbz[k], dtype=np.float64)\n'
     '                                                            for k in _cbz.files if k.startswith("key_")})\n'
     "                    if str(context) not in _RC_CALLSET_BLOCK_CACHE:\n"
     '                        raise SystemExit(f"VCC_RC_CALLSET_BLOCK: no key_{context} in {_cbp}")\n'
     "                    rc_block_key = np.ascontiguousarray(_RC_CALLSET_BLOCK_CACHE[str(context)][rc_control.keep])\n"
     "                    if not np.all(np.isfinite(rc_block_key)):\n"
     '                        raise SystemExit(f"VCC_RC_CALLSET_BLOCK: key_{context} is not finite on the 5-CPM universe")\n'
     "                    _cb_cpm = np.argsort(-np.asarray(rc_control.ref_lin, dtype=np.float64), kind=\"stable\")[:RC_CALLSET_CPM]\n"
     "                    _cb_key = np.argsort(-rc_block_key, kind=\"stable\")[:RC_CALLSET_CPM]\n"
     "                    RC_CALLSET_BLOCK_AUDIT[str(context)] = {\n"
     '                        "n_universe": int(rc_control.keep.sum()),\n'
     '                        "swapped_vs_cpm_no_head": int(np.setdiff1d(_cb_key, _cb_cpm).size),\n'
     '                        "key_sha256": _hashlib.sha256(rc_block_key.tobytes()).hexdigest()}\n'
     '                    print(f"[rc_callset_block] {context}: prop swap {RC_CALLSET_BLOCK[\'swap\']} steer {RC_CALLSET_BLOCK[\'steer\']} (block_keep "\n'
     "                          f\"{RC_CALLSET_BLOCK['block_keep']}); {RC_CALLSET_BLOCK_AUDIT[str(context)]}\", flush=True)\n"),
    # 5. the emitter call: the key reaches the emitter only when set (the unset call is unchanged)
    ("                            call_set=rc_emit.CallSet(RC_CALLSET_HEAD, RC_CALLSET_CPM, RC_CALLSET_Z_CAP),\n"
     "                            classes=de_sign_classes,\n",
     "                            call_set=rc_emit.CallSet(RC_CALLSET_HEAD, RC_CALLSET_CPM, RC_CALLSET_Z_CAP),\n"
     "                            **({} if rc_block_key is None else {\"block_key\": rc_block_key,\n"
     "                                                                \"block_keep\": RC_CALLSET_BLOCK[\"block_keep\"],\n"
     "                                                                \"block_steer\": RC_CALLSET_BLOCK[\"steer\"]}),\n"
     "                            classes=de_sign_classes,\n"),
    # 6. config
    ('                    "rc_callset_z_cap": RC_CALLSET_Z_CAP,\n',
     '                    "rc_callset_z_cap": RC_CALLSET_Z_CAP,\n'
     '                    **({} if RC_CALLSET_BLOCK is None else {"rc_callset_block": {\n'
     '                        **dict(RC_CALLSET_BLOCK), "contexts": dict(RC_CALLSET_BLOCK_AUDIT)}}),\n'),
]
for old, new in edits:
    n = src.count(old)
    if n != 1:
        raise SystemExit(f"anchor found {n} times (want 1): {old[:90]!r}")
    src = src.replace(old, new)
path.write_text(src)
print(f"patched {path}: VCC_RC_CALLSET_BLOCK (emitter pin -> {em_sha[:16]})")
