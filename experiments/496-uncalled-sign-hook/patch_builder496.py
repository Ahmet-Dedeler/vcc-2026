"""496: the builder side of the uncalled-gene sign hook. `VCC_DE_MDE_SIGN="clean4:<tau>"` gives the emitter's MDE add on
complement genes (not head, not block) the CLEAN4 per-gene prior sign where |P_down - 0.5| >= tau (`mde_add_sign`, see
patch_emitter496.py), and re-pins the emitter to the patched file's sha. Unset = no code path runs and no config key is
written; the only always-on change is the emitter pin (config `de_sign_emitter_sha256`), so KA-1 (unset == 446bj2) is
an OUTPUT-sha comparison, and every later arm's expected config diff vs a pre-496 reference includes that key.
usage: uv run python experiments/496-uncalled-sign-hook/patch_builder496.py <builder.py> <patched rc_emit_k.py> [--check]
"""
import hashlib
import sys
from pathlib import Path

path, emitter = Path(sys.argv[1]), Path(sys.argv[2]); src = path.read_text()
em_sha = hashlib.sha256(emitter.read_bytes()).hexdigest()
if "mde_add_sign" not in emitter.read_text():
    raise SystemExit(f"{emitter} is not the 496-patched emitter")
edits = [
    # 1. the emitter pin
    ('DE_SIGN_EMITTER_SHA256 = "edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7"   # 426h: + optional mde_coarse (72be7e16... before)\n',
     f'DE_SIGN_EMITTER_SHA256 = "{em_sha}"   # 496: + optional mde_add_sign (edb5b8dc... before; 426h: + mde_coarse)\n'),
    # 2. the knob, after VCC_DE_MDE's parse
    ('    if DE_SIGN_PRIOR is None:\n        raise SystemExit("VCC_DE_MDE needs VCC_DE_SIGN_PRIOR (the K-class prior-mode emitter)")\n',
     '    if DE_SIGN_PRIOR is None:\n        raise SystemExit("VCC_DE_MDE needs VCC_DE_SIGN_PRIOR (the K-class prior-mode emitter)")\n'
     '# `VCC_DE_MDE_SIGN="clean4:<tau>"` -- 496 (2026-09-27; 478\'s route): the MDE add on complement genes (not head, not block:\n'
     '# the uncalled genes) takes the CLEAN4 per-gene prior sign where |P_down - 0.5| >= tau (the emitter\'s `mde_add_sign`); the\n'
     '# base LFC, the calls and the bulk are unchanged. Needs VCC_DE_MDE. Unset = no code path runs and no config key is written.\n'
     'DE_MDE_SIGN_ENV = os.environ.get("VCC_DE_MDE_SIGN", "").strip()\n'
     'DE_MDE_SIGN = None\n'
     'if DE_MDE_SIGN_ENV:\n'
     '    _ms_src, _, _ms_tau = DE_MDE_SIGN_ENV.partition(":")\n'
     '    try:\n'
     '        _ms_tau_f = float(_ms_tau)\n'
     '    except ValueError:\n'
     '        raise SystemExit(f"VCC_DE_MDE_SIGN must be \'clean4:<tau>\': {DE_MDE_SIGN_ENV!r}")\n'
     '    if _ms_src.strip().lower() != "clean4" or not (0.0 <= _ms_tau_f < 0.5):\n'
     '        raise SystemExit(f"VCC_DE_MDE_SIGN must be \'clean4:<tau>\' with 0 <= tau < 0.5: {DE_MDE_SIGN_ENV!r}")\n'
     '    if DE_MDE is None:\n'
     '        raise SystemExit("VCC_DE_MDE_SIGN needs VCC_DE_MDE (it sets the direction of the MDE add)")\n'
     '    if DE_SIGN_PRIOR not in ("clean4", "clean4q80"):\n'
     '        raise SystemExit("VCC_DE_MDE_SIGN reads the CLEAN4 prior file: it needs VCC_DE_SIGN_PRIOR=clean4 or clean4q80")\n'
     '    DE_MDE_SIGN = {"source": "clean4", "tau": _ms_tau_f}\n'),
    # 3. the per-context init
    ('            de_mde_add = None\n            if DE_SIGN_PRIOR is not None:\n',
     '            de_mde_add = None\n            de_mde_sign = None\n            if DE_SIGN_PRIOR is not None:\n'),
    # 4. the per-context sign vector on the emitter's gated axis, next to the add
    ('                    print(f"[de_mde] {context}: c {DE_MDE[0]} cap {DE_MDE[1]} min_count {DE_MDE[2]}; add median "\n'
     '                          f"{float(np.median(de_mde_add)):.4f}, {int((de_mde_add > 0).sum())}/{de_mde_add.size} genes", flush=True)\n',
     '                    print(f"[de_mde] {context}: c {DE_MDE[0]} cap {DE_MDE[1]} min_count {DE_MDE[2]}; add median "\n'
     '                          f"{float(np.median(de_mde_add)):.4f}, {int((de_mde_add > 0).sum())}/{de_mde_add.size} genes", flush=True)\n'
     '                    if DE_MDE_SIGN is not None:\n'
     '                        import hashlib as _hashlib  # noqa: PLC0415\n'
     '                        if _hashlib.sha256(DE_SIGN_PRIOR_PATH.read_bytes()).hexdigest() != DE_SIGN_PRIOR_SHA256:\n'
     '                            raise SystemExit("VCC_DE_MDE_SIGN: the CLEAN4 prior file does not match its pin")\n'
     '                        with np.load(DE_SIGN_PRIOR_PATH, allow_pickle=False) as _zp:\n'
     '                            if list(map(str, _zp["genes"])) != list(map(str, spec.genes)):\n'
     '                                raise SystemExit("VCC_DE_MDE_SIGN: the prior gene axis differs from the official order")\n'
     '                            _pdn = np.asarray(_zp["P_down"], dtype=np.float64)[rc_control.keep]\n'
     '                        de_mde_sign = np.ascontiguousarray(np.where(\n'
     '                            np.isfinite(_pdn) & (np.abs(_pdn - 0.5) >= DE_MDE_SIGN["tau"]), np.where(_pdn > 0.5, -1.0, 1.0), 0.0))\n'
     '                        print(f"[de_mde_sign] {context}: clean4 tau {DE_MDE_SIGN[\'tau\']}; {int((de_mde_sign != 0).sum())}/"\n'
     '                              f"{de_mde_sign.size} gated genes carry a sign ({int((de_mde_sign < 0).sum())} down)", flush=True)\n'),
    # 5. into the emitter's lfc_target
    ('                                **({"mde_coarse": True} if (de_mde_add is not None and DE_MDE_COARSE) else {}),\n',
     '                                **({"mde_coarse": True} if (de_mde_add is not None and DE_MDE_COARSE) else {}),\n'
     '                                **({} if de_mde_sign is None else {"mde_add_sign": de_mde_sign}),\n'),
    # 6. the config key
    ('                                                            **({"coarse": True} if DE_MDE_COARSE else {})}}),\n',
     '                                                            **({"coarse": True} if DE_MDE_COARSE else {})}}),\n'
     '                    **({} if DE_MDE_SIGN is None else {"de_mde_sign": dict(DE_MDE_SIGN)}),\n'),
]
for old, new in edits:
    if src.count(old) != 1:
        raise SystemExit(f"anchor count {src.count(old)} != 1: {old[:80]!r}")
if "VCC_DE_MDE_SIGN" in src:
    raise SystemExit("already patched")
for old, new in edits:
    src = src.replace(old, new)
compile(src, str(path), "exec")
if "--check" in sys.argv:
    print(f"patch applies cleanly (not written); emitter pin -> {em_sha[:16]}")
else:
    path.write_text(src)
    print(f"patched {path}; emitter pin -> {em_sha[:16]}")
