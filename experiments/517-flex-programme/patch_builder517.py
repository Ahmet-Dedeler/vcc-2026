"""517: the uncalled-gene sign hook (496) sourced from the context's PROGRAMME (Arc's H1 programme, 458d/461), gated like the
block prior. NOT APPLIED: a patch file for the coordinator, to go on top of 496's builder patch (builder 24b95dec... =
HEAD b4975024... + patch_builder496.py; emitter 1658ae2c... = edb5b8dc... + patch_emitter496.py). The emitter needs nothing
new: its `mde_add_sign` is read per block from lfc_target, so a per-target vector is already supported.

`VCC_DE_MDE_SIGN="programme:<tau>[:<tau_c4>]"` (tau_c4 default 0.05):
  - needs VCC_SAMELINE without noprog (the artifact must carry pdown_<c>; every def_sameline.py / artifact458*.py output
    does) and VCC_DE_MDE, like 496's form;
  - in a context whose block prior the programme replaced, the MDE add on complement genes (not head, not block) takes the
    programme's sign where |pdown - 0.5| >= tau; with VCC_DE_PRIOR_GATE set, the gate's weak targets (the ones that fell
    back to the CLEAN4 block prior) take CLEAN4's complement sign at tau_c4 instead (the hook follows the block gate);
  - a context without the programme (not in the artifact, or programme_ok false) takes CLEAN4's at tau_c4 (496's form);
  - config `de_mde_sign` = {"source": "programme", "tau", "tau_c4", "contexts": per-context audit}.
`VCC_DE_MDE_SIGN="file:<npz>:<sha256>:<tau>[:<tau_c4>[:le<N>]]"`: the same, with pdown_<c> read from ANY artifact in the
same-line format (e.g. 458d's `arch1d_programme_abc.npz`, or the day's def_sameline.py output) WITHOUT switching the block
prior: the hook's own source for a CLEAN4-block base (517: the hook's bar is the design's ~0.53 sign accuracy on uncalled
genes, lower than the block's CLEAN4 bar). With VCC_DE_PRIOR_GATE active (only with VCC_SAMELINE) weak targets take
CLEAN4 at tau_c4; a context absent from the file takes CLEAN4 at tau_c4 (recorded). `:le<N>` gives the hook its OWN
strength gate when there is no block gate (468's measure: count(|gated design row| > 0.2) <= N -> CLEAN4 at tau_c4),
refused together with VCC_DE_PRIOR_GATE (which the hook already follows).
`VCC_DE_MDE_SIGN="clean4:<tau>"` parses exactly as 496's (same config dict, no audit key), so KA-496 (clean4:0.05 rebuilds
496a byte for byte) and KA-1 (unset rebuilds the base byte for byte) are the known answers for this patch.
usage: uv run python experiments/517-flex-programme/patch_builder517.py <builder-with-496.py> [--check]
"""
import sys
from pathlib import Path

path = Path(sys.argv[1]); src = path.read_text()
if "VCC_DE_MDE_SIGN" not in src:
    raise SystemExit(f"{path} lacks 496's VCC_DE_MDE_SIGN: apply patch_builder496.py first")
if "programme:<tau>" in src or "file:<npz>" in src:
    raise SystemExit("already patched")
edits = [
    # 1. the parse: clean4:<tau> exactly as 496; programme:<tau>[:<tau_c4>] new
    ('    _ms_src, _, _ms_tau = DE_MDE_SIGN_ENV.partition(":")\n'
     '    try:\n'
     '        _ms_tau_f = float(_ms_tau)\n'
     '    except ValueError:\n'
     '        raise SystemExit(f"VCC_DE_MDE_SIGN must be \'clean4:<tau>\': {DE_MDE_SIGN_ENV!r}")\n'
     '    if _ms_src.strip().lower() != "clean4" or not (0.0 <= _ms_tau_f < 0.5):\n'
     '        raise SystemExit(f"VCC_DE_MDE_SIGN must be \'clean4:<tau>\' with 0 <= tau < 0.5: {DE_MDE_SIGN_ENV!r}")\n',
     '    # 517: `programme:<tau>[:<tau_c4>]` -- the context\'s same-line/lab PROGRAMME (VCC_SAMELINE\'s pdown_<c>) is the sign\n'
     '    # source where it replaced CLEAN4; the prior gate\'s weak targets and programme-less contexts take CLEAN4 at tau_c4.\n'
     '    _ms_parts = [s.strip() for s in DE_MDE_SIGN_ENV.split(":")]\n'
     '    _ms_src = _ms_parts[0].lower()\n'
     '    try:\n'
     '        _ms_tau_f = float(_ms_parts[1]) if (len(_ms_parts) >= 2 and _ms_src != "file") else float("nan")\n'
     '        _ms_tau_c4 = float(_ms_parts[2]) if (len(_ms_parts) == 3 and _ms_src != "file") else 0.05\n'
     '    except ValueError:\n'
     '        raise SystemExit(f"VCC_DE_MDE_SIGN must be \'clean4:<tau>\' or \'programme:<tau>[:<tau_c4>]\': {DE_MDE_SIGN_ENV!r}")\n'
     '    if _ms_src == "clean4":\n'
     '        if len(_ms_parts) != 2 or not (0.0 <= _ms_tau_f < 0.5):\n'
     '            raise SystemExit(f"VCC_DE_MDE_SIGN must be \'clean4:<tau>\' with 0 <= tau < 0.5: {DE_MDE_SIGN_ENV!r}")\n'
     '    elif _ms_src == "programme":\n'
     '        if len(_ms_parts) not in (2, 3) or not (0.0 <= _ms_tau_f < 0.5) or not (0.0 <= _ms_tau_c4 < 0.5):\n'
     '            raise SystemExit(f"VCC_DE_MDE_SIGN must be \'programme:<tau>[:<tau_c4>]\' with 0 <= taus < 0.5: {DE_MDE_SIGN_ENV!r}")\n'
     '        if SAMELINE is None or not SAMELINE["prog"]:\n'
     '            raise SystemExit("VCC_DE_MDE_SIGN=programme reads the programme\'s pdown: it needs VCC_SAMELINE without noprog")\n'
     '    elif _ms_src == "file":\n'
     '        # file:<npz>:<sha256>:<tau>[:<tau_c4>] -- the taus are parts 3 and 4\n'
     '        import re as _re  # noqa: PLC0415\n'
     '        _ms_le = None\n'
     '        if len(_ms_parts) == 6:\n'
     '            _ms_lm = _re.fullmatch(r"le(\\d{1,2})", _ms_parts[5].lower())\n'
     '            if _ms_lm is None:\n'
     '                raise SystemExit(f"VCC_DE_MDE_SIGN file form: the 6th field must be le<N>: {DE_MDE_SIGN_ENV!r}")\n'
     '            _ms_le = int(_ms_lm.group(1))\n'
     '            if DE_PRIOR_GATE is not None:\n'
     '                raise SystemExit("VCC_DE_MDE_SIGN file form with its own le<N> gate: unset VCC_DE_PRIOR_GATE (the hook follows it)")\n'
     '        try:\n'
     '            _ms_tau_f = float(_ms_parts[3]) if len(_ms_parts) in (4, 5, 6) else float("nan")\n'
     '            _ms_tau_c4 = float(_ms_parts[4]) if len(_ms_parts) in (5, 6) else 0.05\n'
     '        except ValueError:\n'
     '            raise SystemExit(f"VCC_DE_MDE_SIGN must be \'file:<npz>:<sha256>:<tau>[:<tau_c4>[:le<N>]]\': {DE_MDE_SIGN_ENV!r}")\n'
     '        if len(_ms_parts) not in (4, 5, 6) or len(_ms_parts[2]) != 64 or not (0.0 <= _ms_tau_f < 0.5) or not (0.0 <= _ms_tau_c4 < 0.5):\n'
     '            raise SystemExit(f"VCC_DE_MDE_SIGN must be \'file:<npz>:<sha256>:<tau>[:<tau_c4>[:le<N>]]\' with 0 <= taus < 0.5: {DE_MDE_SIGN_ENV!r}")\n'
     '    else:\n'
     '        raise SystemExit(f"VCC_DE_MDE_SIGN source must be clean4, programme or file: {DE_MDE_SIGN_ENV!r}")\n'),
    ('    DE_MDE_SIGN = {"source": "clean4", "tau": _ms_tau_f}\n',
     '    DE_MDE_SIGN = ({"source": "clean4", "tau": _ms_tau_f} if _ms_src == "clean4"\n'
     '                   else {"source": "programme", "tau": _ms_tau_f, "tau_c4": _ms_tau_c4} if _ms_src == "programme"\n'
     '                   else {"source": "file", "path": _ms_parts[1], "sha256": _ms_parts[2], "tau": _ms_tau_f,\n'
     '                         "tau_c4": _ms_tau_c4, **({} if _ms_le is None else {"le": _ms_le})})\n'
     'DE_MDE_SIGN_AUDIT: dict[str, dict[str, object]] = {}   # 517: per-context record of the programme form only\n'),
    # 2. sameline_load: also carry the programme's P_down when the artifact has it (nothing else reads the key)
    ('            rec = {k: np.asarray(z[f"{k}_{c}"]) for k in ("programme", "targets", "rows", "signs")}\n',
     '            rec = {k: np.asarray(z[f"{k}_{c}"]) for k in ("programme", "targets", "rows", "signs")}\n'
     '            if f"pdown_{c}" in z.files:   # 517: the programme\'s P_down, read only by VCC_DE_MDE_SIGN=programme\n'
     '                rec["pdown"] = np.asarray(z[f"pdown_{c}"], dtype=np.float64)\n'),
    # 3. the per-context init
    ('            de_mde_add = None\n            de_mde_sign = None\n            if DE_SIGN_PRIOR is not None:\n',
     '            de_mde_add = None\n            de_mde_sign = None\n'
     '            de_mde_sign_c4 = None                  # 517: CLEAN4\'s complement sign for the gate\'s weak targets\n'
     '            de_mde_sign_weak: set[str] = set()\n'
     '            _sl_pdown = None\n'
     '            if DE_SIGN_PRIOR is not None:\n'),
    # 4. the programme's P_down, where it replaced CLEAN4 in this context
    ('                    if _prog:\n                        de_sign_full = np.asarray(_slr["programme"], dtype=np.float64)\n',
     '                    if _prog:\n                        de_sign_full = np.asarray(_slr["programme"], dtype=np.float64)\n'
     '                        _sl_pdown = _slr.get("pdown")   # 517\n'),
    # 5. the hook follows the prior gate: weak targets that fell back to the CLEAN4 block prior
    ('                        for _t in _pg_weak:\n'
     '                            if _t not in sl_prior:\n'
     '                                sl_prior[_t] = _pg_c4\n',
     '                        for _t in _pg_weak:\n'
     '                            if _t not in sl_prior:\n'
     '                                sl_prior[_t] = _pg_c4\n'
     '                        de_mde_sign_weak = {_t for _t in _pg_weak if _t not in _pg_kept}   # 517\n'),
    # 6. the sign vectors
    ('                            _pdn = np.asarray(_zp["P_down"], dtype=np.float64)[rc_control.keep]\n'
     '                        de_mde_sign = np.ascontiguousarray(np.where(\n'
     '                            np.isfinite(_pdn) & (np.abs(_pdn - 0.5) >= DE_MDE_SIGN["tau"]), np.where(_pdn > 0.5, -1.0, 1.0), 0.0))\n',
     '                            _pdn = np.asarray(_zp["P_down"], dtype=np.float64)[rc_control.keep]\n'
     '                        if DE_MDE_SIGN["source"] == "clean4":\n'
     '                            de_mde_sign = np.ascontiguousarray(np.where(\n'
     '                                np.isfinite(_pdn) & (np.abs(_pdn - 0.5) >= DE_MDE_SIGN["tau"]), np.where(_pdn > 0.5, -1.0, 1.0), 0.0))\n'
     '                        else:\n'
     '                            # 517: the programme where it replaced CLEAN4; CLEAN4 at tau_c4 for weak targets / elsewhere\n'
     '                            _c4s = np.ascontiguousarray(np.where(\n'
     '                                np.isfinite(_pdn) & (np.abs(_pdn - 0.5) >= DE_MDE_SIGN["tau_c4"]), np.where(_pdn > 0.5, -1.0, 1.0), 0.0))\n'
     '                            if DE_MDE_SIGN["source"] == "file":\n'
     '                                _msp = Path(DE_MDE_SIGN["path"]) if Path(DE_MDE_SIGN["path"]).is_absolute() else ROOT / DE_MDE_SIGN["path"]\n'
     '                                if _hashlib.sha256(_msp.read_bytes()).hexdigest() != DE_MDE_SIGN["sha256"]:\n'
     '                                    raise SystemExit(f"VCC_DE_MDE_SIGN=file: {_msp} does not match its sha256")\n'
     '                                with np.load(_msp, allow_pickle=False) as _zf:\n'
     '                                    if list(map(str, _zf["genes"])) != list(map(str, spec.genes)):\n'
     '                                        raise SystemExit("VCC_DE_MDE_SIGN=file: the artifact gene axis differs from the official order")\n'
     '                                    _sl_pdown = (np.asarray(_zf[f"pdown_{context}"], dtype=np.float64)\n'
     '                                                 if f"pdown_{context}" in _zf.files else None)\n'
     '                                if "le" in DE_MDE_SIGN:\n'
     '                                    # the hook\'s own strength gate: 468\'s measure on this context\'s design rows\n'
     '                                    _hs = np.empty(len(spec.perts), dtype=np.float64)\n'
     '                                    for _hs_t in range(len(spec.perts)):\n'
     '                                        _hs_g = emission_rank_axis(context_rows[_hs_t], b_ctrl, min_cpm=min_cpm)[0]\n'
     '                                        _hs_k = np.where(_hs_g, np.asarray(context_rows[_hs_t], dtype=np.float64), 0.0)\n'
     '                                        _hs[_hs_t] = float((np.abs(_hs_k) > 0.2).sum())\n'
     '                                    de_mde_sign_weak = {str(spec.perts[_i]) for _i in np.flatnonzero(_hs <= DE_MDE_SIGN["le"])}\n'
     '                            if _sl_pdown is not None:\n'
     '                                if _sl_pdown.shape != (len(spec.genes),):\n'
     '                                    raise SystemExit(f"VCC_DE_MDE_SIGN=programme: {context}\'s pdown is not on the official gene axis")\n'
     '                                _ppd = np.asarray(_sl_pdown, dtype=np.float64)[rc_control.keep]\n'
     '                                de_mde_sign = np.ascontiguousarray(np.where(\n'
     '                                    np.isfinite(_ppd) & (np.abs(_ppd - 0.5) >= DE_MDE_SIGN["tau"]), np.where(_ppd > 0.5, -1.0, 1.0), 0.0))\n'
     '                                de_mde_sign_c4 = _c4s\n'
     '                                _ms_used = DE_MDE_SIGN["source"]\n'
     '                            else:\n'
     '                                de_mde_sign, _ms_used = _c4s, f"clean4 (no {DE_MDE_SIGN[\'source\']} pdown in this context)"\n'
     '                                de_mde_sign_weak = set()\n'
     '                            DE_MDE_SIGN_AUDIT[str(context)] = {\n'
     '                                "source": _ms_used, "n_sign": int((de_mde_sign != 0).sum()), "n_down": int((de_mde_sign < 0).sum()),\n'
     '                                "weak_targets_clean4": int(len(de_mde_sign_weak)), "n_sign_clean4": int((_c4s != 0).sum())}\n'
     '                            print(f"[de_mde_sign] {context}: {_ms_used} tau {DE_MDE_SIGN[\'tau\']}; "\n'
     '                                  f"{DE_MDE_SIGN_AUDIT[str(context)]}", flush=True)\n'),
    # 7. into the emitter's lfc_target, per target
    ('                                **({} if de_mde_sign is None else {"mde_add_sign": de_mde_sign}),\n',
     '                                **({} if de_mde_sign is None else {"mde_add_sign": (\n'
     '                                    de_mde_sign_c4 if (de_mde_sign_c4 is not None and str(target_value) in de_mde_sign_weak)\n'
     '                                    else de_mde_sign)}),\n'),
    # 8. the config key (the clean4 form writes exactly 496's dict)
    ('                    **({} if DE_MDE_SIGN is None else {"de_mde_sign": dict(DE_MDE_SIGN)}),\n',
     '                    **({} if DE_MDE_SIGN is None else {"de_mde_sign": {\n'
     '                        **dict(DE_MDE_SIGN), **({"contexts": dict(DE_MDE_SIGN_AUDIT)} if DE_MDE_SIGN_AUDIT else {})}}),\n'),
]
for old, new in edits:
    if src.count(old) != 1:
        raise SystemExit(f"anchor count {src.count(old)} != 1: {old[:90]!r}")
for old, new in edits:
    src = src.replace(old, new)
compile(src, str(path), "exec")
if "--check" in sys.argv:
    print("patch applies cleanly (not written)")
else:
    path.write_text(src)
    print(f"patched {path}")
