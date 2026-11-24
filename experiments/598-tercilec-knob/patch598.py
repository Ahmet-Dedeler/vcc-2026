#!/usr/bin/env python3
"""598 -- `tercilec:<low>,<mid>,<high>`, 493h's carried-only key tercile, as an UNAPPLIED patch (the coordinator applies
it AFTER the 10-05 bundle: 538 -> 554 -> 590 -> 598).

Why (`docs/october-decisions.json` sharpen_pt, DECIDED 10-02 `tercile:2.5,2.0,1.75`: "UPGRADE QUEUED: 493h's carried-only
cut points"): `tercile` cuts the key ratio r_t = ||bulk_t|| / ||carriers_t|| (non-panel columns; r = c * w_t) at its
1/3 and 2/3 quantiles over ALL targets. Non-carried targets have w = 1 (r = c up to ~1e-14), so on a panel where they are
numerous (efp: 153-156 of ~300 on two folds) they form a tie block that straddles a cut point and empties MID (493g).
493h priced the rule on 438's October-like panels: cut points over carried targets only, non-carried -> MID, floor ->
LOW: efp +0.84 % rel cos (3/3) vs the plain rule's +0.48 % (2/3); on A/B/C it moves 4 targets (98 / 104 / 98).

What it changes (six files; tercile / stercile / unset keep today's computation and today's meta):

1. `scripts/make_submission_pm_cd4.py` (the builder):
   - the per-target parse accepts `tercilec`;
   - `sharpen_q_by_tercile(..., carried=None)`: in `tercilec` the cut points, the floor and the tie blocks are read on the
     carried targets only, and every non-carried target sits in MID. With every target carried it IS `tercile` (tested).
     The audit adds n_carried, noncarried_to_mid and moved_vs_tercile ({target: [plain tier, tercilec tier]});
   - new `dump_key_carried(dump_env, perts, bulk, carriers, genes)`: `carried` comes from the key that scaled THE LOADED
     DUMP, i.e. the dump env's VCC_PERT_SCALE resolved through PERT_SCALE_ARTIFACTS, sha-checked, target order checked.
     Not this build's own VCC_PERT_SCALE: on D/E/F the level builds run the 410 key (DEF_ENV wshare_lo0.3_hi3_def) while
     the dump ran the LODO key (def_day.sh 2b); 479p2g's A/B/C build was the same (wshare_lo0.3_hi3 vs the dump's
     lodo415_lo0.5_hi3). The artifact must reproduce the dump's own ratio, r_t / w_t constant to 1e-6 relative over every
     target with a carrier (A/B/C: ~1e-14), or the build refuses: a wrong key would mis-assign without an error.
     The audit records carried_key, carried_key_sha256, key_ratio_c and key_ratio_maxdev;
   - the LOAD block passes `carried` for tercilec only.
2. `scripts/def_day.sh` (the POST-538+554+590 file): DEF_DAY_SHARPEN_PT accepts `tercilec:<q>,<q>,<q>`; a tercilec value
   needs the builder's dump_key_carried. Doc lines.
3. `scripts/def_day_helpers.py` meta_matches and 4. `scripts/def_pipeline.py` september_rows: both parse the per-target
   forms by an explicit tuple; without `tercilec` in it, the helper would read the value as a scalar q (mismatch, the
   re-level loop stops) and the driver would call float() on it after a good h5ad existed (499's failure mode).
5. `scripts/october_recipe.py`: sharpen_pt's form is checked (mode tercile|stercile|tercilec, three qs in (1, 3.5], as the
   builder parses it); a tercilec value needs the builder's dump_key_carried and def_day.sh's tercilec case.
6. `experiments/524-october-recipe-dr/gate524.py`: a tercilec run must show n_carried, noncarried_to_mid = 300 - n_carried,
   the key that scaled the dump (the LODO id def_day logged, else the build's own key) and key_ratio_maxdev <= 1e-6.

    uv run python experiments/598-tercilec-knob/patch598.py --root <repo> --check     # anchors only, writes nothing
    uv run python experiments/598-tercilec-knob/patch598.py --root <repo>             # apply (refuses a second time)

It refuses a tree where def_day.sh lacks 554's mark (the bundle must be applied first, in order), and a partial apply.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
import tempfile
from pathlib import Path

MARK = "598: tercilec"            # present in every patched file; a second apply is refused
MARK554 = "554: the K-class emitter"
BUILDER = Path("scripts/make_submission_pm_cd4.py")
DAY = Path("scripts/def_day.sh")
HELPERS = Path("scripts/def_day_helpers.py")
DRIVER = Path("scripts/def_pipeline.py")
RECIPE = Path("scripts/october_recipe.py")
GATE = Path("experiments/524-october-recipe-dr/gate524.py")
FILES = (BUILDER, DAY, HELPERS, DRIVER, RECIPE, GATE)

# ------------------------------------------------------------------------------------------------------------- builder
B_PARSE_A = '''# from the LOADED stack dump (sharpen_q_by_tercile). The scalar form below is unchanged.
BULK_SHARPEN_PT: dict[str, object] | None = None
BULK_SHARPEN_PT_STATE: dict[str, object] = {}
if BULK_SHARPEN_Q is not None and BULK_SHARPEN_Q.split(":", 1)[0].strip().lower() in ("tercile", "stercile"):
    _pt_mode, _pt_q = BULK_SHARPEN_Q.split(":", 1)
    try:
        _pt_qv = [float(x) for x in _pt_q.split(",")]
    except ValueError:
        raise SystemExit(f"VCC_BULK_SHARPEN_Q per-target form must be tercile:<q>,<q>,<q> or stercile:...: {BULK_SHARPEN_Q!r}")
'''
B_PARSE_B = '''# from the LOADED stack dump (sharpen_q_by_tercile). The scalar form below is unchanged.
# 598: tercilec (2026-10-02) is `tercile` with 493h's carried-only cut points: the cuts, the floor and the tie blocks are
# read on the targets the dump's OWN key artifact marks carried (dump_key_carried: the dump env's VCC_PERT_SCALE, resolved
# and sha-checked through the registry); every non-carried target sits in MID.
BULK_SHARPEN_PT: dict[str, object] | None = None
BULK_SHARPEN_PT_STATE: dict[str, object] = {}
if BULK_SHARPEN_Q is not None and BULK_SHARPEN_Q.split(":", 1)[0].strip().lower() in ("tercile", "stercile", "tercilec"):
    _pt_mode, _pt_q = BULK_SHARPEN_Q.split(":", 1)
    try:
        _pt_qv = [float(x) for x in _pt_q.split(",")]
    except ValueError:
        raise SystemExit(f"VCC_BULK_SHARPEN_Q per-target form must be tercile:<q>,<q>,<q>, stercile:... or tercilec:...: {BULK_SHARPEN_Q!r}")
'''

B_FUNC_A = '''def sharpen_q_by_tercile(bulk, carriers, genes, perts, pt) -> tuple[np.ndarray, dict[str, object]]:
    """479: a per-target sharpening q from the loaded stack dump's arrays (see VCC_BULK_SHARPEN_Q's per-target forms)."""
    panel = set(map(str, perts))
    nonpanel = np.fromiter((str(g) not in panel for g in genes), dtype=bool, count=len(genes))
    cn = np.linalg.norm(np.asarray(carriers, dtype=np.float64)[:, nonpanel], axis=1)
    floor = np.zeros(len(cn), dtype=bool)
    if pt["mode"] == "stercile":
        v = cn
    else:
        bn = np.linalg.norm(np.asarray(bulk, dtype=np.float64)[:, nonpanel], axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            v = np.where(cn > 0, bn / np.where(cn > 0, cn, 1.0), np.nan)
        if not np.isfinite(v).any():
            raise SystemExit("VCC_BULK_SHARPEN_Q tercile: no finite key ratio in the loaded dump")
        vmin = float(np.nanmin(v))
        floor = np.isfinite(v) & (v <= vmin * (1.0 + 1e-9))
    # 511 (501 section 6): cuts over finite values only; tie blocks (within 1e-9 relative, which also absorbs the ~1e-14
    # float noise of r = c*w) are assigned explicitly instead of by digitize's edge semantics.
    finite = np.isfinite(v)
    if not finite.any():
        raise SystemExit("VCC_BULK_SHARPEN_Q: no finite value in the loaded dump")
    sample = v[finite]
    fill = float(np.median(sample))
    vv = np.where(finite, v, fill)
    edges = np.quantile(sample, [1.0 / 3.0, 2.0 / 3.0])
    t = np.digitize(vv, edges)   # right=False; the tie blocks below override it
    tie_mid = np.zeros(len(vv), dtype=bool)

    def _near(e):
        return np.abs(vv - e) <= abs(float(e)) * 1e-9
    if pt["mode"] == "tercile":
'''
B_FUNC_B = '''def sharpen_q_by_tercile(bulk, carriers, genes, perts, pt, carried=None) -> tuple[np.ndarray, dict[str, object]]:
    """479: a per-target sharpening q from the loaded stack dump's arrays (see VCC_BULK_SHARPEN_Q's per-target forms).
    598: tercilec is `tercile` read on the targets `carried` marks (the dump key artifact's flags, dump_key_carried): the
    cut points, the floor and the tie blocks use carried targets only, and every non-carried target sits in MID (493h).
    With every target carried it equals `tercile` exactly. `tercile` and `stercile` ignore `carried`."""
    panel = set(map(str, perts))
    nonpanel = np.fromiter((str(g) not in panel for g in genes), dtype=bool, count=len(genes))
    cn = np.linalg.norm(np.asarray(carriers, dtype=np.float64)[:, nonpanel], axis=1)
    floor = np.zeros(len(cn), dtype=bool)
    car = np.ones(len(cn), dtype=bool)   # tercile / stercile: every target takes part, i.e. today's computation
    if pt["mode"] == "tercilec":
        car = np.asarray(carried, dtype=bool) if carried is not None else np.zeros(0, dtype=bool)
        if car.shape != cn.shape or not car.any():
            raise SystemExit("VCC_BULK_SHARPEN_Q tercilec: needs the dump key's carried flags, one per target, at least one True")
    if pt["mode"] == "stercile":
        v = cn
    else:
        bn = np.linalg.norm(np.asarray(bulk, dtype=np.float64)[:, nonpanel], axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            v = np.where(cn > 0, bn / np.where(cn > 0, cn, 1.0), np.nan)
        if not (np.isfinite(v) & car).any():
            raise SystemExit(f"VCC_BULK_SHARPEN_Q {pt['mode']}: no finite key ratio in the loaded dump")
        vmin = float(np.nanmin(np.where(car, v, np.nan)))
        floor = np.isfinite(v) & car & (v <= vmin * (1.0 + 1e-9))
    # 511 (501 section 6): cuts over finite values only; tie blocks (within 1e-9 relative, which also absorbs the ~1e-14
    # float noise of r = c*w) are assigned explicitly instead of by digitize's edge semantics.
    finite = np.isfinite(v)
    if not finite.any():
        raise SystemExit("VCC_BULK_SHARPEN_Q: no finite value in the loaded dump")
    sample = v[finite & car]
    fill = float(np.median(sample))
    vv = np.where(finite, v, fill)
    edges = np.quantile(sample, [1.0 / 3.0, 2.0 / 3.0])
    t = np.digitize(vv, edges)   # right=False; the tie blocks below override it
    tie_mid = np.zeros(len(vv), dtype=bool)

    def _near(e):
        return car & (np.abs(vv - e) <= abs(float(e)) * 1e-9)
    if pt["mode"] in ("tercile", "tercilec"):
'''

B_TAIL_A = '''    t = np.where(finite, t, 1)       # a non-finite value does not enter the cuts; it sits in MID
    q = np.asarray(pt["q"], dtype=np.float64)[t]
    audit = {"mode": pt["mode"], "q": list(pt["q"]), "edges": [float(e) for e in edges],
             "counts": [int((t == i).sum()) for i in range(3)], "floor_to_low": int(floor.sum()),
             "tie_to_mid": int(tie_mid.sum()), "nonfinite": int((~np.isfinite(v)).sum())}
    return q, audit
'''
B_TAIL_B = '''    t = np.where(finite, t, 1)       # a non-finite value does not enter the cuts; it sits in MID
    if pt["mode"] == "tercilec":
        t = np.where(car, t, 1)      # 493h: a non-carried target has no key (w = 1); it sits in MID
    q = np.asarray(pt["q"], dtype=np.float64)[t]
    audit = {"mode": pt["mode"], "q": list(pt["q"]), "edges": [float(e) for e in edges],
             "counts": [int((t == i).sum()) for i in range(3)], "floor_to_low": int(floor.sum()),
             "tie_to_mid": int(tie_mid.sum()), "nonfinite": int((~np.isfinite(v)).sum())}
    if pt["mode"] == "tercilec":
        # the plain rule's tiers on the same arrays (q = tier index), so the record shows exactly what the cut moved
        tp = np.rint(sharpen_q_by_tercile(bulk, carriers, genes, perts, {"mode": "tercile", "q": [0.0, 1.0, 2.0]})[0]).astype(int)
        audit.update({"n_carried": int(car.sum()), "noncarried_to_mid": int((~car).sum()),
                      "moved_vs_tercile": {str(perts[i]): [int(tp[i]), int(t[i])] for i in np.flatnonzero(tp != t)}})
    return q, audit


def dump_key_carried(dump_env, perts, bulk, carriers, genes) -> tuple[np.ndarray, dict[str, object]]:
    """598: the `carried` flags of the key artifact that scaled the LOADED stack dump: the dump env's VCC_PERT_SCALE,
    resolved through PERT_SCALE_ARTIFACTS and sha-checked like apply_pert_scale (on D/E/F the day-of LODO id is in the
    registry from def_day.sh's pin until its restore trap). This build's own VCC_PERT_SCALE is not used: the level builds
    run the 410 key while the dump ran the LODO key. The dump's own key ratio must be reproduced, ||bulk_t|| / ||carriers_t||
    = c * w_t over non-panel columns, constant to 1e-6 relative on every target with a carrier (A/B/C ~1e-14), or the
    build refuses: a wrong key would mis-assign the terciles without an error."""
    kid = str(dump_env.get("VCC_PERT_SCALE", "") or "").strip().lower()
    if kid not in PERT_SCALE_ARTIFACTS:
        raise SystemExit(f"VCC_BULK_SHARPEN_Q tercilec: the loaded dump's VCC_PERT_SCALE {kid!r} is not a registered key")
    path, pinned = PERT_SCALE_ARTIFACTS[kid]
    got = file_sha256(path)
    if got != pinned:
        raise SystemExit(f"VCC_BULK_SHARPEN_Q tercilec: key artifact sha256 {got} != pinned {pinned}")
    with np.load(path, allow_pickle=False) as z:
        if "carried" not in z.files:
            raise SystemExit(f"VCC_BULK_SHARPEN_Q tercilec: {path.name} has no carried array")
        targets = [str(t) for t in z["targets"]]
        w = np.asarray(z["w"], dtype=np.float64)
        carried = np.asarray(z["carried"], dtype=bool)
    if targets != [str(t) for t in perts] or w.shape != carried.shape or len(w) != len(targets):
        raise SystemExit("VCC_BULK_SHARPEN_Q tercilec: key artifact target order differs from the panel")
    panel = set(map(str, perts))
    nonpanel = np.fromiter((str(g) not in panel for g in genes), dtype=bool, count=len(genes))
    bn = np.linalg.norm(np.asarray(bulk, dtype=np.float64)[:, nonpanel], axis=1)
    cn = np.linalg.norm(np.asarray(carriers, dtype=np.float64)[:, nonpanel], axis=1)
    ok = (cn > 0) & np.isfinite(w) & (w > 0)
    if not ok.any():
        raise SystemExit("VCC_BULK_SHARPEN_Q tercilec: no target with a carrier and a positive key weight")
    rc = bn[ok] / (cn[ok] * w[ok])
    c = float(np.median(rc))
    dev = float(np.max(np.abs(rc / c - 1.0))) if c > 0 else float("inf")
    if not (np.isfinite(dev) and dev <= 1e-6):
        raise SystemExit(f"VCC_BULK_SHARPEN_Q tercilec: key {kid} does not reproduce the dump's ratio ||bulk||/||carriers|| "
                         f"= c*w (max relative deviation {dev:.3e} > 1e-6)")
    return carried, {"carried_key": kid, "carried_key_sha256": got, "key_ratio_c": c, "key_ratio_maxdev": dev}
'''

B_LOAD_A = '''        if BULK_SHARPEN_PT is not None:
            _qt, _pt_audit = sharpen_q_by_tercile(loaded, _pt_carriers, spec.genes, spec.perts, BULK_SHARPEN_PT)
            BULK_SHARPEN_PT_STATE.update({"qt": _qt, "audit": _pt_audit})
'''
B_LOAD_B = '''        if BULK_SHARPEN_PT is not None:
            # 598: tercilec reads `carried` from the key that scaled THIS dump (dump_env), never this build's VCC_PERT_SCALE
            _pt_carried, _pt_key = ((None, None) if BULK_SHARPEN_PT["mode"] != "tercilec"
                                    else dump_key_carried(dump_env, spec.perts, loaded, _pt_carriers, spec.genes))
            _qt, _pt_audit = sharpen_q_by_tercile(loaded, _pt_carriers, spec.genes, spec.perts, BULK_SHARPEN_PT,
                                                  carried=_pt_carried)
            if _pt_key is not None:
                _pt_audit.update(_pt_key)
            BULK_SHARPEN_PT_STATE.update({"qt": _qt, "audit": _pt_audit})
'''

# --------------------------------------------------------------------------------------------------------- def_day.sh
D_DOC_A = '''#   DEF_DAY_SHARPEN_PT=<tercile|stercile>:<q>,<q>,<q> (497, default off; e.g. stercile:2.5,2.0,1.75 = 479's S_P2): the
#    level builds sharpen each target's bulk with its tercile's q (VCC_BULK_SHARPEN_Q's per-target form, builder >=
#    d09a986c) instead of the scalar q. The stack dump is unchanged (q is a render knob).
'''
D_DOC_B = '''#   DEF_DAY_SHARPEN_PT=<tercile|stercile|tercilec>:<q>,<q>,<q> (497, default off; e.g. stercile:2.5,2.0,1.75 = 479's S_P2): the
#    level builds sharpen each target's bulk with its tercile's q (VCC_BULK_SHARPEN_Q's per-target form, builder >=
#    d09a986c) instead of the scalar q. The stack dump is unchanged (q is a render knob).
#    598: tercilec = tercile with 493h's carried-only cut points (the carried flags of the key that scaled the dump;
#    non-carried targets in MID). Builder >= 598's patch (dump_key_carried).
'''
D_CASE_A = '''  case "$DEF_DAY_SHARPEN_PT" in tercile:*,*,*|stercile:*,*,*) ;; *) die "DEF_DAY_SHARPEN_PT must be tercile:<q>,<q>,<q> or stercile:<q>,<q>,<q> (got $DEF_DAY_SHARPEN_PT)" ;; esac
  grep -q "def sharpen_q_by_tercile" scripts/make_submission_pm_cd4.py || die "the builder has no per-target sharpening (479)"
'''
D_CASE_B = '''  case "$DEF_DAY_SHARPEN_PT" in tercile:*,*,*|stercile:*,*,*|tercilec:*,*,*) ;; *) die "DEF_DAY_SHARPEN_PT must be tercile:<q>,<q>,<q>, stercile:<q>,<q>,<q> or tercilec:<q>,<q>,<q> (got $DEF_DAY_SHARPEN_PT)" ;; esac
  grep -q "def sharpen_q_by_tercile" scripts/make_submission_pm_cd4.py || die "the builder has no per-target sharpening (479)"
  # 598: tercilec reads the carried flags of the dump's own key artifact (the builder's dump_key_carried)
  case "$DEF_DAY_SHARPEN_PT" in tercilec:*) grep -q "def dump_key_carried" scripts/make_submission_pm_cd4.py || die "the builder has no tercilec (598)" ;; esac
'''

# ------------------------------------------------------------------------------------------- def_day_helpers.py / driver
H_DOC_A = "  def_day_helpers.py meta_matches <meta.json> <factors> <q | tercile:q,q,q | stercile:q,q,q>\n"
H_DOC_B = "  def_day_helpers.py meta_matches <meta.json> <factors> <q | tercile:q,q,q | stercile:q,q,q | tercilec:q,q,q>\n"
H_TUPLE_A = '            if qarg.split(":", 1)[0] in ("tercile", "stercile"):\n'
H_TUPLE_B = '            if qarg.split(":", 1)[0] in ("tercile", "stercile", "tercilec"):   # 598: tercilec, the carried-only key tercile\n'
P_TUPLE_A = '    if pt_mode in ("tercile", "stercile"):\n'
P_TUPLE_B = '    if pt_mode in ("tercile", "stercile", "tercilec"):   # 598: tercilec, the carried-only key tercile\n'

# --------------------------------------------------------------------------------------------------- october_recipe.py
R_CALL_A = "    _problems_538(env, problems)\n"
R_CALL_B = "    _problems_538(env, problems)\n    _problems_598(env, problems)\n"
R_FUNC_A = "\n\ndef installed_problems(env: dict, flags: str) -> list:\n"
R_FUNC_B = '''

def _problems_598(env: dict, problems: list) -> None:
    """598: tercilec (and the two older forms): DEF_DAY_SHARPEN_PT as the builder parses it, three qs in (1, 3.5]."""
    if "DEF_DAY_SHARPEN_PT" not in env:
        return
    raw = env["DEF_DAY_SHARPEN_PT"]
    mode, _, qs = raw.partition(":")
    try:
        qv = [float(x) for x in qs.split(",")]
        ok = mode in ("tercile", "stercile", "tercilec") and len(qv) == 3 and all(1.0 < x <= 3.5 for x in qv)
    except ValueError:
        ok = False
    if not ok:
        problems.append(f"sharpen_pt {raw!r} must be tercile:, stercile: or tercilec:<q>,<q>,<q> with every q in (1, 3.5]")


def installed_problems(env: dict, flags: str) -> list:
'''
R_NEEDS_A = ('    out += [f"the builder lacks {need} (needed by {var})" for var, need in BUILDER_NEEDS.items() '
             'if var in env and need not in bu]\n')
R_NEEDS_B = R_NEEDS_A + '''    # 598: tercilec needs the builder's dump_key_carried and def_day.sh's tercilec case
    if env.get("DEF_DAY_SHARPEN_PT", "").startswith("tercilec:"):
        if "def dump_key_carried" not in bu:
            out.append("the builder lacks dump_key_carried (needed by DEF_DAY_SHARPEN_PT=tercilec, 598)")
        if "tercilec:*,*,*" not in dd:
            out.append("scripts/def_day.sh does not accept DEF_DAY_SHARPEN_PT=tercilec (598)")
'''

# ---------------------------------------------------------------------------------------------------------- gate524.py
G_PT_A = '''    chk("per-target sharpening", pt.get("mode") == mode and pt.get("q") == [float(x) for x in qs.split(",")],
        f"meta {json.dumps({k: pt.get(k) for k in ('mode', 'q', 'counts')})}")
'''
G_PT_B = G_PT_A + '''    if mode == "tercilec":
        # 598: tercilec reads the carried flags of the key that scaled the dump: the LODO id def_day logged, else the
        # build's own key (the dump subshell then inherits it)
        lodo_ids = re.findall(r"LODO key ON for the stack dump: (\\S+)", dd)
        want_key = lodo_ids[-1] if lodo_ids else c.get("pert_scale")
        n_car, cnts = pt.get("n_carried"), pt.get("counts") or []
        try:
            dev_ok = float(pt.get("key_ratio_maxdev", 1.0)) <= 1e-6
        except (TypeError, ValueError):
            dev_ok = False
        chk("carried-only cut points (598)",
            isinstance(n_car, int) and 0 < n_car <= sum(cnts) and pt.get("noncarried_to_mid") == sum(cnts) - n_car
            and pt.get("carried_key") == want_key and dev_ok,
            f"n_carried {n_car} noncarried_to_mid {pt.get('noncarried_to_mid')} key {pt.get('carried_key')} (want {want_key}) "
            f"maxdev {pt.get('key_ratio_maxdev')} moved_vs_tercile {len(pt.get('moved_vs_tercile') or {})}")
'''


def edits() -> dict[Path, list[tuple[str, str]]]:
    return {
        BUILDER: [(B_PARSE_A, B_PARSE_B), (B_FUNC_A, B_FUNC_B), (B_TAIL_A, B_TAIL_B), (B_LOAD_A, B_LOAD_B)],
        DAY: [(D_DOC_A, D_DOC_B), (D_CASE_A, D_CASE_B)],
        HELPERS: [(H_DOC_A, H_DOC_B), (H_TUPLE_A, H_TUPLE_B)],
        DRIVER: [(P_TUPLE_A, P_TUPLE_B)],
        RECIPE: [(R_CALL_A, R_CALL_B), (R_FUNC_A, R_FUNC_B), (R_NEEDS_A, R_NEEDS_B)],
        GATE: [(G_PT_A, G_PT_B)],
    }


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def patched_texts(root: Path) -> dict[Path, str]:
    """The six new texts (nothing written). Refuses a pre-bundle tree, a partial apply and any anchor miss."""
    texts = {f: (root / f).read_text() for f in FILES}
    marked = [str(f) for f in FILES if MARK in texts[f]]
    if marked:
        raise SystemExit(f"already patched ({', '.join(marked)}): refusing a second apply" if len(marked) == len(FILES)
                         else f"partially patched ({', '.join(marked)}): restore the six files and re-apply")
    if MARK554 not in texts[DAY]:
        raise SystemExit("scripts/def_day.sh has no 554 mark: apply the 10-05 bundle (538 -> 554 -> 590) first; "
                         "598 anchors on the post-bundle def_day.sh")
    new = {}
    for f, pairs in edits().items():
        t = texts[f]
        for old, rep in pairs:
            n = t.count(old)
            if n != 1:
                raise SystemExit(f"{f}: anchor found {n} times (need 1): {old.splitlines()[0][:90]!r}")
            t = t.replace(old, rep)
        if MARK not in t:
            raise SystemExit(f"{f}: the patched text lacks the mark {MARK!r} (patch bug)")
        new[f] = t
    if new[BUILDER].count("PERT_SCALE_ARTIFACTS = {\n") != 1:   # def_pin_panel_artifacts.py's anchor must stay unique
        raise SystemExit("builder: the registry anchor would not be unique after the patch (patch bug)")
    return new


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    root = a.root.resolve()
    new = patched_texts(root)
    for f in FILES:
        print(f"{f}\n  before {sha(root / f)}")
    if a.check:
        print("patches apply cleanly (not written)")
        return
    for f in FILES:                               # atomic per file: write a sibling temp, then os.replace
        dst = root / f
        fd, tmp = tempfile.mkstemp(dir=dst.parent, prefix=f".{dst.name}.598-")
        with os.fdopen(fd, "w") as fh:
            fh.write(new[f])
        os.chmod(tmp, os.stat(dst).st_mode)
        os.replace(tmp, dst)
        print(f"{f}\n  after  {sha(dst)}")
    print("patched " + ", ".join(str(f) for f in FILES))


if __name__ == "__main__":
    sys.exit(main())
