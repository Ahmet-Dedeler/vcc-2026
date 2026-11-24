#!/usr/bin/env python3
"""696 R-ABC3 slot gates (gate695 with R-ABC2 as the reference; "adopted" = new vs R-ABC2's recipe):
695 R-ABC2 slot gates: gate694 (gate524's checks vs 453dr + the UPRN / DNMAG rows) with the dump-hook row added and the
strict diff now vs R-ABC (694), content-aware via scripts/def_metadiff.py: only the keys the adopted switches act on may differ
(θ 0.25, 707's hook, K), and each must move its own. Tag [rabc3]; arm1015 greps "^[rabc3] PASS" / "^[rabc3] CLIFF OK".
Everything below the line is gate524's docstring.
----
524dr's slot gates: the October recipe AS RESOLVED on 10-05 (`october_recipe.py command --allow-pending`, saved by
recipe_dr524.sh as work/recipe524.json) run by def_day.sh on 180's relabelled A/B/C bundle, vs 453dr (the plain October
path, board entry q3VuTImh4FRGKXH2waCp). Writes work/gates524.log with --write. PASS requires:
  - def_day's postbuild.json present with every check passing, converged (or its nearest safe build), every rho <= 1;
  - the zero-gene audit at 0;
  - the relabelled meta's config differs from 453dr's ONLY in keys the resolved switches (or the known builder changes
    since 453dr) explain, and every switch that is ON shows up in the meta with its resolved value;
  - the same-line record agrees with the recipe (w, no `pt`, programme on iff DEF_DAY_PROGRAMME_SOURCE is set);
  - recipe_dr524's novelty rc 0 and a package sha that matches the file.
usage: uv run python experiments/696-rabc3/gate696.py [--work W] [--stem S] [--chainlog L] [--write]"""
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
X = ROOT / "experiments/696-rabc3"; W = X / "work"; RUN = W / "run"
STEM = ROOT / "submissions/696rabc3-october-recipe-abc"
CHAINLOG = W / "rabc3_696.log"
# test overrides (the smoke test replays 497dr's finished rehearsal): --work <dir> (holds recipe524.json and run/),
# --stem <submissions stem>, --chainlog <log with the package/novelty lines>. Production passes none of them.
for _i, _a in enumerate(sys.argv):
    if _a == "--work":
        W = Path(sys.argv[_i + 1]); RUN = W / "run"
    elif _a == "--stem":
        STEM = Path(sys.argv[_i + 1])
    elif _a == "--chainlog":
        CHAINLOG = Path(sys.argv[_i + 1])
PKG, H5 = Path(str(STEM) + ".vcc"), Path(str(STEM) + ".h5ad")
REF = ROOT / "submissions/453dr-october-path-abc-q1.75.h5ad.meta.json"
# keys any 10-05 rehearsal differs from 453dr in, whatever the recipe: the re-level factors, the stack dump (LODO is
# DECIDED on; 453dr ran without it), the same-line record (453dr ran w auto + pt; the recipe runs w 0.5 without pt), and
# the emitter pin (496p's 1658ae2c landed after 453dr's edb5b8dc; no output effect without VCC_DE_MDE_SIGN, 499 F2).
BASE_ALLOWED = {"bulk_ctx_scale", "bulk_load", "bulk_load_sha256", "sameline", "de_sign_emitter_sha256"}
# DEF_DAY_* export -> the meta keys it may change
SWITCH_KEYS = {
    "DEF_DAY_LODO": {"bulk_load", "bulk_load_sha256"},
    "DEF_DAY_LODO_HI_PCT": {"bulk_load", "bulk_load_sha256"},
    "DEF_DAY_SHARPEN_PT": {"bulk_sharpen_pt", "bulk_sharpen_q"},
    "DEF_DAY_PROGRAMME_SOURCE": {"sameline"},
    "DEF_DAY_PRIOR_GATE": {"de_prior_gate"},
    "DEF_DAY_DE_MAG": {"de_mag"},
    "DEF_DAY_MDE_SIGN": {"de_mde_sign"},
    "DEF_DAY_RC_CALLSET_CPM": {"rc_callset_cpm"},
    "DEF_DAY_RHO_SCALE": {"bulk_ctx_scale"},
    "DEF_DAY_CD4_STACK": {"bulk_load", "bulk_load_sha256"},
    "DEF_DAY_LODO_STRENGTH": {"bulk_load", "bulk_load_sha256", "pert_scale", "pert_scale_sha256"},
    "DEF_DAY_HOOK_ESTIMATOR": {"de_mde_sign"},
    "DEF_DAY_LODO_HI": {"bulk_load", "bulk_load_sha256"},
    "DEF_DAY_CALLSET_BLOCK": {"rc_callset_block"},
    "DEF_DAY_DE_UPRN": {"de_uprn"},      # 630: rank-neutral UP pass (rank_up module)
    "DEF_DAY_DE_DNMAG": {"de_dnmag"},    # 631: DOWN magnitude pass (dn_mag2 module since 693p)
    "DEF_DAY_DUMP_HOOK": {"bulk_load_sha256"},   # 690g: the day-of dump hook (707)
}
# 694: the board-scored 524dr2 (entry dVHLws6KdMxudq0Yn5ap) and the reviewed module shas def_day pins
REF2 = ROOT / "submissions/524dr2-october-recipe-abc.h5ad.meta.json"
UPRN_SHA = "7e6b3e52e2ad374b7b6e21e92327fcf65771433693eb3846c09c0189f3f2959f"
DNMAG_SHA = "b9d6cfcf17ad783c1ea2607f95366173a0fc1f7920c4104604f9ee925ac55c87"
out, ok = [], True
rec = json.loads((W / "recipe696.json").read_text())
env, args = rec["env"], rec["args"]
unknown = sorted(set(env) - set(SWITCH_KEYS))
if unknown or rec.get("problems"):
    print("recipe has unknown exports or problems:", unknown, rec.get("problems")); sys.exit(1)
out.append(f"[rabc3] recipe: env {json.dumps(env, sort_keys=True)} args {args[3:]} pending-off {rec.get('pending_treated_off')}")
dd = (RUN / "def_day.log").read_text()
pb_path = RUN / "postbuild.json"
pb = json.loads(pb_path.read_text()) if pb_path.is_file() else {}
failed = [c.get("check") for c in pb.get("checks", []) if not c.get("pass")]
if not pb or failed:
    print("postbuild.json missing or failing:", failed or "no file"); sys.exit(1)
conv = re.findall(r"CONVERGED at (\S+) after", dd) or re.findall(r"NOT converged .*nearest safe build: (\S+)", dd)
fin = re.findall(r"final: (\S+) factors (\S+) ->", dd)
if not fin:
    print("def_day has no final line"); sys.exit(1)
final_name, factors = fin[-1]
nt = json.loads((RUN / f"{final_name}.ntclevel.json").read_text())
rho = {c: round(v["rho"], 4) for c, v in nt["per_context"].items()}
cliff_ok = all(v <= 1.0 for v in rho.values())
ok &= cliff_ok and bool(conv)
out.append(f"[rabc3] {'CLIFF OK' if cliff_ok else 'CLIFF FAIL'} {json.dumps(rho)} final {final_name} factors {factors} "
           f"({'converged' if 'CONVERGED at' in dd else 'nearest safe build'}); postbuild {len(pb.get('checks', []))}/"
           f"{len(pb.get('checks', []))} PASS")
audit_path = RUN / "audit428.txt"
audit_lines = audit_path.read_text().splitlines() if audit_path.is_file() else []
zero_ok = bool(audit_lines) and re.search(r"ZERO block total \(target excluded\): 0 in 0 blocks", audit_lines[0]) is not None
out.append(f"zero-gene audit: {'0 (OK)' if zero_ok else 'NOT 0 / unparsed'} :: {audit_lines[0][:160] if audit_lines else '<no audit428.txt>'}")
ok &= zero_ok
meta = json.loads(Path(str(H5) + ".meta.json").read_text()); ref = json.loads(REF.read_text())
c, r = dict(meta["config"]), dict(ref["config"])
diff = sorted(k for k in set(c) | set(r) if c.get(k) != r.get(k))
allowed = set(BASE_ALLOWED)
for k in env:
    allowed |= SWITCH_KEYS[k]
extra = [k for k in diff if k not in allowed]
out.append(f"[rabc3] {'PASS' if not extra else 'FAIL'}: config diff vs 453dr {diff}" + (f"  UNEXPECTED {extra}" if extra else ""))
ok &= not extra
# 695: vs R-ABC (694; board 0.2743318, the October recipe before the 10-17/18 reads), content-aware (scripts/def_metadiff.py):
# only the keys the ADOPTED switches act on may differ, and each adopted switch must move its own keys.
import subprocess  # noqa: E402
import os  # noqa: E402
# test-only overrides (test696 points them at R-ABC, since R-ABC2 does not exist yet); production passes none
REF_ABC = Path(os.environ.get("GATE696_REF_META", ROOT / "submissions/695rabc2-october-recipe-abc.h5ad.meta.json"))   # 696: vs R-ABC2
REF_ABC_RUN = Path(os.environ.get("GATE696_REF_RUN", ROOT / "experiments/695-rabc2/work/run"))
ACTS = {   # export -> (allowed keys, keys that must move)
    "DEF_DAY_LODO_STRENGTH": ({"bulk_load", "bulk_load_sha256", "pert_scale", "pert_scale_sha256", "bulk_ctx_scale", "bulk_sharpen_pt"},
                              {"bulk_load"}),   # the dump renamed _s0.25; the level builds keep their own pert_scale
    "DEF_DAY_DUMP_HOOK": ({"bulk_load_sha256", "bulk_ctx_scale"}, {"bulk_load_sha256"}),
    "DEF_DAY_RC_CALLSET_CPM": ({"rc_callset_cpm"}, {"rc_callset_cpm"}),
}
# 696: "adopted" = exports NEW relative to R-ABC2 (its own θ / hook are the reference, not changes)
REF_RECIPE = Path(os.environ.get("GATE696_REF_RECIPE", ROOT / "experiments/695-rabc2/work/recipe695.json"))
_ref_env = json.loads(REF_RECIPE.read_text())["env"] if REF_RECIPE.is_file() else {}
adopted = sorted(k for k in env if k in ACTS and k not in _ref_env)
allow = sorted(set().union(*[ACTS[k][0] for k in adopted])) if adopted else []
need = sorted(set().union(*[ACTS[k][1] for k in adopted])) if adopted else []
md = subprocess.run([sys.executable, str(ROOT / "scripts/def_metadiff.py"), str(H5) + ".meta.json", str(RUN), str(REF_ABC),
                     str(REF_ABC_RUN), "--allow", ",".join(allow), "--need", ",".join(need)], capture_output=True, text=True)
mdl = (md.stdout.strip().splitlines() or ["(no output) " + md.stderr[-200:]])[-1]
ok2 = md.returncode == 0 and bool(adopted)
out.append(f"[rabc3] {'PASS' if ok2 else 'FAIL'}: config diff vs R-ABC2 (adopted {adopted}): {mdl[:300]}")
for n in md.stdout.strip().splitlines()[:-1]:
    out.append("  " + n.strip()[:200])
ok &= ok2


def _dump_name530(stack: str, strength: str | None) -> str:
    """Same filename rule as def_day.sh / def_day_helpers._dump_name530."""
    base = "dump_h1x2_cd4x15" if stack == "x1.5" else "dump_h1x2_cd4x1"
    if strength:
        base += "_s" + str(strength).replace("-", "m")
    return base + ".npz"


def _mde_ok530(spec: str, ms) -> bool:
    """526 compared source and tau after one partition, which is clean4:<tau> only.
    hookbank:<bank>:<tau>[:tau_c4[:leN]] is recorded with source file (the shell rewrites it).
    """
    if not isinstance(ms, dict):
        return False
    p = spec.split(":")
    try:
        if p[0] == "clean4" and len(p) == 2:
            return ms.get("source") == "clean4" and abs(float(ms.get("tau", -1)) - float(p[1])) < 1e-12
        if p[0] == "programme" and len(p) in (2, 3):
            tau_c4 = float(p[2]) if len(p) == 3 else 0.05
            return (ms.get("source") == "programme"
                    and abs(float(ms.get("tau", -1)) - float(p[1])) < 1e-12
                    and abs(float(ms.get("tau_c4", -1)) - tau_c4) < 1e-12)
        if p[0] == "hookbank" and len(p) in (3, 4, 5) and p[1]:
            tau, tau_c4 = float(p[2]), float(p[3]) if len(p) >= 4 else 0.05
            ok = (ms.get("source") == "file" and isinstance(ms.get("sha256"), str) and len(ms["sha256"]) == 64
                  and abs(float(ms.get("tau", -1)) - tau) < 1e-12
                  and abs(float(ms.get("tau_c4", -1)) - tau_c4) < 1e-12)
            if len(p) == 5:
                ok = ok and p[4].startswith("le") and ms.get("le") == int(p[4][2:])
            else:
                ok = ok and "le" not in ms
            return bool(ok)
    except (TypeError, ValueError, IndexError):
        return False
    return False


def chk(name: str, cond: bool, detail: str) -> None:
    global ok
    ok &= bool(cond)
    out.append(f"  {'ok  ' if cond else 'FAIL'} {name}: {detail}")


# every ON switch is visible in the final artifact with its resolved value
if env.get("DEF_DAY_LODO") == "1":
    chk("LODO key on the stack dump", "LODO key ON for the stack dump" in dd and "bulk_load_sha256" in diff,
        "def_day logged it; the dump differs from 453dr's")
if "DEF_DAY_LODO_HI_PCT" in env:
    side = RUN / "pert_scale_lodo_def_lo0.5_pct.json"
    hp = json.loads(side.read_text()).get("hi_pct") if side.is_file() else None
    chk("LODO percentile cap", "lodo_lo0.5_pct_def" in dd and hp is not None and float(hp) == float(env["DEF_DAY_LODO_HI_PCT"]),
        f"sidecar hi_pct {hp}")
if "DEF_DAY_SHARPEN_PT" in env:
    mode, qs = env["DEF_DAY_SHARPEN_PT"].split(":")
    pt = c.get("bulk_sharpen_pt") or {}
    chk("per-target sharpening", pt.get("mode") == mode and pt.get("q") == [float(x) for x in qs.split(",")],
        f"meta {json.dumps({k: pt.get(k) for k in ('mode', 'q', 'counts')})}")
    if mode == "tercilec":
        # 598: tercilec reads the carried flags of the key that scaled the dump: the LODO id def_day logged, else the
        # build's own key (the dump subshell then inherits it)
        lodo_ids = re.findall(r"LODO key ON for the stack dump: (\S+)", dd)
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
if "DEF_DAY_DE_MAG" in env:
    chk("DE magnitude rider", (c.get("de_mag") or {}).get("variant") == env["DEF_DAY_DE_MAG"],
        f"meta de_mag {json.dumps(c.get('de_mag'))[:160]}")
if "DEF_DAY_MDE_SIGN" in env:
    ms = c.get("de_mde_sign")
    tau_ok = _mde_ok530(env["DEF_DAY_MDE_SIGN"], ms)
    chk("uncalled-gene sign hook", tau_ok, f"meta de_mde_sign {json.dumps(ms)[:160]} vs {env['DEF_DAY_MDE_SIGN']}")
if "DEF_DAY_PRIOR_GATE" in env:
    pg = c.get("de_prior_gate")
    spec = env["DEF_DAY_PRIOR_GATE"]
    try:
        if spec.startswith("le"):
            gate_ok = isinstance(pg, dict) and pg.get("le") == int(spec[2:]) and "q" not in pg
        elif spec.startswith("q"):
            gate_ok = isinstance(pg, dict) and abs(float(pg.get("q", -1)) - int(spec[1:]) / 100.0) < 1e-12 and "le" not in pg
        else:
            gate_ok = False
    except (TypeError, ValueError):
        gate_ok = False
    chk("prior gate", gate_ok, f"meta de_prior_gate {json.dumps(pg)[:160]} vs {spec}")
if "DEF_DAY_RC_CALLSET_CPM" in env:
    chk("call-set CPM block", c.get("rc_callset_cpm") == int(env["DEF_DAY_RC_CALLSET_CPM"]), f"meta {c.get('rc_callset_cpm')}")
if "DEF_DAY_CD4_STACK" in env or "DEF_DAY_LODO_STRENGTH" in env:
    stack = env.get("DEF_DAY_CD4_STACK", "x1")
    theta = env.get("DEF_DAY_LODO_STRENGTH")
    want_dump = _dump_name530(stack, theta)
    g = c.get("cd4_target_g")
    try:
        dose_ok = g is not None and abs(float(g) - 0.3443655278649641) < 1e-12
    except (TypeError, ValueError):
        dose_ok = False
    chk("stack dump name", c.get("bulk_load") == want_dump and dose_ok,
        f"bulk_load {c.get('bulk_load')} cd4_target_g {g} vs {want_dump} (level g stays DEF_ENV)")
else:
    chk("CD4 stack x1", c.get("bulk_load") == "dump_h1x2_cd4x1.npz", f"bulk_load {c.get('bulk_load')}")
if "DEF_DAY_RHO_SCALE" in env:
    scale = env["DEF_DAY_RHO_SCALE"]
    chk("rho scale", re.search(rf"rho-scale={re.escape(scale)}(\s|$)", dd) is not None, f"def_day log names rho-scale={scale}")
if "DEF_DAY_LODO_STRENGTH" in env:
    side = RUN / "pert_scale_lodo_strength.json"
    theta = None
    if side.is_file():
        try:
            theta = json.loads(side.read_text()).get("theta")
        except (OSError, json.JSONDecodeError):
            theta = None
    try:
        th_ok = theta is not None and abs(float(theta) - float(env["DEF_DAY_LODO_STRENGTH"])) < 1e-12
    except (TypeError, ValueError):
        th_ok = False
    chk("LODO strength rider", th_ok and f"LODO strength rider: theta {env['DEF_DAY_LODO_STRENGTH']}" in dd,
        f"sidecar theta {theta}")
if "DEF_DAY_HOOK_ESTIMATOR" in env:
    est = env["DEF_DAY_HOOK_ESTIMATOR"]
    chk("hook estimator", est == "meanz_ne" and f"estimator {est}" in dd
        and str(env.get("DEF_DAY_MDE_SIGN", "")).startswith("hookbank:"),
        f"log names estimator {est}; mde {env.get('DEF_DAY_MDE_SIGN')}")
if "DEF_DAY_LODO_HI" in env:
    h = env["DEF_DAY_LODO_HI"]
    side = RUN / f"pert_scale_lodo_def_lo0.5_hi{h}.json"
    got_hi = None
    if side.is_file():
        try:
            got_hi = json.loads(side.read_text()).get("hi")
        except (OSError, json.JSONDecodeError):
            got_hi = None
    try:
        hi_ok = got_hi is not None and abs(float(got_hi) - float(h)) < 1e-12
    except (TypeError, ValueError):
        hi_ok = False
    chk("LODO fixed cap", hi_ok and f"lodo_lo0.5_hi{h}_def" in dd, f"sidecar hi {got_hi}")
if "DEF_DAY_CALLSET_BLOCK" in env:
    raw = env["DEF_DAY_CALLSET_BLOCK"]
    parts = raw.split(":")
    form_ok, want_swap, want_steer = False, None, None
    if raw == "prop":
        form_ok, want_swap, want_steer = True, 1.0, "kept"
    elif parts and parts[0] == "prop" and len(parts) in (2, 3):
        try:
            want_swap = float(parts[1])
            want_steer = parts[2].lower() if len(parts) == 3 else "kept"
            form_ok = 0.0 < want_swap <= 1.0 and want_steer in {"kept", "key", "cpm"}
        except ValueError:
            form_ok = False
    blk = c.get("rc_callset_block") or {}
    sha = blk.get("sha256")
    try:
        swap_ok = want_swap is not None and abs(float(blk.get("swap", -1.0)) - want_swap) < 1e-12
    except (TypeError, ValueError):
        swap_ok = False
    chk("call-propensity block",
        form_ok and blk.get("source") == "prop" and isinstance(sha, str) and len(sha) == 64
        and swap_ok and blk.get("steer") == want_steer,
        f"meta {json.dumps(blk)[:160]} vs {raw}")
if "DEF_DAY_DE_UPRN" in env:
    u = c.get("de_uprn") or {}
    chk("rank-neutral UP rider (630)", u.get("variant") == env["DEF_DAY_DE_UPRN"] and u.get("module_sha256") == UPRN_SHA,
        f"meta de_uprn variant {u.get('variant')} module {str(u.get('module_sha256'))[:16]}")
if "DEF_DAY_DE_DNMAG" in env:
    u = c.get("de_dnmag") or {}
    chk("DOWN magnitude rider (631, dn_mag2)", u.get("variant") == env["DEF_DAY_DE_DNMAG"] and u.get("module_sha256") == DNMAG_SHA,
        f"meta de_dnmag variant {u.get('variant')} module {str(u.get('module_sha256'))[:16]}")
if "DEF_DAY_DUMP_HOOK" in env:
    hk = env["DEF_DAY_DUMP_HOOK"]
    chk("dump hook (690g / 707)", f"690g: dump hook {hk}" in dd and "690g: dump replaced by the hook's" in dd,
        f"def_day log names the hook {hk} and the validated replacement dump")
# ... and every OFF switch is absent
for k, key in (("DEF_DAY_SHARPEN_PT", "bulk_sharpen_pt"), ("DEF_DAY_DE_MAG", "de_mag"), ("DEF_DAY_MDE_SIGN", "de_mde_sign"),
               ("DEF_DAY_PRIOR_GATE", "de_prior_gate"), ("DEF_DAY_DE_UPRN", "de_uprn"), ("DEF_DAY_DE_DNMAG", "de_dnmag")):
    if k not in env:
        chk(f"{key} off", c.get(key) is None, f"meta {json.dumps(c.get(key))[:80]}")
# the same-line record follows the recipe's arguments (w, flags) and the programme switch
sl = c.get("sameline")
w_arg, flags_arg = args[4], args[5]
if isinstance(sl, dict):
    want_prog = "DEF_DAY_PROGRAMME_SOURCE" in env
    ctxs = sl.get("contexts", {})
    ws = sorted({v.get("w") for v in ctxs.values()})
    chk("same-line flags", sl.get("pt") is ("pt" in flags_arg.split(":")) and sl.get("prog") is want_prog,
        f"pt {sl.get('pt')} prog {sl.get('prog')} (args w {w_arg!r} flags {flags_arg!r}; programme {'on' if want_prog else 'off'})")
    if w_arg != "auto":
        chk("same-line w", bool(ctxs) and all(x == float(w_arg) for x in ws), f"per-context w {ws}")
else:
    # Absence is not agreement. 453dr carried a record; BASE_ALLOWED lets the key
    # differ, so a dropped record used to pass (526).
    chk("same-line record", False, "no record (no context matched)")
if "DEF_DAY_PROGRAMME_SOURCE" in env:
    want_bank = env["DEF_DAY_PROGRAMME_SOURCE"]
    chk("programme source bank", f"programme-source {want_bank}" in dd, f"def_day log names {want_bank}")
out.append(f"relabel: {meta.get('relabel453', {}).get('categories_after')} cells {meta.get('relabel453', {}).get('cells')}")
cl = CHAINLOG.read_text()
nov = re.findall(r"novelty rc=(\d+)", cl); pk = re.findall(r"package (\S+) sha256 ([0-9a-f]{64})", cl)
h = hashlib.sha256()
with open(PKG, "rb") as f:
    for chunk in iter(lambda: f.read(1 << 24), b""):
        h.update(chunk)
sha_ok = bool(pk) and pk[-1][1] == h.hexdigest()
ok &= bool(nov) and nov[-1] == "0" and sha_ok
out.append(f"novelty rc={nov[-1] if nov else 'missing'}")
out.append(f"package sha256: {h.hexdigest()}" + ("" if sha_ok else "  (MISMATCH with the chain log)"))
print("\n".join(out)); print("ALL GATES", "PASS" if ok else "FAIL")
if "--write" in sys.argv:
    if not ok:
        print("not writing gates696.log (a gate failed)"); sys.exit(1)
    (W / "gates696.log").write_text("\n".join(out) + "\n")
    print("wrote", W / "gates696.log")
