"""October recipe resolver: turn `docs/october-decisions.json` into the exact `scripts/def_day.sh` call.

The runbook (`docs/DEF-RUNBOOK.md`) carries the reasoning for each October switch; the decisions file carries its value,
status and evidence. This script is the one place that maps switches to `def_day.sh` arguments and `DEF_DAY_*` exports,
so the Oct 22 command (and any A/B/C rehearsal of it) is resolved from the record instead of assembled by hand.

  uv run python scripts/october_recipe.py show
  uv run python scripts/october_recipe.py command --bundle data/challenge --work <dir> --name <stem> [--allow-pending] [--json]
  uv run python scripts/october_recipe.py decide <switch> <value> --evidence "<read, entry id, file>"
  uv run python scripts/october_recipe.py check [--allow-pending]

`command` refuses while any switch is PENDING unless `--allow-pending` (which treats PENDING as off and says so). Both
`command` and `check` verify that the installed `scripts/def_day.sh` and builder support every knob the recipe uses.
Exit codes: 0 ok, 2 pending switches, 3 an installed script lacks a knob or a dependency is violated, 1 usage error.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import shlex
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# The OCTOBER_* overrides exist for tests (tests/test_october_recipe.py); the day uses the defaults.
DECISIONS = Path(os.environ.get("OCTOBER_DECISIONS", ROOT / "docs/october-decisions.json"))
DEF_DAY = Path(os.environ.get("OCTOBER_DEF_DAY", ROOT / "scripts/def_day.sh"))
BUILDER = Path(os.environ.get("OCTOBER_BUILDER", ROOT / "scripts/make_submission_pm_cd4.py"))

# switch -> the DEF_DAY_* export it drives (None: a positional argument or a day-of note)
ENV_OF = {
    "lodo": "DEF_DAY_LODO",
    "lodo_hi_pct": "DEF_DAY_LODO_HI_PCT",
    "sharpen_pt": "DEF_DAY_SHARPEN_PT",
    "programme_source": "DEF_DAY_PROGRAMME_SOURCE",
    "prior_gate": "DEF_DAY_PRIOR_GATE",
    "de_mag": "DEF_DAY_DE_MAG",
    "mde_sign": "DEF_DAY_MDE_SIGN",
    "rho_scale": "DEF_DAY_RHO_SCALE",
    "rc_callset_cpm": "DEF_DAY_RC_CALLSET_CPM",
    "hook_estimator": "DEF_DAY_HOOK_ESTIMATOR",
    "lodo_strength": "DEF_DAY_LODO_STRENGTH",
    "lodo_hi": "DEF_DAY_LODO_HI",
    "callset_block": "DEF_DAY_CALLSET_BLOCK",
    # 690a: the late DE riders (def_day.sh carries all three; PENDING until a board read and a held-out read decide them)
    "de_uprn": "DEF_DAY_DE_UPRN",
    "mde_sign_keygate": "DEF_DAY_MDE_SIGN_KEYGATE",
    "de_dnmag": "DEF_DAY_DE_DNMAG",
    "dump_hook": "DEF_DAY_DUMP_HOOK",   # 690g: a committed script run after the stack dump (the PIE day-of step)
}
# knob -> text the installed builder must contain for the knob to do anything
BUILDER_NEEDS = {"DEF_DAY_DE_MAG": "VCC_DE_MAG", "DEF_DAY_MDE_SIGN": "VCC_DE_MDE_SIGN",
                 "DEF_DAY_PRIOR_GATE": "VCC_DE_PRIOR_GATE", "DEF_DAY_SHARPEN_PT": "stercile",
                 "DEF_DAY_RC_CALLSET_CPM": "VCC_RC_CALLSET_CPM",
                 "DEF_DAY_CALLSET_BLOCK": "VCC_RC_CALLSET_BLOCK",
                 # 690a
                 "DEF_DAY_DE_UPRN": "630: rank-neutral UP", "DEF_DAY_MDE_SIGN_KEYGATE": "VCC_DE_MDE_SIGN_KEYGATE",
                 "DEF_DAY_DE_DNMAG": "VCC_DE_DNMAG"}
SAMELINE = Path(os.environ.get("OCTOBER_SAMELINE", ROOT / "scripts/def_sameline.py"))
EMITTER = Path(os.environ.get("OCTOBER_EMITTER", ROOT / "src/vcc2026/rc_emit_k.py"))
CALLBLOCK = Path(os.environ.get("OCTOBER_CALLBLOCK", ROOT / "experiments/528-callprop-block/callblock528.py"))
OFF = {None, "", "off", "0", "none"}
EMPTY_FLAGS_MARK = 'SLF="${6-pt:noprog}"'  # def_day.sh >= e29e793a (497): "" means no flags, not pt:noprog


def load() -> dict:
    return json.loads(DECISIONS.read_text())


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def is_on(sw: dict) -> bool:
    return sw.get("status") == "DECIDED" and str(sw.get("value")).strip().lower() not in {str(x) for x in OFF}


def resolve(d: dict, allow_pending: bool) -> tuple[dict, list, list, list]:
    """Returns (env, positional args after bundle/work/name, pending switch names, problems)."""
    s = d["switches"]
    pending = [k for k, v in s.items() if v.get("status") == "PENDING"]
    problems: list[str] = []
    env: dict[str, str] = {}
    for k, var in ENV_OF.items():
        sw = s.get(k, {})
        if is_on(sw):
            env[var] = str(sw["value"])
    prog = "DEF_DAY_PROGRAMME_SOURCE" in env
    if "DEF_DAY_PRIOR_GATE" in env and not prog:
        problems.append("prior_gate is on without programme_source (def_day.sh refuses it)")
    if "DEF_DAY_LODO_HI_PCT" in env and env.get("DEF_DAY_LODO") != "1":
        problems.append("lodo_hi_pct needs lodo = 1")
    # 530: x1 is today's dump and is not exported, so a DECIDED x1 leaves the command unchanged.
    cd4 = s.get("cd4_stack", {})
    if cd4.get("status") == "DECIDED":
        cd4_val = str(cd4.get("value"))
        if cd4_val == "x1.5":
            env["DEF_DAY_CD4_STACK"] = "x1.5"
        elif cd4_val not in {"x1", "off", "0", "none"}:
            problems.append(f"cd4_stack {cd4_val!r} is not x1 or x1.5")
    _problems_530(env, problems)
    _problems_538(env, problems)
    _problems_598(env, problems)
    _problems_690a(env, problems)
    _problems_690g(env, problems)
    # is_on treats "0" as off, so a DECIDED theta of 0 would vanish instead of being refused.
    # 0 still rewrites the LODO npz (522); it is not a way to leave the rider off.
    _ls = s.get("lodo_strength", {})
    if _ls.get("status") == "DECIDED" and str(_ls.get("value")).strip() == "0":
        problems.append("lodo_strength '0' must be a non-zero theta in [-1, 1]")
    flags = s["sameline_flags"]["value"]
    if prog:
        flags = s["sameline_flags"].get("with_programme_source", "")
        if flags and "noprog" in flags:
            problems.append("programme_source is on but the flags keep noprog (the builder would ignore the programme)")
    args = [str(s["q"]["value"]), str(s["sameline_w"]["value"]), flags, str(s["max_levels"]["value"])]
    if pending and not allow_pending:
        problems.append("PENDING: " + ", ".join(pending))
    return env, args, pending, problems


def _problems_530(env: dict, problems: list) -> None:
    """Range and dependency checks for 530's knobs. PENDING stays out of env, so this does not run for them."""
    if "DEF_DAY_RHO_SCALE" in env:
        try:
            s = float(env["DEF_DAY_RHO_SCALE"])
            ok = 0.9 <= s <= 1.1
        except ValueError:
            ok = False
        if not ok:
            problems.append(f"rho_scale {env['DEF_DAY_RHO_SCALE']!r} must be in [0.9, 1.1]")
    if "DEF_DAY_RC_CALLSET_CPM" in env:
        n = env["DEF_DAY_RC_CALLSET_CPM"]
        if not (n.isdigit() and 100 <= int(n) <= 5000):
            problems.append(f"rc_callset_cpm {n!r} must be an integer in [100, 5000]")
    if "DEF_DAY_HOOK_ESTIMATOR" in env:
        if env["DEF_DAY_HOOK_ESTIMATOR"] != "meanz_ne":
            problems.append(f"hook_estimator {env['DEF_DAY_HOOK_ESTIMATOR']!r} must be meanz_ne")
        if not env.get("DEF_DAY_MDE_SIGN", "").startswith("hookbank:"):
            problems.append("hook_estimator needs DEF_DAY_MDE_SIGN=hookbank:<bank>:<tau>...")
    if "DEF_DAY_LODO_STRENGTH" in env:
        if env.get("DEF_DAY_LODO") != "1":
            problems.append("lodo_strength needs lodo = 1")
        if "DEF_DAY_LODO_HI_PCT" not in env:
            problems.append("lodo_strength needs the percentile cap (lodo_hi_pct)")
        import re as _re
        th = env["DEF_DAY_LODO_STRENGTH"]
        try:
            v = float(th)
            ok = _re.fullmatch(r"-?(?:0|[1-9]\d*)(?:\.\d+)?", th) is not None and v != 0.0 and abs(v) <= 1.0
        except ValueError:
            ok = False
        if not ok:
            problems.append(f"lodo_strength {th!r} must be a non-zero theta in [-1, 1]")
    mde = env.get("DEF_DAY_MDE_SIGN", "")
    if mde and not (mde.startswith("clean4:") or mde.startswith("programme:") or mde.startswith("hookbank:")):
        problems.append(f"mde_sign {mde!r} must be clean4:, programme: or hookbank:")


def _callset_ok538(raw: str) -> bool:
    if raw == "prop":
        return True
    p = raw.split(":")
    if not p or p[0] != "prop" or len(p) not in (2, 3):
        return False
    try:
        v = float(p[1])
    except ValueError:
        return False
    if not (0.0 < v <= 1.0):
        return False
    if len(p) == 3 and p[2].lower() not in {"kept", "key", "cpm"}:
        return False
    return True


def _problems_538(env: dict, problems: list) -> None:
    """Range and dependency checks for 538's knobs. PENDING stays out of env."""
    if "DEF_DAY_LODO_HI" in env:
        if env.get("DEF_DAY_LODO") != "1":
            problems.append("lodo_hi needs lodo = 1")
        if "DEF_DAY_LODO_HI_PCT" in env:
            problems.append("lodo_hi and lodo_hi_pct are two cap forms; set only one")
        import re as _re
        t = env["DEF_DAY_LODO_HI"]
        try:
            h = float(t)
            ok = _re.fullmatch(r"[1-9]\d*(?:\.\d+)?", t) is not None and 1.0 < h <= 20.0
        except ValueError:
            ok = False
        if not ok:
            problems.append(f"lodo_hi {t!r} must satisfy 1 < h <= 20")
    if "DEF_DAY_CALLSET_BLOCK" in env and not _callset_ok538(env["DEF_DAY_CALLSET_BLOCK"]):
        problems.append(f"callset_block {env['DEF_DAY_CALLSET_BLOCK']!r} must be prop or prop:<swap>[:<steer>]")


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


def _problems_690a(env: dict, problems: list) -> None:
    """690a: dependencies and grammar of the late DE riders, matching def_day.sh's own validation. PENDING stays out of env."""
    import re  # noqa: PLC0415
    if "DEF_DAY_DE_UPRN" in env:
        if "DEF_DAY_DE_MAG" not in env:
            problems.append("de_uprn runs after the zs pass: it needs de_mag")
        m = re.fullmatch(r"(s|b)_k(\d{2,3})(?:_q([1-9]\d*)r([1-9]\d*)p([1-9]\d*))?", env["DEF_DAY_DE_UPRN"])
        if not (m and 50 <= int(m.group(2)) <= 300):
            problems.append(f"de_uprn {env['DEF_DAY_DE_UPRN']!r} must be <s|b>_k<K*100>[_q<n>r<n>p<n>] with K in [0.5, 3]")
    if "DEF_DAY_MDE_SIGN_KEYGATE" in env:
        if "DEF_DAY_MDE_SIGN" not in env:
            problems.append("mde_sign_keygate gates the hook's sign: it needs mde_sign")
        try:
            q = float(env["DEF_DAY_MDE_SIGN_KEYGATE"])
        except ValueError:
            q = -1.0
        if not 0.0 < q < 1.0:
            problems.append(f"mde_sign_keygate {env['DEF_DAY_MDE_SIGN_KEYGATE']!r} must be a quantile 0 < q < 1")
    if "DEF_DAY_DE_DNMAG" in env:
        if "DEF_DAY_DE_MAG" not in env:
            problems.append("de_dnmag runs after the zs pass: it needs de_mag")
        m = re.fullmatch(r"d_m(\d{3})(?:_q([1-9]\d*))?", env["DEF_DAY_DE_DNMAG"])
        if not (m and 125 <= int(m.group(1)) <= 400):
            problems.append(f"de_dnmag {env['DEF_DAY_DE_DNMAG']!r} must be d_m<m*100>[_q<n>] with m in [1.25, 4]")


def _problems_690g(env: dict, problems: list) -> None:
    """690g: the dump hook must name an existing script in the repo (def_day.sh also requires it to be committed)."""
    if "DEF_DAY_DUMP_HOOK" in env and not (ROOT / env["DEF_DAY_DUMP_HOOK"]).is_file():
        problems.append(f"dump_hook {env['DEF_DAY_DUMP_HOOK']!r} is not a file in the repo")


def installed_problems(env: dict, flags: str) -> list:
    dd, bu = DEF_DAY.read_text(), BUILDER.read_text()
    sameline = SAMELINE.read_text() if SAMELINE.is_file() else ""
    # The export, not a comment that mentions the knob (519): a substring match let a comment pass.
    # 517 rewrites the MDE export to $MDS (clean4 keeps the request; hookbank becomes file:<npz>:<sha>:...).
    export_of = {
        "DEF_DAY_DE_MAG": 'export VCC_DE_MAG="$DEF_DAY_DE_MAG"',
        "DEF_DAY_MDE_SIGN": 'export VCC_DE_MDE_SIGN="$MDS"',
        "DEF_DAY_RC_CALLSET_CPM": 'export VCC_RC_CALLSET_CPM="$DEF_DAY_RC_CALLSET_CPM"',
        "DEF_DAY_RHO_SCALE": '--rho-scale "$DEF_DAY_RHO_SCALE"',
        "DEF_DAY_CD4_STACK": 'VCC_CD4_TARGET_G="$CD4G"',
        "DEF_DAY_LODO_STRENGTH": "def_strength_key522.py",
        "DEF_DAY_HOOK_ESTIMATOR": '--programme-estimator "$DEF_DAY_HOOK_ESTIMATOR"',
        "DEF_DAY_LODO_HI": '--hi "$DEF_DAY_LODO_HI"',
        "DEF_DAY_CALLSET_BLOCK": 'export VCC_RC_CALLSET_BLOCK="$CSB"',
    }
    out = [f"scripts/def_day.sh has no {var}" for var in env if export_of.get(var, var) not in dd]
    if "DEF_DAY_CD4_STACK" in env and 'VCC_CD4_TARGET_G="$CD4G"' not in dd:
        out.append("cd4_stack x1.5: scripts/def_day.sh has no DEF_DAY_CD4_STACK dump override (530)")
    out += [f"the builder lacks {need} (needed by {var})" for var, need in BUILDER_NEEDS.items() if var in env and need not in bu]
    # 598: tercilec needs the builder's dump_key_carried and def_day.sh's tercilec case
    if env.get("DEF_DAY_SHARPEN_PT", "").startswith("tercilec:"):
        if "def dump_key_carried" not in bu:
            out.append("the builder lacks dump_key_carried (needed by DEF_DAY_SHARPEN_PT=tercilec, 598)")
        if "tercilec:*,*,*" not in dd:
            out.append("scripts/def_day.sh does not accept DEF_DAY_SHARPEN_PT=tercilec (598)")
    if env.get("DEF_DAY_MDE_SIGN", "").startswith("hookbank:") and "hookbank:" not in dd:
        out.append("scripts/def_day.sh has no hookbank form (517)")
    if env.get("DEF_DAY_MDE_SIGN", "").startswith("hookbank:") and "file:<npz>" not in bu:
        out.append("the builder lacks file:<npz> (needed by hookbank)")
    if "DEF_DAY_HOOK_ESTIMATOR" in env and "--programme-estimator" not in sameline:
        out.append("scripts/def_sameline.py has no --programme-estimator (517)")
    if "DEF_DAY_CALLSET_BLOCK" in env:
        emitter = EMITTER.read_text() if EMITTER.is_file() else ""
        if "block_steer=" not in emitter:
            out.append("the emitter has no block_steer (528)")
        if not CALLBLOCK.is_file():
            out.append("experiments/528-callprop-block/callblock528.py is missing")
    if flags == "" and EMPTY_FLAGS_MARK not in dd:
        out.append("scripts/def_day.sh treats an empty flags argument as pt:noprog (needs >= e29e793a, 497)")
    return out


def cmd_show(d: dict) -> int:
    print(f"october-decisions.json (updated {d.get('updated')})")
    for k, v in d["switches"].items():
        val = v.get("value")
        extra = v.get("evidence") if v.get("status") == "DECIDED" else v.get("decided_by", v.get("evidence", ""))
        print(f"  {k:17s} {v.get('status'):8s} {str(val):28s} {str(extra)[:110]}")
    return 0


def cmd_command(d: dict, a) -> int:
    env, args, pending, problems = resolve(d, a.allow_pending)
    problems += installed_problems(env, args[2])
    if a.json:
        print(json.dumps({"env": env, "args": [a.bundle, a.work, a.name] + args, "pending_treated_off": pending if a.allow_pending else [],
                          "problems": problems}, indent=1))
    else:
        line = " ".join(f"{k}={shlex.quote(v)}" for k, v in env.items())
        line += (" " if line else "") + "bash scripts/def_day.sh " + " ".join(shlex.quote(x) for x in [a.bundle, a.work, a.name] + args)
        if pending and a.allow_pending:
            print(f"# PENDING switches treated as OFF: {', '.join(pending)}")
        if d["switches"].get("cd4_thresholds", {}).get("status") == "DAY-OF":
            print("# day-of: if stage 1 stops with the CD4-gate message, check the mapping and rerun with DEF_DAY_CD4_THRESHOLDS=relaxed")
        # A nonzero exit still used to print a pasteable command. On today's def_day.sh that command
        # sets DEF_DAY_DE_MAG and nothing reads it (519).
        print(("# REFUSED " + line) if problems else line)
    for p in problems:
        print(f"PROBLEM: {p}", file=sys.stderr)
    if any(p.startswith("PENDING") for p in problems):
        return 2
    return 3 if problems else 0


def cmd_check(d: dict, a) -> int:
    env, args, pending, problems = resolve(d, a.allow_pending)
    problems += installed_problems(env, args[2])
    print(f"env {env}  args {args}  pending {pending}")
    for p in problems:
        print(f"PROBLEM: {p}")
    if any(p.startswith("PENDING") for p in problems):
        return 2
    return 3 if problems else 0


def cmd_decide(d: dict, a) -> int:
    s = d["switches"]
    if a.switch not in s:
        print(f"unknown switch {a.switch!r}; known: {', '.join(s)}", file=sys.stderr)
        return 1
    if not a.evidence.strip():
        print("--evidence is required (the read, the entry id, the file)", file=sys.stderr)
        return 1
    sw = s[a.switch]
    old = (sw.get("status"), sw.get("value"))
    sw.update({"status": "DECIDED", "value": a.value, "at": now(), "evidence": a.evidence})
    sw.setdefault("history", []).append({"from": list(old), "to": a.value, "at": sw["at"]})
    d["updated"] = sw["at"]
    tmp = DECISIONS.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(d, indent=1) + "\n")
    tmp.replace(DECISIONS)
    print(f"{a.switch}: {old} -> DECIDED {a.value!r} at {sw['at']}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("show")
    c = sub.add_parser("command")
    c.add_argument("--bundle", required=True)
    c.add_argument("--work", required=True)
    c.add_argument("--name", required=True)
    c.add_argument("--allow-pending", action="store_true")
    c.add_argument("--json", action="store_true")
    k = sub.add_parser("check")
    k.add_argument("--allow-pending", action="store_true")
    e = sub.add_parser("decide")
    e.add_argument("switch")
    e.add_argument("value")
    e.add_argument("--evidence", required=True)
    a = ap.parse_args()
    d = load()
    return {"show": lambda: cmd_show(d), "command": lambda: cmd_command(d, a), "check": lambda: cmd_check(d, a),
            "decide": lambda: cmd_decide(d, a)}[a.cmd]()


if __name__ == "__main__":
    sys.exit(main())
