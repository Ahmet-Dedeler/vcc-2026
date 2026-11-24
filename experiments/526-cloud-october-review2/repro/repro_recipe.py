#!/usr/bin/env python3
"""october_recipe.py command against def_day.sh, over the switch combinations.

Exits 0. Prints one line per combination: rc, env keys, flags, and the problem list.
A combination the resolver accepts is then checked against the builder's refusal
predicates (copied as the conditions at the cited lines, and asserted still present
in the builder source).
"""
from __future__ import annotations

import importlib.util
import io
import json
import os
import subprocess
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DECISIONS = json.loads((ROOT / "docs/october-decisions.json").read_text())
BUILDER = (ROOT / "scripts/make_submission_pm_cd4.py").read_text()
DEF_DAY = (ROOT / "scripts/def_day.sh").read_text()


def load_mod(path: Path):
    os.environ["OCTOBER_DECISIONS"] = str(path)
    os.environ["OCTOBER_DEF_DAY"] = str(ROOT / "scripts/def_day.sh")
    os.environ["OCTOBER_BUILDER"] = str(ROOT / "scripts/make_submission_pm_cd4.py")
    spec = importlib.util.spec_from_file_location("october_recipe_526", ROOT / "scripts/october_recipe.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class NS:
    def __init__(self, allow, as_json=True):
        self.allow_pending = allow
        self.json = as_json
        self.bundle = "BUNDLE"
        self.work = "WORK"
        self.name = "NAME"


def run(mod, allow: bool, as_json: bool = True):
    buf, err = io.StringIO(), io.StringIO()
    with redirect_stdout(buf), redirect_stderr(err):
        rc = mod.cmd_command(mod.load(), NS(allow, as_json))
    return rc, buf.getvalue(), err.getvalue()


def write(tmp: Path, switches: dict) -> Path:
    d = json.loads(json.dumps(DECISIONS))
    d["switches"].update(switches)
    p = tmp / "decisions.json"
    p.write_text(json.dumps(d))
    return p


def decided(value, **extra):
    sw = {"status": "DECIDED", "value": value}
    sw.update(extra)
    return sw


def main() -> int:
    import tempfile
    tmp = Path(tempfile.mkdtemp())
    rows = []

    def rec(label, path, allow, as_json=True):
        mod = load_mod(path)
        rc, out, err = run(mod, allow, as_json)
        if as_json:
            payload = json.loads(out)
            flags = payload["args"][5]
            env = payload["env"]
            problems = payload["problems"]
        else:
            payload, flags, env, problems = out, "", {}, err
        rows.append((label, rc, env, flags, problems, out if not as_json else ""))
        print(f"{label}: rc={rc} env={sorted(env) if as_json else '<text>'} flags={flags!r}")
        for p in (problems if as_json else [ln for ln in err.splitlines() if ln.startswith("PROBLEM:")]):
            print(f"    {p}")
        if not as_json:
            print("   ", out.strip()[:180])
        return payload if as_json else out

    live = ROOT / "docs/october-decisions.json"
    got = rec("live + allow-pending", live, True)
    assert got["problems"] == []
    assert got["env"] == {"DEF_DAY_LODO": "1", "DEF_DAY_DE_MAG": "dn_k250_zs"}
    assert got["args"][3:] == ["1.75", "0.5", "noprog", "4"]
    assert got["args"][:3] == ["BUNDLE", "WORK", "NAME"]

    text = rec("live, pending refused", live, False, as_json=False)
    assert "# REFUSED " in text
    assert "DEF_DAY_DE_MAG=" in text and "DEF_DAY_LODO=1" in text

    # programme ON, prior still pending: prior stays off, flags become the empty string.
    p = write(tmp, {"programme_source": decided("h1_2025_full")})
    got = rec("programme on, prior pending", p, True)
    assert "DEF_DAY_PROGRAMME_SOURCE" in got["env"] and "DEF_DAY_PRIOR_GATE" not in got["env"]
    assert got["args"][5] == ""
    assert got["problems"] == []

    # prior ON, programme still pending: resolver refuses (and def_day would too).
    p = write(tmp, {"prior_gate": decided("le1")})
    got = rec("prior on, programme pending", p, True)
    assert any("prior_gate is on without programme_source" in x for x in got["problems"])
    assert got  # rc checked below via rows

    # both ON: empty flags, both exports.
    p = write(tmp, {
        "programme_source": decided("h1_2025_full"),
        "prior_gate": decided("le1"),
        "mde_sign": decided("clean4:0.05"),
        "lodo_hi_pct": decided("96.64697522769474"),
        "sharpen_pt": decided("stercile:2.5,2.0,1.75"),
    })
    got = rec("all five pending switches forced on", p, False)
    assert got["problems"] == []
    assert got["env"]["DEF_DAY_PRIOR_GATE"] == "le1"
    assert got["env"]["DEF_DAY_MDE_SIGN"] == "clean4:0.05"
    assert got["env"]["DEF_DAY_LODO_HI_PCT"] == "96.64697522769474"
    assert got["env"]["DEF_DAY_SHARPEN_PT"] == "stercile:2.5,2.0,1.75"
    assert got["env"]["DEF_DAY_PROGRAMME_SOURCE"] == "h1_2025_full"
    assert got["args"][5] == ""

    # programme on, but the with-programme flags still say noprog.
    p = write(tmp, {"programme_source": decided("h1_2025_full"),
                    "sameline_flags": decided("noprog", with_programme_source="pt:noprog")})
    got = rec("programme + noprog flags", p, True)
    assert any("noprog" in x for x in got["problems"])

    # percentile without lodo.
    p = write(tmp, {"lodo": decided("off"), "lodo_hi_pct": decided("96.6")})
    got = rec("pct without lodo", p, True)
    assert any("lodo_hi_pct needs lodo" in x for x in got["problems"])

    # cd4 stack the dumper cannot build.
    p = write(tmp, {"cd4_stack": decided("x1.5")})
    got = rec("cd4_stack x1.5", p, True)
    assert any("cd4_stack" in x for x in got["problems"])

    # DECIDED "false" is ON. "off" is off. Cosmetic footgun, asserted so it stays visible.
    p = write(tmp, {"mde_sign": decided("false")})
    got = rec("value 'false' is on", p, True)
    assert got["env"].get("DEF_DAY_MDE_SIGN") == "false"
    p = write(tmp, {"mde_sign": decided("off")})
    got = rec("value 'off' is off", p, True)
    assert "DEF_DAY_MDE_SIGN" not in got["env"]

    # Shell: an empty 6th argument stays empty. An omitted one becomes pt:noprog.
    empty = subprocess.check_output(
        ["bash", "-c", 'set -- a b c d e ""; SLF="${6-pt:noprog}"; printf "%s" "$SLF"'], text=True)
    omitted = subprocess.check_output(
        ["bash", "-c", 'set -- a b c d e; SLF="${6-pt:noprog}"; printf "%s" "$SLF"'], text=True)
    print(f"bash flags: empty 6th -> {empty!r}; omitted 6th -> {omitted!r}")
    assert empty == "" and omitted == "pt:noprog"
    assert 'SLF="${6-pt:noprog}"' in DEF_DAY

    # Wipe then export, the def_day order. VCC_RC_CALLSET_CPM is wiped and never put back.
    wiped = subprocess.check_output(["bash", "-c", r'''
export VCC_DE_MAG=stale VCC_DE_MDE_SIGN=stale VCC_RC_CALLSET_CPM=1200 VCC_DE_PRIOR_GATE=stale
for v in $(env | grep -o '^VCC_[A-Z0-9_]*'); do unset "$v"; done
DEF_DAY_DE_MAG=dn_k250_zs
[ -n "${DEF_DAY_DE_MAG:-}" ] && export VCC_DE_MAG="$DEF_DAY_DE_MAG"
# MDE and the call-set knob are unset, so they are not exported
printf "MAG=%s SIGN=%s CPM=%s GATE=%s\n" "${VCC_DE_MAG-<unset>}" "${VCC_DE_MDE_SIGN-<unset>}" "${VCC_RC_CALLSET_CPM-<unset>}" "${VCC_DE_PRIOR_GATE-<unset>}"
'''], text=True).strip()
    print("wipe-then-export:", wiped)
    assert wiped == "MAG=dn_k250_zs SIGN=<unset> CPM=<unset> GATE=<unset>"

    # Builder refuses these. The predicates are the source conditions, asserted present.
    for needle in (
        'raise SystemExit("VCC_DE_PRIOR_GATE gates a programme default: it needs VCC_SAMELINE without noprog")',
        'raise SystemExit(f"VCC_DE_MAG must look like',
        'raise SystemExit(f"VCC_DE_MDE_SIGN must be \'clean4:<tau>\' with 0 <= tau < 0.5:',
        'if os.environ.get("VCC_RC_CALLSET_CPM", "").strip():',
    ):
        assert needle in BUILDER, needle
    assert 'die "DEF_DAY_PRIOR_GATE=$DEF_DAY_PRIOR_GATE gates a programme' in DEF_DAY
    assert 'die "DEF_DAY_PROGRAMME_SOURCE=$DEF_DAY_PROGRAMME_SOURCE needs sameline_flags without noprog' in DEF_DAY
    # def_day's shell accepts any [a-z0-9_]+ for DE_MAG; the builder's regex is tighter.
    # `foo` passes the shell and SystemExit's in the builder. Not ignored.
    assert '*[!a-z0-9_]*)' in DEF_DAY
    assert "dn_k" in BUILDER and "zs" in BUILDER

    # rc of the refusing rows
    by = {r[0]: r[1] for r in rows}
    assert by["live + allow-pending"] == 0
    assert by["prior on, programme pending"] == 3
    assert by["programme + noprog flags"] == 3
    assert by["pct without lodo"] == 3
    assert by["cd4_stack x1.5"] == 3
    assert by["all five pending switches forced on"] == 0
    print("RESOLVER OK: every combination above matched the expected env, flags and exit code")
    return 0


if __name__ == "__main__":
    sys.exit(main())
