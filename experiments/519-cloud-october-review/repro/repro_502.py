#!/usr/bin/env python3
"""502 pass-through: commit regex, (b) exit status, half-write, riders_match, export order, env merge.

Touches nothing under scripts/ or experiments/502. Helpers are imported from the repro copy.
"""
from __future__ import annotations

import importlib.machinery
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
EXP = REPO / "experiments/502-october-de-riders"
bad = 0


def expect(label, cond, detail=""):
    global bad
    print(f"{'OK' if cond else 'UNEXPECTED'}  {label}" + (f": {detail}" if detail != "" else ""))
    if not cond:
        bad += 1


def load_helpers():
    # .new is not a registered suffix, so spec_from_file_location returns None.
    path = EXP / "repro/def_day_helpers.py.new"
    loader = importlib.machinery.SourceFileLoader("helpers502", str(path))
    spec = importlib.util.spec_from_loader("helpers502", loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod


def spec_reqs() -> list[str]:
    """Same split as scripts/commit_after_chain.sh: read -r, value is everything after the first =."""
    reqs = []
    for line in (EXP / "commit502.spec").read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        k, _, v = line.partition("=")
        if k == "REQ":
            reqs.append(v)
    return reqs


def grep_e(pattern: str, text: str) -> int:
    r = subprocess.run(["grep", "-E", "-c", pattern], input=text, text=True, capture_output=True)
    if r.returncode == 1:
        return 0
    r.check_returncode()
    return int(r.stdout.strip())


def main() -> int:
    new = (EXP / "repro/def_day.sh.new").read_text().splitlines()
    # exports sit after the dump subshell closes
    dump_close = next(i for i, l in enumerate(new) if l.strip() == ")" and i > 190)
    export_mag = next(i for i, l in enumerate(new) if 'export VCC_DE_MAG="$DEF_DAY_DE_MAG"' in l)
    export_bulk = next(i for i, l in enumerate(new) if l.startswith("export VCC_BULK_LOAD="))
    wipe = next(i for i, l in enumerate(new) if l.startswith("for v in $(env | grep -o"))
    expect("wipe is before the dump", wipe < dump_close, f"wipe L{wipe+1} close L{dump_close+1}")
    expect("VCC_DE_MAG export is after the dump subshell", export_mag > dump_close,
           f"export L{export_mag+1} close L{dump_close+1}")
    expect("bulk load export is the line before the rider export", export_bulk + 1 == export_mag,
           f"bulk L{export_bulk+1}")

    # Control flow copied from those lines, run as bash. The dump subshell is forked before the export.
    sim = textwrap.dedent(r"""
        set -uo pipefail
        export VCC_DE_MAG=PRESET_SHOULD_DIE
        export DEF_DAY_DE_MAG=dn_k250_zs
        for v in $(env | grep -o '^VCC_[A-Z0-9_]*'); do unset "$v"; done
        (
          if [ -n "${VCC_DE_MAG:-}" ]; then echo "DUMP_SAW=$VCC_DE_MAG"; else echo DUMP_SAW=; fi
        )
        export VCC_BULK_LOAD=dump:sha
        [ -n "${DEF_DAY_DE_MAG:-}" ] && export VCC_DE_MAG="$DEF_DAY_DE_MAG"
        echo "LEVEL_SAW=${VCC_DE_MAG:-}"
        echo "BULK_SAW=${VCC_BULK_LOAD:-}"
    """)
    r = subprocess.run(["bash", "-c", sim], text=True, capture_output=True, check=True)
    expect("dump does not see VCC_DE_MAG", "DUMP_SAW=\n" in r.stdout or r.stdout.startswith("DUMP_SAW=\n"), r.stdout.strip())
    expect("level build sees the rider", "LEVEL_SAW=dn_k250_zs" in r.stdout, r.stdout.strip())

    # Driver merge: VCC_DE_MAG is not a shipped default, so {**os.environ, **resolved} keeps it.
    # bare_env drops every VCC_* before the overlay, and the overlay does not put it back.
    sys.path.insert(0, str(REPO / "scripts"))
    import def_pipeline as dp
    expect("VCC_DE_MAG not in SHIPPED_ENV or DEF_ENV",
           "VCC_DE_MAG" not in dp.SHIPPED_ENV and "VCC_DE_MAG" not in dp.DEF_ENV, "")
    saved = os.environ.get("VCC_DE_MAG")
    os.environ["VCC_DE_MAG"] = "dn_k250_zs"
    try:
        resolved, _ = dp.resolve_env(dp.day_env(False))
        base = dict(os.environ)
        child = {**base, **resolved}
        bare_base = {k: v for k, v in os.environ.items() if not k.startswith("VCC_")}
        bare_child = {**bare_base, **resolved}
    finally:
        if saved is None:
            os.environ.pop("VCC_DE_MAG", None)
        else:
            os.environ["VCC_DE_MAG"] = saved
    expect("level-build child keeps VCC_DE_MAG", child.get("VCC_DE_MAG") == "dn_k250_zs", child.get("VCC_DE_MAG"))
    expect("bare_env child drops VCC_DE_MAG", "VCC_DE_MAG" not in bare_child, bare_child.get("VCC_DE_MAG"))

    # ALLOW and final.env lines are conditional on the variable being set.
    src = "\n".join(new)
    expect("ALLOW de_mag is conditional", '[ -n "${VCC_DE_MAG:-}" ] && ALLOW="$ALLOW,de_mag"' in src, "")
    expect("final.env de_mag is conditional",
           '[ -n "${VCC_DE_MAG:-}" ] && printf \'VCC_DE_MAG=%s\\n\' "$VCC_DE_MAG" >> "$W/def_day.final.env"' in src, "")

    # commit502.spec as grep -E sees it, against the lines chain502.sh actually logs.
    log = "\n".join([
        "2026-09-28T00:00:00Z applied: def_day 610db78277239f50, helpers 20142142e4574b7f; test_riders502 9/9",
        "2026-09-28T00:00:01Z (a) PASS: the 502 script with neither rider retraces efp relaxed byte for byte (ac49d319f5982507)",
        "2026-09-28T00:00:02Z (b) PASS: rider on -> DONE, postbuild PASS, riders_match, meta de_mag dn_k250_zs, output abcdef != (a)",
        "2026-09-28T00:00:03Z END def_day 610db78277239f50, helpers 20142142e4574b7f",
    ])
    for req in spec_reqs():
        expect(f"REQ matches log: {req}", grep_e(req, log) == 1, grep_e(req, log))
    # the unescaped form is what would have missed; the file does not contain it
    expect("unescaped (a) PASS does not match", grep_e("(a) PASS", log) == 0, "")
    raw = (EXP / "commit502.spec").read_text()
    expect("spec stores the escaped form", r"REQ=\(a\) PASS" in raw and r"REQ=\(b\) PASS" in raw, "")

    # (b) failure falls off the end. Reproduce the tail with set -uo pipefail and no exit in the else.
    tail = textwrap.dedent(r"""
        set -uo pipefail
        GOTA=aaa
        V=nope
        GOTB=aaa
        if false; then
          if [ "$V" = dn_k250_zs ] && [ "$GOTB" != "$GOTA" ]; then echo PASS; else echo "(b) FAIL inner"; fi
        else echo "(b) FAIL outer"; fi
        echo END
    """)
    r = subprocess.run(["bash", "-c", tail], text=True, capture_output=True)
    expect("(b) FAIL tail exits 0", r.returncode == 0 and "(b) FAIL outer" in r.stdout and r.stdout.strip().endswith("END"),
           f"rc={r.returncode} out={r.stdout.strip()}")

    # Half-write: first file replaced, second not, then the release predicate (old sha prefix) is false
    # and nothing restores. chain502.sh:23 requires def_day prefix e29e793a7964df03 before it will run.
    pre = (REPO / "scripts/def_day.sh").read_text()
    post = (EXP / "repro/def_day.sh.new").read_text()
    with tempfile.TemporaryDirectory() as td:
        dd = Path(td) / "def_day.sh"
        hp = Path(td) / "helpers.py"
        dd.write_text(post)          # killed after the first write_text
        hp.write_text("STILL_OLD\n")
        import hashlib
        sha = hashlib.sha256(dd.read_bytes()).hexdigest()[:16]
        released = sha == "e29e793a7964df03"
        expect("half-write is not the pre-502 sha, so a re-arm does not release", not released, sha)
        expect("helpers were not updated", hp.read_text() == "STILL_OLD\n", "")
        expect("pre-502 prefix of the live script is the gate's pin",
               hashlib.sha256(pre.encode()).hexdigest().startswith("e29e793a7964df03"), "")
        expect("repro def_day prefix is commit502's EXPECT",
               hashlib.sha256(post.encode()).hexdigest().startswith("610db78277239f50"), "")
        hp_new = (EXP / "repro/def_day_helpers.py.new").read_bytes()
        expect("repro helpers prefix is commit502's EXPECT",
               hashlib.sha256(hp_new).hexdigest().startswith("20142142e4574b7f"), "")

    # riders_match. ON path fails closed for a missing/empty/non-dict de_mag.
    # OFF path treats a non-dict and an empty dict as absent (fail-open on that path only).
    helpers = load_helpers()
    with tempfile.TemporaryDirectory() as td:
        def meta(cfg):
            p = Path(td) / "m.json"
            p.write_text(json.dumps({"config": cfg}))
            return str(p)

        def call(cfg, mag, sign):
            old = sys.argv
            sys.argv = ["helpers", "riders_match", meta(cfg), mag, sign]
            try:
                return helpers.main()
            finally:
                sys.argv = old

        good = {"de_mag": {"variant": "dn_k250_zs", "module_sha256": "abc"}}
        expect("ON dict matches", call(good, "dn_k250_zs", "-") == 0, "")
        expect("ON missing fails", call({}, "dn_k250_zs", "-") == 1, "")
        expect("ON empty dict fails", call({"de_mag": {}}, "dn_k250_zs", "-") == 1, "")
        expect("ON non-dict fails", call({"de_mag": "dn_k250_zs"}, "dn_k250_zs", "-") == 1, "")
        expect("OFF empty config matches", call({}, "-", "-") == 0, "")
        off_empty = call({"de_mag": {}}, "-", "-")
        off_nonstr = call({"de_mag": "garbage"}, "-", "-")
        expect("OFF empty dict is treated as absent (fail-open)", off_empty == 0, off_empty)
        expect("OFF non-dict is treated as absent (fail-open)", off_nonstr == 0, off_nonstr)
        expect("OFF rejects a real variant", call(good, "-", "-") == 1, "")
        sign_ok = call({"de_mde_sign": {"source": "clean4", "tau": 0.05}}, "-", "clean4:0.05")
        expect("sign dict matches", sign_ok == 0, sign_ok)
        try:
            rc_tau = call({"de_mde_sign": {"source": "clean4", "tau": "foo"}}, "-", "clean4:0.05")
        except ValueError:
            rc_tau = "ValueError"
        expect("non-numeric tau is loud, not a pass", rc_tau != 0, rc_tau)

    # Release awk: interpreter + path, and it does not match chain502.sh itself.
    awk = r"""awk '$1 ~ /^bash$/ && /chain4(7[69]b?|9[5-8])\.sh|chain496p\.sh|chain497e\.sh|chain51[01]\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4/'"""
    ps = "\n".join([
        "bash experiments/502-october-de-riders/chain502.sh",
        "bash scripts/def_day.sh data/challenge work stem",
        "python3 scripts/make_submission_pm_cd4.py",
        "bash experiments/511-tercile-ties/chain511.sh",
    ])
    r = subprocess.run(f"printf '%s\n' \"$PS\" | {awk}", shell=True, text=True,
                       capture_output=True, env={**os.environ, "PS": ps})
    matched = r.stdout.splitlines()
    expect("awk matches def_day.sh and chain511.sh",
           any("def_day.sh" in m for m in matched) and any("chain511.sh" in m for m in matched), matched)
    expect("awk does not match chain502.sh or the python builder",
           not any("chain502" in m or "make_submission" in m for m in matched), matched)

    print("FAILS" if bad else "ALL DEMONSTRATED", bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
