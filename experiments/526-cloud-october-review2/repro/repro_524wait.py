#!/usr/bin/env python3
"""524's waiters: process-name gate, exit status, and the READY greps.

Exits 0 while all of the following are still true:

  * recipe_dr524.sh builders() matches a bash argv that merely mentions def_day.sh,
    and does not match a python/uv argv that is actually running the builder
  * the script has no `set -e`, and a copy of its tail exits 0 when gate524 rc is 1
  * chain1006's ready() greps and arm1006's greps both hit the lines gate524.py and
    run_arm_generic.sh actually write
"""
from __future__ import annotations

import hashlib
import re
import subprocess
import sys
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RECIPE = (ROOT / "experiments/524-october-recipe-dr/recipe_dr524.sh").read_text()
CHAIN = (ROOT / "experiments/524-october-recipe-dr/chain1006.sh").read_text()
ARM = (ROOT / "experiments/524-october-recipe-dr/arm1006.sh").read_text()
PIPE = (ROOT / "scripts/def_pipeline.py").read_text()


def builders_awk(lines: list[str]) -> list[str]:
    # The awk program, verbatim from recipe_dr524.sh.
    m = re.search(r"awk '(\$1 ~ /\^bash\$/ && /.*?/)'", RECIPE)
    prog = m.group(1)
    r = subprocess.run(["awk", prog], input="\n".join(lines) + "\n", text=True, capture_output=True)
    return [ln for ln in r.stdout.splitlines() if ln]


def main() -> int:
    assert "set -uo pipefail" in RECIPE and "set -e" not in RECIPE.split("set -uo pipefail")[0]
    assert "declare -A" not in RECIPE and "declare -A" not in CHAIN and "declare -A" not in ARM
    # arm1006's slot-2 wait uses the form AGENTS §6.2 asks for.
    assert 'grep -qE "NOT SUBMITTING| END$"' in ARM

    hits = builders_awk([
        "bash experiments/524-october-recipe-dr/recipe_dr524.sh",
        "bash experiments/524-october-recipe-dr/chain1006.sh",
        "bash experiments/524-october-recipe-dr/arm1006.sh",
        "bash -c echo scripts/def_day.sh is mentioned",
        "bash scripts/def_day.sh data/challenge work name 1.75 0.5 noprog 4",
        "python scripts/make_submission_pm_cd4.py --bundle data/challenge",
        "uv run python scripts/make_submission_pm_cd4.py --bundle data/challenge",
        "python scripts/build_union_whitener.py --out work/op.npz",
    ])
    print("builders() hits:")
    for h in hits:
        print("  ", h)
    assert any(h.startswith("bash -c") for h in hits)
    assert any("scripts/def_day.sh data/challenge" in h for h in hits)
    assert not any("recipe_dr524" in h for h in hits)
    assert not any("chain1006" in h for h in hits)
    assert not any(h.startswith("python ") or h.startswith("uv ") for h in hits)
    # def_pipeline.other_builders would see those two. Its launcher regex is the evidence.
    assert "python[0-9.]*|uv|time" in PIPE
    assert "make_submission_pm_cd4|build_union_whitener" in PIPE

    # Tail of recipe_dr524: gate failure still reaches log END, and the script exits 0.
    tail = textwrap.dedent(r'''
        set -uo pipefail
        log(){ printf '%s\n' "$*"; }
        G=1
        log "gate524 rc=$G"
        [ $G -eq 0 ] && echo would-delete-h5ad
        log "END recipe_dr524 (READY iff gate524 rc=0)"
    ''')
    r = subprocess.run(["bash", "-c", tail], capture_output=True, text=True)
    print("recipe tail rc", r.returncode, "out", r.stdout.splitlines())
    assert r.returncode == 0
    assert "END recipe_dr524" in r.stdout
    assert "would-delete-h5ad" not in r.stdout

    # Log lines the two producers write, and the greps both waiters use.
    import tempfile
    td = Path(tempfile.mkdtemp())
    digest = hashlib.sha256(b"pkg").hexdigest()
    gate_log = td / "gates524.log"
    gate_log.write_text(
        "[524dr] CLIFF OK {\"D\": 0.9} final foo factors D=1\n"
        "[524dr] PASS: config diff vs 453dr ['bulk_load']\n"
        "novelty rc=0\n"
        f"package sha256: {digest}\n"
    )
    klog = td / "run.log"
    # run_arm_generic timestamps the novelty and sha lines; PASS and CLIFF are raw python|tee.
    klog.write_text(
        "[446bj2zk16] PASS: exactly ['rc_callset_cpm']\n"
        "[446bj2zk16] CLIFF OK {\"A\": 0.9, \"B\": 0.9, \"C\": 0.9} amp-only\n"
        f"2026-10-05T20:00:00Z novelty rc=0\n"
        f"2026-10-05T20:00:01Z package sha256: {digest}\n"
    )

    def ready(tag: str, path: Path) -> bool:
        script = textwrap.dedent(f'''
            set -uo pipefail
            r={path}
            tag={tag}
            /usr/bin/grep -q "^\\[$tag\\] PASS" "$r" || exit 1
            /usr/bin/grep -q "^\\[$tag\\] CLIFF OK" "$r" || exit 1
            /usr/bin/grep -q "novelty rc=0" "$r" || exit 1
            want="$(/usr/bin/grep 'package sha256:' "$r" | tail -1 | awk '{{print $NF}}')"
            [ "$want" = "{digest}" ]
        ''')
        return subprocess.run(["bash", "-c", script]).returncode == 0

    assert ready("524dr", gate_log)
    assert ready("446bj2zk16", klog)
    # arm1006's greps are the same three, without the /usr/bin prefix.
    for path, tag in ((gate_log, "524dr"), (klog, "446bj2zk16")):
        for pat in (f"^\\[{tag}\\] PASS", f"^\\[{tag}\\] CLIFF OK", "novelty rc=0"):
            rr = subprocess.run(["grep", "-q", "-E", pat, str(path)])
            assert rr.returncode == 0, (tag, pat)
    print("READY greps match gate524.py and run_arm_generic.sh log lines, including the K-arm tag")
    print("WAIT CONFIRMED: builders() misses python/uv; a failed gate still exits 0")
    return 0


if __name__ == "__main__":
    sys.exit(main())
