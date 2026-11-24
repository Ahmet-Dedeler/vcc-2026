"""591 — UNAPPLIED patch: a pressure-based alternative for run_arm_generic.sh's PREP gate.

Problem (coordinator, 2026-09-30; experiments/549-cph-arm/RESULTS.md): `gate_ok` requires vm_stat free+inactive ≥ 23 GiB.
With ~18 GiB of other apps' pages in the compressor, free+inactive stays under 23 GiB while memory_pressure reports ~60 %
free and no swap is used, so the gate never opens, `wait_gate` gives up after 2 h (exit 3), and 549 had to be prepped
by hand.

Change (the coordinator's spec, plus one restriction found while reading the callers):
- `gate_ok <mode>` opens on either memory branch, with every other condition unchanged (no other builder, no `vcc prep`,
  disk ≥ 20 GiB):
  - `fi23`: memory_pressure free ≥ 40 % AND free+inactive ≥ 23 GiB. This is today's gate, used in both modes.
  - `pressure55`: memory_pressure free ≥ 55 % AND swap free ≥ 400 MiB. **Prep only.**
- **The build gate keeps `fi23` alone.** Invocation 1 runs `def_pipeline.py --min-free-gib 23`. Before its
  repin-panel-support and build stages, that script re-checks free+inactive+speculative+purgeable ≥ 23 GiB and **refuses
  (SystemExit)** instead of waiting. Opening the shell gate on pressure there would turn today's wait into a failed
  build. Extending the branch to the build needs the same alternative inside `def_pipeline.py`. That is a pinned driver
  (run_arm_generic's `case` list), so it is not done here.
- The gate line keeps its `gate: free …` prefix and field order and appends `swap free … MiB  mode …  memory branch …`.
- `wait_gate` passes the mode; the two call sites become `wait_gate build` (before invocation 1) and `wait_gate prep`
  (before `vcc prep`).

    uv run python experiments/591-prep-gate/patch_gate591.py --check     # writes nothing; bash -n on the result
    uv run python experiments/591-prep-gate/patch_gate591.py --apply     # atomic rename; refuses while a bash runs it
    bash experiments/591-prep-gate/test_gate591.sh                       # the stubbed-reader tests (any time)
"""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TARGET = ROOT / "experiments/428-zero-gene-audit/run_arm_generic.sh"

OLD_GATE = '''gate_ok(){
  local pct fi others preps disk
  pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\\./,"",f);gsub(/\\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  others="$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null || echo 99)"
  preps="$(ps -axo args= | grep -vE 'claude|grep' | grep -cE '(^|/| )vcc prep( |$)' | tr -d ' ')"
  disk="$(free_gib)"
  log "gate: free ${pct}%  free+inactive ${fi} GiB  other builders ${others}  vcc preps ${preps}  disk ${disk} GiB"
  [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=23)}' && [ "$others" = "0" ] \\
    && [ "$preps" = "0" ] && awk -v v="$disk" 'BEGIN{exit !(v>=20)}'
}
wait_gate(){ local n=0; until gate_ok; do n=$((n+1)); [ $n -ge 240 ] && { log "GATE never opened"; exit 3; }; sleep 30; done; }
'''

NEW_GATE = '''gate_ok(){
  # $1 = build | prep. Both need no other builder, no vcc prep, disk >= 20 GiB, and one memory branch:
  #   fi23        memory_pressure free >= 40 % AND vm_stat free+inactive >= 23 GiB (the original gate; both modes)
  #   pressure55  memory_pressure free >= 55 % AND swap free >= 400 MiB (591; PREP ONLY: other apps' pages held in the
  #               compressor keep free+inactive under 23 GiB while the machine has room; 549 was prepped by hand 09-30)
  # The build keeps fi23 alone: def_pipeline.py re-checks free+inactive(+speculative+purgeable) >= --min-free-gib before
  # its build stages and REFUSES instead of waiting, so a looser gate there would turn a wait into a failed build.
  local mode="${1:-build}" pct fi swapfree others preps disk branch
  pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\\./,"",f);gsub(/\\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  swapfree="$(sysctl -n vm.swapusage 2>/dev/null | awk '{for(i=1;i<NF;i++) if($i=="free"){v=$(i+2); u=substr(v,length(v),1); n=substr(v,1,length(v)-1)+0; if(u=="G") n=n*1024; else if(u=="K") n=n/1024; printf "%.0f", n; exit}}')"
  others="$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null || echo 99)"
  preps="$(ps -axo args= | grep -vE 'claude|grep' | grep -cE '(^|/| )vcc prep( |$)' | tr -d ' ')"
  disk="$(free_gib)"
  branch=none
  if [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=23)}'; then branch=fi23
  elif [ "$mode" = prep ] && [ "${pct:-0}" -ge 55 ] && awk -v v="${swapfree:-0}" 'BEGIN{exit !(v>=400)}'; then branch=pressure55
  fi
  log "gate: free ${pct}%  free+inactive ${fi} GiB  other builders ${others}  vcc preps ${preps}  disk ${disk} GiB  swap free ${swapfree:-?} MiB  mode ${mode}  memory branch ${branch}"
  [ "$branch" != none ] && [ "$others" = "0" ] && [ "$preps" = "0" ] && awk -v v="$disk" 'BEGIN{exit !(v>=20)}'
}
wait_gate(){ local n=0; until gate_ok "${1:-build}"; do n=$((n+1)); [ $n -ge 240 ] && { log "GATE never opened ($1)"; exit 3; }; sleep 30; done; }
'''

CALLS = [
    ('wait_gate\nlog "DRIVER invocation 1: preflight ... postbuild  (panel A/B/C, real challenge bundle)"\n',
     'wait_gate build\nlog "DRIVER invocation 1: preflight ... postbuild  (panel A/B/C, real challenge bundle)"\n'),
    ('wait_gate\nlog "INVOCATION 2: vcc prep -> $PKG"\n',
     'wait_gate prep\nlog "INVOCATION 2: vcc prep -> $PKG"\n'),
]


def patched(src: str) -> str:
    if "memory branch" in src:
        raise SystemExit("run_arm_generic.sh already carries 591's gate")
    for old, new in [(OLD_GATE, NEW_GATE), *CALLS]:
        n = src.count(old)
        if n != 1:
            raise SystemExit(f"anchor found {n} times (need 1):\n{old[:160]}")
        src = src.replace(old, new)
    return src


def live_runners() -> list[str]:
    out = subprocess.run(["ps", "-axo", "pid=,comm=,args="], capture_output=True, text=True).stdout
    rows = []
    for line in out.splitlines():
        parts = line.split(None, 2)
        if len(parts) == 3 and parts[1].rsplit("/", 1)[-1] == "bash" and "run_arm_generic.sh" in parts[2]:
            rows.append(line.strip())
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", type=Path, default=TARGET)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--check", action="store_true")
    g.add_argument("--apply", action="store_true")
    ap.add_argument("--force", action="store_true", help="apply even while a bash runs the target (atomic rename)")
    a = ap.parse_args()
    src = a.target.read_text()
    new = patched(src)
    with tempfile.NamedTemporaryFile("w", suffix=".sh", delete=False) as fh:
        fh.write(new)
    r = subprocess.run(["/bin/bash", "-n", fh.name], capture_output=True, text=True)
    os.unlink(fh.name)
    if r.returncode != 0:
        raise SystemExit(f"bash -n failed on the patched script:\n{r.stderr}")
    old_sha, new_sha = hashlib.sha256(src.encode()).hexdigest(), hashlib.sha256(new.encode()).hexdigest()
    print(f"[591] {a.target.name} {old_sha[:16]} -> {new_sha[:16]} (anchors OK, bash -n OK)")
    if a.check:
        print("[591] --check: nothing written")
        return
    runners = live_runners()
    if runners and not a.force:
        raise SystemExit("[591] a bash process is executing run_arm_generic.sh; apply after it exits (or --force):\n  "
                         + "\n  ".join(runners))
    tmp = a.target.with_name(a.target.name + ".591.tmp")
    tmp.write_text(new)
    os.chmod(tmp, os.stat(a.target).st_mode)
    os.replace(tmp, a.target)  # atomic: a running bash keeps reading the old inode
    print(f"[591] applied: {a.target} sha256 {hashlib.sha256(a.target.read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
