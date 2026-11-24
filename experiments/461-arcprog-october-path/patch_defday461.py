"""461: add DEF_DAY_PROGRAMME_SOURCE to scripts/def_day.sh (default off). With it set, stage 2 passes
`--programme-source <bank>` to def_sameline.py, and the script refuses sameline_flags containing `noprog` (the builder would
ignore the programme: a silent no-op). Unset: the def_sameline command line is unchanged (the ${VAR:+...} expansion is
empty). Apply ONLY when no def_day.sh process is running (bash reads a running script incrementally).
usage: uv run python experiments/461-arcprog-october-path/patch_defday461.py [--check]"""
import subprocess
import sys
from pathlib import Path

p = Path(__file__).resolve().parents[2] / "scripts/def_day.sh"
s = p.read_text()
edits = [
    ("""#   DEF_DAY_LODO=1 (460, default off): the stack dump is scaled by 415's LODO-V1 key re-fitted on the day's panel
#    (446bj2's key on A/B/C) instead of 410's; see stage 2b.
""", """#   DEF_DAY_LODO=1 (460, default off): the stack dump is scaled by 415's LODO-V1 key re-fitted on the day's panel
#    (446bj2's key on A/B/C) instead of 410's; see stage 2b.
#   DEF_DAY_PROGRAMME_SOURCE=<bank> (461, default off; e.g. h1_2025_full = Arc's H1 from all public cells, 458d): every
#    context's block prior comes from that bank (def_sameline.py --programme-source). Needs sameline_flags WITHOUT noprog
#    (pass `pt`); the script refuses `noprog` with it, because the builder would then ignore the programme.
"""),
    ("""# 2. same-line
SL=""
""", """# 2. same-line
SL=""
if [ -n "${DEF_DAY_PROGRAMME_SOURCE:-}" ]; then
  case ":$SLF:" in *:noprog:*) die "DEF_DAY_PROGRAMME_SOURCE=$DEF_DAY_PROGRAMME_SOURCE needs sameline_flags without noprog (got $SLF)";; esac
  log "programme-source $DEF_DAY_PROGRAMME_SOURCE: every context's block prior from it (flags $SLF)"
fi
"""),
    ("""    ${DEF_DAY_FORCE:+--force "$DEF_DAY_FORCE"} >> "$LOG" 2>&1 || die "def_sameline.py failed\"""",
     """    ${DEF_DAY_FORCE:+--force "$DEF_DAY_FORCE"} ${DEF_DAY_PROGRAMME_SOURCE:+--programme-source "$DEF_DAY_PROGRAMME_SOURCE"} \\
    >> "$LOG" 2>&1 || die "def_sameline.py failed\""""),
]
for old, new in edits:
    if s.count(old) != 1:
        raise SystemExit(f"anchor count {s.count(old)} != 1: {old[:70]!r}")
for old, new in edits:
    s = s.replace(old, new)
if "--check" in sys.argv:
    print("all anchors found; not writing")
    raise SystemExit(0)
running = subprocess.run(["bash", "-c", "ps -axo comm,args | awk '$1 ~ /^bash$/ && /scripts\\/def_day\\.sh/' | grep -c ."],
                         capture_output=True, text=True).stdout.strip()
if running not in ("", "0"):
    raise SystemExit(f"{running} def_day.sh process(es) running: refusing to edit a live script")
p.write_text(s)
print(f"patched {p}")
