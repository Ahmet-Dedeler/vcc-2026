"""523: make the DE call set's CPM-block size an env knob, VCC_RC_CALLSET_CPM (default unset = 2200, bit-identical).
The head (400) and z_cap (1.90) stay constants. The value is already written to the meta as config.rc_callset_cpm, so a
set knob shows up in the rule-10 diff as exactly that key.
usage: uv run python experiments/523-callset-size/patch_builder523.py <builder> [--check]
"""
import sys
from pathlib import Path

p = Path(sys.argv[1]); s = p.read_text()
old = "RC_CALLSET_CPM = 2200\n"
new = ('RC_CALLSET_CPM = 2200\n'
       '# 523 (2026-09-28): VCC_RC_CALLSET_CPM overrides the CPM block size for a call-set-size dose arm (521 section 8).\n'
       '# Unset = 2200, bit-identical. Validated as an integer in [100, 5000]; the value lands in config.rc_callset_cpm.\n'
       'if os.environ.get("VCC_RC_CALLSET_CPM", "").strip():\n'
       '    _cpm = os.environ["VCC_RC_CALLSET_CPM"].strip()\n'
       '    if not _cpm.isdigit() or not (100 <= int(_cpm) <= 5000):\n'
       '        raise SystemExit(f"VCC_RC_CALLSET_CPM must be an integer in [100, 5000], got {_cpm!r}")\n'
       '    RC_CALLSET_CPM = int(_cpm)\n')
if s.count(old) != 1:
    raise SystemExit(f"anchor count {s.count(old)} != 1")
if "VCC_RC_CALLSET_CPM" in s:
    raise SystemExit("already patched")
head = s[: s.index(old)]
if "\nimport os\n" not in head and not head.startswith("import os\n"):
    raise SystemExit("os is not imported above the anchor")
s2 = s.replace(old, new, 1)
compile(s2, str(p), "exec")
if "--check" in sys.argv:
    print("patch applies cleanly (not written)")
else:
    p.write_text(s2)
    print(f"patched {p}")
