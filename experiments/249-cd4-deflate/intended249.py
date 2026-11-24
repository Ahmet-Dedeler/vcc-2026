"""249 -- reconstruct the INTENDED builder from git HEAD: HEAD + patch249.py + the one follow-up
edit (the deflation log line made safe for a supplied control basis).  Writes
`work/make_submission_pm_cd4.249.py` and `work/builder.intended.sha256`; the heavy stages gate on
the working tree matching this byte for byte (244's rehearsal pins the registry in the same file
mid-run and restores it afterwards).

    uv run python experiments/249-cd4-deflate/intended249.py
"""
import hashlib
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
out = WORK / "make_submission_pm_cd4.249.py"
head = subprocess.run(["git", "-C", str(ROOT), "show", "HEAD:scripts/make_submission_pm_cd4.py"],
                      capture_output=True, text=True, check=True).stdout
pre = WORK / "make_submission_pm_cd4.PRE.py"
assert pre.read_text() == head, "git HEAD moved since the PRE byte copy was taken"
out.write_text(head)
subprocess.run([sys.executable, str(HERE / "patch249.py"), str(out)], check=True)
s = out.read_text()
old = "centred energy {deflation_audit['top_k_energy_frac_of_centred']:.6f}), \"\n"
new = "centred energy {deflation_audit.get('top_k_energy_frac_of_centred', float('nan')):.6f}), \"\n"
assert s.count(old) == 1
out.write_text(s.replace(old, new))
sha = hashlib.sha256(out.read_bytes()).hexdigest()
(WORK / "builder.intended.sha256").write_text(sha + "\n")
cur = hashlib.sha256((ROOT / "scripts/make_submission_pm_cd4.py").read_bytes()).hexdigest()
print(f"intended builder {out}  sha256 {sha}")
print(f"working tree     sha256 {cur}  {'== intended' if cur == sha else '!= intended (244 pins present or edits missing)'}")
