"""254 -- reconstruct the INTENDED builder from git HEAD: HEAD + patch254.py.  Writes
`work/make_submission_pm_cd4.254.py` and `work/builder.intended.sha256`; every heavy stage gates
on the working tree matching this byte for byte (244's rehearsal pins the registry in the same
file mid-run and restores it afterwards; 249 section 12).

    uv run python experiments/254-partial-rows/intended254.py
"""
import hashlib
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
from patch254 import patch  # noqa: E402

out = WORK / "make_submission_pm_cd4.254.py"
head = subprocess.run(["git", "-C", str(ROOT), "show", "HEAD:scripts/make_submission_pm_cd4.py"],
                      capture_output=True, text=True, check=True).stdout
pre = WORK / "make_submission_pm_cd4.PRE.py"
assert pre.read_text() == head, "git HEAD moved since the PRE byte copy was taken"
out.write_text(patch(head))
sha = hashlib.sha256(out.read_bytes()).hexdigest()
(WORK / "builder.intended.sha256").write_text(sha + "\n")
cur = hashlib.sha256((ROOT / "scripts/make_submission_pm_cd4.py").read_bytes()).hexdigest()
pre_sha = hashlib.sha256(pre.read_bytes()).hexdigest()
print(f"PRE (git HEAD) builder sha256 {pre_sha}")
print(f"intended builder {out}  sha256 {sha}")
print(f"working tree     sha256 {cur}  {'== intended' if cur == sha else '!= intended (244 pins present or edits missing)'}")
raise SystemExit(0 if cur == sha else 1)
