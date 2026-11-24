"""Show that arming a sealed R7b surface does not read the ledger.

Imports lib553 only (numpy; no truth file). Does not call open_truth or truth_rows.
"""
from __future__ import annotations

import inspect
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "experiments/553-r7b-harness"))
import lib553 as S  # noqa: E402

src = inspect.getsource(S.arm_access)
assert "LEDGER" not in src and "CONFIRMATION" not in src
spec = S.surface_spec("r7bd5000")
assert spec["sealed"] is True, spec["sealed"]
assert spec["name"] == "r7bd5000"
S.arm_access("r7bd5000")
assert S._ACCESS_ARMED["r7bd5000"] is True
# the draft's sha is not in the ledger; arm_access does not notice
ledger = (ROOT / "docs/CONFIRMATION-LEDGER.md").read_text()
# This branch's working tree does not contain the draft (it landed on origin/master). The snapshot
# is byte-identical to that blob (git hash-object checked in the review).
draft = Path(__file__).resolve().parent / "upstream/CLAIM-R7-DRAFT.md"
import hashlib
h = hashlib.sha256(draft.read_bytes()).hexdigest()
print("r7bd5000 sealed:", spec["sealed"])
print("arm_access source:")
print(src)
print("draft sha256", h)
print("sha in ledger:", h in ledger)
print("arm_access set the flag anyway:", S._ACCESS_ARMED["r7bd5000"])
print('{"owner_go_ignores_ledger": true}')
