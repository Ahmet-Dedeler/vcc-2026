"""Resume and sealing checks against r7a551's real control flow, with guard551 stubbed and every
path redirected under this directory. No T-panel truth file is opened: the only inputs are tiny
synthetic bundles (a pert_counts and a gene_names).
"""
from __future__ import annotations

import importlib.util
import json
import sys
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent / "resume_out"
OUT.mkdir(parents=True, exist_ok=True)

sys.modules.setdefault("guard551", types.ModuleType("guard551"))
sys.path.insert(0, str(ROOT / "experiments/551-r6-c10b"))

src_path = ROOT / "experiments/555-cloud-r7-review/repro/upstream/r7a551.py"
spec = importlib.util.spec_from_file_location("r7a551_under_test", src_path)
# __file__ must be the real 551 path so HERE/ROOT resolve to the repo. We then retarget WORK.
R = importlib.util.module_from_spec(spec)
R.__file__ = str(ROOT / "experiments/551-r6-c10b/r7a551.py")
spec.loader.exec_module(R)

# retarget all writes
R.WORK = OUT / "work"
R.R7A = R.WORK / "r7a"
R.STATE = R.R7A / "state.json"
R.CLAIM = R.R7A / "claim.json"
R.R7A.mkdir(parents=True, exist_ok=True)
LEDGER = OUT / "ledger.md"
R.LEDGER = LEDGER


def reset_state(obj):
    if R.STATE.exists():
        R.STATE.unlink()


def test_second_record_and_ledger_substring():
    print("=== record: substring ledger check, second record refused, config not read from the claim ===")
    claim = OUT / "CLAIM-not-the-real-one.md"
    claim.write_text("this file says F = M5\n")
    h = R.sha256(claim)
    # the hash appears only as a substring of a sentence, not as an R7-section record
    LEDGER.write_text(f"# notes\n\nDo not arm anything with {h} until October.\n")
    ns = types.SimpleNamespace(claim=str(claim.relative_to(ROOT)), config="M6+P2g",
                               window_start="2026-10-08T06:00Z", window_end="2026-10-08T18:00Z",
                               received="2026-10-08T05:00Z", message_file=str(OUT / "msg.txt"))
    (OUT / "msg.txt").write_text("coordinator says go\n")
    # plan() would hash real bundles and git; skip it
    R.plan = lambda write=False: {"instrument": {"blobs": {}}, "external": {}}
    R.record(ns)
    rec = json.loads(R.CLAIM.read_text())
    print(f"  recorded config {rec['config']} while the claim text says M5; sha in ledger by substring: "
          f"{h in LEDGER.read_text()}")
    raised = False
    try:
        R.record(ns)
    except AssertionError as e:
        raised = True
        print(f"  second record refused: {e}")
    # a hash that is NOT in the ledger is refused
    LEDGER.write_text("no hashes here\n")
    R.CLAIM.unlink()
    refused = False
    try:
        R.record(ns)
    except AssertionError as e:
        refused = True
        print(f"  absent hash refused: {e}")
    return raised and refused and rec["config"] == "M6+P2g"


def test_ops_deleted_resume():
    print("=== resume after a crash mid-deletion: npz already gone, fit steps marked done ===")
    reset_state(None)
    p = "T1"
    pw = R.WORK / p
    ops = pw / "ops"
    ops.mkdir(parents=True, exist_ok=True)
    # one operator json kept, npz already unlinked (crash during the deletion loop)
    meta = {"sha256": "abc123", "axis_sha256": "def456", "n_cols": 3, "donors": ["k562_gwps", "orion_hct116"]}
    for tag in R.FITS:
        (ops / f"T1__{tag}.json").write_text(json.dumps(meta))
    # every npz is already gone (the crash was mid-unlink, after the jsons were kept)
    state = {}
    for tag in R.FITS:
        state[f"{p}:fit:{tag}"] = {"done": "t", "rc": 0}
    state[f"{p}:cd4"] = {"done": "t"}
    state[f"{p}:axes"] = {"done": "t"}
    state[f"{p}:full551"] = {"done": "t", "rc": 0}
    state[f"{p}:comp551"] = {"done": "t", "rc": 0}
    R.STATE.write_text(json.dumps(state))
    # comp jsons so the assert before deletion passes
    for f in R.CANCER:
        (pw / f"comp551-{f}.json").write_text(json.dumps({"peak_rss_gib": 1.0, "op_axis_checked": True, "config": "M6"}))
    # the real heavy() returns before G.heavy when the step is already in state. If a done step
    # still reaches the gate, this raises and the test fails for that reason.
    def boom(*a, **k):
        raise AssertionError(f"G.heavy was called for a done step: {a[1] if len(a) > 1 else a}")
    R.G.heavy = boom
    try:
        R.panel(p, "M6")
        print("  panel() returned without error (unexpected)")
        return False
    except FileNotFoundError as e:
        print(f"  FileNotFoundError as predicted: {e}")
        # and the fit steps are still marked done, so a retry cannot rebuild the npz
        s = json.loads(R.STATE.read_text())
        still = all(f"{p}:fit:{tag}" in s for tag in R.FITS) and f"{p}:ops_deleted" not in s
        print(f"  fit steps still done and ops_deleted absent: {still}")
        return still


def test_window_does_not_cover_aggregate_skip():
    print("=== done heavy steps are skipped without consulting the window ===")
    # structural: heavy() returns before G.heavy when the step is in state
    src = src_path.read_text().splitlines()
    # quote the relevant lines (1-indexed in the file, which is the committed r7a551.py)
    body = "\n".join(src[283:294])  # def heavy ... return r   line numbers of the snapshot match the commit
    print(body)
    return "if step in state()" in body and "G.heavy" in body


def main():
    a = test_second_record_and_ledger_substring()
    b = test_ops_deleted_resume()
    c = test_window_does_not_cover_aggregate_skip()
    print("===", json.dumps({"ledger": a, "resume_unlink": b, "heavy_skip": c}))
    if not (a and b and c):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
