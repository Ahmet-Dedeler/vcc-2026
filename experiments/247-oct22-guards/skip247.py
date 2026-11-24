"""247 item B #18 — the whitener stage's skip branch now proves the existing operator's axis.

`verify_existing_operator(wout, wmeta, bundle)` is exercised directly on five cases, and the
driver itself is run (`--stages whitener`, NOT dry-run: the stage is decided before any work)
on the two that matter on the day.  The pre-247 skip branch had no check at all -- its verbatim
lines are recorded as the BLIND row.

    stale-1  A/B/C's real operator (no axis field; legacy table) vs data/challenge  -> SKIP OK
    stale-2  A/B/C's real operator vs 199's def-bundle                             -> REFUSE (different panel)
    stale-3  199's def199 refit (json predates 235: no axis field, sha not legacy) -> REFUSE (no digest)
    fresh-4  a json carrying def-bundle's own axis digest vs def-bundle             -> SKIP OK
    stale-5  a json carrying A/B/C's axis digest vs 238's ef-bundle                -> REFUSE (different panel)

    uv run python experiments/247-oct22-guards/skip247.py
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work" / "skip"
sys.path.insert(0, str(ROOT / "scripts"))

import def_pipeline as P  # noqa: E402

ABC = ROOT / "data/challenge"
DEF199 = ROOT / "experiments/199-def-runbook/work/def-bundle"
EF238 = ROOT / "experiments/238-essential-free-panel/work/ef-bundle"
ABC_OP = ROOT / "data/donors/whitener_union_pruned.npz"
ABC_JSON = ROOT / "data/donors/whitener_union_pruned.json"
DEF199_OP = ROOT / "experiments/199-def-runbook/work/whitener_union_pruned_def199.npz"
DEF199_JSON = ROOT / "experiments/199-def-runbook/work/whitener_union_pruned_def199.json"


def attempt(name: str, wout: Path, wmeta: Path, bundle: Path) -> dict:
    try:
        source = P.verify_existing_operator(wout, wmeta, bundle)
        row = {"case": name, "raised": False, "source": source}
        print(f"[{name}] SKIP ALLOWED  (axis matched via {source})")
    except SystemExit as exc:
        msg = str(exc)
        row = {"case": name, "raised": True, "message": msg}
        print(f"[{name}] RAISED:\n" + "\n".join("    " + ln for ln in msg.splitlines()))
    return row


def main() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    rows = []
    rows.append(attempt("stale-1 A/B/C operator vs A/B/C bundle", ABC_OP, ABC_JSON, ABC))
    rows.append(attempt("stale-2 A/B/C operator vs 199 def-bundle", ABC_OP, ABC_JSON, DEF199))
    rows.append(attempt("stale-3 199 def199 refit (pre-235 json) vs def-bundle", DEF199_OP, DEF199_JSON, DEF199))
    # synthetic jsons: a dummy npz (its sha is in no table) + a json carrying a chosen digest
    dummy = WORK / "dummy.npz"
    dummy.write_bytes(b"not an operator")
    def_axis, def_n = P.bundle_axis(DEF199)
    abc_axis, abc_n = P.bundle_axis(ABC)
    fresh = WORK / "fresh_def199.json"
    fresh.write_text(json.dumps({"n_cols": def_n, "axis_sha256": def_axis, "donors": list(P.PRODUCTION_DONORS)}))
    stale = WORK / "stale_abc_axis.json"
    stale.write_text(json.dumps({"n_cols": abc_n, "axis_sha256": abc_axis, "donors": list(P.PRODUCTION_DONORS)}))
    rows.append(attempt("fresh-4 json with def-bundle's own axis vs def-bundle", dummy, fresh, DEF199))
    rows.append(attempt("stale-5 json with A/B/C's axis vs 238 ef-bundle", dummy, stale, EF238))

    # the driver itself, stage `whitener` only (decided before any work; no fit is started)
    driver = {}
    for tag, bundle, expect_rc in (("driver-abc", ABC, 0), ("driver-def199", DEF199, 1)):
        cmd = ["uv", "run", "python", "scripts/def_pipeline.py", "--bundle", str(bundle),
               "--work", str(WORK / tag), "--name", tag, "--stages", "whitener",
               "--whitener-out", str(ABC_OP), "--whitener-id", "union_pruned"]
        proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
        tail = (proc.stdout + proc.stderr).strip().splitlines()
        driver[tag] = {"rc": proc.returncode, "expect_rc": expect_rc, "tail": tail[-8:]}
        print(f"\n[{tag}] rc={proc.returncode} (expected {expect_rc})")
        for ln in tail[-8:]:
            print(f"    {ln}")

    pre = (HERE / "work/pre247/def_pipeline.py").read_text().splitlines()
    blind = [ln for ln in pre if "exists, skipping" in ln or "if wout.exists() and wmeta.exists()" in ln]
    verdict = {
        "stale-1 skip allowed (legacy table)": not rows[0]["raised"],
        "stale-2 refused: DIFFERENT PANEL": rows[1]["raised"] and "DIFFERENT PANEL" in rows[1]["message"],
        "stale-3 refused: no axis digest": rows[2]["raised"] and "no axis digest" in rows[2]["message"],
        "fresh-4 skip allowed (metadata)": not rows[3]["raised"] and rows[3]["source"].startswith("metadata:"),
        "stale-5 refused: DIFFERENT PANEL": rows[4]["raised"] and "DIFFERENT PANEL" in rows[4]["message"],
        "driver on A/B/C skips the fit (rc 0)": driver["driver-abc"]["rc"] == 0,
        "driver on def-bundle REFUSES the stale operator (rc != 0)": driver["driver-def199"]["rc"] != 0,
    }
    out = {"cases": rows, "driver": driver, "blind_pre247_lines": blind,
           "bundle_axes": {"abc": abc_axis, "def199": def_axis}, "verdict": verdict}
    (WORK / "skip247.json").write_text(json.dumps(out, indent=1))
    print("\nBLIND (pre-247 skip branch, verbatim):")
    for ln in blind:
        print(f"    {ln}")
    print()
    for k, v in verdict.items():
        print(f"[{'PASS' if v else 'FAIL'}] {k}")
    ok = all(verdict.values())
    print(f"\nITEM B#18 {'PASS' if ok else 'FAIL'}")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
