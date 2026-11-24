"""233 Phase 0.2 / Phase 1.2 — re-extract the CD4 banks for each panel.

The CD4 banks hold ONLY the panel's rows (`docs/DEF-RUNBOOK.md` §0.2), so a new panel needs
a new extraction from the pinned 16.8 GB source.  Three runs:

  abc   KNOWN-ANSWER: must reproduce the shipped cloud banks BYTE-FOR-BYTE on every data file
        (DEF-RUNBOOK §4 records 36 of 36, verified 2026-09-03).  Everything downstream is a
        ratio against this control, so if it does not reproduce, stop.
  real  / pess   the two synthetic panels' banks.

047's coverage gates (`--min-targets 275`, `--min-union-targets 290`) are ON for `abc` and
`real`.  On the pessimistic draw they are EXPECTED to fire (§5: 173/173/168, union 174) --
that is the finding, not a threshold to weaken -- so that run is made twice: once with the
gates on (recording the abort) and once with them at 0, which is DEF-RUNBOOK §5's pre-decided
response for exactly this case.

    uv run python cd4_233.py --panels abc real pess
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import panels233 as P  # noqa: E402

SOURCE = ROOT / "experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad"
SHIPPED = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"
EXTRACT = ROOT / "experiments/047-cd4-direct/extract_cd4.py"
CONDITIONS = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(1 << 22):
            h.update(chunk)
    return h.hexdigest()


def data_files(root: Path) -> dict[str, Path]:
    out = {}
    for cond in CONDITIONS:
        for f in sorted((root / cond).iterdir()):
            if f.suffix in (".npy", ".txt", ".csv"):
                out[f"{cond}/{f.name}"] = f
    return out


def run(panel: str, relax: bool) -> dict:
    out = P.cd4_root(panel)
    out.parent.mkdir(parents=True, exist_ok=True)
    bundle = P.BUNDLES[panel]
    cmd = ["uv", "run", "python", str(EXTRACT),
           "--input", str(SOURCE), "--output", str(out),
           "--panel-csv", str(bundle / "pert_counts.csv"),
           "--genes-csv", str(bundle / "gene_names.csv"),
           "--max-rss-gib", "4.0", "--overwrite"]
    if relax:
        cmd += ["--min-targets", "0", "--min-union-targets", "0"]
    t0 = time.time()
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    log = P.WORK / f"cd4-{panel}{'-relaxed' if relax else ''}.log"
    log.write_text(r.stdout + "\n----- stderr -----\n" + r.stderr)
    cov = {}
    for line in (r.stdout + r.stderr).splitlines():
        if "coverage" in line.lower() or "union" in line.lower():
            cov[line.strip()[:160]] = True
    return dict(panel=panel, relaxed=relax, rc=r.returncode,
                seconds=round(time.time() - t0, 1), log=str(log.name),
                tail=(r.stdout.strip().splitlines() or [""])[-6:],
                err=(r.stderr.strip().splitlines() or [""])[-4:])


def coverage_of(panel: str) -> dict:
    root = P.cd4_root(panel)
    man = json.loads((root / "manifest.json").read_text())
    per, union = {}, set()
    for cond in CONDITIONS:
        perts = [x.strip() for x in (root / cond / "perts.txt").read_text().split("\n")
                 if x.strip()]
        per[cond] = len(perts)
        union |= set(perts)
    return dict(per_condition=per, union=len(union),
                genes=len([x for x in (root / CONDITIONS[0] / "genes.txt")
                          .read_text().split("\n") if x.strip()]),
                manifest_keys=sorted(man)[:8])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", nargs="*", default=["abc", "real", "pess"])
    args = ap.parse_args()
    assert SOURCE.exists(), f"pinned CD4 source missing: {SOURCE}"

    results = {}
    for panel in args.panels:
        relax = panel == "pess"
        if relax:                                  # first, with the gates ON: expect an abort
            strict = run(panel, relax=False)
            results[f"{panel}|gates_on"] = strict
            print(f"[{panel}] gates ON -> rc {strict['rc']} "
                  f"({'ABORTED as DEF-RUNBOOK §5 predicts' if strict['rc'] else 'PASSED'})")
            for ln in strict["err"]:
                print("      ", ln[:150])
        r = run(panel, relax=relax)
        results[panel] = r
        print(f"[{panel}] rc {r['rc']} in {r['seconds']}s"
              + ("  (gates relaxed to 0, §5's pre-decided response)" if relax else ""))
        for ln in r["tail"]:
            print("      ", ln[:150])
        if r["rc"] != 0:
            for ln in r["err"]:
                print("   ERR", ln[:200])
            raise SystemExit(f"cd4 extraction failed for {panel}")
        results[f"{panel}|coverage"] = coverage_of(panel)
        print("      coverage:", json.dumps(results[f"{panel}|coverage"]))

    if "abc" in args.panels:                        # KNOWN-ANSWER ROW
        mine, ship = data_files(P.cd4_root("abc")), data_files(SHIPPED)
        common = sorted(set(mine) & set(ship))
        diffs = [k for k in common if sha256(mine[k]) != sha256(ship[k])]
        results["KA_cd4_bytes"] = dict(
            n_files_mine=len(mine), n_files_shipped=len(ship), n_compared=len(common),
            only_mine=sorted(set(mine) - set(ship)), only_shipped=sorted(set(ship) - set(mine)),
            n_identical=len(common) - len(diffs), differing=diffs)
        print(f"\nKA-B  CD4 re-extraction vs the shipped banks: "
              f"{len(common) - len(diffs)}/{len(common)} data files byte-identical"
              + (f"  DIFFER: {diffs}" if diffs else "  -> PASS"))

    (P.WORK / "cd4_233.json").write_text(json.dumps(results, indent=2))
    print("-> work/cd4_233.json")


if __name__ == "__main__":
    main()
