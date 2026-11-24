"""238 Phase 1b — re-extract the CD4 banks for the two essential-free panels.

The CD4 banks hold ONLY the panel's rows (`DEF-RUNBOOK` §0.2), so a new panel needs a new
extraction from the pinned 16.8 GB source.  This is `cd4_233.py` with the panel list changed.

KNOWN-ANSWER ROW.  The `abc` extraction must reproduce the SHIPPED cloud banks byte-for-byte
on every data file (`DEF-RUNBOOK` §4: 36 of 36, verified 2026-09-03; 233 KA-B re-verified it).
Everything downstream is a ratio against that control, so it is re-run here rather than
inherited -- a cheap 2 s check that the extractor and the pinned source are the ones 233 used.

047's coverage gates (`--min-targets 275`, `--min-union-targets 290`) are expected to FIRE on
both essential-free draws, exactly as they fired on 199's pessimistic one: a panel disjoint
from A/B/C simply is not covered by a bank extracted for A/B/C's biology.  `DEF-RUNBOOK` §5's
pre-decided response for that case is to re-run at `--min-targets 0 --min-union-targets 0`
after confirming the mapping, which is what happens here, and the abort is recorded.

    uv run python cd4_238.py --panels abc efp efr
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
import panels238 as P  # noqa: E402

SOURCE = ROOT / "experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad"
SHIPPED = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"
EXTRACT = ROOT / "experiments/047-cd4-direct/extract_cd4.py"
CONDITIONS = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
KA_ABC = dict(per_condition={"cd4_rest": 283, "cd4_stim8hr": 287, "cd4_stim48hr": 284},
              union=291, genes=10113)


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
    return dict(panel=panel, relaxed=relax, rc=r.returncode,
                seconds=round(time.time() - t0, 1), log=log.name,
                tail=(r.stdout.strip().splitlines() or [""])[-6:],
                err=(r.stderr.strip().splitlines() or [""])[-4:])


def coverage_of(panel: str) -> dict:
    root = P.cd4_root(panel)
    per, union = {}, set()
    for cond in CONDITIONS:
        perts = [x.strip() for x in (root / cond / "perts.txt").read_text().splitlines()
                 if x.strip()]
        per[cond] = len(perts)
        union |= set(perts)
    genes = len([x for x in (root / CONDITIONS[0] / "genes.txt").read_text().splitlines()
                 if x.strip()])
    return dict(per_condition=per, union=len(union), genes=genes)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", nargs="*", default=["abc", "efp", "efr"])
    args = ap.parse_args()
    assert SOURCE.exists(), f"pinned CD4 source missing: {SOURCE}"

    results: dict = {}
    for panel in args.panels:
        relax = panel != "abc"
        if relax:                                  # gates ON first: record whether they fire
            strict = run(panel, relax=False)
            results[f"{panel}|gates_on"] = strict
            print(f"[{panel}] gates ON -> rc {strict['rc']} "
                  f"({'ABORTED (DEF-RUNBOOK §5 case)' if strict['rc'] else 'PASSED'})")
            for ln in strict["err"][-2:]:
                print("      ", ln[:150])
        r = run(panel, relax=relax)
        results[panel] = r
        if r["rc"] != 0:
            for ln in r["err"]:
                print("   ERR", ln[:200])
            raise SystemExit(f"cd4 extraction failed for {panel}")
        cov = coverage_of(panel)
        results[f"{panel}|coverage"] = cov
        print(f"[{panel}] rc 0 in {r['seconds']}s"
              + ("  (gates relaxed to 0, §5's pre-decided response)" if relax else "")
              + f"   coverage {cov['per_condition']} union {cov['union']} "
                f"genes {cov['genes']}")

    if "abc" in args.panels:                        # KNOWN-ANSWER ROW
        mine, ship = data_files(P.cd4_root("abc")), data_files(SHIPPED)
        common = sorted(set(mine) & set(ship))
        diffs = [k for k in common if sha256(mine[k]) != sha256(ship[k])]
        cov = results["abc|coverage"]
        cov_ok = (cov["per_condition"] == KA_ABC["per_condition"]
                  and cov["union"] == KA_ABC["union"] and cov["genes"] == KA_ABC["genes"])
        results["KA_cd4"] = dict(n_compared=len(common), n_identical=len(common) - len(diffs),
                                 differing=diffs, coverage_matches_known_answer=cov_ok)
        print(f"\nKA  CD4 abc re-extraction vs the shipped banks: "
              f"{len(common) - len(diffs)}/{len(common)} data files byte-identical; "
              f"coverage {'PASS' if cov_ok else '**FAIL**'}")
        if diffs or not cov_ok:
            raise SystemExit(f"KA FAILED: {len(diffs)} differing files, coverage_ok={cov_ok}")

    (P.WORK / "cd4_238.json").write_text(json.dumps(results, indent=2))
    print(f"-> {P.WORK / 'cd4_238.json'}")


if __name__ == "__main__":
    main()
