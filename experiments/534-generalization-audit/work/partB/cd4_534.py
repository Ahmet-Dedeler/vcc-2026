"""534 Part B -- extract the CD4 banks for one R6 panel: 238's cd4_238.py with the panel changed (047's extract_cd4.py on
the pinned 16.8 GB source, --max-rss-gib 4). Gates ON first (047's coverage gates are expected to fire on a panel
disjoint from A/B/C, exactly as on 238's draws; recorded), then the pre-decided response (DEF-RUNBOOK section 5):
--min-targets 0 --min-union-targets 0. Light (~2 s, < 0.5 GB); not a heavy process.

    uv run python cd4_534.py <panel>          -> work/partB/<panel>/cd4/cd4_compact, work/partB/<panel>/cd4_534.json
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard534  # noqa: E402,F401
import json  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE = ROOT / "experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad"
EXTRACT = ROOT / "experiments/047-cd4-direct/extract_cd4.py"
CONDITIONS = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")


def run(panel: str, relax: bool) -> dict:
    bundle = HERE / panel / "bundle"
    out = HERE / panel / "cd4" / "cd4_compact"
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["uv", "run", "python", str(EXTRACT), "--input", str(SOURCE), "--output", str(out),
           "--panel-csv", str(bundle / "pert_counts.csv"), "--genes-csv", str(bundle / "gene_names.csv"),
           "--max-rss-gib", "4.0", "--overwrite"]
    if relax:
        cmd += ["--min-targets", "0", "--min-union-targets", "0"]
    t0 = time.time()
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    log = HERE / panel / f"cd4-{panel}{'-relaxed' if relax else ''}.log"
    log.write_text(r.stdout + "\n----- stderr -----\n" + r.stderr)
    return dict(relaxed=relax, rc=r.returncode, seconds=round(time.time() - t0, 1), log=log.name,
                err=(r.stderr.strip().splitlines() or [""])[-2:])


def coverage(panel: str) -> dict:
    root = HERE / panel / "cd4" / "cd4_compact"
    per, union = {}, set()
    for cond in CONDITIONS:
        perts = [x.strip() for x in (root / cond / "perts.txt").read_text().splitlines() if x.strip()]
        per[cond] = len(perts)
        union |= set(perts)
    genes = len([x for x in (root / CONDITIONS[0] / "genes.txt").read_text().splitlines() if x.strip()])
    return dict(per_condition=per, union=len(union), genes=genes)


def main(panel: str) -> None:
    res = {"gates_on": run(panel, relax=False)}
    print(f"[{panel}] gates ON -> rc {res['gates_on']['rc']} ({'ABORTED as expected (DEF-RUNBOOK 5)' if res['gates_on']['rc'] else 'PASSED'})",
          flush=True)
    res["relaxed"] = run(panel, relax=True)
    assert res["relaxed"]["rc"] == 0, f"cd4 extraction failed: {res['relaxed']}"
    res["coverage"] = coverage(panel)
    print(f"[{panel}] relaxed rc 0 in {res['relaxed']['seconds']} s; coverage {res['coverage']}", flush=True)
    (HERE / panel / "cd4_534.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
