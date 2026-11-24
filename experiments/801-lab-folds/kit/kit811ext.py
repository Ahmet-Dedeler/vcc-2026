"""801 kit (follow-up 7): a COPY of experiments/811-lolo-model/ext811a.py (git blob 55fbda142ca37e6b127020136bc687455a32443f) with the
panel list as a parameter: env KIT811_PANELS (space-separated kit tags, e.g. "DEF"); HERE points at 811 (its banks, its
ext811 module, its work/), so it writes exactly where ext811a writes, on the EC2 checkout. Nothing else changed. ext811a's
docstring follows.

811 -- extend the compact banks to 801's A-only panels (A1-A6), so the many-line consensus on an A fold has the same lines
as on every other fold. ext811.py's banks cover only the A/B/C + L1-L6 target union; the A panels are disjoint from it.

1. CD4: 047's extract_cd4.py (the extractor 801's cd4_801.py calls, same arguments incl. the relaxed gates) on each A panel
   and, as a KNOWN ANSWER, on L1 -> work/cd4a/<panel>/cd4_compact. L1's delta/genes/perts must be byte-identical (sha256) to
   801's own L1 compact (experiments/801-lab-folds/work/L1/cd4/cd4_compact), else abort.
2. Every extra bank (ext811.SIMPLE + the 3 CD4 conditions): rows for the A panels' targets via ext811.read_bank (unchanged).
3. Merge with work/banks/<bank>.npz -> work/banks_ext/<bank>.npz: old rows kept byte-for-byte, A rows appended; a target in
   both must have identical rows (asserted). work/banks_ext/manifest.json records every sha256.

    bash experiments/602-council-gate/gate602.sh ext811a 3 -- uv run python experiments/811-lolo-model/ext811a.py
"""
from __future__ import annotations

import csv
import hashlib
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[2] / "811-lolo-model"     # 801 kit: 811's dir
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import ext811 as X  # noqa: E402

W801 = ROOT / "experiments/801-lab-folds/work"
SOURCE = ROOT / "experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad"
EXTRACT = ROOT / "experiments/047-cd4-direct/extract_cd4.py"
CD4A = HERE / "work" / "cd4a"
OLD, NEW = HERE / "work" / "banks", HERE / "work" / "banks_ext"
import os  # noqa: E402
APANELS = os.environ["KIT811_PANELS"].split()                      # 801 kit: the kit tags


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def cd4(panel):
    bundle = W801 / panel / "bundle"
    out = CD4A / panel / "cd4_compact"
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["uv", "run", "python", str(EXTRACT), "--input", str(SOURCE), "--output", str(out),
           "--panel-csv", str(bundle / "pert_counts.csv"), "--genes-csv", str(bundle / "gene_names.csv"),
           "--max-rss-gib", "4.0", "--overwrite", "--min-targets", "0", "--min-union-targets", "0"]
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    (CD4A / panel / "extract.log").write_text(r.stdout + "\n----- stderr -----\n" + r.stderr)
    assert r.returncode == 0, f"extract_cd4 {panel} rc {r.returncode}"
    return out


def main():
    t0 = time.time()
    man = dict(cd4={}, banks={})
    # 1. CD4 compacts + the L1 known answer
    l1 = cd4("L1")
    for c in X.CD4_CONDS:
        for f in ("delta.npy", "genes.txt", "perts.txt"):
            a, b = sha(l1 / c / f), sha(W801 / "L1/cd4/cd4_compact" / c / f)
            assert a == b, f"KA FAILED: L1 {c}/{f} {a} != 801's {b}"
    man["ka_L1_cd4"] = "PASS: delta/genes/perts byte-identical to 801's L1 compact (3 conditions)"
    print("KA L1 CD4 compact == 801's: PASS", flush=True)
    dirs = {p: cd4(p) for p in APANELS}
    man["cd4"] = {p: {c: sha(d / c / "delta.npy") for c in X.CD4_CONDS} for p, d in dirs.items()}
    # 2. rows for the A panels' targets
    want = set()
    for p in APANELS:
        with open(W801 / p / "bundle/pert_counts.csv") as f:
            want |= {r["target_gene"] for r in csv.DictReader(f)}
    jobs = {b: [D] for b, D in X.SIMPLE.items()}
    for c in X.CD4_CONDS:
        jobs[c] = [dirs[p] / c for p in APANELS]
    NEW.mkdir(parents=True, exist_ok=True)
    for b, ds in jobs.items():
        rows, measured, src = X.read_bank(ds, want)
        z = np.load(OLD / f"{b}.npz")
        T0, R0, M0 = [str(t) for t in z["targets"]], np.asarray(z["rows"]), np.asarray(z["measured"], bool)
        pos0 = {t: i for i, t in enumerate(T0)}
        add = [t for t in sorted(rows) if t not in pos0]
        for t in sorted(rows):
            if t in pos0:
                assert np.array_equal(R0[pos0[t]], rows[t]), f"{b}: {t} differs between old and new extraction"
        T = T0 + add
        R = np.vstack([R0] + [rows[t][None, :] for t in add]) if add else R0
        M = M0 | measured if b in X.SIMPLE else M0 | measured
        dst = NEW / f"{b}.npz"
        np.savez_compressed(dst, targets=np.asarray(T), rows=R, measured=M)
        man["banks"][b] = dict(old_sha256=sha(OLD / f"{b}.npz"), n_old=len(T0), n_added=len(add), sources=src, out_sha256=sha(dst))
        print(f"{b:16s} old {len(T0):5d} + A {len(add):4d} ({time.time() - t0:.0f}s)", flush=True)
    man["peak_rss_gib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30
    (NEW / "manifest.json").write_text(json.dumps(man, indent=1))
    print(f"done {time.time() - t0:.0f}s peak {man['peak_rss_gib']:.2f} GiB", flush=True)


if __name__ == "__main__":
    main()
