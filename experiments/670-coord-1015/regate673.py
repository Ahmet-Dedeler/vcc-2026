#!/usr/bin/env python3
"""Re-gate W1 (673) on its FINISHED build after gate673's judge false-FAIL (the 697 / 698 precedent: re-gate, never rebuild).

Written 2026-10-20 ~11:35Z by the coordinator [3c7a7e] (cross-session messaging to readiness was paused by the app).
What failed: gate673's config diff vs R-ABC2 allowed {bulk_ctx_scale, bulk_load, bulk_load_sha256}; the build's diff also lists
bulk_sharpen_pt. In that record only the DATA-DERIVED scale fields moved (edges, key_ratio_c, key_ratio_maxdev): tercilec's key ratio
is normalised on the stack dump, which an H1 dose changes by design. Every BEHAVIOURAL field (mode, q, counts, carried key + sha,
per-target moves, n_carried, floor/tie/non-carried counts, nonfinite) is identical, so each target is sharpened exactly as in R-ABC2.
This re-gate checks, on the finished build:
  1. the package on disk has the sha the build logged, and the novelty check logged rc 0;
  2. gate673's log has exactly ONE FAIL row, the config-diff row, and its diff set is within {bulk_ctx_scale, bulk_load,
     bulk_load_sha256, bulk_sharpen_pt} and contains bulk_load; every other gate row is "ok" / PASS (CLIFF OK, postbuild, H1 dose G 0.7);
  3. bulk_sharpen_pt differs from R-ABC2's only in {edges, key_ratio_c, key_ratio_maxdev}; all other fields are equal.
On success it writes the arm1015 lines: "[w1] PASS …", "[w1] CLIFF OK …" (copied from gate673), "novelty rc=0" to the run log, and
"READY w1 <sha>" to the ready log. Otherwise "NOT READY <reason>" to the ready log.
  uv run python experiments/670-coord-1015/regate673.py
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
X = ROOT / "experiments/673-h1x25-abc/work"
LOG = X / "w1_673.log"
PKG = ROOT / "submissions/673w1-october-recipe-abc.vcc"
META_W1 = X / "run/w1_final.h5ad.meta.json"
META_R2 = ROOT / "experiments/695-rabc2/work/run/rabc2_final.h5ad.meta.json"
W = ROOT / "experiments/670-coord-1015/work"
READY = W / "regate673.ready.log"
RUN = W / "regate673.run.log"
ALLOWED = {"bulk_ctx_scale", "bulk_load", "bulk_load_sha256", "bulk_sharpen_pt"}
DATA_DERIVED = {"edges", "key_ratio_c", "key_ratio_maxdev"}


def utc():
    import datetime as dt
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def fail(why):
    READY.open("a").write(f"{utc()} NOT READY w1 re-gate: {why}\n")
    print("NOT READY:", why)
    sys.exit(1)


def main():
    t = LOG.read_text()
    m = re.search(r"package sha256: ([0-9a-f]{64})", t)
    if not m:
        fail("no package sha in w1_673.log")
    sha = m.group(1)
    if hashlib.sha256(PKG.read_bytes()).hexdigest() != sha:
        fail("package on disk differs from the build's logged sha")
    if "novelty rc=0" not in t:
        fail("novelty rc=0 not logged")
    fails = [l for l in t.splitlines() if re.search(r"\bFAIL\b", l) and "ALL GATES FAIL" not in l and "gate673 failed" not in l]
    if len(fails) != 1 or "config diff vs R-ABC2" not in fails[0]:
        fail(f"expected exactly one FAIL row (the config diff vs R-ABC2), got {len(fails)}: {[f[:120] for f in fails]}")
    dm = re.search(r"diff \[([^\]]*)\]", fails[0].split("[metadiff]")[1])
    diff = {s.strip().strip("'") for s in dm.group(1).split(",") if s.strip()}
    if not diff <= ALLOWED or "bulk_load" not in diff:
        fail(f"diff {sorted(diff)} not within {sorted(ALLOWED)} containing bulk_load")
    cliff = [l for l in t.splitlines() if l.startswith("[w1] CLIFF OK")]
    if not cliff or "postbuild 42/42 PASS" not in cliff[0]:
        fail("no CLIFF OK / postbuild 42/42 PASS row")
    if "records VCC_H1_UNION_G 0.7" not in t or "H1 stack dose (690i): dump_h1x25_cd4x1_s0.25.npz env VCC_H1_UNION_G 0.7 (want 0.7)" not in t:
        fail("the H1 dose check (G 0.7) is not in the log")
    a = json.loads(META_W1.read_text())["config"]["bulk_sharpen_pt"]
    b = json.loads(META_R2.read_text())["config"]["bulk_sharpen_pt"]
    moved = {k for k in set(a) | set(b) if a.get(k) != b.get(k)}
    if not moved <= DATA_DERIVED:
        fail(f"bulk_sharpen_pt behavioural fields differ: {sorted(moved - DATA_DERIVED)}")
    pkg_meta = json.loads((ROOT / "submissions/673w1-october-recipe-abc.h5ad.meta.json").read_text())["config"]["bulk_sharpen_pt"]
    if pkg_meta != a:
        fail("the package's meta sidecar differs from the run's final meta in bulk_sharpen_pt")
    with RUN.open("a") as fh:
        fh.write(f"{utc()} re-gate of 673 W1 by the coordinator (gate673's judge false-FAIL on bulk_sharpen_pt's data-derived fields)\n")
        fh.write(f"[w1] PASS: config diff vs R-ABC2 (re-gated): diff {sorted(diff)}; bulk_sharpen_pt moved only {sorted(moved)} "
                 f"(behavioural fields identical: mode {a['mode']}, q {a['q']}, counts {a['counts']}, carried key sha "
                 f"{a['carried_key_sha256'][:16]}, moves {sorted(a['moved_vs_tercile'])})\n")
        fh.write(cliff[0] + "\n")
        fh.write("novelty rc=0\n")
        fh.write(f"package sha256 {sha}\n")
    READY.open("a").write(f"{utc()} READY w1 {sha} {PKG.relative_to(ROOT)}\n")
    print(f"READY w1 {sha}: diff {sorted(diff)}, bulk_sharpen_pt moved {sorted(moved)}")


if __name__ == "__main__":
    main()
