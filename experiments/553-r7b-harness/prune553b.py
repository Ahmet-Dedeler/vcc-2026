#!/usr/bin/env python3
"""553 -- the second prune (2026-10-02, approved by the coordinator for 10-04/10-05 scratch): regenerable ARRAYS only, by
explicit path, each one hashed and recorded before it is deleted.

What goes (arrays >= 5 MB: .npz / .parquet / .h5ad; every json / md / csv / log and every small array stays):
  draw-arm      h1r7[sN]_<arm>, arm != cpm      the priced arms on the five rehearsal draws (results in access553_h1r7*)
  chalf-arm     h1r7c<k><h>_<arm>               arms on the rehearsal's control halves (results in access553_h1r7c*)
  chalf-surf    h1r7c<k><h>                     the rehearsal's control-half prediction sides
  access        access553_<x>, x != h1          truth-side caches of finished accesses (access553.json / md kept)
  r7b-arm       r7b_<arm>                       arms on R7b as built (unfit for the production emitter; not the access)
  sweep         r7bd6000 / 7000 / 8000 (+ arms) the floor sweep's other floors (R7 amendment 3 chose 5,000)
What stays, whatever its size (refused even if a category matched):
  R7b's access inputs   r7b/ (operator, stage, effects, cd4, callblock), r7bd5000/, r7bd5000c<k><h>/ and every arm on them
  re-pricing on h1r7    h1r7, h1r7s2..s5 (controlref + callblock), the base arm h1r7[sN]_cpm, h1/ (bulk553.npz, the cache)
  known answers         h1_<arm> (KA-H1 / KA-S / KA-FINAL), access553_h1 (KA-ACCESS-FINAL's reference, H1-full's --retro),
                        access553_rehearsal5 (R7b's --retro), access553_r7bd5000-xBST2 (the dry-run plan), h1cb*, h1ka, test553

    uv run python experiments/553-r7b-harness/prune553b.py            # dry run: every file that would go, per category
    uv run python experiments/553-r7b-harness/prune553b.py --delete   # hash -> work/PRUNED-1002.tsv -> delete -> PRUNED.md
"""
from __future__ import annotations

import hashlib
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = (HERE / "work").resolve()
MIN_BYTES = 5_000_000
ARRAY_SUFFIXES = {".npz", ".parquet", ".h5ad"}
CATEGORIES = [  # (category, directory-name regex); first match wins
    ("draw-arm", r"h1r7(s\d+)?_(?!cpm$).+"),
    ("chalf-arm", r"h1r7(s\d+)?c\d[AB]_.+"),
    ("chalf-surf", r"h1r7(s\d+)?c\d[AB]"),
    ("access", r"access553_(?!h1$|rehearsal5$|r7bd5000-xBST2$).+"),
    ("r7b-arm", r"r7b_.+"),
    ("sweep", r"r7bd(6000|7000|8000)(_.+)?"),
]
KEEP = [r"r7b", r"r7bd5000.*", r"h1r7(s\d+)?", r"h1r7(s\d+)?_cpm", r"h1", r"h1_.+", r"access553_h1",
        r"access553_rehearsal5", r"access553_r7bd5000-xBST2", r"h1cb.*", r"h1ka", r"test553"]


def category(dirname: str) -> str | None:
    if any(re.fullmatch(k, dirname) for k in KEEP):
        return None
    for cat, rx in CATEGORIES:
        if re.fullmatch(rx, dirname):
            return cat
    return None


def busy() -> list[str]:
    """Running interpreters whose argv names one of 553's scripts (interpreter + path, AGENTS 6.2), other than this one."""
    out = subprocess.run(["ps", "-axo", "comm=,args="], capture_output=True, text=True).stdout.splitlines()
    pat = re.compile(r"(run|access|prep|diag|band|effects|ka|emitdiag)553\.py|drive553\.sh")
    return [ln for ln in out if re.match(r"^\S*(python[0-9.]*|Python|bash)\s", ln.strip()) and pat.search(ln)]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 23), b""):
            h.update(chunk)
    return h.hexdigest()


def plan() -> list[tuple[str, Path, int]]:
    rows = []
    for d in sorted(p for p in WORK.iterdir() if p.is_dir() and not p.is_symlink()):
        cat = category(d.name)
        if cat is None:
            continue
        for f in sorted(d.rglob("*")):
            st = f.lstat()
            if f.is_symlink() or not f.is_file() or f.suffix not in ARRAY_SUFFIXES or st.st_size < MIN_BYTES:
                continue
            r = f.resolve()
            if WORK not in r.parents or category(r.relative_to(WORK).parts[0]) != cat:
                raise SystemExit(f"refused: {f} resolves to {r}")
            rows.append((cat, f, st.st_size))
    return rows


def main() -> None:
    delete = "--delete" in sys.argv[1:]
    rows = plan()
    by = {}
    for cat, f, n in rows:
        by.setdefault(cat, [0, 0])
        by[cat][0] += 1
        by[cat][1] += n
    for cat, (k, n) in by.items():
        print(f"{cat:11s} {k:4d} files {n / 2**30:6.2f} GiB")
    total = sum(n for _, _, n in rows)
    print(f"{'total':11s} {len(rows):4d} files {total / 2**30:6.2f} GiB")
    if not delete:
        for cat, f, n in rows:
            print(f"  {cat:11s} {n / 1e6:9.1f} MB  {f.relative_to(ROOT)}")
        return
    if b := busy():
        raise SystemExit("refused: a 553 script is running:\n" + "\n".join(b))
    t0 = time.time()
    tsv = WORK / "PRUNED-1002.tsv"
    with open(tsv, "w") as out:
        out.write("category\tpath\tbytes\tsha256\n")
        for cat, f, n in rows:
            out.write(f"{cat}\t{f.relative_to(ROOT)}\t{n}\t{sha256(f)}\n")
            out.flush()
    print(f"hashed {len(rows)} files in {time.time() - t0:.0f} s -> {tsv.relative_to(ROOT)}")
    for cat, f, n in rows:
        f.unlink()
    stamp = time.strftime("%Y-%m-%d %H:%MZ", time.gmtime())
    lines = [f"\n## Second prune, {stamp} (`prune553b.py --delete`; approved by the coordinator for 10-04/10-05 scratch)\n",
             f"{len(rows)} arrays ({total / 2**30:.2f} GiB), each hashed before deletion: `work/PRUNED-1002.tsv` (category, path, "
             f"bytes, sha256). Per category: " + "; ".join(f"{c} {k} files {n / 2**30:.2f} GiB" for c, (k, n) in by.items()) + ".",
             "Kept: R7b's access inputs (`r7b/`, `r7bd5000*`), the h1r7 re-pricing inputs (the five draws' controlref and "
             "callblock, `h1r7[sN]_cpm`, `h1/bulk553.npz`), the known answers' files (`h1_*`, `access553_h1`, "
             "`access553_rehearsal5`, `access553_r7bd5000-xBST2`, `h1cb*`, `h1ka`, `test553`) and every json / md / csv / log.",
             "Regenerate: arms with `run553.py <surface> <MODE> base <seed> ...` (emission is deterministic: KA-FINAL); "
             "control halves with `prep553.py chalf h1r7 <k>`; the floor sweep with `prep553.py subset r7bd<floor>` + "
             "`prep553.py controlref r7bd<floor>`; an access's caches by re-running the access (its results json is kept).\n"]
    with open(WORK / "PRUNED.md", "a") as f:
        f.write("\n".join(lines))
    print(f"deleted {len(rows)} files, {total / 2**30:.2f} GiB ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
