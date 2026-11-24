#!/usr/bin/env python3
"""Post-deadline disk cleanup for the VCC repo. Stdlib only.

Without flags it is a dry run: it lists what would go, prints the sizes and writes the reviewed list
to work/delete-list.txt.gz. With --delete it permanently deletes exactly that reviewed list, after
checking the guards below. It only ever touches files git ignores, so nothing tracked can go.

Deleted: ignored files under experiments/, submissions/, data/external/, and the two duplicate
D/E/F download folders.
Kept, on top of everything git tracks: the two final packages, 692's build inputs (not its h5ads),
811's banks, and small text results. data/donors, data/harness, data/challenge*, data/leaderboard,
data/board-entries and .venv are outside the deleted trees.
"""
import argparse
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

SCOPES = ["experiments", "submissions", "data/external", "data/def-download", "data/challenge-def.unzipped"]
FINAL_PACKAGES = ["submissions/x2def.vcc", "submissions/r0def.vcc"]
MUST_SURVIVE = FINAL_PACKAGES + ["data/donors", "data/harness", "data/challenge", "data/challenge-def",
                                 "data/leaderboard", "data/board-entries", ".venv"]
LIST_REL = "experiments/723-post-final-cleanup/work/delete-list.txt.gz"
TEXT_EXT = {".json", ".jsonl", ".md", ".txt", ".tsv", ".csv", ".log", ".py", ".sh", ".yaml", ".yml",
            ".toml", ".ipynb", ".html", ".r", ".nohup"}
SMALL_TEXT_BYTES = 5_000_000
VENDORED = re.compile(r"(^|/)(site-packages|\.venv|venv|node_modules|__pycache__|dist-info|\.cache|lib/python[0-9.]*)(/|$)")


def git(root, *args):
    return subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True).stdout


def z_lines(raw):
    return [p for p in raw.decode("utf8", "surrogateescape").split("\0") if p]


def is_kept(path, size):
    """True for an ignored file the cleanup must leave on disk."""
    if path in FINAL_PACKAGES or path.startswith("experiments/723-post-final-cleanup/"):
        return True
    if path.startswith("experiments/692-def-r0/") and not path.endswith(".h5ad"):
        return True
    if path.startswith("experiments/811-lolo-model/work/banks/"):
        return True
    # Small text outside vendored packages is a result or a record, not a bulk artifact.
    if os.path.splitext(path)[1].lower() in TEXT_EXT and size < SMALL_TEXT_BYTES and not VENDORED.search(path):
        return True
    return False


def build_list(root):
    """Return the sorted delete list; a trailing / marks a nested git repo removed as a whole."""
    tracked = set(z_lines(git(root, "ls-files", "-z")))
    out = []
    for p in z_lines(git(root, "ls-files", "-z", "-o", "-i", "--exclude-standard", "--", *SCOPES)):
        if p.endswith("/"):
            # git reports a nested repo as one directory entry; only data/external holds any.
            if not p.startswith("data/external/") or any(t.startswith(p) for t in tracked):
                raise SystemExit(f"REFUSED: unexpected nested repo {p}")
            out.append(p)
            continue
        if p in tracked:
            raise SystemExit(f"REFUSED: {p} is both ignored and tracked")
        st = os.lstat(root / p)
        if not is_kept(p, st.st_size):
            out.append(p)
    return sorted(out)


def disk_bytes(root, path):
    full = root / path
    if path.endswith("/"):
        return sum(os.lstat(os.path.join(d, f)).st_blocks * 512 for d, _, fs in os.walk(full) for f in fs)
    return os.lstat(full).st_blocks * 512


def encode(paths):
    data = "\n".join(paths).encode("utf8", "surrogateescape") + b"\n"
    return gzip.compress(data, mtime=0)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1 << 24):
            h.update(chunk)
    return h.hexdigest()


def summary(root, paths):
    groups = Counter()
    exps = Counter()
    total = 0
    for p in paths:
        size = disk_bytes(root, p)
        total += size
        parts = p.split("/")
        groups["/".join(parts[:2]) if parts[0] == "data" else parts[0]] += size
        if parts[0] == "experiments":
            exps[parts[1]] += size
    print(f"{len(paths):,} entries, {total / 2**30:.1f} GiB by file size (APFS clones can make the real gain smaller)")
    for k, v in groups.most_common():
        print(f"  {v / 2**30:8.1f} GiB  {k}")
    print("  largest experiments:")
    for k, v in exps.most_common(8):
        print(f"  {v / 2**30:8.1f} GiB  {k}")


def guards(root, paths):
    """Every condition that must hold before a permanent delete. Returns the failures."""
    bad = []
    if git(root, "status", "--porcelain").strip():
        bad.append("git working tree is not clean: commit first")
    git(root, "fetch", "-q", "origin")
    if git(root, "rev-parse", "HEAD") != git(root, "rev-parse", "origin/master"):
        bad.append("HEAD is not origin/master: push first")
    ledger = json.loads((root / "submissions/SUBMITTED_LEDGER.json").read_text())
    for pkg in FINAL_PACKAGES:
        full = root / pkg
        if not full.exists():
            bad.append(f"{pkg} is missing")
        elif sha256(full) not in ledger or ledger[sha256(full)].get("package") != pkg:
            bad.append(f"{pkg} does not match SUBMITTED_LEDGER")
    for keep in MUST_SURVIVE:
        if not (root / keep).exists():
            bad.append(f"{keep} is missing")
        if any(p == keep or p.startswith(keep + "/") for p in paths):
            bad.append(f"the list touches kept path {keep}")
    return bad


def remove_empty_dirs(top):
    for d, _, _ in os.walk(top, topdown=False):
        try:
            os.rmdir(d)
        except OSError:
            pass


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--delete", action="store_true", help="permanently delete the reviewed list")
    ap.add_argument("--check", action="store_true", help="run every guard and the list comparison, delete nothing")
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[2]))
    args = ap.parse_args()
    root = Path(args.root).resolve()
    paths = build_list(root)
    list_path = root / LIST_REL
    summary(root, paths)

    if args.check:
        same = list_path.exists() and list_path.read_bytes() == encode(paths)
        bad = guards(root, paths) + ([] if same else ["the files on disk no longer match the reviewed list"])
        print("CHECK PASS: --delete would run. Nothing deleted." if not bad else "CHECK FAIL:\n  " + "\n  ".join(bad))
        sys.exit(1 if bad else 0)

    if not args.delete:
        list_path.parent.mkdir(parents=True, exist_ok=True)
        list_path.write_bytes(encode(paths))
        print(f"DRY RUN: nothing deleted. List written to {LIST_REL}. Run with --delete to apply it.")
        return

    if not list_path.exists() or list_path.read_bytes() != encode(paths):
        sys.exit("REFUSED: the files on disk no longer match the reviewed list. Re-run the dry run and review it.")
    bad = guards(root, paths)
    if bad:
        sys.exit("REFUSED:\n  " + "\n  ".join(bad))

    before = shutil.disk_usage(root).free
    failed = 0
    for p in paths:
        full = root / p
        try:
            if p.endswith("/"):
                shutil.rmtree(full)
            else:
                full.unlink()
        except OSError as e:
            failed += 1
            print(f"could not remove {p}: {e}", file=sys.stderr)
    for top in ["experiments", "data/external", "data/def-download", "data/challenge-def.unzipped"]:
        if (root / top).exists():
            remove_empty_dirs(root / top)
    freed = shutil.disk_usage(root).free - before
    print(f"DONE: {len(paths) - failed:,} removed, {failed} failed, {freed / 2**30:.1f} GiB freed, "
          f"{shutil.disk_usage(root).free / 2**30:.1f} GiB now free")
    status = git(root, "status", "--porcelain").decode()
    print("git status after delete:", "clean" if not status.strip() else f"{len(status.splitlines())} changes, investigate")


if __name__ == "__main__":
    main()
