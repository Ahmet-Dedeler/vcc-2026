"""Re-pin the builder's PANEL-DEPENDENT frozen inputs for a new round: the whitening operator
and the three CD4 banks.

Experiment 199 found that two of the builder's "ours, never moves" frozen inputs are in fact
functions of the 300-target panel, and that both are pinned by sha256 IN THE BUILDER'S SOURCE
(`scripts/make_submission_pm_cd4.py`: `WHITENER_REGISTRY` and the `cd4/...` rows of
`FROZEN_INPUT_SHA256`), not in `scripts/frozen_inputs.json`:

* the operator's axis is `union(donor coverage) - panel` (`build_union_whitener.py`), so the
  A/B/C operator pairs 98.8% of its coordinates with the WRONG gene on a new panel -- and
  every existing check (n_cols, evec shape, file sha256) passes, because both panels give
  17,942 columns. This is the one silent failure on the Oct-22 path.
* the CD4 compact banks hold only the panel's rows (283/287/284 of A/B/C), so on a new panel
  `cd4_complete_targets` collapses to 0 and `VCC_CD4_TARGET_F` raises.

So on the day both artifacts are re-made for the new panel, and the builder then REJECTS them
("frozen model input hash drift") until their hashes are re-pinned. This script performs that
re-pin as a deliberate, diffed, reviewable source edit:

  * adds (or replaces) one `WHITENER_REGISTRY[<id>]` entry with the new operator's (npz, json)
    sha256 pair -- ADDITIVE, so `--whitener-id union_pruned` builds are untouched;
  * replaces the 15 `cd4/<bank>/<file>` hashes with the new banks' -- these are keyed by bank
    name, not by path, so they cannot be additive.

It never touches anything else, refuses if any anchor is not found exactly once, verifies the
result by re-reading the file, and prints the unified diff. `--restore --record <out.json>` puts
the byte copy saved at pin time back (used by the rehearsal; on the day you COMMIT the pin
instead, so the meta's recorded builder hash is reproducible). It does not restore from git
HEAD by default because HEAD can move under you: during the 2026-09-03 rehearsal a sweep
commit landed while the pin was live and would have made a HEAD restore a no-op.

    # dry: print the diff only
    uv run python scripts/def_pin_panel_artifacts.py \
        --whitener <new>.npz --whitener-metadata <new>.json --whitener-id union_pruned_def \
        --cd4-root <new cd4_compact>
    # apply
    ... --apply
    # rehearsal only: byte-restore the builder to the copy saved at pin time
    uv run python scripts/def_pin_panel_artifacts.py --restore --record <out.json>

Rule 11: every pin written here is a path AND a sha256, printed and saved to --out.
"""

from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts/make_submission_pm_cd4.py"
CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
CD4_FILES = ("delta.npy", "b_ctrl.npy", "genes.txt", "perts.txt", "meta.json")
HEX64 = r"[0-9a-f]{64}"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def head_blob() -> bytes:
    return subprocess.run(["git", "show", f"HEAD:{BUILDER.relative_to(ROOT)}"], cwd=ROOT,
                          capture_output=True, check=True).stdout


def rewrite(src: str, whitener_id: str, npz_sha: str, json_sha: str,
            cd4_sha: dict[str, str], pert_scale: tuple[str, str, str] | None = None) -> str:
    out = src
    # --- 432: 410's per-target key for THIS panel, ADDITIVE in PERT_SCALE_ARTIFACTS -------------
    # (id, repo-relative path, sha256). The artifact is panel- AND operator-specific (DEF-RUNBOOK 410
    # addition) and the builder accepts it only through this source-pinned table.
    if pert_scale is not None:
        ps_id, ps_rel, ps_sha = pert_scale
        anchor = "PERT_SCALE_ARTIFACTS = {\n"
        if out.count(anchor) != 1:
            raise SystemExit(f"PERT_SCALE_ARTIFACTS anchor found {out.count(anchor)} times; refusing")
        ps_entry = (f'    # Panel-dependent 410 key re-fit for this round (scripts/def_pin_panel_artifacts.py).\n'
                    f'    "{ps_id}": (\n        ROOT / "{ps_rel}",\n        "{ps_sha}",\n    ),\n')
        ps_existing = re.compile(r'    # Panel-dependent 410 key[^\n]*\n    "%s": \(\n[^\n]*\n\s+"%s",\n    \),\n' % (
            re.escape(ps_id), HEX64))
        if ps_existing.search(out):
            out = ps_existing.sub(lambda _: ps_entry, out, count=1)
        elif f'"{ps_id}": (' in out:
            raise SystemExit(f"PERT_SCALE_ARTIFACTS already has a hand-written {ps_id!r} entry; refusing")
        else:
            out = out.replace(anchor, anchor + ps_entry, 1)
    # --- CD4 rows: exactly one anchor per key -------------------------------------------
    for key, digest in cd4_sha.items():
        pat = re.compile(r'("%s": \(\n\s+")%s(")' % (re.escape(key), HEX64))
        n = len(pat.findall(out))
        if n != 1:
            raise SystemExit(f"anchor for {key!r} found {n} times in {BUILDER.name}; refusing")
        out = pat.sub(lambda m: m.group(1) + digest + m.group(2), out)
    # --- the registry: replace the id if present, else insert before the closing brace ---
    reg = re.compile(r"WHITENER_REGISTRY: dict\[str, tuple\[str, str\]\] = \{\n(.*?)\n\}\n", re.S)
    m = reg.search(out)
    if not m:
        raise SystemExit("WHITENER_REGISTRY block not found; refusing")
    body = m.group(1)
    entry = (f'    # Panel-dependent operator re-fit for this round (scripts/def_pin_panel_artifacts.py).\n'
             f'    "{whitener_id}": (\n        "{npz_sha}",\n        "{json_sha}",\n    ),')
    existing = re.compile(r'(    # [^\n]*\n)?    "%s": \(\n\s+"%s",\n\s+"%s",\n    \),' % (
        re.escape(whitener_id), HEX64, HEX64))
    if existing.search(body):
        body = existing.sub(lambda _: entry, body, count=1)
    else:
        body = body + "\n" + entry
    out = out[:m.start(1)] + body + out[m.end(1):]
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--whitener", type=Path)
    ap.add_argument("--whitener-metadata", type=Path)
    ap.add_argument("--whitener-id", default="union_pruned_def")
    ap.add_argument("--cd4-root", type=Path)
    ap.add_argument("--apply", action="store_true", help="write the edit (default: diff only)")
    ap.add_argument("--restore", action="store_true",
                    help="byte-restore the builder to the copy saved at pin time (--record) and exit")
    ap.add_argument("--record", type=Path, default=None,
                    help="the --out JSON of the pin to undo; its sibling .builder-before.py is written back. "
                         "Without it, --restore falls back to git HEAD (NOT safe if HEAD moved since the pin)")
    ap.add_argument("--out", type=Path, default=None, help="JSON record of every pin written")
    ap.add_argument("--pert-scale", type=Path, default=None,
                    help="432: the day's 410 key artifact (key432.py on this bundle + the re-fit operator)")
    ap.add_argument("--pert-scale-id", default="wshare_lo0.3_hi3_def",
                    help="registry id the build then selects with VCC_PERT_SCALE (lower case)")
    args = ap.parse_args()

    before = BUILDER.read_bytes()
    if args.restore:
        if args.record:
            rec = json.loads(args.record.read_text())
            copy = args.record.with_suffix(".builder-before.py")
            blob = copy.read_bytes()
            if hashlib.sha256(blob).hexdigest() != rec["builder_sha256_before"]:
                raise SystemExit(f"{copy} does not hash to the recorded builder_sha256_before; refusing")
            source = f"{copy.name} (sha256 {rec['builder_sha256_before'][:16]}…)"
        else:
            blob = head_blob()
            source = "git HEAD"
        BUILDER.write_bytes(blob)
        after = BUILDER.read_bytes()
        print(f"restored {BUILDER.relative_to(ROOT)} from {source}: "
              f"{hashlib.sha256(before).hexdigest()[:16]}… -> {hashlib.sha256(after).hexdigest()[:16]}… "
              f"({'identical' if after == blob else 'MISMATCH'})")
        raise SystemExit(0 if after == blob else 1)

    for name in ("whitener", "whitener_metadata", "cd4_root"):
        if getattr(args, name) is None:
            ap.error(f"--{name.replace('_', '-')} is required unless --restore")
    if not re.fullmatch(r"[A-Za-z0-9_\-]+", args.whitener_id):
        ap.error("--whitener-id must be [A-Za-z0-9_-]+")
    whitener = args.whitener.resolve()
    wmeta = args.whitener_metadata.resolve()
    cd4_root = args.cd4_root.resolve()
    npz_sha, json_sha = sha256(whitener), sha256(wmeta)
    meta = json.loads(wmeta.read_text())
    if meta.get("sha256") not in (None, npz_sha):
        raise SystemExit(f"{wmeta.name} records sha256 {meta.get('sha256')} but the npz hashes to {npz_sha}")
    cd4_sha: dict[str, str] = {}
    for bank in CD4_NAMES:
        for f in CD4_FILES:
            p = cd4_root / bank / f
            if not p.is_file():
                raise SystemExit(f"missing CD4 bank file {p}")
            cd4_sha[f"cd4/{bank}/{f}"] = sha256(p)

    pert_scale = None
    if args.pert_scale is not None:
        ps = args.pert_scale.resolve()
        if not re.fullmatch(r"[a-z0-9_.\-]+", args.pert_scale_id):
            ap.error("--pert-scale-id must be lower case [a-z0-9_.-]+ (the builder lower-cases VCC_PERT_SCALE)")
        pert_scale = (args.pert_scale_id, str(ps.relative_to(ROOT)), sha256(ps))
    src = before.decode()
    new = rewrite(src, args.whitener_id, npz_sha, json_sha, cd4_sha, pert_scale)
    diff = list(difflib.unified_diff(src.splitlines(keepends=True), new.splitlines(keepends=True),
                                     fromfile=f"a/{BUILDER.relative_to(ROOT)}",
                                     tofile=f"b/{BUILDER.relative_to(ROOT)}"))
    sys.stdout.writelines(diff)
    changed_lines = sum(1 for d in diff if d.startswith("+") and not d.startswith("+++"))
    record = {
        "builder": str(BUILDER.relative_to(ROOT)),
        "builder_sha256_before": hashlib.sha256(before).hexdigest(),
        "builder_sha256_after": hashlib.sha256(new.encode()).hexdigest(),
        "whitener_id": args.whitener_id,
        "whitener": {"path": str(whitener), "sha256": npz_sha},
        "whitener_metadata": {"path": str(wmeta), "sha256": json_sha},
        "cd4_root": str(cd4_root),
        "cd4": cd4_sha,
        "pert_scale": None if pert_scale is None else dict(zip(("id", "path", "sha256"), pert_scale)),
        "added_lines": changed_lines,
        "applied": bool(args.apply),
    }
    print(f"\n{changed_lines} lines added/changed; registry id {args.whitener_id!r}; "
          f"operator {npz_sha[:16]}…; {len(cd4_sha)} CD4 hashes")
    # Record and before-image first, then replace. A kill used to land between
    # tmp.replace and these writes: pinned builder, no record, runbook restore
    # cannot run (526). Restoring a record whose replace has not happened yet
    # writes the before-image back onto the still-original file.
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        if args.apply:
            copy = args.out.with_suffix(".builder-before.py")
            copy.write_bytes(before)
        args.out.write_text(json.dumps(record, indent=2) + "\n")
        print(f"wrote {args.out}")
    if args.apply:
        compile(new, str(BUILDER), "exec")                 # must still parse
        tmp = BUILDER.with_suffix(".py.tmp")
        tmp.write_text(new)
        tmp.replace(BUILDER)
        check = BUILDER.read_text()
        for digest in [npz_sha, json_sha, *cd4_sha.values(), *([pert_scale[2]] if pert_scale else [])]:
            if digest not in check:
                raise SystemExit(f"post-write verification failed: {digest[:16]}… not in file")
        print(f"applied -> {BUILDER.relative_to(ROOT)} sha256 {record['builder_sha256_after']}")
        print("next: git diff scripts/make_submission_pm_cd4.py  (review), then COMMIT it on the day")
        if args.out:
            print(f"saved the pre-pin builder to {args.out.with_suffix('.builder-before.py')} (restore with --restore --record {args.out})")
    else:
        print("dry run: nothing written (pass --apply)")


if __name__ == "__main__":
    main()
