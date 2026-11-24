"""432: run the production builder with ONE extra `PERT_SCALE_ARTIFACTS` entry, in memory. The file is never written.

Why this exists. On a new panel, 410's per-target key (`VCC_PERT_SCALE=wshare_lo0.3_hi3`) must be rebuilt for the new
panel and operator (DEF-RUNBOOK, 410 addition), and the builder accepts a key artifact only through its source-pinned
table `PERT_SCALE_ARTIFACTS` (path + sha256). `scripts/def_pin_panel_artifacts.py` pins only the operator and the CD4
banks, so the day needs one more source edit. The 432 mandate says: if the key needs a builder edit, do not make it,
document it. This wrapper documents it by EMULATING exactly that edit in memory, so everything downstream of the key
(ballast, sharpening, dump/load, re-level, MDE, zero repair, zs) can still be rehearsed on D/E/F while
`scripts/make_submission_pm_cd4.py` stays byte-identical on disk.

What it does:
  1. reads the builder bytes and REFUSES unless their sha256 is `--expect-builder-sha` (the driver's pin, i.e. the
     pin record's `builder_sha256_after`);
  2. inserts ONE entry at the top of `PERT_SCALE_ARTIFACTS` (anchor must occur exactly once):
         "<id>": (Path("<artifact>"), "<sha256>"),
     and refuses if the id already exists;
  3. writes the patched source, its sha256 and the unified diff next to `--record` (provenance), and
  4. executes the patched code as `__main__` with `__file__` = the real builder path (so ROOT resolves as in
     production) and `sys.argv` = the builder path + the arguments after `--`.

Caveat, stated in RESULTS.md: the builder hashes `Path(__file__)` for `hashes.implementation`, so the meta records the
on-disk (pinned) sha, not the patched source's. The record written here carries the patched sha.

usage: uv run python run_builder432.py --builder scripts/make_submission_pm_cd4.py --id wshare_lo0.3_hi3_def
           --artifact <npz> --expect-builder-sha <sha> --record <work>/builder432.json -- <builder args...>
"""
from __future__ import annotations

import argparse
import builtins
import difflib
import hashlib
import json
import sys
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILDER = ROOT / "scripts/make_submission_pm_cd4.py"
ANCHOR = "PERT_SCALE_ARTIFACTS = {\n"


def main() -> None:
    argv = sys.argv[1:]
    if "--" not in argv:
        raise SystemExit("usage: run_builder432.py [opts] -- <builder args>")
    cut = argv.index("--")
    ap = argparse.ArgumentParser()
    ap.add_argument("--id", required=True)
    ap.add_argument("--artifact", type=Path, required=True)
    ap.add_argument("--expect-builder-sha", required=True)
    ap.add_argument("--record", type=Path, required=True)
    # the builder path ON THE COMMAND LINE, so def_pipeline.other_builders() (and every gate that uses it) sees this
    # process as a running build; it must name the production builder
    ap.add_argument("--builder", type=Path, required=True)
    a = ap.parse_args(argv[:cut])
    if a.builder.resolve() != BUILDER:
        raise SystemExit(f"--builder must be {BUILDER}")
    builder_args = argv[cut + 1:]

    raw = BUILDER.read_bytes()
    got = hashlib.sha256(raw).hexdigest()
    if got != a.expect_builder_sha:
        raise SystemExit(f"builder on disk {got} != expected {a.expect_builder_sha}; refusing")
    src = raw.decode()
    if src.count(ANCHOR) != 1:
        raise SystemExit(f"anchor {ANCHOR!r} found {src.count(ANCHOR)} times; refusing")
    if f'"{a.id}": (' in src:
        raise SystemExit(f"registry id {a.id!r} already in the builder; refusing")
    art = a.artifact.resolve()
    art_sha = hashlib.sha256(art.read_bytes()).hexdigest()
    entry = (f"    # 432 D/E/F rehearsal: 410's key re-fit for this panel (key432.py). The day-of edit, emulated in memory.\n"
             f"    \"{a.id}\": (\n        Path(\"{art}\"),\n        \"{art_sha}\",\n    ),\n")
    patched = src.replace(ANCHOR, ANCHOR + entry, 1)
    psha = hashlib.sha256(patched.encode()).hexdigest()
    a.record.parent.mkdir(parents=True, exist_ok=True)
    a.record.with_suffix(".patched.py").write_text(patched)
    diff = "".join(difflib.unified_diff(src.splitlines(keepends=True), patched.splitlines(keepends=True),
                                        fromfile="a/scripts/make_submission_pm_cd4.py",
                                        tofile="b/scripts/make_submission_pm_cd4.py"))
    a.record.with_suffix(".diff").write_text(diff)
    rec = dict(builder=str(BUILDER), builder_sha256_on_disk=got, patched_sha256=psha, registry_id=a.id,
               artifact=str(art), artifact_sha256=art_sha, builder_args=builder_args)
    a.record.write_text(json.dumps(rec, indent=1) + "\n")
    print(f"[run_builder432] builder {got[:16]} + registry {a.id!r} -> {art.name} ({art_sha[:16]}); "
          f"patched source sha {psha[:16]} (record {a.record.name})", flush=True)
    sys.stdout.write(diff)
    sys.stdout.flush()

    code = compile(patched, str(BUILDER), "exec")
    sys.argv = [str(BUILDER), *builder_args]
    sys.path[0] = str(BUILDER.parent)            # as when the script is run directly
    mod = types.ModuleType("__main__")           # a real module registered as __main__ (dataclasses/pickle look it up)
    mod.__file__ = str(BUILDER)
    mod.__builtins__ = builtins
    sys.modules["__main__"] = mod
    exec(code, mod.__dict__)


if __name__ == "__main__":
    main()
