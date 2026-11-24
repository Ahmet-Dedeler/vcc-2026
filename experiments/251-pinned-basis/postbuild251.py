"""251 -- `scripts/def_postbuild.py`'s known-answer rows for the pinned-basis change.

Runs the changed script AND the pre-251 byte copy (`work/def_postbuild.PRE.py`) on the metas that
must not move, then (with `--rehearsal <meta>`) on the rehearsal's meta and on tampered copies of
it, each of which must FAIL on the named row:

    KA-1  222's shipped meta                              31/31 under both (nothing about it changed)
    KA-2  249's `cd4svd5` meta, --allow-config-diff cd4_deflate_k   32/32 under both
    KA-3  248's KA rebuild meta                           32/32 under both
    KA-4  the 251 rehearsal meta (D/E/F, pinned)          PASS with the new rows; the PRE script FAILS
          it on rule 10 (`cd4_deflate_basis` unallowed) -- the row the mandate's item 3 is about
    KA-5  tampered rehearsal metas: basis sha moved / mode not pinned / cd4_deflate_k added / name
          unregistered / audit pinned without the config key  -> FAIL, the named row
    KA-6  the rehearsal meta against a pre-251 BUILDER source (`--builder`), i.e. a source with no
          registry: FAIL on 'is in the builder's CD4_DEFLATION_BASIS_REGISTRY'

    uv run python experiments/251-pinned-basis/postbuild251.py [--rehearsal work/rehearsal/<name>.h5ad.meta.json]
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
NEW = ROOT / "scripts/def_postbuild.py"
PRE = WORK / "def_postbuild.PRE.py"
EF_BUNDLE = ROOT / "experiments/238-essential-free-panel/work/ef-bundle"


RUNNER = r'''
import importlib.util, sys
from pathlib import Path
script, pinned, argv = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3:]
spec = importlib.util.spec_from_file_location("postbuild_under_test", script)
m = importlib.util.module_from_spec(spec)
m.__file__ = str(pinned)                 # ROOT / BUILDER / SHIPPED_META resolve from scripts/
sys.modules["postbuild_under_test"] = m
spec.loader.exec_module(m)
sys.argv = [str(pinned)] + argv
m.main()
'''


def run(script: Path, meta: Path, *extra: str) -> tuple[int, str, str]:
    """Run a def_postbuild source (the working tree's or the byte copy) with `__file__` pinned to
    `scripts/def_postbuild.py`, so its ROOT-relative paths resolve wherever the copy lives."""
    proc = subprocess.run(["uv", "run", "python", "-c", RUNNER, str(script), str(NEW),
                           "--meta", str(meta), *extra],
                          capture_output=True, text=True, cwd=ROOT)
    out = proc.stdout + proc.stderr
    tally = next((ln.strip() for ln in out.splitlines() if ln.strip().endswith("checks pass")), "<no tally>")
    fails = [ln.strip() for ln in out.splitlines() if ln.strip().startswith("FAIL")]
    return proc.returncode, tally, "\n".join(fails)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rehearsal", type=Path, default=None)
    ap.add_argument("--cd4-root", type=Path, default=WORK / "rehearsal/cd4_compact")
    args = ap.parse_args()
    rows = []

    def row(name, ok, detail):
        rows.append({"row": name, "ok": bool(ok), "detail": detail})
        print(f"  {'PASS' if ok else 'FAIL'}  {name}: {detail}", flush=True)

    fixed = [
        ("KA-1 222 shipped meta", ROOT / "submissions/222-h1union-g028.h5ad.meta.json", [], "31/31"),
        ("KA-2 249 cd4svd5 meta (--allow-config-diff cd4_deflate_k)",
         ROOT / "submissions/249-cd4svd5-f034.h5ad.meta.json", ["--allow-config-diff", "cd4_deflate_k"], "32/32"),
        ("KA-3 248 KA rebuild meta", ROOT / "submissions/248-ka-h1union-g028.h5ad.meta.json", [], "32/32"),
    ]
    for name, meta, extra, want in fixed:
        for label, script in (("251", NEW), ("pre-251", PRE)):
            rc, tally, fails = run(script, meta, *extra)
            row(f"{name} [{label}]", rc == 0 and tally.startswith(want), f"rc={rc} {tally} {fails[:200]}")

    if args.rehearsal is not None:
        meta = args.rehearsal.resolve()
        common = ["--bundle", str(EF_BUNDLE), "--cd4-root", str(args.cd4_root.resolve()),
                  "--allow-config-diff", "whitener_id"]
        rc, tally, fails = run(NEW, meta, *common)
        row("KA-4 rehearsal meta [251]", rc == 0, f"rc={rc} {tally} {fails[:300]}")
        rc, tally, fails = run(PRE, meta, *common)
        row("KA-4 rehearsal meta [pre-251] FAILS on rule 10 (cd4_deflate_basis unallowed)",
            rc == 1 and "config differences" in fails, f"rc={rc} {tally} {fails[:300]}")
        rc, tally, fails = run(NEW, meta, *common, "--builder", str(WORK / "make_submission_pm_cd4.PRE.py"))
        row("KA-6 rehearsal meta against a pre-251 builder SOURCE (no registry) FAILS on the registry row",
            rc == 1 and "CD4_DEFLATION_BASIS_REGISTRY" in fails, f"rc={rc} {tally} {fails[:300]}")
        base = json.loads(meta.read_text())
        tampers = {
            "basis sha moved": (lambda m: m["effects"]["cd4"]["deflation"].__setitem__("basis_sha256", "0" * 64),
                                "basis_sha256"),
            "mode not pinned": (lambda m: m["effects"]["cd4"]["deflation"].__setitem__("mode", "refit"),
                                "deflation.mode == 'pinned'"),
            "cd4_deflate_k added": (lambda m: m["config"].__setitem__("cd4_deflate_k", 5),
                                    "cd4_deflate_k ABSENT"),
            "name unregistered": (lambda m: (m["config"].__setitem__("cd4_deflate_basis", "V5_nope"),
                                             m["effects"]["cd4"]["deflation"].__setitem__("basis", "V5_nope")),
                                  "CD4_DEFLATION_BASIS_REGISTRY"),
            "audit pinned, config key removed": (lambda m: m["config"].pop("cd4_deflate_basis"),
                                                 "requires config.cd4_deflate_basis"),
            "gene axis not verified": (lambda m: m["effects"]["cd4"]["deflation"].__setitem__("gene_axis_verified", False),
                                       "gene_axis_verified"),
            "V sha moved": (lambda m: m["effects"]["cd4"]["deflation"].__setitem__("basis_v_sha256", "1" * 64),
                            "basis_v_sha256"),
        }
        with tempfile.TemporaryDirectory(prefix="postbuild251-") as tmp:
            for label, (mutate, phrase) in tampers.items():
                m = json.loads(json.dumps(base))
                mutate(m)
                p = Path(tmp) / f"{label.replace(' ', '_').replace(',', '')}.meta.json"
                p.write_text(json.dumps(m))
                rc, tally, fails = run(NEW, p, *common)
                row(f"KA-5 tampered ({label}) FAILS on '{phrase}'", rc == 1 and phrase in fails,
                    f"rc={rc} {tally} {fails[:260]}")
    n_fail = sum(not r["ok"] for r in rows)
    (WORK / ("postbuild251.json" if args.rehearsal else "postbuild251-pre.json")).write_text(json.dumps(rows, indent=1))
    print(f"\nPOSTBUILD KA {len(rows) - n_fail}/{len(rows)} {'PASS' if not n_fail else 'FAIL'}")
    raise SystemExit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
