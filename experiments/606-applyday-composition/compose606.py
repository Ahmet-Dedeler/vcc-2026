#!/usr/bin/env python3
"""606: one composition test for EVERY pending patch, in the coordinator's apply order, on a scratch copy of the tree.

Several UNAPPLIED patches edit the same files (def_day.sh: 641, 598, 612dd, 642, 644; def_sameline.py: 642, 604, 644; the builder:
598, 612, 630). Each patch's own test composes only with its predecessors, so an anchor clash between two of them shows
up only on the apply day. This applies them in order, running each one's --check first:
  10-05 bundle  : 538 -> 554 -> 590  (test598.canon_tree)  then 641 (10-05 checklist step 6b)
  10-06 apply day: 598 -> 612 -> 612dd -> 630 -> 630dd -> 642 (NEXT-SLOT-DECISION "10-06 APPLY DAY"; 612dd per 21dd4c; 630dd 4d0207)
  after         : 604 (identity prior source, council 600) -> 644 (DEF_DAY_SAMELINE_SKIP; anchors on 642, applies LAST)
Then it runs `bash -n` on def_day.sh and py_compile on every touched Python file, and reports the final shas.
A patch whose --check says it is ALREADY in the tree (after the 10-05 / 10-06 applies, HEAD carries it) is skipped, not
failed (B2, 650's tabletop); a PARTIALLY patched tree still fails. Nothing in
the live repo is written (the live files are hashed before and after). Use --order to try another sequence.
"""
from __future__ import annotations

import argparse
import hashlib
import py_compile
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PY = sys.executable
sys.path.insert(0, str(ROOT / "experiments/598-tercilec-knob"))
import test598 as T598  # noqa: E402

PATCHES = {
    "641": ("experiments/641-depth-row/patch_depth641.py", ["--check"], []),
    "598": ("experiments/598-tercilec-knob/patch598.py", ["--check"], []),
    "612": ("experiments/612-key-gated-hook/patch612.py", ["--check"], []),
    "612dd": ("experiments/612-key-gated-hook/patch612dd.py", ["--check"], []),
    "630": ("experiments/630-rank-neutral-up/patch_builder630.py", ["--check"], []),
    "630dd": ("experiments/630-rank-neutral-up/patch630dd.py", ["--check"], []),
    "642": ("experiments/642-hek293-contingency/patch_hek293_642.py", ["--check"], []),
    "604": ("experiments/604-identity-prior/patch604.py", ["--check"], []),
    "644": ("experiments/644-runbook-amendments/patch_skip644.py", ["--check"], []),
    "604f": ("experiments/604-identity-prior/patch604f.py", ["--check"], []),   # 10-08+ (after the freeze); not in DEFAULT_ORDER
}
DEFAULT_ORDER = ["641", "598", "612", "612dd", "630", "630dd", "642", "604", "644"]
WATCH = ["scripts/def_day.sh", "scripts/def_sameline.py", "scripts/def_identity.py", "scripts/def_day_helpers.py",
         "scripts/def_pipeline.py", "scripts/october_recipe.py", "scripts/make_submission_pm_cd4.py",
         "scripts/def_preflight.py", "scripts/def_postbuild.py", "experiments/524-october-recipe-dr/gate524.py"]


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else "-"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--order", default=",".join(DEFAULT_ORDER))
    a = ap.parse_args()
    order = a.order.split(",")
    live = {p: sha(ROOT / p) for p in WATCH}
    ok_all = True
    with tempfile.TemporaryDirectory(prefix="compose606-") as td:
        t = Path(td) / "tree"
        T598.canon_tree(t)
        print(f"canonical post-bundle tree: def_day.sh {sha(t / 'scripts/def_day.sh')[:12]}  def_sameline.py "
              f"{sha(t / 'scripts/def_sameline.py')[:12]}  builder {sha(t / 'scripts/make_submission_pm_cd4.py')[:12]}", flush=True)
        # files some patches need beyond test598's CANON set
        (t / "data/identity_refs").mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / "data/identity_refs/hek293_nonflex_ref642.npz", t / "data/identity_refs/hek293_nonflex_ref642.npz")
        (t / "src/vcc2026").mkdir(parents=True, exist_ok=True)
        for f in (ROOT / "src/vcc2026").glob("*.py"):
            shutil.copy2(f, t / "src/vcc2026" / f.name)
        for name in order:
            script, chk, extra = PATCHES[name]
            c = subprocess.run([PY, str(ROOT / script), "--root", str(t), *chk], capture_output=True, text=True)
            co = (c.stdout + c.stderr).lower()
            if c.returncode != 0 and "already" in co and "partially" not in co:
                # B2 (650's tabletop): after the 10-05 apply (or 10-06's), HEAD itself carries this patch: skip, don't fail
                print(f"SKIP {name}: already in the tree (on HEAD): {(c.stdout + c.stderr).strip().splitlines()[-1][:160]}", flush=True)
                continue
            if c.returncode != 0:
                ok_all = False
                print(f"FAIL {name} --check rc={c.returncode}: {(c.stdout + c.stderr).strip()[-600:]}", flush=True)
                continue
            q = subprocess.run([PY, str(ROOT / script), "--root", str(t), *extra], capture_output=True, text=True)
            if q.returncode != 0:
                ok_all = False
                print(f"FAIL {name} apply rc={q.returncode}: {(q.stdout + q.stderr).strip()[-600:]}", flush=True)
                continue
            print(f"PASS {name} (check + apply): {(q.stdout + q.stderr).strip().splitlines()[-1][:200] if (q.stdout + q.stderr).strip() else ''}", flush=True)
        b = subprocess.run(["bash", "-n", str(t / "scripts/def_day.sh")], capture_output=True, text=True)
        print(f"{'PASS' if b.returncode == 0 else 'FAIL'} bash -n def_day.sh {b.stderr.strip()[:300]}", flush=True)
        ok_all &= b.returncode == 0
        for p in sorted(set(WATCH) | {"scripts/def_identity_prior.py", "src/vcc2026/rank_up.py"}):
            f = t / p
            if f.suffix == ".py" and f.exists():
                try:
                    py_compile.compile(str(f), doraise=True)
                except py_compile.PyCompileError as e:
                    ok_all = False
                    print(f"FAIL py_compile {p}: {e}", flush=True)
        print("final shas: " + "  ".join(f"{p.split('/')[-1]} {sha(t / p)[:12]}" for p in WATCH + ["scripts/def_identity_prior.py", "src/vcc2026/rank_up.py"]), flush=True)
    after = {p: sha(ROOT / p) for p in WATCH}
    if after != live:
        ok_all = False
        print("FAIL the live repo changed during the test", flush=True)
    print("ALL PASS" if ok_all else "SOME FAILED", flush=True)
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
