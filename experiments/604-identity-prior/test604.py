#!/usr/bin/env python3
"""604 tests: patch604 on the def_sameline.py the apply days will actually hold, in two variants:
  BUNDLE  = the canonical post-bundle tree (test598.canon_tree: HEAD + 538 -> 554 -> 590; def_sameline 22c6f0f0);
  FULL    = BUNDLE + 641 (depth row) + 642 (HEK293 contingency, which also edits def_sameline.py), the 10-06 state.
Each run tree is a scratch dir whose data/ and experiments/ are read-only symlinks into the live repo. Nothing in the live
repo is written (the live def_sameline.py and def_day.sh are hashed at start and end).

  T1  --check passes on the canonical file; a second application is refused.
  T2  UNSET path unchanged: the patched def_sameline (no --programme-source) writes the same arrays as the canonical one
      on the relabelled A/B/C bundle (180's def-bundle; D/E/F = A/B/C data; identity = 457's run: JURKAT/HELA/CAL33).
  T3  the existing bank path unchanged: --programme-source h1_2025_full arrays identical before / after the patch.
  T4  KA lam 0: --programme-source identity:0 gives programme_<c> == CLEAN4 psign and pdown_<c> == CLEAN4 P on its
      covered genes, for every context.
  T5  day-of == board: --programme-source identity:0.5 on the relabelled bundle gives programme_D/E/F and pdown_D/E/F
      equal to build604's board artifact's A/B/C (work/ident604v3_flex_lam0.5_ABC.npz): the guard keeps exactly one
      reference per identity call, so D sees {E, F, h1} = A's {B, C, h1}.
  T7  `identity` (bare; lambda 0.5) == identity:0.5, the source def_day's hookbank form passes (hookbank:identity:<tau>).
  T8  `identity@1` (colon-free lambda) == the lambda 1 board artifact.
  T6  a matched context keeps its same-line rows under the identity source (D = Jurkat: rows_D / signs_D / targets_D /
      w_D unchanged vs the no-programme-source run).
usage: uv run python experiments/604-identity-prior/test604.py   (needs work/flex_refs604.npz and the lam 0.5 artifact)
"""
from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PY = sys.executable
BUNDLE = ROOT / "experiments/180-def-readiness/work/def-bundle"
IDENT = ROOT / "experiments/457-cd4-half-def/work/run/identity.json"
BOARD = HERE / "work/ident604v3_flex_lam0.5_ABC.npz"
BOARD1 = HERE / "work/ident604v3_flex_lam1_ABC.npz"
C4 = ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz"
RESULTS: list[tuple[str, bool, str]] = []


def rec(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    print(f"{'PASS' if ok else 'FAIL'} {name}{': ' + detail if detail else ''}", flush=True)


sys.path.insert(0, str(ROOT / "experiments/598-tercilec-knob"))
import test598 as T598  # noqa: E402  (canon_tree: read-only reuse of 598's canonical post-bundle construction)


def variant_samelines() -> dict[str, bytes]:
    """{'BUNDLE': post-bundle def_sameline.py, 'FULL': the same after 641 + 642} built in scratch trees."""
    out = {}
    with tempfile.TemporaryDirectory(prefix="canon604-") as td:
        t = Path(td) / "canon"
        T598.canon_tree(t)
        out["BUNDLE"] = (t / "scripts/def_sameline.py").read_bytes()
        # 642 checks its pinned non-Flex HEK293 reference inside the tree (data/identity_refs/), as its own test does
        (t / "data/identity_refs").mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / "data/identity_refs/hek293_nonflex_ref642.npz", t / "data/identity_refs/hek293_nonflex_ref642.npz")
        for script in ("experiments/641-depth-row/patch_depth641.py", "experiments/642-hek293-contingency/patch_hek293_642.py"):
            c = subprocess.run([PY, str(ROOT / script), "--root", str(t), "--check"], capture_output=True, text=True)
            co = (c.stdout + c.stderr).lower()
            if c.returncode != 0 and "already" in co and "partially" not in co:
                continue   # B2: already on HEAD (after the 10-05 / 10-06 applies)
            q = subprocess.run([PY, str(ROOT / script), "--root", str(t)], capture_output=True, text=True)
            if c.returncode != 0 or q.returncode != 0:
                raise SystemExit(f"{script} failed: {c.stdout[-300:]} {c.stderr[-300:]} {q.stdout[-300:]} {q.stderr[-300:]}")
        out["FULL"] = (t / "scripts/def_sameline.py").read_bytes()
    return out


def run_sameline(tree: Path, out: Path, *extra) -> dict:
    cmd = [PY, str(tree / "scripts/def_sameline.py"), "--bundle", str(BUNDLE), "--identity", str(IDENT), "--out", str(out), *extra]
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tree))
    if p.returncode != 0:
        raise SystemExit(f"def_sameline failed ({' '.join(extra)}): {p.stdout[-600:]} {p.stderr[-900:]}")
    with np.load(out, allow_pickle=False) as z:
        return {k: np.asarray(z[k]) for k in z.files if k != "provenance"}


def same_arrays(a: dict, b: dict) -> tuple[bool, str]:
    if set(a) != set(b):
        return False, f"keys differ: {sorted(set(a) ^ set(b))}"
    for k in a:
        x, y = a[k], b[k]
        if x.shape != y.shape or x.dtype != y.dtype:
            return False, f"{k}: shape/dtype {x.shape}/{x.dtype} vs {y.shape}/{y.dtype}"
        if x.dtype.kind == "f":
            if not np.array_equal(x, y, equal_nan=True):
                return False, f"{k}: values differ"
        elif not np.array_equal(x, y):
            return False, f"{k}: values differ"
    return True, f"{len(a)} arrays identical"


def main() -> None:
    live = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in ("scripts/def_sameline.py", "scripts/def_day.sh")}
    variants = variant_samelines()
    rec("BUNDLE def_sameline.py is the post-bundle 22c6f0f0",
        hashlib.sha256(variants["BUNDLE"]).hexdigest().startswith("22c6f0f0"), hashlib.sha256(variants["BUNDLE"]).hexdigest()[:12])
    rec("FULL def_sameline.py carries 642", b"642: the HEK293 contingency" in variants["FULL"],
        hashlib.sha256(variants["FULL"]).hexdigest()[:12])
    with np.load(C4, allow_pickle=False) as z:
        Pc4, psc4 = np.asarray(z["P_down"], float), np.asarray(z["psign"], float)
    cov = np.isfinite(Pc4)
    with np.load(BOARD, allow_pickle=False) as z:
        bd = {k: np.asarray(z[k]) for k in z.files if k != "provenance"}
    for name, canon in variants.items():
        with tempfile.TemporaryDirectory(prefix=f"t604-{name}-") as td:
            base, pat = Path(td) / "base", Path(td) / "patched"
            for t in (base, pat):
                (t / "scripts").mkdir(parents=True)
                (t / "scripts/def_sameline.py").write_bytes(canon)
                os.symlink(ROOT / "data", t / "data")
                os.symlink(ROOT / "experiments", t / "experiments")
            c = subprocess.run([PY, str(HERE / "patch604.py"), "--root", str(pat), "--check"], capture_output=True, text=True)
            rec(f"[{name}] T1 --check", c.returncode == 0, (c.stdout + c.stderr).strip()[-160:])
            q = subprocess.run([PY, str(HERE / "patch604.py"), "--root", str(pat)], capture_output=True, text=True)
            rec(f"[{name}] T1 apply", q.returncode == 0, (q.stdout + q.stderr).strip()[-200:])
            r = subprocess.run([PY, str(HERE / "patch604.py"), "--root", str(pat)], capture_output=True, text=True)
            rec(f"[{name}] T1 second application refused", r.returncode != 0 and "already" in (r.stdout + r.stderr))
            u0 = run_sameline(base, Path(td) / "u0.npz")
            p0 = run_sameline(pat, Path(td) / "p0.npz")
            ok, d = same_arrays(u0, p0); rec(f"[{name}] T2 unset path unchanged", ok, d)
            u1 = run_sameline(base, Path(td) / "u1.npz", "--programme-source", "h1_2025_full")
            p1 = run_sameline(pat, Path(td) / "p1.npz", "--programme-source", "h1_2025_full")
            ok, d = same_arrays(u1, p1); rec(f"[{name}] T3 bank programme-source path unchanged", ok, d)
            i0 = run_sameline(pat, Path(td) / "i0.npz", "--programme-source", "identity:0")
            ok = all(np.array_equal(i0[f"programme_{c}"], psc4) and np.array_equal(i0[f"pdown_{c}"][cov], Pc4[cov]) for c in "DEF")
            rec(f"[{name}] T4 KA identity:0 == CLEAN4 (psign everywhere, P on covered genes)", ok)
            i5 = run_sameline(pat, Path(td) / "i5.npz", "--programme-source", "identity:0.5")
            ok5 = all(np.array_equal(i5[f"programme_{d}"], bd[f"programme_{a}"])
                      and np.array_equal(i5[f"pdown_{d}"], bd[f"pdown_{a}"], equal_nan=True) for d, a in zip("DEF", "ABC"))
            flips = {d: int((i5[f"programme_{d}"] != psc4).sum()) for d in "DEF"}
            rec(f"[{name}] T5 day-of on the relabelled bundle == the board artifact (programme + pdown)", ok5, f"flips vs CLEAN4 {flips}")
            # T9 (607): an empty / non-numeric dose is a clean SystemExit, and a negative dose is refused
            for bad, msg in (("identity:", "must be a number"), ("identity@", "must be a number"), ("identity:x", "must be a number"),
                             ("identity@-1", "must be >= 0")):
                cmd = [PY, str(pat / "scripts/def_sameline.py"), "--bundle", str(BUNDLE), "--identity", str(IDENT),
                       "--out", str(Path(td) / "bad.npz"), "--programme-source", bad]
                pr = subprocess.run(cmd, capture_output=True, text=True, cwd=str(pat))
                rec(f"[{name}] T9 refuses --programme-source {bad!r}", pr.returncode != 0 and msg in (pr.stdout + pr.stderr)
                    and "Traceback" not in (pr.stdout + pr.stderr))
            # T7/T8: the colon-free names def_day's hookbank form needs (DEF_DAY_MDE_SIGN=hookbank:identity:0)
            ib = run_sameline(pat, Path(td) / "ib.npz", "--programme-source", "identity")
            ok7 = all(np.array_equal(ib[f"programme_{d}"], i5[f"programme_{d}"])
                      and np.array_equal(ib[f"pdown_{d}"], i5[f"pdown_{d}"], equal_nan=True) for d in "DEF")
            rec(f"[{name}] T7 'identity' (bare) == identity:0.5 (the hookbank path's default)", ok7)
            with np.load(BOARD1, allow_pickle=False) as z:
                b1 = {k: np.asarray(z[k]) for k in z.files if k != "provenance"}
            ia = run_sameline(pat, Path(td) / "ia.npz", "--programme-source", "identity@1")
            ok8 = all(np.array_equal(ia[f"programme_{d}"], b1[f"programme_{a}"])
                      and np.array_equal(ia[f"pdown_{d}"], b1[f"pdown_{a}"], equal_nan=True) for d, a in zip("DEF", "ABC"))
            rec(f"[{name}] T8 'identity@1' == the lam 1 board artifact", ok8)
            # T6b: a context matched WITH rows keeps them under the identity source. Nadig Jurkat has no rows for this
            # essential-free panel, so force the HCT116 call (orion_hct116: ~290 rows). D's data is A's, and the pinned A
            # profile carries a DIFFERENT call (JURKAT), so only the 626 sister guard can drop it: programme_D must then
            # equal the board artifact's A (references B, C, h1).
            uf = run_sameline(base, Path(td) / "uf.npz", "--force", "D=HCT116_LARGE_INTESTINE")
            pf = run_sameline(pat, Path(td) / "pf.npz", "--force", "D=HCT116_LARGE_INTESTINE", "--programme-source", "identity:0.5")
            nrows = int(uf["rows_D"].shape[0]) if "rows_D" in uf else -1
            ok_rows = nrows > 0 and all(np.array_equal(pf[f"{k}_D"], uf[f"{k}_D"]) for k in ("rows", "signs", "targets", "w"))
            rec(f"[{name}] T6b forced same-line D (HCT116 bank) keeps its {nrows} rows / signs / targets / w", ok_rows)
            ok_sis = np.array_equal(pf["programme_D"], bd["programme_A"]) and np.array_equal(pf["pdown_D"], bd["pdown_A"], equal_nan=True)
            rec(f"[{name}] T6c the sister guard drops D's own data twin (pinned A, other call): programme_D == board A", ok_sis)
            if "rows_D" in u0:
                ok6 = all(np.array_equal(i5[f"{k}_D"], u0[f"{k}_D"]) for k in ("rows", "signs", "targets", "w"))
                rec(f"[{name}] T6 matched context D keeps its same-line rows / signs / targets / w", ok6, f"D rows {u0['rows_D'].shape[0]}")
            else:
                rec(f"[{name}] T6 matched context D keeps its rows", False, "D not matched in the no-source run")
    after = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in live}
    rec("live def_sameline.py and def_day.sh untouched", after == live)
    n_fail = sum(not ok for _, ok, _ in RESULTS)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} PASS")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
