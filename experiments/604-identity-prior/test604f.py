#!/usr/bin/env python3
"""604f tests: the day-of octile-flip pass-through, applied after the whole 10-06 chain (606's order). ~7 GB peak (the
relabelled bundle's control CPM, like test604); run it through 602's gate.
  T1  --check, apply, second apply refused; refused on a tree without 644 (604f anchors on the end state of the chain)
  T2  bash -n def_day.sh; def_sameline.py compiles
  T3  unset / 00000000: def_sameline's identity:0.5 artifact is array-identical to the unpatched chain's
  T4  day-of == board: identity:0 + 10000000 and identity:0.5 + 10010001 on the relabelled bundle (D/E/F = A/B/C data)
      equal build604 --flip's A/B/C artifacts (programme + pdown)
  T5  def_sameline refuses --flip-octiles with a bank programme source and with none
  T6  def_day's validation block (the inserted text, in a bash harness): unset no-op; identity source -> block flip;
      hookbank:identity -> hook flip; a clean4 hook is refused; nothing to carry the flip is refused; bad masks refused;
      both def_sameline calls carry the flag
usage: uv run python experiments/604-identity-prior/test604f.py
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
sys.path.insert(0, str(ROOT / "experiments/598-tercilec-knob"))
sys.path.insert(0, str(ROOT / "experiments/606-applyday-composition"))
sys.path.insert(0, str(HERE))
import test598 as T598  # noqa: E402
import compose606 as C606  # noqa: E402
import patch604f as P  # noqa: E402
from test604 import BUNDLE, IDENT, same_arrays  # noqa: E402

RES: list[tuple[bool, str]] = []


def rec(ok, name):
    RES.append((bool(ok), name))
    print(f"{'PASS' if ok else 'FAIL'} {name}", flush=True)


def chain_tree(t: Path, order) -> Path:
    T598.canon_tree(t)
    (t / "data/identity_refs").mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / "data/identity_refs/hek293_nonflex_ref642.npz", t / "data/identity_refs/hek293_nonflex_ref642.npz")
    (t / "src/vcc2026").mkdir(parents=True, exist_ok=True)
    for f in (ROOT / "src/vcc2026").glob("*.py"):
        shutil.copy2(f, t / "src/vcc2026" / f.name)
    for name in order:
        c = subprocess.run([PY, str(ROOT / C606.PATCHES[name][0]), "--root", str(t), "--check"], capture_output=True, text=True)
        co = (c.stdout + c.stderr).lower()
        if c.returncode != 0 and "already" in co and "partially" not in co:
            continue   # B2: already on HEAD (after the 10-05 / 10-06 applies)
        r = subprocess.run([PY, str(ROOT / C606.PATCHES[name][0]), "--root", str(t)], capture_output=True, text=True)
        if r.returncode != 0:
            raise SystemExit(f"setup {name}: {(r.stdout + r.stderr)[-400:]}")
    return t


def runnable(t: Path) -> None:
    """def_sameline / def_identity_prior read data/ and experiments/ relative to the tree root. The chain tree already
    holds a few real files under both (e.g. experiments/524-.../gate524.py, data/identity_refs/), so link the rest entry by
    entry instead of linking the whole directory."""
    for top in ("experiments", "data"):
        (t / top).mkdir(exist_ok=True)
        for sub in (ROOT / top).iterdir():
            if not (t / top / sub.name).exists():
                os.symlink(sub, t / top / sub.name)


def sameline(t: Path, out: Path, *extra) -> tuple[int, dict | None, str]:
    cmd = [PY, str(t / "scripts/def_sameline.py"), "--bundle", str(BUNDLE), "--identity", str(IDENT), "--out", str(out), *extra]
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=str(t))
    if p.returncode != 0:
        return p.returncode, None, (p.stdout + p.stderr)[-400:]
    with np.load(out, allow_pickle=False) as z:
        return 0, {k: np.asarray(z[k]) for k in z.files if k != "provenance"}, ""


def main() -> int:
    live = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in ("scripts/def_sameline.py", "scripts/def_day.sh")}
    with tempfile.TemporaryDirectory(prefix="t604f-") as tds:
        td = Path(tds)
        base = chain_tree(td / "base", C606.DEFAULT_ORDER)
        pat = chain_tree(td / "pat", C606.DEFAULT_ORDER)
        c = subprocess.run([PY, str(HERE / "patch604f.py"), "--root", str(pat), "--check"], capture_output=True, text=True)
        rec(c.returncode == 0, f"T1 --check: {c.stdout.strip()[:120]}")
        q = subprocess.run([PY, str(HERE / "patch604f.py"), "--root", str(pat)], capture_output=True, text=True)
        rec(q.returncode == 0, f"T1 apply: {q.stdout.strip()[:120]}")
        r = subprocess.run([PY, str(HERE / "patch604f.py"), "--root", str(pat)], capture_output=True, text=True)
        rec(r.returncode != 0 and "already" in (r.stdout + r.stderr), "T1 second apply refused")
        no644 = chain_tree(td / "no644", [x for x in C606.DEFAULT_ORDER if x != "644"])
        r = subprocess.run([PY, str(HERE / "patch604f.py"), "--root", str(no644), "--check"], capture_output=True, text=True)
        rec(r.returncode != 0 and "after the whole 10-06 chain" in (r.stdout + r.stderr), "T1 refused on a chain without 644")
        b = subprocess.run(["bash", "-n", str(pat / "scripts/def_day.sh")], capture_output=True, text=True)
        rec(b.returncode == 0, f"T2 bash -n def_day.sh {b.stderr.strip()[:120]}")
        try:
            compile((pat / "scripts/def_sameline.py").read_text(), "def_sameline.py", "exec")
            rec(True, "T2 def_sameline.py compiles")
        except SyntaxError as e:
            rec(False, f"T2 def_sameline.py compiles: {e}")
        for t in (base, pat):
            runnable(t)
        rc0, u, e0 = sameline(base, td / "u.npz", "--programme-source", "identity:0.5")
        rc1, p0, e1 = sameline(pat, td / "p0.npz", "--programme-source", "identity:0.5")
        rc2, pz, e2 = sameline(pat, td / "pz.npz", "--programme-source", "identity:0.5", "--flip-octiles", "00000000")
        ok_a, d_a = same_arrays(u, p0) if (rc0 == rc1 == 0) else (False, e0 + e1)
        ok_b, d_b = same_arrays(u, pz) if (rc0 == rc2 == 0) else (False, e2)
        rec(ok_a, f"T3 unset: identity:0.5 artifact unchanged by 604f ({d_a})")
        rec(ok_b, f"T3 00000000: identity:0.5 artifact unchanged ({d_b})")
        for lam, mask in (("0", "10000000"), ("0.5", "10010001")):
            o = td / f"board_{lam}"
            o.mkdir()
            bb = subprocess.run([PY, str(HERE / "build604.py"), "--lam", lam, "--flip", mask, "--out-dir", str(o)],
                                capture_output=True, text=True)
            with np.load(o / f"ident604v3_flex_lam{lam}_flip{mask}_ABC.npz", allow_pickle=False) as z:
                bd = {k: np.asarray(z[k]) for k in z.files if k != "provenance"}
            rc, f, ef = sameline(pat, td / f"f{lam}.npz", "--programme-source", f"identity:{lam}", "--flip-octiles", mask)
            ok = bb.returncode == 0 and rc == 0 and all(
                np.array_equal(f[f"programme_{d}"], bd[f"programme_{a}"])
                and np.array_equal(f[f"pdown_{d}"], bd[f"pdown_{a}"], equal_nan=True) for d, a in zip("DEF", "ABC"))
            rec(ok, f"T4 day-of identity:{lam} + {mask} == build604 --flip board artifact (D/E/F vs A/B/C) {ef[:120]}")
        rc, _, e = sameline(pat, td / "x1.npz", "--programme-source", "h1_2025_full", "--flip-octiles", "10000000")
        rec(rc != 0 and "needs an identity programme source" in e, "T5 refuses a flip with a bank programme source")
        rc, _, e = sameline(pat, td / "x2.npz", "--flip-octiles", "10000000")
        rec(rc != 0 and "needs an identity programme source" in e, "T5 refuses a flip with no programme source")
        # T6: the def_day validation block, cut from the patched file, in a harness
        day = (pat / "scripts/def_day.sh").read_text()
        vblock = P.DAY_V_NEW[: len(P.DAY_V_NEW) - len(P.DAY_V_OLD)]
        rec(day.count(vblock) == 1 and day.count('${FLIP_MAIN:+--flip-octiles "$FLIP_MAIN"}') == 1
            and day.count('${FLIP_HOOK:+--flip-octiles "$FLIP_HOOK"}') == 1, "T6 both def_sameline calls carry the flag")
        w = td / "h"
        (w / "scripts").mkdir(parents=True)
        (w / "scripts/def_sameline.py").write_text("# --flip-octiles\n")

        def harness(env):
            script = ('set -uo pipefail\ndie(){ echo "DIE: $*"; exit 7; }\nlog(){ echo "LOG: $*"; }\n' + vblock
                      + '\necho "MAIN=${FLIP_MAIN:-} HOOK=${FLIP_HOOK:-}"\n')
            return subprocess.run(["bash", "-c", script], cwd=w, capture_output=True, text=True,
                                  env={"PATH": "/usr/bin:/bin", **env})
        h = harness({})
        rec(h.returncode == 0 and "MAIN= HOOK=" in h.stdout and "LOG" not in h.stdout, "T6 unset: no-op")
        h = harness({"DEF_DAY_DE_PRIOR_FLIP": "10000000", "DEF_DAY_PROGRAMME_SOURCE": "identity:0.5"})
        rec(h.returncode == 0 and "MAIN=10000000 HOOK=" in h.stdout, "T6 identity source -> block-prior flip")
        h = harness({"DEF_DAY_DE_PRIOR_FLIP": "10000000", "DEF_DAY_MDE_SIGN": "hookbank:identity:0"})
        rec(h.returncode == 0 and "MAIN= HOOK=10000000" in h.stdout, "T6 hookbank:identity -> hook flip")
        h = harness({"DEF_DAY_DE_PRIOR_FLIP": "10000000", "DEF_DAY_PROGRAMME_SOURCE": "identity", "DEF_DAY_MDE_SIGN": "hookbank:identity@1:0"})
        rec(h.returncode == 0 and "MAIN=10000000 HOOK=10000000" in h.stdout, "T6 both carriers -> both flipped")
        h = harness({"DEF_DAY_DE_PRIOR_FLIP": "10000000", "DEF_DAY_PROGRAMME_SOURCE": "identity", "DEF_DAY_MDE_SIGN": "clean4:0.05"})
        rec(h.returncode == 7 and "reads CLEAN4 directly" in h.stdout, "T6 a clean4 hook is refused with a flip")
        h = harness({"DEF_DAY_DE_PRIOR_FLIP": "10000000", "DEF_DAY_PROGRAMME_SOURCE": "h1_2025_full"})
        rec(h.returncode == 7 and "nothing else can carry the flip" in h.stdout, "T6 a flip with nothing to carry it is refused")
        for bad in ("1000000", "100000000", "1000000x", "10000002"):
            h = harness({"DEF_DAY_DE_PRIOR_FLIP": bad, "DEF_DAY_PROGRAMME_SOURCE": "identity"})
            rec(h.returncode == 7 and "8 characters of 0/1" in h.stdout, f"T6 refuses mask {bad!r}")
    after = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in live}
    rec(after == live, "live def_sameline.py and def_day.sh untouched")
    n = sum(ok for ok, _ in RES)
    print(f"\n{n}/{len(RES)} PASS" if n == len(RES) else f"\n{len(RES) - n} FAILED of {len(RES)}")
    return 0 if n == len(RES) else 1


if __name__ == "__main__":
    sys.exit(main())
