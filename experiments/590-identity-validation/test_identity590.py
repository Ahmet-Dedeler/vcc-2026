"""590 — tests for patch_identity590.py, on COPIES (the live scripts are never written).

A temp root holds patched copies of scripts/def_identity.py and scripts/def_sameline.py; data/ and experiments/ are
symlinks to the real repo (read-only use), so ROOT-relative paths resolve exactly as in production.

T1  patched def_identity --verify PASSES, and its --out JSON equals the LIVE script's --out JSON on A/B/C except for the
    one added key; min_top1_score(A/B/C) equals 590's harness values (0.772 / 0.719 / 0.733 to 1e-3).
T2  patched def_sameline on synthetic identity files (bundle data/challenge, context A):
    a  registered + unanimous + min_top1 0.554 (the Xu HEK293 -> HAP1 case)      -> context unchanged, mode names 590
    b  same with min_top1 0.70                                                     -> same-line ON (rows as unpatched)
    c  no min_top1_score key (a pre-590 identity file)                            -> ON + WARNING (back-compatible)
    d  --min-top1 0 with 0.554                                                     -> ON (gate disabled)
    e  --force A=HAP1_ENGINEERED                                                   -> ON (known-answer path untouched)
    f  not unanimous                                                               -> unchanged by the old rule (as before)
    Every ON case's npz arrays and provenance equal the LIVE script's on the same input (no new field when passing).
T3  --check: anchors hold on the live files and on 530's post-patch def_sameline.py.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import patch_identity590 as P  # noqa: E402

PY = [sys.executable]
HAP1 = "HAP1_ENGINEERED"


def run(cmd: list[str], cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)


def make_root(tmp: Path, patched: bool) -> Path:
    r = tmp / ("patched" if patched else "live")
    (r / "scripts").mkdir(parents=True)
    for name in ("data", "experiments"):
        os.symlink(ROOT / name, r / name)
    ident = (ROOT / "scripts/def_identity.py").read_text()
    same = (ROOT / "scripts/def_sameline.py").read_text()
    if patched:
        ident, same = P.patch_identity(ident), P.patch_sameline(same)
    (r / "scripts/def_identity.py").write_text(ident)
    (r / "scripts/def_sameline.py").write_text(same)
    return r


def npz_equal(a: Path, b: Path) -> bool:
    za, zb = np.load(a, allow_pickle=False), np.load(b, allow_pickle=False)
    if sorted(za.files) != sorted(zb.files):
        return False
    return all(np.array_equal(za[k], zb[k], equal_nan=za[k].dtype.kind == "f") for k in za.files)


def main() -> None:
    fails = []
    with tempfile.TemporaryDirectory(dir=HERE / "work") as td:
        tmp = Path(td)
        live, pat = make_root(tmp, False), make_root(tmp, True)

        # T1
        r_live = run([*PY, "scripts/def_identity.py", "--verify", "--out", str(tmp / "id_live.json")], live)
        r_pat = run([*PY, "scripts/def_identity.py", "--verify", "--out", str(tmp / "id_pat.json")], pat)
        ok_verify = r_pat.returncode == 0 and "VERIFY PASS" in r_pat.stdout and r_live.returncode == 0
        jl, jp = json.loads((tmp / "id_live.json").read_text()), json.loads((tmp / "id_pat.json").read_text())
        extra = set(jp) - set(jl)
        # the two copies live in different temp roots: path-valued keys and timings differ by construction
        diff_keys = [k for k in jl if jl[k] != jp[k]]
        same_rest = set(diff_keys) <= {"timing_seconds", "bundle", "gene_axis_file"}
        if not same_rest:
            print("T1 differing keys:", diff_keys)
        s1 = jp.get("min_top1_score", {})
        want = {"A": 0.772, "B": 0.719, "C": 0.733}
        ok_vals = all(abs(s1.get(c, -1) - v) < 1e-3 for c, v in want.items())
        t1 = ok_verify and extra == {"min_top1_score"} and same_rest and ok_vals
        print(f"T1 verify={ok_verify} added={sorted(extra)} rest_identical={same_rest} min_top1={s1} -> {'PASS' if t1 else 'FAIL'}")
        if not t1:
            fails.append("T1")
            print(r_pat.stdout[-2000:], r_pat.stderr[-2000:])

        # T2
        def ident_file(name: str, call: str, unanimous: bool, s1: float | None) -> Path:
            d = {"call": {"A": call}, "call_unanimous_across_methods": {"A": unanimous}}
            if s1 is not None:
                d["min_top1_score"] = {"A": s1}
            p = tmp / f"{name}.json"
            p.write_text(json.dumps(d))
            return p

        def sameline(root: Path, ident: Path | None, out: str, extra_args: list[str] = ()) -> tuple[subprocess.CompletedProcess, dict]:
            args = [*PY, "scripts/def_sameline.py", "--bundle", str(ROOT / "data/challenge"), "--out", str(tmp / out)]
            args += ["--identity", str(ident)] if ident else []
            r = run(args + list(extra_args), root)
            prov = json.loads(str(np.load(tmp / out)["provenance"])) if (tmp / out).exists() else {}
            return r, prov

        cases = [
            ("a", ident_file("a", HAP1, True, 0.554), [], False, "590"),
            ("b", ident_file("b", HAP1, True, 0.70), [], True, None),
            ("c", ident_file("c", HAP1, True, None), [], True, "WARNING"),
            ("d", ident_file("d", HAP1, True, 0.554), ["--min-top1", "0"], True, None),
            ("e", None, ["--force", f"A={HAP1}"], True, None),
            ("f", ident_file("f", HAP1, False, 0.70), [], False, "not unanimous"),
        ]
        for name, idf, xargs, want_on, want_text in cases:
            rp, prov = sameline(pat, idf, f"s_{name}_pat.npz", xargs)
            xargs_live = [x for x in xargs if x not in ("--min-top1", "0")]
            rl, provl = sameline(live, idf, f"s_{name}_live.npz", xargs_live)
            on = "A" in list(np.load(tmp / f"s_{name}_pat.npz")["contexts"]) if rp.returncode == 0 else None
            text_ok = want_text is None or want_text in (rp.stdout + json.dumps(prov))
            if rp.returncode != 0 or not (tmp / f"s_{name}_pat.npz").exists():
                same_as_live = False
            elif want_on:  # passing / bypass cases must equal the live script's output exactly
                same_as_live = npz_equal(tmp / f"s_{name}_pat.npz", tmp / f"s_{name}_live.npz")
            else:  # blocked case: identical except A's record; f must equal live entirely
                same_as_live = (name != "f") or npz_equal(tmp / f"s_{name}_pat.npz", tmp / f"s_{name}_live.npz")
            ok = rp.returncode == 0 and on == want_on and text_ok and same_as_live
            print(f"T2{name} rc={rp.returncode} on={on} (want {want_on}) text={text_ok} same_as_live={same_as_live} "
                  f"mode={prov.get('contexts', {}).get('A', {}).get('mode')!r} -> {'PASS' if ok else 'FAIL'}")
            if not ok:
                fails.append(f"T2{name}")
                print(rp.stdout[-1500:], rp.stderr[-1500:])

    # T3
    r = run([*PY, str(HERE / "patch_identity590.py"), "--check"], ROOT)
    # pre-530 live file: --check also proves the anchors on the post-530 file; post-530 live file (538 applied): the
    # first line already is the post-530 proof (594 found the old assertion required the pre-530 branch only)
    post530_live = "--programme-estimator" in (ROOT / "scripts/def_sameline.py").read_text()
    t3 = r.returncode == 0 and "nothing written" in r.stdout and (
        post530_live or "after 530's patch_sameline: 590's anchors OK" in r.stdout)
    print(f"T3 --check rc={r.returncode} -> {'PASS' if t3 else 'FAIL'}\n{r.stdout.strip()}")
    if not t3:
        fails.append("T3")
        print(r.stderr[-1500:])
    print("ALL PASS" if not fails else f"FAILED: {fails}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
