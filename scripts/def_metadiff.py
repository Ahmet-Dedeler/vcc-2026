#!/usr/bin/env python3
"""Content-aware config diff between two def_day final metas (694 v2 rule, 2026-10-16).

Every def_day run has its own work dir, and two of its artifacts carry that dir inside a string field: the stack dump's
`env` (VCC_BULK_DUMP) and the same-line artifact's `provenance` (identity_source). Their FILE shas (`bulk_load_sha256`,
`sameline.sha256`) therefore differ between two runs even when every array is the same. A plain meta diff then reports a
change that is not one: on 10-16 it failed R-ABC against 524dr2. Against R0 it would wrongly block a replacement final entry,
and silently keep a worse R0.

Rule: `sameline.path` is ignored; float config values compare to --float-rtol (default 1e-12 relative) with an absolute floor
--float-atol (default 1e-12), so round-off noise in near-zero diagnostics passes; `--float-rtol 0 --float-atol 0` = exact. `bulk_load_sha256` / `sameline` count as unchanged ("provenance-only") only if their npz files
are equal in content: the same keys, the same arrays (NaN-aware, dtype and shape included), and string fields equal once
each run's work dir is replaced by one token. Each npz is found in its run dir and identified by the sha its meta records,
so a dump hook's replacement dump (<run>/dumphook/<name>) is told apart from the stack dump of the same name. Every
other differing key is a real difference.

    uv run python scripts/def_metadiff.py <new.meta.json> <new_run_dir> <ref.meta.json> <ref_run_dir> \
        [--expect k1,k2 | --allow k1,k2 [--need k1,k2]] [--json out.json]
--expect: the real diff must equal this set exactly. --allow / --need: it must be a subset of allow and contain need.
Prints "[metadiff] PASS|FAIL: diff [...]; provenance-only [...]" and exits 0 / 1 (2 = bad input).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def _sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def _tokens(run: Path) -> list[str]:
    r = run.resolve()
    out = {str(r), str(run)}
    try:
        out.add(str(r.relative_to(ROOT)))
    except ValueError:
        pass
    return sorted(out, key=len, reverse=True)   # longest first, so a relative path never eats part of an absolute one


def _norm(s: str, run: Path) -> str:
    for t in _tokens(run):
        s = s.replace(t, "<WORK>")
    return s


def npz_same(fa: Path, run_a: Path, fb: Path, run_b: Path) -> tuple[bool, str]:
    """content equality of two npz files up to each run's work dir inside scalar string fields"""
    a, b = np.load(fa, allow_pickle=False), np.load(fb, allow_pickle=False)
    if sorted(a.files) != sorted(b.files):
        return False, f"keys {sorted(a.files)} vs {sorted(b.files)}"
    texts = []
    for k in sorted(a.files):
        x, y = a[k], b[k]
        if x.dtype.kind == "U" and y.dtype.kind == "U" and x.shape == y.shape == ():
            if _norm(str(x), run_a) != _norm(str(y), run_b):
                return False, f"{k} differs beyond the work dir"
            if str(x) != str(y):
                texts.append(k)
        elif x.shape != y.shape or x.dtype != y.dtype or not np.array_equal(x, y, equal_nan=x.dtype.kind in "fc"):
            return False, f"array {k} differs"
    return True, f"{len(a.files)} arrays equal (NaN-aware); {texts or 'no'} string field(s) differ only by the work dir"


def _find(run: Path, name: str, want_sha: str) -> Path | None:
    for cand in (run / name, run / "dumphook" / name):
        if cand.is_file() and _sha(cand) == want_sha:
            return cand
    return None


def _eq(a, b, rtol: float, atol: float = 1e-12) -> bool:
    """config equality with floats compared to a relative tolerance plus an absolute floor (everything else exact, including
    types of non-floats). 10-17: a dump hook rewrites the dump and the tercilec `edges` moved by 1 ULP (…832 vs …835) with every
    count, flag and assignment identical: float round-off, not a difference. 10-18 (698): tercilec's `key_ratio_maxdev`, itself
    a round-off diagnostic, read 8.7e-15 vs 1.0e-14 (15 % relative): no relative tolerance passes noise near zero, hence atol."""
    if isinstance(a, float) and isinstance(b, float):
        return a == b or abs(a - b) <= max(rtol * max(abs(a), abs(b)), atol)
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_eq(a[k], b[k], rtol, atol) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_eq(x, y, rtol, atol) for x, y in zip(a, b))
    return a == b


def metadiff(new_meta: Path, new_run: Path, ref_meta: Path, ref_run: Path, rtol: float = 1e-12, atol: float = 1e-12) -> dict:
    c = dict(json.loads(new_meta.read_text())["config"])
    r = dict(json.loads(ref_meta.read_text())["config"])
    for d in (c, r):
        if isinstance(d.get("sameline"), dict):
            d["sameline"] = {k: v for k, v in d["sameline"].items() if k != "path"}
    diff = sorted(k for k in set(c) | set(r) if not _eq(c.get(k), r.get(k), rtol, atol))
    prov, notes = [], []
    if "bulk_load_sha256" in diff and c.get("bulk_load") == r.get("bulk_load") and c.get("bulk_load"):
        fa = _find(new_run, str(c["bulk_load"]), str(c["bulk_load_sha256"]))
        fb = _find(ref_run, str(r["bulk_load"]), str(r["bulk_load_sha256"]))
        if fa and fb:
            same, why = npz_same(fa, new_run, fb, ref_run)
        else:
            same, why = False, f"dump not found by its recorded sha (new {fa}, ref {fb})"
        notes.append(f"bulk_load_sha256: {'provenance only' if same else 'REAL'}: {why}")
        if same:
            diff.remove("bulk_load_sha256"); prov.append("bulk_load_sha256")
    if "sameline" in diff:
        sa = {k: v for k, v in (c.get("sameline") or {}).items() if k != "sha256"}
        sb = {k: v for k, v in (r.get("sameline") or {}).items() if k != "sha256"}
        if sa != sb:
            same, why = False, "the same-line record differs beyond its sha256"
        else:
            fa = _find(new_run, "sameline.npz", str((c.get("sameline") or {}).get("sha256")))
            fb = _find(ref_run, "sameline.npz", str((r.get("sameline") or {}).get("sha256")))
            same, why = (npz_same(fa, new_run, fb, ref_run) if fa and fb
                         else (False, f"sameline.npz not found by its recorded sha (new {fa}, ref {fb})"))
        notes.append(f"sameline: {'provenance only' if same else 'REAL'}: {why}")
        if same:
            diff.remove("sameline"); prov.append("sameline.sha256")
    return {"diff": diff, "provenance_only": prov, "notes": notes}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("new_meta", type=Path); ap.add_argument("new_run", type=Path)
    ap.add_argument("ref_meta", type=Path); ap.add_argument("ref_run", type=Path)
    ap.add_argument("--expect"); ap.add_argument("--allow"); ap.add_argument("--need"); ap.add_argument("--json", type=Path)
    ap.add_argument("--float-rtol", type=float, default=1e-12, help="relative tolerance for float config values")
    ap.add_argument("--float-atol", type=float, default=1e-12, help="absolute floor for float config values (with rtol 0: exact)")
    a = ap.parse_args()
    for p in (a.new_meta, a.ref_meta):
        if not p.is_file():
            print(f"[metadiff] bad input: {p} missing", file=sys.stderr); return 2
    for p in (a.new_run, a.ref_run):
        if not p.is_dir():
            print(f"[metadiff] bad input: {p} is not a directory", file=sys.stderr); return 2
    res = metadiff(a.new_meta, a.new_run, a.ref_meta, a.ref_run, rtol=a.float_rtol, atol=a.float_atol)
    split = lambda s: sorted(x for x in (s or "").split(",") if x)
    d = set(res["diff"])
    if a.expect is not None:
        ok = d == set(split(a.expect)); rule = f"== {split(a.expect)}"
    elif a.allow is not None:
        ok = d <= set(split(a.allow)) and set(split(a.need)) <= d; rule = f"within {split(a.allow)}, containing {split(a.need)}"
    else:
        ok = not d; rule = "empty"
    res.update(rule=rule, ok=ok)
    for n in res["notes"]:
        print(f"  {n}")
    print(f"[metadiff] {'PASS' if ok else 'FAIL'}: diff {res['diff']} (want {rule}); provenance-only {res['provenance_only']}")
    if a.json:
        a.json.write_text(json.dumps(res, indent=1))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
