"""813 KA for the UNAPPLIED lib801 patch (patches/lib801-telohaec.patch): existing folds stay bit-identical.

Run by whoever applies the patch (801 or the lab lead), on EC2, never by 813 before a GO.

  uv run python experiments/813-public-data/patches/ka_patch813.py snapshot <out.json>
      sha256 of every existing registry entry (canonical JSON), of every existing dump file it names (where present), of
      lib801.py, and FOLD_SOURCES / PANELS["abc"]["fold_ops"] for every pre-existing key.
  uv run python experiments/813-public-data/patches/ka_patch813.py check <before.json> <after.json>
      PASS iff every pre-existing key is byte-identical (new keys allowed only under "telohaec_gw").
  uv run python experiments/813-public-data/patches/ka_patch813.py rows <stored.npz> <redump.npz>
      a re-dump of an existing fold through the patched lib801 against its stored dump: max |d| per array; PASS iff 0
      (same instance type and thread count as the stored dump) -- otherwise the max |d| is reported against 801's 1e-12.
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
REG = ROOT / "experiments/801-lab-folds/work/registry.json"
NEW = "telohaec_gw"


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def snapshot(out):
    sys.path.insert(0, str(ROOT / "experiments/801-lab-folds"))
    import lib801 as L
    reg = json.loads(REG.read_text())
    rec = {"lib801_sha256": sha((ROOT / "experiments/801-lab-folds/lib801.py").read_bytes()), "entries": {}, "dumps": {},
           "fold_sources": {k: list(v) for k, v in L.FOLD_SOURCES.items()},
           "abc_fold_ops": {k: str(v) for k, v in L.PANELS["abc"]["fold_ops"].items()},
           "abc_folds": list(L.PANELS["abc"]["folds"])}
    for k, v in reg["folds"].items():
        rec["entries"][k] = sha(json.dumps(v, sort_keys=True).encode())
        p = ROOT / v["dump"]
        rec["dumps"][k] = sha(p.read_bytes()) if p.exists() else "absent-locally"
    Path(out).write_text(json.dumps(rec, indent=1))


def check(a, b):
    A, B = json.loads(Path(a).read_text()), json.loads(Path(b).read_text())
    bad = []
    for sec in ("entries", "dumps", "fold_sources", "abc_fold_ops"):
        for k, v in A[sec].items():
            if B[sec].get(k) != v:
                bad.append((sec, k))
        bad += [(sec, k) for k in B[sec] if k not in A[sec] and NEW not in k]
    if A["abc_folds"] != B["abc_folds"]:
        bad.append(("abc_folds", "changed"))
    print("PASS" if not bad else f"FAIL {bad}")
    sys.exit(1 if bad else 0)


def rows(a, b):
    A, B = np.load(a, allow_pickle=False), np.load(b, allow_pickle=False)
    worst = 0.0
    for k in A.files:
        x, y = A[k], B[k]
        if x.dtype.kind in "fc":
            d = float(np.max(np.abs(x.astype(np.float64) - y.astype(np.float64)))) if x.size else 0.0
        else:
            d = 0.0 if np.array_equal(x, y) else float("inf")
        worst = max(worst, d)
        print(f"{k}: max|d| {d:.3e}")
    print("PASS bit-identical" if worst == 0 else f"max|d| {worst:.3e} (801 tolerance 1e-12)")
    sys.exit(0 if worst <= 1e-12 else 1)


if __name__ == "__main__":
    {"snapshot": snapshot, "check": check, "rows": rows}[sys.argv[1]](*sys.argv[2:])
