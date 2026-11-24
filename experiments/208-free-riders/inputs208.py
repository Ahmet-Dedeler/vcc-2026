"""208 — rule 11: a PATH and a SHA256 for every operator and every frozen input.

"A preregistration that says 'essential-pruned held-HCT union operator' is not resolvable to a
file."  Nine near-identically-named whitening operators exist across `experiments/` and
`data/donors/`, two of them agreeing to five significant figures on their largest eigenvalue.
Only the sha separates them.  The shipped one is `data/donors/whitener_union_pruned.npz`,
17,942 x 17,942, sha256 prefix `e2f65ecf`, and both of today's board entries record exactly
that in `meta['operator']['sha256']`.

    uv run python inputs208.py
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)

FROZEN = [
    "submissions/191-rc_hybrid.h5ad",
    "submissions/202-rc_callset.h5ad",
    "submissions/202-rc_callset.h5ad.meta.json",
    "submissions/191-rc_hybrid.h5ad.meta.json",
    "submissions/a15-cd4f010.vcc.meta.json",
    "src/vcc2026/rc_emit.py",
    "src/vcc2026/donor_transfer.py",
    "scripts/make_submission_pm_cd4.py",
    "data/harness/h1-test.h5ad",
    "experiments/178-composed-riders/work/eff.a15.npy",
    "experiments/169-sharpening-six-member/work/h1test-effect.npz",
    "experiments/185-rank-controlled-emission/work/a15_totals.npz",
    "experiments/191-rc-hybrid-build/work/de191.a15.npz",
    "experiments/191-rc-hybrid-build/work/de191.rc_hybrid_seq.npz",
    "experiments/202-callset-rider/work/de202.rc_callset.npz",
    "experiments/202-callset-rider/work/de202.rc_hybrid.npz",
    "experiments/196-emitted-depth/work/chal-totals-a15.npz",
    "experiments/196-emitted-depth/ham196.py",
    "experiments/196-emitted-depth/recon196.py",
    "experiments/188-magnitude-calibration/work/chalframes.npz",
    "experiments/096-mse-floor/work/board_solve.npz",
    "data/challenge/context_A.h5ad",
    "data/challenge/context_B.h5ad",
    "data/challenge/context_C.h5ad",
]


def sha(path: Path, chunk: int = 1 << 24) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def main() -> None:
    out = {"operators": [], "frozen": []}
    print("=== OPERATORS (rule 11: nine near-identically-named files exist) ===")
    print(f"{'path':<62}{'sha256[16]':<18}{'max eig':>10}{'n':>8}")
    for p in sorted(ROOT.glob("**/whitener*.npz")):
        if ".venv" in str(p):
            continue
        rel = str(p.relative_to(ROOT))
        s = sha(p)
        try:
            z = np.load(p, allow_pickle=True)
            k = "eigenvalues" if "eigenvalues" in z else (
                "evals" if "evals" in z else list(z.keys())[0])
            v = np.asarray(z[k]).ravel()
            mx, nn = float(np.max(v)), int(v.size)
        except Exception as exc:                                   # noqa: BLE001
            mx, nn = float("nan"), -1
            print(f"  (could not read eigenvalues from {rel}: {exc})")
        out["operators"].append(dict(path=rel, sha256=s, max_eig=mx, n=nn))
        star = "  <-- SHIPPED" if s.startswith("e2f65ecf") else ""
        print(f"{rel:<62}{s[:16]:<18}{mx:>10.6f}{nn:>8}{star}")

    print("\n=== FROZEN INPUTS ===")
    print(f"{'path':<62}{'sha256[16]':<18}{'bytes':>17}")
    for rel in FROZEN:
        p = ROOT / rel
        if not p.exists():
            print(f"{rel:<62}{'MISSING':<18}")
            out["frozen"].append(dict(path=rel, sha256=None))
            continue
        s = sha(p)
        out["frozen"].append(dict(path=rel, sha256=s, bytes=p.stat().st_size))
        print(f"{rel:<62}{s[:16]:<18}{p.stat().st_size:>17,}")

    print("\n=== the two board entries' recorded operator, straight out of their meta ===")
    for rel in ("submissions/191-rc_hybrid.h5ad.meta.json",
                "submissions/202-rc_callset.h5ad.meta.json"):
        p = ROOT / rel
        if not p.exists():
            continue
        m = json.loads(p.read_text())
        op = m.get("operator", {})
        print(f"  {rel:<52} operator.sha256 {str(op.get('sha256'))[:16]}  "
              f"columns {str(op.get('columns_sha256'))[:16]}  "
              f"effects.final {str(m.get('effects', {}).get('final_sha256'))[:16]}")
        out.setdefault("meta", {})[rel] = dict(
            operator=op, effects_final=m.get("effects", {}).get("final_sha256"),
            whitening=m.get("effects", {}).get("whitening"))

    print("\n=== the two package sha256, against the values the board results record ===")
    want = {"submissions/191-rc_hybrid.h5ad":
            "e0b0677f984a94a1c4a7af2add196908cdab90d7acb463eb1d6e0397c2d870a8",
            "submissions/202-rc_callset.h5ad":
            "2a894aa543db7f574ebc024dbb3ebd5e36f923827c4e9a01f263092fe49692e0"}
    ok = True
    for rel, w in want.items():
        got = next((e["sha256"] for e in out["frozen"] if e["path"] == rel), None)
        hit = got == w
        ok &= hit
        print(f"  {rel:<40} {str(got)[:16]}  expect {w[:16]}  "
              f"{'PASS' if hit else 'FAIL'}")
    out["package_sha_ok"] = ok
    json.dump(out, open(WORK / "inputs208.json", "w"), indent=1, default=float)
    if not ok:
        raise SystemExit("ABORT (rule 11): a board package on disk is not the one scored")


if __name__ == "__main__":
    main()
