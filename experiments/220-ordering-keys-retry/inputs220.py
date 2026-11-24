"""220 — rule 11: a path and a sha256 for every frozen input, never a description.

Eight whitening operators and several near-identically-named DE tables exist across
`experiments/`.  This experiment opens NO donor table and NO operator -- it works entirely
from the frozen a15 h1-test artifact and 169's recovered effect -- so the trap that cost a
whole evening on 2026-08-31 cannot fire here, but the inventory is recorded anyway.

    uv run python inputs220.py
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
WORK = Path(__file__).resolve().parent / "work"

FILES = [
    ("DE surface, OLD instrument (189/190's)",
     "experiments/164-callset-artefact/de164.a15.npz"),
    ("DE surface, NEW instrument (202's own emitted table)",
     "experiments/202-callset-rider/work/de202.rc_callset.npz"),
    ("the recovered design-space effect H (169)",
     "experiments/169-sharpening-six-member/work/h1test-effect.npz"),
    ("189's delivered per-gene prior",
     "experiments/189-per-gene-reliability/work/prior189.npz"),
    ("190's CollecTRI source", "data/external/regulon/collectri_rescued.csv"),
    ("the cached h1-test control block (185)",
     "experiments/185-rank-controlled-emission/work/ctrl_meta.npz"),
    ("...its sorted control matrix",
     "experiments/185-rank-controlled-emission/work/ctrl_cs.npy"),
    ("a15's own delivered h1-test per-gene totals",
     "experiments/185-rank-controlled-emission/work/a15_totals.npz"),
    ("the reference side of every table (191)",
     "experiments/191-rc-hybrid-build/work/de191.a15.npz"),
    ("214's three-seed arm A tables — seed 1",
     "experiments/214-cheap-positives/work/de214.reb_pm00_floor.npz"),
    ("214's three-seed incumbent tables — seed 1",
     "experiments/214-cheap-positives/work/de214.reb_pm50_floor.npz"),
    ("the SHIPPED emitter", "src/vcc2026/rc_emit.py"),
    ("the SHIPPED quasi count law", "src/vcc2026/donor_transfer.py"),
    ("172's instrument", "experiments/172-sign-confidence-key/harness172.py"),
    ("153's scorer", "experiments/153-de-bench/debench.py"),
]


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


def main() -> None:
    out = []
    print(f"{'what':<52}{'sha256[16]':<18}{'MB':>8}  path")
    for what, rel in FILES:
        p = ROOT / rel
        if not p.exists():
            print(f"{what:<52}{'MISSING':<18}{'':>8}  {rel}")
            out.append(dict(what=what, path=rel, sha256=None))
            continue
        d = sha(p)
        print(f"{what:<52}{d[:16]:<18}{p.stat().st_size/1e6:>8.1f}  {rel}")
        out.append(dict(what=what, path=rel, sha256=d, bytes=p.stat().st_size))

    print("\nDONOR AUDIT — this experiment opens NO donor table.")
    print("  `data/donors/h1_2025` (85/85 coverage, pure oracle leakage) is never touched.")
    print("  `k562_gwps` and `orion_hek293t` are not opened either: the effect matrix is")
    print("  169's recovered `H` for the FROZEN a15 artifact, so no donor pool is re-scaled")
    print("  and `Surface(n_real_target=N)` is never called anywhere in this directory.")

    z = np.load(ROOT / "experiments/169-sharpening-six-member/work/h1test-effect.npz",
                allow_pickle=True)
    print(f"\n  H: {z['H'].shape} on {len(z['targets'])} targets x {len(z['genes'])} genes;"
          f" min_cpm {float(z['min_cpm'])}")
    json.dump(out, open(WORK / "inputs220.json", "w"), indent=1)
    print(f"\nwrote {WORK/'inputs220.json'}")


if __name__ == "__main__":
    main()
