#!/usr/bin/env python3
"""693 KA on REAL A/B/C blocks, with the real build as the judge of the input.

Input = 630arm's package (446bj2z + UPRN): DNMAG's exact input in the DN build (the order is ballast → zero repair → zs →
UP → DN, then the writer). Judge = the DN package (630arm + VCC_DE_DNMAG=d_m200, built 10-15 with v1).
For a few context-A targets:
  (1) v1.apply(input block) must equal the DN package's block: proves the input, n_car and controls are reconstructed
      faithfully, so the comparison below is on the real day's object;
  (2) v2.apply(input block) must equal v1.apply(input block), blocks and diagnostics.
n_car = 400 - ballast_K from the DN build's rc_audit dump; the scorer controls = rank_mag.scorer_control_block(context A's
controls, min_cpm 5), as the builder makes them; own_col = the target's kept-axis column.
Restores both packages to work/ (APFS-cheap, deleted by --clean). Peak ~4 GB.
usage: uv run python experiments/693-dnmag-speed/abc693.py [n_targets=3] [--clean]
"""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
from vcc2026 import dn_mag as V1  # noqa: E402
from vcc2026 import rank_mag as RM  # noqa: E402

spec = importlib.util.spec_from_file_location("dn_mag2", HERE / "dn_mag2.py")
V2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(V2)

STEM = "submissions/446bj2z-uprn-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-uprnk100"
PKG_IN, PKG_DN = ROOT / f"{STEM}.vcc", ROOT / f"{STEM}-dnm200.vcc"
AUDIT = ROOT / "experiments/671-dnmag-arm/work/dn631.rc_audit.jsonl"
W = HERE / "work"
n_t = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 3
CTX = "A"


def restore(vcc: Path, out: Path) -> Path:
    if not out.exists():
        subprocess.run(["uv", "run", "python", str(ROOT / "scripts/restore_prediction_h5ad.py"), str(vcc), str(out)],
                       check=True, cwd=ROOT)
    return out


# the restore tool's reader: every string encoding our prediction files use (plain, categorical, nullable values)
_rs = importlib.util.spec_from_file_location("restore_pred", ROOT / "scripts/restore_prediction_h5ad.py")
RP = importlib.util.module_from_spec(_rs)
_rs.loader.exec_module(RP)


def labels(h, key):
    return np.array(RP.strings(h["obs"], key), dtype=object)


def rows(h, idx, n_genes):
    X = h["X"]
    ptr = X["indptr"]
    out = np.zeros((len(idx), n_genes), dtype=np.int64)
    for k, i in enumerate(idx):
        a, b = int(ptr[i]), int(ptr[i + 1])
        out[k, X["indices"][a:b]] = np.rint(X["data"][a:b]).astype(np.int64)
    return out


def main() -> int:
    hin = h5py.File(restore(PKG_IN, W / "pred_630arm.h5ad"), "r")
    hdn = h5py.File(restore(PKG_DN, W / "pred_dn631.h5ad"), "r")
    tg, cx = labels(hin, "target_gene"), labels(hin, "context")
    assert np.array_equal(tg, labels(hdn, "target_gene")) and np.array_equal(cx, labels(hdn, "context")), "row order differs"
    genes = RP.strings(hin["var"], "_index")
    G = len(genes)
    aud = {(r["context"], r["target"]): r for r in map(json.loads, open(AUDIT))}
    ctrl = ad.read_h5ad(ROOT / f"data/challenge/context_{CTX}.h5ad")
    assert list(map(str, ctrl.var_names)) == genes
    cc = RM.scorer_control_block(ctrl.X, min_cpm=5.0)
    del ctrl
    tie = V1.scorer_tie(cc)
    keep = np.asarray(cc.keep, dtype=bool)
    kpos = np.cumsum(keep) - 1
    gpos = {g: i for i, g in enumerate(genes)}
    cands = sorted([k for k in aud if k[0] == CTX], key=lambda k: -aud[k]["dn_moved"])
    # distinct targets: the largest-K block (most ballast rows), then the most-moved ones, then the least-moved one
    picks = [max(cands, key=lambda k: (aud[k]["ballast_K"], aud[k]["dn_moved"]))]
    for k in cands + cands[::-1]:
        if len(picks) >= n_t:
            break
        if k not in picks and (len(picks) < n_t - 1 or k == cands[-1]):
            picks.append(k)
    res = []
    for (c, t) in picks:
        idx = np.flatnonzero((cx == c) & (tg == t))
        assert idx.size == 400, (t, idx.size)
        Xin, Xdn = rows(hin, idx, G), rows(hdn, idx, G)
        n_car = 400 - int(aud[(c, t)]["ballast_K"])
        full = gpos.get(t, -1)
        own = int(kpos[full]) if full >= 0 and keep[full] else -1
        X1, d1 = V1.apply(Xin, n_car, cc, "d_m200", tie=tie, own_col=own)
        X2, d2 = V2.apply(Xin, n_car, cc, "d_m200", tie=tie, own_col=own)
        faithful = bool(np.array_equal(X1, Xdn))
        same = bool(np.array_equal(X1, X2) and d1 == d2)
        r = {"context": c, "target": t, "n_car": n_car, "dn_moved_audit": aud[(c, t)]["dn_moved"],
             "dn_moved_v1": d1["dn_moved"], "v1_equals_DN_package": faithful, "v2_equals_v1": same}
        res.append(r)
        print(json.dumps(r), flush=True)
    ok = all(r["v1_equals_DN_package"] and r["v2_equals_v1"] for r in res)
    json.dump(res, open(W / "abc693.json", "w"), indent=1)
    print(("PASS" if ok else "FAIL") + f": {len(res)} real A/B/C blocks (context {CTX}); v1(630arm block) == the DN build's "
          "block, and v2 == v1")
    hin.close(); hdn.close()
    if "--clean" in sys.argv:
        for f in ("pred_630arm.h5ad", "pred_dn631.h5ad"):
            (W / f).unlink(missing_ok=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
