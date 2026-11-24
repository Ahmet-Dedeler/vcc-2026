"""277 -- V4's own controls: the 272 A' eigenvalue floor (c = 1) applied to DIFFERENT row sets, and to all
rows (272's own arm, the known-answer row).

    PYTHONDONTWRITEBYTECODE=1 uv run python floor277.py --fold <cancer fold>

The expensive object is the transformed base `T_w = ((pre[:, wcols] @ V) * w) @ V^T` for one weight vector
w; every ROW SET is then a splice of `T_flr` into `T_inc` followed by production's single on-axis preserve,
so three matmuls cover every cell.  Arms: the floor on CLASS30 (= V4), on TOPCONS30, on three size-matched
RANDOM sets, on ALL rows (272's `floor:c1`), and -- the mechanism controls -- the floor built from 272's
SHUFFLED and MIRRORED noise vectors on CLASS30.  Each arm re-solves both residual fractions against its own
base (263's `arm_parts`, which is what 272 did).

KA: the ALL-rows cell reproduces 272's `floor:c1` per-fold d pds at 1e-9.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib277 as Z  # noqa: E402

L, plolo, O = Z.L, Z.plolo, Z.O


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", required=True)
    held = ap.parse_args().fold
    t0 = time.time()
    st = L.stage(held)
    targets = [str(t) for t in st["targets"]]
    n = len(targets)
    fold = st["fold"]
    wc = fold["wcols"]
    W = np.asarray(st["whitened"], dtype=np.float64)
    d5, _ = L.deflate(st["cd4"], "svd", 5)
    cl = Z.c_lolo(held, targets)
    C30 = Z.bottom_class(cl, targets, 0.30)
    n_class = int(C30.sum())
    sets = {"CLASS30": C30, "TOPCONS30": Z.bottom_class(cl, targets, 0.30, top=True),
            "ALL": np.ones(n, dtype=bool), "NONE": np.zeros(n, dtype=bool)}
    for i, s in enumerate((20260906, 20260907, 20260908)):
        sets[f"RAND{i+1}"] = Z.random_class(cl, n_class, s)

    spz_path = Z.E272 / (f"work/spectrum272-{held}.npz" if plolo.FOLDS[held]["op_kind"] != "production"
                         else "work/spectrum272-production.npz")
    spz = np.load(spz_path)
    op = plolo.Operator(L.op_path_for(held))
    assert np.array_equal(spz["evals"], op.evals)
    nvec = np.asarray(spz["n"], dtype=np.float64)
    weights = {"true": nvec, "shuffled": O.shuffled_s(nvec, O.SHUF_SEED),
               "mirrored": O.mirrored_s(nvec, op.evals)}
    w_inc = (op.evals + Z.LAM) ** (-Z.GAMMA / 2.0)
    pre = np.asarray(st["pre"], dtype=np.float64) * L.ALPHA
    C = pre[:, wc] @ op.evecs
    T_inc = (C * w_inc) @ op.evecs.T
    before = float(np.median(np.linalg.norm(pre[:, wc], axis=1)))
    print(f"[{held}] {n} targets, CLASS30 {n_class}; base matmul done ({time.time()-t0:.0f}s, peak "
          f"{Z.peak_gib():.2f} GiB)", flush=True)

    def compose_from(T_mixed, tag):
        out = pre.copy()
        after = float(np.median(np.linalg.norm(T_mixed, axis=1)))
        out[:, wc] = T_mixed * (before / after)
        st_arm = dict(st, whitened=out)
        fr, _ = L.solve_fraction(out, d5, Z.CD4_F)
        fx, _ = L.solve_fraction_carrying(out, st["union"], L.H1_G)
        e = L.score(st_arm, L.compose(st_arm, [(d5, fr), (st["union"], fx)]))
        return dict(arm=tag, pds=float(e["pds"]), cos_agg=float(e["cos_agg"]), a23_cos=float(e["a23_cos"]),
                    cd4_fraction=float(fr), union_fraction=float(fx)), \
            np.asarray(e["_pds_per_target_a23"], dtype=np.float64)

    rows, pts = {}, {}
    r0, p0 = compose_from(T_inc, "V0_cd4svd5")
    rows["V0_cd4svd5"], pts["V0_cd4svd5"] = r0, p0
    ka0 = abs(r0["pds"] - json.loads((Z.WORK / f"render277-{held}.json").read_text())
              ["arms"]["V0_cd4svd5"]["pds"])
    print(f"  KA-0 the spliced path with NO row edited == render277's V0: |d pds| {ka0:.1e} "
          f"-> {'PASS' if ka0 < 1e-12 else 'FAIL'}", flush=True)
    for wname, nv in weights.items():
        w_flr = (np.maximum(op.evals, 1.0 * nv) + Z.LAM) ** (-Z.GAMMA / 2.0)
        T_flr = (C * w_flr) @ op.evecs.T
        use = sets if wname == "true" else {"CLASS30": C30}
        for sname, m in use.items():
            if sname == "NONE":
                continue
            tag = f"floor:{wname}:{sname}"
            T = np.where(m[:, None], T_flr, T_inc)
            r, p = compose_from(T, tag)
            r["d_pds"] = r["pds"] - r0["pds"]
            r["rel_agg"] = r["cos_agg"] / r0["cos_agg"] - 1.0
            r["d_pds_on_rows"] = float((p - p0)[m].mean()) if m.any() else 0.0
            rows[tag], pts[tag] = r, p
            print(f"  [{tag:26s}] pds {r['pds']:.6f} (d {r['d_pds']:+.5f}) cos_agg {r['cos_agg']:.6f} "
                  f"({100*r['rel_agg']:+.2f} %) | on the edited rows {r['d_pds_on_rows']:+.5f}", flush=True)
        del T_flr

    # POST-HOC (not in DEFINITION-277): the floor's effect on the class survives shuffling the noise
    # vector, so the mechanism cannot be the noise ESTIMATE.  The natural alternative -- "these rows want
    # LESS whitening" -- is a lower gamma on the class rows only.  Same splice, one matmul per gamma.
    for g in (0.50, 0.35):
        w_g = (op.evals + Z.LAM) ** (-g / 2.0)
        T_g = (C * w_g) @ op.evecs.T
        for sname in ("CLASS30", "TOPCONS30", "RAND1", "ALL"):
            tag = f"gamma{g:g}:{sname}"
            m = sets[sname]
            r, p = compose_from(np.where(m[:, None], T_g, T_inc), tag)
            r["d_pds"] = r["pds"] - r0["pds"]
            r["rel_agg"] = r["cos_agg"] / r0["cos_agg"] - 1.0
            r["d_pds_on_rows"] = float((p - p0)[m].mean()) if m.any() else 0.0
            rows[tag], pts[tag] = r, p
            print(f"  [{tag:26s}] pds {r['pds']:.6f} (d {r['d_pds']:+.5f}) cos_agg {r['cos_agg']:.6f} "
                  f"({100*r['rel_agg']:+.2f} %) | on the edited rows {r['d_pds_on_rows']:+.5f}", flush=True)
        del T_g
    del op, C, T_inc

    # KA: the ALL-rows cell is 272's floor:c1
    p272 = json.loads((Z.E272 / "work/price272.json").read_text())["A"]["floor:c1"]["vs_cd4svd5"]
    i272 = Z.CANCER.index(held) if held in Z.CANCER else None
    ka1 = (abs(rows["floor:true:ALL"]["d_pds"] - p272["per_fold_d_pds"][i272])
           if i272 is not None else float("nan"))
    print(f"  KA-1 floor on ALL rows == 272's `floor:c1` on this fold: |d| {ka1:.1e} "
          f"-> {'PASS' if (ka1 < 1e-9) else 'FAIL'}", flush=True)

    out = dict(fold=held, n_targets=n, n_class=n_class, rows=rows,
               ka=dict(spliced_V0=float(ka0), vs_272_floor_c1=float(ka1)),
               seconds=round(time.time() - t0, 1), peak_rss_gib=Z.peak_gib())
    np.savez_compressed(Z.WORK / f"floor277-{held}.npz", **{f"{k}|pds": v for k, v in pts.items()},
                        targets=np.asarray(targets), **{f"mask|{k}": v for k, v in sets.items()})
    (Z.WORK / f"floor277-{held}.json").write_text(json.dumps(out, indent=1, default=float))
    if not (ka0 < 1e-12 and (np.isnan(ka1) or ka1 < 1e-9)):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
