#!/usr/bin/env python3
"""705: the A/B/C bulk dump with 522's θ 0.25 floor-preserving strength rider, built EXACTLY from 495's cap-6 dump.

Why this is exact. The builder writes `VCC_BULK_DUMP` right after `apply_pert_scale` (make_submission_pm_cd4.py, the
`VCC_PERT_SCALE` block), with nothing in between when VCC_BULK_GENE_GAIN is unset (it is unset in 495's dump env):
    bulk     = carriers * (w * c)[:, None],   c = sqrt(sum(row_e2) / sum(w^2 row_e2)),
    carriers = the unscaled shared effect matrix,  row_e2 = per-row energy over the non-panel columns.
So a dump for any other weight vector w' is `carriers * (w' * c')[:, None]` with the same arithmetic. No builder edit.

Known answers (abort unless exact), BEFORE the rider dump is written:
  KA-sha   : the cap-6 dump is d58750b8… (ALLSW's VCC_BULK_LOAD) and its key is 495's dc374281…;
  KA-cap6  : recomputing bulk from the dump's carriers with dc374281's w reproduces the dump's bulk np.array_equal;
  KA-rider : the rider artifact is def_strength_key522.py's θ 0.25 output on dc374281 (sha 392f1b93…, its own --ka PASS),
             and every array but w equals dc374281's.
Writes work/dump_h1x2_cd4x15_lodo05hi6_s025.npz (bulk, carriers, perts, genes, env byte-copied from the cap-6 dump) and a
provenance sidecar <out>.prov705.json (the rider's per-target multiplier, energy check, shas).
usage: uv run python experiments/705-lodo-strength-safety/make_rider_dump705.py
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CAP6_DUMP = (ROOT / "experiments/495-lodo-cap6/work/dump_h1x2_cd4x15_lodo05hi6.npz",
             "d58750b83915d40d488125926614c83793bd6d7260ceef1192b4b10c19eeb45c")
CAP6_KEY = (ROOT / "experiments/495-lodo-cap6/work/pert_scale_415pp_lo0.5_hi6.npz",
            "dc374281452ed4f7ae3a25b392dbfa4ff0d5cc591bd7b5b71e2c0a2dca86a4c0")
RIDER_KEY = (HERE / "work/pert_scale_415pp_lo0.5_hi6_s025.npz",
             "392f1b9348a7603bff4becf3753c6999844f401ef6cbd0d5b677bad714b4cb55")
OUT = HERE / "work/dump_h1x2_cd4x15_lodo05hi6_s025.npz"


def sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def scale(E: np.ndarray, w: np.ndarray, nonpanel: np.ndarray) -> tuple[np.ndarray, float, np.ndarray]:
    """apply_pert_scale's arithmetic, verbatim (same operations, same order)."""
    row_e2 = np.einsum("ij,ij->i", E[:, nonpanel], E[:, nonpanel])
    c = float(np.sqrt(row_e2.sum() / float((w * w * row_e2).sum())))
    final = w * c
    return E * final[:, None], c, row_e2


def main():
    for p, s in (CAP6_DUMP, CAP6_KEY, RIDER_KEY):
        got = sha256(p)
        assert got == s, f"KA-sha: {p.name} is {got[:16]}, expected {s[:16]}"
    print("KA-sha   PASS: cap-6 dump d58750b8, key dc374281, rider 392f1b93")
    D = np.load(CAP6_DUMP[0], allow_pickle=False)
    perts = [str(x) for x in D["perts"]]
    genes = [str(x) for x in D["genes"]]
    E = np.asarray(D["carriers"], dtype=np.float64)
    panel = set(perts)
    nonpanel = np.asarray([g not in panel for g in genes], dtype=bool)
    k6 = np.load(CAP6_KEY[0], allow_pickle=False)
    kr = np.load(RIDER_KEY[0], allow_pickle=False)
    assert [str(t) for t in k6["targets"]] == perts and [str(t) for t in kr["targets"]] == perts, "target order"
    w6 = np.asarray(k6["w"], dtype=np.float64)
    wr = np.asarray(kr["w"], dtype=np.float64)
    for k in ("targets", "key", "carried", "n_donors"):
        assert np.array_equal(np.asarray(k6[k]), np.asarray(kr[k])), f"KA-rider: {k} differs from dc374281"
    print("KA-rider PASS: every array but w equals dc374281's")
    b6, c6, row_e2 = scale(E, w6, nonpanel)
    ok = np.array_equal(b6, np.asarray(D["bulk"], dtype=np.float64))
    print(f"KA-cap6  {'PASS' if ok else 'FAIL'}: carriers x (w6 c6) == the dump's bulk (c6 {c6:.12f})")
    if not ok:
        raise SystemExit("KA-cap6 failed: the dump is not carriers x (w c); the exact construction does not hold")
    br, cr, _ = scale(E, wr, nonpanel)
    np.savez_compressed(OUT, bulk=br, carriers=np.asarray(D["carriers"]), perts=np.asarray(D["perts"]),
                        genes=np.asarray(D["genes"]), env=np.asarray(D["env"]))
    m = (wr * cr) / (w6 * c6)
    e6 = float(np.sqrt((row_e2 * (w6 * c6) ** 2).sum()))
    er = float(np.sqrt((row_e2 * (wr * cr) ** 2).sum()))
    q = lambda x: [round(float(v), 4) for v in np.quantile(x, [0, .1, .25, .5, .75, .9, 1])]
    prov = {"out": str(OUT.relative_to(ROOT)), "sha256": sha256(OUT),
            "from_dump": [str(CAP6_DUMP[0].relative_to(ROOT)), CAP6_DUMP[1]],
            "cap6_key": [str(CAP6_KEY[0].relative_to(ROOT)), CAP6_KEY[1]],
            "rider_key": [str(RIDER_KEY[0].relative_to(ROOT)), RIDER_KEY[1]],
            "rule": "bulk' = carriers x (w' c')[:, None]; w' = 522's theta 0.25 floor-preserving rider on dc374281",
            "c_cap6": c6, "c_rider": cr, "row_multiplier_quantiles": q(m),
            "rows_up": int((m > 1 + 1e-12).sum()), "rows_down": int((m < 1 - 1e-12).sum()),
            "nonpanel_frobenius_cap6": e6, "nonpanel_frobenius_rider": er,
            "env": "byte copy of the cap-6 dump's env (the dump run's VCC_* settings; the rider is recorded here)"}
    Path(str(OUT) + ".prov705.json").write_text(json.dumps(prov, indent=1))
    print(json.dumps(prov, indent=1))


if __name__ == "__main__":
    main()
