#!/usr/bin/env python3
"""553 — builds the synthetic `fx` surface for an end-to-end access test in seconds (test only): test553's big fixture as
the truth, its control reference, and three fake "arms" x 3 seeds of predicted-block statistics (resampled truth cells
with a per-arm tilt, 400 cells per target) computed with evaluate323.block_stats against the fixture's controls.
usage: uv run python experiments/553-r7b-harness/fx553.py && uv run python experiments/553-r7b-harness/access553.py fx \
         --arms armA,armB --seeds 0,1,2
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lib553 as S  # noqa: E402
import test553  # noqa: E402

TAG = "base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse"


def main():
    import evaluate323 as E
    import scipy.sparse as sp
    spec = S.surface_spec("fx")
    d = Path(spec["truth_h5ad"]).parent
    d.mkdir(parents=True, exist_ok=True)
    a = test553.fixture_big(seed=5531)
    a.write_h5ad(spec["truth_h5ad"])
    cr = S.controlref(Path(spec["control_h5ad"]))
    lab = a.obs.target_gene.to_numpy().astype(str)
    targets = np.array(sorted(set(lab) - {S.CTRL}))
    np.savez(spec["controlref"], targets=targets, **cr)
    S.write_json(Path(spec["controlref_json"]), dict(statistics_sha256=S.sha256(Path(spec["controlref"]))))
    from vcc2026 import rc_emit
    c = rc_emit.RcControlBlock(keep=cr["keep"], cs=cr["cs"], tie=cr["tie"], ref_lin=cr["ref_lin"],
                               n_controls=int(cr["n_controls"]), min_cpm=5)
    X = a.X.tocsr()
    excl = sys.argv[sys.argv.index("--exclude") + 1].split(",") if "--exclude" in sys.argv else []
    if excl:
        targets = np.array([t for t in targets if t not in excl])
    for arm, tilt in (("cpm", 1.0), ("armA", 1.3), ("armB", 0.8)):
        out = S.WORK / (f"fx_{arm}" + (("-x" + "+".join(sorted(excl))) if excl else ""))
        out.mkdir(parents=True, exist_ok=True)
        for seed in (0, 1, 2):
            rng = np.random.default_rng(100 * seed + 7)
            st = []
            for t in targets:
                rows = np.flatnonzero(lab == t)
                blk = np.asarray(X[rng.choice(rows, 400, replace=True)].toarray())
                blk = np.floor(blk * np.where(rng.random(blk.shape[1]) < 0.3, tilt, 1.0) + rng.random(blk.shape)).astype(np.float64)
                st.append(E.block_stats(sp.csr_matrix(blk), c))
            s_, j_, p_, q_, l_ = (list(v) for v in zip(*st))
            for zs in (False, True):
                t_ = TAG.format(seed=seed) + ("_dnk250zs" if zs else "")
                np.savez(out / f"stats-{t_}.npz", targets=targets, sums=np.asarray(s_), jk=np.asarray(j_), p=np.asarray(p_),
                         q=np.asarray(q_), lfc=np.asarray(l_))
    print(json.dumps(dict(targets=len(targets), controls=int(cr["n_controls"]))))


if __name__ == "__main__":
    main()
