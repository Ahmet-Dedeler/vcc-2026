"""811 C1 hook, KA-same-construction (EC2): on engine folds, hook811.transform applied to composite811's B1 rows (the swapped
composition) with the fold's held group out reproduces composite811's C1 rows; and a positive per-row scale (the key) commutes
with it. Also reports the board rule (hold none) against the fold rule (held group out) on the same rows.
    uv run python experiments/811-lolo-model/kasame811.py <fold-id>      # writes $KASAME811_OUT/<fold>.json
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import composite811 as C  # noqa: E402
import hook811 as H  # noqa: E402

OUT = Path(os.environ.get("KASAME811_OUT", "/work/kasame811"))


def main(fold):
    t0 = time.time()
    F = C.E.load(fold, check=True)
    arms, info, (fin, fin707, U, npn) = C.build(F, "day")
    panel = C.M.panel_targets(F)
    held = C.B.HELD_GROUP.get(F.held)
    c1 = np.asarray(arms["C1"].fin)
    t, ti = H.transform(fin707, F.targets, C.B.GENES, panel, lam=0.25, hold=frozenset({held}))
    d_exact = float(np.max(np.abs(t - c1)))
    rng = np.random.default_rng(811)
    s = np.exp(rng.normal(0, 1, F.n))[:, None]
    ts, _ = H.transform(fin707 * s, F.targets, C.B.GENES, panel, lam=0.25, hold=frozenset({held}))
    d_scale = float(np.max(np.abs(ts - c1 * s)) / np.max(np.abs(c1 * s)))
    tb, tbi = H.transform(fin707, F.targets, C.B.GENES, panel, lam=0.25, hold=frozenset())
    ch_fold, ch_board = (c1 - fin707), (tb - fin707)
    nf = np.linalg.norm(ch_fold, axis=1)
    nbd = np.linalg.norm(ch_board, axis=1)
    both = (nf > 0) & (nbd > 0)
    cosd = [float(ch_fold[i] @ ch_board[i] / (nf[i] * nbd[i])) for i in np.flatnonzero(both)]
    rec = dict(fold=fold, held_group=held, n=F.n, ka_engine=F.ka,
               KA_same_exact_max_abs_d=d_exact, KA_same_pass=d_exact <= 1e-12,
               changed_targets_same=bool(np.array_equal(np.linalg.norm(t - fin707, axis=1) > 0, np.linalg.norm(c1 - fin707, axis=1) > 0)),
               KA_scale_max_rel_d=d_scale, KA_scale_pass=d_scale <= 1e-12,
               board_vs_fold=dict(n_moved_fold=int((nf > 0).sum()), n_moved_board=int((nbd > 0).sum()),
                                  change_dir_cos_q=np.round(np.quantile(cosd, [0, .25, .5, .75, 1]), 4).tolist() if cosd else None,
                                  board_groups=tbi["groups"], fold_groups=ti["groups"]),
               seconds=round(time.time() - t0, 1))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{fold.replace('/', '__')}.json").write_text(json.dumps(rec, indent=1, default=float))
    print(json.dumps(rec, default=float))
    assert rec["KA_same_pass"] and rec["KA_scale_pass"] and rec["changed_targets_same"], "KA FAILED"


if __name__ == "__main__":
    main(sys.argv[1])
