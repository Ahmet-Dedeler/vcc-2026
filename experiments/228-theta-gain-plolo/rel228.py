"""228 — the mechanism triangulation 224 ran, re-run on THIS axis.

224 §7's strongest predictor of the fitted gain was its per-gene CROSS-LINE RELIABILITY
(`r = +0.938`): the mean, over line pairs, of the correlation across targets between the two
lines' unit-normalised effects at that gene.  This recomputes exactly that statistic on the
official 18,533-gene axis over the NON-PANEL fit targets, and correlates it with the gain
this experiment fitted, so the "it is 189's `theta` in the direction channel" claim is
tested on my own instrument rather than inherited.

    uv run python rel228.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(HERE))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
from plolo import E062, FOLDS, OFFICIAL, PANEL, PANEL_SET, build_fold  # noqa: E402
from gain228 import HELD_GROUP, POOL, read_perts, unit32  # noqa: E402

HELD = "kolf2"          # the fold whose fit pool is all six lines


def spearman(a, b):
    return float(np.corrcoef(np.argsort(np.argsort(a)), np.argsort(np.argsort(b)))[0, 1])


def main() -> None:
    fold = build_fold(HELD)
    fit_lines = [n for n, (_, g) in POOL.items() if g != HELD_GROUP[HELD]]
    fit_targets = sorted(set.intersection(*[read_perts(POOL[n][0]) for n in fit_lines])
                         - PANEL_SET)
    official_pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([official_pos[g] for g in PANEL if g in official_pos],
                                 dtype=np.int64)
    G = len(OFFICIAL)

    U, cov = {}, {}
    for n in fit_lines:
        bank = E062.Bank.load(n, ROOT / POOL[n][0])
        amap = E062.map_bank_to_axis(bank, OFFICIAL)
        c = np.zeros(G, dtype=bool)
        c[amap.unique_target_positions] = True
        c[panel_positions] = False
        cov[n] = c
        p = E062.prepare_main_source(bank, fit_targets, OFFICIAL, panel_positions,
                                     chunk_rows=256, max_rss_gib=9.0)
        U[n] = unit32(np.asarray(p.effects, dtype=np.float64))
        del bank, p
        print(f"  loaded {n}", flush=True)

    rel = np.zeros(G); npair = np.zeros(G)
    sd_line = np.zeros(G); nline = np.zeros(G)
    for i, a in enumerate(fit_lines):
        A = U[a].astype(np.float64); A = A - A.mean(0, keepdims=True)
        na = np.linalg.norm(A, axis=0)
        sd_line += A.std(0) * cov[a]; nline += cov[a]
        for b in fit_lines[i + 1:]:
            B = U[b].astype(np.float64); B = B - B.mean(0, keepdims=True)
            nb = np.linalg.norm(B, axis=0)
            m = cov[a] & cov[b]
            rel += np.where(m, (A * B).sum(0) / np.maximum(na * nb, 1e-300), 0.0)
            npair += m
    rel = np.where(npair > 0, rel / np.maximum(npair, 1), np.nan)
    sd_line = np.where(nline > 0, sd_line / np.maximum(nline, 1), np.nan)
    cons = np.mean([U[n].astype(np.float64) for n in fit_lines], axis=0)
    sd_cons = cons.std(0)
    b_ctrl = np.asarray(fold["b_ctrl"], dtype=np.float64)

    out = {}
    print("\nr(gain, predictor) on the official axis, non-panel fit targets, "
          f"held {HELD}, {len(fit_lines)} lines, {len(fit_targets)} targets")
    for key in ("gain", "gain_reg"):
        z = np.load(WORK / f"gain-{HELD}.npz", allow_pickle=True)
        g = np.asarray(z[key], dtype=np.float64)
        live = np.asarray(z["live"], bool)
        full = np.asarray(z["full_cov"], bool)
        for nm, v, mask in (("cross_line_reliability", rel, live),
                            ("cross_line_reliability|fullcov", rel, live & full),
                            ("sd_single_line", sd_line, live),
                            ("sd_consensus", sd_cons, live),
                            ("sd_ratio", sd_line / np.maximum(sd_cons, 1e-300), live),
                            ("basal_b_ctrl", b_ctrl, live)):
            ok = mask & np.isfinite(v) & np.isfinite(g)
            r = float(np.corrcoef(g[ok], v[ok])[0, 1])
            s = spearman(g[ok], v[ok])
            out[f"{key}|{nm}"] = dict(pearson=r, spearman=s, n=int(ok.sum()))
            print(f"   {key:9s} {nm:32s} r {r:+.4f}  spearman {s:+.4f}  n {int(ok.sum())}")

    np.savez_compressed(WORK / "rel228.npz", rel=rel, sd_line=sd_line, sd_cons=sd_cons)
    (WORK / "rel228.json").write_text(json.dumps(out, indent=2))
    print(f"\n-> {WORK / 'rel228.json'}")


if __name__ == "__main__":
    main()
