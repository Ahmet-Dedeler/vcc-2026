"""228 §5 — the October property: is the gain selected against TODAY'S panel?

The panel is replaced on 2026-10-22 with 300 targets disjoint from today's, and only the last
entry after that date is prize-eligible (`docs/DEF-RUNBOOK.md` §0, §5).  Every candidate in
225's source ladder is selected against today's panel -- the source, and above all the dose,
are chosen by performance on the 290/300 panel targets a fold can adjudicate.  This gain is
different by construction: it is fitted on NON-PANEL targets only.

Two measurements, because "by construction" is a claim and this file is where it is tested:

  A. SYNTHETIC PANEL.  Refit the gain with a different, disjoint 300-gene set playing the
     panel's role -- excluded from the fit-target pool and zeroed in the gene axis exactly as
     the real panel is.  If the resulting gain is the same vector, then which 300 genes are
     "the panel" does not enter the fit.  Reported as `r(gain_real, gain_synthetic)` and as
     the PLOLO price of the synthetic-panel gain on the REAL panel evaluation.

  B. ZERO-PANEL-COVERAGE FIT POOL.  224 §5's strictest control: fit only on `rpe1`,
     `nadig_jurkat`, `nadig_hepg2`, which carry 0 of the 300 panel targets, so the gain cannot
     have seen a panel response in any line at any point.

WHAT THIS DOES NOT COVER.  The whitening operator's gene axis is `union(donor coverage) -
panel` and MUST be re-fitted on 2026-10-22 (AGENTS.md §8.3, `docs/DEF-RUNBOOK.md`).  That
refit is 437 s and 8.1 GiB and is the runbook's job, not this experiment's; the operator here
is held at today's.

    uv run python october228.py
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
import gain228  # noqa: E402
from plolo import OFFICIAL, PANEL, PANEL_SET  # noqa: E402

SEED = 771103
HELD = ("k562_gwps", "orion_hct116", "orion_hek293t")


class Args:
    seed = gain228.SEED


def synthetic_panel(rng) -> list[str]:
    """300 genes disjoint from the real panel, drawn to match its basal-expression profile.

    The real panel was selected for high basal expression of the target (AGENTS.md §2.1), so
    a uniformly random draw would be an easier panel; the draw is stratified by the challenge
    context-A basal decile of the real panel's own targets.
    """
    b = plolo.challenge_bctrl("A")
    pos = {g: i for i, g in enumerate(OFFICIAL)}
    real = [g for g in PANEL if g in pos]
    pool = [g for g in OFFICIAL if g not in PANEL_SET]
    bp = np.asarray([b[pos[g]] for g in pool])
    edges = np.quantile([b[pos[g]] for g in real], np.linspace(0, 1, 11))
    out: list[str] = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        want = sum(1 for g in real if lo <= b[pos[g]] <= hi)
        cand = [pool[i] for i in np.flatnonzero((bp >= lo) & (bp <= hi))
                if pool[i] not in out]
        if not cand:
            continue
        take = rng.choice(len(cand), size=min(want, len(cand)), replace=False)
        out += [cand[i] for i in take]
    while len(out) < 300:
        g = pool[int(rng.integers(len(pool)))]
        if g not in out:
            out.append(g)
    return out[:300]


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--held", nargs="*", default=list(HELD))
    ap.add_argument("--variant", choices=["synthetic", "zeropanel"], required=True)
    cli = ap.parse_args()

    rng = np.random.default_rng(SEED)
    synth = synthetic_panel(rng)
    assert not (set(synth) & PANEL_SET), "synthetic panel overlaps the real one"
    out_path = WORK / "october228.json"
    out = json.loads(out_path.read_text()) if out_path.exists() else {}
    out["synthetic_panel"], out["seed"] = synth, SEED
    print(f"synthetic panel: {len(synth)} genes, 0 shared with the real 300", flush=True)

    real_panel, real_set, pool_real = gain228.PANEL, gain228.PANEL_SET, gain228.POOL
    for held in cli.held:
        _z = np.load(WORK / f"gain-{held}.npz", allow_pickle=True)
        z = {k: np.array(_z[k]) for k in _z.files}
        _z.close()
        live = np.asarray(z["live"], bool)
        src = WORK / f"gain-{held}.npz"

        if cli.variant == "synthetic":
            # BOTH panels are withheld: the real one because the evaluation is on it and a
            # leak would be fatal, the synthetic one because it is the variable.  What this
            # measures is whether swapping WHICH 300 targets/genes are withheld moves the
            # gain -- which is exactly what happens on 2026-10-22.
            gain228.PANEL = tuple(real_panel) + tuple(synth)
            gain228.PANEL_SET = set(real_set) | set(synth)
            try:
                r = gain228.fit_one(held, Args())
            finally:
                gain228.PANEL, gain228.PANEL_SET = real_panel, real_set
            dst = WORK / f"gain-{held}-synthpanel.npz"
        else:
            gain228.POOL = {k: v for k, v in pool_real.items()
                            if k in ("rpe1", "nadig_jurkat", "nadig_hepg2")}
            try:
                r = gain228.fit_one(held, Args())
            finally:
                gain228.POOL = pool_real
            dst = WORK / f"gain-{held}-zeropanel.npz"

        src.replace(dst)
        np.savez_compressed(src, **z)          # restore the real gain, untouched
        zs = np.load(dst, allow_pickle=True)
        m = live & np.asarray(zs["live"], bool)
        cor = {k: float(np.corrcoef(np.asarray(z[k])[m], np.asarray(zs[k])[m])[0, 1])
               for k in ("gain", "gain_reg")}
        out[f"{held}|{cli.variant}"] = dict(
            fit_lines=r["fit_lines"], n_fit_targets=r["n_fit_targets"],
            shrink=r["shrink"], shrink_pooled=r["shrink_pooled"],
            corr_with_real=cor, n_shared_live=int(m.sum()))
        print(f"[{held}] {cli.variant}: {len(r['fit_lines'])} lines, "
              f"{r['n_fit_targets']} fit targets, r(real, variant) = "
              f"{cor['gain']:+.5f} (gain) / {cor['gain_reg']:+.5f} (gain_reg)", flush=True)
        out_path.write_text(json.dumps(out, indent=2))

    out_path.write_text(json.dumps(out, indent=2))
    print(f"\n-> {out_path}")


if __name__ == "__main__":
    main()
