"""810 -- the composition PARTS for A/B/C's folds under a chosen held-out operator, through 801's dump518p path (EC2 only).

dump518p.py (experiments/801-lab-folds, the ABCp builder) composes A/B/C's PLOLO fold on lib420's path (the board's A/B/C rules:
CD4 x1.5 all-rows sizing, SVD-5 re-fitted, H1 x2) with the held-out operator swapped in-process. This script runs the same
calls, verbatim up to `fST`, and then splits fST into decomp810's four parts with compose_arm's own algebra (lib418):
    fST = s * (W + f_cd4 C + f_h1 H + f_p P)   ->   base = sW, cd4 = s f_cd4 C, h1 = s f_h1 H, part = s f_p P
plus `pre` (the normrel'd donor mean the whitener received, panel zeroed) and 707's U (decomp810.unit_mean over the fold's
sources). Output: $PARTS810/<ABC|ABCp>__<held>-parts.npz in price810's format, and a json record.

  --ops 062   KA (abort unless < 1e-12): parts sum == compose_arm's fST == 518's stored fST; base, cd4, h1, part, U == the GO
              run's ABC parts (810parts, decomp810's lib801 route); pre == them on every target with a donor row.
  --ops day   the DAY-recipe operator (801 work/ABCp/ops/held_<held>.npz, via plolo's production branch, as dump518p);
              KA: parts sum == fST == 801's ABCp dump fST (work/dumps/abcp/dump518-<held>.npz).

    PARTS810=/work/parts810 uv run python experiments/810-strong-direction/decompabcp810.py <held> 062|day
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
from pathlib import Path  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/420-bulk-stack"))
sys.path.insert(0, str(ROOT / "experiments/518-bulk-mse-gap"))
sys.path.insert(0, str(ROOT / "experiments/801-lab-folds"))
sys.path.insert(0, str(HERE))
import guard801  # noqa: E402,F401  (writes inside the repo outside 801 raise; this script writes only outside the repo)
import hashlib  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import lib420  # noqa: E402
from lib420 import (COMPS, E062, FOLDS, Donors, Operator, build_fold, cd4_frozen_and_partial, compose_arm,  # noqa: E402
                    compose_pp, deflate, pre_normrel, source_residual, union_rows)
import lib418  # noqa: E402
import lodo415  # noqa: E402
import plolo  # noqa: E402
import decomp810 as D  # noqa: E402

HELD, OPS = sys.argv[1], sys.argv[2]
assert OPS in ("062", "day"), OPS
PARTS = Path(os.environ.get("PARTS810", "/work/parts810"))
PLOLO = Path(os.environ.get("PLOLO810", "/work/plolo810")) / OPS        # fold caches OUTSIDE the repo
PLOLO.mkdir(parents=True, exist_ok=True)
for c in "ABC":                                                          # frames from 420's caches, as dump518p's runs
    src = ROOT / f"experiments/420-bulk-stack/work/plolo/bctrl-challenge-{c}.npy"
    (PLOLO / src.name).write_bytes(src.read_bytes())
plolo.WORK = PLOLO
DAY_OPS = ROOT / "experiments/801-lab-folds/work/ABCp/ops"
REF = {"062": ROOT / f"experiments/518-bulk-mse-gap/work/dump518-{HELD}.npz",
       "day": ROOT / f"experiments/801-lab-folds/work/dumps/abcp/dump518-{HELD}.npz"}[OPS]
GO_PARTS = ROOT / f"experiments/810-strong-direction/work/parts/ABC__{HELD}-parts.npz"   # the GO run's parts (input mirror)
TOL = 1e-12


def main():
    t0 = time.time()
    if OPS == "day":                                                     # dump518p's operator swap, verbatim
        plolo.FOLDS[HELD] = dict(plolo.FOLDS[HELD], op=DAY_OPS / f"held_{HELD}.npz", op_kind="production")
    build_fold(HELD)
    fold = build_fold(HELD)
    pp_pos, wcols = np.asarray(fold["panel_positions"]), np.asarray(fold["wcols"])
    pre_cached = np.asarray(fold["pre"], dtype=np.float64)
    targets = list(fold["targets"])
    sources = list(FOLDS[HELD]["sources"])
    op = Operator(FOLDS[HELD]["op"] if OPS == "day" else
                  (lib420.PROD_OPS / "whitener_union_pruned.npz" if FOLDS[HELD]["op_kind"] == "production" else FOLDS[HELD]["op"]))
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in lodo415.UNION])
    dn = Donors(targets, sources)
    frozen, partial = cd4_frozen_and_partial(targets)
    cd4, _ = deflate(frozen, 5)
    pre228 = dn.pre(sources, False)
    pre228[:, pp_pos] = 0.0
    nocov = sum(dn.A[s].astype(int) for s in sources) == 0
    pre228[nocov] = pre_cached[nocov]
    pre = pre_normrel(dn, sources, "cap", -1.5, fallback_rows=pre228.copy())
    whit = op.apply(pre * lib420.ALPHA, wcols, lib420.GAMMA, lib420.LAM)
    mc, gh = COMPS["stack"]
    fST, _ = compose_arm(whit, cd4, partial, h1, wcols, pp_pos, m_cd4=mc, g_h1=gh)
    # ---- the parts: compose_arm's algebra (lib418.compose_arm, abc rules), kept apart
    f_cd4, _ = lib418.solve_fraction(whit, cd4, lib418.CD4_F * mc)
    comp = {"cd4": np.where(np.asarray(cd4.gene_votes) > 0, np.asarray(cd4.mean_effect), 0.0) * f_cd4}
    support = np.asarray(cd4.gene_votes) > 0
    frac = {"cd4": float(f_cd4)}
    for tag, res, g in (("h1", h1, gh), ("part", partial, lib418.PARTIAL_G)):
        if not g:
            comp[tag] = np.zeros_like(whit)
            frac[tag] = 0.0
            continue
        f, _ = lib418.solve_fraction_carrying(whit, res, g)
        v = np.asarray(res.gene_votes) > 0
        comp[tag] = np.where(v, np.asarray(res.mean_effect), 0.0) * f
        support |= v
        frac[tag] = float(f)
    cand = np.asarray(whit, np.float64).copy()
    merged = comp["cd4"] + comp["h1"] + comp["part"]
    cand[support] += merged[support]
    s = (float(np.median(np.linalg.norm(whit[:, wcols], axis=1))) or 1.0) / \
        (float(np.median(np.linalg.norm(cand[:, wcols], axis=1))) or 1.0)
    parts = {"base": s * np.asarray(whit, np.float64), **{k: s * v for k, v in comp.items()}}
    for k in parts:
        parts[k][:, pp_pos] = 0.0
    fin = parts["base"] + parts["cd4"] + parts["h1"] + parts["part"]
    pre0 = np.asarray(pre, np.float64).copy()
    pre0[:, pp_pos] = 0.0
    nonpanel = np.ones(len(D.GENES), bool)
    nonpanel[pp_pos] = False
    U = D.unit_mean([D.raw_rows(b, targets)[0] for b in sources], nonpanel)
    ref = np.load(REF)
    assert [str(t) for t in ref["targets"]] == targets
    ka = {"parts sum == compose_arm fST": float(np.max(np.abs(fin - fST))),
          f"fST == stored dump fST ({REF.relative_to(ROOT)})": float(np.max(np.abs(fST - ref["fST"]))),
          "parts sum == stored dump fST": float(np.max(np.abs(fin - ref["fST"])))}
    if OPS == "062":
        G = np.load(GO_PARTS)
        assert [str(t) for t in G["targets"]] == targets
        for k in ("base", "cd4", "h1", "part", "U"):
            ka[f"{k} == GO-run ABC parts"] = float(np.max(np.abs((U if k == "U" else parts[k]) - G[k])))
        have = np.linalg.norm(U, axis=1) > 0
        ka["pre == GO-run ABC parts (targets with a donor row)"] = float(np.max(np.abs(pre0[have] - G["pre"][have])))
    for k, d in ka.items():
        print(f"[{OPS}/{HELD}] KA {k}: {d:.2e}", flush=True)
        assert d < TOL, f"KA FAILED {OPS}/{HELD} {k}: {d}"
    stem = f"{'ABCp' if OPS == 'day' else 'ABC'}__{HELD}"
    PARTS.mkdir(parents=True, exist_ok=True)
    out = PARTS / f"{stem}-parts.npz"
    np.savez_compressed(out, targets=np.asarray(targets), scored=np.asarray(ref["scored"]), pp_pos=pp_pos, **parts, U=U, pre=pre0)
    rec = dict(held=HELD, ops=OPS, op=str(op.path.relative_to(ROOT)), op_sha256=op.sha256, ref_dump=str(REF.relative_to(ROOT)),
               ka=ka, fractions=frac, s=s, fin_sha256=hashlib.sha256(np.ascontiguousarray(fin).tobytes()).hexdigest(),
               peak_rss_gib=E062.peak_rss_gib(), seconds=time.time() - t0)
    (PARTS / f"{stem}-parts.json").write_text(json.dumps(rec, indent=1, default=float))
    print(f"[{OPS}/{HELD}] parts -> {out} ({time.time() - t0:.0f}s, op {op.sha256[:16]})", flush=True)


if __name__ == "__main__":
    main()
