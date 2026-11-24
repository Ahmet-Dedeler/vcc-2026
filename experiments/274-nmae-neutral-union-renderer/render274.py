"""274 -- an `nmae`-NEUTRAL RENDERING of the six-source union's DIRECTIONS on top of `cd4svd5`.
MEASUREMENT ONLY: no build, no `vcc`, no `cell-eval2`, nothing under `submissions/`, nothing written outside this
directory.  258's instrument (`price258.py`) with ONE new stage, the renderer, inserted at the one place production
would need it: AFTER the merged add (`candidate = whitened + merged`, the first line of
`E062.add_post_whiten_residual`) and BEFORE the global preserve scalar.

THE RENDERER (frozen before any fold was measured).  Let X_ref = wh + merged_ref be the INCUMBENT's pre-preserve
composition (CD4 deflated k=5 at f=0.34 + the H1 pair at g=0.28) and X_arm = wh + merged_arm the union arm's
(same CD4, the union at g).  On every row the arm changed (the union's carried rows; all other rows pass through
bit-for-bit), a magnitude-control set S and a matching mode produce X_out:
    set   wald89   the identity's truth pool (top-89 by |lfc| sqrt(cpm))     -- ORACLE, the ceiling; not buildable
          ownK     the incumbent's own top-K entries by |X_ref| (K = 89 / 400 / 2000)  -- prediction-defined, buildable
          supp     the entries the union changed (X_arm != X_ref)            -- buildable
          all      the whole row                                             -- buildable
    mode  ratio    X_out[S] = X_arm[S] * (sum|X_ref[S]| / sum|X_arm[S]|)  -- L1 total on S held at the incumbent's
          norm     X_out[S] = X_arm[S] * (||X_ref[S]|| / ||X_arm[S]||)    -- L2 norm on S held (all+norm = R2a)
          rank     |X_out[S]| = sorted |X_ref[S]| placed at the ranks of |X_arm[S]|, signs of X_arm  (all+rank = R2b:
                   the union's ORDER and SIGNS with the incumbent's magnitude spectrum, per row)
    off-S entries are left as composed.  Then production's ONE global preserve (median row norm over the operator
    columns restored to the whitened base's) and the panel columns zeroed, exactly as the builder does.
R1 = ownK / wald89 sets, R2 = all, R3 = the same on stack3 (pair + dld1) at g 0.28/0.35 and stack6 at g 0.20/0.28,
R4 = 259's construction B (pair kept at 0.28 + union2) rendered.  Controls: the renderer with X_arm := X_ref (must
reproduce the incumbent at 0.0), the union at its own magnitudes (258's rows, bit-identical), and energy-matched
RANDOM directions rendered the same way (the union's per-row contribution d = X_arm - X_ref replaced by a gene
permutation of itself (PERM) or by another carried row's d (SHUF, 218's derangement)).

    uv run python experiments/274-nmae-neutral-union-renderer/render274.py --fold k562_gwps   (x4 folds, one process each)
    uv run python experiments/274-nmae-neutral-union-renderer/render274.py --prod
    uv run python experiments/274-nmae-neutral-union-renderer/render274.py --price
"""
from __future__ import annotations

import argparse
import contextlib
import gc
import importlib.util
import io
import json
import os
import resource
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ["VCC_ALPHA"] = "0.982"
os.environ["VCC_CD4_TARGET_F"] = "0.34"
os.environ["VCC_H1_UNION_G"] = "0.28"
os.environ["VCC_CD4_DEFLATE_K"] = "5"
for _k in ("VCC_CD4_SIZING", "VCC_CD4_TARGET_G", "VCC_UNION_SOURCES", "VCC_UNION2_SOURCES", "VCC_UNION2_G",
           "VCC_CD4_DEFLATE_BASIS", "VCC_CD4_PARTIAL_G", "VCC_CD4_FRACTION", "VCC_SIGN_REPAIR"):
    os.environ.pop(_k, None)

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
E229 = ROOT / "experiments/229-stack4-build"
E243 = ROOT / "experiments/243-transfer-factor"
E246 = ROOT / "experiments/246-direction-hypotheses"
E248 = ROOT / "experiments/248-second-residual"
E249 = ROOT / "experiments/249-cd4-deflate"
E255 = ROOT / "experiments/255-cd4-deflate-k10"
E258 = ROOT / "experiments/258-union-dose-on-cd4svd5"
E259 = ROOT / "experiments/259-union-without-igvf"
E262 = ROOT / "experiments/262-where-pds-is-lost"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/188-magnitude-calibration"))
sys.path.insert(0, str(ROOT / "experiments/203-reopen-125"))
sys.path.insert(0, str(E243))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"        # byte copies of 258's (= 249's) fold caches, SHA256SUMS.txt checked

from plolo import (E062, FOLDS, OFFICIAL, PANEL, PROD_OPS, Operator, axis_sha256,  # noqa: E402
                   build_fold, cd4_for, deliver, endpoints)
from mapfit import to_lfc  # noqa: E402
from phase23 import POOL_N, key_precision, topn, umd  # noqa: E402
from lib243 import CP46, KC, SUMD  # noqa: E402
from cell_eval2.metrics.delta import _across_pert_budget, _row_weights  # noqa: E402

ALPHA, GAMMA, LAM = 0.982, 0.65, 0.0
NMAE_SPAN = 0.60                 # 188/216/229's identity span (overall = -d nmae / 0.60 / 6)
SPAN_MSE = 0.950162566           # the scorer's mse span (243/248/255/258/259)
PDS_SCALE = 2.2138               # overall pds = 2.2138 x k_pds x d pds / 6
SIGMA_PDS = 0.001494             # 162's leave-one-arm-out RMS (a calibration RMS, NOT a sampling sigma)
GATE = 0.0019
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
CTX = ("A", "B", "C")
N_GENES = 18533
POOLS = ("wald89", "lfc89", "h89")
TS = plolo.TS
# the bands (259's convention, from BR237 / BR249 / 243)
K_PDS_BAND = (0.70, 2.48)        # 243's class value (h1union realised) / ARM6-realised A-B-C mean
K_COS_BAND = (0.94, 0.96)        # the mandate: new-source band for the union legs (pessimistic / point)
K_COS_OPT = 1.00                 # the top of the new-source band (243: 0.94-1.00)
DE3_BAND = (0.0425, 0.24)        # ARM6-realised DE-three per unit board d cos / 243's recipe
NMAE_LOSS_X = 1.5                # the identity on a LOSS (BR237 realised 1.497x)
POINT = dict(k_cos=0.96, k_pds=2.48, de3=0.0425, nmae_x=1.0)
K_PDS_A_POINT = 3.79
RATE218 = dict(reach=0.002446, fid=0.000847, jac=0.000100)
DE3_SPLIT = {m: RATE218[m] / sum(RATE218.values()) for m in RATE218}
BOOT_SEED, BOOT_B = 20260906, 20_000
AMP_TOL_PCT = 0.5

PAIR = ("gse295214_h1", "h1_2025")
DLD1, JIANG, IGVF, MYL = ("dld1_lowmoi",), ("jiang_5line_noK562",), ("igvf_cardiomyocyte",), ("myllia_a549",)
SEL = {
    "pair": PAIR,
    "stack3": PAIR + DLD1,
    "stack6": PAIR + DLD1 + JIANG + IGVF + MYL,
    "u2_three": DLD1 + JIANG + MYL,                 # 259's construction B second residual
    "u2_four": DLD1 + JIANG + IGVF + MYL,           # 248's CUM4 order (259's B4+igvf)
}
INC = "cd4svd5 <- INCUMBENT"
OLD = "h1union g=0.28"
ARM6 = "stack6 g=0.35 on h1union"
R0 = "stack6 g=0.35"


# (name, deflate k, cd4 f, g1, sel1, g2, sel2, variant, rset, rmode, K)
def A(name, k, f, g1, sel1, g2=None, sel2=None, variant=None, rset=None, rmode=None, K=None):
    return (name, k, f, g1, sel1, g2, sel2, variant, rset, rmode, K)


def rname(base, rset, rmode, K=None, variant=None):
    s = f"{rset}{K}" if rset == "own" else rset
    v = f" {variant}" if variant else ""
    return f"{base}{v} | {s} {rmode}"


ARMS: list[tuple] = [
    A("a15 f=0.10 SHIPPED", None, 0.10, None, None),
    A("cd4 f=0.34", None, 0.34, None, None),
    A(OLD, None, 0.34, 0.28, "pair"),
    A(ARM6, None, 0.34, 0.35, "stack6"),
    A(INC, 5, 0.34, 0.28, "pair"),
    # R0 -- the unions at their own magnitudes (258 / 259's rows, bit-identical known answers)
    A("stack6 g=0.35", 5, 0.34, 0.35, "stack6"),
    A("stack6 g=0.28", 5, 0.34, 0.28, "stack6"),
    A("stack6 g=0.2", 5, 0.34, 0.20, "stack6"),
    A("stack3 g=0.35", 5, 0.34, 0.35, "stack3"),
    A("stack3 g=0.28", 5, 0.34, 0.28, "stack3"),
    A("B3 g2=0.35", 5, 0.34, 0.28, "pair", 0.35, "u2_three"),
    A("B3 g2=0.45", 5, 0.34, 0.28, "pair", 0.45, "u2_three"),
    A("B4+igvf g2=0.35", 5, 0.34, 0.28, "pair", 0.35, "u2_four"),
]
# R1 / R2 on the six-source union at g = 0.35: the full set x mode grid
_MAIN = ("stack6 g=0.35", 5, 0.34, 0.35, "stack6", None, None)
for rset, K in (("wald89", None), ("own", 89), ("own", 400), ("own", 2000), ("supp", None)):
    for rmode in ("ratio", "rank"):
        ARMS.append(A(rname(_MAIN[0], rset, rmode, K), *_MAIN[1:], None, rset, rmode, K))
for rmode in ("norm", "rank", "ratio"):
    ARMS.append(A(rname(_MAIN[0], "all", rmode), *_MAIN[1:], None, "all", rmode, None))
# R3 -- the cheap subset and the other doses, six renderers each
_R3 = (("stack3 g=0.35", 0.35, "stack3"), ("stack3 g=0.28", 0.28, "stack3"),
       ("stack6 g=0.28", 0.28, "stack6"), ("stack6 g=0.2", 0.20, "stack6"))
_SIX = (("wald89", "ratio", None), ("own", "ratio", 89), ("own", "ratio", 400), ("all", "norm", None),
        ("all", "rank", None), ("supp", "rank", None))
for base, g, sel in _R3:
    for rset, rmode, K in _SIX:
        ARMS.append(A(rname(base, rset, rmode, K), 5, 0.34, g, sel, None, None, None, rset, rmode, K))
# R4 -- 259's construction B rendered
for base, g2, sel2, specs in (("B3 g2=0.35", 0.35, "u2_three", (("all", "norm", None), ("all", "rank", None), ("own", "ratio", 400), ("supp", "rank", None))),
                              ("B3 g2=0.45", 0.45, "u2_three", (("all", "rank", None), ("all", "norm", None))),
                              ("B4+igvf g2=0.35", 0.35, "u2_four", (("all", "rank", None), ("own", "ratio", 400), ("wald89", "ratio", None)))):
    for rset, rmode, K in specs:
        ARMS.append(A(rname(base, rset, rmode, K), 5, 0.34, 0.28, "pair", g2, sel2, None, rset, rmode, K))
# the ORACLE ceiling proper: the truth pool held at the incumbent's POST-preserve magnitude exactly
ARMS.append(A(rname("stack6 g=0.35", "wald89", "exact"), 5, 0.34, 0.35, "stack6", None, None, None, "wald89", "exact", None))
ARMS.append(A(rname("stack6 g=0.28", "wald89", "exact"), 5, 0.34, 0.28, "stack6", None, None, None, "wald89", "exact", None))
ARMS.append(A(rname("B4+igvf g2=0.35", "wald89", "exact"), 5, 0.34, 0.28, "pair", 0.35, "u2_four", None, "wald89", "exact", None))
ARMS.append(A(rname("B3 g2=0.35", "wald89", "exact"), 5, 0.34, 0.28, "pair", 0.35, "u2_three", None, "wald89", "exact", None))
# CONTROLS -- (i) the renderer's identity check: the arm IS the incumbent's composition
for rset, rmode, K in (("wald89", "ratio", None), ("all", "rank", None), ("all", "norm", None), ("own", "ratio", 400)):
    ARMS.append(A(rname("cd4svd5 IDENT", rset, rmode, K), 5, 0.34, 0.28, "pair", None, None, None, rset, rmode, K))
# (ii) energy-matched random directions rendered the same way (the union's per-row contribution randomised)
for variant in ("PERM20260906", "PERM771103", "SHUF20726"):
    ARMS.append(A(f"stack6 g=0.35 {variant}", 5, 0.34, 0.35, "stack6", None, None, variant))
    for rset, rmode, K in (("all", "rank", None), ("wald89", "ratio", None), ("own", "ratio", 400)):
        if variant != "PERM20260906" and rset == "own":
            continue
        ARMS.append(A(rname("stack6 g=0.35", rset, rmode, K, variant), 5, 0.34, 0.35, "stack6", None, None, variant, rset, rmode, K))
# (iii) the same random-direction control for every R3 / R4 family on its two buildable renderers
for base, g, sel, g2, sel2 in (("stack3 g=0.35", 0.35, "stack3", None, None), ("stack3 g=0.28", 0.28, "stack3", None, None),
                               ("stack6 g=0.28", 0.28, "stack6", None, None), ("stack6 g=0.2", 0.20, "stack6", None, None),
                               ("B3 g2=0.35", 0.28, "pair", 0.35, "u2_three"), ("B3 g2=0.45", 0.28, "pair", 0.45, "u2_three"),
                               ("B4+igvf g2=0.35", 0.28, "pair", 0.35, "u2_four")):
    for rset, rmode, K in (("all", "rank", None), ("own", "ratio", 400)):
        ARMS.append(A(rname(base, rset, rmode, K, "PERM20260906"), 5, 0.34, g, sel, g2, sel2, "PERM20260906", rset, rmode, K))
NAMES = [a[0] for a in ARMS]
assert len(NAMES) == len(set(NAMES)), "duplicate arm names"


def peak_gib() -> float:
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024**3


def load_builder():
    path = ROOT / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("render274_builder", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["render274_builder"] = module
    spec.loader.exec_module(module)
    return module, path


def sources_of(MOD, names):
    return tuple((n, MOD.UNION_SOURCE_REGISTRY[n]["path"]) for n in names)


def compose(MOD, wh, cd4, wcols, panel_positions, k, resf, union, g, union2=None, g2=None):
    """Production's own path (`apply_frozen_cd4_residual`; 258/259's `compose`): CD4 deflated at k BEFORE sizing at
    f; union 1 by the carrying-row rule at g; union 2 (248's `extra2`) at g2 after it; ONE merged add, ONE preserve."""
    MOD.CD4_DEFLATE_K = k
    MOD.CD4_TARGET_F = resf
    MOD.H1_UNION_G = g
    MOD.UNION2_G = g2
    audits: dict = {}
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        final, cd4_audit, u_audit = MOD.apply_frozen_cd4_residual(
            wh, cd4, np.asarray(wcols, dtype=np.int64), np.asarray(panel_positions, dtype=np.int64),
            extra=union, extra_audit=None if union is None else {},
            extra2=union2, extra2_audit=None if union2 is None else {}, audits_out=audits)
    u2 = audits.get("union2")
    preserve = float(cd4_audit["post_add_preserve"]["post_add_global_preserve_scale"])
    info = dict(
        cd4_fraction=float(cd4_audit["fraction"]),
        union_fraction=(None if u_audit is None else float(u_audit["residual_frame"]["solved_fraction"])),
        union_carrying=(None if u_audit is None else int(u_audit["residual_frame"]["carrying_rows"])),
        union_coverage=(None if u_audit is None else float(u_audit["residual_frame"]["coverage_frac"])),
        union2_fraction=(None if u2 is None else float(u2["residual_frame"]["solved_fraction"])),
        union2_carrying=(None if u2 is None else int(u2["residual_frame"]["carrying_rows"])),
        preserve_scale=preserve,
        energy_removed=(None if "deflation" not in cd4_audit else float(cd4_audit["deflation"]["energy_removed_frac"])),
        n_carried=(None if "deflation" not in cd4_audit else int(cd4_audit["deflation"]["n_carried"])),
    )
    return final, info


# ------------------------------------------------------------------------------ THE RENDERER
def preserve(X, wh, wcols):
    """`E062.add_post_whiten_residual`'s preserve, verbatim: one scalar restoring the whitened base's median row norm
    over the operator columns."""
    base_norm = float(np.median(np.linalg.norm(wh[:, wcols], axis=1))) or 1.0
    cand_norm = float(np.median(np.linalg.norm(X[:, wcols], axis=1))) or 1.0
    s = base_norm / cand_norm
    return X * s, s


REL_TOL = 1e-9   # one-ulp rounding of `final / s` against the reference is 1e-16; a union contribution is >= 1e-4 of the row


def recover(final_b, s, X_ref):
    """The pre-preserve composition `final / s`, with every row the arm did NOT change snapped to the reference row
    (the two compositions differ there only by the rounding of two different preserve scalars, ~1e-16)."""
    X = final_b / s
    scale = np.maximum(np.abs(X_ref).max(1), 1e-300)
    changed = np.abs(X - X_ref).max(1) > REL_TOL * scale
    X[~changed] = X_ref[~changed]
    return X, changed


def support_mask(a, r, ok_cols):
    """The entries the union changed on one row, above the rounding floor."""
    return ok_cols & (np.abs(a - r) > REL_TOL * max(float(np.abs(r).max()), 1e-300))


def randomise_contribution(X_arm, X_ref, variant, ok_cols):
    """Energy-matched random directions: d = X_arm - X_ref on the changed rows, replaced by (PERM) a gene permutation
    of d over its own support, or (SHUF) another changed row's d (218's derangement, seed 20726)."""
    d = X_arm - X_ref
    changed = np.flatnonzero(np.any(d != 0, axis=1))
    for i in changed:                      # the contribution above the rounding floor only
        d[i] = np.where(support_mask(X_arm[i], X_ref[i], ok_cols), d[i], 0.0)
    kind, seed = variant[:4], int(variant[4:])
    rng = np.random.default_rng(seed)
    out = X_ref.copy()
    if kind == "PERM":
        for i in changed:
            idx = np.flatnonzero(d[i] != 0)
            vals = d[i, idx].copy()
            rng.shuffle(vals)
            out[i, idx] += vals
    elif kind == "SHUF":
        sh = changed.copy()
        while True:
            rng.shuffle(sh)
            if len(changed) < 2 or not np.any(sh == changed):
                break
        for i, j in zip(changed, sh):
            out[i] += d[j]
    else:
        raise ValueError(variant)
    return out


def render_exact_pool(X_arm, final_ref, wald_cols, changed, wh, wcols, iters=12):
    """ORACLE ceiling: hold the truth pool's L1 on every changed row at the INCUMBENT's POST-preserve value exactly
    (the ratio renderer holds it pre-preserve, and the union's off-pool magnitudes then shrink it by the preserve);
    everything off the pool stays as the union composed it.  Fixed-point iteration on the pool scales."""
    X = X_arm.copy()
    rows = np.flatnonzero(changed)
    worst = float("nan")
    for _ in range(iters):
        F, _s = preserve(X, wh, wcols)
        worst = 0.0
        for i in rows:
            idx = wald_cols[i]
            cur = float(np.abs(F[i, idx]).sum())
            tgt = float(np.abs(final_ref[i, idx]).sum())
            if cur > 0:
                c = tgt / cur
                X[i, idx] *= c
                worst = max(worst, abs(c - 1.0))
        if worst < 1e-12:
            break
    return X, dict(n_changed=int(rows.size), scale_mean=float("nan"), scale_median=float("nan"), scale_min=float("nan"),
                   scale_max=float("nan"), rank_rho_mean=float("nan"), set_size_median=89.0, zero_sign_entries=0, fixed_point_residual=worst)


def render(X_arm, X_ref, rset, rmode, K, ok_cols, own_sets, wald_cols):
    """THE RENDERER on the pre-preserve compositions.  Rows the arm did not change pass through untouched.
    Returns X_out and per-row diagnostics (scale factors, rank agreement on S, S size)."""
    X_out = X_arm.copy()
    changed = np.flatnonzero(np.any(X_arm != X_ref, axis=1))
    scales, rho, sizes, zero_sign = [], [], [], 0
    for i in changed:
        a, r = X_arm[i], X_ref[i]
        if rset == "all":
            idx = np.flatnonzero(ok_cols)
        elif rset == "supp":
            idx = np.flatnonzero(support_mask(a, r, ok_cols))
        elif rset == "own":
            idx = own_sets[K][i]
        elif rset == "wald89":
            idx = wald_cols[i]
        else:
            raise ValueError(rset)
        aS, rS = a[idx], r[idx]
        if rmode == "ratio":
            den = float(np.abs(aS).sum())
            c = float(np.abs(rS).sum()) / den if den > 0 else 1.0
            X_out[i, idx] = aS * c
            scales.append(c)
        elif rmode == "norm":
            den = float(np.linalg.norm(aS))
            c = float(np.linalg.norm(rS)) / den if den > 0 else 1.0
            X_out[i, idx] = aS * c
            scales.append(c)
        elif rmode == "rank":
            order = np.argsort(np.abs(aS), kind="stable")
            vals = np.sort(np.abs(rS))
            mag = np.empty_like(vals)
            mag[order] = vals
            sg = np.sign(aS)
            zero_sign += int(np.count_nonzero((sg == 0) & (mag != 0)))
            X_out[i, idx] = sg * mag
            scales.append(float(np.linalg.norm(X_out[i, idx]) / max(np.linalg.norm(aS), 1e-300)))
        else:
            raise ValueError(rmode)
        if idx.size > 2:
            ra = np.argsort(np.argsort(np.abs(aS)))
            rr = np.argsort(np.argsort(np.abs(rS)))
            rho.append(float(np.corrcoef(ra, rr)[0, 1]))
        sizes.append(int(idx.size))
    diag = dict(n_changed=int(changed.size), scale_mean=float(np.mean(scales)) if scales else 1.0,
                scale_median=float(np.median(scales)) if scales else 1.0,
                scale_min=float(np.min(scales)) if scales else 1.0, scale_max=float(np.max(scales)) if scales else 1.0,
                rank_rho_mean=float(np.mean(rho)) if rho else float("nan"), set_size_median=float(np.median(sizes)) if sizes else 0.0,
                zero_sign_entries=zero_sign)
    return X_out, diag, changed


# ------------------------------------------------------------------------------ fold phase
def challenge_contexts():
    out = {}
    for c in CTX:
        b = plolo.challenge_bctrl(c)
        e = np.expm1(b)
        gate = (1e6 * e / float(e.sum()) > plolo.MIN_CPM) & (b > 0)
        out[c] = (b, gate)
    return out


def _per_target(P, T):
    A_ = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    B_ = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    S = A_ @ B_.T
    n = S.shape[0]
    own = np.diag(S).copy()
    pds = ((S < own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1)) / (n - 1)
    return pds, own


def amp_terms(dl_rows, targets, gpos):
    """243's A_c on an effect-level render for one context (the energy-weighted delivered amplitude, target column
    removed), exactly `terms_from_final`'s `A` -- computed on the fold's targets here, on the production panel in --prod."""
    cols = {j: gpos[t] for j, t in enumerate(targets) if t in gpos}
    w = _row_weights(len(targets), N_GENES, cols)
    a = np.einsum("ij,ij->i", dl_rows, dl_rows)
    for j, col in cols.items():
        a[j] -= dl_rows[j, col] ** 2
    return float((a * w).sum())


def score_arm(final, fold, Ts, ctx, LR, pools, gpos):
    """258's two functionals from six renders (2 tags x 3 frames) + the fold-frame delivered amplitude A_c."""
    sc = fold["scored"]
    out: dict[str, object] = {}
    per_pds, per_cos = [], []
    H_A = None
    for tag, mass in (("a23", "minnorm"), ("a15", None)):
        acc = {}
        for c in CTX:
            b, gate = ctx[c]
            dl = deliver(final, b, gate, mass)
            P = dl[:, sc]
            e = endpoints(P, Ts)
            e["cos_agg"] = float(np.sqrt(max(0.0, 1.0 - e["sse"])))
            acc[c] = e
            if tag == "a23":
                p, own = _per_target(P, Ts)
                per_pds.append(p)
                per_cos.append(own)
                out[f"ampA_{c}"] = amp_terms(dl, fold["targets"], gpos)
                if c == "A":
                    out["pds_A"] = float(np.mean(p))
            elif c == "A":
                H_A = dl
        for k in acc["A"]:
            out[f"{tag}_{k}"] = float(np.mean([acc[c][k] for c in CTX]))
    out.update({k: out[f"a23_{k}"] for k in ("pds", "cos", "cos_agg", "sgn", "sgn_w", "sse")})
    out["median_row_norm_final"] = float(np.median(np.linalg.norm(final, axis=1)))
    bA, gA = ctx["A"]
    LP = to_lfc(H_A, bA, gA)[:, sc]
    umd_out = {pk: umd(LP, LR, mask) for pk, mask in pools.items()}
    out["top89_precision"] = float(key_precision(np.abs(H_A[:, sc]), H_A[:, sc], Ts, n=POOL_N))
    out["delivered_norm_A"] = float(np.median(np.linalg.norm(H_A[:, fold["wcols"]], axis=1)))
    arrays = dict(pds_pt=np.mean(per_pds, axis=0), cos_pt=np.mean(per_cos, axis=0))
    for pk, (U, M, D) in umd_out.items():
        arrays[f"U_{pk}"], arrays[f"M_{pk}"], arrays[f"D_{pk}"] = U, M, D
    # effect-level L1 on the wald89 pool per row (the magnitude the renderer controls, before the render)
    arrays["poolL1_wald89"] = np.array([float(np.abs(final[i, sc][pools["wald89"][i]]).sum()) for i in range(final.shape[0])])
    arrays["rownorm_final"] = np.linalg.norm(final, axis=1)
    return out, arrays


def run_fold(held: str) -> None:
    t0 = time.time()
    MOD, builder_path = load_builder()
    assert MOD.CD4_DEFLATE_K == 5 and MOD.H1_UNION_G == 0.28 and MOD.CD4_TARGET_F == 0.34
    for name in ("CD4_SIZING", "UNION2_G", "CD4_DEFLATE_BASIS", "CD4_PARTIAL_G"):
        assert getattr(MOD, name) is None, f"builder global {name} is live"
    print(f"builder {builder_path.relative_to(ROOT)} sha256 {MOD.file_sha256(builder_path)[:16]}", flush=True)

    fold = build_fold(held)
    targets = fold["targets"]
    truth_h = np.asarray(fold["truth"], dtype=np.float64)
    pre = np.asarray(fold["pre"], dtype=np.float64)
    cd4 = cd4_for(targets)
    wcols = np.asarray(fold["wcols"], dtype=np.int64)
    panel_positions = np.asarray(fold["panel_positions"], dtype=np.int64)
    needed = sorted({a[4] for a in ARMS if a[4] is not None} | {a[6] for a in ARMS if a[6] is not None})
    unions = {}
    for sel in needed:
        u, _ = MOD.prepare_h1_union(targets, OFFICIAL, PANEL, max_rss_gib=40.0, sources=sources_of(MOD, SEL[sel]))
        unions[sel] = u
        print(f"[{held}] union {sel:10s} {int(np.sum(np.asarray(u.target_condition_count) > 0)):3d}/{len(targets)} carried ({time.time()-t0:.0f}s)", flush=True)
    op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production" else PROD_OPS / "whitener_union_pruned.npz")
    op = Operator(op_path)
    wh = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    del op
    gc.collect()
    ctx = challenge_contexts()
    bA, gA = ctx["A"]
    sc = np.asarray(fold["scored"])
    Ts = truth_h[:, sc]
    cpm = 1e6 * np.expm1(fold["b_ctrl"][sc]) / float(np.expm1(fold["b_ctrl"]).sum())
    LR = to_lfc(truth_h, fold["b_ctrl"], gA)[:, sc]
    pools = {
        "wald89": topn(np.abs(LR) * np.sqrt(np.maximum(cpm, 0))[None, :], POOL_N),
        "lfc89": topn(np.abs(LR), POOL_N),
        "h89": topn(np.abs(truth_h[:, sc]), POOL_N),
    }
    gpos = {g: i for i, g in enumerate(OFFICIAL)}
    print(f"[{held}] staged {len(targets)} targets ({time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB)", flush=True)

    # ---- the renderer's fixed inputs: the incumbent's pre-preserve composition and its magnitude sets ----
    ok_cols = np.ones(N_GENES, dtype=bool)
    ok_cols[panel_positions] = False
    final_inc, info_inc = compose(MOD, wh, cd4, wcols, panel_positions, 5, 0.34, unions["pair"], 0.28)
    X_ref = final_inc / info_inc["preserve_scale"]
    Xr, s_r = preserve(X_ref, wh, wcols)
    print(f"[{held}] incumbent: preserve {info_inc['preserve_scale']:.12f}; re-preserved {s_r:.12f}; max|d| {np.abs(Xr - final_inc).max():.2e}", flush=True)
    own_sets = {}
    absr = np.where(ok_cols[None, :], np.abs(X_ref), -1.0)
    for K in (89, 400, 2000):
        own_sets[K] = [np.argpartition(-absr[i], K - 1)[:K] for i in range(X_ref.shape[0])]
    wald_cols = [sc[np.flatnonzero(pools["wald89"][i])] for i in range(X_ref.shape[0])]

    rows: dict[str, dict] = {}
    arrays: dict[str, np.ndarray] = {}
    finals_cache: dict[tuple, tuple] = {}
    for name, k, resf, g, sel, g2, sel2, variant, rset, rmode, K in ARMS:
        t1 = time.time()
        key = (k, resf, g, sel, g2, sel2)
        if key not in finals_cache:
            finals_cache[key] = compose(MOD, wh, cd4, wcols, panel_positions, k, resf, unions.get(sel), g,
                                        unions.get(sel2) if sel2 else None, g2)
        final_b, info = finals_cache[key]
        rdiag = None
        X_arm, changed_mask = recover(final_b, info["preserve_scale"], X_ref)
        if variant is None and rset is None:
            final = final_b.copy()
        else:
            if variant is not None:
                X_arm = randomise_contribution(X_arm, X_ref, variant, ok_cols)
            if rset is not None and rmode == "exact":
                X_arm, rdiag = render_exact_pool(X_arm, final_inc, wald_cols, changed_mask, wh, wcols)
            elif rset is not None:
                X_arm, rdiag, _ = render(X_arm, X_ref, rset, rmode, K, ok_cols, own_sets, wald_cols)
            final, s_out = preserve(X_arm, wh, wcols)
            info = {**info, "preserve_scale": float(s_out)}
        final[:, panel_positions] = 0.0
        s, arr = score_arm(final, fold, Ts, ctx, LR, pools, gpos)
        rows[name] = dict(arm=name, k=k, f=resf, g=g, selection=sel, g2=g2, selection2=sel2, variant=variant,
                          rset=rset, rmode=rmode, K=K, **info, **s, renderer=rdiag, seconds=round(time.time() - t1, 1))
        for kk, v in arr.items():
            arrays[f"{name}|{kk}"] = v
        arrays[f"{name}|changed"] = changed_mask
        rd = "" if rdiag is None else f" | rows {rdiag['n_changed']} scale med {rdiag['scale_median']:.4f} rho {rdiag['rank_rho_mean']:.3f}"
        print(f"[{held}] {name:44s} agg {s['cos_agg']:.6f} mac {s['a23_cos']:.6f} pds {s['pds']:.6f} pdsA {s['pds_A']:.6f} "
              f"prec {s['top89_precision']:.5f} pres {info['preserve_scale']:.6f}{rd} ({rows[name]['seconds']}s)", flush=True)
        del final
        gc.collect()

    # ---------------------------------------------------------------- known-answer rows, this fold
    ka = []

    def check(label, got, want, tol=1e-9):
        d = abs(float(got) - float(want))
        ok = d <= tol
        ka.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, ok=bool(ok)))
        print(f"  {'PASS' if ok else 'FAIL'} {label:84s} got {got:.10f} want {want:.10f} |d| {d:.1e}", flush=True)

    print(f"\n[{held}] KNOWN-ANSWER ROWS", flush=True)
    # (a) 258's rows, bit-identical (scalars and per-target arrays)
    f258 = json.loads((E258 / f"work/fold-{held}.json").read_text())["rows"]
    a258 = np.load(E258 / f"work/fold-{held}.npz")
    for n in ("a15 f=0.10 SHIPPED", "cd4 f=0.34", OLD, ARM6, INC, "stack6 g=0.35", "stack6 g=0.28", "stack6 g=0.2", "stack3 g=0.35", "stack3 g=0.28"):
        for key in ("cos_agg", "a23_cos", "pds", "pds_A", "top89_precision", "delivered_norm_A", "preserve_scale", "cd4_fraction"):
            check(f"258 {n}.{key}", rows[n][key], f258[n][key], 0.0)
        for kk in ("pds_pt", "cos_pt", "U_wald89", "M_wald89", "D_wald89"):
            check(f"258 {n}|{kk} max|d|", float(np.abs(arrays[f"{n}|{kk}"] - a258[f"{n}|{kk}"]).max()), 0.0, 0.0)
    f259 = json.loads((E259 / f"work/fold-{held}.json").read_text())["rows"]
    a259 = np.load(E259 / f"work/fold-{held}.npz")
    for n in ("B3 g2=0.35", "B3 g2=0.45", "B4+igvf g2=0.35"):
        for key in ("cos_agg", "a23_cos", "pds", "pds_A", "top89_precision", "preserve_scale", "union2_fraction"):
            check(f"259 {n}.{key}", rows[n][key], f259[n][key], 0.0)
        for kk in ("pds_pt", "U_wald89", "M_wald89"):
            check(f"259 {n}|{kk} max|d|", float(np.abs(arrays[f"{n}|{kk}"] - a259[f"{n}|{kk}"]).max()), 0.0, 0.0)
    # (b) 246 / 249 / 229 rows as 258 checked them
    s246a = {r["arm"]: r for r in json.loads((E246 / f"work/sweep246-{held}.json").read_text())["rows"]}
    s246b = {r["arm"]: r for r in json.loads((E246 / f"work/sweep246b-{held}.json").read_text())["rows"]}
    for key in ("cos_agg", "a23_cos", "pds"):
        check(f"246 sweep246 STACK.{key} == h1union", rows[OLD][key], s246a["STACK"][key])
        check(f"246 sweep246b B2_cd4_svd5_resolved.{key} == cd4svd5", rows[INC][key], s246b["B2_cd4_svd5_resolved"][key])
    r249 = json.loads((E249 / f"work/repro249-{held}.json").read_text())
    for key in ("cos_agg", "a23_cos", "pds"):
        check(f"249 repro249 arm.{key} == cd4svd5", rows[INC][key], r249["arm"][key])
    six229 = json.loads((E229 / "work/six229.json").read_text())["arms"]
    check("229 six229 stack6 g=0.35 <- ARM pds_per_fold (frame A) == ARM6", rows[ARM6]["pds_A"], six229["stack6 g=0.35 <- ARM"]["pds_per_fold"][held])
    check("229 six229 h1union pds_per_fold (frame A) == h1union", rows[OLD]["pds_A"], six229["h1union g=0.28 <- INCUMBENT"]["pds_per_fold"][held])
    # (c) THE RENDERER's own known answers
    check("renderer: re-preserved incumbent composition == builder's final (max|d|)", float(np.abs(Xr - final_inc).max()), 0.0, 1e-12)
    fb, ib = finals_cache[(5, 0.34, 0.35, "stack6", None, None)]
    X6, ch6 = recover(fb, ib["preserve_scale"], X_ref)
    check("recover(): changed rows == the six-source union's carried rows", int(ch6.sum()), int(np.sum(np.asarray(unions["stack6"].target_condition_count) > 0)), 0.0)
    check("recover(): max relative |final/s - X_ref| on the UNchanged rows (rounding floor)", float(np.abs((fb / ib["preserve_scale"])[~ch6] - X_ref[~ch6]).max() / np.abs(X_ref[~ch6]).max()) if (~ch6).any() else 0.0, 0.0, 1e-13)
    Xa, _, _ = render(X6, X_ref, "all", "rank", None, ok_cols, own_sets, wald_cols)
    # the transplant: every changed row's sorted magnitudes == the incumbent's, its norm == the incumbent's
    ch = np.flatnonzero(ch6)
    md = max(float(np.abs(np.sort(np.abs(Xa[i, ok_cols])) - np.sort(np.abs(X_ref[i, ok_cols]))).max()) for i in ch)
    check("renderer all+rank: sorted |X_out| == sorted |X_ref| on every changed row (max|d|)", md, 0.0, 1e-12)
    md = max(abs(np.linalg.norm(Xa[i, ok_cols]) / np.linalg.norm(X_ref[i, ok_cols]) - 1) for i in ch)
    check("renderer all+rank: row norm ratio to the incumbent - 1 (max)", md, 0.0, 1e-12)
    Xw, _, _ = render(X6, X_ref, "wald89", "ratio", None, ok_cols, own_sets, wald_cols)
    md = max(abs(np.abs(Xw[i, wald_cols[i]]).sum() / np.abs(X_ref[i, wald_cols[i]]).sum() - 1) for i in ch)
    check("renderer wald89+ratio: pool L1 ratio to the incumbent - 1 (max over changed rows)", md, 0.0, 1e-12)
    # untouched rows pass through bit for bit (== the reference rows exactly; == the arm's own rows to rounding)
    unt = np.flatnonzero(~ch6)
    check("renderer: untouched rows == the reference rows (max|d|)", float(np.abs(Xa[unt] - X_ref[unt]).max()) if unt.size else 0.0, 0.0, 0.0)
    check("renderer: untouched rows == arm's own rows (max|d|, rounding)", float(np.abs(Xa[unt] - (fb / ib["preserve_scale"])[unt]).max()) if unt.size else 0.0, 0.0, 1e-12)
    # the rendered stack6 g=0.35 | all rank cell's CHANGED entries are exactly the union's carried rows
    ns = np.array([int(support_mask(X6[i], X_ref[i], ok_cols).sum()) for i in ch])
    check("supp set: every changed row has > 0 entries above the rounding floor (min count)", float(ns.min() > 0), 1.0, 0.0)
    check("supp set: no changed row has every column above the floor (max count / N_GENES < 0.99)", float(ns.max() / N_GENES < 0.99), 1.0, 0.0)
    print(f"  supp set sizes over the changed rows: min {ns.min()} median {int(np.median(ns))} max {ns.max()} of {N_GENES}", flush=True)
    # the identity check: every IDENT cell == the incumbent at 0.0 on every scalar and array
    for n in [a[0] for a in ARMS if a[0].startswith("cd4svd5 IDENT")]:
        for key in ("cos_agg", "pds", "pds_A", "top89_precision", "preserve_scale"):
            check(f"IDENT {n}.{key} == cd4svd5", rows[n][key], rows[INC][key], 1e-12)
        check(f"IDENT {n}|pds_pt max|d| (rank statistic, exact)", float(np.abs(arrays[f"{n}|pds_pt"] - arrays[f"{INC}|pds_pt"]).max()), 0.0, 0.0)
        check(f"IDENT {n}|U_wald89 max|d|", float(np.abs(arrays[f"{n}|U_wald89"] - arrays[f"{INC}|U_wald89"]).max()), 0.0, 1e-12)
        check(f"IDENT {n}: no row changed", int(arrays[f"{n}|changed"].sum()), 0, 0.0)
    # the identity renderer on the arm (rset None) is the builder's final: re-preserve of final/s reproduces final
    fr, sr = preserve(X6, wh, wcols)
    check("renderer at identity on stack6 g=0.35: re-preserve(recover(final)) == builder's final (max|d|)", float(np.abs(fr - fb).max()), 0.0, 1e-12)
    Xe, de = render_exact_pool(X6, final_inc, wald_cols, ch6, wh, wcols)
    Fe, _ = preserve(Xe, wh, wcols)
    md = max(abs(np.abs(Fe[i, wald_cols[i]]).sum() / np.abs(final_inc[i, wald_cols[i]]).sum() - 1) for i in ch)
    check("renderer wald89+exact: POST-preserve pool L1 ratio to the incumbent - 1 (max over changed rows)", md, 0.0, 1e-9)
    check("renderer wald89+exact: off-pool entries untouched on changed rows (max|d| vs X6 off the pool)",
          max(float(np.abs(np.delete(Xe[i], wald_cols[i]) - np.delete(X6[i], wald_cols[i])).max()) for i in ch), 0.0, 0.0)
    n_ok = sum(r["ok"] for r in ka)
    print(f"[{held}] known-answer rows {n_ok}/{len(ka)} pass", flush=True)

    # 262's per-target base cosine to the truth, aligned by target (for the group tables)
    z262 = np.load(E262 / f"work/diag262-{held}.npz", allow_pickle=True)
    t262 = {str(t): i for i, t in enumerate(z262["targets"])}
    arrays["cov|base_cos_truth"] = np.array([z262["cov|base_cos_truth"][t262[t]] if t in t262 else np.nan for t in targets])
    out = dict(fold=held, n_targets=len(targets), targets=list(targets), builder_sha256=MOD.file_sha256(builder_path),
               rows=rows, known_answer=ka, union_carried={sel: int(np.sum(np.asarray(u.target_condition_count) > 0)) for sel, u in unions.items()},
               seconds=round(time.time() - t0, 1), peak_rss_gib=peak_gib())
    (WORK / f"fold-{held}.json").write_text(json.dumps(out, indent=1, default=float))
    np.savez_compressed(WORK / f"fold-{held}.npz", **arrays)
    print(f"[{held}] {len(rows)} arms in {out['seconds']}s, peak RSS {out['peak_rss_gib']:.2f} GiB; wrote {WORK/f'fold-{held}.json'}", flush=True)
    if n_ok != len(ka):
        raise SystemExit(1)


# ------------------------------------------------------------------------------ production phase
def terms_from_final(final, targets, b, gate):
    gpos = {g: i for i, g in enumerate(OFFICIAL)}
    per = {}
    for c in CTX:
        dh = deliver(final, b[c], gate[c], "minnorm")
        mu = dh + b[c][None, :]
        cols = {j: gpos[t] for j, t in enumerate(targets) if t in gpos}
        w = _row_weights(len(targets), N_GENES, cols)
        a = np.einsum("ij,ij->i", dh, dh)
        for j, col in cols.items():
            a[j] -= dh[j, col] ** 2
        per[c] = dict(A=float((a * w).sum()), B=float(_across_pert_budget(mu, cols, w)), Kc=KC[c],
                      cap=float(w.sum()) * CP46[c], n=len(targets), resolved=len(cols), wsum=float(w.sum()))
    return per


def amp_ratio(per_arm, per_base):
    return float(np.mean([np.sqrt(per_arm[c]["A"] / per_base[c]["A"]) for c in CTX]))


def run_prod() -> None:
    t0 = time.time()
    MOD, builder_path = load_builder()
    stage = E249 / "work/stage"
    z = np.load(stage / "cd4-off.npz", allow_pickle=True)
    targets = tuple(str(t) for t in z["targets"])
    cd4 = E062.CD4Prepared(mode="complete_mean3", mean_effect=np.asarray(z["mean_effect"], dtype=np.float64),
                           gene_votes=np.asarray(z["gene_votes"], dtype=np.uint8),
                           target_condition_count=np.asarray(z["target_condition_count"], dtype=np.uint8), condition_details=[])
    u = np.load(stage / "union-off.npz", allow_pickle=True)
    union_pair_saved = E062.CD4Prepared(mode="union", mean_effect=np.asarray(u["mean_effect"], dtype=np.float64),
                                        gene_votes=np.asarray(u["gene_votes"], dtype=np.uint8),
                                        target_condition_count=np.asarray(u["target_condition_count"], dtype=np.uint8), condition_details=[])
    final_off = np.load(stage / "effects-off.npy")
    aud_off = json.loads((stage / "stage-off.json").read_text())
    aud_svd5 = json.loads((stage / "stage-svd5.json").read_text())
    s_off = float(aud_off["audit"]["h1_union"]["merged_preserve"]["post_add_global_preserve_scale"])
    fc = float(aud_off["audit"]["cd4"]["fraction"])
    fu = float(aud_off["audit"]["h1_union"]["residual_frame"]["solved_fraction"])
    base_norm = float(aud_off["audit"]["h1_union"]["residual_frame"]["whitened_base_median_row_norm"])
    zk = np.load(WORK / "plolo/fold-kolf2.npz", allow_pickle=True)
    wcols = np.asarray(zk["wcols"], dtype=np.int64)
    axis = axis_sha256(OFFICIAL[i] for i in wcols)
    assert axis == plolo.OPERATOR_REGISTRY["union_pruned"]["axis_sha256"], axis
    gpos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([gpos[t] for t in targets], dtype=np.int64)
    vc = np.asarray(cd4.gene_votes) > 0
    vu = np.asarray(union_pair_saved.gene_votes) > 0
    merged = np.where(vc, np.asarray(cd4.mean_effect), 0.0) * fc + np.where(vu, np.asarray(union_pair_saved.mean_effect), 0.0) * fu
    wh = final_off / s_off - merged
    del merged
    ka = []

    def check(label, got, want, tol):
        d = abs(float(got) - float(want))
        ok = d <= tol
        ka.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, ok=bool(ok)))
        print(f"  {'PASS' if ok else 'FAIL'} {label:84s} got {got:.12g} want {want:.12g} |d| {d:.1e}", flush=True)

    print("\nKNOWN-ANSWER ROWS (production frame)", flush=True)
    check("recovered base: median row norm == stage-off whitened_base_median_row_norm", float(np.median(np.linalg.norm(wh, axis=1))), base_norm, 1e-9)
    unions = {sel: MOD.prepare_h1_union(targets, OFFICIAL, PANEL, max_rss_gib=40.0, sources=sources_of(MOD, SEL[sel]))[0] for sel in SEL}
    check("prepare_h1_union(PAIR) == 249's union-off.npz (max|d|)", float(np.abs(np.asarray(unions["pair"].mean_effect) - np.asarray(union_pair_saved.mean_effect)).max()), 0.0, 0.0)
    f_off, i_off = compose(MOD, wh, cd4, wcols, panel_positions, None, 0.34, unions["pair"], 0.28)
    check("re-composed h1union == 249 effects-off.npy (abeb0749), max|d|", float(np.abs(f_off - final_off).max()), 0.0, 1e-9)
    final_svd5 = np.load(stage / "effects-svd5.npy")
    f_5, i_5 = compose(MOD, wh, cd4, wcols, panel_positions, 5, 0.34, unions["pair"], 0.28)
    check("re-composed cd4svd5 == 249 effects-svd5.npy (a28d5be2), max|d|", float(np.abs(f_5 - final_svd5).max()), 0.0, 1e-9)
    check("re-composed cd4svd5 preserve == stage-svd5", i_5["preserve_scale"], float(aud_svd5["audit"]["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"]), 1e-12)
    check("effects-svd5.npy sha == a28d5be2...", 1.0 if MOD.array_sha256(final_svd5).startswith("a28d5be2") else 0.0, 1.0, 0.0)
    final_arm6 = np.load(E229 / "work/effects-stack6g035.npy")
    f_6, i_6 = compose(MOD, wh, cd4, wcols, panel_positions, None, 0.34, unions["stack6"], 0.35)
    check("re-composed ARM6 == 229 effects-stack6g035.npy (c5721825), max|d|", float(np.abs(f_6 - final_arm6).max()), 0.0, 1e-9)
    check("effects-stack6g035.npy sha == c5721825...", 1.0 if MOD.array_sha256(final_arm6).startswith("c5721825") else 0.0, 1.0, 0.0)
    final_svd10 = np.load(E255 / "work/stage/effects-svd10.npy")
    # the renderer on the production frame: X_ref from the incumbent
    ok_cols = np.ones(N_GENES, dtype=bool)
    ok_cols[panel_positions] = False
    X_ref = f_5 / i_5["preserve_scale"]
    Xr, s_r = preserve(X_ref, wh, wcols)
    check("renderer: re-preserved incumbent == effects-svd5.npy (max|d|)", float(np.abs(Xr - final_svd5).max()), 0.0, 1e-12)
    absr = np.where(ok_cols[None, :], np.abs(X_ref), -1.0)
    own_sets = {K: [np.argpartition(-absr[i], K - 1)[:K] for i in range(X_ref.shape[0])] for K in (89, 400, 2000)}
    n_ok = sum(r["ok"] for r in ka)
    print(f"known-answer rows {n_ok}/{len(ka)} pass  ({time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB)", flush=True)
    if n_ok != len(ka):
        (WORK / "prod274.json").write_text(json.dumps(dict(known_answer=ka, failed=True), indent=1))
        raise SystemExit(1)

    ctx = challenge_contexts()
    b = {c: ctx[c][0] for c in CTX}
    gate = {c: ctx[c][1] for c in CTX}
    bA, gA = ctx["A"]
    terms: dict[str, dict] = {}
    info_all: dict[str, dict] = {}
    rdiags: dict[str, dict] = {}
    finals = {"h1union g=0.28 (file)": final_off, "cd4svd5 (file)": final_svd5, "ARM6 (file)": final_arm6, "cd4svd10 (255 file)": final_svd10}
    cache: dict[tuple, tuple] = {}
    for name, k, resf, g, sel, g2, sel2, variant, rset, rmode, K in ARMS:
        if name in ("a15 f=0.10 SHIPPED", "cd4 f=0.34"):
            continue
        if rset == "wald89":
            continue                      # ORACLE set: not composable without the truth
        key = (k, resf, g, sel, g2, sel2)
        if key not in cache:
            cache[key] = compose(MOD, wh, cd4, wcols, panel_positions, k, resf, unions.get(sel), g, unions.get(sel2) if sel2 else None, g2)
        fb, info = cache[key]
        if variant is None and rset is None:
            final = fb.copy()
        else:
            X_arm, _ = recover(fb, info["preserve_scale"], X_ref)
            if variant is not None:
                X_arm = randomise_contribution(X_arm, X_ref, variant, ok_cols)
            if rset is not None:
                X_arm, rd, _ = render(X_arm, X_ref, rset, rmode, K, ok_cols, own_sets, None)
                rdiags[name] = rd
            final, s_out = preserve(X_arm, wh, wcols)
            info = {**info, "preserve_scale": float(s_out)}
        finals[name] = final
        info_all[name] = info
    for name, final in finals.items():
        t1 = time.time()
        per = terms_from_final(final, targets, b, gate)
        H = deliver(final, bA, gA, None)
        terms[name] = dict(per=per, delivered_norm_A=float(np.median(np.linalg.norm(H[:, wcols], axis=1))),
                           median_row_norm_final=float(np.median(np.linalg.norm(final, axis=1))))
        print(f"  {name:46s} A {per['A']['A']:.6f}/{per['B']['A']:.6f}/{per['C']['A']:.6f}  dlvA {terms[name]['delivered_norm_A']:.5f} ({time.time()-t1:.0f}s)", flush=True)
        del H
    inv243 = json.loads((E243 / "work/invert243.json").read_text())
    amp249 = json.loads((E249 / "work/amp249.json").read_text())
    amp255 = json.loads((E255 / "work/amp255.json").read_text())
    arm6_243 = json.loads((E243 / "work/arm6_243.json").read_text())
    pkg = dict(h1union=inv243["arms"]["h1union"]["per"], cd4svd5=amp249["per_arm"], ARM6=arm6_243["per"], cd4svd10=amp255["per_arm"])
    calib = []
    for arm, base, fa, fbn in (("ARM6", "h1union", "ARM6 (file)", "h1union g=0.28 (file)"), ("cd4svd5", "h1union", "cd4svd5 (file)", "h1union g=0.28 (file)"),
                               ("cd4svd10", "cd4svd5", "cd4svd10 (255 file)", "cd4svd5 (file)")):
        r_pkg = amp_ratio(pkg[arm], pkg[base])
        r_eff = amp_ratio(terms[fa]["per"], terms[fbn]["per"])
        calib.append(dict(arm=arm, base=base, package_pct=100 * (r_pkg - 1), effect_level_pct=100 * (r_eff - 1),
                          ratio_pkg_over_eff=(r_pkg - 1) / (r_eff - 1) if r_eff != 1 else None))
        print(f"  amplitude calibration {arm:8s} vs {base:8s}: package {100*(r_pkg-1):+.4f}%  effect-level {100*(r_eff-1):+.4f}%", flush=True)
    out = dict(known_answer=ka, calibration=calib, terms=terms, compose_info=info_all, renderer=rdiags, package_terms=pkg,
               seconds=round(time.time() - t0, 1), peak_rss_gib=peak_gib())
    (WORK / "prod274.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"wrote {WORK/'prod274.json'}  ({out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB)", flush=True)


# ------------------------------------------------------------------------------ price phase
def raw_for(per, cos):
    num = float(np.mean([per[c]["A"] + per[c]["Kc"] - min(per[c]["cap"], per[c]["B"]) for c in CTX]))
    msa = float(np.mean([np.sqrt(per[c]["A"]) for c in CTX]))
    return 1.0 + (num - 2.0 * cos * np.sqrt(SUMD) * msa) / SUMD


def invert(per, raw):
    num = float(np.mean([per[c]["A"] + per[c]["Kc"] - min(per[c]["cap"], per[c]["B"]) for c in CTX]))
    X = (num - (raw - 1.0) * SUMD) / 2.0
    msa = float(np.mean([np.sqrt(per[c]["A"]) for c in CTX]))
    return X / (np.sqrt(SUMD) * msa)


def scale_terms(per, r):
    return {c: {**per[c], "A": per[c]["A"] * r * r, "B": per[c]["B"] * r * r} for c in CTX}


def board_entry(eid):
    e = json.loads((ROOT / f"data/board-entries/{eid}.json").read_text())["entry"]
    return dict(overall=e["score_avg"],
                raw=dict(pds=e["pds_cosine"], mse=e["expr_mse_unbiased_capped_norm"], nmae=e["de_wilcoxon_lfc_nmae"],
                         fid=e["de_wilcoxon_direction_fidelity_yield_raw"], reach=e["de_wilcoxon_direction_reach_raw"], jac=e["de_wilcoxon_sig_jaccard"]),
                scaled={m: e[f"score_{m}"] for m in ("pds", "mse", "nmae", "fid", "reach", "jac")})


def legs(rel_agg, d_pds, d_nmae, r_amp, base_per, cos_base, raw_base, k_cos, k_pds, de_rate, nmae_x):
    """258/259's six-member legs of one cell on one factor set."""
    dcos = k_cos * rel_agg * cos_base
    per_arm = scale_terms(base_per, r_amp)
    raw = raw_for(per_arm, cos_base + dcos)
    raw_amp_only = raw_for(per_arm, cos_base)
    o_mse = -(raw - raw_base) / SPAN_MSE / 6
    o_amp = -(raw_amp_only - raw_base) / SPAN_MSE / 6
    o_pds = PDS_SCALE * k_pds * d_pds / 6
    o_de3 = de_rate * dcos
    o_nmae = -d_nmae / NMAE_SPAN / 6 * nmae_x
    de = {m: o_de3 * DE3_SPLIT[m] for m in DE3_SPLIT}
    return dict(k_cos=k_cos, k_pds=k_pds, de3_rate=de_rate, nmae_x=nmae_x, dcos=dcos, raw=raw, o_mse=o_mse, o_amp=o_amp,
                o_mse_cos_only=o_mse - o_amp, o_pds=o_pds, o_de3=o_de3, o_nmae=o_nmae, **{f"o_{m}": v for m, v in de.items()},
                total=o_mse + o_pds + o_de3 + o_nmae)


def scenario_factors(which, rel_agg, d_pds, d_nmae):
    """COHERENT scenarios (259 / READ-FIRST s3): per leg, the end of its band chosen by the SIGN of that leg for this
    cell.  Union legs carry the new-source bands: k_cos 0.94 / 0.96 (/ 1.00 optimistic), k_pds 0.7 / 2.48, DE-three
    0.0425 / 0.24, the identity x1.5 on a loss and x0 on a gain (pessimistic), x1.0 (point), x1.0 loss / x1.5 gain (optimistic)."""
    gain_cos, gain_pds, gain_nmae = rel_agg > 0, d_pds > 0, d_nmae < 0
    if which == "POINT":
        return dict(POINT)
    if which == "PESSIMISTIC":
        return dict(k_cos=K_COS_BAND[0] if gain_cos else K_COS_OPT, k_pds=K_PDS_BAND[0] if gain_pds else K_PDS_BAND[1],
                    de3=DE3_BAND[0] if gain_cos else DE3_BAND[1], nmae_x=0.0 if gain_nmae else NMAE_LOSS_X)
    if which == "OPTIMISTIC":
        return dict(k_cos=K_COS_OPT if gain_cos else K_COS_BAND[0], k_pds=K_PDS_BAND[1] if gain_pds else K_PDS_BAND[0],
                    de3=DE3_BAND[1] if gain_cos else DE3_BAND[0], nmae_x=NMAE_LOSS_X if gain_nmae else 1.0)
    raise ValueError(which)


def boot(diff, rng, B=BOOT_B):
    n = diff.size
    idx = rng.integers(0, n, size=(B, n))
    means = diff[idx].mean(1)
    lo, hi = np.percentile(means, [2.5, 97.5])
    return dict(point=float(diff.mean()), lo=float(lo), hi=float(hi), sd=float(means.std(ddof=1)), frac_le0=float((means <= 0).mean()), n=int(n))


def boot_folds(diffs, rng, B=BOOT_B):
    acc = np.zeros(B)
    for d in diffs:
        idx = rng.integers(0, d.size, size=(B, d.size))
        acc += d[idx].mean(1)
    means = acc / len(diffs)
    lo, hi = np.percentile(means, [2.5, 97.5])
    return dict(point=float(np.mean([d.mean() for d in diffs])), lo=float(lo), hi=float(hi), sd=float(means.std(ddof=1)),
                frac_le0=float((means <= 0).mean()), per_fold=[float(d.mean()) for d in diffs])


def _retention(arrs, name, base, kk, d_mean, targets=None):
    loo, share, top = [], [], {}
    for h in CANCER:
        d = arrs[h][f"{name}|pds_pt"] - arrs[h][f"{base}|pds_pt"]
        order = np.argsort(-d)
        keep = np.ones(d.size, dtype=bool)
        keep[order[:kk]] = False
        loo.append(float(d[keep].mean()))
        share.append(float(d[order[:kk]].sum() / d.sum()) if d.sum() != 0 else float("nan"))
        if targets is not None:
            top[h] = [f"{targets[h][i]} {d[i]:+.3f}" for i in order[:kk]]
    return dict(loo_mean=float(np.mean(loo)), retained=(float(np.mean(loo) / d_mean) if d_mean != 0 else float("nan")),
                share=float(np.nanmean(share)), **({"top_names": top} if targets is not None else {}))


def run_price() -> None:
    lines = []

    def P(t=""):
        print(t, flush=True)
        lines.append(t)

    folds = {h: json.loads((WORK / f"fold-{h}.json").read_text()) for h in FOLDS}
    arrs = {h: np.load(WORK / f"fold-{h}.npz") for h in FOLDS}
    prod = json.loads((WORK / "prod274.json").read_text())
    p258 = json.loads((E258 / "work/price258.json").read_text())
    p259 = json.loads((E259 / "work/price259.json").read_text())
    names = NAMES
    spec = {a[0]: a for a in ARMS}
    ka_fold = {h: (sum(r["ok"] for r in folds[h]["known_answer"]), len(folds[h]["known_answer"])) for h in FOLDS}
    ka_prod = (sum(r["ok"] for r in prod["known_answer"]), len(prod["known_answer"]))
    P("=" * 132)
    P("274 -- an nmae-NEUTRAL RENDERING of the union's directions on top of cd4svd5 (PLOLO 162, three cancer folds; kolf2 apart;")
    P("       nmae identity pooled over 4 folds as six229 does; three coherent scenarios with the new-source bands on the union legs)")
    P("=" * 132)
    P("per-fold known-answer rows: " + ", ".join(f"{h} {a}/{b}" for h, (a, b) in ka_fold.items()) + f"; production-frame {ka_prod[0]}/{ka_prod[1]}")
    for h in FOLDS:
        P(f"  {h:14s} {folds[h]['n_targets']} targets, {folds[h]['seconds']}s, peak RSS {folds[h]['peak_rss_gib']:.2f} GiB")
    P(f"  production frame: {prod['seconds']}s, peak RSS {prod['peak_rss_gib']:.2f} GiB")

    agg: dict[str, dict] = {}
    for n in names:
        r = {h: folds[h]["rows"][n] for h in FOLDS}
        a: dict[str, object] = dict(per_fold=r)
        for key in ("cos_agg", "a23_cos", "pds", "pds_A", "a15_sgn_w", "top89_precision", "delivered_norm_A", "cd4_fraction",
                    "union_fraction", "union2_fraction", "preserve_scale"):
            vals = [r[h][key] for h in FOLDS if r[h][key] is not None]
            a[f"{key}_4f"] = float(np.mean(vals)) if vals else None
            vals3 = [r[h][key] for h in CANCER if r[h][key] is not None]
            a[f"{key}_3c"] = float(np.mean(vals3)) if vals3 else None
        for pk in POOLS:
            U = np.concatenate([arrs[h][f"{n}|U_{pk}"] for h in FOLDS])
            M = np.concatenate([arrs[h][f"{n}|M_{pk}"] for h in FOLDS])
            D = np.concatenate([arrs[h][f"{n}|D_{pk}"] for h in FOLDS])
            a[f"nmae_{pk}"] = float(np.mean(1 + (U - 2 * M) / D))
            U3 = np.concatenate([arrs[h][f"{n}|U_{pk}"] for h in CANCER])
            M3 = np.concatenate([arrs[h][f"{n}|M_{pk}"] for h in CANCER])
            D3 = np.concatenate([arrs[h][f"{n}|D_{pk}"] for h in CANCER])
            a[f"nmae_{pk}_3c"] = float(np.mean(1 + (U3 - 2 * M3) / D3))
            if pk == "wald89":
                a["U_sum"], a["M_sum"], a["D_sum"] = float(U.sum()), float(M.sum()), float(D.sum())
                a["U_sum_3c"], a["M_sum_3c"] = float(U3.sum()), float(M3.sum())
        agg[n] = a
    a15n = agg["a15 f=0.10 SHIPPED"]["delivered_norm_A_4f"]
    for n in names:
        agg[n]["delivered_norm_ratio"] = agg[n]["delivered_norm_A_4f"] / a15n

    def delta(a_, b_, key, folds_=CANCER, rel=False):
        vals = [agg[a_]["per_fold"][h][key] / agg[b_]["per_fold"][h][key] - 1 if rel else agg[a_]["per_fold"][h][key] - agg[b_]["per_fold"][h][key]
                for h in folds_]
        return float(np.mean(vals)), [float(v) for v in vals]

    # ---- amplitude: production frame (calibrated) where composable; fold frame A_c for the oracle cells ----
    ratios = [c["ratio_pkg_over_eff"] for c in prod["calibration"] if c["ratio_pkg_over_eff"] is not None]
    amp_cal = float(np.mean(ratios))
    T = prod["terms"]

    def amp_fold(name, base_name):
        return float(np.mean([np.mean([np.sqrt(agg[name]["per_fold"][h][f"ampA_{c}"] / agg[base_name]["per_fold"][h][f"ampA_{c}"]) for c in CTX]) for h in CANCER]))

    def amp_vs(name, base_name):
        if name in T:
            r_eff = amp_ratio(T[name]["per"], T[base_name]["per"])
            return r_eff, 1.0 + amp_cal * (r_eff - 1.0), "production"
        r_eff = amp_fold(name, INC)
        return r_eff, 1.0 + amp_cal * (r_eff - 1.0), "fold (oracle set, not composable on the production panel)"

    # ---- cross-fold known-answer rows: 258 / 259 / 229 / 246 -------------------------------------
    ka = []

    def check(label, got, want, tol=0.0):
        d = abs(float(got) - float(want))
        ok = d <= tol
        ka.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, ok=bool(ok)))
        P(f"  {'PASS' if ok else 'FAIL'} {label:90s} got {got:+.10f} want {want:+.10f} |d| {d:.1e}")

    P("\nCROSS-FOLD KNOWN-ANSWER ROWS")
    for n in ("a15 f=0.10 SHIPPED", "cd4 f=0.34", OLD, ARM6, INC, "stack6 g=0.35", "stack6 g=0.28", "stack6 g=0.2", "stack3 g=0.35", "stack3 g=0.28"):
        for key in ("cos_agg_3c", "pds_3c", "pds_A_3c", "nmae_wald89", "nmae_wald89_3c", "nmae_lfc89", "nmae_h89", "top89_precision_4f", "delivered_norm_ratio", "U_sum", "M_sum"):
            check(f"258 aggregate {n}.{key}", agg[n][key], p258["aggregates"][n][key], 0.0)
    for n in ("B3 g2=0.35", "B3 g2=0.45", "B4+igvf g2=0.35"):
        for key in ("cos_agg_3c", "pds_3c", "nmae_wald89", "U_sum", "M_sum"):
            check(f"259 aggregate {n}.{key}", agg[n][key], p259["aggregates"][n][key], 0.0)
    six229 = json.loads((E229 / "work/six229.json").read_text())["arms"]
    check("six229 stack6 g=0.35 <- ARM nmae_wald89 == ARM6", agg[ARM6]["nmae_wald89"], six229["stack6 g=0.35 <- ARM"]["nmae_wald89"], 5e-6)
    check("six229 h1union nmae_wald89 == h1union", agg[OLD]["nmae_wald89"], six229["h1union g=0.28 <- INCUMBENT"]["nmae_wald89"], 5e-6)
    d_nmae_arm6 = agg[ARM6]["nmae_wald89"] - agg[OLD]["nmae_wald89"]
    check("KA-1 ARM6's identity leg: overall_nmae_vs_incumbent == -0.0047829 (mandate)", -d_nmae_arm6 / NMAE_SPAN / 6, -0.0047829, 5e-8)
    check("KA-1 ARM6's nmae_wald89 == 1.0029112", agg[ARM6]["nmae_wald89"], 1.0029112, 5e-8)
    check("amplitude calibration factor == 258", amp_cal, p258["factors"]["amp_calibration"], 1e-12)
    B_h = board_entry("iUziqE1gewfcPKzYUDaB")
    B_6 = board_entry("ydAemVdWmsTx9xhDWcou")
    B_5 = board_entry("Xv0ayvuoHpsKYWv2iu42")
    pkg = prod["package_terms"]
    cos_h = invert(pkg["h1union"], B_h["raw"]["mse"])
    cos_5 = invert(pkg["cd4svd5"], B_5["raw"]["mse"])
    cos_6 = invert(pkg["ARM6"], B_6["raw"]["mse"])
    raw_5 = B_5["raw"]["mse"]
    base_per = pkg["cd4svd5"]
    check("measured board cosine cd4svd5 == 258", cos_5, p258["factors"]["cos_cd4svd5"], 1e-12)
    for n in ("stack6 g=0.35", "stack3 g=0.35"):
        rel, _ = delta(n, INC, "cos_agg", rel=True)
        dpA, _ = delta(n, INC, "pds_A")
        dpABC, _ = delta(n, INC, "pds")
        d_nmae = agg[n]["nmae_wald89"] - agg[INC]["nmae_wald89"]
        _, r_cal, _ = amp_vs(n, INC)
        c258 = p258["cells"][n]
        check(f"258 cell {n}.amp_calibrated_pct", 100 * (r_cal - 1), c258["amp_calibrated_pct"], 1e-9)
        cons = legs(rel, dpA, d_nmae, r_cal, base_per, cos_5, raw_5, 0.96, 0.7, p258["factors"]["de3_rate_ARM6"], 1.5)
        check(f"258 cell {n}.conservative.total_A", cons["total"], c258["conservative"]["total_A"], 1e-9)
        check(f"258 cell {n}.conc[3].retained", _retention(arrs, n, INC, 3, dpABC)["retained"], c258["conc"]["3"]["retained"], 1e-9)
        c259 = p259["cells"][n]
        sc_ = {}
        for which in ("PESSIMISTIC", "POINT", "OPTIMISTIC"):
            f = scenario_factors(which, rel, dpABC, d_nmae)
            sc_[which] = legs(rel, dpABC, d_nmae, r_cal, base_per, cos_5, raw_5, f["k_cos"], f["k_pds"], f["de3"], f["nmae_x"])
        check(f"259 cell {n}.POINT.total", sc_["POINT"]["total"], c259["scenarios"]["POINT"]["total"], 1e-9)
    n_ok = sum(r["ok"] for r in ka)
    P(f"cross-fold known-answer rows {n_ok}/{len(ka)} pass")
    ka_pass = n_ok == len(ka) and all(a == b for a, b in ka_fold.values()) and ka_prod[0] == ka_prod[1]

    # ---- KA-2: the two board retrodictions (259's R2 form) ----------------------------------------
    P("\nKA-2 -- retrodictions of the two 2026-09-05 board points (a FITTED calibration on pds/mse; the identity's x1.5 and the amplitude leg are the independent content)")
    rel6, _ = delta(ARM6, OLD, "cos_agg", rel=True)
    dp6_ABC, _ = delta(ARM6, OLD, "pds")
    dp6_A, _ = delta(ARM6, OLD, "pds_A")
    rel5, _ = delta(INC, OLD, "cos_agg", rel=True)
    dp5_ABC, _ = delta(INC, OLD, "pds")
    d_nmae6 = agg[ARM6]["nmae_wald89"] - agg[OLD]["nmae_wald89"]
    d_nmae5 = agg[INC]["nmae_wald89"] - agg[OLD]["nmae_wald89"]
    de3_5 = sum(B_5["scaled"][m] - B_h["scaled"][m] for m in ("fid", "reach", "jac")) / 6
    de3_6 = sum(B_6["scaled"][m] - B_h["scaled"][m] for m in ("fid", "reach", "jac")) / 6
    board6 = {m: (B_6["scaled"][m] - B_h["scaled"][m]) / 6 for m in ("pds", "mse", "nmae", "fid", "reach", "jac")}
    board5 = {m: (B_5["scaled"][m] - B_h["scaled"][m]) / 6 for m in ("pds", "mse", "nmae", "fid", "reach", "jac")}
    r6_pkg = amp_ratio(pkg["ARM6"], pkg["h1union"])
    r5_pkg = amp_ratio(pkg["cd4svd5"], pkg["h1union"])
    p_a = legs(rel6, dp6_ABC, d_nmae6, r6_pkg, pkg["h1union"], cos_h, B_h["raw"]["mse"], 0.96, 2.48, 0.0425, 1.5)
    P(f"  (a) BR237 ARM6 by member (k_cos 0.96, k_pds 2.48 ABC, identity x1.5, package amplitude {100*(r6_pkg-1):+.2f} %, DE-three 0.0425): "
      f"pds {p_a['o_pds']:+.5f} (board {board6['pds']:+.5f})  mse {p_a['o_mse']:+.5f} ({board6['mse']:+.5f})  DE3 {p_a['o_de3']:+.5f} ({de3_6:+.5f})  "
      f"nmae {p_a['o_nmae']:+.5f} ({board6['nmae']:+.5f})  TOTAL {p_a['total']:+.5f} (board {B_6['overall']-B_h['overall']:+.5f})")
    ka2a = bool(abs(p_a["total"] - (B_6["overall"] - B_h["overall"])) <= 0.0010 and abs(p_a["o_pds"] - board6["pds"]) <= 0.0008
                and abs(p_a["o_mse"] - board6["mse"]) <= 0.0003 and abs(p_a["o_nmae"] - board6["nmae"]) <= 0.0025)
    p_b = legs(rel5, dp5_ABC, d_nmae5, r5_pkg, pkg["h1union"], cos_h, B_h["raw"]["mse"], 3.76, 1.22, de3_5 / (cos_5 - cos_h), 1.0)
    P(f"  (b) BR249 cd4svd5 vs h1union (k_cos 3.76, k_pds 1.22 ABC, package amplitude {100*(r5_pkg-1):+.2f} %): pds {p_b['o_pds']:+.5f} ({board5['pds']:+.5f})  "
      f"mse {p_b['o_mse']:+.5f} ({board5['mse']:+.5f})  nmae {p_b['o_nmae']:+.5f} ({board5['nmae']:+.5f})  TOTAL {p_b['total']:+.5f} (board {B_5['overall']-B_h['overall']:+.5f}; band +0.0043 +- 0.0008)")
    ka2b = bool(abs(p_b["total"] - 0.0043) <= 0.0008)
    P(f"  KA-2: (a) {'PASS' if ka2a else 'MISS'} (b) {'PASS' if ka2b else 'MISS'}")

    # ---- THE SCREEN ------------------------------------------------------------------------------
    P("\n" + "=" * 132)
    P("THE SCREEN -- every cell vs cd4svd5 (the incumbent), three cancer folds, THREE COHERENT SCENARIOS")
    P("=" * 132)
    chk = raw_for(base_per, cos_5)
    P(f"  197's identity at cd4svd5's own terms reproduces its board raw mse: {chk:.7f} vs {raw_5:.7f}; cos_cd4svd5 {cos_5:.6f}")
    P(f"  bands: k_cos {K_COS_BAND} (opt. {K_COS_OPT}), k_pds {K_PDS_BAND}, DE-three {DE3_BAND}, identity x{NMAE_LOSS_X} on a loss / x0 on a gain (pess.); "
      f"amplitude leg (calibrated effect-level, factor {amp_cal:.3f}) always charged; POINT = {POINT}")
    rng = np.random.default_rng(BOOT_SEED)
    cells = {}
    screen = [a for a in ARMS if a[1] == 5 and a[0] != INC]
    tg = {h: folds[h]["targets"] for h in CANCER}
    base_cos = {h: arrs[h]["cov|base_cos_truth"] for h in CANCER}
    for name, k, resf, g, sel, g2, sel2, variant, rset, rmode, K in screen:
        rel, rel_pf = delta(name, INC, "cos_agg", rel=True)
        relm, _ = delta(name, INC, "a23_cos", rel=True)
        dpA, dpA_pf = delta(name, INC, "pds_A")
        dpABC, dpABC_pf = delta(name, INC, "pds")
        dqw, _ = delta(name, INC, "a15_sgn_w")
        d_nmae = agg[name]["nmae_wald89"] - agg[INC]["nmae_wald89"]
        d_nmae3 = agg[name]["nmae_wald89_3c"] - agg[INC]["nmae_wald89_3c"]
        dU = agg[name]["U_sum"] - agg[INC]["U_sum"]
        dM = agg[name]["M_sum"] - agg[INC]["M_sum"]
        r_eff, r_cal, amp_src = amp_vs(name, INC)
        amp_fold_pct = 100 * (amp_fold(name, INC) - 1)
        kolf = dict(rel_agg=delta(name, INC, "cos_agg", ("kolf2",), rel=True)[0], d_pds=delta(name, INC, "pds", ("kolf2",))[0])
        conc = {kk: _retention(arrs, name, INC, kk, dpABC, targets=tg) for kk in (1, 3, 5)}
        diffs = [arrs[h][f"{name}|pds_pt"] - arrs[h][f"{INC}|pds_pt"] for h in CANCER]
        bt = dict(pooled=boot_folds(diffs, rng), per_fold={h: boot(d, rng) for h, d in zip(CANCER, diffs)})
        # the identity per target (wald89): 1 + (U - 2M)/D, paired over the 4-fold pool and the 3 cancer folds
        nm_d = {h: (1 + (arrs[h][f"{name}|U_wald89"] - 2 * arrs[h][f"{name}|M_wald89"]) / arrs[h][f"{name}|D_wald89"])
                - (1 + (arrs[h][f"{INC}|U_wald89"] - 2 * arrs[h][f"{INC}|M_wald89"]) / arrs[h][f"{INC}|D_wald89"]) for h in FOLDS}
        bt["nmae_4f"] = boot(np.concatenate([nm_d[h] for h in FOLDS]), rng)
        bt["nmae_3c"] = boot(np.concatenate([nm_d[h] for h in CANCER]), rng)
        # 262's groups on the base cosine: (-1, 0] anti-aligned, (0, 0.05] marginal, > 0.05 saturated; targets crossing the 0.02 knee
        groups = {}
        for gname, lo, hi in (("anti", -np.inf, 0.0), ("marginal", 0.0, 0.05), ("saturated", 0.05, np.inf)):
            vals, ns = [], []
            for h, d in zip(CANCER, diffs):
                m = (base_cos[h] > lo) & (base_cos[h] <= hi)
                vals.append(float(d[m].mean()) if m.any() else float("nan"))
                ns.append(int(m.sum()))
            groups[gname] = dict(d_pds=float(np.nanmean(vals)), per_fold=vals, n=ns)
        up = down = 0
        for h in CANCER:
            c0 = arrs[h][f"{INC}|cos_pt"]
            c1 = arrs[h][f"{name}|cos_pt"]
            up += int(np.sum((c0 <= 0.02) & (c1 > 0.02)))
            down += int(np.sum((c0 > 0.02) & (c1 <= 0.02)))
        sc_ = {}
        for which in ("PESSIMISTIC", "POINT", "OPTIMISTIC"):
            f = scenario_factors(which, rel, dpABC, d_nmae)
            sc_[which] = legs(rel, dpABC, d_nmae, r_cal, base_per, cos_5, raw_5, f["k_cos"], f["k_pds"], f["de3"], f["nmae_x"])
        cons258 = legs(rel, dpA, d_nmae, r_cal, base_per, cos_5, raw_5, 0.96, 0.7, DE3_BAND[0], 1.5)
        point_A = legs(rel, dpA, d_nmae, r_cal, base_per, cos_5, raw_5, 0.96, K_PDS_A_POINT, 0.24, 1.0)
        rd = {h: folds[h]["rows"][name].get("renderer") for h in FOLDS}
        cells[name] = dict(arm=name, g=g, g2=g2, selection=sel, selection2=sel2, variant=variant, rset=rset, rmode=rmode, K=K,
                           oracle=bool(rset == "wald89"), control=bool(variant is not None or name.startswith("cd4svd5 IDENT")),
                           carrying1={h: folds[h]["rows"][name]["union_carrying"] for h in FOLDS},
                           carrying2={h: folds[h]["rows"][name]["union2_carrying"] for h in FOLDS},
                           renderer={h: rd[h] for h in FOLDS}, prod_info=prod["compose_info"].get(name), prod_renderer=prod["renderer"].get(name),
                           rel_cos_agg=rel, rel_cos_agg_per_fold=rel_pf, rel_macro=relm,
                           d_pds_A=dpA, d_pds_A_per_fold=dpA_pf, d_pds_ABC=dpABC, d_pds_ABC_per_fold=dpABC_pf,
                           pds_pos_folds=int(sum(v > 0 for v in dpABC_pf)), pds_over_loo_rms=abs(dpABC) / SIGMA_PDS,
                           boot=bt, groups=groups, cross_up=up, cross_down=down,
                           d_qw=dqw, d_nmae_wald89=d_nmae, d_nmae_wald89_3c=d_nmae3, dU=dU, dM=dM, e=(dM / dU if dU != 0 else float("nan")),
                           d_nmae_lfc89=agg[name]["nmae_lfc89"] - agg[INC]["nmae_lfc89"], d_nmae_h89=agg[name]["nmae_h89"] - agg[INC]["nmae_h89"],
                           d_precision89=agg[name]["top89_precision_4f"] - agg[INC]["top89_precision_4f"],
                           amp_effect_level_pct=100 * (r_eff - 1), amp_calibrated_pct=100 * (r_cal - 1), amp_source=amp_src, amp_fold_pct=amp_fold_pct,
                           preserve_3c=agg[name]["preserve_scale_3c"],
                           kolf2=kolf, conc=conc, scenarios=sc_, cons258=cons258, point_A=point_A,
                           clears_pessimistic=bool(sc_["PESSIMISTIC"]["total"] >= GATE),
                           retention_ok=bool(conc[3]["retained"] >= 0.5),
                           pds_ok=bool(all(v > 0 for v in dpABC_pf) and bt["pooled"]["lo"] > 0),
                           amp_ok=bool(abs(100 * (r_cal - 1)) <= AMP_TOL_PCT))
    # controls below half: for every real rendered cell, the same renderer on the random-direction controls
    for n, c in cells.items():
        if c["control"] or c["rset"] is None:
            c["controls_ok"] = None
            c["control_rows"] = []
            continue
        ctrls = [m for m, d in cells.items() if d["variant"] is not None and d["rset"] == c["rset"] and d["rmode"] == c["rmode"] and d["K"] == c["K"]
                 and d["selection"] == c["selection"] and d["g"] == c["g"] and d["g2"] == c["g2"]]
        c["control_rows"] = ctrls
        if not ctrls:
            c["controls_ok"] = None
        else:
            ok = all(cells[m]["d_pds_ABC"] < 0.5 * c["d_pds_ABC"] and cells[m]["scenarios"]["PESSIMISTIC"]["total"] < 0.5 * c["scenarios"]["PESSIMISTIC"]["total"] for m in ctrls) \
                if c["d_pds_ABC"] > 0 and c["scenarios"]["PESSIMISTIC"]["total"] > 0 else False
            c["controls_ok"] = bool(ok)
        c["cleared"] = bool(c["clears_pessimistic"] and c["pds_ok"] and c["retention_ok"] and c["amp_ok"] and c["controls_ok"] is True and not c["oracle"])

    def table(title, subset):
        P(f"\n{title}")
        P(f"  {'cell':46s}{'agg%':>8s}{'mac%':>8s}{'dpds ABC':>10s}{'+':>2s}{'boot lo':>9s}{'boot hi':>9s}{'dpds A':>9s}{'dnmae':>9s}{'dnmae3c':>9s}{'e':>7s}"
          f"{'dprec':>8s}{'amp%':>7s}{'|':>2s}{'PESS':>9s}{'POINT':>9s}{'OPT':>9s}{'|':>2s}{'C258':>8s}{'ret3':>7s}{'kolf dpds':>10s}")
        for n in subset:
            c = cells[n]
            s = c["scenarios"]
            P(f"  {n:46s}{100*c['rel_cos_agg']:>+8.3f}{100*c['rel_macro']:>+8.3f}{c['d_pds_ABC']:>+10.5f}{c['pds_pos_folds']:>2d}"
              f"{c['boot']['pooled']['lo']:>+9.5f}{c['boot']['pooled']['hi']:>+9.5f}{c['d_pds_A']:>+9.5f}{c['d_nmae_wald89']:>+9.5f}{c['d_nmae_wald89_3c']:>+9.5f}{c['e']:>7.3f}"
              f"{c['d_precision89']:>+8.4f}{c['amp_calibrated_pct']:>+7.2f}{'|':>2s}"
              f"{s['PESSIMISTIC']['total']:>+9.5f}{s['POINT']['total']:>+9.5f}{s['OPTIMISTIC']['total']:>+9.5f}{'|':>2s}{c['cons258']['total']:>+8.5f}"
              f"{100*c['conc'][3]['retained']:>6.0f}%{c['kolf2']['d_pds']:>+10.5f}")

    P("\ncolumns: agg%/mac% = rel cos_agg / macro cosine vs cd4svd5 (3 cancer folds); dpds ABC = 246's functional (+ = folds > 0; boot = paired target bootstrap, "
      "stratified by fold, 20,000 draws, 95 %); dpds A = six229's frame-A functional; dnmae = the exact identity (wald89; 4 folds / 3 cancer); e = dM/dU (break-even 0.5); "
      "dprec = top-89 precision; amp% = calibrated delivered amplitude (production frame; fold frame for the oracle set); PESS/POINT/OPT = the coherent six-member "
      "totals; C258 = 258's conservative set (frame A, k_pds 0.7, x1.5); ret3 = leave-top-3-out pds retention")
    main_cells = [n for n in names if n.startswith("stack6 g=0.35 |")]
    table("0. BASELINES -- the unions at their OWN magnitudes (R0; 258 / 259's rows)", ["stack6 g=0.35", "stack6 g=0.28", "stack6 g=0.2", "stack3 g=0.35", "stack3 g=0.28", "B3 g2=0.35", "B3 g2=0.45", "B4+igvf g2=0.35"])
    table("1. R1 / R2 -- the six-source union at g = 0.35, every set x mode (wald89 = ORACLE ceiling)", main_cells)
    for base in ("stack3 g=0.35", "stack3 g=0.28", "stack6 g=0.28", "stack6 g=0.2"):
        table(f"2. R3 -- {base} rendered", [base] + [n for n in names if n.startswith(f"{base} |")])
    table("3. R4 -- 259's construction B (pair kept at 0.28 + union2) rendered", [n for n in names if n.startswith("B3 ") or n.startswith("B4")])
    table("4. CONTROLS -- the identity check (must be 0 everywhere) and energy-matched random directions rendered the same way",
          [n for n in names if n.startswith("cd4svd5 IDENT") or "PERM" in n or "SHUF" in n])

    # ---- the union's magnitude on the pool, before and after rendering (the mechanism) -----------
    P("\n5. THE MECHANISM -- effect-level L1 on the truth's wald89 pool per row (mean over carried rows, 3 cancer folds) and the identity's split")
    P(f"  {'cell':46s}{'poolL1/inc':>11s}{'dU/D':>9s}{'dM/D':>9s}{'dnmae':>9s}{'d own cos':>10s}{'cross up':>9s}{'cross dn':>9s}{'anti':>9s}{'marg':>9s}{'sat':>9s}")
    for n in ["stack6 g=0.35"] + main_cells + ["stack3 g=0.35", "stack3 g=0.35 | all rank", "B3 g2=0.35", "B3 g2=0.35 | all rank", "B4+igvf g2=0.35 | wald89 exact"]:
        c = cells[n]
        ratios_ = []
        dcos = []
        for h in CANCER:
            ch = arrs[h][f"{n}|changed"].astype(bool)
            a_, b_ = arrs[h][f"{n}|poolL1_wald89"], arrs[h][f"{INC}|poolL1_wald89"]
            m = ch & (b_ > 0)
            ratios_.append(float(np.mean(a_[m] / b_[m])) if m.any() else float("nan"))
            dcos.append(float(np.mean(arrs[h][f"{n}|cos_pt"] - arrs[h][f"{INC}|cos_pt"])))
        D = agg[INC]["D_sum"]
        P(f"  {n:46s}{np.nanmean(ratios_):>11.4f}{c['dU']/D:>+9.5f}{c['dM']/D:>+9.5f}{c['d_nmae_wald89']:>+9.5f}{np.mean(dcos):>+10.5f}{c['cross_up']:>9d}{c['cross_down']:>9d}"
          f"{c['groups']['anti']['d_pds']:>+9.4f}{c['groups']['marginal']['d_pds']:>+9.4f}{c['groups']['saturated']['d_pds']:>+9.4f}")
    P("  poolL1/inc = mean over the union's carried rows of (sum of |effect| on the truth's wald89 pool) / the incumbent's; dU/D, dM/D = the identity's unmatched / "
      "matched shifts over the pooled D; cross up/dn = targets crossing own cosine 0.02 upward / downward (3 folds); anti/marg/sat = d pds by 262's base-cosine group")

    # ---- the best cell in full ---------------------------------------------------------------------
    real = [n for n, c in cells.items() if not c["control"] and not c["oracle"] and c["rset"] is not None]
    oracle = [n for n, c in cells.items() if c["oracle"] and not c["control"]]
    best = max(real, key=lambda n: cells[n]["scenarios"]["PESSIMISTIC"]["total"])
    best_or = max(oracle, key=lambda n: cells[n]["scenarios"]["PESSIMISTIC"]["total"])
    P("\n" + "=" * 132)
    P(f"THE BEST BUILDABLE CELL: {best}   (best oracle cell: {best_or} at {cells[best_or]['scenarios']['PESSIMISTIC']['total']:+.5f} pessimistic)")
    P("=" * 132)
    for n in (best, best_or):
        c = cells[n]
        P(f"\n  {n}")
        P(f"  PLOLO vs cd4svd5: rel cos_agg {100*c['rel_cos_agg']:+.4f}% (per fold {['%+.3f' % (100*v) for v in c['rel_cos_agg_per_fold']]}), rel macro {100*c['rel_macro']:+.4f}%")
        P(f"  d pds ABC {c['d_pds_ABC']:+.6f} (per fold {['%+.5f' % v for v in c['d_pds_ABC_per_fold']]}; {c['pds_pos_folds']}/3 > 0; pooled bootstrap [{c['boot']['pooled']['lo']:+.5f}, {c['boot']['pooled']['hi']:+.5f}], "
          f"P(<= 0) {c['boot']['pooled']['frac_le0']:.3f}; per-fold intervals " + "; ".join(f"{h} [{c['boot']['per_fold'][h]['lo']:+.4f}, {c['boot']['per_fold'][h]['hi']:+.4f}]" for h in CANCER) + f"); frame A {c['d_pds_A']:+.6f}")
        P(f"  nmae identity: d wald89 {c['d_nmae_wald89']:+.6f} (3 cancer {c['d_nmae_wald89_3c']:+.6f}; paired bootstrap 4f [{c['boot']['nmae_4f']['lo']:+.5f}, {c['boot']['nmae_4f']['hi']:+.5f}]); "
          f"lfc89 {c['d_nmae_lfc89']:+.6f} h89 {c['d_nmae_h89']:+.6f}; dU {c['dU']:+.2f} dM {c['dM']:+.2f} e {c['e']:.3f}; d precision89 {c['d_precision89']:+.5f}; d q_w {c['d_qw']:+.5f}")
        P(f"  amplitude vs cd4svd5: {c['amp_source']} effect-level {c['amp_effect_level_pct']:+.3f}%, calibrated {c['amp_calibrated_pct']:+.3f}% (fold-frame read {c['amp_fold_pct']:+.3f}%); preserve 3c {c['preserve_3c']:.6f} vs {agg[INC]['preserve_scale_3c']:.6f}")
        P(f"  renderer: " + "; ".join(f"{h} rows {c['renderer'][h]['n_changed']} scale med {c['renderer'][h]['scale_median']:.4f} [{c['renderer'][h]['scale_min']:.3f}, {c['renderer'][h]['scale_max']:.3f}] rank rho {c['renderer'][h]['rank_rho_mean']:.3f}" for h in CANCER if c['renderer'][h]))
        P(f"  262 groups (d pds): anti {c['groups']['anti']['d_pds']:+.4f} (n {c['groups']['anti']['n']}), marginal {c['groups']['marginal']['d_pds']:+.4f} (n {c['groups']['marginal']['n']}), saturated {c['groups']['saturated']['d_pds']:+.4f}; "
          f"targets crossing 0.02 own cosine: up {c['cross_up']} / down {c['cross_down']}")
        P(f"  kolf2 (apart): rel cos_agg {100*c['kolf2']['rel_agg']:+.3f}%, d pds {c['kolf2']['d_pds']:+.5f}")
        P(f"  leave-top-k-out pds retention: k=1 {100*c['conc'][1]['retained']:.1f}% (share {100*c['conc'][1]['share']:.0f}%), k=3 {100*c['conc'][3]['retained']:.1f}% (share {100*c['conc'][3]['share']:.0f}%), k=5 {100*c['conc'][5]['retained']:.1f}%; "
          f"top-3 per fold: " + " | ".join(f"{h}: {', '.join(c['conc'][3]['top_names'][h])}" for h in CANCER))
        P(f"  {'scenario':14s}{'k_cos':>6s}{'k_pds':>6s}{'DE3':>7s}{'nmae x':>7s}{'d cos':>10s}{'raw mse':>9s}{'| pds':>9s}{'mse':>9s}{'nmae':>9s}{'reach':>9s}{'fid':>9s}{'jac':>9s}{'| amp':>9s}{'TOTAL':>10s}")
        for which in ("PESSIMISTIC", "POINT", "OPTIMISTIC"):
            s = c["scenarios"][which]
            P(f"  {which:14s}{s['k_cos']:>6.2f}{s['k_pds']:>6.2f}{s['de3_rate']:>7.4f}{s['nmae_x']:>7.1f}{s['dcos']:>+10.6f}{s['raw']:>9.5f}{s['o_pds']:>+9.5f}{s['o_mse']:>+9.5f}{s['o_nmae']:>+9.5f}"
              f"{s['o_reach']:>+9.5f}{s['o_fid']:>+9.5f}{s['o_jac']:>+9.5f}{s['o_amp']:>+9.5f}{s['total']:>+10.5f}")
        P(f"  258's conservative (frame A, k_pds 0.7, x1.5) {c['cons258']['total']:+.5f}; point (frame A, k_pds 3.79, DE3 0.24) {c['point_A']['total']:+.5f}")
        P(f"  conditions: pessimistic >= +{GATE} {'yes' if c['clears_pessimistic'] else 'NO'}; pds > 0 on 3/3 with the pooled interval above 0 {'yes' if c['pds_ok'] else 'NO'}; "
          f"retention >= 50 % {'yes' if c['retention_ok'] else 'NO'}; amplitude within {AMP_TOL_PCT} % {'yes' if c['amp_ok'] else 'NO'}; controls below half "
          f"{c['controls_ok'] if c['controls_ok'] is not None else 'n/a'} {c['control_rows']}; oracle {c['oracle']}")

    # ---- verdict ----------------------------------------------------------------------------------
    pess_clear = [n for n in real if cells[n]["clears_pessimistic"]]
    point_clear = [n for n in real if cells[n]["scenarios"]["POINT"]["total"] >= GATE]
    opt_clear = [n for n in real if cells[n]["scenarios"]["OPTIMISTIC"]["total"] >= GATE]
    all_clear = [n for n in real if cells[n]["cleared"]]
    or_pess = [n for n in oracle if cells[n]["clears_pessimistic"]]
    P("\n" + "=" * 132)
    P("VERDICT")
    P("=" * 132)
    P(f"  known-answer rows: per fold " + ", ".join(f"{h} {a}/{b}" for h, (a, b) in ka_fold.items()) + f"; production {ka_prod[0]}/{ka_prod[1]}; cross-fold {n_ok}/{len(ka)} -> {'PASS' if ka_pass else 'FAIL'}")
    P(f"  KA-1 (ARM6's identity leg -0.0047829): {'PASS' if any(r['label'].startswith('KA-1') and r['ok'] for r in ka) else 'FAIL'}; KA-2: (a) BR237 {'PASS' if ka2a else 'MISS'} (b) BR249 {'PASS' if ka2b else 'MISS'}")
    P(f"  buildable rendered cells: {len(real)}; PESSIMISTIC >= +{GATE}: {len(pess_clear)} {pess_clear}; POINT >= gate: {len(point_clear)}; OPTIMISTIC >= gate: {len(opt_clear)}")
    P(f"  ALL conditions (pess. gate AND pds 3/3 with interval > 0 AND retention >= 50 % AND amplitude within 0.5 % AND controls below half): {len(all_clear)} {all_clear}")
    P(f"  oracle (wald89) cells: {len(oracle)}; pessimistic >= gate: {len(or_pess)} {or_pess}; best oracle {best_or} {cells[best_or]['scenarios']['PESSIMISTIC']['total']:+.5f}")
    P(f"  best buildable cell {best}: pessimistic {cells[best]['scenarios']['PESSIMISTIC']['total']:+.5f}, point {cells[best]['scenarios']['POINT']['total']:+.5f}, optimistic {cells[best]['scenarios']['OPTIMISTIC']['total']:+.5f}")
    verdict = "CLEARS" if all_clear else "DOES NOT CLEAR"
    P(f"\n{verdict}")
    out = dict(known_answer=ka, ka_pass=ka_pass, ka2a=p_a, ka2a_pass=ka2a, ka2b=p_b, ka2b_pass=ka2b,
               factors=dict(cos_h1union=cos_h, cos_ARM6=cos_6, cos_cd4svd5=cos_5, amp_calibration=amp_cal, rel6=rel6, dp6_ABC=dp6_ABC, dp6_A=dp6_A, d_nmae_arm6=d_nmae6),
               aggregates={n: {k: v for k, v in agg[n].items() if k != "per_fold"} for n in names},
               cells=cells, best=best, best_oracle=best_or, pess_clear=pess_clear, point_clear=point_clear, opt_clear=opt_clear, all_clear=all_clear, verdict=verdict,
               boot=dict(seed=BOOT_SEED, draws=BOOT_B))
    (WORK / "price274.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "price274.log").write_text("\n".join(lines) + "\n")
    P(f"\nwrote {WORK/'price274.json'} and price274.log")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold")
    ap.add_argument("--prod", action="store_true")
    ap.add_argument("--price", action="store_true")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.list:
        for n in NAMES:
            print(n)
        print(len(NAMES), "arms")
    elif a.fold:
        run_fold(a.fold)
    elif a.prod:
        run_prod()
    elif a.price:
        run_price()
    else:
        ap.error("one of --fold <held>, --prod, --price, --list")


if __name__ == "__main__":
    main()
