"""259 -- the union WITHOUT `igvf_cardiomyocyte` on top of `cd4svd5`, in the two constructions 258 did
not test, priced on COHERENT PESSIMISTIC / POINT / OPTIMISTIC scenarios.  MEASUREMENT ONLY: no build, no
`vcc`, no `cell-eval2`, nothing under `submissions/`, nothing written outside this directory.

258's instrument, reused (its code paths copied, its fold caches byte-copied, its numbers the known-answer
rows), plus the builder's SECOND union residual (248's `extra2` / `audits_out["union2"]` path):

  --fold <held>   stage one PLOLO fold (162), build every arm through the PRODUCTION builder's own
                  `apply_frozen_cd4_residual` (deflate k=5 -> size CD4 -> union 1 at g -> union 2 at g2 ->
                  ONE merged add -> ONE preserve; 248/249's route, module globals as the knobs) and score
                  it on (a) 218's `score_effects` functional (cos_agg / macro / pds, frames A,B,C) and
                  (b) 216/229's frame-A functional (188's exact nmae identity on three pools, top-89
                  precision, delivered norm, pds on frame A).  Known-answer rows: 258's own fold rows and
                  per-target arrays (must be bit-identical), and 248's union2 rows (repro248-<fold>.json).
  --prod          production-frame reconstruction off 249's staged artifacts (258 s1): re-compose cd4svd5
                  (a28d5be2...), ARM6 (c5721825...), 248's staged cum1 (6225a582...) and cum4 (2f9eba51...)
                  through the union2 path, then every screen arm on the real 300-target panel for its
                  delivered amplitude (211/243's A_c/B_c terms), calibrated on the three packages.
  --price         R1-R3, the A and B tables, the controls, the best cells, leave-top-3-out, the three
                  scenarios, CLEARS / DOES NOT CLEAR -> work/price259.{json,log}.

    uv run python experiments/259-union-without-igvf/price259.py --fold k562_gwps   (x4 folds, serially)
    uv run python experiments/259-union-without-igvf/price259.py --prod
    uv run python experiments/259-union-without-igvf/price259.py --price
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
# The builder reads its knobs at import.  This is the incumbent's configuration (258's env); every arm
# rebinds the module globals it varies (248's pattern) and nothing else.
os.environ["VCC_ALPHA"] = "0.982"
os.environ["VCC_CD4_TARGET_F"] = "0.34"
os.environ["VCC_H1_UNION_G"] = "0.28"
os.environ["VCC_CD4_DEFLATE_K"] = "5"
for _k in ("VCC_CD4_SIZING", "VCC_CD4_TARGET_G", "VCC_UNION_SOURCES", "VCC_UNION2_SOURCES",
           "VCC_UNION2_G", "VCC_CD4_DEFLATE_BASIS", "VCC_CD4_PARTIAL_G", "VCC_CD4_FRACTION",
           "VCC_SIGN_REPAIR"):
    os.environ.pop(_k, None)

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
E225 = ROOT / "experiments/225-source-sweep"
E229 = ROOT / "experiments/229-stack4-build"
E243 = ROOT / "experiments/243-transfer-factor"
E248 = ROOT / "experiments/248-second-residual"
E249 = ROOT / "experiments/249-cd4-deflate"
E255 = ROOT / "experiments/255-cd4-deflate-k10"
E258 = ROOT / "experiments/258-union-dose-on-cd4svd5"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/188-magnitude-calibration"))
sys.path.insert(0, str(ROOT / "experiments/203-reopen-125"))
sys.path.insert(0, str(E243))

import plolo  # noqa: E402
from plolo import (E062, FOLDS, OFFICIAL, PANEL, PROD_OPS, Operator, axis_sha256,  # noqa: E402
                   build_fold, cd4_for, deliver, endpoints)
from mapfit import to_lfc  # noqa: E402
from phase23 import POOL_N, key_precision, topn, umd  # noqa: E402
from lib243 import CP46, KC, SUMD  # noqa: E402
from cell_eval2.metrics.delta import _across_pert_budget, _row_weights  # noqa: E402

plolo.WORK = WORK / "plolo"        # AFTER every import: 259's OWN byte copies of the fold caches (= 258's = 249's)

ALPHA, GAMMA, LAM = 0.982, 0.65, 0.0
NMAE_SPAN = 0.60                 # 188/216/229's identity span (overall = -d nmae / 0.60 / 6)
SPAN_MSE = 0.950162566           # the scorer's mse span (243/248/255/258)
PDS_SCALE = 2.2138               # overall pds = 2.2138 x k_pds x d pds / 6
SIGMA_PDS = 0.001494             # 162's leave-one-arm-out sigma on raw PLOLO pds
GATE = 0.0019                    # shape-only gate
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
CTX = ("A", "B", "C")
N_GENES = 18533
POOLS = ("wald89", "lfc89", "h89")
TS = plolo.TS

# ---- the measured bands (mandate s"The price"; 258 s3 for the ARM6/cd4svd5 realised values) ----------
K_PDS_BAND = (0.70, 2.48)        # 243's class value (h1union realised) / ARM6-realised A-B-C mean
K_COS_BAND = (0.94, 3.76)        # new-source into empty rows / re-shape of occupied rows (BR249)
DE3_BAND = (0.0425, 0.24)        # ARM6-realised DE-three per unit board d cos / 243's recipe
NMAE_LOSS_BAND = (1.0, 1.5)      # the identity on a LOSS (BR237 realised 1.497x); a GAIN is credited 0 (pess.) / 1.0 (opt.)
POINT = dict(k_cos=0.96, k_pds=2.48, de3=0.0425, nmae_x=1.0)   # ARM6-realised for the union legs
K_PDS_A_POINT = 3.79             # ARM6-realised on the frame-A functional (reported alongside, never the verdict)
RATE218 = dict(reach=0.002446, fid=0.000847, jac=0.000100)      # 218 s7.4 rate card, SHARES only (257's split)
DE3_SPLIT = {m: RATE218[m] / sum(RATE218.values()) for m in RATE218}

PAIR = ("gse295214_h1", "h1_2025")
DLD1, JIANG, IGVF, MYL = ("dld1_lowmoi",), ("jiang_5line_noK562",), ("igvf_cardiomyocyte",), ("myllia_a549",)
SEL = {
    "pair": PAIR,
    "stack3": PAIR + DLD1,
    "pair+igvf": PAIR + IGVF,
    "stack6": PAIR + DLD1 + JIANG + IGVF + MYL,
    "union5": PAIR + DLD1 + JIANG + MYL,            # construction A: the single union without igvf
    "u2_dld1": DLD1, "u2_jiang": JIANG, "u2_myllia": MYL,
    "u2_three": DLD1 + JIANG + MYL,                 # construction B: the second residual without igvf
    "u2_four": DLD1 + JIANG + IGVF + MYL,           # control (i), igvf back, in 248's CUM4 order (its KA row)
    "u2_four_mandate": DLD1 + JIANG + MYL + IGVF,   # the mandate's spelling of the same set (difference measured)
}
DOSES_A = (0.20, 0.28, 0.35, 0.45)
DOSES_B = (0.15, 0.28, 0.35, 0.45)
DOSES_CTRL = (0.15, 0.28, 0.35)
# (name, deflate k, cd4 f, union-1 g, union-1 selection, union-2 g2, union-2 selection, variant)
ARMS: list[tuple] = [
    ("a15 f=0.10 SHIPPED", None, 0.10, None, None, None, None, None),
    ("h1union g=0.28", None, 0.34, 0.28, "pair", None, None, None),
    ("stack6 g=0.35 on h1union", None, 0.34, 0.35, "stack6", None, None, None),
    ("cd4svd5 <- INCUMBENT", 5, 0.34, 0.28, "pair", None, None, None),
    ("stack6 g=0.35", 5, 0.34, 0.35, "stack6", None, None, None),
    ("stack3 g=0.35", 5, 0.34, 0.35, "stack3", None, None, None),
    ("pair+igvf g=0.28", 5, 0.34, 0.28, "pair+igvf", None, None, None),
    # R3 -- 248's second-residual rows on the h1union base (deflation off)
    ("248 CUM1 dld1 g2=0.28 on h1union", None, 0.34, 0.28, "pair", 0.28, "u2_dld1", None),
    ("248 CUM1 dld1 g2=0.35 on h1union", None, 0.34, 0.28, "pair", 0.35, "u2_dld1", None),
    ("248 CUM4 four g2=0.35 on h1union", None, 0.34, 0.28, "pair", 0.35, "u2_four", None),
]
ARMS += [(f"A5 g={g:g}", 5, 0.34, g, "union5", None, None, None) for g in DOSES_A]
ARMS += [(f"B3 g2={g:g}", 5, 0.34, 0.28, "pair", g, "u2_three", None) for g in DOSES_B]
ARMS += [(f"B4+igvf g2={g:g}", 5, 0.34, 0.28, "pair", g, "u2_four", None) for g in DOSES_CTRL]
ARMS += [(f"B1 dld1 g2={g:g}", 5, 0.34, 0.28, "pair", g, "u2_dld1", None) for g in DOSES_CTRL]
ARMS += [(f"B1 jiang g2={g:g}", 5, 0.34, 0.28, "pair", g, "u2_jiang", None) for g in (0.28, 0.35)]
ARMS += [(f"B1 myllia g2={g:g}", 5, 0.34, 0.28, "pair", g, "u2_myllia", None) for g in (0.28, 0.35)]
for _g in DOSES_CTRL:
    ARMS += [(f"B3 g2={_g:g} SHUF20726", 5, 0.34, 0.28, "pair", _g, "u2_three", "SHUF20726"),
             (f"B3 g2={_g:g} SHUF771103", 5, 0.34, 0.28, "pair", _g, "u2_three", "SHUF771103"),
             (f"B3 g2={_g:g} CD4ROWS", 5, 0.34, 0.28, "pair", _g, "u2_three", "CD4ROWS")]
INC = "cd4svd5 <- INCUMBENT"
OLD = "h1union g=0.28"
ARM6 = "stack6 g=0.35 on h1union"
SHARED_258 = ("a15 f=0.10 SHIPPED", OLD, ARM6, INC, "stack6 g=0.35", "stack3 g=0.35", "pair+igvf g=0.28")
R3_ROWS = {  # my name -> (248 arm, 248 g)
    "248 CUM1 dld1 g2=0.28 on h1union": ("CUM1", 0.28),
    "248 CUM1 dld1 g2=0.35 on h1union": ("CUM1", 0.35),
    "248 CUM4 four g2=0.35 on h1union": ("CUM4", 0.35),
}


def peak_gib() -> float:
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024**3


def load_builder():
    path = ROOT / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("price259_builder", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["price259_builder"] = module
    spec.loader.exec_module(module)
    return module, path


def sources_of(MOD, names):
    return tuple((n, MOD.UNION_SOURCE_REGISTRY[n]["path"]) for n in names)


def shuffled(res, seed: int):
    """218's target-shuffle control (`resid218.shuffled`, verbatim): a derangement of the carried rows."""
    rng = np.random.default_rng(seed)
    n = res.mean_effect.shape[0]
    carry = np.flatnonzero(np.asarray(res.target_condition_count) > 0)
    perm = np.arange(n)
    sh = carry.copy()
    while True:
        rng.shuffle(sh)
        if len(carry) < 2 or not np.any(sh == carry):
            break
    perm[carry] = sh
    return E062.CD4Prepared(mode=res.mode + "+shuffled", mean_effect=np.asarray(res.mean_effect)[perm].copy(),
                            gene_votes=np.asarray(res.gene_votes)[perm].copy(),
                            target_condition_count=np.asarray(res.target_condition_count)[perm].copy(),
                            condition_details=res.condition_details)


def restrict_rows(src, like):
    """225's F3 row-matched control (`compose218.restrict_rows`, verbatim): `src` kept on the rows `like` carries."""
    keep = np.asarray(like.target_condition_count) > 0
    eff = np.asarray(src.mean_effect).copy()
    votes = np.asarray(src.gene_votes).copy()
    counts = np.asarray(src.target_condition_count).copy()
    eff[~keep] = 0.0
    votes[~keep] = 0
    counts[~keep] = 0
    return E062.CD4Prepared(mode=f"{src.mode}|rows_of_{like.mode}", mean_effect=eff, gene_votes=votes,
                            target_condition_count=counts, condition_details=[])


def compose(MOD, wh, cd4, wcols, panel_positions, k, resf, union, g, union2=None, g2=None):
    """Production's own path (`apply_frozen_cd4_residual`): CD4 deflated at k (None = off) BEFORE sizing at
    f = resf; union 1 sized by the carrying-row rule at g; union 2 (248's `extra2`) sized by the same rule
    at g2 and accumulated after union 1; ONE merged add, ONE preserve.  The module globals ARE the builder's
    knobs (248's pattern); `UNION2_G` is read at call time exactly as `H1_UNION_G` is."""
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
    if u2 is not None:
        assert float(u2["merged_preserve"]["post_add_global_preserve_scale"]) == preserve
    if u_audit is not None and u2 is not None:
        assert float(u_audit["merged_preserve"]["post_add_global_preserve_scale"]) == preserve
    info = dict(
        cd4_fraction=float(cd4_audit["fraction"]),
        union_fraction=(None if u_audit is None else float(u_audit["residual_frame"]["solved_fraction"])),
        union_carrying=(None if u_audit is None else int(u_audit["residual_frame"]["carrying_rows"])),
        union_coverage=(None if u_audit is None else float(u_audit["residual_frame"]["coverage_frac"])),
        union2_fraction=(None if u2 is None else float(u2["residual_frame"]["solved_fraction"])),
        union2_carrying=(None if u2 is None else int(u2["residual_frame"]["carrying_rows"])),
        union2_coverage=(None if u2 is None else float(u2["residual_frame"]["coverage_frac"])),
        union2_shared_targets=(None if u2 is None else u2["residual_frame"].get("targets_shared_with_h1_union")),
        union2_outside_cd4=(None if u2 is None else u2["residual_frame"].get("targets_outside_cd4_support")),
        preserve_scale=preserve,
        energy_removed=(None if "deflation" not in cd4_audit else float(cd4_audit["deflation"]["energy_removed_frac"])),
        n_carried=(None if "deflation" not in cd4_audit else int(cd4_audit["deflation"]["n_carried"])),
    )
    return final, info


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
    A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    B = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S).copy()
    pds = ((S < own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1)) / (n - 1)
    return pds, own


def score_arm(final, fold, Ts, ctx, LR, pools):
    """Both functionals from six renders (2 tags x 3 frames) -- 258's `score_arm`, verbatim."""
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
    return out, arrays


def run_fold(held: str) -> None:
    t0 = time.time()
    assert plolo.WORK == WORK / "plolo", plolo.WORK
    MOD, builder_path = load_builder()
    assert MOD.CD4_DEFLATE_K == 5 and MOD.H1_UNION_G == 0.28 and MOD.CD4_TARGET_F == 0.34
    for name in ("CD4_SIZING", "UNION2_G", "UNION2_SOURCES", "CD4_DEFLATE_BASIS", "CD4_PARTIAL_G"):
        assert getattr(MOD, name) is None, f"builder global {name} is live"
    print(f"builder {builder_path.relative_to(ROOT)} sha256 {MOD.file_sha256(builder_path)[:16]}", flush=True)

    fold = build_fold(held)
    targets = fold["targets"]
    truth_h = np.asarray(fold["truth"], dtype=np.float64)
    pre = np.asarray(fold["pre"], dtype=np.float64)
    cd4 = cd4_for(targets)
    needed = sorted({a[4] for a in ARMS if a[4] is not None} | {a[6] for a in ARMS if a[6] is not None} | {"u2_four_mandate"})
    unions = {}
    for sel in needed:
        u, _details = MOD.prepare_h1_union(targets, OFFICIAL, PANEL, max_rss_gib=40.0, sources=sources_of(MOD, SEL[sel]))
        unions[sel] = u
        print(f"[{held}] union {sel:16s} {int(np.sum(np.asarray(u.target_condition_count) > 0)):3d}/{len(targets)} carried  ({time.time()-t0:.0f}s)", flush=True)
    # the two spellings of the four-source second union: a per-gene mean, so only the float accumulation order differs
    d_order = float(np.abs(np.asarray(unions["u2_four"].mean_effect) - np.asarray(unions["u2_four_mandate"].mean_effect)).max())
    same_votes = bool(np.array_equal(np.asarray(unions["u2_four"].gene_votes), np.asarray(unions["u2_four_mandate"].gene_votes)))
    print(f"[{held}] u2_four (248 order) vs u2_four_mandate (mandate order): max|d mean_effect| {d_order:.3e}, votes equal {same_votes}", flush=True)
    variants = {}

    def variant_of(sel, variant):
        key = (sel, variant)
        if key not in variants:
            u = unions[sel]
            if variant is None:
                variants[key] = u
            elif variant.startswith("SHUF"):
                variants[key] = shuffled(u, int(variant[4:]))
            elif variant == "CD4ROWS":
                variants[key] = restrict_rows(cd4, u)
            else:
                raise ValueError(variant)
        return variants[key]

    op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production" else PROD_OPS / "whitener_union_pruned.npz")
    op = Operator(op_path)
    wh = op.apply(pre * ALPHA, fold["wcols"], GAMMA, LAM)
    del op
    gc.collect()
    ctx = challenge_contexts()
    bA, gA = ctx["A"]
    sc = fold["scored"]
    Ts = truth_h[:, sc]
    cpm = 1e6 * np.expm1(fold["b_ctrl"][sc]) / float(np.expm1(fold["b_ctrl"]).sum())
    LR = to_lfc(truth_h, fold["b_ctrl"], gA)[:, sc]
    pools = {
        "wald89": topn(np.abs(LR) * np.sqrt(np.maximum(cpm, 0))[None, :], POOL_N),
        "lfc89": topn(np.abs(LR), POOL_N),
        "h89": topn(np.abs(truth_h[:, sc]), POOL_N),
    }
    print(f"[{held}] staged {len(targets)} targets ({time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB)", flush=True)

    rows: dict[str, dict] = {}
    arrays: dict[str, np.ndarray] = {}
    for name, k, resf, g, sel, g2, sel2, variant in ARMS:
        t1 = time.time()
        u2 = None if sel2 is None else variant_of(sel2, variant)
        final, info = compose(MOD, wh, cd4, fold["wcols"], fold["panel_positions"], k, resf, unions.get(sel), g, u2, g2)
        final[:, fold["panel_positions"]] = 0.0
        s, arr = score_arm(final, fold, Ts, ctx, LR, pools)
        rows[name] = dict(arm=name, k=k, f=resf, g=g, selection=sel, g2=g2, selection2=sel2, variant=variant,
                          **info, **s, seconds=round(time.time() - t1, 1))
        for kk, v in arr.items():
            arrays[f"{name}|{kk}"] = v
        print(f"[{held}] {name:34s} agg {s['cos_agg']:.6f} macro {s['a23_cos']:.6f} pds {s['pds']:.6f} pdsA {s['pds_A']:.6f} "
              f"prec {s['top89_precision']:.5f} fc {info['cd4_fraction']:.5f} "
              f"fu {'-' if info['union_fraction'] is None else round(info['union_fraction'], 5)} "
              f"fu2 {'-' if info['union2_fraction'] is None else round(info['union2_fraction'], 5)} "
              f"c2 {'-' if info['union2_carrying'] is None else info['union2_carrying']} ({rows[name]['seconds']}s)", flush=True)
        del final
        gc.collect()

    # ---------------------------------------------------------------- known-answer rows, this fold
    ka = []

    def check(label, got, want, tol=0.0):
        d = abs(float(got) - float(want))
        ok = d <= tol
        ka.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, ok=bool(ok)))
        if not ok or tol > 0:
            print(f"  {'PASS' if ok else 'FAIL'} {label:80s} got {got:.12g} want {want:.12g} |d| {d:.1e}", flush=True)

    print(f"\n[{held}] KNOWN-ANSWER ROWS", flush=True)
    # R1 (per fold): 258's rows and per-target arrays for the seven shared arms, at 0.0e+00
    f258 = json.loads((E258 / f"work/fold-{held}.json").read_text())["rows"]
    z258 = np.load(E258 / f"work/fold-{held}.npz")
    n_scalar = n_arr = 0
    for name in SHARED_258:
        theirs, mine = f258[name], rows[name]
        for key, want in theirs.items():
            if key in ("arm", "selection", "seconds", "k", "f", "g") or want is None or isinstance(want, str):
                continue
            check(f"258 fold row {name}.{key}", mine[key], want, 0.0)
            n_scalar += 1
        for kk in ("pds_pt", "cos_pt") + tuple(f"{p}_{pk}" for pk in POOLS for p in ("U", "M", "D")):
            same = bool(np.array_equal(z258[f"{name}|{kk}"], arrays[f"{name}|{kk}"]))
            ka.append(dict(label=f"258 fold array {name}|{kk} bit-identical", ok=same))
            n_arr += 1
            if not same:
                print(f"  FAIL 258 fold array {name}|{kk} differs (max|d| {np.abs(z258[f'{name}|{kk}'] - arrays[f'{name}|{kk}']).max():.3e})", flush=True)
    n_ok_258 = sum(r["ok"] for r in ka)
    print(f"[{held}] R1 per-fold: {n_ok_258}/{len(ka)} pass ({n_scalar} scalars at 0.0e+00, {n_arr} arrays bit-identical)", flush=True)
    # R3 (per fold): 248's union2 rows through my union2 path
    r248 = json.loads((E248 / f"work/repro248-{held}.json").read_text())
    for mine, (arm248, g248) in R3_ROWS.items():
        r = next(x for x in r248["arms"] if x["arm"] == arm248 and x["variant"] == "" and abs(x["g"] - g248) < 1e-12)
        m = rows[mine]
        for key_m, key_t in (("pds", "arm_pds"), ("cos_agg", "arm_cos_agg"), ("a23_cos", "arm_cos"), ("a15_sgn_w", "arm_a15_sgn_w"),
                             ("union2_fraction", "fraction"), ("preserve_scale", "preserve_scale"),
                             ("cd4_fraction", "frac_cd4"), ("union_fraction", "frac_h1"),
                             ("median_row_norm_final", "arm_median_row_norm_final")):
            check(f"248 repro248 {arm248}@{g248:g}.{key_t} == {mine}.{key_m}", m[key_m], r[key_t], 1e-9)
        check(f"248 repro248 {arm248}@{g248:g}.carrying_rows == {mine}", m["union2_carrying"], r["carrying_rows"], 0.0)
        # and 248's base (KA-A == 218's) is my h1union row
        for key_m, key_t in (("pds", "base_pds"), ("cos_agg", "base_cos_agg"), ("a23_cos", "base_cos")):
            check(f"248 repro248 base.{key_t} == {OLD}.{key_m}", rows[OLD][key_m], r[key_t], 1e-9)
    n_ok = sum(r["ok"] for r in ka)
    print(f"[{held}] known-answer rows {n_ok}/{len(ka)} pass", flush=True)

    out = dict(fold=held, n_targets=len(targets), targets=list(targets), builder_sha256=MOD.file_sha256(builder_path),
               rows=rows, known_answer=ka, n_ka_258=len(ka) - 0,
               union_carried={sel: int(np.sum(np.asarray(u.target_condition_count) > 0)) for sel, u in unions.items()},
               u2_four_order=dict(max_abs_diff_mean_effect=d_order, votes_equal=same_votes),
               seconds=round(time.time() - t0, 1), peak_rss_gib=peak_gib())
    (WORK / f"fold-{held}.json").write_text(json.dumps(out, indent=1, default=float))
    np.savez_compressed(WORK / f"fold-{held}.npz", **arrays)
    print(f"[{held}] {len(rows)} arms in {out['seconds']}s, peak RSS {out['peak_rss_gib']:.2f} GiB; wrote {WORK/f'fold-{held}.json'}", flush=True)
    if n_ok != len(ka):
        raise SystemExit(1)


# ------------------------------------------------------------------------------ production phase
def terms_from_final(final, targets, b, gate):
    """211/243's per-context (A_c, B_c, Kc_c, cap_c) on the EFFECT-LEVEL render (258's, verbatim)."""
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
    same_order = tuple(targets) == tuple(PANEL)
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
        print(f"  {'PASS' if ok else 'FAIL'} {label:90s} got {got:.12g} want {want:.12g} |d| {d:.1e}", flush=True)

    print("\nKNOWN-ANSWER ROWS (production frame)", flush=True)
    check("recovered base: median row norm == stage-off whitened_base_median_row_norm", float(np.median(np.linalg.norm(wh, axis=1))), base_norm, 1e-9)
    unions = {}
    for sel in sorted({a[4] for a in ARMS if a[4] is not None} | {a[6] for a in ARMS if a[6] is not None}):
        unions[sel], _ = MOD.prepare_h1_union(targets, OFFICIAL, PANEL, max_rss_gib=40.0, sources=sources_of(MOD, SEL[sel]))
    d_u = float(np.abs(np.asarray(unions["pair"].mean_effect) - np.asarray(union_pair_saved.mean_effect)).max())
    check("prepare_h1_union(PAIR) on the production panel == 249's union-off.npz (max|d|)", d_u, 0.0, 0.0)
    f_off, i_off = compose(MOD, wh, cd4, wcols, panel_positions, None, 0.34, unions["pair"], 0.28)
    check("re-composed h1union == 249 effects-off.npy (abeb0749), max|d|", float(np.abs(f_off - final_off).max()), 0.0, 1e-9)
    final_svd5 = np.load(stage / "effects-svd5.npy")
    f_5, i_5 = compose(MOD, wh, cd4, wcols, panel_positions, 5, 0.34, unions["pair"], 0.28)
    check("re-composed cd4svd5 == 249 effects-svd5.npy (a28d5be2), max|d|", float(np.abs(f_5 - final_svd5).max()), 0.0, 1e-9)
    check("re-composed cd4svd5 cd4 fraction == stage-svd5", i_5["cd4_fraction"], float(aud_svd5["audit"]["cd4"]["fraction"]), 1e-12)
    check("re-composed cd4svd5 preserve == stage-svd5", i_5["preserve_scale"], float(aud_svd5["audit"]["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"]), 1e-12)
    check("effects-svd5.npy array sha == a28d5be2... (249 record)", 1.0 if MOD.array_sha256(final_svd5).startswith("a28d5be2") else 0.0, 1.0, 0.0)
    final_arm6 = np.load(E229 / "work/effects-stack6g035.npy")
    f_6, i_6 = compose(MOD, wh, cd4, wcols, panel_positions, None, 0.34, unions["stack6"], 0.35)
    check("re-composed ARM6 == 229 effects-stack6g035.npy (c5721825), max|d|", float(np.abs(f_6 - final_arm6).max()), 0.0, 1e-9)
    check("effects-stack6g035.npy array sha == c5721825... (229 record)", 1.0 if MOD.array_sha256(final_arm6).startswith("c5721825") else 0.0, 1.0, 0.0)
    # R3 at the production level: 248's staged second-residual arrays through the union2 path
    for tag, sel2, sha in (("cum1", "u2_dld1", "6225a582"), ("cum4", "u2_four", "2f9eba51")):
        f248 = np.load(E248 / f"work/stage/effects-{tag}.npy")
        a248 = json.loads((E248 / f"work/stage/stage-{tag}.json").read_text())
        f_x, i_x = compose(MOD, wh, cd4, wcols, panel_positions, None, 0.34, unions["pair"], 0.28, unions[sel2], 0.35)
        check(f"re-composed h1union + union2 {tag}@0.35 == 248 effects-{tag}.npy ({sha}), max|d|", float(np.abs(f_x - f248).max()), 0.0, 1e-9)
        check(f"effects-{tag}.npy array sha == {sha}... (248 record)", 1.0 if MOD.array_sha256(f248).startswith(sha) else 0.0, 1.0, 0.0)
        fr = a248["audit"]["union2"]["residual_frame"]
        check(f"re-composed {tag} union2 solved fraction == 248 stage-{tag}", i_x["union2_fraction"], float(fr["solved_fraction"]), 1e-12)
        check(f"re-composed {tag} union2 carrying rows == 248 stage-{tag}", i_x["union2_carrying"], int(fr["carrying_rows"]), 0.0)
        check(f"re-composed {tag} preserve == 248 stage-{tag}", i_x["preserve_scale"], float(a248["audit"]["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"]), 1e-12)
        check(f"re-composed {tag} h1 fraction == 248 stage-{tag}", i_x["union_fraction"], float(a248["audit"]["h1_union"]["residual_frame"]["solved_fraction"]), 1e-12)
        del f_x, f248
    final_svd10 = np.load(E255 / "work/stage/effects-svd10.npy")
    n_ok = sum(r["ok"] for r in ka)
    print(f"known-answer rows {n_ok}/{len(ka)} pass  ({time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB)", flush=True)
    if n_ok != len(ka):
        (WORK / "prod259.json").write_text(json.dumps(dict(known_answer=ka, failed=True), indent=1))
        raise SystemExit(1)

    ctx = challenge_contexts()
    b = {c: ctx[c][0] for c in CTX}
    gate = {c: ctx[c][1] for c in CTX}
    bA, gA = ctx["A"]
    variants = {}

    def variant_of(sel, variant):
        key = (sel, variant)
        if key not in variants:
            uu = unions[sel]
            variants[key] = uu if variant is None else (shuffled(uu, int(variant[4:])) if variant.startswith("SHUF") else restrict_rows(cd4, uu))
        return variants[key]

    terms: dict[str, dict] = {}
    info_all: dict[str, dict] = {}
    finals = {"h1union g=0.28 (file)": final_off, "cd4svd5 (file)": final_svd5, "ARM6 (file)": final_arm6, "cd4svd10 (255 file)": final_svd10}
    for name, k, resf, g, sel, g2, sel2, variant in ARMS:
        if name == "a15 f=0.10 SHIPPED":
            continue
        u2 = None if sel2 is None else variant_of(sel2, variant)
        final, info = compose(MOD, wh, cd4, wcols, panel_positions, k, resf, unions.get(sel), g, u2, g2)
        finals[name] = final
        info_all[name] = info
    for name, final in finals.items():
        t1 = time.time()
        per = terms_from_final(final, targets, b, gate)
        H = deliver(final, bA, gA, None)
        terms[name] = dict(per=per, delivered_norm_A=float(np.median(np.linalg.norm(H[:, wcols], axis=1))),
                           median_row_norm_final=float(np.median(np.linalg.norm(final, axis=1))))
        print(f"  {name:36s} A {per['A']['A']:.6f}/{per['B']['A']:.6f}/{per['C']['A']:.6f}  dlvA {terms[name]['delivered_norm_A']:.5f} ({time.time()-t1:.0f}s)", flush=True)
        del H
    inv243 = json.loads((E243 / "work/invert243.json").read_text())
    amp249 = json.loads((E249 / "work/amp249.json").read_text())
    amp255 = json.loads((E255 / "work/amp255.json").read_text())
    arm6_243 = json.loads((E243 / "work/arm6_243.json").read_text())
    pkg = dict(h1union=inv243["arms"]["h1union"]["per"], cd4svd5=amp249["per_arm"], ARM6=arm6_243["per"], cd4svd10=amp255["per_arm"])
    calib = []
    for arm, base, fa, fb in (("ARM6", "h1union", "ARM6 (file)", "h1union g=0.28 (file)"),
                              ("cd4svd5", "h1union", "cd4svd5 (file)", "h1union g=0.28 (file)"),
                              ("cd4svd10", "cd4svd5", "cd4svd10 (255 file)", "cd4svd5 (file)")):
        r_pkg = amp_ratio(pkg[arm], pkg[base])
        r_eff = amp_ratio(terms[fa]["per"], terms[fb]["per"])
        r_med = terms[fa]["delivered_norm_A"] / terms[fb]["delivered_norm_A"]
        calib.append(dict(arm=arm, base=base, package_pct=100 * (r_pkg - 1), effect_level_pct=100 * (r_eff - 1),
                          median_norm_A_pct=100 * (r_med - 1), ratio_pkg_over_eff=(r_pkg - 1) / (r_eff - 1) if r_eff != 1 else None))
        print(f"  amplitude calibration {arm:8s} vs {base:8s}: package {100*(r_pkg-1):+.4f}%  effect-level A_c {100*(r_eff-1):+.4f}%  216 median-norm A {100*(r_med-1):+.4f}%", flush=True)
    out = dict(known_answer=ka, targets_in_panel_order=same_order, calibration=calib, terms=terms, compose_info=info_all,
               package_terms=pkg, seconds=round(time.time() - t0, 1), peak_rss_gib=peak_gib())
    (WORK / "prod259.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"wrote {WORK/'prod259.json'}  ({out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB)", flush=True)


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
    """258's `price()`: the six-member legs of one cell on one factor set (the pds leg on ONE functional)."""
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
    """COHERENT scenarios: per leg, the end of its measured band chosen by the SIGN of that leg for this cell.
    A cosine gain feeds `mse` and DE-three the same way, so one k_cos choice is coherent for both."""
    gain_cos, gain_pds, gain_nmae = rel_agg > 0, d_pds > 0, d_nmae < 0
    if which == "POINT":
        return dict(POINT)
    if which == "PESSIMISTIC":
        return dict(k_cos=K_COS_BAND[0] if gain_cos else K_COS_BAND[1], k_pds=K_PDS_BAND[0] if gain_pds else K_PDS_BAND[1],
                    de3=DE3_BAND[0] if gain_cos else DE3_BAND[1], nmae_x=0.0 if gain_nmae else NMAE_LOSS_BAND[1])
    if which == "OPTIMISTIC":
        return dict(k_cos=K_COS_BAND[1] if gain_cos else K_COS_BAND[0], k_pds=K_PDS_BAND[1] if gain_pds else K_PDS_BAND[0],
                    de3=DE3_BAND[1] if gain_cos else DE3_BAND[0], nmae_x=1.0)
    raise ValueError(which)


def run_price() -> None:
    lines = []

    def P(t=""):
        print(t, flush=True)
        lines.append(t)

    folds = {h: json.loads((WORK / f"fold-{h}.json").read_text()) for h in FOLDS}
    arrs = {h: np.load(WORK / f"fold-{h}.npz") for h in FOLDS}
    prod = json.loads((WORK / "prod259.json").read_text())
    p258 = json.loads((E258 / "work/price258.json").read_text())
    names = [a[0] for a in ARMS]
    ka_fold = {h: (sum(r["ok"] for r in folds[h]["known_answer"]), len(folds[h]["known_answer"])) for h in FOLDS}
    ka_prod = (sum(r["ok"] for r in prod["known_answer"]), len(prod["known_answer"]))
    P("=" * 130)
    P("259 -- the union WITHOUT igvf_cardiomyocyte on top of cd4svd5: construction A (single union replacing the pair) and")
    P("       construction B (the pair kept at 0.28 + a SECOND residual, 248's union2 path); three coherent scenarios")
    P("       (PLOLO 162, three cancer folds; kolf2 apart; nmae identity pooled over 4 folds as six229 does)")
    P("=" * 130)
    P("per-fold known-answer rows: " + ", ".join(f"{h} {a}/{b}" for h, (a, b) in ka_fold.items()) + f"; production-frame {ka_prod[0]}/{ka_prod[1]}")
    for h in FOLDS:
        P(f"  {h:14s} {folds[h]['n_targets']} targets, {folds[h]['seconds']}s, peak RSS {folds[h]['peak_rss_gib']:.2f} GiB; "
          f"u2_four order test max|d| {folds[h]['u2_four_order']['max_abs_diff_mean_effect']:.2e}, votes equal {folds[h]['u2_four_order']['votes_equal']}")
    P(f"  production frame: {prod['seconds']}s, peak RSS {prod['peak_rss_gib']:.2f} GiB")

    # ---- cross-fold aggregates of both functionals (258's convention) --------------------------
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
                a["U_sum"], a["M_sum"] = float(U.sum()), float(M.sum())
        agg[n] = a
    a15n = agg["a15 f=0.10 SHIPPED"]["delivered_norm_A_4f"]
    for n in names:
        agg[n]["delivered_norm_ratio"] = agg[n]["delivered_norm_A_4f"] / a15n

    def delta(a_, b_, key, folds_=CANCER, rel=False):
        vals = [agg[a_]["per_fold"][h][key] / agg[b_]["per_fold"][h][key] - 1 if rel else agg[a_]["per_fold"][h][key] - agg[b_]["per_fold"][h][key]
                for h in folds_]
        return float(np.mean(vals)), [float(v) for v in vals]

    # ---- amplitude ------------------------------------------------------------------------------
    ratios = [c["ratio_pkg_over_eff"] for c in prod["calibration"] if c["ratio_pkg_over_eff"] is not None]
    amp_cal = float(np.mean(ratios))
    T = prod["terms"]

    def amp_vs(name, base_name):
        r_eff = amp_ratio(T[name]["per"], T[base_name]["per"])
        return r_eff, 1.0 + amp_cal * (r_eff - 1.0)

    # ---- R1: cross-fold, against 258's price258.json ---------------------------------------------
    ka = []

    def check(label, got, want, tol=0.0):
        d = abs(float(got) - float(want))
        ok = d <= tol
        ka.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, ok=bool(ok)))
        P(f"  {'PASS' if ok else 'FAIL'} {label:86s} got {got:+.10f} want {want:+.10f} |d| {d:.1e}")

    P("\nR1 -- 258's rows reproduce from this runner (cross-fold aggregates at 0.0e+00; per-fold rows and arrays in fold-<held>.json)")
    for n in SHARED_258:
        for key in ("cos_agg_3c", "a23_cos_3c", "pds_3c", "pds_A_3c", "nmae_wald89", "nmae_wald89_3c", "nmae_lfc89", "nmae_h89",
                    "top89_precision_4f", "delivered_norm_ratio", "U_sum", "M_sum"):
            check(f"258 aggregate {n}.{key}", agg[n][key], p258["aggregates"][n][key], 0.0)
    B_h = board_entry("iUziqE1gewfcPKzYUDaB")      # h1union
    B_6 = board_entry("ydAemVdWmsTx9xhDWcou")      # ARM6
    B_5 = board_entry("Xv0ayvuoHpsKYWv2iu42")      # cd4svd5
    pkg = prod["package_terms"]
    cos_h = invert(pkg["h1union"], B_h["raw"]["mse"])
    cos_5 = invert(pkg["cd4svd5"], B_5["raw"]["mse"])
    cos_6 = invert(pkg["ARM6"], B_6["raw"]["mse"])
    raw_5 = B_5["raw"]["mse"]
    base_per = pkg["cd4svd5"]
    check("amplitude calibration factor == 258", amp_cal, p258["factors"]["amp_calibration"], 1e-12)
    check("measured board cosine cd4svd5 == 258", cos_5, p258["factors"]["cos_cd4svd5"], 1e-12)
    # 258's own conservative / point totals for the three screen KA cells, through my legs()
    for n in ("stack6 g=0.35", "stack3 g=0.35", "pair+igvf g=0.28"):
        rel, _ = delta(n, INC, "cos_agg", rel=True)
        dpA, _ = delta(n, INC, "pds_A")
        dpABC, _ = delta(n, INC, "pds")
        d_nmae = agg[n]["nmae_wald89"] - agg[INC]["nmae_wald89"]
        _, r_cal = amp_vs(n, INC)
        c258 = p258["cells"][n]
        check(f"258 cell {n}.amp_calibrated_pct", 100 * (r_cal - 1), c258["amp_calibrated_pct"], 1e-9)
        cons = legs(rel, dpA, d_nmae, r_cal, base_per, cos_5, raw_5, 0.96, 0.7, p258["factors"]["de3_rate_ARM6"], 1.5)
        point = legs(rel, dpA, d_nmae, r_cal, base_per, cos_5, raw_5, 0.96, p258["factors"]["k_pds_ARM6_A"], 0.24, 1.0)
        check(f"258 cell {n}.conservative.total_A (k_pds 0.7 frame A, nmae x1.5)", cons["total"], c258["conservative"]["total_A"], 1e-9)
        check(f"258 cell {n}.point.total_A (k_pds 3.79 frame A, DE3 0.24)", point["total"], c258["point"]["total_A"], 1e-9)
        check(f"258 cell {n}.conc[3].retained", _retention(arrs, n, INC, 3, dpABC)["retained"], c258["conc"]["3"]["retained"], 1e-9)
    r1_n = sum(r["ok"] for r in ka)
    P(f"R1: {r1_n}/{len(ka)} cross-fold rows pass; per fold: " + ", ".join(f"{h} {a}/{b}" for h, (a, b) in ka_fold.items()))
    r1_pass = r1_n == len(ka) and all(a == b for a, b in ka_fold.values())

    # ---- R3: 248's second-residual row, cross-fold ------------------------------------------------
    P("\nR3 -- 248's second residual through my union2 path (dld1_lowmoi @0.35 on the h1union base; repro248-summary.json)")
    s248 = json.loads((E248 / "work/repro248-summary.json").read_text())
    ka3 = []

    def check3(label, got, want, tol):
        d = abs(float(got) - float(want))
        ok = d <= tol
        ka3.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, ok=bool(ok)))
        P(f"  {'PASS' if ok else 'FAIL'} {label:86s} got {got:+.10f} want {want:+.10f} |d| {d:.1e}")

    for mine, key248 in (("248 CUM1 dld1 g2=0.35 on h1union", "CUM1@0.35"), ("248 CUM1 dld1 g2=0.28 on h1union", "CUM1@0.28"),
                         ("248 CUM4 four g2=0.35 on h1union", "CUM4@0.35")):
        ref = s248["arms"][key248]
        rel_m, _ = delta(mine, OLD, "a23_cos", rel=True)
        rel_a, _ = delta(mine, OLD, "cos_agg", rel=True)
        dp, _ = delta(mine, OLD, "pds")
        check3(f"248 {key248} rel macro cosine (3 cancer) == {mine}", rel_m, ref["rel_macro"], 1e-9)
        check3(f"248 {key248} rel cos_agg (3 cancer)", rel_a, ref["rel_agg"], 1e-9)
        check3(f"248 {key248} d pds (3 cancer)", dp, ref["d_pds"], 1e-9)
        check3(f"248 {key248} kolf2 d pds", delta(mine, OLD, "pds", ("kolf2",))[0], ref["kolf2"]["d_pds"], 1e-9)
    conc225 = json.loads((E225 / "work/conc225.json").read_text())["summary"]["dld1_lowmoi"]
    dp_cum1, _ = delta("248 CUM1 dld1 g2=0.35 on h1union", OLD, "pds")
    ret_cum1 = _retention(arrs, "248 CUM1 dld1 g2=0.35 on h1union", OLD, 3, dp_cum1)
    check3("248 KA-E: leave-top-3-out pds retention of dld1@0.35 == 225 conc225 (25.2 %)", ret_cum1["retained"], conc225["k"]["3"]["pds_retained"], 1e-6)
    check3("248 KA-E: top-3 pds share == 225 conc225 (81.9 %)", ret_cum1["share"], conc225["k"]["3"]["top_share_pds"], 1e-6)
    P(f"  the mandate's digits: rel macro +1.405 %, d pds +0.00256 3/3 folds, leave-top-3-out 25.2 % (got {100*rel_m:+.3f} % / {dp:+.5f} / {100*ret_cum1['retained']:.1f} %)")
    r3_pass = all(r["ok"] for r in ka3)
    P(f"R3: {sum(r['ok'] for r in ka3)}/{len(ka3)} pass")

    # ---- R2: the two board retrodictions ------------------------------------------------------------
    P("\nR2 -- retrodictions of the two 2026-09-05 board points from THIS runner's rows")
    rel6, _ = delta(ARM6, OLD, "cos_agg", rel=True)
    dp6_ABC, _ = delta(ARM6, OLD, "pds")
    dp6_A, _ = delta(ARM6, OLD, "pds_A")
    rel5, _ = delta(INC, OLD, "cos_agg", rel=True)
    dp5_ABC, _ = delta(INC, OLD, "pds")
    dp5_A, _ = delta(INC, OLD, "pds_A")
    d_nmae6 = agg[ARM6]["nmae_wald89"] - agg[OLD]["nmae_wald89"]
    d_nmae5 = agg[INC]["nmae_wald89"] - agg[OLD]["nmae_wald89"]
    k_cos6 = (cos_6 - cos_h) / (rel6 * cos_h)
    k_cos5 = (cos_5 - cos_h) / (rel5 * cos_h)
    k_pds6_ABC = (B_6["raw"]["pds"] - B_h["raw"]["pds"]) / dp6_ABC
    k_pds6_A = (B_6["raw"]["pds"] - B_h["raw"]["pds"]) / dp6_A
    k_pds5_ABC = (B_5["raw"]["pds"] - B_h["raw"]["pds"]) / dp5_ABC
    de3_6 = sum(B_6["scaled"][m] - B_h["scaled"][m] for m in ("fid", "reach", "jac")) / 6
    de3_5 = sum(B_5["scaled"][m] - B_h["scaled"][m] for m in ("fid", "reach", "jac")) / 6
    de3_rate6 = de3_6 / (cos_6 - cos_h)
    board6 = {m: (B_6["scaled"][m] - B_h["scaled"][m]) / 6 for m in ("pds", "mse", "nmae", "fid", "reach", "jac")}
    board5 = {m: (B_5["scaled"][m] - B_h["scaled"][m]) / 6 for m in ("pds", "mse", "nmae", "fid", "reach", "jac")}
    P(f"  board cosines (197's identity at own terms): h1union {cos_h:.6f}  ARM6 {cos_6:.6f}  cd4svd5 {cos_5:.6f}")
    P(f"  realised factors from these rows: ARM6 k_cos {k_cos6:.3f}, k_pds {k_pds6_ABC:.2f} (ABC) / {k_pds6_A:.2f} (A), DE-three {de3_rate6:.4f}/unit d cos, "
      f"nmae identity {-d_nmae6/NMAE_SPAN/6:+.5f} vs board {board6['nmae']:+.5f} (ratio {board6['nmae']/(-d_nmae6/NMAE_SPAN/6):.3f}); "
      f"cd4svd5 k_cos {k_cos5:.3f}, k_pds {k_pds5_ABC:.2f} (ABC)")
    r6_pkg = amp_ratio(pkg["ARM6"], pkg["h1union"])
    r5_pkg = amp_ratio(pkg["cd4svd5"], pkg["h1union"])
    # (a) BR237's ARM6 by member at the factor set 258 s0 found: k_cos 0.96, k_pds ARM6-realised, identity x1.5, amplitude leg
    P("  (a) BR237's ARM6 (stack6 g=0.35 on h1union) by member -- factor set: k_cos 0.96, k_pds 2.48 (ABC) [3.79 on frame A], "
      "DE-three at the ARM6-realised rate and at 0.24, nmae identity x1.5, the package amplitude leg (-1.24 %)")
    r2a = {}
    for de_label, de in (("DE3 0.0425", 0.0425), ("DE3 0.24", 0.24)):
        for fr_label, dp6, kp in (("ABC", dp6_ABC, 2.48), ("A", dp6_A, 3.79)):
            p = legs(rel6, dp6, d_nmae6, r6_pkg, pkg["h1union"], cos_h, B_h["raw"]["mse"], 0.96, kp, de, 1.5)
            r2a[f"{de_label}|{fr_label}"] = p
            P(f"      {de_label:10s} frame {fr_label:3s}: pds {p['o_pds']:+.5f} (board {board6['pds']:+.5f})  mse {p['o_mse']:+.5f} (board {board6['mse']:+.5f})  "
              f"DE3 {p['o_de3']:+.5f} (board {de3_6:+.5f})  nmae {p['o_nmae']:+.5f} (board {board6['nmae']:+.5f})  TOTAL {p['total']:+.5f} (board {B_6['overall']-B_h['overall']:+.5f})")
    p_a = r2a["DE3 0.0425|ABC"]
    r2a_pass = bool(abs(p_a["total"] - (B_6["overall"] - B_h["overall"])) <= 0.0010 and abs(p_a["o_pds"] - board6["pds"]) <= 0.0008
                    and abs(p_a["o_mse"] - board6["mse"]) <= 0.0003 and abs(p_a["o_nmae"] - board6["nmae"]) <= 0.0025)
    # (b) BR249: cd4svd5 vs h1union at k_cos 3.76 / k_pds 1.22
    p_b = legs(rel5, dp5_ABC, d_nmae5, r5_pkg, pkg["h1union"], cos_h, B_h["raw"]["mse"], 3.76, 1.22, de3_5 / (cos_5 - cos_h), 1.0)
    p_b_de = legs(rel5, dp5_ABC, d_nmae5, r5_pkg, pkg["h1union"], cos_h, B_h["raw"]["mse"], 3.76, 1.22, 0.24, 1.0)
    P(f"  (b) BR249's cd4svd5 vs h1union at k_cos 3.76 / k_pds 1.22 (ABC), package amplitude (-0.41 %), identity x1.0: "
      f"pds {p_b['o_pds']:+.5f} (board {board5['pds']:+.5f})  mse {p_b['o_mse']:+.5f} (board {board5['mse']:+.5f})  DE3 {p_b['o_de3']:+.5f} / {p_b_de['o_de3']:+.5f} at 0.24 (board {de3_5:+.5f})  "
      f"nmae {p_b['o_nmae']:+.5f} (board {board5['nmae']:+.5f})  TOTAL {p_b['total']:+.5f} / {p_b_de['total']:+.5f} (board {B_5['overall']-B_h['overall']:+.5f}; band +0.0043 +- 0.0008)")
    r2b_pass = bool(abs(p_b["total"] - 0.0043) <= 0.0008)
    P(f"  R2: (a) {'PASS' if r2a_pass else 'MISS'} within BR237's bands (total +-0.0010, pds +-0.0008, mse +-0.0003, nmae +-0.0025); (b) {'PASS' if r2b_pass else 'MISS'}")
    P("  honesty: k_cos and k_pds in (a) and (b) are the values those two board points themselves fixed, so those two legs pass by construction; "
      "the non-circular content is the identity at x1.5, the amplitude leg and the DE-three rate (258 s2).")

    # ---- THE SCREEN --------------------------------------------------------------------------------
    P("\n" + "=" * 130)
    P("THE SCREEN -- every cell vs cd4svd5 (the incumbent), three cancer folds, THREE COHERENT SCENARIOS")
    P("=" * 130)
    chk = raw_for(base_per, cos_5)
    P(f"  197's identity at cd4svd5's own terms reproduces its board raw mse: {chk:.7f} vs {raw_5:.7f} (|d| {abs(chk-raw_5):.1e}); cos_cd4svd5 {cos_5:.6f}")
    P(f"  bands: k_pds {K_PDS_BAND} (ABC functional), k_cos {K_COS_BAND}, DE-three {DE3_BAND} per unit board d cos, nmae identity x{NMAE_LOSS_BAND} on a loss / "
      f"x0 (pess.) or x1.0 (opt.) credit on a gain; amplitude leg (calibrated effect-level, factor {amp_cal:.3f}) always charged; POINT = {POINT}")
    P("  PESSIMISTIC picks, per leg, the band end that LOWERS the total for that cell's sign; OPTIMISTIC the other end; the frame-A pds is reported, never priced")
    cells = {}
    screen = [a for a in ARMS if a[1] == 5 and a[0] != INC]
    for name, k, resf, g, sel, g2, sel2, variant in screen:
        rel, rel_pf = delta(name, INC, "cos_agg", rel=True)
        relm, _ = delta(name, INC, "a23_cos", rel=True)
        dpA, dpA_pf = delta(name, INC, "pds_A")
        dpABC, dpABC_pf = delta(name, INC, "pds")
        dqw, _ = delta(name, INC, "a15_sgn_w")
        d_nmae = agg[name]["nmae_wald89"] - agg[INC]["nmae_wald89"]
        d_nmae3 = agg[name]["nmae_wald89_3c"] - agg[INC]["nmae_wald89_3c"]
        dU = agg[name]["U_sum"] - agg[INC]["U_sum"]
        dM = agg[name]["M_sum"] - agg[INC]["M_sum"]
        r_eff, r_cal = amp_vs(name, INC)
        kolf = dict(rel_agg=delta(name, INC, "cos_agg", ("kolf2",), rel=True)[0], d_pds=delta(name, INC, "pds", ("kolf2",))[0])
        conc = {kk: _retention(arrs, name, INC, kk, dpABC, targets={h: folds[h]["targets"] for h in CANCER}) for kk in (1, 3, 5)}
        sc = {}
        for which in ("PESSIMISTIC", "POINT", "OPTIMISTIC"):
            f = scenario_factors(which, rel, dpABC, d_nmae)
            sc[which] = legs(rel, dpABC, d_nmae, r_cal, base_per, cos_5, raw_5, f["k_cos"], f["k_pds"], f["de3"], f["nmae_x"])
        cons258 = legs(rel, dpA, d_nmae, r_cal, base_per, cos_5, raw_5, 0.96, 0.7, DE3_BAND[0], 1.5)
        point_A = legs(rel, dpA, d_nmae, r_cal, base_per, cos_5, raw_5, 0.96, K_PDS_A_POINT, 0.24, 1.0)
        sig = abs(dpABC) / SIGMA_PDS
        cells[name] = dict(arm=name, g=g, g2=g2, selection=sel, selection2=sel2, variant=variant,
                           sources1=list(SEL[sel]) if sel else None, sources2=list(SEL[sel2]) if sel2 else None,
                           carrying1={h: folds[h]["rows"][name]["union_carrying"] for h in FOLDS},
                           carrying2={h: folds[h]["rows"][name]["union2_carrying"] for h in FOLDS},
                           shared2={h: folds[h]["rows"][name]["union2_shared_targets"] for h in FOLDS},
                           prod_info=prod["compose_info"].get(name),
                           rel_cos_agg=rel, rel_cos_agg_per_fold=rel_pf, rel_macro=relm,
                           d_pds_A=dpA, d_pds_A_per_fold=dpA_pf, d_pds_ABC=dpABC, d_pds_ABC_per_fold=dpABC_pf,
                           pds_pos_folds=int(sum(v > 0 for v in dpABC_pf)), pds_sigma=sig, pds_fold_se=float(np.std(dpABC_pf, ddof=1) / np.sqrt(3)),
                           d_qw=dqw, d_nmae_wald89=d_nmae, d_nmae_wald89_3c=d_nmae3, dU=dU, dM=dM, e=(dM / dU if dU != 0 else float("nan")),
                           d_precision89=agg[name]["top89_precision_4f"] - agg[INC]["top89_precision_4f"],
                           amp_effect_level_pct=100 * (r_eff - 1), amp_calibrated_pct=100 * (r_cal - 1),
                           kolf2=kolf, conc=conc, scenarios=sc, cons258=cons258, point_A=point_A,
                           clears_pessimistic=bool(sc["PESSIMISTIC"]["total"] >= GATE),
                           retention_ok=bool(conc[3]["retained"] >= 0.5), pds_ok=bool(all(v > 0 for v in dpABC_pf) and sig >= 2.0),
                           cleared=bool(sc["PESSIMISTIC"]["total"] >= GATE and conc[3]["retained"] >= 0.5 and all(v > 0 for v in dpABC_pf) and sig >= 2.0))

    def table(title, subset):
        P(f"\n{title}")
        P(f"  {'cell':26s}{'c1':>4s}{'c2':>4s}{'agg%':>8s}{'mac%':>8s}{'dpds ABC':>10s}{'+':>2s}{'sig':>5s}{'dpds A':>9s}{'dnmae':>9s}{'e':>7s}{'amp%':>7s}"
          f"{'|':>2s}{'PESS':>9s}{'POINT':>9s}{'OPT':>9s}{'|':>2s}{'C258':>8s}{'ret3':>7s}{'kolf agg%':>10s}{'kolf dpds':>10s}")
        for n in subset:
            c = cells[n]
            s = c["scenarios"]
            c1 = c["carrying1"]["orion_hek293t"]
            c2 = c["carrying2"]["orion_hek293t"]
            P(f"  {n:26s}{c1 if c1 is not None else 0:>4d}{c2 if c2 is not None else 0:>4d}{100*c['rel_cos_agg']:>+8.3f}{100*c['rel_macro']:>+8.3f}{c['d_pds_ABC']:>+10.5f}"
              f"{c['pds_pos_folds']:>2d}{c['pds_sigma']:>5.1f}{c['d_pds_A']:>+9.5f}{c['d_nmae_wald89']:>+9.5f}{c['e']:>7.3f}{c['amp_calibrated_pct']:>+7.2f}{'|':>2s}"
              f"{s['PESSIMISTIC']['total']:>+9.5f}{s['POINT']['total']:>+9.5f}{s['OPTIMISTIC']['total']:>+9.5f}{'|':>2s}{c['cons258']['total']:>+8.5f}"
              f"{100*c['conc'][3]['retained']:>6.0f}%{100*c['kolf2']['rel_agg']:>+10.3f}{c['kolf2']['d_pds']:>+10.5f}")

    P("\ncolumns: c1/c2 = union-1 / union-2 carrying rows on the hek293t fold (/294); agg%/mac% = rel cos_agg / macro cosine vs cd4svd5; dpds ABC = 246's functional "
      "(+ = folds > 0, sig = |dpds|/0.001494); dpds A = six229's frame-A functional; dnmae = exact identity (4 folds, wald89); e = dM/dU (break-even 0.5); "
      "amp% = calibrated delivered amplitude; PESS/POINT/OPT = the three coherent six-member totals; C258 = 258's conservative set (k_pds 0.7 on frame A, "
      "k_cos 0.96, DE3 0.0425, nmae x1.5) for continuity; ret3 = leave-top-3-out pds retention; kolf2 apart")
    A_cells = [f"A5 g={g:g}" for g in DOSES_A]
    B_cells = [f"B3 g2={g:g}" for g in DOSES_B]
    table("A. the single union WITHOUT igvf replacing the H1 pair (gse295214_h1,h1_2025,dld1_lowmoi,jiang_5line_noK562,myllia_a549 at g), on cd4svd5", A_cells)
    table("B. the pair kept at 0.28 + the SECOND residual (dld1_lowmoi,jiang_5line_noK562,myllia_a549 at g2), 248's union2 path, on cd4svd5", B_cells)
    table("(i) B with igvf back (union2 = dld1,jiang,igvf,myllia -- 248's CUM4 set/order) at g2", [f"B4+igvf g2={g:g}" for g in DOSES_CTRL])
    table("(ii) B with dld1_lowmoi alone (248's rider, re-priced on cd4svd5) at g2", [f"B1 dld1 g2={g:g}" for g in DOSES_CTRL])
    table("attribution inside B: jiang alone / myllia alone as union2", [f"B1 jiang g2={g:g}" for g in (0.28, 0.35)] + [f"B1 myllia g2={g:g}" for g in (0.28, 0.35)])
    table("controls for B (energy-matched): target-shuffled union2 (two seeds) and CD4's own information on the rows it carries (row-matched)",
          [f"B3 g2={g:g} {v}" for g in DOSES_CTRL for v in ("SHUF20726", "SHUF771103", "CD4ROWS")])
    table("258's cells re-run (the same numbers): stack6 g=0.35, stack3 g=0.35, pair+igvf g=0.28", ["stack6 g=0.35", "stack3 g=0.35", "pair+igvf g=0.28"])
    bestA = max(A_cells, key=lambda n: cells[n]["scenarios"]["PESSIMISTIC"]["total"])
    bestB = max(B_cells, key=lambda n: cells[n]["scenarios"]["PESSIMISTIC"]["total"])
    real_cells = [n for n in cells if cells[n]["variant"] is None]
    best_all = max(real_cells, key=lambda n: cells[n]["scenarios"]["PESSIMISTIC"]["total"])
    P(f"\n  best A cell by the PESSIMISTIC total: {bestA} ({cells[bestA]['scenarios']['PESSIMISTIC']['total']:+.5f}); best B cell: {bestB} "
      f"({cells[bestB]['scenarios']['PESSIMISTIC']['total']:+.5f}); best real cell overall: {best_all} ({cells[best_all]['scenarios']['PESSIMISTIC']['total']:+.5f})")

    # ---- the best cells in full --------------------------------------------------------------------
    def full(name):
        c = cells[name]
        P("\n" + "=" * 130)
        P(f"IN FULL: {name}  (union 1 {c['sources1']} at g {c['g']}; union 2 {c['sources2']} at g2 {c['g2']})")
        P("=" * 130)
        P(f"  PLOLO vs cd4svd5: rel cos_agg {100*c['rel_cos_agg']:+.4f}% (per fold {['%+.3f' % (100*v) for v in c['rel_cos_agg_per_fold']]}), rel macro {100*c['rel_macro']:+.4f}%")
        P(f"  d pds frames ABC {c['d_pds_ABC']:+.6f} (per fold {['%+.5f' % v for v in c['d_pds_ABC_per_fold']]}; {c['pds_pos_folds']}/3 positive; {c['pds_sigma']:.1f} LOO sigma; "
          f"fold se {c['pds_fold_se']:.6f}); frame A {c['d_pds_A']:+.6f} (per fold {['%+.5f' % v for v in c['d_pds_A_per_fold']]})")
        P(f"  nmae identity: d wald89 {c['d_nmae_wald89']:+.6f} (3 cancer folds {c['d_nmae_wald89_3c']:+.6f}); dU {c['dU']:+.2f} dM {c['dM']:+.2f} e {c['e']:.3f}; "
          f"d precision89 {c['d_precision89']:+.5f}; d q_w {c['d_qw']:+.5f}")
        P(f"  amplitude vs cd4svd5: effect-level {c['amp_effect_level_pct']:+.3f}%, calibrated {c['amp_calibrated_pct']:+.3f}% (charged in every scenario)")
        P(f"  kolf2 (apart): rel cos_agg {100*c['kolf2']['rel_agg']:+.3f}%, d pds {c['kolf2']['d_pds']:+.5f}")
        for kk in (1, 3, 5):
            r = c["conc"][kk]
            P(f"  leave-top-{kk}-out pds retention vs cd4svd5: {100*r['retained']:.1f}% (top-{kk} share {100*r['share']:.1f}%)" + (f"; top-3 per fold: {r['top_names']}" if kk == 3 else ""))
        pi = c["prod_info"]
        if pi:
            P(f"  production-frame composition: union 1 {pi['union_carrying']}/300 carrying, solved fraction {pi['union_fraction']}; "
              f"union 2 {pi['union2_carrying']}/300 carrying (shared with union 1 {pi['union2_shared_targets']}), solved fraction {pi['union2_fraction']}; "
              f"CD4 fraction {pi['cd4_fraction']}; preserve {pi['preserve_scale']} (cd4svd5 0.9115882187626339); deflation energy {pi['energy_removed']}")
        P(f"  {'scenario':14s}{'k_cos':>6s}{'k_pds':>6s}{'DE3':>7s}{'nmae x':>7s}{'d cos':>10s}{'raw mse':>9s}{'| pds':>9s}{'mse':>9s}{'nmae':>9s}{'reach':>9s}{'fid':>9s}{'jac':>9s}{'| amp':>9s}{'TOTAL':>10s}")
        for which in ("PESSIMISTIC", "POINT", "OPTIMISTIC"):
            s = c["scenarios"][which]
            P(f"  {which:14s}{s['k_cos']:>6.2f}{s['k_pds']:>6.2f}{s['de3_rate']:>7.4f}{s['nmae_x']:>7.1f}{s['dcos']:>+10.6f}{s['raw']:>9.5f}{s['o_pds']:>+9.5f}{s['o_mse']:>+9.5f}"
              f"{s['o_nmae']:>+9.5f}{s['o_reach']:>+9.5f}{s['o_fid']:>+9.5f}{s['o_jac']:>+9.5f}{s['o_amp']:>+9.5f}{s['total']:>+10.5f}")
        s = c["cons258"]
        P(f"  {'258-conservative':14s} (k_pds 0.7 on FRAME A {c['d_pds_A']:+.5f}, k_cos 0.96, DE3 0.0425, nmae x1.5): pds {s['o_pds']:+.5f} mse {s['o_mse']:+.5f} DE3 {s['o_de3']:+.5f} nmae {s['o_nmae']:+.5f} TOTAL {s['total']:+.5f}")
        s = c["point_A"]
        P(f"  {'258-point':14s} (k_pds 3.79 on FRAME A, DE3 0.24, nmae x1.0): TOTAL {s['total']:+.5f}")
        P(f"  gate: PESSIMISTIC >= +{GATE} -> {'yes' if c['clears_pessimistic'] else 'NO'}; leave-top-3-out >= 50 % -> {'yes' if c['retention_ok'] else 'NO'}; "
          f"pds > 0 on 3/3 folds at >= 2 LOO-sigma -> {'yes' if c['pds_ok'] else 'NO'}  =>  {'CLEARED' if c['cleared'] else 'NOT CLEARED'}")

    for n in dict.fromkeys([bestA, bestB, best_all]):
        full(n)

    # ---- verdict -------------------------------------------------------------------------------------
    cleared = [n for n in real_cells if cells[n]["cleared"]]
    pess_clear = [n for n in real_cells if cells[n]["clears_pessimistic"]]
    point_clear = [n for n in real_cells if cells[n]["scenarios"]["POINT"]["total"] >= GATE]
    opt_clear = [n for n in real_cells if cells[n]["scenarios"]["OPTIMISTIC"]["total"] >= GATE]
    P("\n" + "=" * 130)
    P("VERDICT")
    P("=" * 130)
    P(f"  R1 (258's rows reproduce at 0.0e+00): {'PASS' if r1_pass else 'FAIL'}")
    P(f"  R2 (a) BR237 ARM6 by member at 258 s0's factor set: {'PASS' if r2a_pass else 'MISS'}; (b) BR249 cd4svd5 vs h1union +0.0043 +- 0.0008: {'PASS' if r2b_pass else 'MISS'} ({p_b['total']:+.5f})")
    P(f"  R3 (248's second residual through the union2 path): {'PASS' if r3_pass else 'FAIL'}")
    P(f"  real cells (no controls): {len(real_cells)}; PESSIMISTIC >= +{GATE}: {len(pess_clear)} {pess_clear}; POINT >= gate: {len(point_clear)} {point_clear}; "
      f"OPTIMISTIC >= gate: {len(opt_clear)} {opt_clear}")
    P(f"  cells meeting ALL THREE conditions (pessimistic gate, retention >= 50 %, pds 3/3 at >= 2 sigma): {len(cleared)} {cleared}")
    verdict = "CLEARS" if cleared else "DOES NOT CLEAR"
    P(f"\n{verdict}")
    out = dict(known_answer_crossfold=ka, known_answer_r3=ka3, r1_pass=r1_pass, r2a=r2a, r2a_pass=r2a_pass, r2b=p_b, r2b_de24=p_b_de, r2b_pass=r2b_pass, r3_pass=r3_pass,
               factors=dict(cos_h1union=cos_h, cos_ARM6=cos_6, cos_cd4svd5=cos_5, k_cos_ARM6=k_cos6, k_cos_cd4svd5=k_cos5, k_pds_ARM6_ABC=k_pds6_ABC,
                            k_pds_ARM6_A=k_pds6_A, k_pds_cd4svd5_ABC=k_pds5_ABC, de3_rate_ARM6=de3_rate6, de3_rate_cd4svd5=de3_5 / (cos_5 - cos_h),
                            nmae_ratio_ARM6=board6["nmae"] / (-d_nmae6 / NMAE_SPAN / 6), amp_calibration=amp_cal, board_delta_ARM6=board6, board_delta_cd4svd5=board5,
                            bands=dict(k_pds=K_PDS_BAND, k_cos=K_COS_BAND, de3=DE3_BAND, nmae_loss=NMAE_LOSS_BAND, point=POINT), de3_split=DE3_SPLIT),
               aggregates={n: {k: v for k, v in agg[n].items() if k != "per_fold"} for n in names},
               cells=cells, best_A=bestA, best_B=bestB, best_all=best_all, cleared=cleared, pessimistic_clear=pess_clear, point_clear=point_clear,
               optimistic_clear=opt_clear, verdict=verdict)
    (WORK / "price259.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "price259.log").write_text("\n".join(lines) + "\n")
    P(f"\nwrote {WORK/'price259.json'} and price259.log")


def _retention(arrs, name, base, kk, d_mean, targets=None):
    """222/225/248/258's leave-top-k-out on pds: per fold, rank targets by their own d pds, drop the top k, mean the rest;
    retained = mean of per-fold LOO means / the full 3-fold mean; share = the top k's part of the summed delta."""
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold")
    ap.add_argument("--prod", action="store_true")
    ap.add_argument("--price", action="store_true")
    a = ap.parse_args()
    if a.fold:
        run_fold(a.fold)
    elif a.prod:
        run_prod()
    elif a.price:
        run_price()
    else:
        ap.error("one of --fold <held>, --prod, --price")


if __name__ == "__main__":
    main()
