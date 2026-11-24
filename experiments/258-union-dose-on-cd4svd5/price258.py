"""258 -- the six-source union on top of `cd4svd5`: dose and subsets, priced on the instruments the
2026-09-05 board points calibrated (BR237: k_cos 0.96 for a new source into empty rows, the exact
`nmae` identity VALIDATED at 1.5x; BR249: k_cos 3.76 / k_pds 1.22 for a re-shape of occupied rows).
MEASUREMENT ONLY: no build, no `vcc`, no `cell-eval2`, nothing under `submissions/`.

Three phases, one file, every number an artifact on disk can be checked against:

  --fold <held>   stage one PLOLO fold (162), build every arm through the PRODUCTION builder's own
                  `apply_frozen_cd4_residual` (deflate k=5 -> size -> merge the union -> ONE preserve;
                  248/249's route, module globals as the knobs), and score it on
                    (a) 218's `score_effects` functional: `cos_agg` / macro / `pds` on the a23 tag,
                        mean over challenge frames A, B, C -- the frame 243's class rule and 246/249's
                        Delta pds +0.012263 live in; and
                    (b) 216/229's frame-A functional: 188's exact `nmae` identity on three pools,
                        top-89 precision, delivered norm ratio, `pds` on frame A -- the frame
                        `six229.json` lives in.
                  Writes work/fold-<held>.json and work/fold-<held>.npz (per-target `pds`/cosine, and
                  the identity's per-row U/M/D arrays).
  --prod          PRODUCTION-FRAME reconstruction, no 20 GiB staging: the whitened base is recovered
                  from 249's staged `effects-off.npy` (one invertible preserve scalar; `wh = final/s -
                  merged`), the incumbent `cd4svd5` is re-composed (KA: 249's `a28d5be2...` to numeric
                  identity) and so is 229's ARM6 (KA: `c5721825...`); then every screen arm is composed
                  on the real 300-target panel and its delivered amplitude is read with 211/243's
                  (A_c, B_c) terms on the effect-level render -- calibrated against the three packages
                  on disk (ARM6 -1.24 %, cd4svd5 -0.41 %, cd4svd10 -0.456 % vs cd4svd5).
  --price         the retrodiction rows R1/R2/R3, the dose and subset tables, the best cell,
                  leave-top-3-out retention, CLEARS / DOES NOT CLEAR -> work/price258.{json,log}.

    uv run python experiments/258-union-dose-on-cd4svd5/price258.py --fold k562_gwps   (x4 folds)
    uv run python experiments/258-union-dose-on-cd4svd5/price258.py --prod
    uv run python experiments/258-union-dose-on-cd4svd5/price258.py --price
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
# The builder reads its knobs at import.  This is the incumbent's configuration; every arm rebinds
# the module globals it varies (248's pattern) and nothing else.
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
E229 = ROOT / "experiments/229-stack4-build"
E243 = ROOT / "experiments/243-transfer-factor"
E246 = ROOT / "experiments/246-direction-hypotheses"
E249 = ROOT / "experiments/249-cd4-deflate"
E255 = ROOT / "experiments/255-cd4-deflate-k10"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/188-magnitude-calibration"))
sys.path.insert(0, str(ROOT / "experiments/203-reopen-125"))
sys.path.insert(0, str(E243))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"        # 258's OWN copies of the fold caches (sha-identical to 249's)

from plolo import (E062, FOLDS, OFFICIAL, PANEL, PROD_OPS, Operator, axis_sha256,  # noqa: E402
                   build_fold, cd4_for, deliver, endpoints)
from mapfit import to_lfc  # noqa: E402
from phase23 import POOL_N, key_precision, topn, umd  # noqa: E402
from lib243 import CP46, KC, SUMD  # noqa: E402
from cell_eval2.metrics.delta import _across_pert_budget, _row_weights  # noqa: E402

ALPHA, GAMMA, LAM = 0.982, 0.65, 0.0
NMAE_SPAN = 0.60                 # 188/216/229's identity span (overall = -d nmae / 0.60 / 6)
SPAN_MSE = 0.950162566           # the scorer's mse span (243/248/255)
PDS_SCALE = 2.2138               # mandate: overall pds = 2.2138 x k_pds x d pds / 6  (= 1/0.4517)
SIGMA_PDS = 0.001494             # 162's leave-one-arm-out sigma on raw PLOLO pds
GATE = 0.0019                    # shape-only gate
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
CTX = ("A", "B", "C")
N_GENES = 18533
POOLS = ("wald89", "lfc89", "h89")
TS = plolo.TS

PAIR = ("gse295214_h1", "h1_2025")
DLD1 = ("dld1_lowmoi",)
JIANG = ("jiang_5line_noK562",)
IGVF = ("igvf_cardiomyocyte",)
MYL = ("myllia_a549",)
SELECTIONS = {
    "pair": PAIR,
    "stack3": PAIR + DLD1,
    "pair+jiang": PAIR + JIANG,
    "pair+igvf": PAIR + IGVF,
    "pair+myllia": PAIR + MYL,
    "stack4": DLD1 + JIANG + IGVF + MYL,
    "stack5": PAIR + DLD1 + JIANG + IGVF,
    "stack6": PAIR + DLD1 + JIANG + IGVF + MYL,
}
DOSES6 = (0.10, 0.15, 0.20, 0.28, 0.35)
DOSES_SUB = (0.20, 0.28, 0.35)
# (name, deflate k, cd4 f, union g, selection).  The first seven are known-answer / retrodiction rows.
ARMS: list[tuple[str, int | None, float, float | None, str | None]] = [
    ("a15 f=0.10 SHIPPED", None, 0.10, None, None),
    ("cd4 f=0.34", None, 0.34, None, None),
    ("h1union g=0.28", None, 0.34, 0.28, "pair"),
    ("stack6 g=0.35 on h1union", None, 0.34, 0.35, "stack6"),
    ("stack3 g=0.35 on h1union", None, 0.34, 0.35, "stack3"),
    ("stack4 g=0.35 on h1union", None, 0.34, 0.35, "stack4"),
    ("cd4svd5 <- INCUMBENT", 5, 0.34, 0.28, "pair"),
]
ARMS += [(f"stack6 g={g:g}", 5, 0.34, g, "stack6") for g in DOSES6]
ARMS += [(f"stack3 g={g:g}", 5, 0.34, g, "stack3") for g in DOSES6]
ARMS += [(f"{sel} g={g:g}", 5, 0.34, g, sel)
         for sel in ("pair+jiang", "pair+igvf", "pair+myllia", "stack4", "stack5") for g in DOSES_SUB]
ARMS += [(f"pair g={g:g}", 5, 0.34, g, "pair") for g in (0.20, 0.35, 0.45)]
# beyond the mandate's grid, because the best cell of the first pass sat at its top (g = 0.35)
ARMS += [("stack3 g=0.45", 5, 0.34, 0.45, "stack3"), ("stack3 g=0.6", 5, 0.34, 0.60, "stack3"),
         ("stack6 g=0.45", 5, 0.34, 0.45, "stack6")]
INC = "cd4svd5 <- INCUMBENT"
OLD = "h1union g=0.28"
ARM6 = "stack6 g=0.35 on h1union"


def peak_gib() -> float:
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024**3


def load_builder():
    path = ROOT / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("price258_builder", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["price258_builder"] = module
    spec.loader.exec_module(module)
    return module, path


def sources_of(MOD, names):
    return tuple((n, MOD.UNION_SOURCE_REGISTRY[n]["path"]) for n in names)


def compose(MOD, wh, cd4, wcols, panel_positions, k, resf, union, g):
    """Production's own path (`apply_frozen_cd4_residual`): CD4 deflated at k (None = off) BEFORE
    sizing at f = resf, the union sized by the carrying-row rule at g, ONE merged add, ONE preserve.
    The module globals ARE the builder's knobs (248's pattern)."""
    MOD.CD4_DEFLATE_K = k
    MOD.CD4_TARGET_F = resf
    MOD.H1_UNION_G = g
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        final, cd4_audit, u_audit = MOD.apply_frozen_cd4_residual(
            wh, cd4, np.asarray(wcols, dtype=np.int64), np.asarray(panel_positions, dtype=np.int64),
            extra=union, extra_audit=None if union is None else {})
    info = dict(
        cd4_fraction=float(cd4_audit["fraction"]),
        union_fraction=(None if u_audit is None else float(u_audit["residual_frame"]["solved_fraction"])),
        union_carrying=(None if u_audit is None else int(u_audit["residual_frame"]["carrying_rows"])),
        union_coverage=(None if u_audit is None else float(u_audit["residual_frame"]["coverage_frac"])),
        preserve_scale=float(cd4_audit["post_add_preserve"]["post_add_global_preserve_scale"]),
        energy_removed=(None if "deflation" not in cd4_audit
                        else float(cd4_audit["deflation"]["energy_removed_frac"])),
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
    """Both functionals from six renders (2 tags x 3 frames)."""
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
                    out["pds_A"] = float(np.mean(p))          # six229's `pds` (frame A, a23 render)
            elif c == "A":
                H_A = dl                                         # six229's frame-A no-mass-fix render
        for k in acc["A"]:
            out[f"{tag}_{k}"] = float(np.mean([acc[c][k] for c in CTX]))
    out.update({k: out[f"a23_{k}"] for k in ("pds", "cos", "cos_agg", "sgn", "sgn_w", "sse")})
    out["median_row_norm_final"] = float(np.median(np.linalg.norm(final, axis=1)))
    # ---- six229's frame-A block: the exact nmae identity, precision, delivered amplitude ----
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
    needed = sorted({a[4] for a in ARMS if a[4] is not None})
    unions = {}
    for sel in needed:
        u, details = MOD.prepare_h1_union(targets, OFFICIAL, PANEL, max_rss_gib=40.0,
                                          sources=sources_of(MOD, SELECTIONS[sel]))
        unions[sel] = u
        print(f"[{held}] union {sel:12s} {int(np.sum(np.asarray(u.target_condition_count) > 0)):3d}/"
              f"{len(targets)} carried  ({time.time()-t0:.0f}s)", flush=True)
    op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
               else PROD_OPS / "whitener_union_pruned.npz")
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
    for name, k, resf, g, sel in ARMS:
        t1 = time.time()
        final, info = compose(MOD, wh, cd4, fold["wcols"], fold["panel_positions"], k, resf,
                              unions.get(sel), g)
        final[:, fold["panel_positions"]] = 0.0
        s, arr = score_arm(final, fold, Ts, ctx, LR, pools)
        rows[name] = dict(arm=name, k=k, f=resf, g=g, selection=sel, **info, **s,
                          seconds=round(time.time() - t1, 1))
        for kk, v in arr.items():
            arrays[f"{name}|{kk}"] = v
        print(f"[{held}] {name:28s} agg {s['cos_agg']:.6f} macro {s['a23_cos']:.6f} pds {s['pds']:.6f} "
              f"pdsA {s['pds_A']:.6f} prec {s['top89_precision']:.5f} fc {info['cd4_fraction']:.5f} "
              f"fu {info['union_fraction'] if info['union_fraction'] is None else round(info['union_fraction'], 5)} "
              f"({rows[name]['seconds']}s)", flush=True)
        del final
        gc.collect()

    # ---------------------------------------------------------------- known-answer rows, this fold
    ka = []

    def check(label, got, want, tol=1e-9):
        d = abs(float(got) - float(want))
        ok = d <= tol
        ka.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, ok=bool(ok)))
        print(f"  {'PASS' if ok else 'FAIL'} {label:70s} got {got:.10f} want {want:.10f} |d| {d:.1e}", flush=True)

    print(f"\n[{held}] KNOWN-ANSWER ROWS", flush=True)
    s246a = {r["arm"]: r for r in json.loads((E246 / f"work/sweep246-{held}.json").read_text())["rows"]}
    s246b = {r["arm"]: r for r in json.loads((E246 / f"work/sweep246b-{held}.json").read_text())["rows"]}
    for key in ("cos_agg", "a23_cos", "pds", "a15_pds", "a15_cos_agg", "a15_sgn_w"):
        check(f"246 sweep246 STACK.{key} == h1union", rows[OLD][key], s246a["STACK"][key])
    for key in ("cos_agg", "a23_cos", "pds", "a15_pds"):
        check(f"246 sweep246b B2_cd4_svd5_resolved.{key} == cd4svd5", rows[INC][key], s246b["B2_cd4_svd5_resolved"][key])
    check("246 B2_cd4_svd5_resolved solved fraction == cd4svd5", rows[INC]["cd4_fraction"],
          s246b["B2_cd4_svd5_resolved"]["audit"]["fraction"], 1e-12)
    r249 = json.loads((E249 / f"work/repro249-{held}.json").read_text())
    for key in ("cos_agg", "a23_cos", "pds"):
        check(f"249 repro249 arm.{key} == cd4svd5", rows[INC][key], r249["arm"][key])
    check("249 repro249 union fraction == cd4svd5 union fraction", rows[INC]["union_fraction"],
          r249["union_frame_solved_fraction"], 1e-12)
    p229 = json.loads((E229 / "work/price229.json").read_text())["arms"]
    for tag, mine in (("ARM6", ARM6), ("ARM3", "stack3 g=0.35 on h1union"), ("ARM4", "stack4 g=0.35 on h1union")):
        r = next(x for x in p229 if x["fold"] == held and x["arm"] == tag and x["on_top"] is None
                 and abs(x["g"] - 0.35) < 1e-12)
        check(f"229 price229 {tag}@0.35 cos_agg == {mine}", rows[mine]["cos_agg"], r["cos_agg"])
        check(f"229 price229 {tag}@0.35 pds (ABC) == {mine}", rows[mine]["pds"], r["pds"])
        check(f"229 price229 {tag}@0.35 macro == {mine}", rows[mine]["a23_cos"], r["cos"])
        check(f"229 price229 {tag}@0.35 union fraction == {mine}", rows[mine]["union_fraction"], r["fraction"], 1e-12)
    six229 = json.loads((E229 / "work/six229.json").read_text())["arms"]
    for theirs, mine in (("h1union g=0.28 <- INCUMBENT", OLD), ("stack6 g=0.35 <- ARM", ARM6),
                         ("stack4 g=0.35", "stack4 g=0.35 on h1union"), ("stack3 g=0.35", "stack3 g=0.35 on h1union"),
                         ("a15 f=0.10 SHIPPED", "a15 f=0.10 SHIPPED"), ("cd4 f=0.34", "cd4 f=0.34")):
        check(f"229 six229 {theirs} pds_per_fold (frame A) == {mine}", rows[mine]["pds_A"],
              six229[theirs]["pds_per_fold"][held])
    n_ok = sum(r["ok"] for r in ka)
    print(f"[{held}] known-answer rows {n_ok}/{len(ka)} pass", flush=True)

    out = dict(fold=held, n_targets=len(targets), targets=list(targets),
               builder_sha256=MOD.file_sha256(builder_path), rows=rows, known_answer=ka,
               union_carried={sel: int(np.sum(np.asarray(u.target_condition_count) > 0)) for sel, u in unions.items()},
               seconds=round(time.time() - t0, 1), peak_rss_gib=peak_gib())
    (WORK / f"fold-{held}.json").write_text(json.dumps(out, indent=1, default=float))
    np.savez_compressed(WORK / f"fold-{held}.npz", **arrays)
    print(f"[{held}] {len(rows)} arms in {out['seconds']}s, peak RSS {out['peak_rss_gib']:.2f} GiB; "
          f"wrote {WORK/f'fold-{held}.json'}", flush=True)
    if n_ok != len(ka):
        raise SystemExit(1)


# ------------------------------------------------------------------------------ production phase
def terms_from_final(final, targets, b, gate):
    """211/243's per-context (A_c, B_c, Kc_c, cap_c) on the EFFECT-LEVEL render (no emitter):
    `deliver` returns dh = log1p(TS p) - b exactly as 243's `terms` forms it from block totals."""
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
                           target_condition_count=np.asarray(z["target_condition_count"], dtype=np.uint8),
                           condition_details=[])
    u = np.load(stage / "union-off.npz", allow_pickle=True)
    union_pair_saved = E062.CD4Prepared(mode="union", mean_effect=np.asarray(u["mean_effect"], dtype=np.float64),
                                        gene_votes=np.asarray(u["gene_votes"], dtype=np.uint8),
                                        target_condition_count=np.asarray(u["target_condition_count"], dtype=np.uint8),
                                        condition_details=[])
    final_off = np.load(stage / "effects-off.npy")
    aud_off = json.loads((stage / "stage-off.json").read_text())
    aud_svd5 = json.loads((stage / "stage-svd5.json").read_text())
    s_off = float(aud_off["audit"]["h1_union"]["merged_preserve"]["post_add_global_preserve_scale"])
    fc = float(aud_off["audit"]["cd4"]["fraction"])
    fu = float(aud_off["audit"]["h1_union"]["residual_frame"]["solved_fraction"])
    base_norm = float(aud_off["audit"]["h1_union"]["residual_frame"]["whitened_base_median_row_norm"])
    print(f"249 stage-off: preserve {s_off:.16f} cd4 fraction {fc:.16f} union fraction {fu:.16f} base median {base_norm}", flush=True)
    # the production operator's columns = the kolf2 fold's (production operator, all three donors)
    zk = np.load(WORK / "plolo/fold-kolf2.npz", allow_pickle=True)
    wcols = np.asarray(zk["wcols"], dtype=np.int64)
    axis = axis_sha256(OFFICIAL[i] for i in wcols)
    assert axis == plolo.OPERATOR_REGISTRY["union_pruned"]["axis_sha256"], axis
    gpos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([gpos[t] for t in targets], dtype=np.int64)
    same_order = tuple(targets) == tuple(PANEL)
    print(f"production target order == plolo.PANEL: {same_order}; wcols {wcols.size} axis {axis[:12]} (registry pin)", flush=True)

    # ---- recover the whitened base: final = s (wh + merged)  ->  wh = final / s - merged ----------
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
        print(f"  {'PASS' if ok else 'FAIL'} {label:78s} got {got:.12g} want {want:.12g} |d| {d:.1e}", flush=True)

    print("\nKNOWN-ANSWER ROWS (production frame)", flush=True)
    check("recovered base: median row norm == stage-off whitened_base_median_row_norm", float(np.median(np.linalg.norm(wh, axis=1))), base_norm, 1e-9)
    # KA-0: the union the builder prepares for the production panel == the saved one, bit for bit
    union_pair, _ = MOD.prepare_h1_union(targets, OFFICIAL, PANEL, max_rss_gib=40.0, sources=sources_of(MOD, PAIR))
    d_u = float(np.abs(np.asarray(union_pair.mean_effect) - np.asarray(union_pair_saved.mean_effect)).max())
    same_votes = np.array_equal(np.asarray(union_pair.gene_votes), np.asarray(union_pair_saved.gene_votes))
    check("prepare_h1_union(PAIR) on the production panel == 249's union-off.npz (max|d|)", d_u, 0.0, 0.0)
    ka.append(dict(label="prepare_h1_union(PAIR) votes bit-identical", ok=bool(same_votes)))
    # KA-1: re-compose the shipped h1union (k off, pair @0.28) -> effects-off.npy
    f_off, i_off = compose(MOD, wh, cd4, wcols, panel_positions, None, 0.34, union_pair, 0.28)
    check("re-composed h1union == 249 effects-off.npy (abeb0749), max|d|", float(np.abs(f_off - final_off).max()), 0.0, 1e-9)
    check("re-composed h1union cd4 fraction", i_off["cd4_fraction"], fc, 1e-12)
    check("re-composed h1union union fraction", i_off["union_fraction"], fu, 1e-12)
    # KA-2: the incumbent cd4svd5 (k=5, pair @0.28) -> 249's effects-svd5.npy (a28d5be2...)
    final_svd5 = np.load(stage / "effects-svd5.npy")
    f_5, i_5 = compose(MOD, wh, cd4, wcols, panel_positions, 5, 0.34, union_pair, 0.28)
    check("re-composed cd4svd5 == 249 effects-svd5.npy (a28d5be2), max|d|", float(np.abs(f_5 - final_svd5).max()), 0.0, 1e-9)
    check("re-composed cd4svd5 cd4 fraction == stage-svd5", i_5["cd4_fraction"], float(aud_svd5["audit"]["cd4"]["fraction"]), 1e-12)
    check("re-composed cd4svd5 preserve == stage-svd5", i_5["preserve_scale"],
          float(aud_svd5["audit"]["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"]), 1e-12)
    check("re-composed cd4svd5 energy removed == stage-svd5", i_5["energy_removed"],
          float(aud_svd5["audit"]["cd4"]["deflation"]["energy_removed_frac"]), 1e-12)
    check("effects-svd5.npy array sha == a28d5be2... (249 record)", 1.0 if MOD.array_sha256(final_svd5).startswith("a28d5be2") else 0.0, 1.0, 0.0)
    # KA-3: 229's ARM6 (k off, six @0.35) -> 229's effects-stack6g035.npy (c5721825...)
    final_arm6 = np.load(E229 / "work/effects-stack6g035.npy")
    union6, _ = MOD.prepare_h1_union(targets, OFFICIAL, PANEL, max_rss_gib=40.0, sources=sources_of(MOD, SELECTIONS["stack6"]))
    f_6, i_6 = compose(MOD, wh, cd4, wcols, panel_positions, None, 0.34, union6, 0.35)
    check("re-composed ARM6 == 229 effects-stack6g035.npy (c5721825), max|d|", float(np.abs(f_6 - final_arm6).max()), 0.0, 1e-9)
    aud6 = json.loads((E229 / "work/stage-stack6g035.json").read_text())["audit"]
    check("re-composed ARM6 union fraction == 229 stage", i_6["union_fraction"], float(aud6["h1_union"]["residual_frame"]["solved_fraction"]), 1e-12)
    check("effects-stack6g035.npy array sha == c5721825... (229 record)", 1.0 if MOD.array_sha256(final_arm6).startswith("c5721825") else 0.0, 1.0, 0.0)
    final_svd10 = np.load(E255 / "work/stage/effects-svd10.npy")
    n_ok = sum(r["ok"] for r in ka)
    print(f"known-answer rows {n_ok}/{len(ka)} pass  ({time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB)", flush=True)
    if n_ok != len(ka):
        (WORK / "prod258.json").write_text(json.dumps(dict(known_answer=ka, failed=True), indent=1))
        raise SystemExit(1)

    # ---- effect-level delivered terms for every arm on the production panel ------------------
    ctx = challenge_contexts()
    b = {c: ctx[c][0] for c in CTX}
    gate = {c: ctx[c][1] for c in CTX}
    bA, gA = ctx["A"]
    unions = {"pair": union_pair, "stack6": union6}
    for sel in sorted({a[4] for a in ARMS if a[4] is not None} - set(unions)):
        unions[sel], _ = MOD.prepare_h1_union(targets, OFFICIAL, PANEL, max_rss_gib=40.0, sources=sources_of(MOD, SELECTIONS[sel]))
    terms: dict[str, dict] = {}
    info_all: dict[str, dict] = {}
    finals = {"h1union g=0.28 (file)": final_off, "cd4svd5 (file)": final_svd5, "ARM6 (file)": final_arm6,
              "cd4svd10 (255 file)": final_svd10}
    for name, k, resf, g, sel in ARMS:
        if name in ("a15 f=0.10 SHIPPED", "cd4 f=0.34"):
            continue
        final, info = compose(MOD, wh, cd4, wcols, panel_positions, k, resf, unions.get(sel), g)
        finals[name] = final
        info_all[name] = info
    for name, final in finals.items():
        t1 = time.time()
        per = terms_from_final(final, targets, b, gate)
        H = deliver(final, bA, gA, None)
        terms[name] = dict(per=per, delivered_norm_A=float(np.median(np.linalg.norm(H[:, wcols], axis=1))),
                           median_row_norm_final=float(np.median(np.linalg.norm(final, axis=1))))
        print(f"  {name:32s} A {per['A']['A']:.6f}/{per['B']['A']:.6f}/{per['C']['A']:.6f}  B/A {np.mean([per[c]['B']/per[c]['A'] for c in CTX]):.5f} "
              f"dlvA {terms[name]['delivered_norm_A']:.5f} ({time.time()-t1:.0f}s)", flush=True)
        del H
    # ---- calibration of the effect-level amplitude against the three packages on disk -----------
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
        print(f"  amplitude calibration {arm:8s} vs {base:8s}: package {100*(r_pkg-1):+.4f}%  effect-level A_c {100*(r_eff-1):+.4f}%  "
              f"216 median-norm frame A {100*(r_med-1):+.4f}%", flush=True)
    out = dict(known_answer=ka, targets_in_panel_order=same_order, calibration=calib,
               terms={k: v for k, v in terms.items()}, compose_info=info_all,
               package_terms=pkg, seconds=round(time.time() - t0, 1), peak_rss_gib=peak_gib())
    (WORK / "prod258.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"wrote {WORK/'prod258.json'}  ({out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB)", flush=True)


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
                         fid=e["de_wilcoxon_direction_fidelity_yield_raw"], reach=e["de_wilcoxon_direction_reach_raw"],
                         jac=e["de_wilcoxon_sig_jaccard"]),
                scaled={m: e[f"score_{m}"] for m in ("pds", "mse", "nmae", "fid", "reach", "jac")})


def run_price() -> None:
    lines = []

    def P(t=""):
        print(t, flush=True)
        lines.append(t)

    folds = {h: json.loads((WORK / f"fold-{h}.json").read_text()) for h in FOLDS}
    arrs = {h: np.load(WORK / f"fold-{h}.npz") for h in FOLDS}
    prod = json.loads((WORK / "prod258.json").read_text())
    names = [a[0] for a in ARMS]
    ka_fold = {h: (sum(r["ok"] for r in folds[h]["known_answer"]), len(folds[h]["known_answer"])) for h in FOLDS}
    ka_prod = (sum(r["ok"] for r in prod["known_answer"]), len(prod["known_answer"]))
    P("=" * 120)
    P("258 -- the six-source union on top of cd4svd5, priced on BR237/BR249's factors  (PLOLO 162, three cancer folds; "
      "kolf2 pooled ONLY into the 4-fold nmae identity as six229 does)")
    P("=" * 120)
    P(f"per-fold known-answer rows: " + ", ".join(f"{h} {a}/{b}" for h, (a, b) in ka_fold.items()) + f"; production-frame {ka_prod[0]}/{ka_prod[1]}")
    for h in FOLDS:
        P(f"  {h:14s} {folds[h]['n_targets']} targets, {folds[h]['seconds']}s, peak RSS {folds[h]['peak_rss_gib']:.2f} GiB")
    P(f"  production frame: {prod['seconds']}s, peak RSS {prod['peak_rss_gib']:.2f} GiB")

    # ---- the cross-fold aggregates of both functionals -------------------------------------------
    agg: dict[str, dict] = {}
    for n in names:
        r = {h: folds[h]["rows"][n] for h in FOLDS}
        a: dict[str, object] = dict(per_fold=r)
        for key in ("cos_agg", "a23_cos", "pds", "pds_A", "a15_sgn_w", "top89_precision", "delivered_norm_A", "cd4_fraction",
                    "union_fraction", "preserve_scale"):
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
        agg[n]["delivered_norm_ratio"] = agg[n]["delivered_norm_A_4f"] / a15n     # six229's `dlv/a15`

    # ---- six229 / six216 / six246 cross-fold known-answer rows ----------------------------------
    P("\nCROSS-FOLD KNOWN-ANSWER ROWS (the frame-A functional, 4 folds pooled as six229 does)")
    six229 = json.loads((E229 / "work/six229.json").read_text())["arms"]
    six246 = json.loads((E246 / "work/six246.json").read_text())["arms"]
    six216 = json.loads((ROOT / "experiments/216-cd4-hard/work/six216.json").read_text())
    ka = []

    def check(label, got, want, tol=5e-6):
        d = abs(float(got) - float(want))
        ok = d <= tol
        ka.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, ok=bool(ok)))
        P(f"  {'PASS' if ok else 'FAIL'} {label:74s} got {got:.8f} want {want:.8f} |d| {d:.1e}")

    pairs = [("a15 f=0.10 SHIPPED", "a15 f=0.10 SHIPPED", six229), ("cd4 f=0.34", "cd4 f=0.34", six229),
             (OLD, "h1union g=0.28 <- INCUMBENT", six229), (ARM6, "stack6 g=0.35 <- ARM", six229),
             ("stack4 g=0.35 on h1union", "stack4 g=0.35", six229), ("stack3 g=0.35 on h1union", "stack3 g=0.35", six229),
             (INC, "B2_cd4_svd5_resolved", six246)]
    for mine, theirs, src in pairs:
        ref = src[theirs]
        tag = "six229" if src is six229 else "six246"
        check(f"{tag} {theirs} nmae_wald89 == {mine}", agg[mine]["nmae_wald89"], ref["nmae_wald89"])
        check(f"{tag} {theirs} top89_precision == {mine}", agg[mine]["top89_precision_4f"], ref["top89_precision"])
        check(f"{tag} {theirs} delivered_norm_ratio == {mine}", agg[mine]["delivered_norm_ratio"], ref["delivered_norm_ratio"])
        check(f"{tag} {theirs} pds_3cancer (frame A) == {mine}", agg[mine]["pds_A_3c"], ref["pds_3cancer"])
        check(f"{tag} {theirs} nmae_lfc89 == {mine}", agg[mine]["nmae_lfc89"], ref["nmae_lfc89"])
        check(f"{tag} {theirs} nmae_h89 == {mine}", agg[mine]["nmae_h89"], ref["nmae_h89"])
    check("six216 a15 nmae_wald89 == a15 f=0.10", agg["a15 f=0.10 SHIPPED"]["nmae_wald89"], six216["a15 f=0.10 SHIPPED"]["nmae_wald89"])
    check("six216 res f=0.34 nmae_wald89 == cd4 f=0.34", agg["cd4 f=0.34"]["nmae_wald89"], six216["res f=0.34 <- INCUMBENT"]["nmae_wald89"])
    # R1 by the mandate's digits
    d_nmae_arm6 = agg[ARM6]["nmae_wald89"] - agg[OLD]["nmae_wald89"]
    check("R1 stack6@0.35 on h1union: nmae_wald89 == 1.0029112 (mandate)", agg[ARM6]["nmae_wald89"], 1.0029112, 5e-8)
    check("R1 stack6@0.35 on h1union: overall_nmae_vs_incumbent == -0.0047829 (mandate)", -d_nmae_arm6 / NMAE_SPAN / 6, -0.0047829, 5e-8)
    n_ok = sum(r["ok"] for r in ka)
    P(f"cross-fold known-answer rows {n_ok}/{len(ka)} pass")
    r1_pass = n_ok == len(ka)

    # ---- board rows and the calibrated factors ---------------------------------------------------
    B_h = board_entry("iUziqE1gewfcPKzYUDaB")      # h1union
    B_6 = board_entry("ydAemVdWmsTx9xhDWcou")      # ARM6
    B_5 = board_entry("Xv0ayvuoHpsKYWv2iu42")      # cd4svd5
    pkg = prod["package_terms"]
    cos_h = invert(pkg["h1union"], B_h["raw"]["mse"])
    cos_5 = invert(pkg["cd4svd5"], B_5["raw"]["mse"])
    cos_6 = invert(pkg["ARM6"], B_6["raw"]["mse"])
    P("\nBOARD-CALIBRATED FACTORS")
    P(f"  measured board cosines (197's identity at each package's own terms): h1union {cos_h:.6f}  ARM6 {cos_6:.6f}  cd4svd5 {cos_5:.6f}")

    def delta(a, b, key, folds_=CANCER, rel=False):
        vals = [agg[a]["per_fold"][h][key] / agg[b]["per_fold"][h][key] - 1 if rel else agg[a]["per_fold"][h][key] - agg[b]["per_fold"][h][key]
                for h in folds_]
        return float(np.mean(vals)), [float(v) for v in vals]

    rel6, rel6_pf = delta(ARM6, OLD, "cos_agg", rel=True)
    dp6_ABC, _ = delta(ARM6, OLD, "pds")
    dp6_A, _ = delta(ARM6, OLD, "pds_A")
    rel5, rel5_pf = delta(INC, OLD, "cos_agg", rel=True)
    dp5_ABC, _ = delta(INC, OLD, "pds")
    dp5_A, _ = delta(INC, OLD, "pds_A")
    k_cos6 = (cos_6 - cos_h) / (rel6 * cos_h)
    k_cos5 = (cos_5 - cos_h) / (rel5 * cos_h)
    dpds6_board = B_6["raw"]["pds"] - B_h["raw"]["pds"]
    dpds5_board = B_5["raw"]["pds"] - B_h["raw"]["pds"]
    k_pds6_A, k_pds6_ABC = dpds6_board / dp6_A, dpds6_board / dp6_ABC
    k_pds5_A, k_pds5_ABC = dpds5_board / dp5_A, dpds5_board / dp5_ABC
    dcos6_board, dcos5_board = cos_6 - cos_h, cos_5 - cos_h
    de3_6 = sum(B_6["scaled"][m] - B_h["scaled"][m] for m in ("fid", "reach", "jac")) / 6
    de3_5 = sum(B_5["scaled"][m] - B_h["scaled"][m] for m in ("fid", "reach", "jac")) / 6
    nm6_board = (B_6["scaled"]["nmae"] - B_h["scaled"]["nmae"]) / 6
    nm6_ident = -d_nmae_arm6 / NMAE_SPAN / 6
    P(f"  ARM6  (new source into empty rows, BR237): PLOLO rel cos_agg {100*rel6:+.4f}% (per fold {['%+.3f' % (100*v) for v in rel6_pf]}), "
      f"board d cos {dcos6_board:+.6f} -> k_cos {k_cos6:.3f}")
    P(f"        d pds PLOLO frame A {dp6_A:+.6f} / frames ABC {dp6_ABC:+.6f}; board raw d pds {dpds6_board:+.6f} -> "
      f"k_pds {k_pds6_A:.3f} (frame A, six229's functional) / {k_pds6_ABC:.3f} (frames ABC, 246's functional)")
    P(f"        DE-three realised {de3_6:+.6f} overall = {de3_6/dcos6_board:.4f} per unit board d cos (the 0.24 recipe would give {0.24*dcos6_board:+.6f}; "
      f"243 s4.3's DE-FOUR rate 0.239 carries nmae at 43 %, so its DE-THREE share is 0.136 -> {0.136*dcos6_board:+.6f})")
    P(f"        nmae: identity {nm6_ident:+.6f} vs board {nm6_board:+.6f} -> realised ratio {nm6_board/nm6_ident:.3f} (the mandate's x1.5)")
    P(f"  cd4svd5 (re-shape of occupied rows, BR249): PLOLO rel cos_agg {100*rel5:+.4f}%, board d cos {dcos5_board:+.6f} -> k_cos {k_cos5:.3f}; "
      f"d pds A {dp5_A:+.6f} / ABC {dp5_ABC:+.6f}, board {dpds5_board:+.6f} -> k_pds {k_pds5_A:.3f} (A) / {k_pds5_ABC:.3f} (ABC); "
      f"DE-three realised {de3_5:+.6f} = {de3_5/dcos5_board:.4f} per unit d cos")
    K_COS_NEW = 0.96
    K_DE3_ARM6 = de3_6 / dcos6_board
    K_DE3_RECIPE = 0.24
    K_DE3_SHARE = 0.136
    NMAE_X_CONS = 1.5

    # ---- amplitude: effect-level A_c on the production panel, calibrated ----------------------------
    P("\nDELIVERED AMPLITUDE: effect-level (A_c, B_c) on the production 300-target panel vs the packages on disk")
    for c in prod["calibration"]:
        P(f"  {c['arm']:8s} vs {c['base']:8s}: package {c['package_pct']:+.4f}%   effect-level {c['effect_level_pct']:+.4f}%   "
          f"216 median-norm frame A {c['median_norm_A_pct']:+.4f}%   package/effect-level {c['ratio_pkg_over_eff'] if c['ratio_pkg_over_eff'] is None else round(c['ratio_pkg_over_eff'], 3)}")
    ratios = [c["ratio_pkg_over_eff"] for c in prod["calibration"] if c["ratio_pkg_over_eff"] is not None]
    amp_cal = float(np.mean(ratios))
    P(f"  calibration factor package/effect-level, mean over the three pairs: {amp_cal:.3f} (applied to every screen cell's effect-level "
      f"amplitude below; the un-calibrated read is carried alongside)")
    T = prod["terms"]

    def amp_vs(name, base_name):
        r_eff = amp_ratio(T[name]["per"], T[base_name]["per"])
        return r_eff, 1.0 + amp_cal * (r_eff - 1.0)

    # ---- the price function --------------------------------------------------------------------
    def price(rel_agg, d_pds_A, d_pds_ABC, d_nmae, r_amp, base_per, cos_base, raw_base,
              k_cos, k_pds_A, k_pds_ABC, de_rate, nmae_x):
        dcos = k_cos * rel_agg * cos_base
        per_arm = scale_terms(base_per, r_amp)
        raw = raw_for(per_arm, cos_base + dcos)
        raw_amp_only = raw_for(per_arm, cos_base)
        o_mse = -(raw - raw_base) / SPAN_MSE / 6
        o_mse_cos_only = -(raw_for(base_per, cos_base + dcos) - raw_base) / SPAN_MSE / 6
        o_pds_A = PDS_SCALE * k_pds_A * d_pds_A / 6
        o_pds_ABC = PDS_SCALE * k_pds_ABC * d_pds_ABC / 6
        o_de3 = de_rate * dcos
        o_nmae = -d_nmae / NMAE_SPAN / 6 * nmae_x
        return dict(dcos=dcos, raw=raw, raw_amp_only=raw_amp_only, o_mse=o_mse, o_mse_cos_only=o_mse_cos_only,
                    o_amp=-(raw_amp_only - raw_base) / SPAN_MSE / 6, o_pds_A=o_pds_A, o_pds_ABC=o_pds_ABC, o_de3=o_de3,
                    o_nmae=o_nmae, total_A=o_mse + o_pds_A + o_de3 + o_nmae, total_ABC=o_mse + o_pds_ABC + o_de3 + o_nmae)

    # ---- R2: ARM6 six-member price vs the board ---------------------------------------------------
    P("\nR2 -- ARM6 (stack6 g=0.35 on h1union) six-member price vs BR237's board point (overall -0.0036403; "
      "pds +0.00248, mse +0.00088, nmae -0.00716)")
    r6_pkg = amp_ratio(pkg["ARM6"], pkg["h1union"])
    r6_eff, r6_cal = amp_vs("ARM6 (file)", "h1union g=0.28 (file)")
    P(f"  ARM6 delivered amplitude vs h1union: package {100*(r6_pkg-1):+.3f}%  effect-level {100*(r6_eff-1):+.3f}%  calibrated {100*(r6_cal-1):+.3f}%")
    grid = []
    r2_table = []
    for amp_label, r_amp in (("package amplitude", r6_pkg), ("calibrated effect-level", r6_cal), ("no amplitude leg", 1.0)):
        for kp_label, kpA, kpABC in (("k_pds 0.7", 0.7, 0.7), (f"k_pds ARM6-realised ({k_pds6_A:.2f} A / {k_pds6_ABC:.2f} ABC)", k_pds6_A, k_pds6_ABC)):
            for de_label, de in (("DE3 0.24 recipe", K_DE3_RECIPE), ("DE3 0.136 share", K_DE3_SHARE), (f"DE3 ARM6-realised {K_DE3_ARM6:.3f}", K_DE3_ARM6)):
                for nx in (1.0, 1.5):
                    p = price(rel6, dp6_A, dp6_ABC, d_nmae_arm6, r_amp, pkg["h1union"], cos_h, B_h["raw"]["mse"],
                              K_COS_NEW, kpA, kpABC, de, nx)
                    ok_tot = abs(p["total_A"] - (-0.0036)) <= 0.0010
                    ok_pds = abs(p["o_pds_A"] - 0.0025) <= 0.0008
                    ok_mse = abs(p["o_mse"] - 0.0009) <= 0.0003
                    ok_nm = abs(p["o_nmae"] - (-0.0072)) <= 0.0025
                    row = dict(amplitude=amp_label, k_pds=kp_label, de3=de_label, nmae_x=nx, **p,
                               pass_total=ok_tot, pass_pds=ok_pds, pass_mse=ok_mse, pass_nmae=ok_nm,
                               pass_all=bool(ok_tot and ok_pds and ok_mse and ok_nm))
                    r2_table.append(row)
    for row in r2_table:
        P(f"  [{'ALL PASS' if row['pass_all'] else 'total ' + ('ok' if row['pass_total'] else 'MISS') + ', pds ' + ('ok' if row['pass_pds'] else 'MISS') + ', mse ' + ('ok' if row['pass_mse'] else 'MISS') + ', nmae ' + ('ok' if row['pass_nmae'] else 'MISS')}] "
          f"{row['amplitude']:24s} {row['k_pds']:44s} {row['de3']:26s} nmae x{row['nmae_x']:.1f}: pds {row['o_pds_A']:+.5f} mse {row['o_mse']:+.5f} "
          f"DE3 {row['o_de3']:+.5f} nmae {row['o_nmae']:+.5f} TOTAL {row['total_A']:+.5f}")
    r2_pass_sets = [r for r in r2_table if r["pass_all"]]
    P(f"  -> {len(r2_pass_sets)}/{len(r2_table)} factor sets land inside every band; the sets that do: "
      + ("; ".join(f"[{r['amplitude']} | {r['k_pds']} | {r['de3']} | nmae x{r['nmae_x']}]" for r in r2_pass_sets) if r2_pass_sets else "NONE"))

    # ---- R3: cd4svd5 vs h1union on the same recipe (k_cos 3.76, k_pds 1.22) ------------------------
    P("\nR3 -- cd4svd5 vs h1union on the same recipe at the re-shape factors (k_cos 3.76 / k_pds 1.22; BR249 +0.0042634)")
    d_nmae_5 = agg[INC]["nmae_wald89"] - agg[OLD]["nmae_wald89"]
    r5_pkg = amp_ratio(pkg["cd4svd5"], pkg["h1union"])
    r5_eff, r5_cal = amp_vs("cd4svd5 (file)", "h1union g=0.28 (file)")
    r3_rows = []
    for amp_label, r_amp in (("package amplitude", r5_pkg), ("calibrated effect-level", r5_cal)):
        for de_label, de in (("DE3 0.24 recipe", K_DE3_RECIPE), ("DE3 0.136 share", K_DE3_SHARE), (f"DE3 cd4svd5-realised {de3_5/dcos5_board:.3f}", de3_5 / dcos5_board)):
            for nx in (1.0, 1.5):
                p = price(rel5, dp5_A, dp5_ABC, d_nmae_5, r_amp, pkg["h1union"], cos_h, B_h["raw"]["mse"], 3.76, 1.22 * dp5_ABC / dp5_A, 1.22, de, nx)
                ok = abs(p["total_ABC"] - 0.0043) <= 0.0008
                r3_rows.append(dict(amplitude=amp_label, de3=de_label, nmae_x=nx, **p, pass_total=ok))
                P(f"  [{'PASS' if ok else 'MISS'}] {amp_label:24s} {de_label:30s} nmae x{nx:.1f}: pds {p['o_pds_ABC']:+.5f} (k 1.22 on the ABC frame) "
                  f"mse {p['o_mse']:+.5f} DE3 {p['o_de3']:+.5f} nmae {p['o_nmae']:+.5f} TOTAL {p['total_ABC']:+.5f}")
    r3_pass = any(r["pass_total"] for r in r3_rows)
    P(f"  (board: pds +0.005452, mse -0.000799, nmae -0.000045, DE-three -0.000344; identity nmae leg for cd4svd5 vs h1union {-d_nmae_5/NMAE_SPAN/6:+.6f})")

    # ---- THE SCREEN: every cell vs cd4svd5 --------------------------------------------------------
    P("\n" + "=" * 120)
    P("THE SCREEN -- every cell vs cd4svd5 (the incumbent), three cancer folds; nmae identity pooled over 4 folds (six229's convention)")
    P("=" * 120)
    base_per = pkg["cd4svd5"]
    raw_5 = B_5["raw"]["mse"]
    chk = raw_for(base_per, cos_5)
    P(f"  197's identity at cd4svd5's own terms reproduces its board raw mse: {chk:.7f} vs {raw_5:.7f} (|d| {abs(chk-raw_5):.1e}); cos_cd4svd5 {cos_5:.6f}")
    P(f"  factors: k_cos {K_COS_NEW} (band 0.94-1.00; ARM6 realised {k_cos6:.3f}); k_pds conservative 0.7, point = ARM6-realised "
      f"{k_pds6_A:.2f} (frame A) / {k_pds6_ABC:.2f} (frames ABC), sensitivity 1.22 (BR249); DE-three conservative {K_DE3_ARM6:.3f}/unit d cos (ARM6 realised), "
      f"point 0.24 (mandate recipe); nmae identity x1.0 point, x{NMAE_X_CONS} conservative; mse from 197's identity at the arm's own terms "
      f"(amplitude leg from the calibrated effect-level read)")
    cells = {}
    screen = [a for a in ARMS if a[1] == 5 and a[0] != INC]
    for name, k, resf, g, sel in screen:
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
        # per-target retention (222's definition), vs cd4svd5
        conc = {}
        for kk in (1, 3, 5):
            loo, share = [], []
            for h in CANCER:
                d = arrs[h][f"{name}|pds_pt"] - arrs[h][f"{INC}|pds_pt"]
                order = np.argsort(-d)
                keep = np.ones(d.size, dtype=bool)
                keep[order[:kk]] = False
                loo.append(float(d[keep].mean()))
                share.append(float(d[order[:kk]].sum() / d.sum()) if d.sum() != 0 else float("nan"))
            conc[kk] = dict(loo_mean=float(np.mean(loo)), retained=(float(np.mean(loo) / dpABC) if dpABC != 0 else float("nan")),
                            share=float(np.nanmean(share)))
        cons = price(rel, dpA, dpABC, d_nmae, r_cal, base_per, cos_5, raw_5, K_COS_NEW, 0.7, 0.7, K_DE3_ARM6, NMAE_X_CONS)
        point = price(rel, dpA, dpABC, d_nmae, r_cal, base_per, cos_5, raw_5, K_COS_NEW, k_pds6_A, k_pds6_ABC, K_DE3_RECIPE, 1.0)
        mid = price(rel, dpA, dpABC, d_nmae, r_cal, base_per, cos_5, raw_5, K_COS_NEW, 1.22, 1.22, K_DE3_SHARE, 1.0)
        best_case = price(rel, dpA, dpABC, d_nmae, r_cal, base_per, cos_5, raw_5, 1.36, k_pds6_A, k_pds6_ABC, K_DE3_RECIPE, 1.0)
        cells[name] = dict(arm=name, g=g, selection=sel, sources=list(SELECTIONS[sel]),
                           carrying={h: folds[h]["rows"][name]["union_carrying"] for h in FOLDS},
                           union_fraction={h: folds[h]["rows"][name]["union_fraction"] for h in FOLDS},
                           prod_info=prod["compose_info"].get(name),
                           rel_cos_agg=rel, rel_cos_agg_per_fold=rel_pf, rel_macro=relm,
                           d_pds_A=dpA, d_pds_A_per_fold=dpA_pf, d_pds_ABC=dpABC, d_pds_ABC_per_fold=dpABC_pf,
                           pds_pos_folds=int(sum(v > 0 for v in dpABC_pf)), pds_sigma=abs(dpABC) / SIGMA_PDS,
                           pds_fold_se=float(np.std(dpABC_pf, ddof=1) / np.sqrt(3)),
                           d_qw=dqw, d_nmae_wald89=d_nmae, d_nmae_wald89_3c=d_nmae3, dU=dU, dM=dM,
                           e=(dM / dU if dU != 0 else float("nan")),
                           d_precision89=agg[name]["top89_precision_4f"] - agg[INC]["top89_precision_4f"],
                           amp_effect_level_pct=100 * (r_eff - 1), amp_calibrated_pct=100 * (r_cal - 1),
                           amp_plolo_median_pct=100 * (agg[name]["delivered_norm_ratio"] / agg[INC]["delivered_norm_ratio"] - 1),
                           kolf2=kolf, conc=conc, conservative=cons, point=point, mid=mid, best_case=best_case,
                           clears_conservative=bool(cons["total_A"] > GATE and cons["total_ABC"] > GATE),
                           clears_point=bool(point["total_A"] > GATE and point["total_ABC"] > GATE))

    def table(title, subset):
        P(f"\n{title}")
        P(f"  {'cell':22s}{'carry':>6s}{'agg%':>8s}{'mac%':>8s}{'dpds A':>10s}{'dpds ABC':>10s}{'+':>3s}{'sig':>5s}{'dq_w':>9s}{'d nmae':>9s}{'e':>7s}"
          f"{'amp%':>7s}{'|':>2s}{'pds.7':>8s}{'mse':>8s}{'DE3':>8s}{'nmae1.5':>9s}{'CONS':>9s}{'|':>2s}{'pds pt':>8s}{'DE3.24':>8s}{'nmae1':>8s}{'POINT':>9s}{'ret3':>7s}")
        for n in subset:
            c = cells[n]
            cs, pt = c["conservative"], c["point"]
            P(f"  {n:22s}{c['carrying']['orion_hek293t']:>6d}{100*c['rel_cos_agg']:>+8.3f}{100*c['rel_macro']:>+8.3f}{c['d_pds_A']:>+10.5f}{c['d_pds_ABC']:>+10.5f}"
              f"{c['pds_pos_folds']:>3d}{c['pds_sigma']:>5.1f}{c['d_qw']:>+9.5f}{c['d_nmae_wald89']:>+9.5f}{c['e']:>7.3f}{c['amp_calibrated_pct']:>+7.2f}{'|':>2s}"
              f"{cs['o_pds_A']:>+8.5f}{cs['o_mse']:>+8.5f}{cs['o_de3']:>+8.5f}{cs['o_nmae']:>+9.5f}{cs['total_A']:>+9.5f}{'|':>2s}"
              f"{pt['o_pds_A']:>+8.5f}{pt['o_de3']:>+8.5f}{pt['o_nmae']:>+8.5f}{pt['total_A']:>+9.5f}{100*c['conc'][3]['retained']:>6.0f}%")

    P("\ncolumns: carry = union carrying rows on the hek293t fold (/294); agg%/mac% = rel cos_agg / macro cosine vs cd4svd5; dpds A / ABC = the two "
      "pds functionals; + = folds with d pds > 0; sig = |d pds ABC| / 0.001494; d nmae = exact identity (4 folds, wald89 pool); e = dM/dU (break-even 0.5); "
      "amp% = calibrated delivered amplitude vs cd4svd5; CONS = conservative six-member total (k_pds 0.7, mse at k_cos 0.96 with the amplitude leg, DE-three "
      f"at ARM6's realised {K_DE3_ARM6:.3f}, nmae x1.5); POINT = k_pds {k_pds6_A:.2f} (ARM6-realised, frame A), DE-three 0.24, nmae x1.0; ret3 = leave-top-3-out pds retention")
    table("1. DOSE -- the six-source union replacing the H1 pair, on top of cd4svd5", [f"stack6 g={g:g}" for g in DOSES6] + ["stack6 g=0.45"])
    table("1b. DOSE -- the H1 pair + dld1_lowmoi (229's ARM3), on top of cd4svd5", [f"stack3 g={g:g}" for g in DOSES6] + ["stack3 g=0.45", "stack3 g=0.6"])
    table("1c. DOSE -- the H1 pair alone re-dosed (the incumbent's own union; g = 0.28 is cd4svd5 itself)", ["pair g=0.2", "pair g=0.35", "pair g=0.45"])
    # best dose by the conservative total among stack6 cells
    best6 = max([f"stack6 g={g:g}" for g in DOSES6] + ["stack6 g=0.45"], key=lambda n: cells[n]["conservative"]["total_A"])
    best_all = max(cells, key=lambda n: cells[n]["conservative"]["total_A"])
    P(f"\n  best six-source dose by the conservative total: {best6} ({cells[best6]['conservative']['total_A']:+.5f}); "
      f"best cell overall: {best_all} ({cells[best_all]['conservative']['total_A']:+.5f} conservative, {cells[best_all]['point']['total_A']:+.5f} point)")
    for g in DOSES_SUB:
        table(f"2. SUBSETS at g = {g:g}", [n for n in cells if cells[n]["g"] == g and not n.startswith("stack6 ") and not n.startswith("stack3 ") and not n.startswith("pair g=")]
              + [f"stack3 g={g:g}", f"stack6 g={g:g}"])
    P("\n  the zero row: cd4svd5 itself is the H1 pair at g = 0.28 (every delta above is against it); `pair g=0.20/0.35` re-dose the incumbent's own union.")

    # ---- the best cell in full -------------------------------------------------------------------
    c = cells[best_all]
    P("\n" + "=" * 120)
    P(f"THE BEST CELL: {best_all}  (sources {c['sources']}, g = {c['g']})")
    P("=" * 120)
    P(f"  PLOLO vs cd4svd5: rel cos_agg {100*c['rel_cos_agg']:+.4f}% (per fold {['%+.3f' % (100*v) for v in c['rel_cos_agg_per_fold']]}), rel macro {100*c['rel_macro']:+.4f}%")
    P(f"  d pds frames ABC {c['d_pds_ABC']:+.6f} (per fold {['%+.5f' % v for v in c['d_pds_ABC_per_fold']]}; {c['pds_pos_folds']}/3 positive; "
      f"{c['pds_sigma']:.1f} LOO sigma; fold se {c['pds_fold_se']:.6f}); frame A {c['d_pds_A']:+.6f}")
    P(f"  nmae identity: d wald89 {c['d_nmae_wald89']:+.6f} (3 cancer folds {c['d_nmae_wald89_3c']:+.6f}); dU {c['dU']:+.2f} dM {c['dM']:+.2f} e {c['e']:.3f}; "
      f"d precision89 {c['d_precision89']:+.5f}; d q_w {c['d_qw']:+.5f}")
    P(f"  amplitude vs cd4svd5: effect-level {c['amp_effect_level_pct']:+.3f}%, calibrated {c['amp_calibrated_pct']:+.3f}%, PLOLO median-norm {c['amp_plolo_median_pct']:+.3f}%"
      f"  ({'> 0.2 %: the amplitude leg is charged' if abs(c['amp_calibrated_pct']) > 0.2 else 'within 0.2 %: the amplitude leg is negligible, charged anyway'})")
    P(f"  kolf2 (apart): rel cos_agg {100*c['kolf2']['rel_agg']:+.3f}%, d pds {c['kolf2']['d_pds']:+.5f}")
    P(f"  leave-top-k-out pds retention vs cd4svd5: k=1 {100*c['conc'][1]['retained']:.1f}% (share {100*c['conc'][1]['share']:.1f}%), "
      f"k=3 {100*c['conc'][3]['retained']:.1f}% (share {100*c['conc'][3]['share']:.1f}%), k=5 {100*c['conc'][5]['retained']:.1f}%")
    P(f"  {'factor set':58s}{'d cos':>10s}{'raw mse':>9s}{'mse':>9s}{'pds A':>9s}{'pds ABC':>9s}{'DE3':>9s}{'nmae':>9s}{'TOTAL A':>9s}{'TOTAL ABC':>10s}")
    sens = {}
    for label, kc, kpA, kpABC, de, nx in (
            ("CONSERVATIVE: k_cos 0.96, k_pds 0.7, DE3 ARM6-realised, nmae x1.5", K_COS_NEW, 0.7, 0.7, K_DE3_ARM6, 1.5),
            ("conservative but nmae x1.0", K_COS_NEW, 0.7, 0.7, K_DE3_ARM6, 1.0),
            ("k_pds 1.22 (BR249's re-shape value), DE3 0.136, nmae x1.0", K_COS_NEW, 1.22, 1.22, K_DE3_SHARE, 1.0),
            (f"POINT: k_pds ARM6-realised ({k_pds6_A:.2f} A / {k_pds6_ABC:.2f} ABC), DE3 0.24, nmae x1.0", K_COS_NEW, k_pds6_A, k_pds6_ABC, K_DE3_RECIPE, 1.0),
            ("point but nmae x1.5", K_COS_NEW, k_pds6_A, k_pds6_ABC, K_DE3_RECIPE, 1.5),
            ("point but DE3 ARM6-realised", K_COS_NEW, k_pds6_A, k_pds6_ABC, K_DE3_ARM6, 1.0),
            ("k_cos 1.36 (re-dose class), point otherwise", 1.36, k_pds6_A, k_pds6_ABC, K_DE3_RECIPE, 1.0),
            ("k_cos 0.94 (band low), conservative otherwise", 0.94, 0.7, 0.7, K_DE3_ARM6, 1.5),
            ("nmae at ZERO (an nmae-neutral rendering, BR237 s3 -- NOT this builder), point otherwise", K_COS_NEW, k_pds6_A, k_pds6_ABC, K_DE3_RECIPE, 0.0),
            ("nmae at ZERO, conservative otherwise", K_COS_NEW, 0.7, 0.7, K_DE3_ARM6, 0.0)):
        p = price(c["rel_cos_agg"], c["d_pds_A"], c["d_pds_ABC"], c["d_nmae_wald89"], 1.0 + amp_cal * (amp_ratio(T[best_all]["per"], T["cd4svd5 (file)"]["per"]) - 1),
                  base_per, cos_5, raw_5, kc, kpA, kpABC, de, nx)
        sens[label] = p
        P(f"  {label:58s}{p['dcos']:>+10.6f}{p['raw']:>9.5f}{p['o_mse']:>+9.5f}{p['o_pds_A']:>+9.5f}{p['o_pds_ABC']:>+9.5f}{p['o_de3']:>+9.5f}{p['o_nmae']:>+9.5f}{p['total_A']:>+9.5f}{p['total_ABC']:>+10.5f}")
    P(f"  gate +{GATE} shape-only; ordering-channel noise sd 0.004502 overall (230) -- a margin must be read against that too")

    # ---- verdict ----------------------------------------------------------------------------------
    clears_any = [n for n in cells if cells[n]["clears_conservative"]]
    clears_point = [n for n in cells if cells[n]["clears_point"]]
    P("\n" + "=" * 120)
    P("VERDICT")
    P("=" * 120)
    P(f"  R1 (stack6@0.35 on h1union reproduces six229.json): {'PASS' if r1_pass else 'FAIL'}")
    P(f"  R2 (ARM6's six-member price lands on the board within the bands): {len(r2_pass_sets)}/{len(r2_table)} factor sets pass all four bands")
    P(f"  R3 (cd4svd5 vs h1union at k_cos 3.76 / k_pds 1.22 within +0.0043 +- 0.0008): {'PASS' if r3_pass else 'FAIL'}")
    P(f"  cells clearing +{GATE} over cd4svd5 at the CONSERVATIVE factors: {len(clears_any)}/{len(cells)}" + (f" -> {clears_any}" if clears_any else ""))
    P(f"  cells clearing +{GATE} at the POINT factors: {len(clears_point)}/{len(cells)}" + (f" -> {clears_point}" if clears_point else ""))
    P(f"  best cell {best_all}: conservative {cells[best_all]['conservative']['total_A']:+.5f}, point {cells[best_all]['point']['total_A']:+.5f}")
    verdict = "CLEARS" if clears_any else "DOES NOT CLEAR"
    P(f"\n{verdict}: {'a cell of the six-source union clears the gate at the conservative factors' if clears_any else 'no dose and no source subset of the six-source union clears +0.0019 over cd4svd5 on the validated nmae identity'}")
    out = dict(known_answer=ka, r1_pass=r1_pass, r2_table=r2_table, r2_pass_sets=r2_pass_sets, r3_rows=r3_rows, r3_pass=r3_pass,
               factors=dict(cos_h1union=cos_h, cos_ARM6=cos_6, cos_cd4svd5=cos_5, k_cos_ARM6=k_cos6, k_cos_cd4svd5=k_cos5,
                            k_pds_ARM6_A=k_pds6_A, k_pds_ARM6_ABC=k_pds6_ABC, k_pds_cd4svd5_A=k_pds5_A, k_pds_cd4svd5_ABC=k_pds5_ABC,
                            de3_rate_ARM6=K_DE3_ARM6, de3_rate_cd4svd5=de3_5 / dcos5_board, nmae_ratio_ARM6=nm6_board / nm6_ident,
                            amp_calibration=amp_cal, rel6=rel6, dp6_A=dp6_A, dp6_ABC=dp6_ABC, rel5=rel5, dp5_A=dp5_A, dp5_ABC=dp5_ABC,
                            d_nmae_arm6=d_nmae_arm6, d_nmae_cd4svd5=d_nmae_5),
               aggregates={n: {k: v for k, v in agg[n].items() if k != "per_fold"} for n in names},
               cells=cells, best_cell=best_all, best_dose6=best6, best_cell_sensitivity=sens,
               clears_conservative=clears_any, clears_point=clears_point, verdict=verdict)
    (WORK / "price258.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "price258.log").write_text("\n".join(lines) + "\n")
    P(f"\nwrote {WORK/'price258.json'} and price258.log")


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
