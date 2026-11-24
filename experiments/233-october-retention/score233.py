"""233 Phase 1/2 — all six members, on every (panel, fold, arm), in the challenge frame.

MEMBERS.  Each is the local expression of a board member that 162/188/203/228 already
validated, computed after `profile_from_effect` + the scorer's read-back (156: an effect-level
score measures an object the board never sees).

  cos_agg  energy-weighted aggregate transfer cosine, sqrt(1 - sse)          [DIRECTION]
  pds      exact rank statistic over the fold's decoy set, ties at half      [DIRECTION]
  mse*     `1 - cos^2` = min_alpha expr_mse (AGENTS.md §4.1); credit = cos^2 [DIRECTION]
  nmae     188's exact identity, `|lfc|*sqrt(cpm)` pool at depth 89          [DIRECTION]
  fid      sign precision of the SHIPPED `rc_callset` call set (400 + 2,200);
           the yield term `min(1, n_pred/n_real)` is pinned at 1            [ORDERING]
  reach    the metric's own definition: deepest prefix of the reference's DE
           pool, ordered by our confidence, that stays >= 90% sign-pure, / n_real  [ORDERING]
  jac      |call n ref| / |call u ref| at n_real in {89, 200, 400}          [ORDERING]

DECOY-COUNT CONTROL.  `pds` is a retrieval statistic and its value depends on how many decoys
there are.  The pessimistic panel gives the k562 fold 150 evaluable targets against today's
269, so a raw `pds` ratio would confound the panel with the decoy count.  `pds_at` reports
`pds` re-computed on random subsets of a fixed size, averaged over draws, so panels can be
compared at MATCHED n.  Every retention ratio for `pds` is quoted at matched n.

ARMS (the Phase-2 attribution).  BASE is the shipped stack.  Each ablation removes exactly one
component, so `BASE - ablation` is that component's contribution ON THAT PANEL, and the
question "which coverage change carries the score change" is the comparison of those
contributions across panels.

  BASE      consensus(+PerturbMatch) -> operator -> + CD4 @ f=0.34 + H1 union @ g=0.28
  cd4_off   CD4 dose 0                      (CD4 union coverage 291 -> 300 / 174)
  h1_off    H1 union dose 0                 (gse295214_h1 64 -> 57 / 34; h1_2025 20 -> 7 / 3)
  resid_off both residuals off
  pm_off    no PerturbMatch blend in the K562 slot
  src_<n>   consensus from ONE source only  (prices the second carrier)
  null      BASE with the predicted rows DERANGED -- the empirical zero of every member on
            that panel and fold, so a retention ratio is taken above a measured null and not
            above an assumed one.

    uv run python score233.py --panels abc real pess --folds k562_gwps ...
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/188-magnitude-calibration"))
sys.path.insert(0, str(ROOT / "experiments/203-reopen-125"))
sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))
sys.path.insert(0, str(ROOT / "experiments/228-theta-gain-plolo"))

import panels233 as P  # noqa: E402
import plolo  # noqa: E402
from plolo import E062, Operator, deliver, endpoints  # noqa: E402
from mapfit import to_lfc  # noqa: E402
from phase23 import POOL_N, topn, umd  # noqa: E402

import fold233  # noqa: E402

ALPHA, GAMMA, LAM = 0.982, 0.65, 0.0
CD4_F, UNION_G = 0.34, 0.28
UNION = ("gse295214_h1", "h1_2025")
RC_HEAD, RC_CPM = 400, 2200
N_REALS = (89, 200, 400)
PURITY = 0.9
NULL_SEED = 233233
PDS_MATCH_REPS = 64
CONTEXTS = ("A", "B", "C")
SHARED = P.WORK / "plolo" / "shared"
SHARED.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------------ residual plumbing
def _resid_module():
    """Import 218's residual builders, then take `plolo.WORK` back.

    `resid218` and `compose218` both re-point `plolo.WORK` at 218's directory on import, and
    this experiment must not write outside its own.  `resid218.PANEL_SET` is also bound at
    IMPORT time and is used to drop panel coordinates from every residual -- so it has to be
    re-bound for each panel or the H1 residual would keep zeroing TODAY'S 300 genes while
    scoring a disjoint panel.  `bind_panel` below does that; forgetting it is a silent leak.
    """
    import resid218
    import compose218
    plolo.WORK = SHARED
    return resid218, compose218


def bind_panel(panel: str, resid218) -> dict:
    audit = P.bind(panel)
    resid218.PANEL_SET = set(P.read_panel(panel))
    resid218._PREP_CACHE.clear()
    plolo.WORK = SHARED                     # challenge_bctrl's cache is panel-independent
    return audit


def cd4_prepared(panel: str, targets):
    """Production `prepare_cd4_mean` on the PANEL'S OWN re-extracted banks."""
    banks = [E062.Bank.load(n, P.cd4_root(panel) / n, require_cd4_scale=True)
             for n in ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")]
    prep = E062.prepare_cd4_mean(banks, targets, plolo.OFFICIAL, set(P.read_panel(panel)),
                                 "gene_norm", chunk_rows=256, max_rss_gib=10.0)
    complete = np.asarray(prep.target_condition_count) == 3
    eff = np.asarray(prep.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(prep.gene_votes, dtype=np.uint8).copy()
    eff[~complete] = 0.0
    votes[~complete] = 0
    return E062.CD4Prepared(mode=prep.mode, mean_effect=eff, gene_votes=votes,
                            target_condition_count=np.asarray(
                                prep.target_condition_count, dtype=np.uint8).copy(),
                            condition_details=prep.condition_details)


# ------------------------------------------------------------------ the six members
def call_set(effect_row, gate10, keep, ref_lin):
    score = np.where(gate10 & keep, np.abs(np.nan_to_num(effect_row)), -np.inf)
    head = np.argsort(-score, kind="stable")[:RC_HEAD]
    rest = np.where(keep, ref_lin, -np.inf).copy()
    rest[head] = -np.inf
    cpm = np.argsort(-rest, kind="stable")[:RC_CPM]
    call = np.zeros(len(score), dtype=bool)
    call[head] = True
    call[cpm] = True
    return call


def reach_depth(conf, s_pred, s_true, pool_idx) -> float:
    """`de_wilcoxon_direction_reach_raw` for one row, on `pool_idx` (the reference's DE pool).

    Genes are walked in DESCENDING confidence (the submission's own ranking); `k*` is the
    DEEPEST prefix whose sign purity is >= 0.9 (the spec is explicit that purity is not
    monotone, so a prefix that dips and recovers counts at its deeper crossing).
    """
    order = pool_idx[np.argsort(-conf[pool_idx], kind="stable")]
    ok = (s_pred[order] != 0) & (s_pred[order] == s_true[order])
    adj = s_true[order] != 0
    keep = np.flatnonzero(adj)
    if keep.size == 0:
        return float("nan")
    hits = np.cumsum(ok[keep].astype(np.float64))
    k = np.arange(1, keep.size + 1)
    good = np.flatnonzero(hits / k >= PURITY)
    return float(k[good[-1]] / keep.size) if good.size else 0.0


def pds_at(P_mat, T_mat, n: int, reps: int, rng) -> float:
    """`pds` re-computed on random subsets of size `n` -- the decoy-count control."""
    m = P_mat.shape[0]
    if n >= m:
        return float("nan")
    A = P_mat / np.maximum(np.linalg.norm(P_mat, axis=1, keepdims=True), 1e-300)
    B = T_mat / np.maximum(np.linalg.norm(T_mat, axis=1, keepdims=True), 1e-300)
    S = A @ B.T
    vals = []
    for _ in range(reps):
        idx = rng.choice(m, size=n, replace=False)
        s = S[np.ix_(idx, idx)]
        own = np.diag(s).copy()
        beaten = (s < own[:, None]).sum(1)
        ties = (s == own[:, None]).sum(1) - 1
        vals.append(float(np.mean((beaten + 0.5 * ties) / (n - 1))))
    return float(np.mean(vals))


def score_arm(final, fold, bctx, gate, ctrl, Ts, LR, pools, real_de, rng, match_n,
              full: bool = True):
    """All six members.  `full=False` computes only the DIRECTION members (cos/pds/mse*),
    which is what a dose sweep needs and is ~20x cheaper than the per-row call-set loop."""
    sc = np.asarray(fold["scored"])
    out: dict[str, object] = {}

    # ---- DIRECTION: cos_agg / pds / mse*, a23 render, mean over the three challenge frames
    acc = {}
    for c in CONTEXTS:
        dl = deliver(final, bctx[c], gate[c], "minnorm")
        e = endpoints(dl[:, sc], Ts)
        e["cos_agg"] = float(np.sqrt(max(0.0, 1.0 - e["sse"])))
        acc[c] = e
    for k in ("pds", "cos", "cos_agg", "sse", "sgn", "sgn_w"):
        out[k] = float(np.mean([acc[c][k] for c in CONTEXTS]))
    out["mse_star"] = out["sse"]
    if match_n:
        dlA = deliver(final, bctx["A"], gate["A"], "minnorm")
        for n in match_n:
            out[f"pds_at{n}"] = pds_at(dlA[:, sc], Ts, n, PDS_MATCH_REPS,
                                       np.random.default_rng(9090 + n))

    if not full:
        return out

    # ---- nmae: 188's identity, a15 render, challenge frame A (six228's convention)
    H_A = deliver(final, bctx["A"], gate["A"], None)
    LP = to_lfc(H_A, bctx["A"], gate["A"])[:, sc]
    out["nmae_pools"] = {}
    for k, pool in pools.items():
        u, m, d = umd(LP, LR, pool)
        out["nmae_pools"][k] = float(np.mean(1 + (u - 2 * m) / d))
    out["nmae"] = out["nmae_pools"]["wald89"]

    # ---- ORDERING: the call set, reach depth and jaccard, per challenge context
    fid, reach, jac = [], [], {n: [] for n in N_REALS}
    for c in CONTEXTS:
        keep, ref_lin = ctrl[c]
        ok = tot = 0
        calls = np.zeros((final.shape[0], len(sc)), dtype=bool)
        depths = []
        for i in range(final.shape[0]):
            p = plolo.profile_from_effect(final[i], bctx[c], gate[c], "minnorm")
            h_row = np.log1p(plolo.TS * p) - bctx[c]
            call = call_set(final[i], gate[c], keep, ref_lin)
            cs = call[sc]
            calls[i] = cs
            s_pred = np.sign(h_row[sc])
            s_true = np.sign(Ts[i])
            live = cs & (s_true != 0) & (s_pred != 0)
            ok += int((s_pred[live] == s_true[live]).sum())
            tot += int(live.sum())
            # reach: the reference's DE pool, ordered by our confidence
            conf = np.abs(h_row[sc]) + cs.astype(np.float64) * 1e6
            depths.append(reach_depth(conf, s_pred, s_true,
                                      np.flatnonzero(real_de[POOL_N][i])))
        fid.append(ok / max(tot, 1))
        reach.append(float(np.nanmean(depths)))
        for n in N_REALS:
            inter = (calls & real_de[n]).sum(1)
            union = (calls | real_de[n]).sum(1)
            jac[n].append(float(np.mean(inter / np.maximum(union, 1))))
    out["fid"] = float(np.mean(fid))
    out["reach"] = float(np.mean(reach))
    out["jac"] = {str(n): float(np.mean(v)) for n, v in jac.items()}
    out["n_pred"] = RC_HEAD + RC_CPM
    out["median_row_norm"] = float(np.median(np.linalg.norm(final, axis=1)))
    return out


# ------------------------------------------------------------------ the arms
def build_final(fold, op, pre, cd4, h1, wcols, cd4_f, union_g, compose,
                cd4_sizing: str = "all_rows", audit: dict | None = None):
    """The shipped stack.  `cd4_sizing` selects production's ALL-ROWS median solve or the
    CARRYING-ROW one; DEF-RUNBOOK §5.0 measures the two diverging by ~1.24x at the
    pessimistic draw's coverage, and the all-rows rule returns an exact zero below 50%."""
    resid218, compose218 = compose
    whitened = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    solve = (compose218.solve_fraction_carrying if cd4_sizing == "carrying"
             else plolo.solve_fraction)
    frac_cd4, pu_cd4 = solve(whitened, cd4, cd4_f) if cd4_f else (0.0, 0.0)
    frac_u, pu_u = (compose218.solve_fraction_carrying(whitened, h1, union_g)
                    if union_g else (0.0, 0.0))
    if audit is not None:
        _, pu_cd4_all = plolo.solve_fraction(whitened, cd4, 1.0)
        _, pu_cd4_car = compose218.solve_fraction_carrying(whitened, cd4, 1.0)
        _, pu_u_all = plolo.solve_fraction(whitened, h1, 1.0)
        _, pu_u_car = compose218.solve_fraction_carrying(whitened, h1, 1.0)
        audit.update(frac_cd4=float(frac_cd4), frac_union=float(frac_u),
                     cd4_per_unit_all_rows=float(pu_cd4_all),
                     cd4_per_unit_carrying=float(pu_cd4_car),
                     h1_per_unit_all_rows=float(pu_u_all),
                     h1_per_unit_carrying=float(pu_u_car),
                     cd4_carrying_rows=int((np.asarray(cd4.target_condition_count) > 0).sum()),
                     h1_carrying_rows=int((np.asarray(h1.target_condition_count) > 0).sum()),
                     cd4_sizing=cd4_sizing)
    if cd4_f == 0.0 and union_g == 0.0:
        final = np.asarray(whitened, dtype=np.float64).copy()
    else:
        final, _ = E062.add_post_whiten_residual(
            whitened, compose218.merged(cd4, frac_cd4, h1, frac_u), 1.0, wcols)
    final[:, fold["panel_positions"]] = 0.0
    return final


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", nargs="*", default=["abc", "real", "pess"])
    ap.add_argument("--folds", nargs="*", default=list(P.CANCER_FOLDS))
    ap.add_argument("--arms", nargs="*", default=None)
    ap.add_argument("--out", default=str(P.WORK / "score233.json"))
    args = ap.parse_args()

    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["panel"], r["fold"], r["arm"]) for r in rows}
    compose = _resid_module()
    resid218, compose218 = compose

    ctrl_mod = __import__("ctrl228")
    ctrl = {c: ctrl_mod.control_block(c) for c in CONTEXTS}
    bctx = {c: plolo.challenge_bctrl(c) for c in CONTEXTS}
    gate = {}
    for c in CONTEXTS:
        e = np.expm1(bctx[c])
        gate[c] = (1e6 * e / float(e.sum()) > plolo.MIN_CPM) & (bctx[c] > 0)

    # decoy-count control: compare every panel at the smallest fold size seen
    match_n = [150, 250]

    for panel in args.panels:
        bind_panel(panel, resid218)
        for held in args.folds:
            t0 = time.time()
            fold = fold233.build(panel, held)
            print(f"[{panel}/{held}] operator {fold233.verify_op(panel, held)}", flush=True)
            sc = np.asarray(fold["scored"])
            wcols = np.asarray(fold["wcols"])
            Ts = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
            targets = fold["targets"]
            cd4 = cd4_prepared(panel, targets)
            h1 = compose218.union_rows([resid218.source_residual(n, targets, "gene_norm")
                                        for n in UNION])
            op = Operator(fold["op_path"])

            cpm_sc = 1e6 * np.expm1(fold["b_ctrl"][sc]) / float(
                np.expm1(fold["b_ctrl"]).sum())
            LR = to_lfc(np.asarray(fold["truth"], dtype=np.float64),
                        fold["b_ctrl"], gate["A"])[:, sc]
            pools = {"wald89": topn(np.abs(LR) * np.sqrt(np.maximum(cpm_sc, 0)), POOL_N),
                     "lfc89": topn(np.abs(LR), POOL_N),
                     "h89": topn(np.abs(Ts), POOL_N)}
            real_de = {n: topn(np.abs(Ts), n) for n in N_REALS}
            real_de[POOL_N] = real_de.get(POOL_N, topn(np.abs(Ts), POOL_N))

            specs = {
                "BASE": dict(pm=True, only=None, cd4_f=CD4_F, g=UNION_G),
                "cd4_off": dict(pm=True, only=None, cd4_f=0.0, g=UNION_G),
                "h1_off": dict(pm=True, only=None, cd4_f=CD4_F, g=0.0),
                "resid_off": dict(pm=True, only=None, cd4_f=0.0, g=0.0),
                "pm_off": dict(pm=False, only=None, cd4_f=CD4_F, g=UNION_G),
                "null": dict(pm=True, only=None, cd4_f=CD4_F, g=UNION_G, derange=True),
                "cd4_carrying": dict(pm=True, only=None, cd4_f=CD4_F, g=UNION_G,
                                     cd4_sizing="carrying"),
            }
            for s in fold["sources"]:
                specs[f"src_{s}"] = dict(pm=True, only=s, cd4_f=CD4_F, g=UNION_G)
            names = args.arms or list(specs)

            pre_cache: dict[tuple, np.ndarray] = {}
            for arm in names:
                if (panel, held, arm) in done:
                    continue
                spec = specs[arm]
                key = (spec["pm"], spec["only"])
                if key not in pre_cache:
                    pre_cache[key] = fold233.consensus(fold, pm=spec["pm"],
                                                       only=spec["only"])[0]
                aud: dict = {}
                final = build_final(fold, op, pre_cache[key], cd4, h1, wcols,
                                    spec["cd4_f"], spec["g"], compose,
                                    cd4_sizing=spec.get("cd4_sizing", "all_rows"),
                                    audit=aud)
                if spec.get("derange"):
                    rng = np.random.default_rng(NULL_SEED)
                    n = final.shape[0]
                    perm = rng.permutation(n)
                    for _ in range(50):
                        if not (perm == np.arange(n)).any():
                            break
                        perm = rng.permutation(n)
                    final = final[perm]
                m = score_arm(final, fold, bctx, gate, ctrl, Ts, LR, pools, real_de,
                              np.random.default_rng(7), match_n)
                row = dict(panel=panel, fold=held, arm=arm, n_targets=int(final.shape[0]),
                           audit=fold["audit"], sizing=aud,
                           seconds=round(time.time() - t0, 1), **m)
                rows.append(row)
                out_path.write_text(json.dumps(rows, indent=2))
                print(f"[{panel}/{held}/{arm}] cos {row['cos_agg']:.5f} pds {row['pds']:.5f} "
                      f"nmae {row['nmae']:.5f} fid {row['fid']:.5f} reach {row['reach']:.5f} "
                      f"jac89 {row['jac']['89']:.5f} ({row['seconds']:.0f}s)", flush=True)
                del final
            del op, pre_cache
    out_path.write_text(json.dumps(rows, indent=2))
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
