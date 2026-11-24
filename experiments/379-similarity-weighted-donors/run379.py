#!/usr/bin/env python3
"""379 -- similarity-weighted donor consensus on the six-line frame.

Reads 334/341/360's outputs read-only (sys.path import, exactly like 341 imports 334).
Writes only under experiments/379-similarity-weighted-donors/work/. No vcc, no build,
no submission, nothing outside this experiment's own directory.
"""
import os
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
E334 = ROOT / "experiments/334-text-conditioned-transfer"
E341 = ROOT / "experiments/341-known-line-calibration"
E360 = ROOT / "experiments/360-depmap-three-lines"
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)

sys.path.insert(0, str(E334))
sys.path.insert(0, str(E341))
import run334 as B334          # noqa: E402
import core334 as C334         # noqa: E402
import core341 as C341         # noqa: E402

SEED = 3790910
TAUS = (0.05, 0.1, 0.2)
LINE_ACH = {
    "k562_gwps": "ACH-000551",
    "orion_hct116": "ACH-000971",
    "orion_hek293t": None,     # not a DepMap model (0 hits in Model.csv)
    "rpe1": None,              # only engineered RPE1SS* clones in Model.csv
    "nadig_jurkat": "ACH-000995",
    "nadig_hepg2": "ACH-000739",
}


# ---------------------------------------------------------------- frame ----

def load_frame():
    selection = json.loads((E334 / "selection.json").read_text())
    targets, genes = selection["targets"], selection["genes"]
    lines = list(B334.LINES)
    assert lines == list(selection["lines"])
    loaded = [B334.load(b, targets, genes) for b in B334.BANKS]
    raw = np.array([a for a, _ in loaded])
    basal = np.array([b for _, b in loaded])
    direction, rms = C334.rms_direction(raw)
    return selection, targets, genes, lines, direction, basal, rms


def two_blocks(targets):
    order = C334.order(targets, "334-split|")
    return [np.array(order[::2]), np.array(order[1::2])]


# ------------------------------------------------------------ similarity ---

def bank_similarity(basal):
    n = basal.shape[0]
    S = np.eye(n)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            S[i, j] = np.corrcoef(basal[i], basal[j])[0, 1]
    return S


def depmap_similarity(lines, genes_512, ids_field, mat_field, bank_fallback):
    d = np.load(E360 / "work/three_lines.npz", allow_pickle=True)
    ids = list(d[ids_field])
    mat = np.asarray(d[mat_field], dtype=np.float64)
    axis_genes = list(d["genes"])
    pos = {g: i for i, g in enumerate(axis_genes)}
    common_idx = [pos[g] for g in genes_512 if g in pos]
    covered = {name: ids.index(ach) for name, ach in LINE_ACH.items()
               if ach is not None and ach in ids}
    sub = {name: mat[row][common_idx] for name, row in covered.items()}
    if sub:
        finite = np.all([np.isfinite(v) for v in sub.values()], axis=0)
    else:
        finite = np.zeros(0, bool)
    n = len(lines)
    S = bank_fallback.copy()
    covered_idx = {lines.index(name) for name in covered}
    for a, ia in enumerate(lines):
        for b, ib in enumerate(lines):
            if a == b:
                continue
            if lines[a] in covered and lines[b] in covered:
                va, vb = sub[lines[a]][finite], sub[lines[b]][finite]
                S[a, b] = np.corrcoef(va, vb)[0, 1]
    meta = {"ids_field": ids_field, "n_common_genes_axis": len(common_idx),
            "n_common_genes_finite": int(finite.sum()) if len(finite) else 0,
            "covered_lines": sorted(covered), "covered_idx": sorted(covered_idx)}
    return S, meta


# ------------------------------------------------------------- weighting ---

def softmax_weights(sim_row, donors, tau):
    s = np.array([sim_row[d] for d in donors], dtype=np.float64)
    z = s / tau
    z = z - z.max()
    w = np.exp(z)
    w = w / w.sum()
    return dict(zip(donors, w))


def uniform_weights(donors):
    w = 1.0 / len(donors)
    return {d: w for d in donors}


def top2_weights(sim_row, donors):
    ordered = sorted(donors, key=lambda d: -sim_row[d])
    keep = set(ordered[:2])
    return {d: (0.5 if d in keep else 0.0) for d in donors}


def wmean(direction, weights, donors, idx):
    acc = None
    for d in donors:
        w = weights[d]
        if w == 0.0:
            continue
        v = direction[d][idx] * w
        acc = v if acc is None else acc + v
    return acc


def diagonal_gain(xs, ys):
    num = sum((x * y).sum(0) for x, y in zip(xs, ys))
    den = sum((x * x).sum(0) for x, y in zip(xs, ys))
    prior = float(np.median(den))
    return (num + prior) / np.maximum(den + prior, 1e-12)


# --------------------------------------------------------------- arm run ---

def run_arm(direction, blocks, scheme, S=None, tau=None, wrong_map=None, rng=None):
    """scheme in {'uniform','softmax','top2','random'}."""
    n_lines, n_targets, n_genes = direction.shape
    pred_cons = np.full_like(direction, np.nan)
    pred_diag = np.full_like(direction, np.nan)

    def get_weights(recipient, donors):
        if scheme == "uniform":
            return uniform_weights(donors)
        if scheme == "random":
            a = rng.dirichlet(np.ones(len(donors)))
            return dict(zip(donors, a))
        eff = wrong_map[recipient] if wrong_map is not None else recipient
        sim_row = {d: S[eff, d] for d in donors}
        if scheme == "softmax":
            return softmax_weights(sim_row, donors, tau)
        if scheme == "top2":
            return top2_weights(sim_row, donors)
        raise ValueError(scheme)

    for outer in range(n_lines):
        donors = [i for i in range(n_lines) if i != outer]
        for block_id in (0, 1):
            train, test = blocks[1 - block_id], blocks[block_id]
            w_outer = get_weights(outer, donors)
            wbase_test = wmean(direction, w_outer, donors, test)
            xs, ys = [], []
            for R in donors:
                inner_donors = [d for d in donors if d != R]
                w_R = get_weights(R, inner_donors)
                x_R = wmean(direction, w_R, inner_donors, train)
                xs.append(x_R)
                ys.append(direction[R][train])
            gain = diagonal_gain(xs, ys)
            pred_cons[outer, test] = wbase_test
            pred_diag[outer, test] = wbase_test * gain
    return pred_cons, pred_diag


# ---------------------------------------------------------------- scoring --

def evaluate(pred, truth):
    n_lines = truth.shape[0]
    sign_rows, retrieval_rows, cos_rows = [], [], []
    for L in range(n_lines):
        se = C341.sign_endpoints(pred[L], truth[L])
        ce = C334.endpoints(pred[L], truth[L])
        sign_rows.append(se["weighted_sign"])
        retrieval_rows.append(ce["retrieval"])
        cos_rows.append(ce["cosine"])
    return (np.array(sign_rows), np.array(retrieval_rows), np.array(cos_rows))


def pool(rows):
    return float(np.mean(rows))


def main():
    t0 = time.monotonic()
    selection, targets, genes, lines, direction, basal, rms = load_frame()
    blocks = two_blocks(targets)
    n_lines = len(lines)
    print(f"lines={lines}")
    print(f"direction shape={direction.shape} basal shape={basal.shape}")

    # ---------------- known-answer reproduction (rule 3) ----------------
    ka = {}
    pred_cons_u, pred_diag_u = run_arm(direction, blocks, "uniform")
    sign_u, retr_u, cos_u = evaluate(pred_diag_u, direction)
    sign_u0, retr_u0, cos_u0 = evaluate(pred_cons_u, direction)
    ka["source_mean_weighted_sign"] = pool(sign_u0)
    ka["source_mean_retrieval"] = pool(retr_u0)
    ka["crossline_diag_weighted_sign"] = pool(sign_u)
    ka["crossline_diag_retrieval"] = pool(retr_u)
    published = {"source_mean_weighted_sign": 0.589398, "source_mean_retrieval": 0.790095,
                 "crossline_diag_weighted_sign": 0.589398, "crossline_diag_retrieval": 0.794614}
    ka["published"] = published
    ka["max_abs_diff"] = max(abs(ka[k] - published[k]) for k in published)
    print("KNOWN ANSWER:", json.dumps(ka, indent=2))

    # cross-check my hand-rolled uniform diagonal against core334's own builtin,
    # on one fold, as an internal consistency check of the reimplemented formula.
    outer, block_id = 0, 0
    train, test = blocks[1 - block_id], blocks[block_id]
    inner_ref = C334.source_blocks(direction, basal, rms, outer, train)
    gain_ref = C334.diagonal(inner_ref)
    donors = [i for i in range(n_lines) if i != outer]
    xs, ys = [], []
    for R in donors:
        inner_donors = [d for d in donors if d != R]
        w_R = uniform_weights(inner_donors)
        x_R = wmean(direction, w_R, inner_donors, train)
        xs.append(x_R)
        ys.append(direction[R][train])
    gain_mine = diagonal_gain(xs, ys)
    max_diff_gain = float(np.abs(gain_ref - gain_mine).max())
    print(f"internal check: hand-rolled uniform diagonal gain vs core334.diagonal, max|diff|={max_diff_gain:.3e}")
    ka["internal_gain_check_max_abs_diff"] = max_diff_gain

    if ka["max_abs_diff"] > 1e-5 or max_diff_gain > 1e-9:
        (WORK / "known_answer.json").write_text(json.dumps(ka, indent=2))
        print("KNOWN ANSWER MISMATCH -- STOPPING PER PROTOCOL RULE 3")
        sys.exit(1)

    (WORK / "known_answer.json").write_text(json.dumps(ka, indent=2))

    # ---------------------------- similarity sources ----------------------------
    S_bank = bank_similarity(basal)
    S_ccle, meta_ccle = depmap_similarity(lines, genes, "expr_ids", "expr_logtpm1p", S_bank)
    S_chronos, meta_dep = depmap_similarity(lines, genes, "dep_ids", "dep_chronos", S_bank)
    sims = {"bank": S_bank, "ccle": S_ccle, "chronos": S_chronos}
    print("S_bank=\n", np.round(S_bank, 3))
    print("S_ccle=\n", np.round(S_ccle, 3), meta_ccle)
    print("S_chronos=\n", np.round(S_chronos, 3), meta_dep)

    wrong_map = dict(enumerate(C334.derangement(n_lines, SEED + 7)))
    print("wrong-recipient derangement:", {lines[k]: lines[v] for k, v in wrong_map.items()})

    results = {"lines": lines, "known_answer": ka,
               "similarity": {k: v.tolist() for k, v in sims.items()},
               "similarity_meta": {"ccle": meta_ccle, "chronos": meta_dep},
               "wrong_recipient_map": {lines[k]: lines[v] for k, v in wrong_map.items()},
               "arms": {}}

    def score_and_store(name, pred_cons, pred_diag):
        sign_c, retr_c, cos_c = evaluate(pred_cons, direction)
        sign_d, retr_d, cos_d = evaluate(pred_diag, direction)
        delta_sign = sign_d - sign_u
        delta_retr = retr_d - retr_u
        delta_cos = cos_d - cos_u
        interval_sign = C341.interval(delta_sign)
        interval_retr = C341.interval(delta_retr)
        results["arms"][name] = {
            "consensus_only": {"weighted_sign": pool(sign_c), "retrieval": pool(retr_c), "cosine": pool(cos_c)},
            "plus_diagonal": {"weighted_sign": pool(sign_d), "retrieval": pool(retr_d), "cosine": pool(cos_d)},
            "delta_vs_uniform_plus_diagonal": {
                "weighted_sign": interval_sign, "retrieval": interval_retr,
                "cosine_mean": pool(delta_cos)},
        }
        print(f"{name:32s} sign={pool(sign_d):.6f} (d={pool(delta_sign):+.6f}) "
              f"retr={pool(retr_d):.6f} (d={pool(delta_retr):+.6f}) cos={pool(cos_d):.6f}")

    results["arms"]["uniform"] = {
        "consensus_only": {"weighted_sign": pool(sign_u0), "retrieval": pool(retr_u0), "cosine": pool(cos_u0)},
        "plus_diagonal": {"weighted_sign": pool(sign_u), "retrieval": pool(retr_u), "cosine": pool(cos_u)},
        "delta_vs_uniform_plus_diagonal": {"weighted_sign": None, "retrieval": None, "cosine_mean": 0.0},
    }
    print(f"{'uniform':32s} sign={pool(sign_u):.6f} retr={pool(retr_u):.6f} cos={pool(cos_u):.6f}  [KNOWN ANSWER]")

    # ---------------------------- honest arms ----------------------------
    for sname, S in sims.items():
        for tau in TAUS:
            pc, pd_ = run_arm(direction, blocks, "softmax", S=S, tau=tau)
            score_and_store(f"{sname}_tau{tau}", pc, pd_)
        pc, pd_ = run_arm(direction, blocks, "top2", S=S)
        score_and_store(f"{sname}_top2", pc, pd_)

    # ---------------------------- wrong-recipient control ----------------------------
    for sname, S in sims.items():
        for tau in TAUS:
            pc, pd_ = run_arm(direction, blocks, "softmax", S=S, tau=tau, wrong_map=wrong_map)
            score_and_store(f"WRONG_{sname}_tau{tau}", pc, pd_)
        pc, pd_ = run_arm(direction, blocks, "top2", S=S, wrong_map=wrong_map)
        score_and_store(f"WRONG_{sname}_top2", pc, pd_)

    # ---------------------------- random-weight control ----------------------------
    random_signs, random_retr = [], []
    for draw in range(20):
        rng = np.random.default_rng(SEED + 1000 + draw)
        pc, pd_ = run_arm(direction, blocks, "random", rng=rng)
        sign_d, retr_d, cos_d = evaluate(pd_, direction)
        random_signs.append(pool(sign_d))
        random_retr.append(pool(retr_d))
    results["arms"]["random_20draws"] = {
        "weighted_sign_mean": float(np.mean(random_signs)),
        "weighted_sign_min": float(np.min(random_signs)), "weighted_sign_max": float(np.max(random_signs)),
        "retrieval_mean": float(np.mean(random_retr)),
        "retrieval_min": float(np.min(random_retr)), "retrieval_max": float(np.max(random_retr)),
        "delta_sign_vs_uniform_mean": float(np.mean(random_signs) - pool(sign_u)),
    }
    print(f"{'random x20':32s} sign mean={np.mean(random_signs):.6f} "
          f"[{np.min(random_signs):.6f},{np.max(random_signs):.6f}] "
          f"retr mean={np.mean(random_retr):.6f}")

    # ---------------------------- per-line breakdown for the best honest arm ----------------------------
    best_name = max((n for n in results["arms"] if not n.startswith("WRONG_") and n not in ("uniform", "random_20draws")),
                     key=lambda n: results["arms"][n]["delta_vs_uniform_plus_diagonal"]["weighted_sign"]["mean"])
    sname, rest = best_name.split("_", 1)
    if rest.startswith("tau"):
        tau = float(rest[3:])
        pc, pd_ = run_arm(direction, blocks, "softmax", S=sims[sname], tau=tau)
    else:
        pc, pd_ = run_arm(direction, blocks, "top2", S=sims[sname])
    sign_d, retr_d, cos_d = evaluate(pd_, direction)
    per_line = {lines[i]: {"weighted_sign": float(sign_d[i].mean()), "retrieval": float(retr_d[i].mean()),
                            "delta_sign_vs_uniform": float(sign_d[i].mean() - sign_u[i].mean())}
                for i in range(n_lines)}
    results["best_honest_arm"] = best_name
    results["best_honest_arm_per_line"] = per_line
    print("best honest arm:", best_name)
    print(json.dumps(per_line, indent=2))

    results["wall_seconds"] = time.monotonic() - t0
    (WORK / "results379.json").write_text(json.dumps(results, indent=2))
    print(f"done in {results['wall_seconds']:.1f}s")


if __name__ == "__main__":
    main()
