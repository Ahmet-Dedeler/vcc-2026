"""359 main arms: expression-gated propagation on the six-line source frame.

See PLAN.md. All arms run in 341's split regime (complementary target
halves); gain/baseline identical across arms; only the test-target
prediction differs. Writes work/endpoints359.npz + work/results359.json.
Nothing outside this directory.
"""

import csv
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIOR = ROOT / "experiments/334-text-conditioned-transfer"
EXP341 = ROOT / "experiments/341-known-line-calibration"
EXP340 = ROOT / "experiments/340-neighborhood-transfer"
sys.path.insert(0, str(PRIOR))
sys.path.insert(0, str(EXP341))

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import run334 as B  # noqa: E402
import core334 as P  # noqa: E402
import core341 as C  # noqa: E402

SEED = 3590909
LAMBDAS = (0.25, 0.5, 1.0)
SOFT_S = 0.5


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def blocks_for(direction, basal, rms, recipients, targets):
    recips = list(recipients)
    return [(r, P.features(direction, basal, rms,
                           [i for i in recips if i != r], targets),
             direction[r][targets]) for r in recips]


def build_P(W, node_gate_row, node_gate_col=None):
    """Row-normalized gated propagation. W signed or unsigned (G,G).

    G[i,j] = gr[i]*gc[j]; D_ii = sum_j |W∘G|; isolated rows -> self-loop.
    """
    W = np.asarray(W, dtype=float)
    gr = np.asarray(node_gate_row, dtype=float)
    gc = gr if node_gate_col is None else np.asarray(node_gate_col, dtype=float)
    M = W * (gr[:, None] * gc[None, :])
    denom = np.abs(M).sum(1)
    Pn = np.zeros_like(M)
    ok = denom > 0
    Pn[ok] = M[ok] / denom[ok, None]
    Pn[~ok, :] = 0.0
    np.fill_diagonal(Pn, np.where(~ok, 1.0, np.diag(Pn)))
    return Pn


def main():
    t0 = time.monotonic()
    selection = json.loads((PRIOR / "selection.json").read_text())
    targets, genes = selection["targets"], selection["genes"]
    T, G = len(targets), len(genes)
    lines = list(B.LINES)
    n_lines = len(lines)

    loaded = [B.load(bank, targets, genes) for bank in B.BANKS]
    raw = np.array([a for a, _ in loaded])
    basalG = np.array([b for _, b in loaded])
    del loaded
    # basal for the 1024 target symbols (seeded arm node gates)
    basalT = np.full((n_lines, T), np.nan)
    missing = {}
    for li, bank in enumerate(B.BANKS):
        b_all = np.load(ROOT / bank / "b_ctrl.npy")
        glist = (ROOT / bank / "genes.txt").read_text().splitlines()
        gi = {g: i for i, g in enumerate(glist)}
        miss = 0
        for ti, t in enumerate(targets):
            if t in gi:
                basalT[li, ti] = float(b_all[gi[t]])
            else:
                miss += 1
        missing[lines[li]] = miss
    print("target-basal missing per line:", missing)

    direction, rms = P.rms_direction(raw)
    del raw

    order = P.order(targets, "334-split|")
    blocks = [np.array(order[::2]), np.array(order[1::2])]

    # --- prior adjacencies ---
    sel340 = json.loads((EXP340 / "selection.json").read_text())
    assert sel340["genes"] == genes and sel340["targets"] == targets
    Wdir = np.zeros((G, G))
    for i, nbrs in enumerate(sel340["graph_neighbors"]):
        for j in nbrs:
            if int(j) >= 0 and int(j) != i:
                Wdir[i, int(j)] = 1.0
    Wsym = ((Wdir + Wdir.T) > 0).astype(float)
    print(f"STRING directed edges {int(Wdir.sum())}, symmetrized {int(Wsym.sum())}")

    Wct = np.zeros((G, G))
    gidx = {g: i for i, g in enumerate(genes)}
    tidx = {t: i for i, t in enumerate(targets)}
    txg = []  # (t_idx, g_idx, mor) for the seeded arm
    n_ct = 0
    with open(ROOT / "data/external/regulon/collectri_rescued.csv") as f:
        for row in csv.DictReader(f):
            s, t, mor = row["source"], row["target"], int(row["weight"])
            assert mor in (1, -1)
            n_ct += 1
            if s in gidx and t in gidx and s != t:
                # smoothing sign: co-response (+mor); duplicates: sum then sign
                Wct[gidx[t], gidx[s]] += mor
            if s in tidx and t in gidx:
                txg.append((tidx[s], gidx[t], mor))
    Wct = np.sign(Wct)
    cov_t = sorted(set(t for t, _, _ in txg))
    print(f"CollecTRI rows {n_ct}; GxG edges {int(np.abs(Wct).sum())} "
          f"genes w/ in-edge {int((np.abs(Wct).sum(1) > 0).sum())}; "
          f"TxG edges {len(txg)} covering {len(cov_t)}/1024 targets")

    # --- gates ---
    tau = np.median(basalG, axis=1)
    iqr = np.quantile(basalG, 0.75, axis=1) - np.quantile(basalG, 0.25, axis=1)
    hard = basalG >= tau[:, None]
    soft = sigmoid((basalG - tau[:, None]) / SOFT_S)
    rng = np.random.default_rng(SEED)
    rands = [np.array([rng.permutation(h) for h in hard]) for _ in range(3)]
    tauT = np.nanmedian(basalT, axis=1)
    hardT = np.where(np.isnan(basalT), True, basalT >= tauT[:, None])
    print("line tau:", [round(float(v), 4) for v in tau], " IQR:",
          [round(float(v), 4) for v in iqr])
    print("hard-gate density per line:",
          [round(float(h.mean()), 4) for h in hard])

    Pcache = {}

    def get_P(Wkey, W, gate_vec):
        key = (Wkey, gate_vec.tobytes())
        if key not in Pcache:
            Pcache[key] = build_P(W, gate_vec)
        return Pcache[key]

    ones = np.ones(G)
    # gate variants per recipient line L: list of (gkey, gatevec)
    def gate_variants(L):
        out = [("hard", hard[L].astype(float)), ("soft", soft[L]),
               ("ungated", ones)]
        for s in range(3):
            out.append((f"rand{s}", rands[s][L].astype(float)))
        for m in range(n_lines):
            if m != L:
                out.append((f"wrong-{lines[m]}", hard[m].astype(float)))
        return out

    # --- donor-only lambda selection (train targets, inner LOLO) ---
    Wprim = Wdir
    configs_sel = {"hard.s1": ("hard", 1), "soft.s1": ("soft", 1),
                   "hard.s2": ("hard", 2), "ungated.s1": ("ungated", 1)}
    chosen = {}
    for cfg, (gk, steps) in configs_sel.items():
        lam_score = {l: [] for l in LAMBDAS}
        for L in range(n_lines):
            donors = [i for i in range(n_lines) if i != L]
            for d in donors:
                others = [i for i in donors if i != d]
                base_d = direction[others].mean(0)  # (T,G)
                for block in (0, 1):
                    train = blocks[1 - block]
                    ib = blocks_for(direction, basalG, rms, others, train)
                    gain_d = P.diagonal(ib)
                    s0 = base_d[train] * gain_d
                    gv = {"hard": hard[d].astype(float), "soft": soft[d],
                          "ungated": ones}[gk]
                    Pm = build_P(Wprim, gv)
                    sk = s0 @ Pm.T
                    if steps == 2:
                        sk = sk @ Pm.T
                    truth = direction[d][train]
                    for l in LAMBDAS:
                        pred = (1 - l) * s0 + l * sk
                        lam_score[l].append(
                            float(C.sign_endpoints(pred, truth)["weighted_sign"].mean()))
        mean = {l: float(np.mean(v)) for l, v in lam_score.items()}
        best = max(LAMBDAS, key=lambda l: mean[l])
        chosen[cfg] = {"lambda": best, "means": mean}
        print(f"lambda-select {cfg}: " +
              " ".join(f"{l}:{mean[l]:.6f}" for l in LAMBDAS) + f" -> {best}")

    # --- arm prediction ---
    # We build per-line arm lists because wrong-recipient names depend on L.
    # wrong-<X> is undefined on line X itself (left NaN; pooled with nanmean,
    # best-wrong computed per line over its 5 valid sources).
    def arms_for_line(L):
        arms = [("diagonal", None, None, None, 0, 0.0)]
        glist = gate_variants(L)
        gmap = dict(glist)
        for gk, _ in glist:
            for lam in LAMBDAS:
                arms.append((f"str.dir.{gk}.s1.l{lam:g}", "dir", Wdir, gmap[gk], 1, lam))
        for gk in ("hard", "soft", "ungated"):
            for lam in LAMBDAS:
                arms.append((f"str.dir.{gk}.s2.l{lam:g}", "dir", Wdir, gmap[gk], 2, lam))
        for gk in ("hard", "ungated"):
            for lam in LAMBDAS:
                arms.append((f"str.sym.{gk}.s1.l{lam:g}", "sym", Wsym, gmap[gk], 1, lam))
        for gk in [g for g, _ in glist if g == "hard" or g == "ungated" or g.startswith("wrong")]:
            for lam in LAMBDAS:
                arms.append((f"ct.{gk}.s1.l{lam:g}", "ct", Wct, gmap[gk], 1, lam))
        return arms

    # union of arm names across lines (wrong-* names differ per line; unify)
    all_names = ["diagonal"]
    for gk in (["hard", "soft", "ungated"] + [f"rand{s}" for s in range(3)]
               + [f"wrong-{l}" for l in lines]):
        for lam in LAMBDAS:
            all_names.append(f"str.dir.{gk}.s1.l{lam:g}")
    for gk in ("hard", "soft", "ungated"):
        for lam in LAMBDAS:
            all_names.append(f"str.dir.{gk}.s2.l{lam:g}")
    for gk in ("hard", "ungated"):
        for lam in LAMBDAS:
            all_names.append(f"str.sym.{gk}.s1.l{lam:g}")
    for gk in (["hard", "ungated"] + [f"wrong-{l}" for l in lines]):
        for lam in LAMBDAS:
            all_names.append(f"ct.{gk}.s1.l{lam:g}")
    all_names.append("seeded-ct")
    A = len(all_names)
    aidx = {n: i for i, n in enumerate(all_names)}
    Ws = np.full((A, n_lines, T), np.nan)
    Wt = np.full((A, n_lines, T), np.nan)
    Wr = np.full((A, n_lines, T), np.nan)

    for L in range(n_lines):
        donors = [i for i in range(n_lines) if i != L]
        baseline = direction[donors].mean(0)
        truth = direction[L]
        stitched = {n: np.full((T, G), np.nan) for n in all_names}
        for block in (0, 1):
            train, test = blocks[1 - block], blocks[block]
            inner = P.source_blocks(direction, basalG, rms, L, train)
            gain = P.diagonal(inner)
            s0 = baseline[test] * gain
            stitched["diagonal"][test] = s0
            for name, Wkey, W, gv, steps, lam in arms_for_line(L):
                if Wkey is None:  # diagonal
                    stitched[name][test] = s0
                    continue
                Pmat = get_P(Wkey, W, gv)
                if Wkey == "ct" and not np.abs(W).sum():
                    pred = s0
                else:
                    sk = s0 @ Pmat.T
                    if steps == 2:
                        sk = sk @ Pmat.T
                    pred = (1 - lam) * s0 + lam * sk
                # map per-line wrong names to unified names
                uname = name
                if ".wrong-" in name and Wkey == "dir":
                    # name already carries the source line: str.dir.wrong-<src>.s1.l..
                    uname = name
                stitched[uname][test] = pred
            # seeded-ct on test targets (fixed lam 0.5)
            seed = np.array(stitched["diagonal"][test])
            pred = seed.copy()
            ghard = hard[L].astype(float)
            for tpos, t in enumerate(test):
                if t in cov_t:
                    u = np.zeros(G)
                    for (tt, g, mor) in txg:
                        if tt == t:
                            u[g] += -mor * -1.0 * (float(hardT[L, t]) * ghard[g])
                    # +mor*gate (knockdown seed -1 through -mor edge)
                    nrm = np.sqrt(np.mean(seed[tpos] ** 2)) / max(
                        np.sqrt(np.mean(u ** 2)), 1e-12)
                    u = u * nrm
                    pred[tpos] = 0.5 * seed[tpos] + 0.5 * u
            stitched["seeded-ct"][test] = pred
        for name in all_names:
            full = stitched[name]
            if not np.isfinite(full).all():
                continue  # wrong-<X> on its own line X: undefined by design
            Ws[aidx[name], L] = C.sign_endpoints(full, truth)["weighted_sign"]
            Wt[aidx[name], L] = C.sign_endpoints(full, truth)["top64_weighted_sign"]
            Wr[aidx[name], L] = P.endpoints(full, truth)["retrieval"]
        print(f"{lines[L]} done ({time.monotonic() - t0:.0f}s)", flush=True)

    (HERE / "work").mkdir(parents=True, exist_ok=True)
    np.savez_compressed(HERE / "work" / "endpoints359.npz",
                        arms=np.array(all_names), lines=np.array(lines),
                        wsign=Ws, top64=Wt, retrieval=Wr)
    # summary + intervals for key contrasts (NaN-aware: wrong-<X> valid off X)
    diag = aidx["diagonal"]
    summ = {}
    with np.errstate(invalid="ignore"):
        for name in all_names:
            i = aidx[name]
            valid = np.isfinite(Ws[i]).all(1)
            d = Ws[i][valid] - Ws[diag][valid]
            iv = C.interval(d)
            summ[name] = {"pooled_wsign": float(np.nanmean(Ws[i])),
                          "valid_lines": int(valid.sum()),
                          "delta_vs_diag": iv["mean"],
                          "target_ci95": iv["target_ci95"],
                          "crossed_ci95": iv["crossed_ci95"],
                          "line_means_delta": iv["line_means"],
                          "positive_lines": iv["positive_lines"],
                          "pooled_top64": float(np.nanmean(Wt[i])),
                          "pooled_retrieval": float(np.nanmean(Wr[i]))}
    # best-wrong controls at matched lambda (per line, best of its 5 valid
    # sources; pooled over all 6 lines x 1024 targets)
    for Wtag in ("str.dir", "ct"):
        for lam in LAMBDAS:
            stack = np.array([Ws[aidx[f"{Wtag}.wrong-{m}.s1.l{lam:g}"]]
                              for m in lines])  # (6 src, 6 lines, T)
            best = np.nanmax(stack, axis=0)  # own-line slot is NaN -> ignored
            d = best - Ws[diag]
            iv = C.interval(d)
            key = f"{Wtag}.WRONGBEST.s1.l{lam:g}"
            summ[key] = {"pooled_wsign": float(best.mean()),
                         "valid_lines": 6,
                         "delta_vs_diag": iv["mean"],
                         "target_ci95": iv["target_ci95"],
                         "crossed_ci95": iv["crossed_ci95"],
                         "line_means_delta": iv["line_means"],
                         "positive_lines": iv["positive_lines"],
                         "pooled_top64": float("nan"),
                         "pooled_retrieval": float("nan")}
    out = {"arms": all_names, "chosen_lambda": chosen,
           "diagonal_pooled": {"wsign": float(np.nanmean(Ws[diag])),
                               "top64": float(np.nanmean(Wt[diag])),
                               "retrieval": float(np.nanmean(Wr[diag]))},
           "summary": summ,
           "seconds": time.monotonic() - t0}
    (HERE / "work" / "results359.json").write_text(json.dumps(out, indent=2))
    # headline print: primary + controls at chosen lambda
    ch = chosen["hard.s1"]["lambda"]
    print(f"\nchosen primary lambda = {ch:g}")
    for name in ([f"str.dir.hard.s1.l{ch:g}", f"str.dir.soft.s1.l{ch:g}",
                  f"str.dir.ungated.s1.l{ch:g}"]
                 + [f"str.dir.rand{s}.s1.l{ch:g}" for s in range(3)]
                 + [f"str.dir.WRONGBEST.s1.l{ch:g}"]
                 + [f"ct.hard.s1.l{ch:g}", f"ct.ungated.s1.l{ch:g}",
                    f"ct.WRONGBEST.s1.l{ch:g}", "seeded-ct", "diagonal"]):
        s = summ[name]
        print(f"{name:42s} wsign {s['pooled_wsign']:.6f} "
              f"d {s['delta_vs_diag']:+.6f} [{s['target_ci95'][0]:+.5f},{s['target_ci95'][1]:+.5f}] "
              f"retr {s['pooled_retrieval']:.6f}")


if __name__ == "__main__":
    main()
