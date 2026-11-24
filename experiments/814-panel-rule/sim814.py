"""814 -- simulate candidate draw rules v2 from the feature table (no data beyond work/run1/features814.csv.gz) and report
A/B/C's percentile on every panel statistic among 1,000 redraws per rule (seeds 814000 + i; never 800107/800108).

The draw is draw801's: cells = (coverage class) x (CD4 Rest strength bin, A/B/C quintile edges), A/B/C's quota per cell,
cells in sorted order, rng.choice over the cell's sorted members. A v2 rule changes only
  pool    : the arc800 null pool (feat814: pool0 minus the sealed L7/L8 lists), minus EXCL, minus L1..L6 if --fresh;
  weights : within a cell, members are drawn with probability proportional to w(class) = A/B/C share / pool share of the
            target's class (rng.choice(..., p=w / sum w)), so a panel's expected class mix equals A/B/C's.
RULES (frozen 2026-10-16 before this script first ran; each is one try):
  R0    arc800 (KNOWN ANSWER: must reproduce feat814 run1's null814-R0.npz statistics exactly; abort otherwise)
  R1    R0 - essential screens (RPE1 u Jurkat)                                     (= feat814's R1; KA: matches its npz)
  V2r   R1 + weights on resp5 (production-bank responsiveness quintile, flags814 --classes)
  V2rk  R1 + weights on resp5 x kd2 (responsiveness quintile x weak-knockdown class)
  V2p   R1 + weights on prop10 (pool decile of flags814's six-feature propensity score; appended 2026-10-16 ~11:35Z after
        R1/V2r/V2rk had run: try 5)
  V2pc  R1 + per-target weights exp(s) from the same propensity score s (the odds, i.e. the density ratio A/B/C : pool),
        capped at the pool's 99th percentile (appended ~11:40Z after V2p had run: try 6)
--seq6: also simulate the validation draw itself (six panels in order, seeds 814101..814106, disjoint, deficits moved
by draw801's order) and read each panel's PRIMARY statistics against the rule's null.
--indep6: the proposed validation design instead: six INDEPENDENT draws (seeds 814101..814106) from the rule's pool;
each panel's PRIMARY percentiles, its overlap with L1..L6 and the pairwise overlaps (no panel list is written, only sha256).
--fresh: the pool also drops L1..L6 (the v2 validation folds would be disjoint from the worn L panels); the script reports
whether six disjoint v2 panels fit (per-cell capacity) in that pool.
PRIMARY statistics (named before the run; the v2 bar is A/B/C inside [0.1, 0.9] on each): in_essscreen:count,
prod_logcnt_mean:med, prod_logcnt_mean:q90, k562_own_lfc:med, hct_own_lfc:med, pan_cov:count, cd4_guide_anysig:count,
in_vip:count. Every other statistic is reported as the count outside [0.1, 0.9] for A/B/C and for each L panel (the
L panels show how many a typical draw has).
    uv run python experiments/814-panel-rule/sim814.py --rules R0,R1,V2r,V2rk [--fresh] --n 1000 --out <dir>
Light: < 1 GB, < 1 min.
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import argparse  # noqa: E402
import collections  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import feat814 as FT  # noqa: E402

RUN1 = HERE / "work/run1"
FLAGS = HERE / "work/flags"
SEED0 = 814000
PRIMARY = ["in_essscreen:count", "prod_logcnt_mean:med", "prod_logcnt_mean:q90", "k562_own_lfc:med", "hct_own_lfc:med",
           "pan_cov:count", "cd4_guide_anysig:count", "in_vip:count"]
CLASSES = ("C3", "K2", "O2", "O1")


def sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read_classes(p: Path) -> dict[str, str]:
    return dict(ln.split("\t") for ln in p.read_text().splitlines() if ln.strip())


def main(rules, n, out: Path, fresh: bool, seq6: bool = False, indep6: bool = False):
    t0 = time.time()
    out.mkdir(parents=True, exist_ok=True)
    F = pd.read_csv(RUN1 / "features814.csv.gz", index_col=0, float_precision="round_trip")   # bit-exact floats
    genes = list(F.index)
    row = {g: i for i, g in enumerate(genes)}
    abc = list(F.index[F.in_abc > 0])
    Ls = {f"L{i}": sorted(F.index[F[f"in_L{i}"] > 0]) for i in range(1, 7)}
    # cells exactly as draw801 (coverage from the three donors' perts lists = k562/hct/hek cov; CD4 Rest n_downstream bins)
    nd = F["cd4_Rest_nds"]
    edges = np.quantile(np.asarray(nd[F.in_abc > 0].dropna(), float), [0.2, 0.4, 0.6, 0.8])
    run1 = json.loads((RUN1 / "feat814.json").read_text())

    def cell(g):
        k, h, e = (F.at[g, c] > 0 for c in ("k562_cov", "hct_cov", "hek_cov"))
        cov = ("C3" if (h and e) else "K2") if k else ("O2" if (h and e) else ("O1" if (h or e) else None))
        x = nd[g]
        return (cov, "NA" if not np.isfinite(x) else f"Q{int(np.digitize(x, edges)) + 1}")

    cells = {g: cell(g) for g in genes}
    quota = dict(collections.Counter(cells[g] for g in Ls["L1"]))   # arc800's quotas (A/B/C's counts after draw801's moves)
    null_pool = sorted(F.index[F.in_null_pool > 0])
    ess = set(F.index[F.in_essscreen > 0])
    Lall = set().union(*map(set, Ls.values()))
    cls = {k: read_classes(FLAGS / f"{k}.tsv") for k in ("resp5", "resp5kd", "prop10")}
    propc = {g: float(v) for g, v in read_classes(FLAGS / "propc.tsv").items()}
    spec = {"R0": (set(), None), "R1": (ess, None), "V2r": (ess, "resp5"), "V2rk": (ess, "resp5kd"), "V2p": (ess, "prop10"),
            "V2pc": (ess, "propc")}
    flags, cont = FT.stat_defs(F.drop(columns=["in_null_pool", "in_abc"] + [f"in_L{i}" for i in range(1, 7)]))
    names = [f"{c}:count" for c in flags] + [f"{c}:{s}" for s in ("med", "q90", "q10", "mean") for c in cont]
    assert names == run1["stat_names"], "feature/statistic layout differs from run1"
    M = F[flags + cont].to_numpy(np.float64)
    fi, ci = np.arange(len(flags)), np.arange(len(flags), len(flags) + len(cont))
    ref = {"ABC": abc, **Ls}
    ref_stats = {k: FT.panel_stats(M, fi, ci, np.asarray([[row[t] for t in v]]))[0] for k, v in ref.items()}
    res = dict(script_sha256=sha256(Path(__file__)), features_sha256=sha256(RUN1 / "features814.csv.gz"),
               class_files={k: sha256(FLAGS / f"{k}.tsv") for k in cls}, n=n, seeds=[SEED0, SEED0 + n - 1], fresh=fresh,
               primary=PRIMARY, rules={})
    for R in rules:
        excl, wname = spec[R]
        pool = [g for g in null_pool if g not in excl and not (fresh and g in Lall)]
        byc = collections.defaultdict(list)
        for g in pool:
            byc[cells[g]].append(g)
        short = {f"{k[0]}|{k[1]}": (q, len(byc[k])) for k, q in quota.items() if len(byc[k]) < q}
        cap6 = {f"{k[0]}|{k[1]}": (q, len(byc[k])) for k, q in quota.items() if len(byc[k]) < 6 * q}
        if short:
            print(f"== {R}: cells short of their quota {short}: not drawable as specified", flush=True)
            res["rules"][R] = dict(pool=len(pool), short=short)
            continue
        W = None
        if wname == "propc":                                   # per-target odds, capped at the pool's 99th percentile
            raw = np.exp(np.asarray([propc[g] for g in pool]))
            cap = float(np.quantile(raw, 0.99))
            C = {g: g for g in pool}                           # every target is its own class
            W = {g: min(float(np.exp(propc[g])), cap) for g in pool}
        elif wname:
            C = cls[wname]
            pa = collections.Counter(C[g] for g in abc)
            pp = collections.Counter(C[g] for g in pool)
            W = {c: (pa.get(c, 0) / len(abc)) / (pp[c] / len(pool)) for c in pp}
        members = {k: np.asarray(sorted(byc[k])) for k in quota}
        probs = {}
        if W:
            for k, m in members.items():
                w = np.asarray([W[C[g]] for g in m], float)
                probs[k] = w / w.sum() if w.sum() > 0 else None
        idx = np.empty((n, 300), np.int64)
        for i in range(n):
            rng = np.random.default_rng(SEED0 + i)
            panel = []
            for k in sorted(quota):
                p = probs.get(k) if W else None
                panel += rng.choice(members[k], size=quota[k], replace=False, p=p).tolist()
            assert len(set(panel)) == 300
            idx[i] = sorted(row[g] for g in panel)            # draw801 sorts the panel; statistics are order-free
        S = FT.panel_stats(M, fi, ci, idx)
        if R in ("R0", "R1") and not fresh:                  # KNOWN ANSWER: feat814's null for the same rule and seeds
            z = np.load(RUN1 / f"null814-{R}.npz", allow_pickle=True)
            same = np.array_equal(np.nan_to_num(z["stats"][:n], nan=-9e99), np.nan_to_num(S, nan=-9e99))
            assert same, f"KA FAIL: {R} null statistics differ from feat814 run1"
            print(f"KA PASS: {R} reproduces feat814 run1's null statistics bit-equal", flush=True)
        pct = {}
        for k, v in ref_stats.items():
            lt, eq = np.sum(S < v[None, :], 0), np.sum(S == v[None, :], 0)
            pct[k] = np.where(np.isfinite(v), (lt + 0.5 * eq) / n, np.nan)
        outside = {k: int(np.sum((pct[k] < 0.1) | (pct[k] > 0.9))) for k in pct}
        prim = {nm: float(pct["ABC"][names.index(nm)]) for nm in PRIMARY}
        res["rules"][R] = dict(pool=len(pool), weights=(W if wname != "propc" else dict(cap=cap, n=len(W))), cells_short_for_6_disjoint=cap6, primary_abc_pct=prim,
                               outside_01_09=outside, n_stats=int(np.sum(np.isfinite(pct["ABC"]))),
                               abc_pct=dict(zip(names, pct["ABC"].tolist())),
                               L_pct={k: dict(zip(names, pct[k].tolist())) for k in Ls},
                               null_mean=dict(zip(names, np.nanmean(S, 0).tolist())), abc_value=dict(zip(names, ref_stats["ABC"].tolist())))
        np.savez_compressed(out / f"sim814-{R}{'-fresh' if fresh else ''}.npz", stats=S, names=np.asarray(names))
        print(f"== {R}{' (fresh pool)' if fresh else ''}: pool {len(pool)}; cells short for 6 disjoint panels: {cap6 or 'none'}")
        if W and wname != "propc":
            print("   weights " + " ".join(f"{c}:{w:.2f}" for c, w in sorted(W.items())))
        print("   PRIMARY  " + "  ".join(f"{nm} {p:.3f}" for nm, p in prim.items()))
        print(f"   stats outside [0.1, 0.9]: ABC {outside['ABC']} | L " + " ".join(f"{outside[k]}" for k in Ls)
              + f"  (of {res['rules'][R]['n_stats']})")
        ext = sorted(((abs(p - 0.5), nm, p) for nm, p in zip(names, pct["ABC"]) if np.isfinite(p)), reverse=True)[:12]
        print("   most extreme: " + "; ".join(f"{nm} {p:.3f}" for _, nm, p in ext), flush=True)
    def wfun(R, wname, g):
        if wname == "propc":
            return min(float(np.exp(propc[g])), res["rules"][R]["weights"]["cap"])
        return res["rules"][R]["weights"][cls[wname][g]]

    if seq6:
        # the validation draw itself, simulated: six panels drawn in order (seeds 814101..814106), each disjoint from the
        # earlier ones, from the rule's pool with its weights; a cell short of its quota takes what is left and moves the
        # deficit by draw801's order (sibling coverage class, then the next strength bin, stronger first). Each panel's
        # PRIMARY statistics are read against the rule's 1,000-draw null (single-panel draws), so drift shows.
        for R in rules:
            if R not in res["rules"] or "short" in res["rules"][R]:
                continue
            excl, wname = spec[R]
            S = np.load(out / f"sim814-{R}{'-fresh' if fresh else ''}.npz")["stats"]
            used: set[str] = set()
            seqrec = {}
            for j, seed in enumerate(range(814101, 814107)):
                pool = [g for g in null_pool if g not in excl and g not in used and not (fresh and g in Lall)]
                byc = collections.defaultdict(list)
                for g in pool:
                    byc[cells[g]].append(g)
                q = dict(quota)
                moves = 0
                for k in sorted(quota):
                    d = q[k] - len(byc[k])
                    if d <= 0:
                        continue
                    q[k] = len(byc[k])
                    c, b = k
                    order = [("K2" if c == "C3" else "C3" if c == "K2" else "O1" if c == "O2" else "O2", b)]
                    if b != "NA":
                        i = int(b[1]) - 1
                        order += [(c, f"Q{x + 1}") for x in list(range(i + 1, 5)) + list(range(i - 1, -1, -1))]
                    for _ in range(d):
                        dest = next(o for o in order if q.get(o, 0) < len(byc[o]))
                        q[dest] = q.get(dest, 0) + 1
                        moves += 1
                rng = np.random.default_rng(seed)
                panel = []
                for k in sorted(q):
                    if not q[k]:
                        continue
                    m = np.asarray(sorted(byc[k]))
                    p = None
                    if wname:
                        w = np.asarray([wfun(R, wname, g) for g in m], float)
                        p = w / w.sum()
                    if p is not None and int(np.sum(p > 0)) < q[k]:
                        panel = None                           # the weighted cell is exhausted: this draw is infeasible
                        break
                    panel += rng.choice(m, size=q[k], replace=False, p=p).tolist()
                if panel is None:
                    seqrec[f"V{j + 1}"] = dict(seed=seed, infeasible=True)
                    print(f"   seq {R} V{j + 1} seed {seed}: INFEASIBLE (a weighted cell has fewer non-zero-weight members than its quota)")
                    break
                assert len(set(panel)) == 300 and not (set(panel) & used)
                used |= set(panel)
                v = FT.panel_stats(M, fi, ci, np.asarray([[row[t] for t in sorted(panel)]]))[0]
                pc = {nm: float((np.sum(S[:, names.index(nm)] < v[names.index(nm)]) + 0.5 * np.sum(S[:, names.index(nm)] == v[names.index(nm)])) / S.shape[0])
                      for nm in PRIMARY}
                abcv = ref_stats["ABC"]
                seqrec[f"V{j + 1}"] = dict(seed=seed, moves=moves, pool=len(pool), primary_pct=pc,
                                         abc_minus_panel={nm: float(abcv[names.index(nm)] - v[names.index(nm)]) for nm in PRIMARY})
                print(f"   seq {R} V{j + 1} seed {seed}: pool {len(pool)}, moves {moves}; PRIMARY pct " +
                      " ".join(f"{pc[nm]:.2f}" for nm in PRIMARY), flush=True)
            res["rules"][R]["seq6"] = seqrec
    if indep6:
        # the proposed validation design: six INDEPENDENT draws (seeds 814101..814106) by the rule, each from the full rule
        # pool (overlap between v2 panels allowed and counted; overlap with L1..L6 counted; L7/L8 never in the pool)
        for R in rules:
            if R not in res["rules"] or "short" in res["rules"][R]:
                continue
            excl, wname = spec[R]
            S = np.load(out / f"sim814-{R}{'-fresh' if fresh else ''}.npz")["stats"]
            pool = [g for g in null_pool if g not in excl and not (fresh and g in Lall)]
            byc = collections.defaultdict(list)
            for g in pool:
                byc[cells[g]].append(g)
            members = {k: np.asarray(sorted(byc[k])) for k in quota}
            panels = {}
            for j, seed in enumerate(range(814101, 814107)):
                rng = np.random.default_rng(seed)
                panel = []
                for k in sorted(quota):
                    p = None
                    if wname:
                        w = np.asarray([wfun(R, wname, g) for g in members[k]], float)
                        p = w / w.sum()
                    panel += rng.choice(members[k], size=quota[k], replace=False, p=p).tolist()
                panels[f"V{j + 1}"] = sorted(panel)
            rec = {}
            for key, panel in panels.items():
                v = FT.panel_stats(M, fi, ci, np.asarray([[row[t] for t in panel]]))[0]
                pc = {nm: float((np.sum(S[:, names.index(nm)] < v[names.index(nm)]) + 0.5 * np.sum(S[:, names.index(nm)] == v[names.index(nm)])) / S.shape[0])
                      for nm in PRIMARY}
                rec[key] = dict(primary_pct=pc, overlap_L=int(len(set(panel) & Lall)),
                                sha256=hashlib.sha256("\n".join(panel).encode()).hexdigest())
                print(f"   indep {R} {key}: overlap with L1-L6 {rec[key]['overlap_L']}; PRIMARY pct " + " ".join(f"{pc[nm]:.2f}" for nm in PRIMARY))
            ov = [len(set(panels[a]) & set(panels[b])) for i, a in enumerate(panels) for b in list(panels)[i + 1:]]
            print(f"   indep {R}: pairwise overlap between v2 panels mean {np.mean(ov):.1f} (min {min(ov)}, max {max(ov)}) of 300")
            res["rules"][R]["indep6"] = dict(panels=rec, pairwise_overlap=ov)
    res["seconds"] = round(time.time() - t0, 1)
    (out / f"sim814{'-fresh' if fresh else ''}.json").write_text(json.dumps(res, indent=1, default=float))
    print(f"done {res['seconds']} s")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--rules", default="R0,R1,V2r,V2rk")
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--fresh", action="store_true")
    ap.add_argument("--indep6", action="store_true", help="also simulate six independent validation draws")
    ap.add_argument("--seq6", action="store_true", help="also simulate the six sequential disjoint validation draws")
    ap.add_argument("--out", default=str(HERE / "work" / "sim1"))
    a = ap.parse_args()
    main(a.rules.split(","), a.n, Path(a.out), a.fresh, a.seq6, a.indep6)
