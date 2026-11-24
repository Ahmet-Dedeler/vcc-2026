"""233 Phase 2 — the causal attribution: is the October loss COVERAGE, or the GENE SET?

Comparing BASE on two panels tells you how much moves.  It does not tell you WHY, because two
things change at once: the sources cover fewer of the new targets (`k562_gwps` 269 -> 150, CD4
union 291 -> 174, the H1 union 84 -> 37), and the targets are different genes with a different
response distribution.

This separates them.  On TODAY'S panel, where the gene set is held fixed, each component's
coverage is THINNED at random down to the rate the pessimistic panel produces, and the stack
is re-scored:

    thin_k562   the K562 slot is marked unavailable on random targets until the fold's
                two-carrier fraction matches the pessimistic fold's
    thin_cd4    CD4 rows dropped at random to the pessimistic carrying-row count
    thin_h1     H1-union rows dropped likewise
    thin_all    all three together

If `thin_all` on today's panel lands on the pessimistic panel's level, **the loss is coverage**
and is fixable by adding sources.  If today's panel stays materially higher after thinning,
the residue is the GENE SET -- a different, harder 300 -- and no amount of coverage fixes it.

Each thinning is run at three seeds, and the mean and spread are reported.

    uv run python thin233.py
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

import panels233 as P  # noqa: E402
import score233 as S  # noqa: E402
import fold233  # noqa: E402
import plolo  # noqa: E402
from plolo import E062, Operator  # noqa: E402

SEEDS = (2331, 2332, 2333)
PAIR = {"orion_hct116": "orion_hct116", "orion_hek293t": "orion_hek293t",
        "k562_gwps": "k562_gwps"}


def thin_residual(res, keep_n: int, rng):
    """Keep `keep_n` of the carrying rows at random; zero the rest (a coverage drop)."""
    counts = np.asarray(res.target_condition_count)
    carrying = np.flatnonzero(counts > 0)
    if keep_n >= len(carrying):
        return res, len(carrying)
    keep = set(rng.choice(carrying, size=keep_n, replace=False).tolist())
    mask = np.array([i in keep for i in range(len(counts))], dtype=bool)
    eff = np.asarray(res.mean_effect).copy()
    votes = np.asarray(res.gene_votes).copy()
    cnt = counts.copy()
    eff[~mask] = 0.0
    votes[~mask] = 0
    cnt[~mask] = 0
    return E062.CD4Prepared(mode=res.mode + f"|thin{keep_n}", mean_effect=eff,
                            gene_votes=votes, target_condition_count=cnt,
                            condition_details=[]), keep_n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", nargs="*", default=["orion_hct116", "orion_hek293t"])
    ap.add_argument("--out", default=str(P.WORK / "thin233.json"))
    args = ap.parse_args()

    score = json.loads((P.WORK / "score233.json").read_text())
    R = {(r["panel"], r["fold"], r["arm"]): r for r in score}
    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"], r["seed"]) for r in rows}
    resid218, compose218 = S._resid_module()

    ctrl_mod = __import__("ctrl228")
    ctrl = {c: ctrl_mod.control_block(c) for c in S.CONTEXTS}
    bctx = {c: plolo.challenge_bctrl(c) for c in S.CONTEXTS}
    gate = {}
    for c in S.CONTEXTS:
        e = np.expm1(bctx[c])
        gate[c] = (1e6 * e / float(e.sum()) > plolo.MIN_CPM) & (bctx[c] > 0)

    from mapfit import to_lfc
    from phase23 import POOL_N, topn

    for held in args.folds:
        # ---- the pessimistic panel's coverage, measured, as the thinning target
        pess = fold233.build("pess", held)
        p_src = pess["audit"]["source_coverage"]
        p_n = pess["audit"]["n_targets"]
        pess_row = R.get(("pess", held, "BASE"), {})
        p_cd4 = pess_row.get("sizing", {}).get("cd4_carrying_rows")
        p_h1 = pess_row.get("sizing", {}).get("h1_carrying_rows")
        thin_source = min(p_src, key=lambda k: p_src[k] / p_n)      # the one that collapses
        frac_target = p_src[thin_source] / p_n
        cd4_frac = (p_cd4 / p_n) if p_cd4 is not None else None
        h1_frac = (p_h1 / p_n) if p_h1 is not None else None
        print(f"\n[{held}] pessimistic target rates: {thin_source} {frac_target:.3f}, "
              f"CD4 {cd4_frac}, H1 {h1_frac}")

        S.bind_panel("abc", resid218)
        fold = fold233.build("abc", held)
        sc = np.asarray(fold["scored"])
        wcols = np.asarray(fold["wcols"])
        Ts = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
        targets = fold["targets"]
        n = len(targets)
        cd4 = S.cd4_prepared("abc", targets)
        h1 = compose218.union_rows(
            [resid218.source_residual(x, targets, "gene_norm") for x in S.UNION])
        op = Operator(fold["op_path"])
        cpm_sc = 1e6 * np.expm1(fold["b_ctrl"][sc]) / float(np.expm1(fold["b_ctrl"]).sum())
        LR = to_lfc(np.asarray(fold["truth"], dtype=np.float64),
                    fold["b_ctrl"], gate["A"])[:, sc]
        pools = {"wald89": topn(np.abs(LR) * np.sqrt(np.maximum(cpm_sc, 0)), POOL_N),
                 "lfc89": topn(np.abs(LR), POOL_N), "h89": topn(np.abs(Ts), POOL_N)}
        real_de = {k: topn(np.abs(Ts), k) for k in S.N_REALS}

        for arm in ("thin_k562", "thin_cd4", "thin_h1", "thin_all"):
            for seed in SEEDS:
                if (held, arm, seed) in done:
                    continue
                rng = np.random.default_rng(seed)
                f = dict(fold)
                if arm in ("thin_k562", "thin_all"):
                    avail = np.asarray(fold[f"avail__{thin_source}"], bool).copy()
                    have = np.flatnonzero(avail)
                    want = int(round(frac_target * n))
                    if want < len(have):
                        drop = rng.choice(have, size=len(have) - want, replace=False)
                        avail[drop] = False
                    f = {**fold, f"avail__{thin_source}": avail}
                cd4_u, cd4_kept = (thin_residual(cd4, int(round(cd4_frac * n)), rng)
                                   if arm in ("thin_cd4", "thin_all") and cd4_frac
                                   else (cd4, int((np.asarray(cd4.target_condition_count) > 0).sum())))
                h1_u, h1_kept = (thin_residual(h1, int(round(h1_frac * n)), rng)
                                 if arm in ("thin_h1", "thin_all") and h1_frac
                                 else (h1, int((np.asarray(h1.target_condition_count) > 0).sum())))
                pre, _ = fold233.consensus(f, pm=True, only=None)
                aud: dict = {}
                final = S.build_final(fold, op, pre, cd4_u, h1_u, wcols, S.CD4_F, S.UNION_G,
                                      (resid218, compose218), audit=aud)
                t0 = time.time()
                m = S.score_arm(final, fold, bctx, gate, ctrl, Ts, LR, pools, real_de,
                                np.random.default_rng(7), [150, 250], full=True)
                rows.append(dict(fold=held, arm=arm, seed=seed, panel="abc_thinned",
                                 n_targets=int(final.shape[0]),
                                 thinned_source=thin_source,
                                 source_kept=int(np.asarray(f[f"avail__{thin_source}"]).sum()),
                                 cd4_kept=cd4_kept, h1_kept=h1_kept,
                                 sizing=aud, seconds=round(time.time() - t0, 1), **m))
                out_path.write_text(json.dumps(rows, indent=2))
                print(f"[{held}/{arm}/s{seed}] cos {m['cos_agg']:.5f} pds {m['pds']:.5f} "
                      f"nmae {m['nmae']:.5f} fid {m['fid']:.5f} reach {m['reach']:.5f} "
                      f"jac89 {m['jac']['89']:.5f}", flush=True)
                del final
        del op
    out_path.write_text(json.dumps(rows, indent=2))
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
