"""322 Phase 3 — `kolf2` as an EXTRA BASE DONOR (264's construction), on all three panels.

264 measured this on today's panel with 162's production operator.  The mandate asks whether
the arm that loses today degrades BETTER tomorrow, so the same construction has to run on
233's two disjoint draws.  233's folds are panel-parameterised and carry the per-donor rows,
so the construction transplants exactly:

    (c)  the incumbent two-donor consensus                 = 233's BASE
    (a)  `E062.build_baseline([d1_pm50, d2, kolf2])`       = 264's `mean3`

264's `lib264.kolf2_as_donor` (lines 144-192) is reproduced operation for operation:
`E062.prepare_main_source` with the fold's targets and panel positions -> restrict the
effects and `mean_on_axis` to the operator's columns -> multiply by the fold weight ->
`build_baseline`.  Two weights, 264's two constructions:

    `equal`   weight 1.0                                     (264 §0.1, "by the letter")
    `matched` mean_j(whitened-unpreserved median norm of production donor j)
              / (whitened-unpreserved median norm of kolf2)   (264 `stage3:209`)

`whitening_audit` (264 lines 116-141) is reproduced verbatim for the amplification /
bottom-eigen-decile columns, because 264's central limitation -- kolf2's rows are pushed
through an operator fitted without them -- applies here too and has to be visible.

WHAT IS *NOT* A KNOWN-ANSWER ROW.  233 re-fits every operator with its own recipe (no
essential prune, no 12,000-row cap), and 233's own "Bridge" paragraph measures that recipe
change at -3.5 / -5.8 / -1.8 % `cos_agg` and +0.0018 / -0.0029 / -0.0046 `pds` on today's
panel.  So the `abc` numbers here are the SAME CONSTRUCTION ON A DIFFERENT OPERATOR than
264's and are expected to agree in sign and order, not to the digit.  They are labelled that
way and 264's numbers are quoted separately.

    uv run python basedonor322.py --panels abc real pess
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib322 as L  # noqa: E402

import panels233 as P  # noqa: E402
import score233 as S  # noqa: E402
import fold233  # noqa: E402
import plolo  # noqa: E402
from plolo import E062, Operator, OFFICIAL  # noqa: E402

THIRD = "kolf2"
GAMMA, LAM = S.GAMMA, S.LAM


def whitening_audit(op, wcols, rows, avail, bA, gA):
    """264 `lib264.whitening_audit`, verbatim in its computed quantities."""
    R = np.asarray(rows)[np.asarray(avail, bool)][:, wcols]
    if R.shape[0] == 0:
        return dict(n=0, raw_median=0.0, wht_unpreserved_median=0.0, amplification=float("nan"),
                    energy_bottom_eigen_decile_wht=float("nan"))
    nr = np.linalg.norm(R, axis=1)
    w = (op.evals + LAM) ** (-GAMMA / 2.0)
    order = np.argsort(op.evals)[::-1]
    T = (R @ op.evecs) * w
    nt = np.linalg.norm(T, axis=1)
    e_T = (T ** 2).sum(0)
    bottom = order[9 * len(order) // 10:]
    return dict(n=int(np.asarray(avail).sum()), raw_median=float(np.median(nr)),
                wht_unpreserved_median=float(np.median(nt)),
                amplification=float(np.median(nt) / max(np.median(nr), 1e-300)),
                energy_bottom_eigen_decile_wht=float(e_T[bottom].sum() / max(e_T.sum(), 1e-300)))


def prepared_list(fold):
    """`fold233.consensus`'s prepared donors (PerturbMatch blend included), exposed."""
    out = []
    for n in fold["sources"]:
        eff = np.asarray(fold[f"eff__{n}"], dtype=np.float64).copy()
        name = n
        if n == "k562_gwps":
            rows = np.flatnonzero(np.asarray(fold["pm_both"], bool))
            sh = np.asarray(fold["pm_shared"], dtype=np.int64)
            if len(rows) and len(sh):
                pm_eff = np.asarray(fold["pm_effects"], dtype=np.float64)
                eff[np.ix_(rows, sh)] = ((1 - fold233.PM_BLEND) * eff[np.ix_(rows, sh)]
                                         + fold233.PM_BLEND * pm_eff[np.ix_(rows, sh)])
            name = "k562_gwps_pm50"
        out.append(fold233._prepared(name, eff, fold[f"avail__{n}"], fold[f"mean__{n}"],
                                     float(fold[f"rmn__{n}"]), 0, ()))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", nargs="*", default=["abc", "real", "pess"])
    ap.add_argument("--folds", nargs="*", default=list(P.CANCER_FOLDS))
    ap.add_argument("--out", default=str(L.WORK / "basedonor322.json"))
    ap.add_argument("--kolf2-root", default=None,
                    help="path template for the kolf2 bank, formatted with {panel}; "
                         "the default is the frozen `data/external/kolf2`, which was "
                         "range-extracted against TODAY'S panel and therefore carries "
                         "nothing on a disjoint one")
    ap.add_argument("--arm-suffix", default="",
                    help="appended to every arm name so a re-extracted run can share a file")
    args = ap.parse_args()

    ok, msg = L.memory_ok()
    print(f"[mem] {msg} -> {'OK' if ok else 'HOLD'}", flush=True)
    assert ok, "memory gate: hold"

    SUF = args.arm_suffix
    out_path = Path(args.out)
    rows_out = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["panel"], r["fold"], r["arm"]) for r in rows_out}
    resid218, compose218 = S._resid_module()

    ctrl_mod = __import__("ctrl228")
    ctrl = {c: ctrl_mod.control_block(c) for c in S.CONTEXTS}
    bctx = {c: plolo.challenge_bctrl(c) for c in S.CONTEXTS}
    gate = {}
    for c in S.CONTEXTS:
        e = np.expm1(bctx[c])
        gate[c] = (1e6 * e / float(e.sum()) > plolo.MIN_CPM) & (bctx[c] > 0)

    for panel in args.panels:
        S.bind_panel(panel, resid218)
        kolf_root = (L.HERE / args.kolf2_root.format(panel=panel)).resolve() \
            if args.kolf2_root else (L.ROOT / "data/external/kolf2")
        kolf_bank = E062.Bank.load(THIRD, kolf_root)
        print(f"[bank] kolf2 <- {kolf_root} "
              f"({len(kolf_bank.pert_pos)} perturbations)", flush=True)
        for held in args.folds:
            t0 = time.time()
            fold = fold233.build(panel, held)
            sc = np.asarray(fold["scored"])
            wcols = np.asarray(fold["wcols"])
            Ts = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
            targets = fold["targets"]
            cd4 = S.cd4_prepared(panel, targets)
            h1 = compose218.union_rows(
                [resid218.source_residual(n, targets, "gene_norm") for n in S.UNION])

            # `prepare_main_source` enforces a PEAK-RSS cap, and the operator's eigenvector
            # matrix is 17,942^2 x 8 B = 2.6 GiB resident, so kolf2 is prepared BEFORE the
            # operator is loaded.  One (panel, fold) per process on top of that (264's rule).
            panel_positions = np.asarray(fold["panel_positions"], dtype=np.int64)
            kp = E062.prepare_main_source(kolf_bank, targets, OFFICIAL, panel_positions,
                                          chunk_rows=256, max_rss_gib=8.0)
            op = Operator(fold["op_path"])

            from mapfit import to_lfc
            from phase23 import POOL_N, topn
            cpm_sc = 1e6 * np.expm1(fold["b_ctrl"][sc]) / float(np.expm1(fold["b_ctrl"]).sum())
            LR = to_lfc(np.asarray(fold["truth"], dtype=np.float64), fold["b_ctrl"],
                        gate["A"])[:, sc]
            pools = {"wald89": topn(np.abs(LR) * np.sqrt(np.maximum(cpm_sc, 0)), POOL_N),
                     "lfc89": topn(np.abs(LR), POOL_N), "h89": topn(np.abs(Ts), POOL_N)}
            real_de = {n: topn(np.abs(Ts), n) for n in S.N_REALS}

            # ---- the incumbent two-donor base (c) --------------------------------------
            prod = prepared_list(fold)
            pre_c, _ = E062.build_baseline(prod)

            # ---- kolf2 prepared as a donor for THIS panel's targets (264 lines 144-192) --
            eff_k = np.asarray(kp.effects, dtype=np.float64)
            avail_k = np.asarray(kp.target_available, bool)
            inside = np.zeros(eff_k.shape[1], dtype=bool)
            inside[wcols] = True
            eff_kr = eff_k.copy()
            eff_kr[:, ~inside] = 0.0
            moa_k = np.asarray(kp.mean_on_axis, dtype=np.float64).copy()
            moa_k[~inside] = 0.0

            aud_prod = {p.name: whitening_audit(op, wcols, p.effects, p.target_available,
                                                bctx["A"], gate["A"]) for p in prod}
            aud_k = whitening_audit(op, wcols, eff_kr, avail_k, bctx["A"], gate["A"])
            if aud_k["n"] and aud_k["wht_unpreserved_median"] > 0:
                w_matched = float(np.mean([aud_prod[p.name]["wht_unpreserved_median"]
                                           for p in prod]) / aud_k["wht_unpreserved_median"])
            else:
                w_matched = float("nan")

            prod_amp = ", ".join(f"{aud_prod[p.name]['amplification']:.1f}x" for p in prod)
            print(f"[{panel}/{held}] kolf2 carries {int(avail_k.sum())}/{len(targets)} fold "
                  f"targets; amplification {aud_k['amplification']:.1f}x "
                  f"(production {prod_amp})", flush=True)

            def score(arm, pre, extra_meta):
                if (panel, held, arm) in done:
                    return None
                whitened = op.apply(pre * S.ALPHA, wcols, S.GAMMA, S.LAM)
                frac_cd4, _ = plolo.solve_fraction(whitened, cd4, S.CD4_F)
                frac_u, _ = compose218.solve_fraction_carrying(whitened, h1, S.UNION_G)
                add = compose218.merged(cd4, frac_cd4, h1, frac_u)
                final, _ = E062.add_post_whiten_residual(whitened, add, 1.0, wcols)
                final[:, fold["panel_positions"]] = 0.0
                m = S.score_arm(final, fold, bctx, gate, ctrl, Ts, LR, pools, real_de,
                                np.random.default_rng(7), None, full=True)
                rows_out.append(dict(panel=panel, fold=held, arm=arm,
                                     n_targets=int(final.shape[0]),
                                     seconds=round(time.time() - t0, 1),
                                     **extra_meta, **m))
                out_path.write_text(json.dumps(rows_out, indent=2))
                print(f"[{panel}/{held}/{arm}] cos {m['cos_agg']:.6f} pds {m['pds']:.6f}",
                      flush=True)
                del final, whitened
                return m

            meta_common = dict(kolf2_carrying_rows=int(avail_k.sum()),
                               kolf2_amplification=aud_k["amplification"],
                               kolf2_bottom_decile_energy=aud_k["energy_bottom_eigen_decile_wht"],
                               production_amplification={k: v["amplification"]
                                                         for k, v in aud_prod.items()})
            score(f"BASE_2donor{SUF}", pre_c, dict(weight=None, rows_changed=0,
                                             max_abs_diff=0.0, **meta_common))

            for tag, w in (("equal", 1.0), ("matched", w_matched)):
                if not np.isfinite(w):
                    # kolf2 carries nothing: build_baseline over [d1, d2, kolf2] is the
                    # two-donor mean EXACTLY -- a zero-row donor contributes no term and no
                    # count.  Demonstrated, not asserted.
                    src = fold233._prepared(THIRD, eff_kr * 0.0, avail_k, moa_k * 0.0,
                                            float(kp.raw_median_norm), 0, ())
                    pre_a, _ = E062.build_baseline(prod + [src])
                    d = float(np.abs(pre_a - pre_c).max())
                    changed = int((np.abs(pre_a - pre_c) > 0).any(1).sum())
                    rows_out.append(dict(panel=panel, fold=held, arm=f"mean3:{tag}{SUF}",
                                         n_targets=len(targets), weight=None,
                                         rows_changed=changed, max_abs_diff=d,
                                         seconds=round(time.time() - t0, 1),
                                         note="kolf2 carries 0 rows on this panel — the "
                                              "three-donor base IS the two-donor base",
                                         **meta_common))
                    out_path.write_text(json.dumps(rows_out, indent=2))
                    print(f"[{panel}/{held}/mean3:{tag}{SUF}] IDENTITY — rows changed {changed}, "
                          f"max|d| {d:.3e}", flush=True)
                    continue
                src = fold233._prepared(THIRD, eff_kr * w, avail_k, moa_k * w,
                                        float(kp.raw_median_norm), 0, ())
                pre_a, _ = E062.build_baseline(prod + [src])
                d = float(np.abs(pre_a - pre_c).max())
                changed = int((np.abs(pre_a - pre_c) > 0).any(1).sum())
                score(f"mean3:{tag}{SUF}", pre_a,
                      dict(weight=w, rows_changed=changed, max_abs_diff=d, **meta_common))
            del op
    out_path.write_text(json.dumps(rows_out, indent=2))
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
