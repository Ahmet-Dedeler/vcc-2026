"""233 — what the Oct-22 operator trap COSTS, measured on the six members.

`trap233.py` shows every production-side check passing with the A/B/C operator on a new
panel.  This scores it: the pessimistic panel, all three cancer folds, the shipped stack,
once with the operator RE-FITTED for that panel and once with the A/B/C operator applied to
the new panel's axis -- the exact silent failure DEF-RUNBOOK §0.1 describes.  212 priced it at
-0.0429 overall on `pds` alone from a different route.
"""
import json, os, sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import panels233 as P, score233 as S, fold233, plolo
from plolo import Operator

resid218, compose218 = S._resid_module()
ctrl_mod = __import__("ctrl228")
ctrl = {c: ctrl_mod.control_block(c) for c in S.CONTEXTS}
bctx = {c: plolo.challenge_bctrl(c) for c in S.CONTEXTS}
gate = {}
for c in S.CONTEXTS:
    e = np.expm1(bctx[c]); gate[c] = (1e6*e/float(e.sum()) > plolo.MIN_CPM) & (bctx[c] > 0)
from mapfit import to_lfc
from phase23 import POOL_N, topn

out, rows = {}, []
for panel in ("pess", "real"):
    S.bind_panel(panel, resid218)
    for held in P.CANCER_FOLDS:
        fold = fold233.build(panel, held)
        sc = np.asarray(fold["scored"]); wcols = np.asarray(fold["wcols"])
        Ts = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
        t = fold["targets"]
        cd4 = S.cd4_prepared(panel, t)
        h1 = compose218.union_rows([resid218.source_residual(n, t, "gene_norm") for n in S.UNION])
        pre, _ = fold233.consensus(fold, pm=True, only=None)
        cpm_sc = 1e6*np.expm1(fold["b_ctrl"][sc])/float(np.expm1(fold["b_ctrl"]).sum())
        LR = to_lfc(np.asarray(fold["truth"], dtype=np.float64), fold["b_ctrl"], gate["A"])[:, sc]
        pools = {"wald89": topn(np.abs(LR)*np.sqrt(np.maximum(cpm_sc, 0)), POOL_N),
                 "lfc89": topn(np.abs(LR), POOL_N), "h89": topn(np.abs(Ts), POOL_N)}
        real_de = {n: topn(np.abs(Ts), n) for n in S.N_REALS}
        for tag, path in (("correct", P.op_path(panel, held)),
                          ("abc_operator", P.op_path("abc", held))):
            op = Operator(path)             # loaded directly: the assertion is bypassed ON PURPOSE
            final = S.build_final(fold, op, pre, cd4, h1, wcols, S.CD4_F, S.UNION_G,
                                  (resid218, compose218))
            m = S.score_arm(final, fold, bctx, gate, ctrl, Ts, LR, pools, real_de,
                            np.random.default_rng(7), None, full=True)
            rows.append(dict(panel=panel, fold=held, operator=tag, **m))
            print(f"[{panel}/{held}/{tag}] cos {m['cos_agg']:.5f} pds {m['pds']:.5f} "
                  f"nmae {m['nmae']:.5f} fid {m['fid']:.5f} reach {m['reach']:.5f} "
                  f"jac89 {m['jac']['89']:.5f}", flush=True)
            del op, final
R = {(r["panel"], r["fold"], r["operator"]): r for r in rows}
print(f"\n{'panel':6}{'fold':15}{'d cos %':>10}{'d pds':>10}{'d nmae':>10}{'d fid':>9}"
      f"{'d reach':>9}{'d jac':>9}")
for panel in ("pess", "real"):
    for held in P.CANCER_FOLDS:
        a, b = R[(panel, held, "correct")], R[(panel, held, "abc_operator")]
        d = dict(rel_cos=100*(b["cos_agg"]/a["cos_agg"]-1), pds=b["pds"]-a["pds"],
                 nmae=b["nmae"]-a["nmae"], fid=b["fid"]-a["fid"],
                 reach=b["reach"]-a["reach"], jac=b["jac"]["89"]-a["jac"]["89"])
        out[f"{panel}|{held}"] = d
        print(f"{panel:6}{held:15}{d['rel_cos']:9.2f}%{d['pds']:10.5f}{d['nmae']:10.5f}"
              f"{d['fid']:9.5f}{d['reach']:9.5f}{d['jac']:9.5f}")
(P.WORK/"trapprice233.json").write_text(json.dumps(dict(rows=rows, delta=out), indent=2))
print("-> work/trapprice233.json")
