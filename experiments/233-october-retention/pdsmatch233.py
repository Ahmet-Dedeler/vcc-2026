"""233 — `pds` at a decoy count every panel can supply.

The evaluable-target count differs by panel and fold (150 to 300), and `pds` is a retrieval
statistic whose value depends on how many decoys there are.  140 is below every fold's n, so
`pds@140` is the one matched statistic all nine (panel, fold) cells can report.  Computed for
the shipped stack and for its deranged null, so the retention ratio is a ratio of credits.
"""
import json, os, sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import panels233 as P, score233 as S, fold233, plolo
from plolo import Operator, deliver

NS = (140, 250)
resid218, compose218 = S._resid_module()
bctx = {c: plolo.challenge_bctrl(c) for c in S.CONTEXTS}
gate = {}
for c in S.CONTEXTS:
    e = np.expm1(bctx[c]); gate[c] = (1e6 * e / float(e.sum()) > plolo.MIN_CPM) & (bctx[c] > 0)

out = {}
for panel in ("abc", "real", "pess"):
    S.bind_panel(panel, resid218)
    for held in P.CANCER_FOLDS:
        fold = fold233.build(panel, held)
        sc = np.asarray(fold["scored"]); wcols = np.asarray(fold["wcols"])
        Ts = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
        t = fold["targets"]
        cd4 = S.cd4_prepared(panel, t)
        h1 = compose218.union_rows([resid218.source_residual(n, t, "gene_norm") for n in S.UNION])
        op = Operator(fold["op_path"])
        pre, _ = fold233.consensus(fold, pm=True, only=None)
        final = S.build_final(fold, op, pre, cd4, h1, wcols, S.CD4_F, S.UNION_G,
                              (resid218, compose218))
        rng = np.random.default_rng(S.NULL_SEED)
        n = final.shape[0]
        perm = rng.permutation(n)
        while (perm == np.arange(n)).any():
            perm = rng.permutation(n)
        for arm, M in (("BASE", final), ("null", final[perm])):
            dl = deliver(M, bctx["A"], gate["A"], "minnorm")[:, sc]
            rec = {f"pds_at{k}": S.pds_at(dl, Ts, k, 128, np.random.default_rng(9090 + k))
                   for k in NS if k < n}
            out[f"{panel}|{held}|{arm}"] = rec | dict(n=int(n))
            print(panel, held, arm, rec, flush=True)
        del op, final
(P.WORK / "pdsmatch233.json").write_text(json.dumps(out, indent=2))
print("-> work/pdsmatch233.json")
