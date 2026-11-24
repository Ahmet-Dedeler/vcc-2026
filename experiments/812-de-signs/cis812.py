"""812 H-CISFRONT — does reach gain if the target's expressed cis neighbours (TSS within 5 kb, GENCODE v44; dCas9-KRAB
spreads to bidirectional partners) are put at the FRONT of our call order with a DOWN sign? 425 closed cis for nmae and
showed the bulk rows carry most neighbours; this asks only where they sit in the DE ORDER, which reach reads.
Analytic on the realized stats: the cis genes get a p below every other call and lfc = min(ours, -1); nothing else moves.
KA (abort on failure): base reach and fid recomputed from the stats equal the arm's recorded raw values (1e-12).

usage: uv run python experiments/812-de-signs/cis812.py   -> work/cis812.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import surf812 as SF  # noqa: E402

MET = {"h1test": SF.ROOT / "experiments/435-zs-on-coarse/work/h1/metrics-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.json",
       "h1tv": SF.ROOT / "experiments/618-h1-trainval-surface/work/g1/metrics-base_s{s}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.json"}
INV = dict(fid=3.32482, reach=1.12927)
tss = pd.read_csv(SF.ROOT / "experiments/425-crispri-cis-neighbours/work/tss_v44.tsv", sep="\t")
tss = tss.drop_duplicates("gene_name").set_index("gene_name")
out = {}
for surf in ("h1test", "h1tv"):
    for s in (0, 1, 2):
        S = SF.load(surf, s)
        adj, sig, call = SF.masks(S)
        rl, ol, oq, op, genes, targets = S["rl"], S["ol"], S["oq"], S["op"], S["genes"], S["targets"]
        T, G = rl.shape
        rec = json.load(open(str(MET[surf]).format(s=s)))["raw"]
        rr0 = np.array([SF.reach_row(SF.realized_order_pos(oq[i], op[i], ol[i], genes), sig[i], adj[i], ol[i], rl[i]) for i in range(T)])
        f0 = float(np.nanmean(SF.fid_rows(call, ol, sig, adj, rl)))
        assert abs(np.nanmean(rr0) - rec["reach"]) < 1e-12 and abs(f0 - rec["fid"]) < 1e-12, f"KA FAIL {surf} s{s}"
        gch = tss.reindex(genes)
        res = {"reach0": float(np.nanmean(rr0)), "fid0": f0}
        for dmax in (1000, 5000):
            oq1, op1, ol1, call1 = oq.copy(), op.copy(), ol.copy(), call.copy()
            ncis, nsig, nright_before, front_before = 0, 0, 0, []
            for i, t in enumerate(targets):
                if t not in tss.index:
                    continue
                c, x = tss.loc[t, "chrom"], tss.loc[t, "tss"]
                d = np.where(gch["chrom"].to_numpy() == c, np.abs(gch["tss"].to_numpy() - x), np.inf)
                cis = (d <= dmax) & ~S["own"][i]
                if not cis.any():
                    continue
                pos = SF.realized_order_pos(oq[i], op[i], ol[i], genes)
                for j in np.flatnonzero(cis):
                    ncis += 1
                    if sig[i, j]:
                        nsig += 1
                        nright_before += int(np.sign(ol[i, j]) == np.sign(rl[i, j]))
                        front_before.append(int(pos[j]))
                # front: below every other p, by distance
                k = np.flatnonzero(cis)[np.argsort(d[cis])]
                oq1[i, k] = 1e-300 * (1 + np.arange(k.size))
                op1[i, k] = 1e-300 * (1 + np.arange(k.size))
                ol1[i, k] = np.minimum(ol[i, k], -1.0)
                call1[i, k] = True
            rr1 = np.array([SF.reach_row(SF.realized_order_pos(oq1[i], op1[i], ol1[i], genes), sig[i], adj[i], ol1[i], rl[i]) for i in range(T)])
            f1 = float(np.nanmean(SF.fid_rows(call1, ol1, sig, adj, rl)))
            dr, df = float(np.nanmean(rr1) - np.nanmean(rr0)), f1 - f0
            res[f"cis{dmax}"] = dict(pairs=ncis, pairs_sig=nsig, our_sign_right_on_sig_before=nright_before,
                                     median_rank_before=float(np.median(front_before)) if front_before else None,
                                     dreach=round(dr, 6), dfid=round(df, 6), d_overall=round((dr * INV["reach"] + df * INV["fid"]) / 6, 6),
                                     targets_reach_changed=int((np.abs(np.nan_to_num(rr1) - np.nan_to_num(rr0)) > 1e-12).sum()))
        out[f"{surf}/s{s}"] = res
        print(surf, s, json.dumps(res), flush=True)
json.dump(out, open(SF.HERE / "work/cis812.json", "w"), indent=1)
