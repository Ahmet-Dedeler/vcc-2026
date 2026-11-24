"""812 step 1 — where a0 (sign accuracy on OUR called genes that the truth does not call) lives, by gene class, on both
DE-H1 surfaces, and what an oracle per class would buy in fid (f 1).

Known answer (aborts if it fails): fid recomputed from the stats equals the arm's recorded raw fid
(h1test s0 0.5089291277653674, 435 metrics; h1tv s0 0.5052857989716864, 618 g1 metrics) to 1e-12.

usage: uv run python experiments/812-de-signs/a0class812.py [seeds=0,1,2]
writes experiments/812-de-signs/work/a0class812.json and .txt
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import surf812 as SF  # noqa: E402

HERE = SF.HERE
KA = {("h1test", 0): 0.5089291277653674, ("h1tv", 0): 0.5052857989716864}
MET = {"h1test": SF.ROOT / "experiments/435-zs-on-coarse/work/h1/metrics-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.json",
       "h1tv": SF.ROOT / "experiments/618-h1-trainval-surface/work/g1/metrics-base_s{s}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.json"}
INV_FID = 3.32482  # board 1/span (192)
seeds = [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else "0,1,2").split(",")]
out = {}
lines = []


def acc(sg, rl, M):
    n = int(M.sum())
    return (round(float((np.sign(sg[M]) == np.sign(rl[M])).mean()), 4) if n else None), n


for surf in ("h1test", "h1tv"):
    for s in seeds:
        S = SF.load(surf, s)
        adj, sig, call = SF.masks(S)
        rl, ol = S["rl"], S["ol"]
        T, G = rl.shape
        fr = SF.fid_rows(call, ol, sig, adj, rl)
        fid = float(np.nanmean(fr))
        rec = json.load(open(str(MET[surf]).format(s=s)))["raw"]["fid"]
        assert abs(fid - rec) < 1e-12, f"KA FAIL {surf} s{s}: fid {fid} vs recorded {rec}"
        if (surf, s) in KA:
            assert abs(fid - KA[(surf, s)]) < 1e-12
        R = dict(fid=fid, ka_recorded_fid=rec)
        M0 = call & adj & ~sig
        M1 = call & adj & sig
        R["a0"], R["n0"] = acc(ol, rl, M0)
        R["a1"], R["n1"] = acc(ol, rl, M1)
        R["share_nonsig_of_calls"] = round(M0.sum() / (M0.sum() + M1.sum()), 4)
        R["K_mean"] = float((call & adj).sum(1).mean())
        R["nreal_median"] = float(np.median(sig.sum(1)))
        # call-set parts
        other = call & ~S["head"] & ~S["block"]
        R["calls_in_prescribed"] = round(float((call & (S["head"] | S["block"])).sum() / call.sum()), 4)
        R["prescribed_called"] = round(float((call & (S["head"] | S["block"])).sum() / (S["head"] | S["block"]).sum()), 4)
        for nm, part in (("head", S["head"]), ("block", S["block"]), ("other", other)):
            R[f"a0_{nm}"], R[f"n0_{nm}"] = acc(ol, rl, M0 & part)
        # CLEAN4, all-down, key on the same genes
        c4 = S["c4"][None, :].repeat(T, 0)
        cov = np.isfinite(c4)
        c4s = np.where(c4 > 0.5, -1.0, 1.0)
        R["a0_emitted_on_c4cov"], _ = acc(ol, rl, M0 & cov)
        R["a0_c4_on_c4cov"], R["n0_c4cov"] = acc(c4s, rl, M0 & cov)
        R["a0_alldown"], _ = acc(-np.ones_like(rl), rl, M0)
        R["a0_key"], _ = acc(S["key"], rl, M0 & (S["key"] != 0))
        R["emitted_down_share_M0"] = round(float((ol[M0] < 0).mean()), 4)
        R["real_down_share_M0"] = round(float((rl[M0] < 0).mean()), 4)
        R["block_delivery_c4sign"] = round(float((np.sign(ol) == c4s)[call & S["block"] & cov].mean()), 4)
        # by control CPM bins
        cpm = S["cpm"][None, :].repeat(T, 0)
        cb = {}
        for lo, hi in ((5, 10), (10, 20), (20, 50), (50, 100), (100, 300), (300, 1000), (1000, 1e9)):
            m = M0 & (cpm >= lo) & (cpm < hi)
            a, n = acc(ol, rl, m)
            ac4, _ = acc(c4s, rl, m & cov)
            ad, _ = acc(-np.ones_like(rl), rl, m)
            cb[f"{lo}-{hi:g}"] = dict(a0=a, n=n, a0_c4=ac4, a0_down=ad, real_down=round(float((rl[m] < 0).mean()), 4) if n else None)
        R["by_cpm"] = cb
        # by CLEAN4 confidence
        cc = {}
        conf = np.abs(c4 - 0.5)
        for lo, hi in ((0, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.51)):
            m = M0 & cov & (conf >= lo) & (conf < hi)
            a, n = acc(ol, rl, m)
            ac4, _ = acc(c4s, rl, m)
            cc[f"{lo}-{hi}"] = dict(a0=a, a0_c4=ac4, n=n)
        m = M0 & ~cov
        cc["uncovered"] = dict(a0=acc(ol, rl, m)[0], n=int(m.sum()), real_down=round(float((rl[m] < 0).mean()), 4))
        R["by_c4conf"] = cc
        # by the truth's own nominal p (how much true signal a non-sig gene carries)
        pb = {}
        rp = S["rp"]
        for lo, hi in ((0, 0.01), (0.01, 0.05), (0.05, 0.2), (0.2, 0.5), (0.5, 1.01)):
            m = M0 & (rp >= lo) & (rp < hi)
            a, n = acc(ol, rl, m)
            ac4, _ = acc(c4s, rl, m & cov)
            pb[f"{lo}-{hi}"] = dict(a0=a, a0_c4=ac4, n=n, share=round(n / max(M0.sum(), 1), 4))
        R["by_truth_p"] = pb
        # by target strength (n_real terciles)
        nr = sig.sum(1)
        qs = np.quantile(nr, [1 / 3, 2 / 3])
        tb = {}
        for nm, tm in (("weak", nr <= qs[0]), ("mid", (nr > qs[0]) & (nr <= qs[1])), ("strong", nr > qs[1])):
            a, n = acc(ol, rl, M0 & tm[:, None])
            ac4, _ = acc(c4s, rl, M0 & tm[:, None] & cov)
            tb[nm] = dict(a0=a, a0_c4=ac4, n=n, nreal_range=[int(nr[tm].min()), int(nr[tm].max())])
        R["by_target_strength"] = tb
        # context-programme oracle: leave-one-target-out mean real lfc sign per gene (uses the truth; an upper bound)
        rlz = np.where(np.isfinite(rl), rl, 0.0)
        tot = rlz.sum(0, keepdims=True)
        loo = (tot - rlz) / max(T - 1, 1)
        R["a0_oracle_ctxmean_loo"], _ = acc(loo, rl, M0 & (loo != 0))
        rln = np.where(sig, 0.0, rlz)  # non-sig parts only
        loon = (rln.sum(0, keepdims=True) - rln) / max(T - 1, 1)
        R["a0_oracle_ctxmean_loo_nonsig"], _ = acc(loon, rl, M0 & (loon != 0))
        # oracle buys in fid: true sign on a class of called non-sig genes
        def buy(cls):
            sg = np.where(cls, np.sign(rl), ol)
            return round(float(np.nanmean(SF.fid_rows(call, sg, sig, adj, rl))) - fid, 5)
        R["oracle_fid_buy"] = {
            "all_M0": buy(M0), "block_M0": buy(M0 & S["block"]), "head_M0": buy(M0 & S["head"]),
            "M0_truth_p<0.05": buy(M0 & (rp < 0.05)), "M0_truth_p<0.2": buy(M0 & (rp < 0.2)),
            "ctxmean_oracle_on_M0": round(float(np.nanmean(SF.fid_rows(call, np.where(M0 & (loo != 0), loo, ol), sig, adj, rl))) - fid, 5),
        }
        R["oracle_fid_buy_overall_eq"] = {k: round(v * INV_FID / 6, 5) for k, v in R["oracle_fid_buy"].items()}
        out[f"{surf}/s{s}"] = R
        print(f"{surf} s{s}: fid {fid:.6f} (KA ok) a0 {R['a0']} n0 {R['n0']} a1 {R['a1']} n1 {R['n1']} K {R['K_mean']:.0f}", flush=True)

out["inputs"] = {k: {kk: (str(vv), SF.sha(vv)) for kk, vv in v.items() if kk != "st"} for k, v in SF.SURF.items()}
out["inputs"]["clean4"] = (str(SF.C4), SF.sha(SF.C4))
json.dump(out, open(HERE / "work/a0class812.json", "w"), indent=1)
print("wrote work/a0class812.json")
