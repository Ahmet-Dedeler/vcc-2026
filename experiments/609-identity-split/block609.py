#!/usr/bin/env python3
"""609 Part A ([4d0207] for [7772dc]): block622.py (622k/622m, 7772dc) VERBATIM except ONE switch, ZC_MODE (argv[2]):
`centred` reproduces block622 (zc = z(P - 1/2), mean-centred; KA: equals 622's block622_<mode>.json); `uncentred` uses
639d's and production v3's rule, zc = (P - 1/2) / sd (the sign of zc is CLEAN4's own sign). Is 622m's block deficit a
consequence of the centring?

Original docstring follows.
622k: the identity rider's gain restricted to the CPM BLOCK on 622g's Broad-16 held-out lines.

This is the hostile-line tail for 639's block prereg ([21dd4c]'s review item 4, asked by [d4cd5f] 2026-10-02 ~13:10Z). It also
checks whether the identity gain is an EXPRESSION-LEVEL effect ("highly expressed genes go down") rather than identity
(the partial-out rows, 7772dc).

Per held-out working line L (12 lines; broad622.load_matrix never parses the 4 reserved ones), per target t:
  head(L, t)  = the 400 genes with the largest |consensus lfc| for t over the OTHER 11 lines, t excluded. Two key analogues
                (argv[1]): `gated` (default) = mean over the other lines of lfc * [fdr < 0.05], i.e. a consensus of their DE
                calls, closest to production's donor-stack |key|; `raw` = NaN-safe mean of their raw lfc, which favours
                noisy low-count genes with extreme log ratios (its head barely meets the gate: 16-71 gate pairs per line).
                If t is screened in < 2 other lines there is no head; block = top 2,200 by expression, and the count is
                reported.
  block(L, t) = the top 2,200 genes by L's control CPM after removing head(L, t) and t. This is an APPROXIMATE production
                mask. Production's head is 400 by |key| from the donor stack and its block is 2,200 by control expression;
                639d's ds612c construction does the same on the folds.
  gate        = the authors' fdr < 0.05, finite lfc != 0, gene != t (622g).
Sign sources (positive means predict down), 622g's construction:
  zc = z(P_C4 - 0.5): CLEAN4, all four CRISPRi banks.
  zi = z(I_L), with I_L = log1p CPM_L - mean of log1p CPM over the other 11 lines (same-lab, Flex-analogue refs),
       on genes > 5 CPM in L.
  rider = zc + 0.5 zi, only where both are finite. There is NO identity-only extension: genes without CLEAN4 drop out
       of the pairing.
Reads, per line: d = acc(rider) - acc(CLEAN4) on gate genes that both sources cover, separately on block, head and all
gate genes; CI from a 2,000-draw target bootstrap.
Partial-out, on all gate genes and on block gate genes:
  E_L = log1p CPM_L (absolute expression);
  rI  = I with E regressed out (OLS on L's expressed genes); rE = E with I regressed out;
  riders zc + 0.5 z(E) and zc + 0.5 z(rI), against the same CLEAN4 pairing.
Analysis only, < 4 GB.
usage: uv run python experiments/622-identity-track/block622.py [gated|raw]   -> work/block622_<mode>.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "622-identity-track"))
import broad622 as B  # noqa: E402

RNG = np.random.default_rng(62211)
HEAD, BLOCK, LAM = 400, 2200, 0.5
MODE = sys.argv[1] if len(sys.argv) > 1 else "gated"
assert MODE in ("gated", "raw"), MODE
ZC_MODE = sys.argv[2] if len(sys.argv) > 2 else "uncentred"
assert ZC_MODE in ("centred", "uncentred"), ZC_MODE


def z(v):
    return (v - np.nanmean(v)) / np.nanstd(v)


def resid(y, x):
    """OLS residual of y on x over the genes where both are finite (NaN elsewhere)."""
    ok = np.isfinite(y) & np.isfinite(x)
    b = np.polyfit(x[ok], y[ok], 1)
    return np.where(ok, y - np.polyval(b, x), np.nan)


def main() -> None:
    gl, LFC, cl = B.load_matrix(B.P / "lfc_matrix.csv", 1)
    gf, FDR, cf = B.load_matrix(B.P / "fdr_matrix.csv", 1)
    assert list(gl) == list(gf) and cl == cf
    gp, PB, cp = B.load_matrix(B.P / "pseudobulk_sum_by_ko.csv", 0)
    lines = sorted({ln for ln, _ in cl})
    assert not (set(lines) & B.RESERVED) and len(lines) == 12, lines
    all_t = set(k for _, k in cl)
    ctrl_cols = {ln: [j for j, (l2, k) in enumerate(cp) if l2 == ln and k not in all_t] for ln in lines}
    gpos = {g: i for i, g in enumerate(gp)}
    rowmap = np.array([gpos.get(g, -1) for g in gl])
    cpm = {}
    for ln in lines:                                   # 622g's control CPM, absent genes NaN, then onto the lfc axis
        s = np.nansum(PB[:, ctrl_cols[ln]], axis=1)
        absent = np.all(~np.isfinite(PB[:, ctrl_cols[ln]]), axis=1)
        v = 1e6 * s / np.nansum(s); v[absent] = np.nan
        cpm[ln] = np.where(rowmap >= 0, v[np.maximum(rowmap, 0)], np.nan)
    c4 = np.load(B.ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(c4["bank_names"].astype(str))
    Pc = c4["P_banks"].astype(float)[[names.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    c4map = dict(zip(c4["genes"].astype(str), Pc))
    _p = np.array([c4map.get(g, np.nan) for g in gl]) - 0.5
    zc = z(_p) if ZC_MODE == "centred" else _p / np.nanstd(_p)
    colline = np.array([ln for ln, _ in cl]); coltgt = np.array([k for _, k in cl])
    out = {"config": dict(head=HEAD, block=BLOCK, lam=LAM, head_mode=MODE, mask="approximate production mask (consensus-|lfc| head from the other "
                          "11 lines, block by the held-out line's control CPM)")}
    for ln in lines:
        others = [o for o in lines if o != ln]
        expressed = cpm[ln] > 5
        I = np.where(expressed, np.log1p(cpm[ln]) - np.nanmean([np.log1p(cpm[o]) for o in others], 0), np.nan)
        E = np.where(expressed, np.log1p(cpm[ln]), np.nan)
        zi, zE = z(I), z(E)
        zrI, zrE = z(resid(I, E)), z(resid(E, I))
        okIE = np.isfinite(I) & np.isfinite(E)
        rIE = float(np.corrcoef(I[okIE], E[okIE])[0, 1])
        srcs = {"rider": zc + LAM * zi, "identity": zi, "expr": zE, "rider_expr": zc + LAM * zE, "rider_I_perp_E": zc + LAM * zrI,
                "I_perp_E": zrI, "E_perp_I": zrE}
        jl = np.flatnonzero(colline == ln)
        rank_expr = np.argsort(-np.nan_to_num(cpm[ln], nan=-1.0), kind="stable")
        regions = {r: np.zeros((len(gl), len(jl)), bool) for r in ("block", "head")}
        no_head = 0
        for c, j in enumerate(jl):
            t = coltgt[j]
            oj = [jj for jj in np.flatnonzero((coltgt == t) & (colline != ln))]
            notg = gl != t
            head = np.zeros(len(gl), bool)
            if len(oj) >= 2:
                if MODE == "gated":
                    cons = np.abs(np.mean(np.where((FDR[:, oj] < 0.05) & np.isfinite(LFC[:, oj]), LFC[:, oj], 0.0), axis=1))
                else:
                    cons = np.abs(np.nanmean(LFC[:, oj], axis=1))
                cons = np.where(notg & np.isfinite(cons), cons, -np.inf)
                head[np.argsort(-cons, kind="stable")[:HEAD]] = True
            else:
                no_head += 1
            cand = [g for g in rank_expr if not head[g] and notg[g] and np.isfinite(cpm[ln][g])][:BLOCK]
            regions["block"][cand, c] = True
            regions["head"][:, c] = head
        gate = (FDR[:, jl] < 0.05) & np.isfinite(LFC[:, jl]) & (LFC[:, jl] != 0) & (gl[:, None] != coltgt[jl][None, :])
        truth_dn = LFC[:, jl] < 0
        res = {"targets": int(len(jl)), "targets_without_head": int(no_head), "corr_I_E": rIE}
        regions["uncalled"] = ~regions["block"] & ~regions["head"]
        for reg, R in (("all", np.ones_like(gate)), ("block", regions["block"]), ("head", regions["head"]), ("uncalled", regions["uncalled"])):
            rr = {}
            for name, sc in srcs.items():
                both = np.isfinite(sc) & np.isfinite(zc)
                m = gate & R & both[:, None]
                okS = (((sc > 0)[:, None] == truth_dn) & m).sum(0).astype(float)
                okC = (((zc > 0)[:, None] == truth_dn) & m).sum(0).astype(float)
                n = m.sum(0).astype(float)
                idx = np.flatnonzero(n > 0)
                bs = [(okS[b].sum() - okC[b].sum()) / n[b].sum() for b in (RNG.choice(idx, idx.size) for _ in range(2000))] if idx.size else [np.nan]
                rr[name] = dict(acc=float(okS.sum() / max(n.sum(), 1)), clean4=float(okC.sum() / max(n.sum(), 1)), dnet=float(okS.sum() - okC.sum()),
                                d=float((okS.sum() - okC.sum()) / max(n.sum(), 1)),
                                ci=[float(np.nanquantile(bs, .025)), float(np.nanquantile(bs, .975))], n=int(n.sum()))
            rr["gate_genes_in_region"] = int((gate & R).sum())
            res[reg] = rr
        for reg in ("block", "head", "uncalled"):    # each region's share of the total change in (right - wrong) on all gate genes
            tot = res["all"]["rider"]["dnet"]
            res[reg]["rider"]["share"] = float(res[reg]["rider"]["dnet"] / tot) if tot else float("nan")
        out[ln] = res
        b, a = res["block"], res["all"]
        print(f"{ln:9s} no-head {no_head:2d}/{len(jl)} corr(I,E) {rIE:+.2f} | BLOCK n {b['rider']['n']:5d} C4 {b['rider']['clean4']:.3f} "
              f"rider d {b['rider']['d']:+.4f} [{b['rider']['ci'][0]:+.4f}, {b['rider']['ci'][1]:+.4f}] | rider_expr {b['rider_expr']['d']:+.4f} "
              f"rider_I_perp_E {b['rider_I_perp_E']['d']:+.4f} | ALL rider {a['rider']['d']:+.4f} expr {a['rider_expr']['d']:+.4f} "
              f"I_perp_E {a['rider_I_perp_E']['d']:+.4f} | alone acc: I {a['identity']['acc']:.3f} E {a['expr']['acc']:.3f} "
              f"I_perp_E {a['I_perp_E']['acc']:.3f} E_perp_I {a['E_perp_I']['acc']:.3f}", flush=True)
    summ = {}
    for reg in ("block", "head", "uncalled", "all"):
        for name in ("rider", "rider_expr", "rider_I_perp_E", "identity", "expr", "I_perp_E", "E_perp_I"):
            ds = np.array([out[ln][reg][name]["d"] for ln in lines])
            summ[f"{reg}/{name}"] = dict(mean_d=float(ds.mean()), n_negative=int((ds < 0).sum()), min=float(ds.min()),
                                         argmin=lines[int(ds.argmin())], ci_below0=int(sum(out[ln][reg][name]["ci"][1] < 0 for ln in lines)),
                                         ci_above0=int(sum(out[ln][reg][name]["ci"][0] > 0 for ln in lines)))
    out["_summary"] = summ
    print()
    for ln in lines:
        print(f"{ln:9s} n(all/block/head/uncalled) " + "/".join(str(out[ln][r]["rider"]["n"]) for r in ("all", "block", "head", "uncalled"))
              + "  d " + "/".join(f"{out[ln][r]['rider']['d']:+.4f}" for r in ("all", "block", "head", "uncalled"))
              + "  share of net " + "/".join(f"{out[ln][r]['rider']['share']:+.2f}" for r in ("block", "head", "uncalled")))
    pooled = {r: sum(out[ln][r]["rider"]["dnet"] for ln in lines) for r in ("all", "block", "head", "uncalled")}
    out["_summary"]["pooled_net_share"] = {r: float(pooled[r] / pooled["all"]) for r in ("block", "head", "uncalled")}
    print("pooled share of the net (right - wrong) gain:", {r: round(v, 2) for r, v in out["_summary"]["pooled_net_share"].items()})
    for k in ("block/rider", "head/rider", "uncalled/rider", "all/rider", "block/rider_expr", "block/rider_I_perp_E", "all/rider_expr", "all/rider_I_perp_E"):
        s = summ[k]
        print(f"{k:22s} mean d {s['mean_d']:+.4f}; negative lines {s['n_negative']}/12; min {s['min']:+.4f} ({s['argmin']}); "
              f"CI<0 {s['ci_below0']}, CI>0 {s['ci_above0']}")
    (HERE / f"work/block609_{MODE}_{ZC_MODE}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
