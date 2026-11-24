#!/usr/bin/env python3
"""622h — H-633's decisive test (d4cd5f): does a SAME-LAB programme learned from only 3 lines transfer to new lines?

A/B/C -> D/E/F analogue on the Broad-16 (12 working lines; reserved 4 never read). For each held-out working line L and
40 random triples of OTHER lines: programme P_down(g) = fraction down among the triple's gated responses (>= 3 gated);
score its gate-sign accuracy on L vs CLEAN4 on shared genes. Also a COARSE version, what a K-module board readout could
recover: genes binned into K = 16 modules by (CLEAN4 P_down octile x control-expression half, defined without truth);
each module gets the triple's majority sign; then (a) module sign alone, (b) CLEAN4 sign flipped to the module sign only
in modules whose triple majority disagrees with CLEAN4's majority. Analysis only, < 4 GB.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import broad622 as B  # noqa: E402

RNG = np.random.default_rng(6228)


def main() -> None:
    gl, LFC, cl = B.load_matrix(B.P / "lfc_matrix.csv", 1)
    gf, FDR, cf = B.load_matrix(B.P / "fdr_matrix.csv", 1)
    gp, PB, cp = B.load_matrix(B.P / "pseudobulk_sum_by_ko.csv", 0)
    lines = sorted({ln for ln, _ in cl})
    assert not (set(lines) & B.RESERVED) and len(lines) == 12
    c4 = np.load(B.ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(c4["bank_names"].astype(str))
    Pc = c4["P_banks"].astype(float)[[names.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    c4map = dict(zip(c4["genes"].astype(str), Pc))
    C4 = np.array([c4map.get(g, np.nan) for g in gl])
    gated = (FDR < 0.05) & np.isfinite(LFC) & (LFC != 0)
    dn = (LFC < 0) & gated
    colline = np.array([ln for ln, _ in cl])
    per_line_dn = {ln: dn[:, colline == ln].sum(1) for ln in lines}
    per_line_g = {ln: gated[:, colline == ln].sum(1) for ln in lines}
    # pooled control expression across all 12 lines (truth-free) for module bins
    all_t = set(k for _, k in cl)
    ctrl = np.nansum(PB[:, [j for j, (l2, k) in enumerate(cp) if k not in all_t]], axis=1)
    gpos = {g: i for i, g in enumerate(gp)}
    ex = np.array([ctrl[gpos[g]] if g in gpos else np.nan for g in gl])
    oct_ = np.digitize(C4, np.nanquantile(C4, np.linspace(0.125, 0.875, 7)))
    half = (ex > np.nanmedian(ex)).astype(int)
    module = np.where(np.isfinite(C4) & np.isfinite(ex), oct_ * 2 + half, -1)
    BATCH = {"HCT15": 31, "KYSE140": 31, "KYSE450": 31, "RPMI7951": 31, "SKGII": 32, "U343": 32, "C4I": 33,
             "ONE58": 34, "T3M5": 34, "KMRC20": 35, "SNU8": 35, "HKGZCC": 38}
    CROSS = "--cross-batch" in sys.argv
    out = {}
    for L in lines:
        others = [o for o in lines if o != L and (not CROSS or BATCH[o] != BATCH[L])]
        jl = np.flatnonzero(colline == L)
        tl = [cl[j][1] for j in jl]
        gate = gated[:, jl] & (gl[:, None] != np.array(tl)[None, :])
        tdn = LFC[:, jl] < 0
        okc4 = np.isfinite(C4)
        rows = {"full11": [], "tri": [], "tri_mod": [], "tri_modflip": []}

        def acc(sign_dn, mask):
            m = gate & mask[:, None]
            return float((((sign_dn)[:, None] == tdn) & m).sum() / max(m.sum(), 1)), int(m.sum())
        c4dn = C4 > 0.5
        # full 11-line programme
        sd = sum(per_line_dn[o] for o in others); sg = sum(per_line_g[o] for o in others)
        pr = (sd + 0.5) / (sg + 1.0); ok = okc4 & (sg >= 3)
        a, n = acc(pr > 0.5, ok); a0, _ = acc(c4dn, ok)
        res = {"full11_d": a - a0}
        dtri, dmod, dflip = [], [], []
        for _ in range(40):
            tri = list(RNG.choice(others, 3, replace=False))
            sd = sum(per_line_dn[o] for o in tri); sg = sum(per_line_g[o] for o in tri)
            pr = (sd + 0.5) / (sg + 1.0); ok = okc4 & (sg >= 3)
            a, _ = acc(pr > 0.5, ok); a0, _ = acc(c4dn, ok)
            dtri.append(a - a0)
            # module readout: majority sign per module from the triple (count-weighted)
            msign = np.full(module.max() + 1, np.nan); c4maj = np.full(module.max() + 1, np.nan)
            for k in range(module.max() + 1):
                mk = module == k
                if sg[mk].sum() >= 20:
                    msign[k] = sd[mk].sum() / sg[mk].sum() > 0.5
                    c4maj[k] = np.nanmean(C4[mk]) > 0.5
            mod_dn = np.where(module >= 0, msign[np.maximum(module, 0)], np.nan)
            okm = okc4 & np.isfinite(mod_dn)
            a, _ = acc(mod_dn > 0.5, okm); a0, _ = acc(c4dn, okm)
            dmod.append(a - a0)
            flip = np.where(module >= 0, (msign != c4maj)[np.maximum(module, 0)], False) & np.isfinite(mod_dn)
            f_dn = np.where(flip, mod_dn > 0.5, c4dn)
            a, _ = acc(f_dn, okm); a0, _ = acc(c4dn, okm)
            dflip.append(a - a0)
        for k, v in (("tri", dtri), ("tri_module16", dmod), ("tri_module16_flip", dflip)):
            res[k] = dict(mean=float(np.mean(v)), p05=float(np.quantile(v, .05)), p95=float(np.quantile(v, .95)))
        out[L] = res
        print(f"{L:9s} full-11 {res['full11_d']:+.3f} | 3-line {res['tri']['mean']:+.3f} [{res['tri']['p05']:+.3f},{res['tri']['p95']:+.3f}] "
              f"| 16-module {res['tri_module16']['mean']:+.3f} | C4 flipped by module {res['tri_module16_flip']['mean']:+.3f}", flush=True)
    summ = {k: float(np.mean([out[L][k]["mean"] if isinstance(out[L][k], dict) else out[L][k] for L in lines]))
            for k in ("full11_d", "tri", "tri_module16", "tri_module16_flip")}
    out["_summary_mean_over_lines"] = summ
    print("\nmean over 12 lines:", {k: round(v, 4) for k, v in summ.items()})
    (HERE / ("work/lab3_622_xbatch.json" if CROSS else "work/lab3_622.json")).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
