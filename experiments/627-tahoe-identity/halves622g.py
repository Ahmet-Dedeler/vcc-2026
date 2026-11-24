#!/usr/bin/env python3
"""627d — does 622g's Broad-16 identity sign result survive CONTROL HALVES? (shared-control-noise check)

622g used the authors' lfc (vs ALL 10 OR controls) and I from the SAME controls. Here: gate = authors' fdr < .05 (the gate
set only); the SIGN truth is recomputed against OR half A (5 KOs); I from OR half B (other 5) vs the other 11 lines' half B.
Reads gate-sign accuracy of identity alone, CLEAN4, and z(C4) + 0.5 z(I), per line, against the same reads with I from the
FULL controls (622g's construction) on the same genes. If identity-alone falls to ~0.5 under halves, 622g was shared noise.
Analysis only, < 4 GB.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/622-identity-track"))
import broad622 as B  # noqa: E402


def main() -> None:
    gl, LFC, cl = B.load_matrix(B.P / "lfc_matrix.csv", 1)
    gf, FDR, cf = B.load_matrix(B.P / "fdr_matrix.csv", 1)
    gp, PB, cp = B.load_matrix(B.P / "pseudobulk_sum_by_ko.csv", 0)
    lines = sorted({ln for ln, _ in cl}); assert len(lines) == 12 and not (set(lines) & B.RESERVED)
    all_t = set(k for _, k in cl)
    ors = sorted({k for _, k in cp if k not in all_t}); hA, hB = ors[0::2], ors[1::2]
    col = {(l2, k): j for j, (l2, k) in enumerate(cp)}
    gpos = {g: i for i, g in enumerate(gp)}; rowmap = np.array([gpos.get(g, -1) for g in gl]); okr = rowmap >= 0

    def cpm(ln, kos):
        s = np.nansum(PB[:, [col[(ln, k)] for k in kos if (ln, k) in col]], 1); v = 1e6 * s / np.nansum(s)
        return np.where(okr, v[np.maximum(rowmap, 0)], np.nan)
    cB = {ln: cpm(ln, hB) for ln in lines}; cA = {ln: cpm(ln, hA) for ln in lines}; cF = {ln: cpm(ln, ors) for ln in lines}
    c4 = np.load(B.ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    nm = list(c4["bank_names"].astype(str))
    Pc = c4["P_banks"].astype(float)[[nm.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    cmap = dict(zip(c4["genes"].astype(str), Pc)); C4 = np.array([cmap.get(g, np.nan) for g in gl])
    zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
    colline = np.array([ln for ln, _ in cl]); tgt = np.array([k for _, k in cl])
    out = {}
    for ln in lines:
        others = [o for o in lines if o != ln]
        res = {}
        for name, ctl in (("halves", cB), ("full(622g)", cF)):
            I = np.log1p(ctl[ln]) - np.nanmean([np.log1p(ctl[o]) for o in others], 0)
            I = np.where(ctl[ln] > 5, I, np.nan); zi = (I - np.nanmean(I)) / np.nanstd(I)
            jl = np.flatnonzero(colline == ln)
            gate = (FDR[:, jl] < 0.05) & np.isfinite(LFC[:, jl]) & (LFC[:, jl] != 0) & (gl[:, None] != tgt[jl][None, :])
            if name == "halves":   # sign truth against half A of the controls (independent of I's half B)
                pert = np.stack([cpm(ln, [t]) for t in tgt[jl]], 1)
                tdn = np.log2((pert + 1) / (cA[ln][:, None] + 1)) < 0
            else:
                tdn = LFC[:, jl] < 0
            both = np.isfinite(zc) & np.isfinite(zi)
            m = gate & both[:, None]
            acc = lambda s: float((((s > 0)[:, None] == tdn) & m).sum() / max(m.sum(), 1))
            res[name] = dict(identity=acc(zi), clean4=acc(zc), rider=acc(zc + 0.5 * zi), n=int(m.sum()))
        out[ln] = res
        h, f = res["halves"], res["full(622g)"]
        print(f"{ln:9s} halves: I-alone {h['identity']:.3f} C4 {h['clean4']:.3f} rider {h['rider']:.3f} (Δ {h['rider'] - h['clean4']:+.4f}) | "
              f"full: I-alone {f['identity']:.3f} C4 {f['clean4']:.3f} rider {f['rider']:.3f} (Δ {f['rider'] - f['clean4']:+.4f})")
    summ = {k: float(np.mean([out[l][k]["identity"] for l in lines])) for k in ("halves", "full(622g)")}
    d = {k: float(np.mean([out[l][k]["rider"] - out[l][k]["clean4"] for l in lines])) for k in ("halves", "full(622g)")}
    print(f"\nmean identity-alone acc: halves {summ['halves']:.4f} vs full {summ['full(622g)']:.4f}; mean rider-C4: halves {d['halves']:+.4f} vs full {d['full(622g)']:+.4f}")
    out["_summary"] = {"identity_alone": summ, "rider_minus_c4": d}
    (HERE / "work/halves622g.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
