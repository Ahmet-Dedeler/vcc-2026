#!/usr/bin/env python3
"""625 — preregistration support: simulate the A/B/C MODULE READOUT and its D/E/F transfer rule on the Broad-16.

What the board readout measures (nmae-primary decode): per module M, the net |l|-weighted agreement of CLEAN4's sign with
the reference sign on GATE genes, pooled over the contexts (A/B/C) and targets: a_M = sum (agree - disagree) |l| / sum |l|.
Simulation (Broad-16, 12 working lines; reserved 4 never read): "A/B/C" = 3 random other lines; modules = 16 FIXED gene sets
(CLEAN4 octile x control-expression half, the expression pooled over the 3 training lines' controls: truth-free); compute
a_M on the triple's gated genes; the transfer rule flips CLEAN4 on module M for the held-out line iff a_M < -bar; read the
held-out line's gate-sign accuracy vs CLEAN4 (plain and |l|-weighted). 40 triples per held-out line. Bars 0 ... 0.2.
Also: the noise floor of a_M on a single "control" (sign-randomised truth) to see where a bar must sit.
Analysis only, < 4 GB.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/622-identity-track"))
import broad622 as B  # noqa: E402

RNG = np.random.default_rng(6250)
BARS = (0.0, 0.02, 0.05, 0.1, 0.15, 0.2)


def main() -> None:
    gl, LFC, cl = B.load_matrix(B.P / "lfc_matrix.csv", 1)
    gf, FDR, cf = B.load_matrix(B.P / "fdr_matrix.csv", 1)
    gp, PB, cp = B.load_matrix(B.P / "pseudobulk_sum_by_ko.csv", 0)
    lines = sorted({ln for ln, _ in cl}); assert len(lines) == 12 and not (set(lines) & B.RESERVED)
    c4 = np.load(B.ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    nm = list(c4["bank_names"].astype(str))
    Pc = c4["P_banks"].astype(float)[[nm.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    cmap = dict(zip(c4["genes"].astype(str), Pc)); C4 = np.array([cmap.get(g, np.nan) for g in gl])
    c4dn = C4 > 0.5
    all_t = set(k for _, k in cl)
    gpos = {g: i for i, g in enumerate(gp)}
    rowmap = np.array([gpos.get(g, -1) for g in gl])
    ctrl = {}
    for ln in lines:
        cc = [j for j, (l2, k) in enumerate(cp) if l2 == ln and k not in all_t]
        s = np.nansum(PB[:, cc], 1); v = 1e6 * s / np.nansum(s)
        ctrl[ln] = np.where(rowmap >= 0, v[np.maximum(rowmap, 0)], np.nan)
    octile = np.digitize(C4, np.nanquantile(C4, np.linspace(0.125, 0.875, 7)))
    gated = (FDR < 0.05) & np.isfinite(LFC) & (LFC != 0)
    colline = np.array([ln for ln, _ in cl])
    tgt = np.array([k for _, k in cl])
    notown = gl[:, None] != tgt[None, :]
    agree = np.where(gated & notown, np.where((LFC < 0) == c4dn[:, None], 1.0, -1.0), 0.0)
    absl = np.where(gated & notown, np.abs(LFC), 0.0)
    out = {"bars": {}, "noise": {}}
    res = {b: {"plain": [], "w": []} for b in BARS}
    flips_per_bar = {b: [] for b in BARS}
    noise_aM = []
    for L in lines:
        others = [o for o in lines if o != L]
        jl = np.flatnonzero(colline == L)
        gL = gated[:, jl] & notown[:, jl]; tdnL = LFC[:, jl] < 0; aL = np.abs(np.nan_to_num(LFC[:, jl]))
        okc = np.isfinite(C4)
        for _ in range(40):
            tri = list(RNG.choice(others, 3, replace=False))
            ex = np.nanmean([np.log1p(ctrl[o]) for o in tri], 0)
            half = (ex > np.nanmedian(ex)).astype(int)
            module = np.where(okc & np.isfinite(ex), octile * 2 + half, -1)
            jt = np.flatnonzero(np.isin(colline, tri))
            num = (agree[:, jt] * absl[:, jt]).sum(1); den = absl[:, jt].sum(1)
            # sign-randomised control (noise floor): random signs on the same genes/weights
            rnd = (np.where(RNG.random(agree[:, jt].shape) < 0.5, 1.0, -1.0) * (absl[:, jt] > 0) * absl[:, jt]).sum(1)
            aM = np.full(16, np.nan); aN = np.full(16, np.nan)
            for k in range(16):
                mk = module == k
                if den[mk].sum() > 0:
                    aM[k] = num[mk].sum() / den[mk].sum(); aN[k] = rnd[mk].sum() / den[mk].sum()
            noise_aM.extend(list(aN[np.isfinite(aN)]))
            for b in BARS:
                flipmod = np.isfinite(aM) & (aM < -b)
                flipg = np.where(module >= 0, flipmod[np.maximum(module, 0)], False)
                sgn = np.where(flipg, ~c4dn, c4dn)
                m = gL & okc[:, None]
                base = ((c4dn[:, None] == tdnL) & m); new = ((sgn[:, None] == tdnL) & m)
                res[b]["plain"].append(float((new.sum() - base.sum()) / max(m.sum(), 1)))
                res[b]["w"].append(float(((new.astype(float) - base) * aL * m).sum() / max((aL * m).sum(), 1e-12)))
                flips_per_bar[b].append(int(flipmod.sum()))
    for b in BARS:
        v = np.array(res[b]["plain"]); w = np.array(res[b]["w"])
        per_line = v.reshape(12, 40).mean(1)
        out["bars"][str(b)] = dict(mean=float(v.mean()), p05=float(np.quantile(v, .05)), mean_w=float(w.mean()),
                                   min_line=float(per_line.min()), lines_pos=int((per_line > 0).sum()),
                                   frac_triples_neg=float((v < 0).mean()), mean_modules_flipped=float(np.mean(flips_per_bar[b])))
        r = out["bars"][str(b)]
        print(f"bar {b:4.2f}: Δ gate-sign {r['mean']:+.4f} (|l|-w {r['mean_w']:+.4f}); worst line {r['min_line']:+.4f}; lines>0 {r['lines_pos']}/12; "
              f"triples<0 {r['frac_triples_neg']:.2f}; modules flipped {r['mean_modules_flipped']:.1f}/16", flush=True)
    nz = np.abs(np.array(noise_aM))
    out["noise"] = dict(abs_aM_noise_p50=float(np.median(nz)), p95=float(np.quantile(nz, .95)), p99=float(np.quantile(nz, .99)))
    print("noise floor |a_M| (random signs):", {k: round(v, 4) for k, v in out["noise"].items()})
    (HERE / "work").mkdir(exist_ok=True)
    (HERE / "work/sim625.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
