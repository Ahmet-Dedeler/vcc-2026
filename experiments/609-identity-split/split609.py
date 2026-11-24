#!/usr/bin/env python3
"""609 ([4d0207] for [7772dc], council 600, 2026-10-02 ~14:10Z): is the 639d (block) vs 622m (uncalled) conflict over WHERE
the identity rider's sign gain sits definitional or real?

Part B (this file): 639's four folds (h1, hct116, hek293t, k562) with 639d's OWN masks (decomp639.callset_for: head 400
by the 3-donor consensus key + the 2,200-gene CPM block; h1 = the real design key) and 639d's identity I (Flex refs).
Per region (head / block / uncalled), the paired gain over the matching CLEAN4 baseline, net (right - wrong) on gate
genes, for these sources:
  rider        zc + 0.5 zi, zc = (P - 1/2)/sd UNCENTRED (639d's and production v3's rule) -- KA: equals decomp639.json
  rider_c      622m's convention: zc_c = z(P - 1/2) mean-CENTRED; the baseline is sign(zc_c) too (as block622 does)
  rider_expr   zc + 0.5 z(E), E = log1p(own control CPM) on the genes with a finite I (an expression-level control)
  rider_IperpE zc + 0.5 z(resid(I | E)) (identity with expression partialled out, OLS on the finite pairs)
plus each region's gate-gene DOWN fraction (truth lfc < 0): the Orion folds' truth carries a technical blob (AGENTS 6.5;
docs/HARNESS-DEPTH-CONFOUND.md: a depth mismatch skews DE calls negative), so a large DOWN skew in the block would make
"identity genes go down" partly an artifact there.
usage: uv run python experiments/609-identity-split/split609.py   -> work/split609.json"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/639-identity-axis-lolo"))
import decomp639 as D  # noqa: E402

LAM = 0.5


def resid(y, x):
    ok = np.isfinite(y) & np.isfinite(x)
    b = np.polyfit(x[ok], y[ok], 1)
    return np.where(ok, y - np.polyval(b, x), np.nan)


def z(v):
    return (v - np.nanmean(v)) / np.nanstd(v)


def main() -> None:
    I = D.Inst()
    prof = D.flex_profiles()
    ka = json.loads((ROOT / "experiments/639-identity-axis-lolo/work/decomp639.json").read_text())
    out = {}
    for s in I.S:
        X = I.S[s]
        T, G, lfc, gate = X["T"], X["G"], X["lfc"], X["gate"]
        if s == "h1":
            own = prof["h1"]
            refs = [prof[k] for k in ("ctxA", "ctxB", "ctxC")]
        else:
            zz = np.load(D.G480 / D.SURF[s], allow_pickle=True)
            own = dict(zip(zz["genes"].astype(str), zz["ctrl_cpm"].astype(float)))
            refs = [prof[k] for k in ("ctxA", "ctxB", "ctxC", "h1")]
        Iv = D.identity(own, refs, list(G))
        zi = (Iv - np.nanmean(Iv)) / np.nanstd(Iv)
        Ev = np.array([np.log1p(own[g]) if (g in own and np.isfinite(Iv[j])) else np.nan for j, g in enumerate(G)])
        rIE = float(np.corrcoef(Iv[np.isfinite(Iv) & np.isfinite(Ev)], Ev[np.isfinite(Iv) & np.isfinite(Ev)])[0, 1])
        C4, _ = D.clean4_lolo(I, s)
        zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
        zc_c = z(C4 - 0.5)
        head, block = D.callset_for(I, s)
        unc = ~head & ~block
        dn = lfc < 0
        srcs = {
            "rider": (zc + LAM * zi, zc),
            "rider_c": (zc_c + LAM * zi, zc_c),
            "rider_expr": (zc + LAM * z(Ev), zc),
            "rider_IperpE": (zc + LAM * z(resid(Iv, Ev)), zc),
        }
        rec = {"corr_I_E": rIE, "regions": {}}
        for rname, R in (("head", head), ("block", block), ("uncalled", unc)):
            g = gate & R
            rr = {"gate_pairs": int(g.sum()), "gate_down_frac": float(dn[g].mean()) if g.any() else float("nan")}
            for name, (sc, base) in srcs.items():
                both = np.isfinite(sc) & np.isfinite(base)
                m = g & both[None, :]
                oks = ((sc > 0)[None, :] == dn) & m
                okb = ((base > 0)[None, :] == dn) & m
                n = int(m.sum())
                rr[name] = dict(n=n, acc_base=float(okb.sum() / max(n, 1)), acc=float(oks.sum() / max(n, 1)),
                                gain=float((oks.sum() - okb.sum()) / max(n, 1)),
                                d_net=float(2 * (int(oks.sum()) - int(okb.sum()))))
            rec["regions"][rname] = rr
        for name in srcs:
            tot = sum(rec["regions"][r][name]["d_net"] for r in ("head", "block", "uncalled"))
            for r in ("head", "block", "uncalled"):
                rec["regions"][r][name]["share"] = rec["regions"][r][name]["d_net"] / tot if tot else float("nan")
            rec[f"total_d_net_{name}"] = tot
        # KA: our `rider` must reproduce 639d's lam0.5 d_net per region exactly
        for r in ("head", "block", "uncalled"):
            want = ka[s]["regions"]["lam0.5"][r]["d_net"]
            got = rec["regions"][r]["rider"]["d_net"]
            if abs(want - got) > 1e-9:
                raise SystemExit(f"KA FAILED {s} {r}: d_net {got} vs 639d {want}")
        out[s] = rec
        print(f"\n== {s}: corr(I, E) {rIE:+.2f}; KA vs 639d: PASS", flush=True)
        print(f"  {'region':9s} {'gate':>6s} {'down':>5s} | " + " | ".join(f"{n:>22s}" for n in srcs))
        for r in ("head", "block", "uncalled"):
            q = rec["regions"][r]
            print(f"  {r:9s} {q['gate_pairs']:6d} {q['gate_down_frac']:5.2f} | " + " | ".join(
                f"gain {q[n]['gain']:+.4f} sh {q[n]['share']:+.2f}" for n in srcs), flush=True)
    (HERE / "work/split609.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
