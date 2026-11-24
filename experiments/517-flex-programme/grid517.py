#!/usr/bin/env python3
"""517 — the H1 dose / protection / attenuation grid for Arc's programme, analytic (price517: validated against 462l / 467 /
470's real runs by blockval517.py: nmae to 5e-5, fid to ~10 %, reach to ~0.002; jac (leaked calls) is NOT seen).

From the shipped c4 base (the incumbent's DE construction), per variant:
  block prior  = sign(w z(Arc_r) + (1 - w) z(CLEAN4)), gated (Arc-family on 467's strong 57, CLEAN4 on the weak 28) or not;
  protect_q    in {0.6, 0.7, 0.8, 0.9} (the shipped 0.8);
  Arc_r        = Arc's P_down with a random fraction r of genes' signs flipped (P -> 1 - P; one fixed mask), an emulation of
                 a cross-line loss of programme quality: on A/B/C the programme is another line's, so the board's 458d read
                 places A/B/C at some r;
  + optionally the complement hook from the same source (gated like the block), tau on |P - 0.5|.
Output work/grid517.json, work/grid517.log. usage: uv run python experiments/517-flex-programme/grid517.py [--quick]
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import price517 as P  # noqa: E402
RC = P.RC


def main():
    quick = "--quick" in sys.argv
    t0 = time.monotonic()
    zr = np.load(RC.REF, allow_pickle=False)
    genes_full = zr["genes"]; keep = np.asarray(zr["keep"], dtype=bool); T = zr["targets"].astype(str)
    zr.close()
    SRC, _ = P.load_sources(genes_full, keep)
    Pc4, Parc = SRC["c4"], SRC["arc"]
    c4 = RC.prior_signs_full(genes_full)[keep]
    strong, thr = P.strong_set(T)
    mask_rng = np.random.default_rng(51717)
    u = mask_rng.random(Parc.size)            # one fixed uniform per gene: flip where u < r (nested masks across r)
    R_LIST = (0.0, 0.1, 0.2, 0.3, 0.4) if not quick else (0.0, 0.3)
    W_LIST = (1.0, 0.75, 0.5, 0.25)
    Q_LIST = (0.5, 0.6, 0.7, 0.8, 0.9)
    TAUS = (None, 0.0, 0.05, 0.1, 0.15)

    def src(r, w):
        Pa = np.where(u < r, 1.0 - Parc, Parc)
        return Pa if w == 1.0 else P.zblend(Pa, Pc4, w)

    def hooksig(Pv, tau):
        return np.where(np.isfinite(Pv) & (np.abs(Pv - 0.5) >= tau), np.where(Pv > 0.5, -1.0, 1.0), 0.0)
    out = {"gate_thr": thr, "n_strong": int(strong.sum()), "rows": []}
    for seed in (0, 1, 2):
        ctx = P.Ctx("c4", seed, lambda i, t: c4)
        rng = np.random.default_rng(5171 + seed)
        # accuracy of each source on H1's gate genes (block: fid-like territory; complement: the hook's territory)
        for r, w in itertools.product(R_LIST, W_LIST):
            Pv = src(r, w)
            blk_sign = P.sign_from_P(Pv)
            for gated, q in itertools.product((True, False), Q_LIST):
                if w != 1.0 and q not in (0.7, 0.8):
                    continue
                pr = (lambda i, t, b=blk_sign: b if strong[i] else c4) if gated else (lambda i, t, b=blk_sign: b)
                for tau in TAUS:
                    if tau is not None and (q != 0.8 or w not in (1.0, 0.5)):
                        continue
                    hk = None
                    if tau is not None:
                        hs = hooksig(Pv, tau)
                        hc = hooksig(Pc4, 0.05)
                        hk = (lambda i, t, hs=hs, hc=hc: hs if strong[i] else hc) if gated else (lambda i, t, hs=hs: hs)
                    lfc_new, rec = ctx.block(pr, protect_q=q, hook_for_target=hk)
                    sc, _ = ctx.score(lfc_new, rng, n_boot=400)
                    row = {"seed": seed, "r": r, "w": w, "gated": gated, "q": q, "tau": tau,
                           **{m: sc[m] for m in P.SP}, "ov": sc["overall_equiv"], "ci95": sc["ci95"],
                           "ne_ov": sc["ne_overall_equiv"], "ne_nmae": sc["ne_nmae"], "ne_fid": sc["ne_fid"],
                           "without_four": sc["without_four"], "flips_vs_base": float(np.median(rec["flips_vs_base"]))}
                    out["rows"].append(row)
                    print(f"s{seed} r{r:.1f} w{w:.2f} {'g33' if gated else 'all'} q{q:.1f} tau {tau} | fid {sc['fid']:+.5f} "
                          f"reach {sc['reach']:+.5f} nmae {sc['nmae']:+.5f} | ov {sc['overall_equiv']:+.5f} ne {sc['ne_overall_equiv']:+.5f}"
                          f" | flips {row['flips_vs_base']:.0f} ({time.monotonic() - t0:.0f}s)", flush=True)
        del ctx
    # accuracies (seed-independent): block genes (fid-like, all block genes vs the reference lfc sign) and gate genes
    (HERE / "work/grid517.json").write_text(json.dumps(out, indent=1))
    print(f"done {time.monotonic() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
