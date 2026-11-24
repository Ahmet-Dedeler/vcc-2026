#!/usr/bin/env python3
"""517 — the complement hook ALONE on the incumbent's block (base c4 = CLEAN4 q80), with Arc's programme attenuated by a
random sign-flip fraction r (grid517's emulation of a cross-line quality loss; the same fixed mask), against 496's CLEAN4
hook. Question: how much programme quality does the hook need to beat the CLEAN4 hook, and does it need the block to have
switched? (The hook replaces the DESIGN's sign on uncalled genes, ~0.53 accurate, a lower bar than CLEAN4's on the block.)
Also records each source's accuracy on the uncalled (complement) gate genes, the hook's own territory.
usage: uv run python experiments/517-flex-programme/hookgrid517.py      (< 2 GB)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import price517 as P  # noqa: E402
RC = P.RC


def main():
    t0 = time.monotonic()
    zr = np.load(RC.REF, allow_pickle=False)
    genes_full = zr["genes"]; keep = np.asarray(zr["keep"], dtype=bool); T = zr["targets"].astype(str)
    zr.close()
    SRC, _ = P.load_sources(genes_full, keep)
    Pc4, Parc = SRC["c4"], SRC["arc"]
    c4 = RC.prior_signs_full(genes_full)[keep]
    strong, _ = P.strong_set(T)
    u = np.random.default_rng(51717).random(Parc.size)       # grid517's mask (same seed, same draw order)

    def hs(Pv, tau):
        return np.where(np.isfinite(Pv) & (np.abs(Pv - 0.5) >= tau), np.where(Pv > 0.5, -1.0, 1.0), 0.0)
    hc = hs(Pc4, 0.05)
    rows = []
    for seed in (0, 1, 2):
        ctx = P.Ctx("c4", seed, lambda i, t: c4)
        rng = np.random.default_rng(5173 + seed)
        lfc, _ = ctx.hook(lambda i, t: hc)
        sc, _ = ctx.score(lfc, rng, n_boot=500)
        acc = ctx.gate_accuracy(lambda i, t: hc, "comp")
        rows.append({"seed": seed, "src": "clean4", "r": None, "tau": 0.05, "gated": False, **sc, "acc": acc})
        print(f"s{seed} clean4 tau .05        nmae {sc['nmae']:+.5f} reach {sc['reach']:+.5f} fid {sc['fid']:+.5f} ov {sc['overall_equiv']:+.5f} "
              f"| acc {acc['acc']:.4f} n {acc['n']}", flush=True)
        for r in (0.0, 0.1, 0.2, 0.3, 0.4):
            Pa = np.where(u < r, 1.0 - Parc, Parc)
            for tau in (0.05, 0.1, 0.15, 0.2):
                ha = hs(Pa, tau)
                for gated in (False, True):
                    fn = (lambda i, t, ha=ha: ha if strong[i] else hc) if gated else (lambda i, t, ha=ha: ha)
                    lfc, rec = ctx.hook(fn)
                    sc, _ = ctx.score(lfc, rng, n_boot=500)
                    acc = ctx.gate_accuracy(fn, "comp")
                    rows.append({"seed": seed, "src": "arc", "r": r, "tau": tau, "gated": gated, **sc, "acc": acc,
                                 "flips": float(np.median(rec["flips"]))})
                    print(f"s{seed} arc r{r:.1f} tau {tau:<4} {'g33' if gated else 'all'} nmae {sc['nmae']:+.5f} reach {sc['reach']:+.5f} "
                          f"fid {sc['fid']:+.5f} ov {sc['overall_equiv']:+.5f} ne {sc['ne_overall_equiv']:+.5f} | acc {acc['acc']:.4f} n {acc['n']} "
                          f"| flips {np.median(rec['flips']):.0f} ({time.monotonic() - t0:.0f}s)", flush=True)
        del ctx
    (HERE / "work/hookgrid517.json").write_text(json.dumps(rows, indent=1, default=float))
    print(f"done {time.monotonic() - t0:.0f}s")


if __name__ == "__main__":
    main()
