#!/usr/bin/env python3
"""657 — the prereg's remaining light measures, read from the stats npz files only (no lib521 state, a few MB per file):
  1. zero-total gated cells (lfc < -30, zerocount654.py's rule), split into own-target and non-own cells, per arm and seed;
  2. mean |delivered lfc| on uncalled (q >= 0.05) octile-4 cells, non-own: over all such cells, and over those inside the
     reference's own call set (ref q < 0.05; a proxy for analyze657's f 1 pool, which needs lib521 state).
Gated axis = reference genes[keep] (618's ref_tv); own = the target's own gene on that axis.
usage: uv run python experiments/657-octile4-up-callneutral/zero657.py"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REF = np.load(ROOT / "experiments/618-h1-trainval-surface/work/ref_tv/reference.npz", allow_pickle=True)
G = REF["genes"].astype(str)[REF["keep"].astype(bool)]
TG = REF["targets"].astype(str)
RQ = np.asarray(REF["q"], float)
V8 = np.load(ROOT / "experiments/625-module-readout/work/modules625_v8.npz", allow_pickle=True)
OCT = dict(zip(V8["genes"].astype(str), V8["octile"].astype(int)))
O4 = np.array([OCT.get(g, -1) == 4 for g in G])
EM = "stats-base_s{s}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz"
ST = HERE / "work/stacks"


def arms(s):
    st = lambda nm, v: ST / f"stats-rnup_{nm}_s{s}_{v}.npz"  # noqa: E731
    return {"em_c4h": HERE / "work/c4h" / EM.format(s=s), "A0": HERE / "work/c4h_o4" / EM.format(s=s),
            "base": st("c4h", "b_k100_q20r3p10"), "A": st("o4", "b_k100_q20r3p10"),
            "B": st("c4h", "b_k100_q20r3p10_o4"), "AB": st("o4", "b_k100_q20r3p10_o4"),
            "base_dn": st("c4h", "b_k100_q20r3p10_D200"), "A_dn": st("o4", "b_k100_q20r3p10_D200"),
            "B_dn": st("c4h", "b_k100_q20r3p10_D200_o4"), "AB_dn": st("o4", "b_k100_q20r3p10_D200_o4")}


def main():
    out, lines = {}, []
    say = lambda x="": (print(x), lines.append(x))  # noqa: E731
    say(f"gated genes {len(G)}, octile-4 on the gated axis {int(O4.sum())}, targets {len(TG)}")
    say("arm | per seed: zero-total cells non-own (own) | mean |lfc| uncalled o4 non-own: all / in ref call set")
    for arm in arms(0):
        rec = {}
        for s in (0, 1, 2):
            z = np.load(arms(s)[arm], allow_pickle=True)
            assert np.array_equal(z["targets"].astype(str), TG), (arm, s)
            lfc, q = np.asarray(z["lfc"], float), np.asarray(z["q"], float)
            assert lfc.shape == (len(TG), len(G)), (arm, s, lfc.shape)
            own = G[None, :] == TG[:, None]
            zero = lfc < -30
            m = (q >= 0.05) & O4[None, :] & ~own & np.isfinite(lfc)
            mc = m & (RQ < 0.05)
            rec[s] = dict(zero_nonown=int((zero & ~own).sum()), zero_own=int((zero & own).sum()),
                          o4_abs_lfc_all=float(np.abs(lfc[m]).mean()), o4_abs_lfc_refcalled=float(np.abs(lfc[mc]).mean()),
                          n_all=int(m.sum()), n_refcalled=int(mc.sum()))
        out[arm] = rec
        say(f"{arm:7s} | " + " | ".join(f"s{s}: {r['zero_nonown']} ({r['zero_own']}) {r['o4_abs_lfc_all']:.4f} / "
                                        f"{r['o4_abs_lfc_refcalled']:.4f}" for s, r in rec.items()))
    (HERE / "work/zero657.json").write_text(json.dumps(out, indent=1))
    (HERE / "work/zero657.out").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
