#!/usr/bin/env python3
"""657 post-hoc (written after the 618 read, NOT preregistered): is octile 4's Broad-16 "UP is right" a property of cancer
lines, or a calibration offset of the Broad KO data?

1. Gate-cell level (625's octflip625b.json): per CLEAN4 octile, the fraction of gate cells where UP is right, on the Broad-16
   KO lines (mean of 12) and on H1 (CRISPRi Flex, Arc's chemistry). UP-right = agreement on UP-signed octiles 0-2,
   1 - agreement on DOWN-signed octiles 4-7. Octile 3 straddles P_down 0.5 (mixed signs), so only its agreement is shown.
2. Gene level (371's consensus_prior.npz, P_banks = each bank's per-gene P(down)): in each octile band of a prior that
   EXCLUDES the bank, the bank's mean P_down and the fraction of genes that lean UP in it (P_down < 0.5). CLEAN4's own
   banks (K562-GWPS, K562-ess, Nadig Jurkat, Nadig HepG2) are read against the mean of the other three (LOLO-3); the
   three non-CLEAN4 banks (H1 2025, Orion HCT116, Orion HEK293T) against CLEAN4. Bands use 625's v8 octile edges.
   Every one of these banks is CRISPRi; K562, Jurkat, HepG2 and HCT116 are cancer lines.
Light (two small npz/json files). usage: uv run python experiments/657-octile4-up-callneutral/calib657.py"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EDGES = [-np.inf, 0.4403, 0.465, 0.4855, 0.5058, 0.5251, 0.5481, 0.583, np.inf]   # 625 v8 octile edges (CLEAN4 P_down)
CLEAN4 = ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")


def main():
    out, lines = {}, []
    say = lambda x="": (print(x), lines.append(x))  # noqa: E731
    fl = json.loads((ROOT / "experiments/625-module-readout/work/octflip625b.json").read_text())
    say("1. gate cells (octflip625b): UP-right fraction per CLEAN4 octile (octile 3: CLEAN4-sign agreement, mixed signs)")
    say("   octile | Broad-16 KO mean of 12 (min..max over lines) | H1 (CRISPRi Flex)")
    gate = {}
    for k in range(8):
        r = fl[str(k)]
        per = np.array([v[0] for v in r["per_line"].values()])
        h1a = r["h1"][0]                                   # [count agreement, |lfc|-weighted agreement, n cells]
        if k <= 3:
            b, lo, hi, h = r["mean_agree"], per.min(), per.max(), h1a
        else:
            b, lo, hi, h = 1 - r["mean_agree"], 1 - per.max(), 1 - per.min(), 1 - h1a
        gate[k] = dict(broad=float(b), broad_min=float(lo), broad_max=float(hi), h1=float(h))
        say(f"   o{k}{'*' if k == 3 else ' '}    | {b:.3f} ({lo:.3f}..{hi:.3f}) | {h:.3f}")
    out["gate_cells"] = gate
    c4 = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(c4["bank_names"].astype(str))
    P = c4["P_banks"].astype(float)
    say("")
    say("2. gene level (371 P_banks): bank's mean P_down [fraction of genes leaning UP] in each octile band of a prior without it")
    gene = {}
    for b in names:
        prior = [x for x in CLEAN4 if x != b]
        ref = P[[names.index(x) for x in prior]].mean(0)
        pb = P[names.index(b)]
        rows = []
        for k in range(8):
            m = (ref >= EDGES[k]) & (ref < EDGES[k + 1]) & np.isfinite(pb) & np.isfinite(ref)
            rows.append(dict(mean_p_down=float(pb[m].mean()), frac_up=float((pb[m] < 0.5).mean()), n=int(m.sum())))
        gene[b] = dict(prior="LOLO-3" if b in CLEAN4 else "CLEAN4", global_mean_p_down=float(np.nanmean(pb)), bands=rows)
        say(f"   {b:15s} vs {gene[b]['prior']:6s} (global {gene[b]['global_mean_p_down']:.3f}) | "
            + " ".join(f"o{k}:{r['mean_p_down']:.3f}[{r['frac_up']:.2f}]" for k, r in enumerate(rows)))
    out["gene_level"] = gene
    o4 = {b: gene[b]["bands"][4]["frac_up"] for b in names}
    say("")
    say("   octile 4, fraction of genes leaning UP: " + ", ".join(f"{b} {v:.2f}" for b, v in o4.items()))
    (HERE / "work/calib657.json").write_text(json.dumps(out, indent=1))
    (HERE / "work/calib657.out").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
