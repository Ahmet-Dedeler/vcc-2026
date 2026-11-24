#!/usr/bin/env python3
"""537 — held-out check for the block + CLEAN4 hook synergy (P1): on public cancer-line truths, how often does CLEAN4's
sign (tau 0.05) match the truth on the genes the call-propensity block swaps IN (the genes the hook reaches under `kept`)?

On H1 the hook's gain on the block comes from these called-but-unsteered genes (design sign 0.508 all / 0.533 gate). If
CLEAN4 carries sign information on them in cancer lines too, the fid part of the synergy should transfer.

Folds (480's gates; RELATIVE DIRECTION ONLY for the Orion truths, which carry a target-independent blob, AGENTS 6.5):
k562_gwps (CLEAN4 minus its two K562 banks = Nadig Jurkat + HepG2: the held line's own banks never vote), orion_hct116,
orion_hek293t (CLEAN4's four banks; no Orion bank is in CLEAN4). Blocks: 528's fold528b `norm` key (lc + 22 x median-rescaled
to H1) vs the CPM block, 2,200 genes among non-panel genes. Gate = q < .05, finite non-zero lfc, not the target.
Reported per fold: CLEAN4 tau .05 coverage and sign accuracy on gate genes that are swapped in / kept block / swapped out /
the rest, the truth's down fraction on those genes (the all-down baseline), and H1's same numbers (h1-test reference).
Known answer: the swap counts equal fold528b's `norm22` rows (work/fold528b.json). Writes work/foldsign537.json. < 1 GB.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
E528 = ROOT / "experiments/528-callprop-block"
GATES = ROOT / "experiments/480-nmae-instrument/work"
BETA = 22.0


def block(score, excl, n=2200):
    b = np.zeros(score.size, bool)
    b[np.argsort(-np.where(excl, -np.inf, score), kind="stable")[:n]] = True
    return b


def clean4_P(genes, banks):
    cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(cp["bank_names"].astype(str))
    P = cp["P_banks"].astype(float)[[names.index(b) for b in banks]]
    with np.errstate(invalid="ignore"):
        Pm = np.nanmean(np.where(np.isfinite(P), P, np.nan), axis=0)
    pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
    col = np.array([pos.get(str(g), -1) for g in genes])
    out = np.full(len(genes), np.nan)
    out[col >= 0] = Pm[col[col >= 0]]
    return out


def stats(gate, lfc, P, cls, adj=None):
    ok = np.isfinite(P) & (np.abs(P - 0.5) >= 0.05)
    sgn = np.where(P > 0.5, -1.0, 1.0)
    out = {}
    for name, M in cls.items():
        G = gate & M[None, :]
        cov = G & ok[None, :]
        n, nc = int(G.sum()), int(cov.sum())
        tr = np.sign(lfc)
        acc = float((tr[cov] == np.broadcast_to(sgn, G.shape)[cov]).mean()) if nc else float("nan")
        down = float((tr[G] < 0).mean()) if n else float("nan")
        down_cov = float((tr[cov] < 0).mean()) if nc else float("nan")
        out[name] = dict(gate_pairs=n, covered=nc, coverage=nc / max(n, 1), clean4_acc=acc, truth_down=down,
                         truth_down_on_covered=down_cov)
        if adj is not None:     # fid's population: every adjudicable pair, significant or not
            A = adj & M[None, :]
            out[name]["up_all_adjudicable"] = float((lfc[A] > 0).mean()) if A.any() else float("nan")
    return out


def main():
    H = dict(np.load(E528 / "work/feat528_H1.npz", allow_pickle=True))
    med_h1 = float(np.median(H["x"][H["keep"]]))
    ka = json.loads((E528 / "work/fold528b.json").read_text())
    res = {"median_x_H1": med_h1}
    for fold in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        F = dict(np.load(E528 / f"work/feat528_{fold}.npz", allow_pickle=True))
        z = np.load(GATES / f"gate-{fold}.npz", allow_pickle=True)
        gg = z["genes"].astype(str)
        pos = {g: i for i, g in enumerate(F["genes"].astype(str))}
        col = np.array([pos[g] for g in gg])
        x = F["x"][col]
        lc = np.log10(np.maximum(F["mean_cpm"][col], 1e-3))
        xn = x * (med_h1 / float(np.median(F["x"][F["keep"]])))
        targets = z["targets"].astype(str)
        q, lfc = np.asarray(z["q"], float), np.asarray(z["lfc"], float)
        gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & ~np.equal.outer(targets, gg)
        excl = np.isin(gg, targets)
        b0, b1 = block(lc, excl), block(lc + BETA * xn, excl)
        sw_in, sw_out, kept = b1 & ~b0, b0 & ~b1, b0 & b1
        rest = ~b0 & ~b1 & ~excl
        # known answer: the swap count equals fold528b's norm row
        rec = ka.get(fold, {}).get("rows", {}).get("norm22", {})
        ka_ok = (rec.get("swapped") == int(sw_in.sum())) if rec else None
        banks = ("nadig_jurkat", "nadig_hepg2") if fold == "k562_gwps" else ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")
        P = clean4_P(gg, banks)
        adj = np.isfinite(lfc) & (lfc != 0) & ~np.equal.outer(targets, gg)
        st = stats(gate, lfc, P, {"swapped_in": sw_in, "kept_block": kept, "swapped_out": sw_out, "rest": rest}, adj)
        res[fold] = dict(banks=list(banks), swapped_in=int(sw_in.sum()), ka_swap_matches_fold528b=ka_ok, classes=st)
        print(f"{fold} (CLEAN4 = {'+'.join(banks)}; swapped in {int(sw_in.sum())}; KA swap == fold528b norm22: {ka_ok})")
        for k, v in st.items():
            print(f"   {k:12s} gate pairs {v['gate_pairs']:7d} coverage {v['coverage']:.3f} CLEAN4 acc {v['clean4_acc']:.3f} "
                  f"| truth down {v['truth_down']:.3f} (on covered {v['truth_down_on_covered']:.3f}) | up over all adjudicable "
                  f"{v['up_all_adjudicable']:.3f}")
    # H1 (h1-test reference; the same block construction on its own controls; CLEAN4's four banks)
    ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
    keep = ref["keep"].astype(bool)
    gk = ref["genes"].astype(str)[keep]
    pos = {g: i for i, g in enumerate(H["genes"].astype(str))}
    colh = np.array([pos[g] for g in gk])
    lc = np.log10(np.maximum(H["mean_cpm"][colh], 1e-3))
    x = H["x"][colh]
    targets = ref["targets"].astype(str)
    q, lfc = np.asarray(ref["q"], float), np.asarray(ref["lfc"], float)
    gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & ~np.equal.outer(targets, gk)
    excl = np.isin(gk, targets)
    b0, b1 = block(lc, excl), block(lc + BETA * x, excl)
    P = clean4_P(gk, ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2"))
    adj = np.isfinite(lfc) & (lfc != 0) & ~np.equal.outer(targets, gk)
    st = stats(gate, lfc, P, {"swapped_in": b1 & ~b0, "kept_block": b0 & b1, "swapped_out": b0 & ~b1, "rest": ~b0 & ~b1 & ~excl}, adj)
    res["H1"] = dict(swapped_in=int((b1 & ~b0).sum()), classes=st)
    print(f"H1 (swapped in {int((b1 & ~b0).sum())}; no head exclusion here)")
    for k, v in st.items():
        print(f"   {k:12s} gate pairs {v['gate_pairs']:7d} coverage {v['coverage']:.3f} CLEAN4 acc {v['clean4_acc']:.3f} "
              f"| truth down {v['truth_down']:.3f} (on covered {v['truth_down_on_covered']:.3f}) | up over all adjudicable "
              f"{v['up_all_adjudicable']:.3f}")
    (HERE / "work/foldsign537.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
