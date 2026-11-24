#!/usr/bin/env python3
"""634 — DE-conditional direction prior vs CLEAN4's unconditional P_down, on H1's real reference gate.

Question (council 600, d4cd5f; lead from chill's day-1 design "separate DE membership/direction ranking ... signed
conditional response estimates", fid 0.5225 / nmae 0.835 on its first board entry): CLEAN4's per-gene P_down is the
fraction of ALL perturbations in a bank whose delta < 0 (371: P = neg / n_t), i.e. an unconditional programme sign
dominated by non-responding (target, gene) pairs. nmae and reach are scored only on the reference's DE genes, where the
sign that matters is P(down | the gene responds). Does a direction prior estimated from DE EVENTS in public banks predict
H1's gate-gene signs better than CLEAN4?

DE event in a bank: robust z = (delta - median_g) / (1.4826 * MAD_g) over the bank's perturbations, |z| >= Z.
Banks: cancer/other-line public CRISPRi screens on disk (never h1_2025*: it covers h1-test's targets = oracle leakage).
Truth: 323's h1-test reference (85 targets x 10,779 gated genes); gate = q < 0.05, finite non-zero lfc, gene != target.
Analysis only: mmap + chunks, < 4 GB. Writes work/cond634.json and prints the tables.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
BANKS = ["k562_gwps", "orion_hct116", "orion_hek293t", "rpe1", "nadig_jurkat", "k562_essential", "viperturb_k562"]
ZS = (3.0, 4.0, 6.0)
CHUNK = 1500
RNG = np.random.default_rng(634)

# ---- truth: h1-test reference --------------------------------------------------------------------------------------
REF = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
keep = np.asarray(REF["keep"], bool)
G = REF["genes"].astype(str)[keep]
T = REF["targets"].astype(str)
q, lfc = REF["q"], REF["lfc"]
gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & (G[None, :] != T[:, None])
ess = {ln.strip() for ln in open(ROOT / "data/donors/k562_essential/perts.txt") if ln.strip()}
noness_t = np.array([t not in ess for t in T])
gpos = {g: i for i, g in enumerate(G)}
print(f"truth: {len(T)} targets ({noness_t.sum()} non-essential), {len(G)} genes, gate pairs {int(gate.sum())}", flush=True)


def bank_stats(name: str) -> dict:
    d = ROOT / "data/donors" / name
    D = np.load(d / "delta.npy", mmap_mode="r")
    genes = [ln.strip() for ln in open(d / "genes.txt")]
    perts = [ln.strip() for ln in open(d / "perts.txt")]
    n_t, n_g = D.shape
    col = np.array([gpos.get(g, -1) for g in genes])
    use = col >= 0
    cidx = np.flatnonzero(use)
    # robust per-gene centre/scale from a row subsample
    sub = np.sort(RNG.choice(n_t, size=min(n_t, 3000), replace=False))
    S = np.asarray(D[sub][:, cidx], dtype=np.float32)
    med = np.median(S, axis=0)
    mad = np.median(np.abs(S - med), axis=0) * 1.4826
    mad[mad <= 0] = np.inf
    del S
    noness = np.array([p.split("_")[0] not in ess for p in perts])
    out = {"n_perts": n_t, "n_noness": int(noness.sum()), "n_genes_on_axis": int(use.sum())}
    neg_all = np.zeros(cidx.size, np.int64)
    for scope in ("all", "noness"):
        out[scope] = {Z: [np.zeros(cidx.size, np.int64), np.zeros(cidx.size, np.int64)] for Z in ZS}  # [down, up]
    for s in range(0, n_t, CHUNK):
        e = min(s + CHUNK, n_t)
        C = np.asarray(D[s:e][:, cidx], dtype=np.float32)
        neg_all += (C < 0).sum(0)
        Zm = (C - med) / mad
        ne = noness[s:e]
        for Z in ZS:
            dn, up = (Zm <= -Z), (Zm >= Z)
            out["all"][Z][0] += dn.sum(0)
            out["all"][Z][1] += up.sum(0)
            out["noness"][Z][0] += dn[ne].sum(0)
            out["noness"][Z][1] += up[ne].sum(0)
        del C, Zm
    # map to the truth axis
    def full(v, fill=np.nan):
        o = np.full(len(G), fill, dtype=float)
        o[col[cidx]] = v
        return o
    res = {"meta": {k: out[k] for k in ("n_perts", "n_noness", "n_genes_on_axis")},
           "P_uncond": full(neg_all / n_t)}
    for scope in ("all", "noness"):
        for Z in ZS:
            dn, up = out[scope][Z]
            res[f"dn_{scope}_{Z}"] = full(dn, 0)
            res[f"up_{scope}_{Z}"] = full(up, 0)
    print(f"bank {name}: {n_t} perts ({out['n_noness']} non-ess), {int(use.sum())} genes on axis; "
          f"events/gene z4 median {np.median(out['all'][4.0][0] + out['all'][4.0][1]):.0f}", flush=True)
    return res


def clean4_P() -> np.ndarray:
    cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(cp["bank_names"].astype(str))
    P = cp["P_banks"].astype(float)[[names.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
    c = np.array([pos.get(g, -1) for g in G])
    o = np.full(len(G), np.nan)
    o[c >= 0] = P[c[c >= 0]]
    return o


def score(P: np.ndarray, label: str, tmask=None, conf: float = 0.0) -> dict:
    """Sign accuracy of 'down if P > 0.5' on gate pairs where P is defined and |P - .5| >= conf."""
    tm = np.ones(len(T), bool) if tmask is None else tmask
    defined = np.isfinite(P) & (np.abs(P - 0.5) >= conf) & (P != 0.5)
    m = gate & tm[:, None] & defined[None, :]
    pred_dn = (P > 0.5)[None, :]
    truth_dn = lfc < 0
    corr = (pred_dn == truth_dn) & m
    w = np.abs(np.nan_to_num(lfc)) * m
    n_gate = int((gate & tm[:, None]).sum())
    up_m = m & ~truth_dn
    dn_m = m & truth_dn
    r = dict(label=label, coverage=float(m.sum() / max(n_gate, 1)), n=int(m.sum()),
             acc=float(corr.sum() / max(m.sum(), 1)), acc_w=float((w * corr).sum() / max(w.sum(), 1e-12)),
             acc_on_truth_up=float((corr & up_m).sum() / max(up_m.sum(), 1)),
             acc_on_truth_dn=float((corr & dn_m).sum() / max(dn_m.sum(), 1)),
             frac_pred_dn=float((pred_dn & m).sum() / max(m.sum(), 1)))
    return r


def main():
    stats = {b: bank_stats(b) for b in BANKS}
    P4 = clean4_P()
    rows = []
    truth_up_frac = float((gate & (lfc > 0)).sum() / gate.sum())
    print(f"truth gate up-fraction {truth_up_frac:.4f}", flush=True)
    for tm, tlabel in ((None, "all85"), (noness_t, "noness")):
        rows.append(score(np.full(len(G), 0.9), "all-down", tm) | {"targets": tlabel})
        rows.append(score(P4, "CLEAN4 (uncond, 4 banks)", tm) | {"targets": tlabel})
        Pu = np.nanmean(np.stack([stats[b]["P_uncond"] for b in BANKS]), 0)
        rows.append(score(Pu, "uncond mean, 7 banks", tm) | {"targets": tlabel})
        for scope in ("all", "noness"):
            for Z in ZS:
                for minev in (3, 10):
                    per = []
                    for b in BANKS:
                        dn, up = stats[b][f"dn_{scope}_{Z}"], stats[b][f"up_{scope}_{Z}"]
                        n = dn + up
                        per.append(np.where(n >= minev, (dn + 0.5) / (n + 1.0), np.nan))
                    Pf = np.nanmean(np.stack(per), 0)  # family-balanced: equal weight per bank where defined
                    DN = sum(stats[b][f"dn_{scope}_{Z}"] for b in BANKS)
                    UP = sum(stats[b][f"up_{scope}_{Z}"] for b in BANKS)
                    Pp = np.where(DN + UP >= minev, (DN + 0.5) / (DN + UP + 1.0), np.nan)
                    for P, lab in ((Pf, "COND family-bal"), (Pp, "COND pooled")):
                        rows.append(score(P, f"{lab} perts={scope} z{Z:g} minev{minev}", tm) | {"targets": tlabel})
                        if Z == 4.0 and minev == 3:
                            rows.append(score(P, f"{lab} perts={scope} z4 minev3 conf>=0.1", tm, 0.1) | {"targets": tlabel})
                    # fallback: conditional where defined, CLEAN4 elsewhere (full coverage comparison)
                    if Z == 4.0 and minev == 3:
                        Pc = np.where(np.isfinite(Pf), Pf, P4)
                        rows.append(score(Pc, f"COND fam-bal z4 -> CLEAN4 fallback perts={scope}", tm) | {"targets": tlabel})
        # head-to-head on the SAME genes: where both CLEAN4 and COND(fam-bal, all, z4, minev3) are defined
        per = []
        for b in BANKS:
            n = stats[b]["dn_all_4.0"] + stats[b]["up_all_4.0"]
            per.append(np.where(n >= 3, (stats[b]["dn_all_4.0"] + 0.5) / (n + 1.0), np.nan))
        Pf = np.nanmean(np.stack(per), 0)
        both = np.isfinite(Pf) & np.isfinite(P4)
        rows.append(score(np.where(both, P4, np.nan), "H2H CLEAN4 on shared genes", tm) | {"targets": tlabel})
        rows.append(score(np.where(both, Pf, np.nan), "H2H COND fam-bal z4 on shared genes", tm) | {"targets": tlabel})
        agree = both & ((Pf > 0.5) == (P4 > 0.5))
        rows.append(score(np.where(both & ~agree, Pf, np.nan), "H2H COND where it DISAGREES with CLEAN4", tm) | {"targets": tlabel})
    print(f"\n{'targets':7s} {'label':58s} {'cov':>6s} {'n':>7s} {'acc':>6s} {'acc_w':>6s} {'up':>6s} {'dn':>6s} {'%dn':>5s}")
    for r in rows:
        print(f"{r['targets']:7s} {r['label'][:58]:58s} {r['coverage']:6.3f} {r['n']:7d} {r['acc']:6.4f} {r['acc_w']:6.4f} "
              f"{r['acc_on_truth_up']:6.4f} {r['acc_on_truth_dn']:6.4f} {r['frac_pred_dn']:5.2f}")
    (WORK / "cond634.json").write_text(json.dumps(dict(rows=rows, truth_up_frac=truth_up_frac,
                                                        banks={b: stats[b]["meta"] for b in BANKS}), indent=1))


if __name__ == "__main__":
    sys.exit(main())
