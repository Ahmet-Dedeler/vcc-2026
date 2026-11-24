#!/usr/bin/env python3
"""623 — does a target's IDENTITY-EROSION strength transfer across lines? (council 600, [7772dc], exploratory)

622/639: under perturbation, genes a line expresses above other lines tend to go DOWN (the identity axis I). A gene-level
transfer from donors cannot carry this to a recipient: the recipient's identity genes are not the donors'. But if HOW
STRONGLY a knockdown erodes identity is a property of the TARGET (conserved across lines), then
    delta_hat_L(t) = shared(t) + kappa * beta_hat(t) * I_L_unit,
with beta_hat(t) = mean over donor lines d of the projection of delta_d(t) on d's own identity axis, adds a recipient-
specific bulk component that no donor row contains.

Lines with real control profiles and per-target log-pseudobulk effects (TS 5e4): K562 (k562_gwps), HCT116, HEK293T
(Orion), H1 (Arc, h1_2025_full), RPE1, Jurkat, HepG2 (Weissman essential screens). Identity axes from the banks' b_ctrl
renormalised over the common genes, references = every other line + A/B/C Flex controls (622's lesson: Flex references).
Tests, per held-out line L with donors = the other lines covering each target:
  (1) conservation: Spearman across shared targets of the target-centred loading beta_c(t) between line pairs;
  (2) prediction: mean per-target cosine of delta_hat with the truth, baseline (donor mean) vs + kappa * beta_hat * I_L,
      kappa in {0.5, 1, 2}, target-centred and raw loadings; energy-weighted cosine (X) too.
Same-study pairs (HCT116<->HEK293T; the Weissman lines) are flagged; K562-held-out and H1-held-out are cross-study.
Analysis only, < 4 GB.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
TS = 5e4
BANKS = {"K562": ("data/donors/k562_gwps", "weissman"), "HCT116": ("data/donors/orion_hct116", "xaira"),
         "HEK293T": ("data/donors/orion_hek293t", "xaira"), "H1": ("data/donors/h1_2025_full", "arc"),
         "RPE1": ("data/donors/rpe1", "weissman"), "Jurkat": ("data/donors/nadig_jurkat", "weissman"),
         "HepG2": ("experiments/166-donor-acquisition/donors/nadig_hepg2", "weissman")}


def rl(p: Path) -> list[str]:
    return [x.strip() for x in p.read_text().split("\n") if x.strip()]


def main() -> None:
    meta = {}
    for L, (d, fam) in BANKS.items():
        d = ROOT / d
        meta[L] = dict(genes=rl(d / "genes.txt"), perts=rl(d / "perts.txt"), b=np.load(d / "b_ctrl.npy").astype(float),
                       D=np.load(d / "delta.npy", mmap_mode="r"), fam=fam)
    abc = np.load(ROOT / "experiments/622-identity-track/work/abc_cpm622.npz", allow_pickle=True)
    abc_genes = abc["genes"].astype(str)
    common = set(abc_genes)
    for L in meta:
        common &= set(meta[L]["genes"])
    common = sorted(common)
    G = len(common)
    print(f"common genes {G}", flush=True)
    # control log-profiles renormalised over the common genes
    prof = {}
    for L, m in meta.items():
        pos = {g: i for i, g in enumerate(m["genes"])}
        v = np.expm1(m["b"][[pos[g] for g in common]])
        prof[L] = np.log1p(TS * v / v.sum())
    apos = {g: i for i, g in enumerate(abc_genes)}
    for c in ("ctxA", "ctxB", "ctxC"):
        v = abc[c].astype(float)[[apos[g] for g in common]]
        prof[c] = np.log1p(TS * v / v.sum())
    ident = {}
    for L in meta:
        refs = [k for k in prof if k != L]
        I = prof[L] - np.mean([prof[k] for k in refs], 0)
        I = I - I.mean()
        ident[L] = I / np.linalg.norm(I)
    # per-target effect rows on the common genes, for the union of targets that >= 2 lines carry
    from collections import Counter
    cnt = Counter(t for m in meta.values() for t in set(m["perts"]))
    targets = sorted(t for t, k in cnt.items() if k >= 2)
    rows = {}
    for L, m in meta.items():
        gpos = {g: i for i, g in enumerate(m["genes"])}
        cols = np.array([gpos[g] for g in common])
        ppos = {p: i for i, p in enumerate(m["perts"])}
        ts = [t for t in targets if t in ppos]
        ri = np.array([ppos[t] for t in ts])
        order = np.argsort(ri)
        # chunked: never materialise full-width rows (the first two runs peaked at 12.5 / 14.8 GB doing that)
        R = np.concatenate([np.asarray(m["D"][ri[order][a:a + 1000]], dtype=np.float32)[:, cols]
                            for a in range(0, len(ri), 1000)]) if len(ri) else np.zeros((0, len(cols)), np.float32)
        rows[L] = dict(zip([ts[k] for k in order], R))
        print(f"  {L}: {len(ts)} targets carried", flush=True)
    # loadings: beta[L][t] = <delta_L(t), I_L>; target-centred within line
    beta = {L: {t: float(v @ ident[L]) for t, v in rows[L].items()} for L in meta}
    betac = {L: {t: b - np.mean(list(beta[L].values())) for t, b in beta[L].items()} for L in meta}
    out = {"common_genes": G, "conservation": {}, "prediction": {}}
    from scipy.stats import spearmanr
    Ls = list(meta)
    for i, a in enumerate(Ls):
        for b2 in Ls[i + 1:]:
            sh = sorted(set(betac[a]) & set(betac[b2]))
            if len(sh) < 30:
                continue
            r = spearmanr([betac[a][t] for t in sh], [betac[b2][t] for t in sh]).correlation
            same = meta[a]["fam"] == meta[b2]["fam"]
            out["conservation"][f"{a}~{b2}"] = dict(n=len(sh), spearman=float(r), same_study=bool(same))
            print(f"[1] {a:8s}~{b2:8s} n {len(sh):5d}  Spearman(beta_c) {r:+.3f}{'  (same study)' if same else ''}", flush=True)
    for L in Ls:
        donors = [d for d in Ls if d != L]
        ts = [t for t in rows[L] if any(t in rows[d] for d in donors)]
        if len(ts) < 30:
            continue
        res = {"n_targets": len(ts)}
        cos_base, cos_add, e_base, e_add = [], {}, [], {}
        for t in ts:
            ds = [d for d in donors if t in rows[d]]
            base = np.mean([rows[d][t] for d in ds], 0)
            bhat = np.mean([betac[d][t] for d in ds]) + np.mean(list(beta[L].values())) * 0  # centred loading
            truth = rows[L][t]
            nt = np.linalg.norm(truth)
            if nt == 0 or np.linalg.norm(base) == 0:
                continue
            cos_base.append(float(base @ truth / (np.linalg.norm(base) * nt))); e_base.append((base, truth))
            for kappa in (0.5, 1.0, 2.0):
                p = base + kappa * bhat * ident[L]
                cos_add.setdefault(kappa, []).append(float(p @ truth / (np.linalg.norm(p) * nt)))
                e_add.setdefault(kappa, []).append((p, truth))

        def X(pairs):  # energy-weighted cosine: sum <p, y> / sqrt(sum |p|^2 sum |y|^2) over targets (per-target scaled)
            num = sum(float(p @ y) / np.linalg.norm(p) * np.linalg.norm(y) for p, y in pairs)
            den = sum(np.linalg.norm(y) ** 2 for _, y in pairs)
            return num / den
        res["cos_base"] = float(np.mean(cos_base)); res["X_base"] = X(e_base)
        for kappa in cos_add:
            d = np.array(cos_add[kappa]) - np.array(cos_base)
            idx = np.arange(d.size); rng = np.random.default_rng(6230)
            bs = [d[rng.choice(idx, idx.size)].mean() for _ in range(2000)]
            res[f"k{kappa}"] = dict(dcos=float(d.mean()), ci=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))],
                                    dX=float(X(e_add[kappa]) - res["X_base"]))
        out["prediction"][L] = res
        print(f"[2] held-out {L:8s} n {res['n_targets']:5d}  cos_base {res['cos_base']:.4f}  X_base {res['X_base']:.4f} | " +
              "  ".join(f"κ{k}: Δcos {res[f'k{k}']['dcos']:+.4f} [{res[f'k{k}']['ci'][0]:+.4f},{res[f'k{k}']['ci'][1]:+.4f}] ΔX {res[f'k{k}']['dX']:+.4f}"
                        for k in (0.5, 1.0, 2.0)), flush=True)
    (WORK / "erode623.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
