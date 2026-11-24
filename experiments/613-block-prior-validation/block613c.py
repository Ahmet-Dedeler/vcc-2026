#!/usr/bin/env python3
"""613c: which block SIGN SOURCE maximises null-block precision (fid's filler) on held-out cancer lines? A screen on 613's
restricted-LOLO instrument (validated on 458d), committed before its first run.

The question. Our board fid (0.5144) implies a null-call sign precision of ~0.507-0.510 on A/B/C. The pack sits at fid
0.528-0.530, which needs ~0.52 (610). The block (~2,200 of 2,600 calls, ~95 % null genes) sets that precision, and today
it carries CLEAN4's sign (clean4q80). Does any public sign source beat CLEAN4 on NULL block genes of held-out lines?

Retrodiction rows, fixed before any number exists. The screen is read only if both hold:
  R1 (458d, board HARM): 458's H1 programme − CLEAN4 < 0 on the null block, on >= 2 of 3 cancer folds (613 already shows
     3/3; re-checked here).
  R2 (392, board fid +0.0082 raw when CLEAN4 replaced the design's own block signs): CLEAN4 − the KEY's own sign > 0 on
     the null block, on >= 2 of 3 cancer folds. The key is the 3-donor consensus on the folds and the design key on h1.
Candidates, all compared with CLEAN4 on paired cells (both signed), null block AND gate block:
  all-down; key sign; each CLEAN4 bank alone (where allowed on the surface); Replogle RPE1 (a public line outside CLEAN4,
  371's unconditional P_down); identity alone (z(I) > 0 = down, Flex refs).
Primary: the pooled paired Δ precision vs CLEAN4 on the null block (613's statistic), on the 3 cancer folds; h1 reported
beside. A candidate is a LEAD only if it is > +0.005 on all 3 cancer folds AND not below CLEAN4 on the gate block on any
cancer fold by more than 0.01. Try count: 7 candidates x 4 surfaces, descriptive; any lead needs its own prereg.
Output: work/block613c.out.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import block613 as B  # noqa: E402
import block613b as BB  # noqa: E402
import decomp639 as DC  # noqa: E402
from inst636 import BANK_FORBIDDEN, G480, ROOT, SURF, Inst, bank_uncond, clean4_lolo  # noqa: E402

BANKS = ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")


def main():
    I = Inst()
    prof = DC.flex_profiles()
    cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names, cg = list(cp["bank_names"].astype(str)), cp["genes"].astype(str)
    rg, rp = bank_uncond("rpe1")
    z458 = np.load(B.A458, allow_pickle=True)
    rng = np.random.default_rng(6131)
    lines, R = [], {}
    for s in ("h1",) + B.CANCER:
        D = I.S[s]
        T, G, lfc, gate, tested = D["T"], D["G"], D["lfc"], D["gate"], D["tested"]
        q = np.load(G480 / SURF[s], allow_pickle=True)["q"].astype(float)
        null = tested & np.isfinite(q) & (q >= 0.05) & (lfc != 0)
        _, block = DC.callset_for(I, s)
        C4, _ = clean4_lolo(I, s)
        c4dn, c4def = C4 > 0.5, np.isfinite(C4)
        dn = lfc < 0
        if s == "h1":
            own, refs = prof["h1"], [prof[k] for k in ("ctxA", "ctxB", "ctxC")]
        else:
            zz = np.load(G480 / SURF[s], allow_pickle=True)
            own = dict(zip(zz["genes"].astype(str), zz["ctrl_cpm"].astype(float)))
            refs = [prof[k] for k in ("ctxA", "ctxB", "ctxC", "h1")]
        Iv = DC.identity(own, refs, list(G))
        E = BB.e0_for(I, s)          # |key|; its SIGN is needed too
        key_signed = signed_key(I, s)
        cands = {"458 programme (R1)": (I.gene_vec(s, z458["genes"].astype(str), np.asarray(z458["pdown_A"], float)) > 0.5,
                                        np.isfinite(I.gene_vec(s, z458["genes"].astype(str), np.asarray(z458["pdown_A"], float))), False),
                 "key sign (R2: CLEAN4 - this)": (key_signed < 0, key_signed != 0, True),
                 "all-down": (np.ones(len(G), bool), np.ones(len(G), bool), False),
                 "identity alone": (Iv > np.nanmean(Iv), np.isfinite(Iv), False)}
        for b in BANKS:
            if s in BANK_FORBIDDEN.get(b, set()):
                continue
            P = I.gene_vec(s, cg, cp["P_banks"].astype(float)[names.index(b)])
            cands[f"bank {b}"] = (P > 0.5, np.isfinite(P) & (P != 0.5), False)
        if rg is not None:
            P = I.gene_vec(s, np.array(rg), rp)
            cands["rpe1 (outside CLEAN4)"] = (P > 0.5, np.isfinite(P) & (P != 0.5), False)
        R[s] = {}
        lines.append(f"== {s}: CLEAN4 null-block precision {float((c4dn[None, :] == dn)[null & block & c4def[None, :]].mean()):.4f}")
        for lab, (src, sdef, per_target) in cands.items():
            if per_target:      # the key is per target: a (T, G) prediction
                res = {}
                for rl, M in (("null", null & block), ("gate", gate & block)):
                    m = M & sdef & c4def[None, :]
                    s_ok = ((src == dn) & m).sum(1).astype(float)
                    c_ok = ((c4dn[None, :] == dn) & m).sum(1).astype(float)
                    n = m.sum(1).astype(float)
                    res[rl] = dict(pooled=float((s_ok.sum() - c_ok.sum()) / max(n.sum(), 1)), n=int(n.sum()))
                # R2 is CLEAN4 - key: flip the sign of the difference
                R[s][lab] = {k: dict(v, pooled=-v["pooled"]) for k, v in res.items()}
            else:
                R[s][lab] = {rl: B.paired(src, c4dn, sdef & c4def, dn, M, rng)
                             for rl, M in (("null", null & block), ("gate", gate & block))}
            v = R[s][lab]
            lines.append(f"  {lab:30s} null {v['null']['pooled']:+.4f} (n {v['null']['n']})   gate {v['gate']['pooled']:+.4f} "
                         f"(n {v['gate']['n']})")
    r1 = sum(R[s]["458 programme (R1)"]["null"]["pooled"] < 0 for s in B.CANCER) >= 2
    r2 = sum(R[s]["key sign (R2: CLEAN4 - this)"]["null"]["pooled"] > 0 for s in B.CANCER) >= 2
    lines.append(f"\nRetrodiction R1 (458d): {'PASS' if r1 else 'FAIL'}; R2 (392): {'PASS' if r2 else 'FAIL'}"
                 + ("" if (r1 and r2) else "  -> the screen is NOT read"))
    if r1 and r2:
        leads = []
        for lab in R["hct116"]:
            if lab.startswith(("458", "key")):
                continue
            ok = all(lab in R[s] and R[s][lab]["null"]["pooled"] > 0.005 and R[s][lab]["gate"]["pooled"] >= -0.01 for s in B.CANCER)
            if ok:
                leads.append(lab)
        lines.append("LEADS (> +0.005 null on all 3 cancer folds, gate not below −0.01): " + (", ".join(leads) if leads else "none"))
    print("\n".join(lines))
    (HERE / "work/block613c.out").write_text("\n".join(lines) + "\n")


def signed_key(I, s):
    """The signed key on the instrument's axes: the design key on h1 (lib521), the 3-donor consensus on the folds."""
    D = I.S[s]
    T, G = D["T"], D["G"]
    if s == "h1":
        st = B.L.load_state(0)
        pos = {g: i for i, g in enumerate(st["genes"])}
        cols = np.array([pos.get(g, -1) for g in G])
        K = np.zeros((len(T), len(G)))
        K[:, cols >= 0] = st["key"][:, cols[cols >= 0]]
        return K
    import ds612 as DS
    Sx = DS.F.load_surface(DC.FOLD[s])
    genes_f, targets_f = [str(g) for g in Sx["genes"]], list(Sx["targets"])
    cons0 = np.nan_to_num(np.nanmean(np.stack([DS.load_donor_rows(d, targets_f, np.array(genes_f), pc_cpm=1.0)[0]
                                               for d in Sx["donors"]]), 0))
    tpos, gpos = {t: i for i, t in enumerate(targets_f)}, {g: i for i, g in enumerate(genes_f)}
    ti, gi = np.array([tpos.get(t, -1) for t in T]), np.array([gpos.get(g, -1) for g in G])
    K = np.zeros((len(T), len(G)))
    for a, t in enumerate(ti):
        if t >= 0:
            K[a, gi >= 0] = cons0[t, gi[gi >= 0]]
    return K


if __name__ == "__main__":
    main()
