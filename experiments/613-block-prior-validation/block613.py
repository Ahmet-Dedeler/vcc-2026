#!/usr/bin/env python3
"""613 — two PRE-REGISTERED readings for 639's block-prior arm, handed to [21dd4c] by [d4cd5f] (2026-10-02 ~13:35Z).
The rules are fixed in experiments/639-identity-axis-lolo/PREREGISTRATION-BLOCK-DRAFT.md (commit 187ee529), and this
file was committed before its first run.

(A) 458d RETRODICTION of the restricted-LOLO instrument.
    Source: 458d's SHIPPED block prior, Arc's H1 programme (experiments/458-arc-lab-programme/work/arch1d_programme_abc.npz,
    pdown_A, identical for every context; sign DOWN where P_down > 0.5). It is compared against CLEAN4 on block genes.
    Rule: PASS iff (458 − CLEAN4) < 0 on BOTH the gate-gene and the null-gene block readings on >= 2 of the 3 cancer
    folds (hct116 / hek293t / k562). The 458d board read was HARM: fid −0.0055 raw, reach −0.0267 raw.
    If it FAILS, the prereg's WEAK bin becomes "October stays CLEAN4", and only TRANSFERS ships.
(B) NULL-BLOCK (fid filler) reading.
    Source: identity in v3's side rule, sign(z(P_C4 − 1/2) + 0.5·z(I)), on CLEAN4-covered genes only. I is computed with
    Flex-only references, as decomp639 builds it. It is compared against CLEAN4 on NULL block genes.
    Rule: if identity's null-block precision is below CLEAN4's on >= 2 of 3 cancer folds, the arm ships only on TRANSFERS.
    Beside it, NOT used by any rule: the H1 control-halves version (635's construction). Identity comes from control
    half A, the truth from half B, so no control noise is shared.

Definitions, fixed before any number exists:
- surfaces = 636's Inst (inst636.py); CLEAN4 = clean4_lolo (banks forbidden on a surface dropped);
- block = decomp639.callset_for's block (head excluded): lib521's real call set on h1, ds612c.callset(cons0, cnt, 400,
  2200) on the folds;
- gate genes = Inst's gate (q < 0.05, finite lfc != 0, own target excluded);
- null genes = Inst's tested (finite lfc, own target excluded) with q >= 0.05 (finite) and lfc != 0;
- prediction DOWN where z > 0 (CLEAN4: P > 1/2; identity: the combined score > 0; 458: P_down > 1/2), as decomp639;
- paired cells only: both the source and CLEAN4 carry a sign;
- PRIMARY statistic = the pooled paired accuracy difference (source − CLEAN4) over the masked cells, which is 639d's
  "gain". Reported beside it: the per-target mean difference (targets with >= 10 paired cells) and a 95 % target-
  bootstrap CI of the pooled difference (4,000 draws, seed 613);
- KA: identity λ 0.5's gate-block gains must reproduce decomp639.json (h1 / hct116 / hek293t / k562) to 1e-12, or the
  run aborts.
Analysis only; peak < 4 GB (the H1 halves part streams 323's sorted control matrix in gene chunks, as 635).
Output: work/block613.out, work/block613.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
(HERE / "work").mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/639-identity-axis-lolo"))
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
sys.path.insert(0, str(ROOT / "experiments/612-key-gated-hook"))
import decomp639 as DC  # noqa: E402
from inst636 import G480, SURF, Inst, clean4_lolo  # noqa: E402
import lib521 as L  # noqa: E402

CANCER = ("hct116", "hek293t", "k562")
A458 = ROOT / "experiments/458-arc-lab-programme/work/arch1d_programme_abc.npz"
LAM = 0.5
EPS = 1e-9
LINES = []


def say(s=""):
    print(s, flush=True)
    LINES.append(s)


def paired(src_dn, c4_dn, defined, truth_dn, mask, rng):
    """Pooled paired accuracy difference (source - CLEAN4) over mask & defined, plus per-target mean and target bootstrap."""
    m = mask & defined[None, :]
    n = m.sum(1).astype(float)
    s = ((src_dn[None, :] == truth_dn) & m).sum(1).astype(float)
    c = ((c4_dn[None, :] == truth_dn) & m).sum(1).astype(float)
    N = n.sum()
    pooled = float((s.sum() - c.sum()) / max(N, 1))
    ok = n >= 10
    per_t = float(np.mean((s[ok] - c[ok]) / n[ok])) if ok.any() else float("nan")
    idx = rng.integers(0, len(n), size=(4000, len(n)))
    bs = (s[idx].sum(1) - c[idx].sum(1)) / np.maximum(n[idx].sum(1), 1)
    return dict(n=int(N), acc_src=float(s.sum() / max(N, 1)), acc_c4=float(c.sum() / max(N, 1)), pooled=pooled,
                per_target=per_t, lo=float(np.quantile(bs, 0.025)), hi=float(np.quantile(bs, 0.975)))


def main() -> None:
    I = Inst()
    prof = DC.flex_profiles()
    z458 = np.load(A458, allow_pickle=True)
    g458, p458 = z458["genes"].astype(str), np.asarray(z458["pdown_A"], float)
    assert np.array_equal(np.asarray(z458["pdown_A"], float), np.asarray(z458["pdown_B"], float), equal_nan=True)
    ka = json.loads((ROOT / "experiments/639-identity-axis-lolo/work/decomp639.json").read_text())
    rng = np.random.default_rng(613)
    res = {}
    for s in ("h1",) + CANCER:
        D = I.S[s]
        T, G, lfc, gate, tested = D["T"], D["G"], D["lfc"], D["gate"], D["tested"]
        q = np.load(G480 / SURF[s], allow_pickle=True)["q"].astype(float)
        null = tested & np.isfinite(q) & (q >= 0.05) & (lfc != 0)
        assert not (null & gate).any()
        if s == "h1":
            own, refs = prof["h1"], [prof[k] for k in ("ctxA", "ctxB", "ctxC")]
        else:
            zz = np.load(G480 / SURF[s], allow_pickle=True)
            own = dict(zip(zz["genes"].astype(str), zz["ctrl_cpm"].astype(float)))
            refs = [prof[k] for k in ("ctxA", "ctxB", "ctxC", "h1")]
        Iv = DC.identity(own, refs, list(G))
        zi = (Iv - np.nanmean(Iv)) / np.nanstd(Iv)
        C4, banks = clean4_lolo(I, s)
        zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
        sc = zc + LAM * zi
        head, block = DC.callset_for(I, s)
        dn = lfc < 0
        both_id = np.isfinite(zc) & np.isfinite(zi)
        P458 = I.gene_vec(s, g458, p458)
        both_458 = np.isfinite(zc) & np.isfinite(P458)
        r = {"banks": banks, "n_targets": int(len(T))}
        # KA: identity's gate-block gain == decomp639 (639d)
        g = paired(sc > 0, zc > 0, both_id, dn, gate & block, rng)
        want = ka[s]["regions"]["lam0.5"]["block"]["gain"]
        if abs(g["pooled"] - want) > 1e-12:
            raise SystemExit(f"KA FAIL {s}: identity gate-block gain {g['pooled']!r} vs decomp639 {want!r}")
        r["identity_gate"] = g
        r["identity_null"] = paired(sc > 0, zc > 0, both_id, dn, null & block, rng)
        r["h1prog458_gate"] = paired(P458 > 0.5, zc > 0, both_458, dn, gate & block, rng)
        r["h1prog458_null"] = paired(P458 > 0.5, zc > 0, both_458, dn, null & block, rng)
        # descriptive: 458's extension (CLEAN4 uncovered) accuracy on block cells
        ext = ~np.isfinite(zc) & np.isfinite(P458)
        for lab, M in (("gate", gate & block), ("null", null & block)):
            m = M & ext[None, :]
            r[f"h1prog458_ext_{lab}"] = dict(n=int(m.sum()), acc=float((((P458 > 0.5)[None, :] == dn) & m).sum() / max(m.sum(), 1)))
        res[s] = r
    say("KA: identity λ 0.5 gate-block gains reproduce decomp639.json on h1 / hct116 / hek293t / k562: PASS")
    say("\nPooled paired accuracy difference vs CLEAN4 on BLOCK genes (+ = the source is right more often), [95% target "
        "bootstrap], per-target mean, n paired cells")
    for s, r in res.items():
        say(f"== {s} (CLEAN4 banks: {', '.join(r['banks'])})")
        for k in ("identity_gate", "identity_null", "h1prog458_gate", "h1prog458_null"):
            v = r[k]
            say(f"  {k:16s} {v['pooled']:+.4f} [{v['lo']:+.4f},{v['hi']:+.4f}]  per-target {v['per_target']:+.4f}  "
                f"acc {v['acc_src']:.4f} vs CLEAN4 {v['acc_c4']:.4f}  n {v['n']}")
        say(f"  (458 on CLEAN4-uncovered block cells, descriptive: gate acc {r['h1prog458_ext_gate']['acc']:.4f} n "
            f"{r['h1prog458_ext_gate']['n']}; null acc {r['h1prog458_ext_null']['acc']:.4f} n {r['h1prog458_ext_null']['n']})")
    # (A) the 458d retrodiction verdict
    neg_both = {s: (res[s]["h1prog458_gate"]["pooled"] < 0) and (res[s]["h1prog458_null"]["pooled"] < 0) for s in CANCER}
    passA = sum(neg_both.values()) >= 2
    say(f"\n(A) 458d RETRODICTION: 458 − CLEAN4 negative on BOTH readings: " + ", ".join(f"{s} {v}" for s, v in neg_both.items())
        + f" -> {'PASS' if passA else 'FAIL'} (needs >= 2 of 3)")
    # (B) the null-block reading
    below = {s: res[s]["identity_null"]["pooled"] < 0 for s in CANCER}
    nB = sum(below.values())
    say(f"(B) NULL-BLOCK: identity below CLEAN4 on: " + ", ".join(f"{s} {v}" for s, v in below.items())
        + f" -> {nB} of 3 -> " + ("the arm ships ONLY on TRANSFERS" if nB >= 2 else "no restriction from the null reading"))
    out = {"A_pass": passA, "A_neg_both": neg_both, "B_below": below, "B_restricts": nB >= 2, "surfaces": res}
    # (B') H1 control halves (635's construction), reported beside; not used by any rule
    out["h1_halves"] = h1_halves(I, prof)
    (HERE / "work/block613.json").write_text(json.dumps(out, indent=1, default=float))
    (HERE / "work/block613.out").write_text("\n".join(LINES) + "\n")
    say("wrote work/block613.out, work/block613.json")


def h1_halves(I, prof):
    """635's control halves on h1-test: identity from control half A, null truth vs half B (no shared control noise),
    against the full-pool version (identity and truth from the same pool). Block = lib521's real block."""
    z = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
    keep = np.asarray(z["keep"], bool)
    G = z["genes"].astype(str)[keep]
    T = z["targets"].astype(str)
    q, lfc, ref_lin = z["q"], z["lfc"], np.asarray(z["ref_lin"], float)
    cs = z["cs"]
    assert cs.shape[0] == len(G)
    rng = np.random.default_rng(635)          # 635's seed and per-gene split: the same halves as 635
    h = cs.shape[1] // 2
    mA, mB = np.empty(cs.shape[0]), np.empty(cs.shape[0])
    for i in range(0, cs.shape[0], 1000):
        X = np.expm1(cs[i:i + 1000])
        idx = np.argsort(rng.random(X.shape), axis=1)
        Xs = np.take_along_axis(X, idx, axis=1)
        mA[i:i + 1000], mB[i:i + 1000] = Xs[:, :h].mean(1), Xs[:, h:].mean(1)
        del X, idx, Xs
    del cs
    ka = float(np.nanmax(np.abs((mA + mB) / 2 - ref_lin) / np.maximum(ref_lin, 1e-12)))
    mp = np.where(np.isfinite(lfc), ref_lin[None, :] * np.exp2(lfc), np.nan)
    LB = np.log2((mp + EPS) / (mB[None, :] + EPS))
    own = G[None, :] == T[:, None]
    st = L.load_state(0)
    assert list(map(str, st["genes"])) == list(G) and list(map(str, st["targets"])) == list(T)
    block = st["block"]
    gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & ~own
    refs = [prof[k] for k in ("ctxA", "ctxB", "ctxC")]
    zfull = (lambda v: (v - np.nanmean(v)) / np.nanstd(v))
    I_A = zfull(DC.identity(dict(zip(G, mA)), refs, list(G)))
    I_F = zfull(DC.identity(dict(zip(G, ref_lin)), refs, list(G)))
    C4h, _ = clean4_lolo(I, "h1")
    C4 = np.full(len(G), np.nan)
    pos = {g: i for i, g in enumerate(I.S["h1"]["G"])}
    col = np.array([pos.get(g, -1) for g in G])
    C4[col >= 0] = C4h[col[col >= 0]]
    zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
    rng2 = np.random.default_rng(613)
    out = {"ka_halves_mean_rel_err": ka}
    for tname, TL in (("truth_full", lfc), ("truth_B", LB)):
        nullm = block & ~gate & ~own & np.isfinite(TL) & (TL != 0)
        tdn = TL < 0
        for iname, Iz in (("identity_full", I_F), ("identity_A", I_A)):
            sc = zc + LAM * Iz
            out[f"{iname}|{tname}"] = paired(sc > 0, zc > 0, np.isfinite(zc) & np.isfinite(Iz), tdn, nullm, rng2)
    say(f"\n(B') H1 control halves (635), null block genes, identity vs CLEAN4 (pooled; not used by a rule); halves KA "
        f"max rel err {ka:.1e}")
    for k, v in out.items():
        if k == "ka_halves_mean_rel_err":
            continue
        say(f"  {k:28s} {v['pooled']:+.4f} [{v['lo']:+.4f},{v['hi']:+.4f}]  acc {v['acc_src']:.4f} vs CLEAN4 {v['acc_c4']:.4f}  n {v['n']}")
    return out


if __name__ == "__main__":
    main()
