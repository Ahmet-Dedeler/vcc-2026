#!/usr/bin/env python3
"""632 — analytic H1 price of confidence-graded signs + magnitudes on UNCALLED genes (the hook's territory), i.e. what a
rank-neutral delivery could earn if it removed the leak cap (426/435), with identity-boosted signs (639).

L1 theory: with sign accuracy a and the gate-conditional |l| distribution F, the expected-L1-optimal magnitude is
m* = F^{-1}((2a - 1) / (2a)) (0 if a <= 0.5). At a ~ 0.55 that is ~the MDE (why 480 found nothing beats it); at a 0.65-0.75
it is the 23rd-33rd percentile of |l|, too large for today's delivery (it leaks calls) but free for a rank-neutral one.

Policy on uncalled genes of 446bj2z's H1 construction (435 zs stats, seeds 0/1/2): score S = z(CLEAN4 - .5) + lam z(I)
(639; identity alone where CLEAN4 is undefined); direction = down if S > 0; magnitude m_bin by |S| decile, with a_bin and F_bin
CROSS-FITTED on the other half of the targets (odd/even split, both ways). Called genes keep the delivered lfc (their sign
and size are tied to the call design). Reported: Δnmae (raw, BRIEF §7 per-target mean, targets with gate >= 10) vs the
delivered lfc; split into UP-direction and DOWN-direction genes (DOWN capacity is bounded in a rank-neutral delivery);
and the same with signs from CLEAN4 alone (lam 0), and with the delivered sign kept (magnitude-only).
Analysis only.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
sys.path.insert(0, str(ROOT / "experiments/639-identity-axis-lolo"))
from inst636 import G480, SURF, Inst, clean4_lolo  # noqa: E402
import ident639 as M  # noqa: E402

WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
NB = 10

ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
keep = np.asarray(ref["keep"], bool)
G = ref["genes"].astype(str)[keep]
T = ref["targets"].astype(str)
q, lfc = ref["q"], ref["lfc"]
gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & (G[None, :] != T[:, None])
absl = np.abs(np.nan_to_num(lfc))
scored = gate.sum(1) >= 10


def nmae(pred):
    num = (np.abs(np.nan_to_num(pred) - np.nan_to_num(lfc)) * gate).sum(1)
    den = (absl * gate).sum(1)
    return float(np.mean((num / np.maximum(den, 1e-12))[scored]))


def identity_h1():
    prof = {}
    for s, f in SURF.items():
        if s == "h1":
            continue
        z = np.load(G480 / f, allow_pickle=True)
        prof[s] = dict(zip(z["genes"].astype(str), z["ctrl_cpm"].astype(float)))
    prof.update(M.abc_profiles())
    gi = {g: i for i, g in enumerate(G)}
    cg = sorted(set(G) & set.intersection(*(set(p) for p in prof.values())))
    cols = np.array([gi[g] for g in cg])
    others = np.mean([np.log1p(np.array([prof[k][g] for g in cg])) for k in prof], 0)
    iv = np.full(len(G), np.nan)
    iv[cols] = np.log1p(np.asarray(ref["ref_lin"], float)[cols]) - others
    return iv


def main():
    I = Inst()
    C4, _ = clean4_lolo(I, "h1")
    zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
    iv = identity_h1()
    zi = (iv - np.nanmean(iv)) / np.nanstd(iv)
    halves = [np.arange(len(T)) % 2 == 0, np.arange(len(T)) % 2 == 1]
    out = {}
    for seed in (0, 1, 2):
        st = np.load(ROOT / f"experiments/435-zs-on-coarse/work/h1/stats-base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz",
                     allow_pickle=True)
        d = np.nan_to_num(np.asarray(st["lfc"], float))
        called = np.asarray(st["q"], float) < 0.05
        unc = ~called & np.isfinite(lfc) & (G[None, :] != T[:, None])
        base = nmae(d)
        res = {"base_nmae": base}
        for lab, lam, keep_sign in (("C4+1.0I signs, graded mag", 1.0, False), ("C4+0.5I signs, graded mag", 0.5, False),
                                    ("CLEAN4 signs, graded mag", 0.0, False), ("delivered signs, graded mag", None, True)):
            if keep_sign:
                S = -np.sign(d)                    # positive -> down, as for the scores
                conf_src = np.where(np.isfinite(zc), np.abs(zc), 0.0)[None, :] * np.ones_like(d)
            else:
                Sg = np.where(np.isfinite(zc), zc + lam * np.nan_to_num(zi), zi if lam > 0 else np.nan)
                S = np.broadcast_to(Sg[None, :], d.shape)
                conf_src = np.abs(S)
            dirn = np.where(S > 0, -1.0, np.where(S < 0, 1.0, 0.0))   # predicted lfc sign
            pred = d.copy()
            parts = {"up": d.copy(), "down": d.copy()}
            for hi in range(2):
                fit, ev = halves[hi], halves[1 - hi]
                cf = conf_src[fit][unc[fit] & gate[fit] & np.isfinite(conf_src[fit])]
                edges = np.quantile(cf, np.linspace(0, 1, NB + 1)) if cf.size else np.array([0, np.inf])
                for b in range(NB):
                    lo, hb = edges[b], edges[b + 1] if b < NB - 1 else np.inf
                    inb = (conf_src >= lo) & (conf_src < hb) & unc & np.isfinite(conf_src) & (dirn != 0)
                    mf = inb & fit[:, None] & gate
                    if mf.sum() < 50:
                        continue
                    a = float(((np.sign(np.nan_to_num(lfc)) == dirn) & mf).sum() / mf.sum())
                    if a <= 0.5:
                        m = 0.0
                    else:
                        m = float(np.quantile(absl[mf], (2 * a - 1) / (2 * a)))
                    me = inb & ev[:, None]
                    pred[me] = dirn[me] * m
                    for k, sel in (("up", dirn > 0), ("down", dirn < 0)):
                        mm = me & sel
                        parts[k][mm] = dirn[mm] * m
            res[lab] = dict(d_nmae=nmae(pred) - base, d_nmae_up_only=nmae(parts["up"]) - base,
                            d_nmae_down_only=nmae(parts["down"]) - base)
        out[f"seed{seed}"] = res
        print(f"seed {seed}: base nmae {base:.5f}", flush=True)
        for k, v in res.items():
            if k != "base_nmae":
                print(f"   {k:30s} dnmae {v['d_nmae']:+.5f}  (UP-direction genes only {v['d_nmae_up_only']:+.5f}, "
                      f"DOWN only {v['d_nmae_down_only']:+.5f})", flush=True)
    (WORK / "grade632.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    sys.exit(main())
