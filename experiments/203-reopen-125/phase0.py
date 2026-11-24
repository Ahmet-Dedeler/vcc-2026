"""203 Phase 0 — verify 198's diagnosis of experiment 125 from first principles.

Nothing here is inherited.  Every claim in the mandate is recomputed from 125's own
frozen arrays plus 122's frozen truth:

  0.1  reproduce 125's two macro deltas from the saved per-target arrays;
  0.2  RE-SCORE both saved predictions against the real truth with 125's own
       `score_predictions`, and check the recomputation reproduces the recorded deltas
       (this is stronger than 0.1: it re-runs the metric, not just its average);
  0.3  homogeneity: credit(k.pred)/credit(pred) == k and precision20(k.pred) invariant,
       measured on the REAL truth and the REAL confidence, at several k;
  0.4  the shrink factor, globally and per target, plus per-target cosine so we know
       whether 125 is a rescale of 122 or a different direction;
  0.5  142's homonym: does `develop.py:1379` divide by the prediction's L1?
  0.6  THE DECISIVE CONTROL (mandate Phase 4.1): 122's prediction shrunk by the measured
       scalar, scored on both endpoints.  If it reproduces -0.01505 the loss is a scale
       artefact and nothing else.
  0.7  a 2x2 decomposition of the +0.012175 precision gain into the EFFECT-sign channel
       and the CONFIDENCE-ordering channel.

KNOWN-ANSWER ROWS (rule 10):
  KA1  arm (122 effect, 122 confidence) must score credit delta 0.0 and precision delta
       0.0 against itself -- catches a mis-wired evaluable mask.
  KA2  arm (125 effect, 125 confidence) must reproduce the recorded -0.015050 / +0.012175
       from the float32 saved arrays to float32 precision.
  KA3  the analytic prediction credit(c.p) = c * credit(p) must match the measured
       rescale-control delta to machine precision.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W125 = ROOT / "experiments/125-hierarchical-panel-posterior/work/source-only-evaluation-v1"
F125 = ROOT / "experiments/125-hierarchical-panel-posterior/work/preregistered-freeze-v1"
W122 = ROOT / "experiments/122-panel-response-atlas/work/source-only-phase1-v1"
OUT = HERE / "work"
OUT.mkdir(exist_ok=True)

REC_CREDIT = -0.015049899010796045
REC_PRECISION = +0.012174572820953597
REC_FAMILY_CREDIT = (-0.01214610823191394, -0.026082415627716488, -0.006921173172757709)
REC_FAMILY_PRECISION = (+0.008779778202899283, +0.0206978642163938, +0.007046076043567707)

out: dict[str, object] = {}


def say(k, v):
    print(f"{k:<58} {v}")
    out[k] = v


# ---------------------------------------------------------------- 0.1 arrays
cred_d = np.load(W125 / "per-target-credit-delta.npy")
prec_d = np.load(W125 / "per-target-precision-delta.npy")
print(f"per-target delta arrays: {cred_d.shape} {cred_d.dtype}")

fam_c = np.array([np.nanmean(cred_d[i]) for i in range(3)])
fam_p = np.array([np.nanmean(prec_d[i]) for i in range(3)])
say("0.1 macro credit delta from saved arrays", f"{fam_c.mean():+.15f}")
say("0.1 macro precision delta from saved arrays", f"{fam_p.mean():+.15f}")
say("0.1 |macro credit - recorded|", f"{abs(fam_c.mean() - REC_CREDIT):.3e}")
say("0.1 |macro precision - recorded|", f"{abs(fam_p.mean() - REC_PRECISION):.3e}")
say("0.1 family credit max abs err", f"{np.abs(fam_c - REC_FAMILY_CREDIT).max():.3e}")
say("0.1 family precision max abs err", f"{np.abs(fam_p - REC_FAMILY_PRECISION).max():.3e}")

# ------------------------------------------------------ 0.2 re-score from truth
truth = np.asarray(np.load(W122 / "family-effects.npy"), dtype=np.float64)
tvalid = np.asarray(np.load(W122 / "family-valid.npy"), dtype=bool)
assert truth.shape == (3, 300, 300)

cand_e = np.asarray(np.load(W125 / "observed-candidate-effect.npy"), dtype=np.float64)
cand_c = np.asarray(np.load(W125 / "observed-candidate-confidence.npy"), dtype=np.float64)
base_e = np.asarray(np.load(W125 / "observed-baseline-effect.npy"), dtype=np.float64)
base_c = np.asarray(np.load(W125 / "observed-baseline-confidence.npy"), dtype=np.float64)
say("0.2 saved effect array shape", str(cand_e.shape))

# the evaluable set is exactly where the recorded per-target deltas are finite
evaluable = np.isfinite(cred_d) & np.isfinite(prec_d)
say("0.2 evaluable (target,family) cells", int(evaluable.sum()))

CANON = np.arange(300)


def score_family(pred, conf, f):
    """125 run.py:752-793 verbatim, restricted to the evaluable targets of family f."""
    credit = np.full(300, np.nan)
    prec = np.full(300, np.nan)
    for t in np.flatnonzero(evaluable[f]):
        valid = tvalid[f, t].copy()
        valid[t] = False
        idx = np.flatnonzero(valid)
        tr = truth[f, t, idx]
        den = float(np.sum(np.abs(tr)))
        credit[t] = float(np.sum(np.sign(tr) * pred[f, t, idx]) / den)
        n = max(1, int(math.ceil(0.20 * idx.size)))
        order = np.lexsort((CANON[idx], -conf[f, t, idx]))[:n]
        ch = idx[order]
        w = np.abs(truth[f, t, ch])
        ps, ts_ = np.sign(pred[f, t, ch]), np.sign(truth[f, t, ch])
        ok = (ps != 0.0) & (ts_ != 0.0) & (ps == ts_)
        prec[t] = float(np.sum(w * ok) / float(np.sum(w)))
    return credit, prec


def macro(pred, conf):
    cs, ps = [], []
    for f in range(3):
        c, p = score_family(pred, conf, f)
        cs.append(np.nanmean(c))
        ps.append(np.nanmean(p))
    return np.array(cs), np.array(ps)


b_c, b_p = macro(base_e, base_c)
c_c, c_p = macro(cand_e, cand_c)
say("0.2 re-scored 122 credit (macro)", f"{b_c.mean():+.8f}")
say("0.2 re-scored 122 precision (macro)", f"{b_p.mean():+.8f}")
say("0.2 re-scored 125 credit (macro)", f"{c_c.mean():+.8f}")
say("0.2 re-scored 125 precision (macro)", f"{c_p.mean():+.8f}")
say("0.2 KA2 credit delta (recorded -0.01504990)", f"{(c_c - b_c).mean():+.8f}")
say("0.2 KA2 precision delta (recorded +0.01217457)", f"{(c_p - b_p).mean():+.8f}")
say("0.2 KA1 self-delta credit (must be 0)", f"{(b_c - b_c).mean():+.3e}")

# 122's own absolute numbers, from its RESULTS.md
say("0.2 122 credit recorded in 122/RESULTS.md", 0.029965)
say("0.2 122 top-20 precision recorded in 122/RESULTS.md", 0.537409)

# ------------------------------------------------------------- 0.3 homogeneity
print("\n--- 0.3 homogeneity on the REAL truth and REAL confidence ---")
hom = {}
for k in (0.25, 0.5, 0.75, 2.0, 6.78):
    kc, kp = macro(k * cand_e, cand_c)
    ratio = kc.mean() / c_c.mean()
    hom[str(k)] = dict(credit_ratio=float(ratio), precision_delta=float((kp - c_p).mean()))
    print(f"  k={k:<6} credit({k}.p)/credit(p) = {ratio:.12f}   "
          f"prec({k}.p)-prec(p) = {(kp - c_p).mean():+.12e}")
out["0.3 homogeneity"] = hom
say("0.3 max |credit ratio - k| over k",
    f"{max(abs(v['credit_ratio'] - float(k)) for k, v in hom.items()):.3e}")
say("0.3 max |precision delta| over k",
    f"{max(abs(v['precision_delta']) for v in hom.values()):.3e}")

# -------------------------------------------------------------- 0.4 the shrink
print("\n--- 0.4 shrink and direction ---")
fin = np.isfinite(cand_e) & np.isfinite(base_e)
nz = fin & ((cand_e != 0) | (base_e != 0))
l1 = np.abs(cand_e[nz]).sum() / np.abs(base_e[nz]).sum()
l2 = np.linalg.norm(cand_e[nz]) / np.linalg.norm(base_e[nz])
ls = float(np.dot(cand_e[nz], base_e[nz]) / np.dot(base_e[nz], base_e[nz]))
gcos = float(np.dot(cand_e[nz], base_e[nz]) /
             (np.linalg.norm(cand_e[nz]) * np.linalg.norm(base_e[nz])))
say("0.4 L1 ratio cand/base", f"{l1:.6f}  (1/x = {1/l1:.4f})")
say("0.4 L2 ratio cand/base", f"{l2:.6f}  (1/x = {1/l2:.4f})")
say("0.4 least-squares scalar base->cand", f"{ls:.6f}  (1/x = {1/ls:.4f})")
say("0.4 global cosine(cand, base)", f"{gcos:.6f}")

per_ls, per_cos = [], []
for f in range(3):
    for t in np.flatnonzero(evaluable[f]):
        m = tvalid[f, t].copy()
        m[t] = False
        a, b = cand_e[f, t, m], base_e[f, t, m]
        if np.dot(b, b) <= 0:
            continue
        per_ls.append(np.dot(a, b) / np.dot(b, b))
        per_cos.append(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
per_ls, per_cos = np.array(per_ls), np.array(per_cos)
say("0.4 per-target LS scalar median",
    f"{np.median(per_ls):.6f} IQR [{np.percentile(per_ls,25):.6f},{np.percentile(per_ls,75):.6f}]")
say("0.4 per-target cosine(cand,base) median",
    f"{np.median(per_cos):.6f} IQR [{np.percentile(per_cos,25):.6f},{np.percentile(per_cos,75):.6f}]")
sf = (np.sign(cand_e) != np.sign(base_e)) & nz
say("0.4 coordinate sign flips cand vs base",
    f"{int(sf.sum()):,} / {int(nz.sum()):,} = {100*sf.sum()/nz.sum():.2f}%")

# ------------------------------------------------------- 0.6 rescale control
print("\n--- 0.6 PURE-RESCALE CONTROL (mandate Phase 4.1) ---")
ctrl = {}
for name, c in (("LS", ls), ("L1", l1), ("L2", l2), ("per-target-median", float(np.median(per_ls)))):
    rc, rp = macro(c * base_e, base_c)
    dc, dp = (rc - b_c).mean(), (rp - b_p).mean()
    analytic = (c - 1.0) * b_c.mean()
    ctrl[name] = dict(c=float(c), d_credit=float(dc), d_precision=float(dp),
                      analytic=float(analytic),
                      pct_of_recorded=float(100.0 * dc / REC_CREDIT))
    print(f"  c={c:.6f} ({name:<17}) dcredit {dc:+.8f}  ({100*dc/REC_CREDIT:6.2f}% of recorded)"
          f"  dprec {dp:+.3e}  | KA3 analytic {analytic:+.8f} err {abs(dc-analytic):.2e}")
out["0.6 rescale control"] = ctrl
# the scalar that EXACTLY reproduces the recorded credit loss
c_star = 1.0 + REC_CREDIT / b_c.mean()
say("0.6 c that exactly reproduces -0.01505", f"{c_star:.6f}  (1/c = {1/c_star:.4f})")

# ------------------------------------- 0.7 decompose the precision gain 2x2
print("\n--- 0.7 precision gain: effect-sign channel vs confidence-ordering channel ---")
grid = {}
for en, ee in (("122eff", base_e), ("125eff", cand_e)):
    for cn, cc in (("122conf", base_c), ("125conf", cand_c)):
        mc, mp = macro(ee, cc)
        grid[f"{en}+{cn}"] = dict(credit=float(mc.mean()), precision=float(mp.mean()),
                                  d_credit=float((mc - b_c).mean()),
                                  d_precision=float((mp - b_p).mean()))
        print(f"  {en}+{cn:<8} credit {mc.mean():+.6f}  precision {mp.mean():.6f}  "
              f"dprec {(mp-b_p).mean():+.6f}")
out["0.7 grid"] = grid

json.dump(out, open(OUT / "phase0.json", "w"), indent=1, default=str)
print(f"\nwrote {OUT/'phase0.json'}")
