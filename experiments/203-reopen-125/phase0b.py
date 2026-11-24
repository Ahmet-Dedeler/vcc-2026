"""203 Phase 0b — is 125's DIRECTION better than 122's, in 125's own frame?

Phase 0 proved the recorded credit loss is scale.  That leaves the question the mandate
actually cares about: after the scale is divided out, does the hierarchical posterior
point in a better direction?  Measured here on 125's own frozen 300x300 panel frame with
the four scale-free instruments the competition reads through:

  cos       per-target cosine(prediction, truth)              -> feeds `pds` and `mse`
  pds*      panel-frame retrieval rank (own truth vs decoys)   -> the `pds` shape
  q_w       |truth|-weighted sign accuracy over ALL coords     -> `nmae` moves only here (187 5)
  credit142 142's homonym: signed credit / L1(PREDICTION)      -> the scale-invariant twin

plus the direction-only variant (122's effect, 125's confidence), which by construction is
identical to 122 on every one of those four and can only move the ranked endpoints.

KNOWN-ANSWER ROWS:
  KA1  credit142 must be EXACTLY invariant to rescaling the prediction (it is 142's
       scale-invariant definition); the 125-frame credit must scale by exactly k.
  KA2  the rescale control (c * 122) must reproduce 122's cos / pds* / q_w to 0.0 --
       all three are scale-free, so a nonzero delta means the instrument is wired wrong.
  KA3  q_w of a sign-randomised prediction must sit at ~0.5.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W125 = ROOT / "experiments/125-hierarchical-panel-posterior/work/source-only-evaluation-v1"
W122 = ROOT / "experiments/122-panel-response-atlas/work/source-only-phase1-v1"
OUT = HERE / "work"

truth = np.asarray(np.load(W122 / "family-effects.npy"), dtype=np.float64)
tvalid = np.asarray(np.load(W122 / "family-valid.npy"), dtype=bool)
cand_e = np.asarray(np.load(W125 / "observed-candidate-effect.npy"), dtype=np.float64)
base_e = np.asarray(np.load(W125 / "observed-baseline-effect.npy"), dtype=np.float64)
cred_d = np.load(W125 / "per-target-credit-delta.npy")
prec_d = np.load(W125 / "per-target-precision-delta.npy")
evaluable = np.isfinite(cred_d) & np.isfinite(prec_d)
FAMS = ("replogle_k562", "orion", "cd4")


def instruments(pred, f):
    """cos, q_w, credit142 per target; and the panel-frame retrieval rank."""
    rows = np.flatnonzero(evaluable[f])
    cos, qw, c142, craw = [], [], [], []
    P = np.zeros((len(rows), 300))
    T = np.zeros((len(rows), 300))
    for i, t in enumerate(rows):
        m = tvalid[f, t].copy()
        m[t] = False
        p, tr = pred[f, t], truth[f, t]
        idx = np.flatnonzero(m)
        pv, tv = p[idx], tr[idx]
        np_, nt = np.linalg.norm(pv), np.linalg.norm(tv)
        cos.append(float(pv @ tv / (np_ * nt)) if np_ > 0 else 0.0)
        live = tv != 0
        w = np.abs(tv) * live
        ok = (np.sign(pv) == np.sign(tv)) & live
        qw.append(float(np.where(ok, w, 0.0).sum() / max(w.sum(), 1e-300)))
        l1p = float(np.abs(pv).sum())
        sc = float(np.sum(np.sign(tv) * pv))
        c142.append(sc / l1p if l1p > 0 else 0.0)
        craw.append(sc / float(np.abs(tv).sum()))
        # retrieval frame: shared coordinate support across the family's evaluable rows
        P[i] = np.where(m, p, 0.0)
        T[i] = np.where(m, tr, 0.0)
    A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    B = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S)
    beaten = (S < own[:, None]).sum(1)
    ties = (S == own[:, None]).sum(1) - 1
    pds = float(np.mean((beaten + 0.5 * ties) / (n - 1)))
    return dict(cos=float(np.mean(cos)), qw=float(np.mean(qw)),
                credit142=float(np.mean(c142)), credit125=float(np.mean(craw)),
                pds=pds, n=n)


def macro(pred):
    per = [instruments(pred, f) for f in range(3)]
    keys = ("cos", "qw", "credit142", "credit125", "pds")
    return {k: float(np.mean([p[k] for p in per])) for k in keys}, per


ARMS = {
    "122 equal_mean_agreement": base_e,
    "125 hierarchical posterior": cand_e,
    "122 x 0.147491 (rescale control)": 0.147491 * base_e,
    "125 x 6.7801 (unshrunk 125)": 6.7801 * cand_e,
}
res, per_fam = {}, {}
for name, arr in ARMS.items():
    res[name], per_fam[name] = macro(arr)

hdr = f"{'arm':<34}{'cos':>10}{'pds*':>10}{'q_w':>10}{'credit142':>12}{'credit125':>12}"
print(hdr)
print("-" * len(hdr))
for name in ARMS:
    r = res[name]
    print(f"{name:<34}{r['cos']:>10.5f}{r['pds']:>10.5f}{r['qw']:>10.5f}"
          f"{r['credit142']:>12.6f}{r['credit125']:>12.6f}")

b, c = res["122 equal_mean_agreement"], res["125 hierarchical posterior"]
print()
print(f"{'DELTA 125 - 122':<34}{c['cos']-b['cos']:>+10.5f}{c['pds']-b['pds']:>+10.5f}"
      f"{c['qw']-b['qw']:>+10.5f}{c['credit142']-b['credit142']:>+12.6f}"
      f"{c['credit125']-b['credit125']:>+12.6f}")

rc = res["122 x 0.147491 (rescale control)"]
print("\nKA2 rescale control vs 122 (all four scale-free instruments must be 0.0):")
for k in ("cos", "pds", "qw", "credit142"):
    print(f"   d{k:<10} {rc[k]-b[k]:+.3e}")
print(f"   KA1 credit125 ratio  {rc['credit125']/b['credit125']:.12f}  (must be 0.147491)")

rng = np.random.default_rng(0)
rand = np.abs(base_e) * rng.choice([-1.0, 1.0], size=base_e.shape)
print(f"KA3 q_w of sign-randomised prediction  {macro(rand)[0]['qw']:.5f}   (must be ~0.5)")

print("\nper family (delta 125 - 122):")
print(f"{'family':<16}{'d cos':>10}{'d pds*':>10}{'d q_w':>10}{'d credit142':>13}")
for f in range(3):
    pb, pc = per_fam["122 equal_mean_agreement"][f], per_fam["125 hierarchical posterior"][f]
    print(f"{FAMS[f]:<16}{pc['cos']-pb['cos']:>+10.5f}{pc['pds']-pb['pds']:>+10.5f}"
          f"{pc['qw']-pb['qw']:>+10.5f}{pc['credit142']-pb['credit142']:>+13.6f}")

json.dump(dict(macro=res, per_family={k: v for k, v in per_fam.items()}),
          open(OUT / "phase0b.json", "w"), indent=1)
print(f"\nwrote {OUT/'phase0b.json'}")
