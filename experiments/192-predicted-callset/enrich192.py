"""192 Phase 1 — THE ENRICHMENT CURVE.  The number this project has never had.

For a key `k` over (perturbation, gene), the call set "top n of row p by k" has

    prec(n)   = mean_p |top_n(p) ∩ SIG_p| / n            <- what `fid`'s precision rides on
    rec(n)    = mean_p |top_n(p) ∩ SIG_p| / n_real_p     <- what `jac`'s numerator rides on
    fold(n)   = prec(n) / base(n),   base(n) = mean_p n_real_p / G

`fold` is the membership enrichment INBOX-5 §5 says must clear ~2x before cutting the call set
can raise `jac`.  Rows with n_real_p = 0 are dropped from every row-average (stated, not hidden).

The keys are split into three classes and NEVER mixed in a headline:
  BUILDABLE  -- computable from our prediction + the challenge's own control cells
  ORACLE     -- uses the reference's own lfc; bounds the family (188's 0.7219 is one of these)
  CONTROL    -- shuffled, the Phase-3 falsifier 1 row, evaluated at every n
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
E169 = os.path.join(REPO, "experiments/169-sharpening-six-member")
E185 = os.path.join(REPO, "experiments/185-rank-controlled-emission")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
sys.path.insert(0, HERE)
from debench import Surface                                    # noqa: E402
from keys192 import auc_key, gene_stats, h_to_lfc, shift_auc_table  # noqa: E402

NS = [25, 50, 89, 150, 300, 600, 1200, 2400, 4800]
SHIFTS = np.linspace(-4.0, 4.0, 33)


def topn_mask(key, n_per_row):
    out = np.zeros(key.shape, bool)
    for i in range(key.shape[0]):
        n = int(n_per_row[i])
        if n <= 0:
            continue
        n = min(n, key.shape[1])
        out[i, np.argpartition(-key[i], n - 1)[:n]] = True
    return out


def curve(key, sig, nreal, ns=NS):
    """prec / rec / fold at each fixed n, plus the n = n_real_p row (188's statistic)."""
    G = key.shape[1]
    rows = np.flatnonzero(nreal > 0)
    out = {}
    for n in ns:
        t = topn_mask(key, np.full(key.shape[0], n))
        inter = np.array([(t[i] & sig[i]).sum() for i in rows], float)
        out[n] = dict(prec=float((inter / min(n, G)).mean()),
                      rec=float((inter / nreal[rows]).mean()),
                      base=float((nreal[rows] / G).mean()))
        out[n]["fold"] = out[n]["prec"] / out[n]["base"]
    t = topn_mask(key, nreal)
    inter = np.array([(t[i] & sig[i]).sum() for i in rows], float)
    out["nreal"] = dict(prec=float((inter / nreal[rows]).mean()),
                        rec=float((inter / nreal[rows]).mean()),
                        base=float((nreal[rows] / G).mean()))
    out["nreal"]["fold"] = out["nreal"]["prec"] / out["nreal"]["base"]
    return out


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    E = np.load(os.path.join(E169, "work/h1test-effect.npz"), allow_pickle=True)
    pos = {str(g): i for i, g in enumerate(E["genes"])}
    tpos = {str(t): i for i, t in enumerate(E["targets"])}
    s = Surface(os.path.join(E185, "work/de185.a15.npz"))
    z15 = np.load(os.path.join(E185, "work/de185.a15.npz"), allow_pickle=True)
    gi = np.array([pos[g] for g in s.genes])
    ti = np.array([tpos[t] for t in s.perts])
    cpm = E["cpm"][gi].astype(np.float64)                 # pooled control CPM, exact
    Hm = E["H"][np.ix_(ti, gi)].astype(np.float64)        # our shipped effect, h space
    sig = s.real_sig & np.isfinite(s.real_lfc)
    nreal = sig.sum(1)
    G = len(s.genes)

    C = np.load(os.path.join(E185, "work/ctrl_meta.npz"), allow_pickle=True)
    ref_lin = C["ref_lin"].astype(np.float64)             # mean per-cell control CPM
    gs = gene_stats(os.path.join(E185, "work/ctrl_cs.npy"),
                    os.path.join(WORK, "genestats.npz"))
    print(f"control gene stats: detection rate median {np.median(gs['nz']):.3f}, "
          f"cv2 median {np.median(gs['cv2']):.2f}", flush=True)
    A = shift_auc_table(os.path.join(E185, "work/ctrl_cs.npy"), gs["grid"], SHIFTS,
                        os.path.join(WORK, "aucshift.npz"))
    print(f"shift-AUC table {A.shape} ready", flush=True)

    lfcP = h_to_lfc(Hm, cpm)                              # our predicted lfc, exact from H
    lfcR = np.nan_to_num(s.real_lfc)                      # the reference's own lfc (oracle)
    Rh = np.nan_to_num(z15["absH"]).astype(np.float64)    # == |H| on the gated axis

    keys = {}
    # ---------------- BUILDABLE
    keys["B |H|"] = np.abs(Hm)
    keys["B |lfc_pred|"] = np.abs(lfcP)
    keys["B |lfc_pred|*sqrt(cpm)"] = np.abs(lfcP) * np.sqrt(cpm)[None, :]
    keys["B |lfc_pred|*sqrt(reflin)"] = np.abs(lfcP) * np.sqrt(ref_lin)[None, :]
    keys["B shiftAUC(lfc_pred)"] = auc_key(lfcP, A, SHIFTS)
    keys["B sqrt(cpm) only"] = np.broadcast_to(np.sqrt(cpm)[None, :], Hm.shape).copy()
    keys["B detection rate only"] = np.broadcast_to(gs["nz"][None, :], Hm.shape).copy()
    keys["B a15 emitted p-order"] = -np.log(np.maximum(z15["pred_p"], 1e-300))

    # ---- gene-level DE propensity.  The reference's DE set is a DETECTABILITY set, so the
    # strongest buildable signal may be a per-gene prior.  Two versions:
    #   LOTO  -- estimated on the OTHER 84 targets of the same line.  NOT buildable on the
    #            challenge (we never see its reference); it bounds the gene-level family.
    #   Orion -- estimated on a DIFFERENT cell line's public reference.  Buildable.
    rate_loto = np.zeros(Hm.shape)
    tot = sig.sum(0).astype(float)
    for p in range(len(s.perts)):
        rate_loto[p] = (tot - sig[p]) / (len(s.perts) - 1)
    keys["B gene DE-rate LOTO (same line)"] = rate_loto
    oz = np.load(os.path.join(REPO, "experiments/095-clean-harness/de_thinned.npz"),
                 allow_pickle=True)
    osig = (oz["real_padj"] < 0.05) & np.isfinite(oz["real_padj"])
    orate = dict(zip([str(g) for g in oz["genes"]], osig.mean(0)))
    ov = np.array([orate.get(g, np.nan) for g in s.genes])
    print(f"    Orion HCT116 DE-rate covers {int(np.isfinite(ov).sum())}/{G} h1-test genes; "
          f"missing filled at the observed median")
    ov = np.where(np.isfinite(ov), ov, np.nanmedian(ov))
    keys["B gene DE-rate Orion (cross-line)"] = np.broadcast_to(ov[None, :], Hm.shape).copy()
    keys["B Orion rate x |lfc_pred|"] = ov[None, :] * np.abs(lfcP)
    # ---------------- ORACLE
    keys["O |lfc_real|"] = np.abs(lfcR)
    keys["O |h_real|"] = np.abs(Rh) * 0 + np.abs(
        np.log1p(50000 * 1e-6 * cpm[None, :] * np.exp2(lfcR))
        - np.log1p(50000 * 1e-6 * cpm)[None, :])
    keys["O |lfc_real|*sqrt(cpm)"] = np.abs(lfcR) * np.sqrt(cpm)[None, :]
    keys["O shiftAUC(lfc_real)"] = auc_key(lfcR, A, SHIFTS)
    # ---- fitted membership predictor, leave-one-TARGET-out (never fits on the scored row)
    from sklearn.linear_model import LogisticRegression
    lo = np.log
    X = np.stack([np.broadcast_to(lo(cpm + 1)[None, :], Hm.shape),
                  np.broadcast_to(lo(gs["nz"][:] + 1e-3)[None, :], Hm.shape),
                  np.broadcast_to(lo(gs["cv2"][:] + 1e-3)[None, :], Hm.shape),
                  lo(np.abs(lfcP) + 1e-3), lo(np.abs(Hm) + 1e-6)], -1)
    Xg = X[:, :, :3]                                       # gene-only sub-model
    mu, sd = X.reshape(-1, X.shape[-1]).mean(0), X.reshape(-1, X.shape[-1]).std(0)
    Xz = (X - mu) / sd
    rngf = np.random.default_rng(7)
    for lab, cols in (("B fit LOTO (gene+pred)", [0, 1, 2, 3, 4]),
                      ("B fit LOTO (gene only)", [0, 1, 2])):
        out = np.zeros(Hm.shape)
        for p in range(len(s.perts)):
            tr = np.setdiff1d(np.arange(len(s.perts)), [p])
            Xt = Xz[tr][:, :, cols].reshape(-1, len(cols))
            yt = sig[tr].reshape(-1)
            idx = rngf.choice(Xt.shape[0], size=min(300_000, Xt.shape[0]), replace=False)
            m = LogisticRegression(max_iter=400, C=1.0).fit(Xt[idx], yt[idx])
            out[p] = m.decision_function(Xz[p][:, cols])
        keys[lab] = out
        print(f"    fitted {lab}", flush=True)
    del X, Xz, Xg

    # ---------------- CONTROL
    rng = np.random.default_rng(192)
    keys["C shuffled |H| within row"] = np.apply_along_axis(rng.permutation, 1, np.abs(Hm))
    keys["C uniform random"] = rng.random(Hm.shape)

    print(f"\nh1-test: {len(s.perts)} perturbations x {G} gated genes; "
          f"n_real median {np.median(nreal):.0f} mean {nreal.mean():.0f} "
          f"({int((nreal == 0).sum())} rows with n_real = 0, dropped)")
    res = {}
    hdr = "".join(f"{n:>9}" for n in NS) + f"{'n_real_p':>10}"
    for lab, tag in (("PRECISION p(n)  = mean_p |top_n ∩ SIG_p| / n", "prec"),
                     ("FOLD-ENRICHMENT p(n) / (mean_p n_real_p / G)", "fold"),
                     ("RECALL  r(n)    = mean_p |top_n ∩ SIG_p| / n_real_p", "rec")):
        print(f"\n{lab}")
        print(f"  {'key':<28}{hdr}")
        for k, kv in keys.items():
            if k not in res:
                res[k] = curve(kv, sig, nreal)
            c = res[k]
            fmt = "{:>9.4f}" if tag != "fold" else "{:>9.2f}"
            print(f"  {k:<28}" + "".join(fmt.format(c[n][tag]) for n in NS)
                  + ("{:>10.4f}".format(c["nreal"][tag]) if tag != "fold"
                     else "{:>10.2f}".format(c["nreal"][tag])))

    with open(os.path.join(WORK, "enrich192.json"), "w") as fh:
        json.dump({k: {str(n): v for n, v in c.items()} for k, c in res.items()}, fh,
                  indent=1, default=float)
    np.savez_compressed(os.path.join(WORK, "keycache.npz"),
                        cpm=cpm, ref_lin=ref_lin, lfcP=lfcP, H=Hm, nreal=nreal,
                        nz=gs["nz"][:], cv2=gs["cv2"][:], genes=s.genes, perts=s.perts)
    print(f"\nwrote {WORK}/enrich192.json and keycache.npz")


if __name__ == "__main__":
    main()
