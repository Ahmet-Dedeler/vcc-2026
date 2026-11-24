"""161 SS3 — the cross-line RESIDUAL predictor.

THE OBJECT
----------
Production emits ONE context-independent effect matrix.  159 measured that only ~30% of a
target's target-specific response energy is shared across cell lines, so ~70% is line-specific
and structurally inexpressible by that matrix.  160 showed the line-specific part is not
readable off the recipient's control-cell geometry.  Every one of this project's ~40 NO-GO
experiments tried to improve the CONSENSUS.  This one predicts the RESIDUAL.

    R_L  = E_L - mean_t(E_L)            the recipient's target-specific truth
    C    = mean over donor LINES of R_d, equal weights (AGENTS.md §4.2)
    y    = R_L - C                      the line-specific residual -- what we fit

ARCHITECTURE CONSTRAINT (hard)
------------------------------
A universal function of measurable per-(target, gene, context) SCALARS.  No gene identity, no
target identity, no line identity.  Anything else cannot generalise to the challenge's unseen
lines.  `AgingLab_BI` shipped this exact shape with four TARGET-side scalars and LOST on the
board (0.1746 -> 0.1673).  159's claim is that the missing information is CONTEXT, so the
features are split into families and ablated.

    T  target-side   consensus magnitude and sign, cross-donor disagreement and sign
                     agreement, the target's consensus norm, the gene's donor
                     responsiveness, saturation.
    C  context       recipient basal, donor-mean basal, their difference and spread, the
                     target gene's own basal in the recipient (knockdown headroom), the
                     recipient/donor control-profile divergence, AND the two scalars a donor
                     table cannot express -- per-gene control-cell DETECTION RATE and
                     DISPERSION in the recipient (`ctrlstats.py`).
    X  interaction   consensus x context.  159's claim has to live HERE: see below.

WHY THE INTERACTION IS THE WHOLE EXPERIMENT
-------------------------------------------
`yhat` is centred across targets before use, because a target-independent addition is the
scorer's own zero AND a `pds` hub -- that is 160's entire result.  A feature that is a function
of the GENE alone therefore predicts a value constant across targets, and centring annihilates
it EXACTLY.  So a context-only model is identically zero and 159's "context is worth 3x target"
cannot be tested through context main effects.  It can only appear as `consensus x context`:
a context-dependent RE-WEIGHTING of the consensus.  (The two `btg` features are functions of
the TARGET alone, so they survive centring as a constant-per-row shift; the C family is
therefore not exactly zero, and the printed C row shows how little is left.)

SCORING IS IN THE DELIVERED SPACE (156)
---------------------------------------
gate -> production `mass_consistent_effect` -> render -> 156 `endpoints`.  An arm scored on the
pre-render effect measures a different object (+0.0032/+0.0038 raw pds of difference).

EVERY ENDPOINT IS SCALE-INVARIANT, so arms are compared at RESTORED ROW NORMS: `w` changes the
DIRECTION only.  Without that, a larger `w` both rotates and inflates the effect, and inflation
alone moves the delivered geometry (156: the render residual is ~4% of the row norm and is 75x
more mutually-aligned than the signal).  That would be two variables sold as one.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "156-delivered-geometry"))
from screen import endpoints, render                                    # noqa: E402

from common import SEED, WORK                                           # noqa: E402

GATE_CPM = 10.0
_PROD = None


def prod():
    global _PROD
    if _PROD is None:
        import importlib.util
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "..", "..", "scripts/make_submission_pm_cd4.py")
        spec = importlib.util.spec_from_file_location("prod_builder_161", path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules["prod_builder_161"] = mod
        spec.loader.exec_module(mod)
        assert mod.MASS_FIX is None, "production default changed; re-read the builder"
        _PROD = mod
    return _PROD


def mass_fix(D, b, gate):
    p = prod()
    out = np.empty_like(D)
    for i in range(D.shape[0]):
        out[i] = p.mass_consistent_effect(D[i], b, gate, mode="minnorm")
    return out


def deliver(D, b, gate):
    """gate -> production mass fix -> render.  Exactly the production path (156)."""
    Dg = np.where(gate[None, :], D, 0.0)
    return render(mass_fix(Dg, b, gate), b)


def restore_norms(Q, P):
    """Put each row of Q on the row norm it had in P, so an arm changes DIRECTION only."""
    nq = np.maximum(np.linalg.norm(Q, axis=1, keepdims=True), 1e-300)
    return Q * (np.linalg.norm(P, axis=1, keepdims=True) / nq)


def split_par_orth(Y, C):
    """Split each row of Y into its component along C and the rest.

    THIS IS THE LOAD-BEARING PARAMETERISATION, and the first run of this screen got it
    wrong.  The fitted residual comes out at cos(yhat, C) = -0.99: the model overwhelmingly
    learns the SHRINKAGE `y = R_L - C ~= -C`, because only ~30% of the consensus transfers
    (159).  Every endpoint here unit-normalises each row, so a per-target rescale of the
    prediction is EXACTLY inert -- the parallel component cannot move `pds` or `cos` at all.

    Left in, it wrecks both the arm and its null:
      * the arm `C + w*yhat` becomes `(1 - 0.94w) C + w*u`, which passes through zero near
        w = 1.06 and then FLIPS THE SIGN of the consensus.  That is what the -0.4 pds at
        w = 2.0 was; `w` was mixing "add new direction" with "shrink and then invert".
      * the permutation null then reads hugely negative for a trivial reason -- it is
        splicing in ANOTHER target's consensus with a negative sign -- so it would look
        like a passed control even if the model carried no target-specific information at
        all.  That is precisely the failure mode AGENTS.md flags: ask what the control
        would read if the thing being tested were absent.

    So the arm is built from the ORTHOGONAL component only, and the parallel component is
    scored separately as a known-answer row that must read 0.
    """
    d = np.maximum((C * C).sum(1, keepdims=True), 1e-300)
    a = (Y * C).sum(1, keepdims=True) / d
    par = a * C
    return par, Y - par


# ------------------------------------------------------------------ the surface
class Surface:
    def __init__(self, name, ctrl=True):
        z = np.load(os.path.join(WORK, f"surface-{name}.npz"), allow_pickle=True)
        self.name = name
        self.targets = [str(t) for t in z["targets"]]
        self.genes = [str(g) for g in z["genes"]]
        self.members = [str(m) for m in z["members"]]
        self.mline = {m: str(l) for m, l in zip(self.members, z["lines"])}
        self.lines = sorted(set(self.mline.values()))
        self.E, self.b, self.gate, self.btg = {}, {}, {}, {}
        for m in self.members:
            self.E[m] = z[f"E_{m}"].astype(np.float64)
            self.b[m] = z[f"b_{m}"]
            self.gate[m] = z[f"gate_{m}"]
            self.btg[m] = z[f"btg_{m}"]
        # per LINE: average the tables of that line (two K562 libraries are ONE line)
        self.R, self.bL, self.gateL, self.btgL, self.EL = {}, {}, {}, {}, {}
        for L in self.lines:
            mem = [m for m in self.members if self.mline[m] == L]
            E = np.mean([self.E[m] for m in mem], 0)
            self.EL[L] = E
            self.R[L] = E - E.mean(0, keepdims=True)
            self.bL[L] = np.mean([self.b[m] for m in mem], 0)
            self.gateL[L] = np.all([self.gate[m] for m in mem], 0)
            self.btgL[L] = np.nanmean(np.stack([self.btg[m] for m in mem]), 0)
        self.nT, self.nG = len(self.targets), len(self.genes)
        self.ctrl = self._load_ctrl() if ctrl else {}

    def _load_ctrl(self):
        """Per-gene control-cell detection rate and dispersion, on THIS surface's gene axis."""
        out = {}
        for L in self.lines:
            p = os.path.join(WORK, f"ctrl-{L}.npz")
            if not os.path.exists(p):
                print(f"  [ctrl] {L}: absent -- context cell-level features disabled")
                return {}
            z = np.load(p, allow_pickle=True)
            pos = {str(g): i for i, g in enumerate(z["genes"])}
            idx = np.asarray([pos.get(g, -1) for g in self.genes])
            hit = idx >= 0
            d = {}
            for k in ("det", "vlog", "fano", "mlog"):
                v = np.full(self.nG, np.nan)
                v[hit] = np.asarray(z[k])[idx[hit]]
                d[k] = np.nan_to_num(v, nan=float(np.nanmedian(v)))
            # dispersion BEYOND the mean-variance trend: regress vlog on mlog (quadratic)
            # and keep the residual, so the feature is not a repackaging of the pseudobulk.
            A = np.stack([np.ones(self.nG), d["mlog"], d["mlog"] ** 2], 1)
            coef, *_ = np.linalg.lstsq(A, d["vlog"], rcond=None)
            d["xvar"] = d["vlog"] - A @ coef
            d["n_cells"] = int(z["n_cells"])
            d["cover"] = float(hit.mean())
            out[L] = d
        cov = ", ".join(f"{L} {out[L]['n_cells']:,}c/{out[L]['cover']:.0%}g" for L in self.lines)
        print(f"  [ctrl] cell-level control stats loaded: {cov}")
        return out

    def consensus(self, donors):
        Rs = np.stack([self.R[d] for d in donors])
        return Rs.mean(0), Rs.std(0), np.sign(Rs).mean(0)


# ---------------------------------------------------------------------- features
# Families are defined by WHAT THEY REQUIRE, which is the only split that answers the
# mandate's question ("do CONTEXT features beat target-side ones?"):
#   T  computable WITHOUT the recipient's basal profile -- donor consensus statistics and
#      donor basal.  This is exactly the shape `AgingLab_BI` shipped, and it LOST.
#   C  requires the recipient's context, entering additively.
#   X  consensus x recipient-context interaction.
# `sat` = C/(1+b_recipient) was labelled T in the first draft.  It is definitionally an
# interaction with the recipient's basal profile, and mislabelling it would have moved
# credit from the context families to the target family -- i.e. corrupted the headline.
BASE_FAMILY = {
    "cons": "T", "acons": "T", "dsd": "T", "sign": "T", "tmag": "T", "gresp": "T",
    "bD": "T", "bsd": "T", "x_cons_gresp": "T",
    "bB": "C", "bdiff": "C", "btgB": "C", "btgdiff": "C", "div": "C",
    "sat": "X", "x_cons_bB": "X", "x_cons_bdiff": "X", "x_cons_btgB": "X",
}
CTRL_FAMILY = {                                     # only when ctrl-*.npz exist
    "cdet": "C", "cxvar": "C", "cfano": "C", "cdetdiff": "C", "cxvardiff": "C",
    "x_cons_cdet": "X", "x_cons_cxvar": "X", "x_cons_cdetdiff": "X", "x_cons_cxvardiff": "X",
}
CTX_CELL = set(CTRL_FAMILY)                          # the cell-level-only block


def feature_names(has_ctrl: bool):
    f = dict(BASE_FAMILY)
    if has_ctrl:
        f.update(CTRL_FAMILY)
    return list(f), f


def features(S: Surface, recipient: str, donors: list[str], names: list[str]):
    """(nT*nG, F) float32.  Every column is a measurable scalar; none is an identity."""
    nT, nG = S.nT, S.nG
    C, DSD, SGN = S.consensus(donors)
    bB = S.bL[recipient]
    bD = np.mean([S.bL[d] for d in donors], 0)
    bSD = np.std([S.bL[d] for d in donors], 0)
    btgB = np.nan_to_num(S.btgL[recipient], nan=float(np.nanmedian(S.btgL[recipient])))
    btgD = np.nanmean(np.stack([S.btgL[d] for d in donors]), 0)
    btgD = np.nan_to_num(btgD, nan=float(np.nanmedian(btgD)))
    gresp = np.mean([np.abs(S.R[d]).mean(0) for d in donors], 0)
    tmag = np.linalg.norm(C, axis=1)
    # recipient/donor control-profile divergence: ONE scalar per (recipient, fold).
    # It is line identity in disguise -- flagged, in family C, ablated by `ALL-div`.
    div = float(1.0 - np.dot(bB, bD) / (np.linalg.norm(bB) * np.linalg.norm(bD)))

    one_t, one_g = np.ones((nT, 1)), np.ones((1, nG))
    cols = {
        "bB": one_t * bB[None, :], "bD": one_t * bD[None, :],
        "bdiff": one_t * (bB - bD)[None, :], "bsd": one_t * bSD[None, :],
        "btgB": btgB[:, None] * one_g, "btgdiff": (btgB - btgD)[:, None] * one_g,
        "div": np.full((nT, nG), div),
        "cons": C, "acons": np.abs(C), "dsd": DSD, "sign": SGN,
        "tmag": tmag[:, None] * one_g, "gresp": one_t * gresp[None, :],
        "sat": C / (1.0 + one_t * bB[None, :]),
        "x_cons_bdiff": C * (bB - bD)[None, :], "x_cons_bB": C * bB[None, :],
        "x_cons_gresp": C * gresp[None, :], "x_cons_btgB": C * btgB[:, None],
    }
    if S.ctrl:
        det, xv, fa = (S.ctrl[recipient][k] for k in ("det", "xvar", "fano"))
        detD = np.mean([S.ctrl[d]["det"] for d in donors], 0)
        xvD = np.mean([S.ctrl[d]["xvar"] for d in donors], 0)
        cols.update({
            "cdet": one_t * det[None, :], "cxvar": one_t * xv[None, :],
            "cfano": one_t * np.log1p(np.maximum(fa, 0))[None, :],
            "cdetdiff": one_t * (det - detD)[None, :],
            "cxvardiff": one_t * (xv - xvD)[None, :],
            "x_cons_cdet": C * det[None, :], "x_cons_cxvar": C * xv[None, :],
            "x_cons_cdetdiff": C * (det - detD)[None, :],
            "x_cons_cxvardiff": C * (xv - xvD)[None, :],
        })
    X = np.empty((nT * nG, len(names)), dtype=np.float32)
    for j, n in enumerate(names):
        X[:, j] = np.asarray(cols[n], dtype=np.float32).reshape(-1)
    y = (S.R[recipient] - C).astype(np.float32).reshape(-1)
    return X, y, C


def fit(X, y, cols, seed=SEED, max_iter=150):
    from sklearn.ensemble import HistGradientBoostingRegressor
    m = HistGradientBoostingRegressor(
        max_iter=max_iter, learning_rate=0.08, max_leaf_nodes=31, min_samples_leaf=200,
        l2_regularization=1.0, early_stopping=True, validation_fraction=0.1,
        n_iter_no_change=15, random_state=seed)
    m.fit(X[:, cols], y)
    return m


def clusters(C, k=40, seed=SEED):
    from sklearn.cluster import KMeans
    U = C / np.maximum(np.linalg.norm(C, axis=1, keepdims=True), 1e-12)
    return KMeans(n_clusters=min(k, len(U)), n_init=4, random_state=seed).fit_predict(U)


def _pds_contrib(P, Ptruth):
    A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    B = Ptruth / np.maximum(np.linalg.norm(Ptruth, axis=1, keepdims=True), 1e-300)
    Sm = A @ B.T
    own = np.diag(Sm).copy()
    beaten = (Sm < own[:, None]).sum(1)
    ties = (Sm == own[:, None]).sum(1) - 1
    return (beaten + 0.5 * ties) / (Sm.shape[0] - 1)


def boot_pds(Ppred, base_contrib, Ptruth, lab, n=400, seed=SEED):
    """Target-cluster bootstrap on the paired per-target pds contribution."""
    d = _pds_contrib(Ppred, Ptruth) - base_contrib
    rng = np.random.default_rng(seed)
    cl = np.unique(lab)
    idx = {c: np.where(lab == c)[0] for c in cl}
    draws = np.empty(n)
    for i in range(n):
        pick = rng.choice(cl, len(cl), replace=True)
        draws[i] = np.mean(np.concatenate([d[idx[c]] for c in pick]))
    return float(d.mean()), float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


# ------------------------------------------------------------------------- folds
def known_answers(S, L, Ce, alpha, b, gate, Etruth):
    """Rule 10's corollary: rows whose answer is known before the run.

    A screen with no known-answer row cannot detect a whole class of bug, because the wrong
    numbers look exactly as reasonable as the right ones.  Each of these CAN fail.
    """
    rows = []
    mass = np.abs(np.expm1(b[None, :] + Etruth).sum(1) / 5.0e4 - 1.0).max()
    rt = np.abs(render(Etruth, b) - Etruth).max()
    print(f"  known-answer  {'truth is on the mass manifold':<38} "
          f"max|mass-1| {mass:.2e}   render(E)==E to {rt:.2e}")
    # Tolerances are float64 round-off through expm1/log1p summed over ~5-7k genes, not a
    # tuned gate: a wrong b_ctrl, a wrong gene axis or a broken sub-frame moves these to
    # 1e-2 or worse, i.e. six orders above the bar.
    assert mass < 1e-8 and rt < 1e-6, "the surface truth is not renderer-invariant"

    # A delivered ZERO effect is zero only to float64 round-off (~1e-12), and that round-off
    # is the SAME vector for every target, so unit-normalising it yields one arbitrary
    # direction shared by all rows -> a uniformly random ranking, mean 0.5, sd sqrt(1/(12 n)).
    # The check is therefore "at chance", not "exactly 0.5" -- and it still fails loudly if
    # anything truth-dependent leaks into deliver().
    # The MEANINGFUL assertion is on the delivered zero effect's NORM, not its direction: a
    # ~1e-15 vector has an arbitrary direction, so its cosine to the truth is noise (core4
    # read -0.0264, core5 +0.0016, both from a ratio of ~2e-15).  Asserting on the norm ratio
    # still fails loudly if anything leaks into deliver(); asserting on the cosine was
    # measuring round-off.
    z = deliver(np.zeros_like(Ce), b, gate)
    zr = float(np.linalg.norm(z) / np.linalg.norm(Etruth))
    print(f"  known-answer  {'delivered zero effect is zero':<38} "
          f"||deliver(0)||/||truth|| {zr:.2e}")
    assert zr < 1e-9, f"deliver() of a zero effect is not zero ({zr:.2e})"

    sd0 = float(np.sqrt(1.0 / (12.0 * S.nT)))
    checks = [
        ("zero effect (chance)", lambda: endpoints(z, Etruth), 0.5, 5 * sd0, None, None, True),
        ("ORACLE truth, ungated, alpha=1", lambda: endpoints(render(Etruth, b), Etruth),
         1.0, 1e-6, 1.0, 1e-6, True),
        ("ORACLE truth through deliver()",
         lambda: endpoints(deliver(alpha * Etruth, b, gate), Etruth), 1.0, 2e-3, None, None, True),
        ("ORACLE residual R_L (this expt's ceiling)",
         lambda: endpoints(deliver(alpha * S.R[L], b, gate), Etruth), 1.0, 5e-2, None, None, False),
    ]
    for nm, fn, ep, tp, ec, tc, hard in checks:
        e = fn()
        bad = (ep is not None and abs(e["pds"] - ep) > tp) or \
              (ec is not None and abs(e["cos"] - ec) > tc)
        print(f"  known-answer  {nm:<38} pds {e['pds']:.5f}  cos {e['cos']:.5f}"
              f"{'  <-- FAIL' if bad else '  ok'}")
        if bad and hard:
            raise AssertionError(f"known-answer row {nm!r}: pds {e['pds']}, cos {e['cos']}")
        rows.append(dict(line=L, regime="known", fset=nm, w=0.0, **e))
    return rows


def run(surface_name, regimes, n_train, fset_spec, out_path, n_tgroups=3, wsweep=None,
        max_iter=150, use_ctrl=True, base_power=1.0):
    S = Surface(surface_name, ctrl=use_ctrl)
    NAMES, FAMILY = feature_names(bool(S.ctrl))
    rng = np.random.default_rng(SEED)
    wsweep = wsweep or [0.5, 1.0, 2.0, 4.0, 8.0, 16.0]
    boot_w = {0.5, 1.0, 4.0}          # bootstrap + null only where an arm might be shipped
    print(f"surface {surface_name}: {S.nT} targets x {S.nG} genes, lines "
          f"{', '.join(S.lines)}  (tables: {', '.join(S.members)})")
    print(f"  {len(NAMES)} features: "
          + ", ".join(f"{f}={sum(1 for n in NAMES if FAMILY[n] == f)}" for f in "TCX"))

    idx = {n: i for i, n in enumerate(NAMES)}
    fsets = {}
    for tag in fset_spec:
        if tag == "T":
            cols = [i for i, n in enumerate(NAMES) if FAMILY[n] == "T"]
        elif tag == "C":
            cols = [i for i, n in enumerate(NAMES) if FAMILY[n] == "C"]
        elif tag == "X":
            cols = [i for i, n in enumerate(NAMES) if FAMILY[n] == "X"]
        elif tag == "T+C":
            cols = [i for i, n in enumerate(NAMES) if FAMILY[n] in "TC"]
        elif tag == "T+X":
            cols = [i for i, n in enumerate(NAMES) if FAMILY[n] in "TX"]
        elif tag == "ALL":
            cols = list(range(len(NAMES)))
        elif tag == "ALL-div":
            cols = [i for i, n in enumerate(NAMES) if n != "div"]
        elif tag == "ALL-cell":                 # everything except the cell-level context
            cols = [i for i, n in enumerate(NAMES) if n not in CTX_CELL]
        else:
            raise ValueError(tag)
        if cols:
            fsets[tag] = cols
    print(f"  feature sets: " + "; ".join(f"{k}({len(v)})" for k, v in fsets.items()))

    tg = rng.integers(0, n_tgroups, S.nT)          # target groups for the LOTO axis
    results = []

    for L in S.lines:
        pool = [x for x in S.lines if x != L]
        t0 = time.time()
        Xe, ye, Ce = features(S, L, pool, NAMES)
        b, gate, Etruth = S.bL[L], S.gateL[L], S.EL[L]

        # alpha: match the consensus row norm to the donors' median row norm.  Deployment-
        # available, ONE constant per fold, IDENTICAL across arms.  Every endpoint is
        # scale-invariant, so alpha enters only through deliver()'s nonlinearity.
        dn = np.median([np.median(np.linalg.norm(S.R[d], axis=1)) for d in pool])
        alpha = dn / np.median(np.linalg.norm(Ce, axis=1))

        # `base_power` replaces the comparison baseline with the context-free per-entry
        # re-shaping sign(C)|C|^p (`recal.py`).  The FEATURES still come from the raw C,
        # which is what deployment has.  This asks the only question that decides the
        # verdict: does the learned residual add anything ON TOP of the trivial sharpening
        # that already beats it 3x?
        if base_power != 1.0:
            Cb = np.sign(Ce) * np.abs(Ce) ** base_power
            Cb = restore_norms(Cb, Ce)
        else:
            Cb = Ce
        base_eff = alpha * Cb
        Pbase = deliver(base_eff, b, gate)
        base = endpoints(Pbase, Etruth)
        Ptruth = Etruth
        lab = clusters(Ce)
        base_contrib = _pds_contrib(Pbase, Ptruth)     # reused by every bootstrap

        print(f"\n{'=' * 118}\nHELD LINE {L}   donors {pool}   alpha {alpha:.4f}")
        print(f"  CONSENSUS BASELINE (equal weights, delivered): pds {base['pds']:.5f}  "
              f"cos {base['cos']:.5f}  row-row {base['rr']:+.5f}  sse* {base['sse']:.6f}")
        results.append(dict(line=L, regime="baseline", fset="-", w=0.0, **base))
        results += known_answers(S, L, Ce, alpha, b, gate, Etruth)

        for regime in regimes:
            groups = [None] if regime == "R1" else list(range(n_tgroups))
            for g in groups:
                # ---------- training rows.  R1/R3 never see L; R2 trains on L's OTHER targets.
                train_lines = pool if regime in ("R1", "R3") else S.lines
                pairs = []
                for B in train_lines:
                    others = [x for x in train_lines if x != B]
                    if len(others) < 2:
                        continue
                    pairs.append((B, others))
                    if len(others) > 2:            # donor-count augmentation: eval uses
                        sub = list(rng.choice(others, 2, replace=False))   # |pool| donors,
                        pairs.append((B, sub))     # training recipients only |pool|-1
                per = max(n_train // max(len(pairs), 1), 1)
                Xs, ys = [], []
                for B, dn_B in pairs:
                    Xb, yb, _ = features(S, B, dn_B, NAMES)
                    keep = np.ones(S.nT, dtype=bool)
                    if g is not None:
                        keep &= (tg != g)          # held targets never appear in training
                    ridx = np.where(np.repeat(keep, S.nG))[0]
                    take = rng.choice(ridx, min(per, len(ridx)), replace=False)
                    Xs.append(Xb[take]); ys.append(yb[take])
                    del Xb, yb
                Xtr = np.concatenate(Xs); ytr = np.concatenate(ys)
                del Xs, ys

                emask = np.ones(S.nT, dtype=bool) if g is None else (tg == g)
                bs = base if g is None else endpoints(Pbase[emask], Ptruth[emask])
                for fs_name, cols in fsets.items():
                    if regime != "R1" and fs_name not in ("ALL", "T"):
                        continue
                    m = fit(Xtr, ytr, cols, max_iter=max_iter)
                    yh = m.predict(Xe[:, cols]).reshape(S.nT, S.nG).astype(np.float64)
                    yh_c = yh - yh.mean(0, keepdims=True)      # mandatory: 160's trap
                    par, orth = split_par_orth(yh_c, Cb)
                    aln = float(np.mean((yh_c * Cb).sum(1) / np.maximum(
                        np.linalg.norm(yh_c, axis=1) * np.linalg.norm(Cb, axis=1), 1e-300)))
                    r2 = 1.0 - float(np.mean((ye.reshape(S.nT, S.nG) - yh) ** 2)
                                     / np.mean(ye ** 2))
                    frac = float(np.median(np.linalg.norm(orth, axis=1)
                                           / np.linalg.norm(Cb, axis=1)))
                    # how much of the residual the model COULD have found is orthogonal:
                    _, u_true = split_par_orth(ye.reshape(S.nT, S.nG).astype(np.float64), Cb)
                    ucos = float(np.mean((orth * u_true).sum(1) / np.maximum(
                        np.linalg.norm(orth, axis=1) * np.linalg.norm(u_true, axis=1), 1e-300)))
                    print(f"  {regime}{'' if g is None else f'/g{g}'} {fs_name:<9} "
                          f"[R^2(y) {r2:+.4f}  cos(yhat,C) {aln:+.3f}  |orth|/|C| {frac:.3f}"
                          f"  cos(orth, true orth) {ucos:+.4f}]")
                    for w in wsweep:
                        for cen, tagc in ((orth, "orth"), (par, "par-CONTROL")):
                            if tagc.startswith("par") and w != 1.0:
                                continue
                            eff = restore_norms(Cb + w * cen, Cb) * alpha
                            P = deliver(eff, b, gate)
                            e = (endpoints(P[emask], Ptruth[emask]) if g is not None
                                 else endpoints(P, Ptruth))
                            row = dict(line=L, regime=regime, tgroup=g, fset=fs_name, w=w,
                                       arm=tagc, eval_r2=r2, cos_yhat_C=aln, relnorm=frac,
                                       cos_orth_true=ucos, **e,
                                       d_pds=e["pds"] - bs["pds"], d_cos=e["cos"] - bs["cos"],
                                       base_pds=bs["pds"], base_cos=bs["cos"])
                            if g is None and tagc == "orth" and w in boot_w:
                                mu, lo, hi = boot_pds(P, base_contrib, Ptruth, lab)
                                row.update(boot_lo=lo, boot_hi=hi, boot_mean=mu)
                                # PERMUTATION NULL: the same orthogonal residual rows dealt to
                                # the WRONG targets.  Now that the parallel (inert, but
                                # sign-flipping) component is gone, this null reads ~0 if the
                                # prediction is not target-specific and collapses only if it is.
                                pm = rng.permutation(S.nT)
                                en = endpoints(deliver(restore_norms(Cb + w * cen[pm], Cb)
                                                       * alpha, b, gate), Ptruth)
                                row["null_d_pds"] = en["pds"] - bs["pds"]
                                row["null_d_cos"] = en["cos"] - bs["cos"]
                            results.append(row)
                            print(f"      w={w:<5} {tagc:<12} pds {e['pds']:.5f} "
                                  f"({row['d_pds']:+.5f})  cos {e['cos']:.5f} "
                                  f"({row['d_cos']:+.5f})  rr {e['rr']:+.5f}"
                                  + (f"  boot[{row['boot_lo']:+.5f},{row['boot_hi']:+.5f}]"
                                     f"  null {row['null_d_pds']:+.5f}"
                                     if "boot_lo" in row else ""))
                del Xtr, ytr
        del Xe, ye
        print(f"  fold {L} done in {time.time() - t0:.0f}s", flush=True)

    with open(out_path, "w") as fh:
        json.dump(dict(surface=surface_name, seed=SEED, rows=results,
                       features=NAMES, family=FAMILY), fh, indent=1, default=float)
    print(f"\nwrote {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--surface", default="core5")
    ap.add_argument("--regimes", default="R1")
    ap.add_argument("--fsets", default="T,C,X,T+C,T+X,ALL,ALL-div,ALL-cell")
    ap.add_argument("--n-train", type=int, default=1_500_000)
    ap.add_argument("--max-iter", type=int, default=150)
    ap.add_argument("--no-ctrl", action="store_true")
    ap.add_argument("--base-power", type=float, default=1.0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out = a.out or os.path.join(WORK, f"residual-{a.surface}-{a.regimes.replace(',', '')}.json")
    run(a.surface, a.regimes.split(","), a.n_train, a.fsets.split(","), out,
        max_iter=a.max_iter, use_ctrl=not a.no_ctrl, base_power=a.base_power)


if __name__ == "__main__":
    main()
