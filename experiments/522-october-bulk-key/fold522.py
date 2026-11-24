"""522 — exact per-target-KEY arms from the dumps, no operator (~1-2 GB), on two kinds of surface:

  Oct(panel, fold)  an October-like panel's PLOLO fold (dump522-<panel>-<fold>.npz; 546o's path: lib438, the October
                    stack h2c1, the fold's own LODO fit): key_scale(fin, w) -> row sharpening q 1.75 on each frame's gate
                    -> render (minnorm) -> ballast realisation at KMAX 16 against the h1c1 carriers -> L.row_score.
                    KNOWN ANSWER on load (abort unless |d| < 1e-12, L.cmp): 546o's `base q1.75`, `fixed 6`, `p96.6` rows.
  Abc(fold)         an A/B/C PLOLO fold through 518's score518.Fold (dump518-<fold>.npz; 540's path: the A/B/C stack,
                    q 1.75, realised 16). Its own KA on load: arms540's V1_0.5_3 and V1_0.5_6 rows.

Both expose: targets, kv (LODO key vector), carried, sources, held, score(w, mask=None) -> (summary, arrays).
`mask` restricts the surface to a SUB-PANEL (key_scale, carriers, truth and pds all over the kept targets only).
The strength feature z (518's "src": centred percentile rank of the robust strength, mean over the fold's own sources
that carry the target, 0 where none) is `zsrc(surface)`; `zsrc(surface, extra=("cd4_robust",))` adds CD4.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/438-october-stack-repricing"))
sys.path.insert(0, str(ROOT / "experiments/518-bulk-mse-gap"))
import json  # noqa: E402

import numpy as np  # noqa: E402

import lib438 as L  # noqa: E402
from lib438 import SLOTS  # noqa: E402
from plolo import TS, profile_from_effect  # noqa: E402

WORK = HERE / "work"
W546 = ROOT / "experiments/546-october-cap-pertarget/work"
TOL = 1e-12
SRC = ("k562_gwps", "orion_hct116", "orion_hek293t")
P_PROD = 96.64697522769474          # 497's unrounded P (== cap 6 on the A/B/C artifact)
P_546 = 96.647                      # the rounded P 543/546 priced (known-answer rows)


def crank(v):
    """518's lofo518.crank: centred percentile rank in [-0.5, 0.5] over the finite entries; 0 elsewhere."""
    v = np.asarray(v, float)
    ok = np.isfinite(v)
    out = np.zeros(v.size)
    if ok.sum() > 1:
        r = np.argsort(np.argsort(v[ok]))
        out[ok] = r / (ok.sum() - 1) - 0.5
    return out


def pct_hi(kv, carried, p, mask=None):
    """543's percentile cap: the p-th percentile (numpy linear) of key / median over the carried targets."""
    c = carried if mask is None else (carried & mask)
    kn = kv[c] / float(np.median(kv[c]))
    return float(np.percentile(kn, p))


def wclip(kv, carried, lo, hi, mask=None):
    """key410.weights (clip(key / median_carried, lo, hi) on carried, 1 elsewhere), median over the kept targets."""
    c = carried if mask is None else (carried & mask)
    med = float(np.median(kv[c]))
    w = np.ones(len(kv))
    w[carried] = np.clip(kv[carried] / med, lo, hi)
    return w


class Oct:
    kind = "oct"

    def __init__(self, panel: str, held: str, check: bool = True):
        self.panel, self.held = panel, held
        z = np.load(WORK / f"dump522-{panel}-{held}.npz")
        self.targets = [str(t) for t in z["targets"]]
        self.sources = [str(s) for s in z["sources"]]
        self.fin, self.fPP = z["fin"], z["fPP"]
        self.kv, self.carried = z["kv"], z["carried"]
        self.avail = z["avail"]
        self.V, self.dnorm, self.pnorm, self.cosPd = z["V"], z["dnorm"], z["pnorm"], z["cosPd"]
        self.cos_dd, self.Fnorm = z["cos_dd"], z["Fnorm"]
        self.scored, self.pp_pos = z["scored"], z["panel_positions"]
        self.truth_full = np.asarray(z["truth"], np.float64)
        self.truth = self.truth_full[:, self.scored]
        self.betas = {c: z[f"beta_{c}"] for c in SLOTS}
        self.frames = {c: (z[f"b_{c}"], z[f"gate_{c}"]) for c in SLOTS}
        self.nonpanel = np.ones(self.fin.shape[1], bool)
        self.nonpanel[self.pp_pos] = False
        self.n = len(self.targets)
        self.P_DE = {c: [profile_from_effect(self.fPP[i], b, g, "minnorm") for i in range(self.n)]
                     for c, (b, g) in self.frames.items()}
        self.d = (self.truth * self.truth).sum(1)
        if check:
            ref = json.loads((W546 / f"comp546o-{panel}-{held}.json").read_text())["arms"]
            hi = pct_hi(self.kv, self.carried, P_546)
            self.ka = {}
            for tag, w in (("base q1.75", wclip(self.kv, self.carried, 0.5, 3.0)),
                           ("fixed 6", wclip(self.kv, self.carried, 0.5, 6.0)),
                           ("p96.6", wclip(self.kv, self.carried, 0.5, hi))):
                s, _ = self.score(w)
                dd = L.cmp(s, ref[tag]["realised"]["16"])
                self.ka[tag] = dd
                assert dd < TOL, f"KA FAILED {panel}/{held} {tag}: {dd}"

    def score(self, w, mask=None, q: float = 1.75):
        idx = np.arange(self.n) if mask is None else np.flatnonzero(mask)
        bulk = L.key_scale(self.fin[idx], np.asarray(w, float)[idx])[0]
        rea, K = {}, []
        for c, (b, gate) in self.frames.items():
            F = L.sharpen(bulk, float(q), gate, mode="row")
            R = np.empty_like(F)
            for r, i in enumerate(idx):
                p_bk = profile_from_effect(F[r], b, gate, "minnorm")
                p_re, _, k, _, _ = L.realise_k(self.P_DE[c][i], p_bk, 16)
                R[r] = np.log1p(TS * p_re) - b
                K.append(k)
            rea[c] = R
        s, arr = L.row_score(rea, self.scored, self.truth[idx], self.betas, self.nonpanel)
        s["K_mean"] = float(np.mean(K))
        arr["d"] = self.d[idx]
        return s, arr


class Abc:
    kind = "abc"

    def __init__(self, held: str, check: bool = True):
        import plolo
        # score518.Fold reads the challenge frames through arms414.frames() -> plolo.challenge_bctrl, which caches under
        # plolo.WORK. lib438's import points that at an empty dir; point it at 420's read-only cache (the files exist, so
        # nothing is written anywhere).
        plolo.WORK = ROOT / "experiments/420-bulk-stack/work/plolo"
        for c in "ABC":
            assert (plolo.WORK / f"bctrl-challenge-{c}.npy").exists()
        from score518 import Fold
        self.F = Fold(held, check=check)
        self.panel, self.held = "abc", held
        self.targets = self.F.targets
        self.sources = [s for s in SRC if s != held]
        self.kv, self.carried = self.F.kv1, self.F.carried
        self.n = self.F.n
        self.d = (self.F.truth * self.F.truth).sum(1)
        self.ka = {"arms540 V1_0.5_3 / V1_0.5_6 (score518.Fold)": "PASS" if check else "unchecked"}

    def score(self, w, mask=None, q: float = 1.75):
        F = self.F
        if mask is None:
            s, arr = F.score(w=np.asarray(w, float), q=q)
            arr["d"] = self.d
            return s, arr
        # sub-panel: score518.Fold's rows() on the kept targets (its arithmetic, restricted)
        from lib420 import TS as TS4, key_scale, profile_from_effect as pfe, realise_k, sharpen
        from comp420 import score_rows
        idx = np.flatnonzero(mask)
        bulk = key_scale(F.fST[idx], np.asarray(w, float)[idx])[0]
        rea, K = {}, []
        for c, (b, g) in F.fr.items():
            X = sharpen(bulk, q, g, mode="row")
            R = np.empty_like(X)
            for r, i in enumerate(idx):
                p_bk = pfe(X[r], b, g, "minnorm")
                p_re, _, k, _, _ = realise_k(F.P_DE[c][i], p_bk, 16)
                R[r] = np.log1p(TS4 * p_re) - b
                K.append(k)
            rea[c] = R
        s, arr = score_rows(rea, F.scored, F.truth[idx], F.betas, F.nonpanel)
        s["K_mean"] = float(np.mean(K))
        arr["d"] = self.d[idx]
        return s, arr


def zsrc(S, extra=(), strength=None):
    """518's src feature for surface S: mean over the fold's own sources (+ `extra` screens) of crank(strength over
    S's targets), NaN-skipped per target, 0 where no screen carries it. On A/B/C this is skey518.zfeat(F, 'src')."""
    st = np.load(WORK / "strength522.npz") if strength is None else strength
    panel = [str(t) for t in st[f"{S.panel}__targets"]]
    idx = np.asarray([panel.index(t) for t in S.targets])
    cols = []
    for k in list(S.sources) + list(extra):
        v = np.asarray(st[f"{S.panel}__{k}"], float)[idx]
        r = crank(v)
        r[~np.isfinite(v)] = np.nan
        cols.append(r)
    with np.errstate(all="ignore"):
        m = np.nanmean(np.vstack(cols), 0)
    return np.where(np.isfinite(m), m, 0.0)


def surfaces(which=("abc", "efp", "efr", "real"), check=True):
    out = []
    for p in which:
        if p == "abc":
            for f in ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"):
                out.append((p, f))
        else:
            for f in SRC:
                out.append((p, f))
    return out


def load(panel, held, check=True):
    return Abc(held, check) if panel == "abc" else Oct(panel, held, check)
