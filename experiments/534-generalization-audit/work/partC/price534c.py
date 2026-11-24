"""534 Part C pricer: 438's price438.Panel (levelled route) on Part C's comp534c / full534c (or fold-route) outputs, the
step values Delta_s = B0 - B_s (price(arm = B0, ref = B_s), levelled, k = 1), and the paired target bootstrap.

Energies: the panel's full-panel route (full534c.json) on abc / real / efp; the fold route on efr (foldlevel534c's
fold_energies over the final, levelled fold jsons), exactly as 438 priced efr. The X anchors (h1c1|K3, h1c1|none) and
ret_f are the production base's (designs534c). Bootstrap (PREREGISTRATION-C section 4.1): targets resampled within each
fold, one draw per fold shared by every design (paired); the pooled statistic is the mean over the October-like panels'
nine folds; 4,000 draws; one Generator seeded 534, consumed panel by panel (real, efp, efr), fold by fold.
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard534  # noqa: E402,F401
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

import designs534c as D  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
W438 = ROOT / "experiments/438-october-stack-repricing/work"
sys.path.insert(0, str(ROOT / "experiments/418-bulk-composition"))
import analyse418 as A8  # noqa: E402

SLOTS = ("A", "B", "C")
K_KEY = A8.K_KEY


class Panel:
    def __init__(self, tag, cliff=None, sfx=""):
        self.tag, self.sfx = tag, sfx
        self.R, self.Z = {}, {}
        for f in D.CANCER:
            d = json.loads((HERE / tag / f"comp534c{sfx}-{f}.json").read_text())
            assert "peak_rss_gib" in d and d.get("level_file"), (tag, f, "not a finished levelled pass")
            self.R[f] = d["arms"]
            self.Z[f] = HERE / tag / f"arrays534c{sfx}-{f}.npz"
        self.folds = list(D.CANCER)
        if tag in D.FOLD_ROUTE:
            import foldlevel534c
            self.E = foldlevel534c.fold_energies(tag, sfx)
            self.claims = json.loads((HERE / tag / f"level534c{sfx}.json").read_text())["claims"]
            self.route = "fold"
        else:
            full = json.loads((HERE / tag / f"full534c{sfx}.json").read_text())
            self.E, self.claims, self.route = full["designs"], full["claims"], "full"
        self.cl = cliff or A8.Cliff()
        self.A410 = float(np.mean([self.cl.dec["410b"][c]["A"] for c in SLOTS]))
        self.ret = {}
        for f in self.folds:
            abc = json.loads((W438 / "abc" / f"comp438-{f}-held.json").read_text())["arms"]
            self.ret[f] = self.R[f]["h1c1|K3"]["realised"]["16"]["cos_agg"] / abc["h1c1|K3"]["realised"]["16"]["cos_agg"]
        self._arr = {}

    def raw(self, AB, X, sx=None):
        sx = sx or {c: 1.0 for c in SLOTS}
        return float(np.mean([A8.mse_c(AB[c]["A"], AB[c]["B"], self.cl.Kc[c], self.claims[c], sx[c] * X, self.cl.sumD)
                              for c in SLOTS]))

    def stats(self, f, name):
        e = self.R[f][name]
        s = e["realised"]["16"] if e.get("realised") else e["design"]
        return dict(cos=s["cos_agg"], pdsl=s["pds_lam0.5"], pds=s["pds"])

    def price(self, f, arm, k=1.0, ref="B0", level=True, k_pds=1.0, S=None):
        """price438.Panel.price, verbatim (arm vs ref on fold f; S overrides the PLOLO statistics for the bootstrap)."""
        L = ":L" if level else ""
        g = (lambda n: S[n]) if S is not None else (lambda n: self.stats(f, n))
        sa, sr, sK, s0 = g(arm + L), g(ref + L), g("h1c1|K3"), g("h1c1|none")
        en = lambda n: self.E[n]
        Ea, Er = en(arm)["AB"], en(ref)["AB"]
        rel_ref, relK3 = sr["cos"] / s0["cos"] - 1, sK["cos"] / s0["cos"] - 1
        Xr = self.cl.X410b * self.ret[f] * (1 + K_KEY * rel_ref) / (1 + K_KEY * relK3) * \
            np.sqrt(np.mean([Er[c]["A"] for c in SLOTS]) / self.A410)
        r = sa["cos"] / sr["cos"] - 1
        amp = np.sqrt(np.mean([Ea[c]["A"] for c in SLOTS]) / np.mean([Er[c]["A"] for c in SLOTS]))
        Xa = Xr * (1 + k * r) * amp
        if level:
            La, Lr = en(arm)["AB_level"], en(ref)["AB_level"]
            sxa = {c: np.sqrt(La[c]["A"] / Ea[c]["A"]) for c in SLOTS}
            sxr = {c: np.sqrt(Lr[c]["A"] / Er[c]["A"]) for c in SLOTS}
            dm = self.raw(La, Xa, sxa) - self.raw(Lr, Xr, sxr)
        else:
            dm = self.raw(Ea, Xa) - self.raw(Er, Xr)
        dpb = k_pds * (sa["pdsl"] - sr["pdsl"])
        return dict(ov=(A8.PDS_SC * dpb + A8.MSE_SC * dm) / 6, d_rel=r, dpds_bd=dpb, dpds_lit=sa["pds"] - sr["pds"], dmse=dm)

    def agg_fn(self, f, names):
        if f not in self._arr:
            self._arr[f] = np.load(self.Z[f])
        Z = self._arr[f]
        d = Z["d"]
        X = {}
        for nm in names:
            t = f"{nm}@R16" if f"{nm}@R16__x" in Z.files else f"{nm}@D"
            X[nm] = (Z[f"{t}__x"], Z[f"{t}__a"], Z[f"{t}__pdsl"], Z[f"{t}__pds"])

        def agg(i):
            out = {}
            di = d[i].sum()
            for nm, (x, a, pl, pp) in X.items():
                cos = float(np.mean([abs(x[c, i].sum()) / np.sqrt(a[c, i].sum() * di) for c in range(3)]))
                out[nm] = dict(cos=cos, pdsl=float(pl[:, i].mean()), pds=float(pp[:, i].mean()))
            return out
        S0 = agg(np.arange(d.size))
        for nm in names:                                  # KA-BOOT: the identity draw reproduces the json statistics
            ref = self.stats(f, nm)
            assert abs(S0[nm]["cos"] - ref["cos"]) < 1e-12 and abs(S0[nm]["pdsl"] - ref["pdsl"]) < 1e-12, (self.tag, f, nm)
        return d.size, agg


# secondary (--pk): the toggled base with the PRODUCTION base's LODO key (what a builder-only October knob would deliver)
STEPS_PK = tuple((f"{sid}_pk", arm, f"{ref}pk", *rest) for sid, arm, ref, *rest in D.STEPS)


def need(steps):
    return sorted({"h1c1|K3", "h1c1|none", *[f"{s[1]}:L" for s in steps], *[f"{s[2]}:L" for s in steps]})


def step_prices(P, f, S=None, level=True, k=1.0, steps=D.STEPS):
    return {sid: P.price(f, arm, k, ref, level, S=S)["ov"] for sid, arm, ref, *_ in steps}


def point(P, level=True, steps=D.STEPS):
    per = {f: step_prices(P, f, level=level, steps=steps) for f in P.folds}
    return {sid: dict(per_fold={f: per[f][sid] for f in P.folds}, mean=float(np.mean([per[f][sid] for f in P.folds])),
                      n_pos=int(sum(per[f][sid] > 0 for f in P.folds))) for sid, *_ in steps}


def bootstrap(panels, n=4000, seed=534, steps=D.STEPS):
    """Paired target bootstrap over the given panels' folds (in order); returns {step: draws of the mean over folds}."""
    rng = np.random.default_rng(seed)
    fold_draws = []
    names = need(steps)
    for P in panels:
        for f in P.folds:
            nt, agg = P.agg_fn(f, names)
            v = np.empty((n, len(steps)))
            for b in range(n):
                S = agg(rng.integers(0, nt, nt))
                pr = step_prices(P, f, S=S, steps=steps)
                v[b] = [pr[s[0]] for s in steps]
            fold_draws.append(v)
    M = np.mean(fold_draws, axis=0)
    return {s[0]: M[:, j] for j, s in enumerate(steps)}
