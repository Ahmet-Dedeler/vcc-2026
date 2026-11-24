"""534 Part B pricer: 438's price438.Panel (420's r417a route with the panel as an argument), reading this part's own
comp534 / full534 outputs; the step table; and the pooled paired target bootstrap (PREREGISTRATION-B section 6).

    overall = (2.2138 dpds_bd - 1.051 dmse_raw) / 6                    (420's constants; board-emulated pds at lam 0.5)
    mse_c   = 1 + (A_c + Kc_c - min(claim_c, B_c) - 2 s_c X) / sumD     (187's per-context identity)
    X_ref   = X_410b x ret_f x (1 + k_key rel_ref)/(1 + k_key rel_PP|K3) x sqrt(mean A_ref / A_410b)
    X_arm   = X_ref (1 + k r) sqrt(mean A_arm / mean A_ref),   r = cos_agg(arm)/cos_agg(ref) - 1
LEVELLED (the D/E/F default): each design at its own ntclevel factor (full534), pds/cos re-scored at that factor
(comp534 --level). `Panel.price` is price438.Panel.price verbatim; only the file layout and the control's energy name
(designs534.energy_name) differ. ret_f = cos_agg(h1c1|K3 on this panel's fold f) / the same on A/B/C's fold f, A/B/C's
read from 438's stored abc jsons (read-only).

Bootstrap: targets resampled within each fold (one draw per fold shared by every design: paired), cos_agg / pds
re-aggregated from the per-target arrays, the price recomputed with energies held fixed (438's Panel.boot), the pooled
statistic = the mean over every fold of every panel; 4,000 draws, one Generator seeded 534, consumed panel by panel in
the order given, fold by fold (k562, hct116, hek293t), draw by draw.
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard534  # noqa: E402,F401
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

import designs534 as D  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
WORK = HERE
W438 = ROOT / "experiments/438-october-stack-repricing/work"
sys.path.insert(0, str(ROOT / "experiments/418-bulk-composition"))
import analyse418 as A8  # noqa: E402

SLOTS = ("A", "B", "C")
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
K_KEY = A8.K_KEY


class Panel:
    def __init__(self, tag, cliff=None):
        self.tag = tag
        self.R, self.Z = {}, {}
        for f in CANCER:
            p = WORK / tag / f"comp534-{f}.json"
            if p.exists():
                d = json.loads(p.read_text())
                if "peak_rss_gib" in d:
                    self.R[f] = d["arms"]
                    self.Z[f] = WORK / tag / f"arrays534-{f}.npz"
        self.folds = [f for f in CANCER if f in self.R]
        self.full = json.loads((WORK / tag / "full534.json").read_text())
        self.E = self.full["designs"]
        self.claims = self.full["claims"]
        self.cl = cliff or A8.Cliff()
        self.A410 = float(np.mean([self.cl.dec["410b"][c]["A"] for c in SLOTS]))
        self.ret = {}
        for f in self.R:
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

    def price(self, f, arm, k=1.0, ref=D.M["M1"], level=False, k_pds=1.0, S=None):
        """price438.Panel.price, verbatim (arm vs ref on fold f; S overrides the PLOLO statistics for the bootstrap)."""
        L = ":L" if level else ""
        g = (lambda n: S[n]) if S is not None else (lambda n: self.stats(f, n))
        sa, sr, sK, s0 = g(arm + L), g(ref + L), g("h1c1|K3"), g("h1c1|none")
        en = lambda n: self.E[D.energy_name(n)]
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
            rho = {c: La[c]["B"] / self.claims[c] for c in SLOTS}
        else:
            dm = self.raw(Ea, Xa) - self.raw(Er, Xr)
            rho = {c: Ea[c]["B"] / self.claims[c] for c in SLOTS}
        dpb = k_pds * (sa["pdsl"] - sr["pdsl"])
        return dict(ov=(A8.PDS_SC * dpb + A8.MSE_SC * dm) / 6, d_rel=r, dpds_bd=dpb, dpds_lit=sa["pds"] - sr["pds"],
                    dmse=dm, rho=rho, X=Xa)

    # ------------------------------------------------------------------ per-target arrays (the bootstrap)
    def arrays(self, f):
        if f not in self._arr:
            self._arr[f] = np.load(self.Z[f])
        return self._arr[f]

    def agg_fn(self, f, names):
        """-> (n, agg(i) -> {name: dict(cos, pdsl, pds)}) on target index vector i; identity-checked against the json."""
        Z = self.arrays(f)
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
        full = np.arange(d.size)
        S0 = agg(full)
        for nm in names:                          # KA-BOOT: the identity draw reproduces the json statistics
            ref = self.stats(f, nm)
            assert abs(S0[nm]["cos"] - ref["cos"]) < 1e-12 and abs(S0[nm]["pdsl"] - ref["pdsl"]) < 1e-12, (self.tag, f, nm)
        return d.size, agg


def step_names(step):
    sid, arm, ref, what, cls = step
    return D.M[arm], D.M[ref]


def need_names(level=True):
    L = ":L" if level else ""
    out = {"h1c1|K3", "h1c1|none"}
    for st in D.STEPS:
        a, r = step_names(st)
        out |= {a + L, r + L}
    return sorted(out)


def step_prices(P, f, S=None, level=True, cal=False):
    """{step id: overall} on fold f, levelled; cal=True applies the calibrated transfer (designs534.K_CAL)."""
    out = {}
    for st in D.STEPS:
        a, r = step_names(st)
        k, kp = D.K_CAL[st[4]] if cal else (1.0, 1.0)
        out[st[0]] = P.price(f, a, k, r, level, k_pds=kp, S=S)["ov"]
    return out


def point_table(panels, level=True, cal=False):
    """{panel: {step: dict(per_fold, mean, n_pos)}} plus the pooled mean over every fold."""
    tab, pooled = {}, {st[0]: [] for st in D.STEPS}
    for t, P in panels.items():
        per = {f: step_prices(P, f, level=level, cal=cal) for f in P.folds}
        tab[t] = {}
        for st in D.STEPS:
            v = {f: per[f][st[0]] for f in P.folds}
            tab[t][st[0]] = dict(per_fold=v, mean=float(np.mean(list(v.values()))), n_pos=int(sum(x > 0 for x in v.values())),
                                 n_folds=len(v))
            pooled[st[0]] += list(v.values())
    tab["pooled"] = {s: dict(mean=float(np.mean(v)), n_pos=int(sum(x > 0 for x in v)), n_folds=len(v), per_fold=v)
                     for s, v in pooled.items()}
    return tab


def pooled_bootstrap(panels, n=4000, seed=534, level=True, cal=False):
    """Paired target bootstrap; returns {step: draws of the pooled mean} and {panel: {step: draws of its panel mean}},
    plus 'cum' = the sum of the six milestone steps."""
    rng = np.random.default_rng(seed)
    names = need_names(level)
    per_panel = {}
    all_fold_draws = []
    for t, P in panels.items():
        fold_draws = []
        for f in P.folds:
            nt, agg = P.agg_fn(f, names)
            vals = np.empty((n, len(D.STEPS)))
            for b in range(n):
                i = rng.integers(0, nt, nt)
                S = agg(i)
                pr = step_prices(P, f, S=S, level=level, cal=cal)
                vals[b] = [pr[st[0]] for st in D.STEPS]
            fold_draws.append(vals)
        per_panel[t] = np.mean(fold_draws, axis=0)
        all_fold_draws += fold_draws
    pooled = np.mean(all_fold_draws, axis=0)

    def pack(M):
        o = {st[0]: M[:, j] for j, st in enumerate(D.STEPS)}
        o["cum"] = sum(M[:, j] for j, st in enumerate(D.STEPS) if st in D.MAIN_STEPS)
        return o
    return pack(pooled), {t: pack(v) for t, v in per_panel.items()}


def summ(x):
    x = np.asarray(x)
    return dict(ci95=[float(np.quantile(x, 0.025)), float(np.quantile(x, 0.975))], p_le0=float(np.mean(x <= 0)),
                sd=float(np.std(x)))
