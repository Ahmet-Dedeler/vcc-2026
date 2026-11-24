"""438 pricer: 420's r417a route (analyse420.Pricer.price) with the panel as an argument.

    overall = (2.2138 dpds_bd - 1.051 dmse_raw) / 6                    (420's constants; board-emulated pds at lam 0.5)
    mse_c   = 1 + (A_c + Kc_c - min(claim_c, B_c) - 2 s_c X) / sumD     (187's per-context identity)
    A_c/B_c = 417a's MEASURED energy x (this panel's realised energy of the arm / A/B/C's realised energy of 417a)
    X_ref   = X_410b x ret_f x (1 + k_key rel_ref)/(1 + k_key rel_PP|K3) x sqrt(mean A_ref / A_410b)
    X_arm   = X_ref (1 + k r) sqrt(mean A_arm / mean A_ref),   r = cos_agg(arm)/cos_agg(ref) - 1
Panel dependence: claims (187's on A/B/C; ntclevel's f = 1 claims from the bundle's controls elsewhere), energies (the
panel's own full-panel render), and ret_f = cos_agg(PP|K3 on this panel's fold f) / cos_agg(PP|K3 on A/B/C's fold f):
the board's alignment term X scales with the PLOLO transfer cosine of the same design on the new panel (233/238's
retention assumption; ret = 1 on A/B/C, and `--no-ret` prices with ret = 1 everywhere as a sensitivity).
LEVELLED prices (the D/E/F default: the re-level is mandatory, 432): each design at ntclevel's factor, its pds/cos
re-scored on PLOLO at that factor (comp438 --level), its energies re-rendered at that factor (full438).

KA-P: on A/B/C this pricer reproduces analyse420.json's r417a one-step prices (k = 1 and 1.407, four folds) exactly.

    uv run python price438.py              -> work/price438.json, work/price438.md
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
W420 = ROOT / "experiments/420-bulk-stack/work"
sys.path.insert(0, str(ROOT / "experiments/418-bulk-composition"))
import analyse418 as A8  # noqa: E402

SLOTS = ("A", "B", "C")
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
K_KEY = A8.K_KEY
REF = "h1c1|K3+sh1.25"
NBOOT = 2000


def october_union():
    """The union of the four essential screens' perturbation lists (238 / 436: n 2,389). Arc's A/B/C panel carries 0 of
    them, so an Arc-like October panel should carry ~0 too (436: point estimate 0, 95% upper bound 3)."""
    out = set()
    for p in ("experiments/166-donor-acquisition/donors/nadig_hepg2", "data/donors/nadig_jurkat", "data/donors/rpe1",
              "data/donors/k562_essential"):
        out |= {l.strip() for l in (ROOT / p / "perts.txt").read_text().splitlines() if l.strip()}
    return out


def fold_energies(tag, fold_ops, claims):
    """The FOLD-ROUTE energy surrogate (418's 'PLOLO-energy' column): each design's realised render energy averaged over the
    three cancer folds, anchored at 417a's measured energy through A/B/C's own fold-route 417a energy. Used only where a
    panel has no 3-donor operator; validated against the full-panel route on abc/real/efp (surr438)."""
    a420 = json.loads((W420 / "abc420.json").read_text())
    M0 = {c: a420["measured"]["417a"][c] for c in SLOTS}
    ld = lambda t: {f: json.loads((WORK / t / f"comp438-{f}-{fold_ops if t == tag else 'held'}.json").read_text())["arms"]
                    for f in CANCER}
    R, R0 = ld(tag), ld("abc")
    mE = lambda RR, n, k: {c: {x: float(np.mean([RR[f][n]["realised"]["16"]["AB"][c][x] for f in CANCER])) for x in ("A", "B")}
                           for c in SLOTS}
    E0 = mE(R0, REF, "16")
    out = {}
    for n in R[CANCER[0]]:
        if n.endswith(":L") or "#" in n or not R[CANCER[0]][n].get("realised"):
            continue
        e = mE(R, n, "16")
        d = dict(AB={c: {x: M0[c][x] * e[c][x] / E0[c][x] for x in ("A", "B")} for c in SLOTS})
        if all(f"{n}:L" in R[f] for f in CANCER):
            el = mE(R, f"{n}:L", "16")
            d["AB_level"] = {c: {x: M0[c][x] * el[c][x] / E0[c][x] for x in ("A", "B")} for c in SLOTS}
            d["rho_level"] = {c: d["AB_level"][c]["B"] / claims[c] for c in SLOTS}
        d["rho"] = {c: d["AB"][c]["B"] / claims[c] for c in SLOTS}
        out[n] = d
    return out


class Panel:
    def __init__(self, tag, fold_ops="held", full_tag=None, ret=True, cliff=None, exclude=None, energy="full"):
        """exclude: target symbols dropped from every fold's PLOLO statistics (re-aggregated from the per-target arrays;
        energies are the full panel's)."""
        self.tag, self.fold_ops = tag, fold_ops
        self.R, self.Z = {}, {}
        for f in (*CANCER, "kolf2"):
            p = WORK / tag / f"comp438-{f}-{fold_ops}.json"
            if p.exists():
                d = json.loads(p.read_text())
                if "peak_rss_gib" in d:
                    self.R[f] = d["arms"]
                    self.Z[f] = WORK / tag / f"arrays438-{f}-{fold_ops}.npz"
        self.folds = [f for f in CANCER if f in self.R]
        self.keep, self._st = {}, {}
        for f in self.R:
            self.keep[f] = None
            if exclude:
                tg = [str(t) for t in np.load(WORK / tag / f"fold-{f}-{fold_ops}.npz", allow_pickle=True)["targets"]]
                m = np.asarray([t not in exclude for t in tg])
                self.keep[f] = m if not m.all() else None
        self.full = json.loads((WORK / (full_tag or tag) / "full438.json").read_text())
        self.E = self.full["designs"]
        self.claims = self.full["claims"]
        if energy == "fold":
            self.E = fold_energies(tag, fold_ops, self.claims)
        self.cl = cliff or A8.Cliff()
        self.A410 = float(np.mean([self.cl.dec["410b"][c]["A"] for c in SLOTS]))
        self.ret = {}
        for f in self.R:
            abc = json.loads((WORK / "abc" / f"comp438-{f}-held.json").read_text())["arms"]
            self.ret[f] = (self.R[f]["h1c1|K3"]["realised"]["16"]["cos_agg"] / abc["h1c1|K3"]["realised"]["16"]["cos_agg"]
                           if ret else 1.0)
        self._boot = {}

    def raw(self, AB, X, sx=None):
        sx = sx or {c: 1.0 for c in SLOTS}
        return float(np.mean([A8.mse_c(AB[c]["A"], AB[c]["B"], self.cl.Kc[c], self.claims[c], sx[c] * X, self.cl.sumD)
                              for c in SLOTS]))

    def stats(self, f, name):
        if self.keep.get(f) is not None:
            if (f, name) not in self._st:
                Z = self._arrays(f)
                t = f"{name}@R16" if f"{name}@R16__x" in Z.files else f"{name}@D"
                x, a, pl, pp, d = Z[f"{t}__x"], Z[f"{t}__a"], Z[f"{t}__pdsl"], Z[f"{t}__pds"], Z["d"]
                i = np.flatnonzero(self.keep[f])
                self._st[(f, name)] = dict(cos=float(np.mean([abs(x[c, i].sum()) / np.sqrt(a[c, i].sum() * d[i].sum())
                                                             for c in range(3)])),
                                           pdsl=float(pl[:, i].mean()), pds=float(pp[:, i].mean()))
            return self._st[(f, name)]
        e = self.R[f][name]
        s = e["realised"]["16"] if e.get("realised") else e["design"]
        return dict(cos=s["cos_agg"], pdsl=s["pds_lam0.5"], pds=s["pds"])

    def price(self, f, arm, k=1.0, ref=REF, level=False, k_pds=1.0, S=None):
        """arm vs ref on fold f. S overrides the PLOLO statistics (the bootstrap)."""
        L = ":L" if level else ""
        g = (lambda n: S[n]) if S is not None else (lambda n: self.stats(f, n))
        sa, sr, sK, s0 = g(arm + L), g(ref + L), g("h1c1|K3"), g("h1c1|none")
        en = lambda n: self.E["h1c1|V1+sh1.25" if n.startswith("CTRL:") else n]   # 420's convention for the control
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

    # ------------------------------------------------------------------ bootstrap over targets, per fold
    def _arrays(self, f):
        if f not in self._boot:
            self._boot[f] = np.load(self.Z[f])
        return self._boot[f]

    def boot(self, arm, ref=REF, level=False, k=1.0, n=NBOOT, seed=438):
        """CI of the cancer-fold MEAN overall: targets resampled within each fold (the same draw for every design of
        that fold), cos_agg and pds re-aggregated, the price recomputed with energies held fixed."""
        L = ":L" if level else ""
        names = {arm + L, ref + L, "h1c1|K3", "h1c1|none"}
        rng = np.random.default_rng(seed)
        per = {}
        for f in self.folds:
            Z = self._arrays(f)
            d = Z["d"]
            tag = {nm: (f"{nm}@R16" if f"{nm}@R16__x" in Z.files else f"{nm}@D") for nm in names}
            X = {nm: (Z[f"{tag[nm]}__x"], Z[f"{tag[nm]}__a"], Z[f"{tag[nm]}__pdsl"], Z[f"{tag[nm]}__pds"]) for nm in names}

            def agg(nm, i):
                x, a, pl, pp = X[nm]
                cos = float(np.mean([abs(x[c, i].sum()) / np.sqrt(a[c, i].sum() * d[i].sum()) for c in range(3)]))
                return dict(cos=cos, pdsl=float(pl[:, i].mean()), pds=float(pp[:, i].mean()))
            full = np.arange(d.size) if self.keep.get(f) is None else np.flatnonzero(self.keep[f])
            S0 = {nm: agg(nm, full) for nm in names}
            for nm in names:            # the identity draw must reproduce the json's statistics (KA-BOOT)
                ref_s = self.stats(f, nm)
                assert abs(S0[nm]["cos"] - ref_s["cos"]) < 1e-12 and abs(S0[nm]["pdsl"] - ref_s["pdsl"]) < 1e-12, (f, nm)
            vals = []
            for _ in range(n):
                i = full[rng.integers(0, full.size, full.size)]
                vals.append(self.price(f, arm, k, ref, level, S={nm: agg(nm, i) for nm in names})["ov"])
            per[f] = np.asarray(vals)
        mean = np.mean([per[f] for f in self.folds], axis=0)
        return dict(ci=np.quantile(mean, [0.025, 0.975]).tolist(), p_pos=float(np.mean(mean > 0)),
                    per_fold_p_pos={f: float(np.mean(per[f] > 0)) for f in self.folds})

    def row(self, arm, ref=REF, level=False, boot=True):
        out = dict(arm=arm, ref=ref, level=level)
        for kk, k in (("k1", 1.0), ("kK", K_KEY)):
            p = {f: self.price(f, arm, k, ref, level) for f in self.folds}
            out[kk] = dict(per_fold={f: p[f]["ov"] for f in self.folds}, mean=float(np.mean([p[f]["ov"] for f in self.folds])),
                           n_pos=int(sum(p[f]["ov"] > 0 for f in self.folds)),
                           d_rel=float(np.mean([p[f]["d_rel"] for f in self.folds])),
                           dpds_bd=float(np.mean([p[f]["dpds_bd"] for f in self.folds])),
                           dpds_lit=float(np.mean([p[f]["dpds_lit"] for f in self.folds])),
                           dmse=float(np.mean([p[f]["dmse"] for f in self.folds])))
            out["rho"] = p[self.folds[0]]["rho"]
        if "kolf2" in self.R:
            out["kolf2_k1"] = self.price("kolf2", arm, 1.0, ref, level)["ov"] if (arm + (":L" if level else "")) in self.R["kolf2"] else None
        if boot:
            out["boot"] = self.boot(arm, ref, level)
        return out


def ka_p():
    """KA-P: my pricer on my A/B/C comp jsons + abc energies == analyse420.json's r417a one-step prices."""
    P = Panel("abc")
    a420 = json.loads((W420 / "analyse420.json").read_text())["rows"]
    worst = 0.0
    for mine, theirs in (("h2c1|K3+sh1.25", "h1x2|K3+sh1.25"), ("h1c1.5|K3+sh1.25", "cd4x1.5|K3+sh1.25"),
                         ("h2c1.5|K3+sh1.25", "stack|K3+sh1.25"), ("h1c1|K3+sh1.5", "PP|K3+sh1.5"),
                         ("h1c1|V1_5_3+sh1.25", "PP|V1_5_3+sh1.25"), ("h2c1.5|V1_5_3+sh1.5", "stack|V1_5_3+sh1.5")):
        for f in P.R:
            for kk, k in (("r417a_k1", 1.0), ("r417a_k1.41", K_KEY)):
                got = P.price(f, mine, k)["ov"]
                exp = a420[theirs][kk][f]["ov"]
                worst = max(worst, abs(got - exp))
    assert worst < 1e-12, f"KA-P FAILED: {worst}"
    return worst


if __name__ == "__main__":
    print(f"KA-P: analyse420 r417a one-step prices reproduced, max |d| {ka_p():.1e}")
