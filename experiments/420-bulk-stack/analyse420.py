"""420 analysis: every PLOLO bulk design -> board prediction vs 417a's bulk (PP + 410 key clip(0.3,3) + sharpening q 1.25).

Two mse energy routes (418's cliff-aware model, 187's per-context identity, X solved from 410b's board raw):
  d418   418's route exactly: A_c/B_c = 410b's MEASURED x (DESIGN render energy of the arm / of PP|K3) on the real panel;
         pds/cos from the DESIGNED rows. Reproduces 418's table (KA-A below) and is kept for continuity.
  r417a  PRIMARY: A_c/B_c = 417a's MEASURED artifact energies x (ballast-REALISED render energy of the arm / of 417a's
         design, KMAX 16) on the real panel; pds/cos from the REALISED rows. Validated on the four measured artifacts
         (abc420: 410b/395, 410c/410b, 417a/410b), where it is 2-5x closer than the design route (see VALIDATION).
         The KMAX 24 / 32 columns use the same route with the arm realised at that KMAX (the reference stays 417a, K 16).
X_ref = X_410b (1 + k_key rel_ref)/(1 + k_key rel_PP|K3) sqrt(mean A_ref / mean A_410b); the arm's own direction change r
(cos_agg vs the reference) enters as X_arm = X_ref (1 + k r) sqrt(A ratio). k = 1 and 1.407 (the key class) are reported.
overall = (2.2138 dpds_bd + (-1.051) dmse_raw) / 6 with dpds_bd the board-emulated pds (lam 0.5).
Re-levelling (d): per-context EV-optimal bulk scale on the arm's own A_c/B_c/X under 187's cliff sd (ctxlevel417's
model), the builder factor read off the real-panel realised energy map (abc420 --level); its pds from comp420 --level.

    uv run python analyse420.py            # -> work/analyse420.{json,log}
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
W418 = ROOT / "experiments/418-bulk-composition/work"
sys.path.insert(0, str(ROOT / "experiments/418-bulk-composition"))
import analyse418 as A8  # noqa: E402  (adds 187/096/src to the path; Cliff asserts X_410b == analyse414)

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
CTX = ("A", "B", "C")
K_KEY = A8.K_KEY
KS = (1.0, K_KEY)
REF = "PP|K3+sh1.25"
KOLF2_REL, KOLF2_PDS = 0.371, 0.332   # 414: kolf2 / cancer-mean sharpening transfer at q 1.25 (rel +0.0401 vs +0.1083; pds +0.0063 vs +0.0190)
SD187 = {"A": 0.0025, "B": 0.0004, "C": 0.0013}      # ctxlevel417 / 414 section 6: 187's cliff position sd (log scale)
BASE_ENV = ("VCC_NORMREL_EXP=-1.5 VCC_NORMREL_MODE=cap VCC_CD4_PARTIAL_G=0.5 VCC_DE_SIGN_PRIOR=clean4q80 "
            "VCC_PERT_SCALE_BALLAST=1")          # 417a's env apart from the three bulk knobs below
LEVEL_DESIGNS = ("PP|K3+sh1.25", "PP|K3+sh1.35", "PP|K3+sh1.5", "PP|K3+sh1.75", "PP|K3+sh2", "PP|V1_5_3+sh1.25",
                 "PP|V1_5_3+sh1.5", "PP|V1_5_3+sh1.75", "stack|V1_5_3+sh1.25", "stack|V1_5_3+sh1.5", "stack|V1_5_3+sh1.75",
                 "stack|K3+sh1.25", "stack|K3+sh1.5")
sys.path.insert(0, str(HERE))
from lib420 import COMP_ENV, KEY_ENV  # noqa: E402


def load_folds(tag="comp420"):
    RS = {}
    for f in (*CANCER, "kolf2"):
        p = WORK / f"{tag}-{f}.json"
        if p.exists():
            d = json.loads(p.read_text())
            if "peak_rss_gib" in d:
                RS[f] = d["arms"]
    return RS


class Pricer:
    def __init__(self):
        self.cl = A8.Cliff()
        a = json.loads((WORK / "abc420.json").read_text())
        self.abc, self.meas = a["designs"], a["measured"]
        self.none = {f: json.loads((W418 / f"comp418-{f}.json").read_text())["arms"]["PP|none"]["cos_agg"]
                     for f in (*CANCER, "kolf2") if (W418 / f"comp418-{f}.json").exists()}
        self.A410 = float(np.mean([self.cl.dec["410b"][c]["A"] for c in CTX]))
        lv = WORK / "abc420_level.json"
        self.level = json.loads(lv.read_text()) if lv.exists() else {}

    @staticmethod
    def abc_name(name):
        if name.startswith("CTRL:"):
            return "PP|V1+sh1.25"
        base = name.split(":")[0]           # level variants "design:A=..,B=..,C=.."
        return base.replace("V1own", "V1")  # own-key arms: energy of the fixed-key design (418's convention)

    def E(self, name, src):
        d = self.abc[self.abc_name(name)]
        return d["design"] if src == "design" else d["realised"][src[1:]]

    def price(self, R, f, arm, k, route="r417a", kmax=16, ref=REF, k_pds=1.0, s=None, E_arm=None, s_ref=None):
        """Arm vs ref on fold f. Returns d_rel, dpds_bd, dmse, ov, A ratio vs ref per context, rho per context, X."""
        cl = self.cl
        base = self.none[f]
        if route == "d418":
            sa, sr, sK = R[arm]["design"], R[ref]["design"], R["PP|K3"]["design"]
            Ea, Er, E0, M0 = self.E(arm, "design"), self.E(ref, "design"), self.E("PP|K3", "design"), cl.dec["410b"]
        else:
            sa, sr, sK = R[arm]["realised"][str(kmax)], R[ref]["realised"]["16"], R["PP|K3"]["realised"]["16"]
            Ea, Er = self.E(arm, f"R{kmax}"), self.E(ref, "R16")
            E0, M0 = self.E(REF, "R16"), self.meas["417a"]     # anchor: 417a's design realised == 417a's measured energy
        if E_arm is not None:
            Ea = E_arm
        AB = lambda E: {c: dict(A=M0[c]["A"] * E[c]["A"] / E0[c]["A"], B=M0[c]["B"] * E[c]["B"] / E0[c]["B"]) for c in CTX}
        ABa, ABr = AB(Ea), AB(Er)
        rel_ref, relK3 = sr["cos_agg"] / base - 1, sK["cos_agg"] / base - 1
        Xr = cl.X410b * (1 + K_KEY * rel_ref) / (1 + K_KEY * relK3) * np.sqrt(np.mean([ABr[c]["A"] for c in CTX]) / self.A410)
        rr = cl.raw(ABr, Xr, s_ref)
        r = sa["cos_agg"] / sr["cos_agg"] - 1
        amp = np.sqrt(np.mean([ABa[c]["A"] for c in CTX]) / np.mean([ABr[c]["A"] for c in CTX]))
        Xa = Xr * (1 + k * r) * amp
        ra = cl.raw(ABa, Xa, s)
        dpb = k_pds * (sa["pds_lam0.5"] - sr["pds_lam0.5"])
        dm = ra - rr
        return dict(d_rel=r, dpds_bd=dpb, dpds_lit=sa["pds"] - sr["pds"], dmse=dm, ov=(A8.PDS_SC * dpb + A8.MSE_SC * dm) / 6,
                    Ar={c: ABa[c]["A"] / ABr[c]["A"] for c in CTX}, rho={c: ABa[c]["B"] / cl.claim[c] for c in CTX},
                    X=Xa, ABa=ABa)

    def level_opt(self, ABa, X, design, rng_seed=420):
        """EV-optimal per-context builder factor f_c for one arm (its own A_c/B_c/X), under 187's cliff sd. Uses the
        measured realised energy map A(f)/A(1), B(f)/B(1) of `design` on the real panel (abc420 --level); falls back to
        the render-linear map (s^2) when the design has no map. Returns {c: f}, {c: d raw EV}, {c: d raw nominal}."""
        cl = self.cl
        mp = self.level.get(design)
        grid = np.round(np.arange(0.90, 1.2001, 0.001), 3)   # the builder accepts factors in [0.8, 1.2]
        if mp is not None:
            fs = sorted(float(x) for x in mp)
            assert max(fs) >= grid[-1] - 1e-9 and min(fs) <= grid[0] + 0.04, f"level map of {design} does not span the grid"
            grid = grid[grid >= min(fs)]
            ra = {c: np.interp(grid, fs, [mp[str(x) if str(x) in mp else f"{x:g}"][c]["realised"]["A"] /
                                          mp["1.0"][c]["realised"]["A"] for x in fs]) for c in CTX}
            rb = {c: np.interp(grid, fs, [mp[str(x) if str(x) in mp else f"{x:g}"][c]["realised"]["B"] /
                                          mp["1.0"][c]["realised"]["B"] for x in fs]) for c in CTX}
        else:
            ra = rb = {c: grid ** 2 for c in CTX}
        rng = np.random.default_rng(rng_seed)
        best, dev, dnom = {}, {}, {}
        for c in CTX:
            A, B, Kc, claim = ABa[c]["A"], ABa[c]["B"], cl.Kc[c], cl.claim[c]
            draws = claim * np.exp(2 * rng.normal(0, SD187[c], 400))
            base = (A - np.minimum(draws, B) - 2 * X)
            ev = np.array([np.mean((A * ra[c][i] - np.minimum(draws, B * rb[c][i]) - 2 * X * np.sqrt(ra[c][i])) - base)
                           for i in range(grid.size)]) / cl.sumD
            nom = np.array([(A * ra[c][i] - min(claim, B * rb[c][i]) - 2 * X * np.sqrt(ra[c][i])) - (A - min(claim, B) - 2 * X)
                            for i in range(grid.size)]) / cl.sumD
            i = int(np.argmin(ev))
            best[c], dev[c], dnom[c] = float(grid[i]), float(ev[i]), float(nom[i])
        return best, dev, dnom


def fmt_fold(v):
    return " / ".join(f"{x:+.4f}" for x in v)


def main():
    P = Pricer()
    RS = load_folds()
    folds = [f for f in (*CANCER, "kolf2") if f in RS]
    cancer = [f for f in CANCER if f in RS]
    L = []
    out = dict(folds=folds, ka={}, validation={}, rows={})
    say = lambda *a: (print(*a), L.append(" ".join(str(x) for x in a)))

    # ---------------------------------------------------------------- KA-A: 418's table through my pricer (d418 route)
    a418 = json.loads((W418 / "analyse418.json").read_text())["rows"]
    worst = 0.0
    for mine, theirs in (("h1x2|K3+sh1.25", "b:h1x2"), ("cd4x1.5|K3+sh1.25", "a:cd4x1.5")):
        for f in folds:
            for k, key in ((1.0, "abc_k1"), (K_KEY, "abc_k1.41")):
                got = P.price(RS[f], f, mine, k, route="d418")["ov"]
                exp = a418[theirs]["K3+sh"][f][key]["ov"]
                worst = max(worst, abs(got - exp))
    say(f"KA-A: 418's A/B/C-energy prices of H1 x2 and CD4 x1.5 vs 417a reproduced through my pricer: max |d| {worst:.1e}")
    assert worst < 1e-12, "KA-A FAILED"
    out["ka"]["analysis_vs_418"] = worst

    # ---------------------------------------------------------------- VALIDATION of the energy routes (real panel)
    M = P.meas
    say("\nVALIDATION: render-energy ratios on the real A/B/C panel, predicted vs MEASURED on the shipped artifacts (A/B/C)")
    rows = (("410b/395", ("PP|K3", "PP|none"), ("410b", "395")), ("410c/410b", ("PP|Kinf", "PP|K3"), ("410c", "410b")),
            ("417a/410b", ("PP|K3+sh1.25", "PP|K3"), ("417a", "410b")))
    cc = json.loads((ROOT / "experiments/414-bulk-reprice/work/cliffcheck414.json").read_text())
    B417 = "my-thesis-advisor-is-a-cat.h5ad"          # 417b (q 1.5) = submissions/417b-...-sharpen150; measured after abc420 ran
    if B417 in cc:
        M = dict(M, **{"417b": {c: dict(A=cc[B417]["per_context"][c]["A"], B=cc[B417]["per_context"][c]["B"]) for c in CTX}})
        rows = rows + (("417b/417a", ("PP|K3+sh1.5", "PP|K3+sh1.25"), ("417b", "417a")),
                       ("417b/410b", ("PP|K3+sh1.5", "PP|K3"), ("417b", "410b")))
        out["measured_417b"] = M["417b"]
    for lab, (a, b), (ma, mb) in rows:
        meas = [M[ma][c]["A"] / M[mb][c]["A"] for c in CTX]
        des = [P.abc[a]["design"][c]["A"] / P.abc[b]["design"][c]["A"] for c in CTX]
        rb = "design" if b == "PP|none" else "R16"         # 395 has no ballast: its delivered bulk is the design
        rea = [P.abc[a]["realised"]["16"][c]["A"] / (P.abc[b]["design"][c]["A"] if rb == "design"
                                                     else P.abc[b]["realised"]["16"][c]["A"]) for c in CTX]
        ed = [d - m for d, m in zip(des, meas)]
        er = [d - m for d, m in zip(rea, meas)]
        out["validation"][lab] = dict(measured=meas, design=des, realised=rea)
        say(f"  {lab:10s} measured {fmt_fold(meas)} | design {fmt_fold(des)} (err {fmt_fold(ed)}) | realised K16 "
            f"{fmt_fold(rea)} (err {fmt_fold(er)})")

    # ---------------------------------------------------------------- rows
    names = [n for n in RS[folds[0]] if "realised" in RS[folds[0]][n]]
    for n in names:
        row = {}
        for route in ("r417a", "d418"):
            for k in KS:
                row[f"{route}_k{k:.3g}"] = {f: P.price(RS[f], f, n, k, route) for f in folds}
        for km in (24, 32):
            row[f"r417a_K{km}_k1"] = {f: P.price(RS[f], f, n, 1.0, "r417a", kmax=km) for f in folds}
            row[f"r417a_K{km}_k1.41"] = {f: P.price(RS[f], f, n, K_KEY, "r417a", kmax=km) for f in folds}
        out["rows"][n] = row

    def m(n, key, fld="ov", fl=None):
        v = [out["rows"][n][key][f][fld] for f in (fl or cancer)]
        return float(np.mean(v)), sum(x > 0 for x in v), v

    def energy_ratio(n):
        """realised K16 A/B/C render energy of the arm vs 417a's (the r417a route's ratio)."""
        return [P.E(n, "R16")[c]["A"] / P.E(REF, "R16")[c]["A"] for c in CTX]

    def builder(n):
        comp, rest = n.split("|")
        key, q = (rest.split("+sh") + [None])[:2]
        parts = [BASE_ENV]
        if key == "V1own":
            parts.append("VCC_PERT_SCALE=<NEW artifact: LODO-V1 recomputed on this composition>")
        else:
            parts.append(KEY_ENV[key])
        parts.append(f"VCC_BULK_SHARPEN_Q={q}" if q else "(no VCC_BULK_SHARPEN_Q)")
        if COMP_ENV[comp]:
            parts.append(f"VCC_BULK_LOAD=<dump made with {COMP_ENV[comp]} + the same VCC_PERT_SCALE>")
        return " ".join(parts[1:])

    # ---------------------------------------------------------------- the full ranked table
    say(f"\n==== ALL DESIGNS vs 417a ({REF}); PRIMARY route r417a (realised, K16); cancer folds k562/hct116/hek293t; "
        "kolf2 no vote")
    hdr = (f"{'design':26s} {'d rel':>7s} {'dpds_bd':>8s} {'dmse k1':>8s} {'dmse kK':>8s} | {'ov k=1 per fold':>26s} "
           f"{'mean':>8s} n+ | {'ov k=1.407 per fold':>26s} {'mean':>8s} n+ | {'d418 k1':>8s} {'d418 kK':>8s} | "
           f"{'E ratio vs 417a':>20s} flag | kolf2 k1/kK")
    say(hdr)
    ranked = sorted(names, key=lambda n: -m(n, "r417a_k1")[0])
    out["ranked"] = ranked
    for n in ranked:
        k1, kk = m(n, "r417a_k1"), m(n, "r417a_k1.41")
        er = energy_ratio(n)
        flag = "HI" if max(er) > 1.03 else ("LO" if min(er) < 0.90 else "")
        mm = lambda fld, key="r417a_k1": float(np.mean([out["rows"][n][key][f][fld] for f in cancer]))
        kf = (f"{out['rows'][n]['r417a_k1']['kolf2']['ov']:+.4f}/{out['rows'][n]['r417a_k1.41']['kolf2']['ov']:+.4f}"
              if "kolf2" in folds else "-")
        say(f"{n:26s} {mm('d_rel'):+7.4f} {mm('dpds_bd'):+8.5f} {mm('dmse'):+8.4f} {mm('dmse', 'r417a_k1.41'):+8.4f} | "
            f"{fmt_fold(k1[2]):>26s} {k1[0]:+8.5f} {k1[1]}/{len(cancer)} | {fmt_fold(kk[2]):>26s} {kk[0]:+8.5f} {kk[1]}/{len(cancer)} | "
            f"{m(n, 'd418_k1')[0]:+8.5f} {m(n, 'd418_k1.41')[0]:+8.5f} | {'/'.join(f'{x:.3f}' for x in er):>20s} {flag:4s} | {kf}")

    # ---------------------------------------------------------------- KMAX columns for the sharpening doses
    say("\n==== (b) sharpening dose on PP|K3 and PP|V1: ballast realisability and the KMAX 24/32 emulation (r417a, cancer mean)")
    say(f"{'design':22s} {'capped K16':>10s} {'cos(real,des) q1':>16s} {'ov K16 k1/kK':>18s} {'ov K24 k1/kK':>18s} "
        f"{'ov K32 k1/kK':>18s} {'capped K24/K32':>14s} {'E ratio K16 vs 417a':>22s}")
    for key in ("K3", "V1"):
        for q in ("1.25", "1.35", "1.5", "1.75", "2"):
            n = f"PP|{key}+sh{q}"
            if n not in out["rows"]:
                continue
            b = P.abc[n]["ballast"]
            ov = lambda kk: m(n, kk)[0]
            say(f"{n:22s} {100 * b['16']['frac_capped']:9.1f}% {b['16']['cos_realised_vs_design_q'][0]:16.4f} "
                f"{ov('r417a_k1'):+8.5f}/{ov('r417a_k1.41'):+8.5f} {ov('r417a_K24_k1'):+8.5f}/{ov('r417a_K24_k1.41'):+8.5f} "
                f"{ov('r417a_K32_k1'):+8.5f}/{ov('r417a_K32_k1.41'):+8.5f} {100 * b['24']['frac_capped']:6.1f}/{100 * b['32']['frac_capped']:4.1f}% "
                f"{'/'.join(f'{x:.3f}' for x in energy_ratio(n)):>22s}")

    out["builder"] = {n: builder(n) for n in names if not n.startswith("CTRL")}
    H = json.loads((W418 / "analyse418.json").read_text())["hist_constants"]["abc"]
    out["kolf2_like_sharp"] = dict(k=KOLF2_REL, k_pds=KOLF2_PDS)

    def step(arm, ref, cls, f):
        """One step arm-vs-ref on fold f under three transfer readings: k=1, k=1.407, and the class-calibrated one
        (composition: 418's history constants; key: k_key 1.407; sharpening: kolf2-like 0.371 / 0.332; level: 1)."""
        e = {"k1": P.price(RS[f], f, arm, 1.0, ref=ref), "kK": P.price(RS[f], f, arm, K_KEY, ref=ref)}
        if cls in ("a", "b"):
            e["cls"] = P.price(RS[f], f, arm, H[cls]["k_mse"], ref=ref, k_pds=H[cls]["k_pds"])
        elif cls == "key":
            e["cls"] = e["kK"]
        elif cls == "sharp":
            e["cls"] = P.price(RS[f], f, arm, KOLF2_REL, ref=ref, k_pds=KOLF2_PDS)
        else:
            e["cls"] = e["k1"]
        return e

    def boot(f, arm, ref, n=2000, seed=420):
        Z = np.load(WORK / f"arrays420-{f}.npz")
        d = Z["d"]
        rng = np.random.default_rng(seed)
        g = lambda nm, fld: Z[f"{nm}@R16__{fld}"]
        xa, aa, pa, xr, ar, pr = g(arm, "x"), g(arm, "a"), g(arm, "pdsl"), g(ref, "x"), g(ref, "a"), g(ref, "pdsl")
        agg = lambda x, a, i: np.mean([abs(x[c, i].sum()) / np.sqrt(a[c, i].sum() * d[i].sum()) for c in range(3)])
        rel, dp = [], []
        for _ in range(n):
            i = rng.integers(0, d.size, d.size)
            rel.append(agg(xa, aa, i) / agg(xr, ar, i) - 1)
            dp.append(pa[:, i].mean() - pr[:, i].mean())
        return dict(p_rel=float(np.mean(np.asarray(rel) > 0)), p_pds=float(np.mean(np.asarray(dp) > 0)),
                    pds_ci=np.quantile(dp, [0.025, 0.975]).tolist())

    def show_step(lab, arm, ref, cls, store):
        per = {f: step(arm, ref, cls, f) for f in folds}
        bt = {f: boot(f, arm, ref) for f in cancer}
        mean = lambda key, fld="ov": float(np.mean([per[f][key][fld] for f in cancer]))
        npos = lambda key: sum(per[f][key]["ov"] > 0 for f in cancer)
        er = [P.E(arm, "R16")[c]["A"] / P.E(ref, "R16")[c]["A"] for c in CTX]
        er417 = energy_ratio(arm)
        flag = "HI" if max(er417) > 1.03 else ("LO" if min(er417) < 0.90 else "")
        store[lab] = dict(arm=arm, ref=ref, cls=cls, per=per, boot=bt, energy_vs_ref=er, energy_vs_417a=er417,
                          builder=builder(arm) if not arm.startswith("CTRL") else "-")
        kf = f"{per['kolf2']['k1']['ov']:+.4f}" if "kolf2" in per else "-"
        say(f"{lab:34s} {mean('k1', 'd_rel'):+7.4f} {mean('k1', 'dpds_bd'):+8.5f} {mean('k1', 'dmse'):+8.4f}/{mean('kK', 'dmse'):+8.4f} | "
            f"{fmt_fold([per[f]['k1']['ov'] for f in cancer])} {mean('k1'):+8.5f} {npos('k1')}/{len(cancer)} | "
            f"{fmt_fold([per[f]['kK']['ov'] for f in cancer])} {mean('kK'):+8.5f} {npos('kK')}/{len(cancer)} | cls {mean('cls'):+8.5f} {npos('cls')}/{len(cancer)} | "
            f"{'/'.join(f'{x:.3f}' for x in er417)} {flag:2s} | boot P(rel>0)/P(pds>0) "
            + " ".join(f"{bt[f]['p_rel']:.2f}/{bt[f]['p_pds']:.2f}" for f in cancer) + f" | kolf2 {kf}")

    HDR = (f"{'step':34s} {'d rel':>7s} {'dpds_bd':>8s} {'dmse k1/kK':>17s} | {'ov k=1 per fold':>26s} {'mean':>8s} n+ | "
           f"{'ov k=1.407 per fold':>26s} {'mean':>8s} n+ | {'class-calibrated':>16s} | E vs 417a | bootstrap | kolf2 k1")
    out["one_step"] = {}
    say("\n==== ONE-VARIABLE STEPS FROM 417a (primary route r417a); class: comp = 418 history constants, key = k 1.407, "
        "sharpening = kolf2-like transfer (rel x0.371, pds x0.332)")
    say(HDR)
    ONE = [("q 1.25 -> 1.35", "PP|K3+sh1.35", "sharp"), ("q 1.25 -> 1.5 (= 417b)", "PP|K3+sh1.5", "sharp"),
           ("q 1.25 -> 1.75", "PP|K3+sh1.75", "sharp"), ("q 1.25 -> 2.0", "PP|K3+sh2", "sharp"),
           ("key K3 -> V1 (lodo415_lo0.3_hi3)", "PP|V1+sh1.25", "key"), ("key K3 -> K5_3 (wshare_lo0.5_hi3)", "PP|K5_3+sh1.25", "key"),
           ("key K3 -> V1_5_3 (new artifact)", "PP|V1_5_3+sh1.25", "key"),
           ("comp H1 x2", "h1x2|K3+sh1.25", "b"), ("comp CD4 x1.5", "cd4x1.5|K3+sh1.25", "a"),
           ("comp stack H1x2+CD4x1.5", "stack|K3+sh1.25", "stack"), ("comp stack H1x1.5+CD4x1.5", "stack_h1x1.5|K3+sh1.25", "stack"),
           ("comp stack H1x2+CD4x1.25", "stack_cd4x1.25|K3+sh1.25", "stack"),
           ("CONTROL V1 permuted", "CTRL:PP|V1shuf+sh1.25", "key")]
    for lab, arm, cls in ONE:
        show_step(lab, arm, REF, cls, out["one_step"])
    # the stack's class-calibrated price = the H1 step (class b) then the CD4 step on top of it (class a), summed
    for lab in ("comp stack H1x2+CD4x1.5",):
        st = {f: step("h1x2|K3+sh1.25", REF, "b", f)["cls"]["ov"] + step("stack|K3+sh1.25", "h1x2|K3+sh1.25", "a", f)["cls"]["ov"]
              for f in folds}
        out["one_step"][lab]["cls_path_h1_then_cd4"] = st
        say(f"   stack, class-calibrated along H1 then CD4: {fmt_fold([st[f] for f in cancer])} mean {np.mean([st[f] for f in cancer]):+.5f}")

    # ---------------------------------------------------------------- interactions: does each lever survive a higher q?
    say("\n==== INTERACTIONS: each non-q lever priced vs the SAME q without it (r417a, k=1 / k=1.407 cancer means, n+)")
    out["interactions"] = {}
    for lab, fmt_arm, fmt_ref in (("V1 vs K3", "PP|V1+sh{q}", "PP|K3+sh{q}"), ("K5_3 vs K3", "PP|K5_3+sh{q}", "PP|K3+sh{q}"),
                                  ("V1_5_3 vs K3", "PP|V1_5_3+sh{q}", "PP|K3+sh{q}"), ("V1_5_3 vs V1", "PP|V1_5_3+sh{q}", "PP|V1+sh{q}"),
                                  ("stack vs PP (K3)", "stack|K3+sh{q}", "PP|K3+sh{q}"), ("stack vs PP (V1_5_3)", "stack|V1_5_3+sh{q}", "PP|V1_5_3+sh{q}"),
                                  ("H1x2 vs PP (K3)", "h1x2|K3+sh{q}", "PP|K3+sh{q}"), ("CD4x1.5 vs PP (K3)", "cd4x1.5|K3+sh{q}", "PP|K3+sh{q}"),
                                  ("stack vs H1x2 (K3)", "stack|K3+sh{q}", "h1x2|K3+sh{q}")):
        cells = []
        for q in ("1.25", "1.35", "1.5", "1.75", "2"):
            a, r = fmt_arm.format(q=q), fmt_ref.format(q=q)
            v1 = [P.price(RS[f], f, a, 1.0, ref=r)["ov"] for f in cancer]
            vk = [P.price(RS[f], f, a, K_KEY, ref=r)["ov"] for f in cancer]
            out["interactions"].setdefault(lab, {})[q] = dict(k1=v1, kK=vk)
            cells.append(f"q{q}: {np.mean(v1):+.5f}/{np.mean(vk):+.5f} {sum(x > 0 for x in v1)}/{len(cancer)}")
        say(f"  {lab:22s} " + " | ".join(cells))

    # ---------------------------------------------------------------- stepwise paths (each step one variable from the previous)
    say("\n==== PATHS: each step priced vs the previous arm (r417a; class-calibrated as above); cumulative = sum of steps")
    out["paths"] = {}
    PATHS = {
        "P1 q-first": [("q 1.5", "PP|K3+sh1.5", "sharp"), ("key V1_5_3", "PP|V1_5_3+sh1.5", "key"),
                       ("stack", "stack|V1_5_3+sh1.5", "stack")],
        "P2 key-first": [("key V1_5_3", "PP|V1_5_3+sh1.25", "key"), ("q 1.5", "PP|V1_5_3+sh1.5", "sharp"),
                         ("stack", "stack|V1_5_3+sh1.5", "stack")],
        "P3 key V1": [("key V1", "PP|V1+sh1.25", "key"), ("q 1.5", "PP|V1+sh1.5", "sharp"), ("stack", "stack|V1+sh1.5", "stack")],
        "P4 stack-first": [("stack", "stack|K3+sh1.25", "stack"), ("q 1.5", "stack|K3+sh1.5", "sharp"),
                           ("key V1_5_3", "stack|V1_5_3+sh1.5", "key")],
        "P5 q-first to 1.75": [("q 1.5", "PP|K3+sh1.5", "sharp"), ("q 1.75", "PP|K3+sh1.75", "sharp"),
                               ("key V1_5_3", "PP|V1_5_3+sh1.75", "key")],
    }
    for pname, steps in PATHS.items():
        prev, cum = REF, {"k1": 0.0, "kK": 0.0, "cls": 0.0}
        rec = []
        say(f" {pname}")
        for lab, arm, cls in steps:
            store = {}
            show_step(f"   {lab}", arm, prev, cls if cls != "stack" else "k1", store)
            e = store[f"   {lab}"]
            for kk in cum:
                cum[kk] += float(np.mean([e["per"][f][kk]["ov"] for f in cancer]))
            rec.append(dict(step=lab, arm=arm, ref=prev, cls=cls, **{f"mean_{kk}": float(np.mean([e["per"][f][kk]["ov"] for f in cancer]))
                                                                     for kk in ("k1", "kK", "cls")},
                            n_pos_k1=int(sum(e["per"][f]["k1"]["ov"] > 0 for f in cancer)), builder=e["builder"],
                            energy_vs_417a=e["energy_vs_417a"]))
            prev = arm
        direct = {kk: float(np.mean([P.price(RS[f], f, prev, k, ref=REF)["ov"] for f in cancer]))
                  for kk, k in (("k1", 1.0), ("kK", K_KEY))}
        say(f"   cumulative k1 {cum['k1']:+.5f} kK {cum['kK']:+.5f} cls {cum['cls']:+.5f} (direct end-vs-417a k1 "
            f"{direct['k1']:+.5f} kK {direct['kK']:+.5f})")
        out["paths"][pname] = dict(steps=rec, cumulative=cum, direct=direct)
    # ---------------------------------------------------------------- (d) per-context re-levelling on top
    say("\n==== (d) PER-CONTEXT RE-LEVELLING (VCC_BULK_CTX_SCALE), EV-optimal under 187's cliff sd, on the arm's own "
        "A_c/B_c (r417a) and X; factor read off the real-panel realised energy map")
    out["level"] = {}
    RL = load_folds("comp420L")
    for design in LEVEL_DESIGNS:
        if design not in out["rows"]:
            continue
        for kname, k in (("k1", 1.0), ("kK", K_KEY)):
            pr = {f: P.price(RS[f], f, design, k) for f in folds}
            ABa = pr[cancer[0]]["ABa"]
            Xm = float(np.mean([pr[f]["X"] for f in cancer]))
            fac, dev, dnom = P.level_opt(ABa, Xm, design)
            mp = P.level.get(design)
            if mp is not None:
                fs = sorted(float(x) for x in mp)
                sA = {c: float(np.sqrt(np.interp(fac[c], fs, [mp[f"{x:g}" if f"{x:g}" in mp else str(x)][c]["realised"]["A"] /
                                                              mp["1.0"][c]["realised"]["A"] for x in fs]))) for c in CTX}
            else:
                sA = dict(fac)
            spec = ",".join(f"{c}={fac[c]:.3f}" for c in CTX)
            lv_name = f"{design}:{spec}"
            per = {}
            for f in folds:
                R = dict(RS[f])
                have = f in RL and lv_name in RL[f]
                if have:
                    R[lv_name] = RL[f][lv_name]
                    e = P.price(R, f, lv_name, k, s=sA)
                else:
                    e = P.price(R, f, design, k, s=sA)
                e0 = pr[f]
                per[f] = dict(ov_total=e["ov"], ov_level_only=e["ov"] - e0["ov"], dmse_level=e["dmse"] - e0["dmse"],
                              dpds_level=e["dpds_bd"] - e0["dpds_bd"], pds_measured=have)
            mean = lambda fld: float(np.mean([per[f][fld] for f in cancer]))
            out["level"][f"{design}|{kname}"] = dict(factors=fac, render_scale=sA, ev_draw=dev, nominal=dnom, per=per,
                                                     spec=f"VCC_BULK_CTX_SCALE={spec}")
            say(f"  {design:22s} {kname}: VCC_BULK_CTX_SCALE={spec} (render scale {'/'.join(f'{sA[c]:.4f}' for c in CTX)}) | "
                f"level-only d ov {fmt_fold([per[f]['ov_level_only'] for f in cancer])} mean {mean('ov_level_only'):+.5f} "
                f"(dmse {mean('dmse_level'):+.4f}, dpds {mean('dpds_level'):+.5f}{'' if all(per[f]['pds_measured'] for f in cancer) else ' NOT MEASURED'}) | "
                f"arm+level vs 417a {fmt_fold([per[f]['ov_total'] for f in cancer])} mean {mean('ov_total'):+.5f} "
                f"{sum(per[f]['ov_total'] > 0 for f in cancer)}/{len(cancer)}")
    if "measured_417b" in out:
        say("\n  417b (q 1.5) re-levelling from its MEASURED artifact energy (the factors to build with):")
        for kname, k in (("k1", 1.0), ("kK", K_KEY)):
            pr = {f: P.price(RS[f], f, "PP|K3+sh1.5", k) for f in cancer}
            Xm = float(np.mean([pr[f]["X"] for f in cancer]))
            ABm = {c: dict(A=out["measured_417b"][c]["A"], B=out["measured_417b"][c]["B"]) for c in CTX}
            fac, dev, dnom = P.level_opt(ABm, Xm, "PP|K3+sh1.5")
            say(f"    {kname}: VCC_BULK_CTX_SCALE=" + ",".join(f"{c}={fac[c]:.3f}" for c in CTX) +
                f" | EV d raw per ctx {'/'.join(f'{dev[c]:+.4f}' for c in CTX)} -> mean {np.mean(list(dev.values())):+.4f} raw, "
                f"{A8.MSE_SC * np.mean(list(dev.values())) / 6:+.5f} overall (mse part only) | rho before "
                + "/".join(f"{ABm[c]['B'] / P.cl.claim[c]:.4f}" for c in CTX))
            out.setdefault("level_417b_measured", {})[kname] = dict(factors=fac, ev=dev, nominal=dnom)

    # ---------------------------------------------------------------- recommended sequences, re-levelling as a real step
    say("\n==== SEQUENCES (each arm one variable from the previous; levelling steps use the level variant's measured PLOLO "
        "cos/pds and its factors' render scale; r417a route, cancer folds)")
    lvl = {k[:-3]: v for k, v in out["level"].items() if k.endswith("|k1")}
    RL = load_folds("comp420L")

    def lv(design):
        v = lvl[design]
        return f"{design}:{v['spec'].split('=', 1)[1]}", v["render_scale"], v["spec"]

    def seq_price(f, arm, ref, k, k_pds=1.0):
        R = {**RS[f], **RL.get(f, {})}
        sa = sr = None
        if ":" in arm:
            sa = lvl[arm.split(":")[0]]["render_scale"]
        if ":" in ref:
            sr = lvl[ref.split(":")[0]]["render_scale"]
        return P.price(R, f, arm, k, ref=ref, k_pds=k_pds, s=sa, s_ref=sr)

    SEQS = {
        "STRONG (417a reads k_sh >= 0.75)": [("417b: q 1.5", "PP|K3+sh1.5", "sharp"),
                                             ("+ levelling", lv("PP|K3+sh1.5")[0], "level"),
                                             ("+ stack (LOAD, factors re-derived)", lv("stack|K3+sh1.5")[0], "stack")],
        "WEAK (417a reads k_sh < 0.75)": [("stack (LOAD)", "stack|K3+sh1.25", "stack"),
                                          ("+ key V1_5_3 (new dump)", "stack|V1_5_3+sh1.25", "key"),
                                          ("+ levelling", lv("stack|V1_5_3+sh1.25")[0], "level")],
        "WEAK-alt (existing key only)": [("stack (LOAD)", "stack|K3+sh1.25", "stack"),
                                         ("+ q 1.5", "stack|K3+sh1.5", "sharp"),
                                         ("+ levelling", lv("stack|K3+sh1.5")[0], "level")],
        "LEVEL-FIRST": [("levelling on 417a", lv("PP|K3+sh1.25")[0], "level"),
                        ("+ q 1.5 (factors re-derived)", lv("PP|K3+sh1.5")[0], "sharp"),
                        ("+ stack (LOAD, factors re-derived)", lv("stack|K3+sh1.5")[0], "stack")],
    }
    out["sequences"] = {}
    for sname, steps in SEQS.items():
        say(f" {sname}")
        prev, cum, rows = REF, {"k1": 0.0, "kK": 0.0, "cls": 0.0}, []
        for lab, arm, cls in steps:
            per = {}
            for f in folds:
                e = {"k1": seq_price(f, arm, prev, 1.0), "kK": seq_price(f, arm, prev, K_KEY)}
                if cls == "sharp":
                    e["cls"] = seq_price(f, arm, prev, KOLF2_REL, KOLF2_PDS)
                elif cls == "stack":
                    e["cls"] = e["k1"]
                elif cls == "key":
                    e["cls"] = e["kK"]
                else:
                    e["cls"] = e["k1"]
                per[f] = e
            mean = lambda kk, fld="ov": float(np.mean([per[f][kk][fld] for f in cancer]))
            for kk in cum:
                cum[kk] += mean(kk)
            npos = sum(per[f]["k1"]["ov"] > 0 for f in cancer)
            design = arm.split(":")[0]
            spec = lvl[design]["spec"] if ":" in arm else ""
            b = builder(design) + (f" {spec}" if spec else "")
            rows.append(dict(step=lab, arm=arm, ref=prev, cls=cls, per=per, k1=mean("k1"), kK=mean("kK"), cls_ov=mean("cls"),
                             n_pos_k1=int(npos), d_rel=mean("k1", "d_rel"), dpds_bd=mean("k1", "dpds_bd"),
                             dmse_k1=mean("k1", "dmse"), dmse_kK=mean("kK", "dmse"), builder=b,
                             per_fold_k1=[per[f]["k1"]["ov"] for f in cancer], per_fold_kK=[per[f]["kK"]["ov"] for f in cancer],
                             kolf2_k1=per["kolf2"]["k1"]["ov"] if "kolf2" in per else None))
            say(f"   {lab:36s} rel {mean('k1', 'd_rel'):+.4f} dpds_bd {mean('k1', 'dpds_bd'):+.5f} dmse {mean('k1', 'dmse'):+.4f}/"
                f"{mean('kK', 'dmse'):+.4f} | ov k1 {fmt_fold(rows[-1]['per_fold_k1'])} {mean('k1'):+.5f} {npos}/3 | kK "
                f"{mean('kK'):+.5f} | cls {mean('cls'):+.5f} | cum k1 {cum['k1']:+.5f} kK {cum['kK']:+.5f} cls {cum['cls']:+.5f}")
            say(f"     builder: {b}")
            prev = arm
        out["sequences"][sname] = dict(steps=rows, cumulative=dict(cum))

    # ---------------------------------------------------------------- the 417a read -> sharpening transfer -> q steps
    say("\n==== READING 417a's BOARD ROW: implied sharpening transfer k_sh (mse) and k_pds_sh (pds), and what the q steps are "
        "then worth vs 417a (r417a route, cancer means)")
    cl = P.cl
    AB0 = {c: dict(A=P.meas["410b"][c]["A"], B=P.meas["410b"][c]["B"]) for c in CTX}
    AB1 = {c: dict(A=P.meas["417a"][c]["A"], B=P.meas["417a"][c]["B"]) for c in CTX}
    amp = float(np.sqrt(np.mean([AB1[c]["A"] for c in CTX]) / np.mean([AB0[c]["A"] for c in CTX])))
    r_sh = float(np.mean([RS[f][REF]["realised"]["16"]["cos_agg"] / RS[f]["PP|K3"]["realised"]["16"]["cos_agg"] - 1 for f in cancer]))
    dp_sh = float(np.mean([RS[f][REF]["realised"]["16"]["pds_lam0.5"] - RS[f]["PP|K3"]["realised"]["16"]["pds_lam0.5"] for f in cancer]))
    PDS_410B, MSE_410B = 0.8511557041991825, A8.BOARD_410B_RAW
    say(f"  417a vs 410b on PLOLO (realised, cancer mean): rel cos {r_sh:+.4f}, pds_bd {dp_sh:+.5f}; measured energy amp "
        f"{amp:.4f}; board 410b raw pds {PDS_410B:.6f} mse {MSE_410B:.6f}")
    out["read417a"] = dict(r_sh=r_sh, dp_sh=dp_sh, amp=amp, rows=[])
    say(f"  {'k_sh':>5s} {'pred 417a raw mse':>17s} {'pred raw pds (k_pds=k_sh*.9)':>28s} {'pred d overall vs 410b':>22s} | "
        f"{'q1.35':>8s} {'q1.5':>8s} {'q1.75':>8s} {'q2':>8s} (step vs 417a, ov)")
    for ksh in (0.0, 0.2, KOLF2_REL, 0.55, 0.75, 1.0, 1.2, K_KEY):
        kp = ksh * KOLF2_PDS / KOLF2_REL
        raw1 = cl.raw(AB1, cl.X410b * (1 + ksh * r_sh) * amp)
        pds1 = PDS_410B + kp * dp_sh
        dov = (A8.PDS_SC * (pds1 - PDS_410B) + A8.MSE_SC * (raw1 - MSE_410B)) / 6
        qs = {q: float(np.mean([P.price(RS[f], f, f"PP|K3+sh{q}", ksh, k_pds=kp)["ov"] for f in cancer]))
              for q in ("1.35", "1.5", "1.75", "2")}
        out["read417a"]["rows"].append(dict(k_sh=ksh, k_pds=kp, raw_mse=raw1, raw_pds=pds1, d_overall=dov, q_steps=qs))
        say(f"  {ksh:5.3f} {raw1:17.4f} {pds1:28.4f} {dov:+22.4f} | " + " ".join(f"{qs[q]:+8.5f}" for q in qs))
    # the same read for the RE-LEVELLED q ladder: (PP|K3+shq + its EV-optimal levelling) vs (417a + its levelling)
    say(f"  RE-LEVELLED q ladder (each design at its own EV-optimal VCC_BULK_CTX_SCALE) vs 417a re-levelled, same k_sh grid:")
    out["read417a"]["levelled"] = []
    ref_l = lv("PP|K3+sh1.25")[0]
    for ksh in (0.0, 0.2, KOLF2_REL, 0.55, 0.75, 1.0, 1.2, K_KEY):
        kp = ksh * KOLF2_PDS / KOLF2_REL
        cells = {}
        for q in ("1.35", "1.5", "1.75", "2"):
            d = f"PP|K3+sh{q}"
            if d not in lvl or not all(lv(d)[0] in RL.get(f, {}) for f in cancer):
                continue
            cells[q] = float(np.mean([seq_price(f, lv(d)[0], ref_l, ksh, kp)["ov"] for f in cancer]))
        out["read417a"]["levelled"].append(dict(k_sh=ksh, q_steps=cells))
        say(f"  {ksh:5.3f} " + " ".join(f"q{q} {v:+.5f}" for q, v in cells.items()))
    say("  factors: " + "; ".join(f"q{q}: {lvl[f'PP|K3+sh{q}']['spec']}" for q in ("1.25", "1.35", "1.5", "1.75", "2")
                                   if f"PP|K3+sh{q}" in lvl))
    say("  inverse: k_sh = (solve_X(417a measured A/B, board raw mse) / (X_410b * amp) - 1) / r_sh; "
        "k_pds_sh = (board raw pds - 0.851156) / pds_bd step")
    (WORK / "analyse420.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "analyse420.log").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
