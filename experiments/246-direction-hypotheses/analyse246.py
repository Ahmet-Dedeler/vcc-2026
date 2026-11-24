"""246 — aggregate the per-fold sweeps: cancer-fold means (`kolf2` apart), every falsifier,
concentration (leave-top-3-out on cosine and `pds` separately), and the six-member price
through 218 §7.4's rate card (`price218.price`, uncalibrated k = 0 as the headline).

    uv run python analyse246.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))
from price218 import GATE_SHAPE, GATE_AMPL, price  # noqa: E402

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
FOLDS = CANCER + ("kolf2",)
BASE = "STACK"


def load():
    rows, pt, meta = {}, {}, {}
    for held in FOLDS:
        for tag in ("sweep246", "sweep246b", "sweep246c"):
            p = WORK / f"{tag}-{held}.json"
            if not p.exists():
                continue
            d = json.loads(p.read_text())
            for r in d["rows"]:
                rows[(held, r["arm"])] = r
            z = np.load(WORK / f"pertarget246{tag[len('sweep246'):]}-{held}.npz")
            for k in z.files:
                pt[(held, k)] = z[k]
            if tag == "sweep246":
                meta[held] = d["meta"]
    return rows, pt, meta


def rel(r, b, key):
    return r[key] / b[key] - 1.0


def summarise(rows, arms, folds):
    out = {}
    for arm in arms:
        sel = [(f, rows[(f, arm)], rows[(f, BASE)]) for f in folds if (f, arm) in rows]
        if not sel:
            continue
        agg = [rel(r, b, "cos_agg") for _, r, b in sel]
        mac = [rel(r, b, "a23_cos") for _, r, b in sel]
        dp = [r["pds"] - b["pds"] for _, r, b in sel]
        dq = [r["a15_sgn_w"] - b["a15_sgn_w"] for _, r, b in sel]
        out[arm] = dict(n=len(sel), rel_agg=float(np.mean(agg)), rel_macro=float(np.mean(mac)),
                        d_pds=float(np.mean(dp)), d_qw=float(np.mean(dq)),
                        pds_pos=int(sum(x > 0 for x in dp)), cos_pos=int(sum(x > 0 for x in mac)),
                        per_fold={f: dict(rel_agg=a, rel_macro=m, d_pds=p, d_qw=q)
                                  for (f, _, _), a, m, p, q in zip(sel, agg, mac, dp, dq)},
                        audit=sel[0][1].get("audit", {}))
    return out


def concentration(pt, arm, folds):
    """Leave-top-k-out retention, cosine and pds SEPARATELY (each on its own ordering), plus
    222's version (cosine retention on the pds ordering)."""
    res = {}
    for k in (1, 3, 5):
        rp, rc, rc222 = [], [], []
        for f in folds:
            if (f, f"{arm}|pds") not in pt:
                continue
            dp = pt[(f, f"{arm}|pds")] - pt[(f, f"{BASE}|pds")]
            dc = pt[(f, f"{arm}|cos")] - pt[(f, f"{BASE}|cos")]
            op = np.argsort(-dp)[:k]
            oc = np.argsort(-dc)[:k]
            keep_p = np.ones(dp.size, bool); keep_p[op] = False
            keep_c = np.ones(dc.size, bool); keep_c[oc] = False
            rp.append((dp[keep_p].mean(), dp.mean()))
            rc.append((dc[keep_c].mean(), dc.mean()))
            rc222.append((dc[keep_p].mean(), dc.mean()))
        if not rp:
            continue
        full_p = float(np.mean([t for _, t in rp])); loo_p = float(np.mean([l for l, _ in rp]))
        full_c = float(np.mean([t for _, t in rc])); loo_c = float(np.mean([l for l, _ in rc]))
        loo_c222 = float(np.mean([l for l, _ in rc222]))
        res[k] = dict(pds_full=full_p, pds_loo=loo_p,
                      pds_retained=(loo_p / full_p) if full_p else float("nan"),
                      cos_full=full_c, cos_loo=loo_c,
                      cos_retained=(loo_c / full_c) if full_c else float("nan"),
                      cos_loo_on_pds_order=loo_c222)
    return res


def noncarrying(pt, meta, arm, folds, key):
    """rel change of the mean per-target cosine on rows the residual does NOT carry."""
    out = {}
    for f in folds:
        if (f, f"{arm}|cos") not in pt:
            continue
        carried = np.asarray(meta[f][key], bool)
        c = pt[(f, f"{arm}|cos")]; b = pt[(f, f"{BASE}|cos")]
        nc = ~carried
        out[f] = dict(n_noncarry=int(nc.sum()),
                      rel_out=float(c[nc].mean() / b[nc].mean() - 1) if nc.any() else 0.0,
                      rel_in=float(c[carried].mean() / b[carried].mean() - 1))
    return out


def fmt_table(summary, title):
    print("\n" + "=" * 118)
    print(title)
    print("=" * 118)
    print(f"{'arm':30s} {'nf':>2s} {'rel agg%':>9s} {'rel macro%':>10s} {'d pds':>10s} "
          f"{'d q_w':>10s} {'+cos':>4s} {'+pds':>4s} {'E rm':>6s} {'price(macro)':>12s} "
          f"{'price(agg)':>11s}")
    for arm, s in summary.items():
        pm = price(s["rel_macro"], s["d_pds"], k=0)["total_pds_from_plolo"]
        pa = price(s["rel_agg"], s["d_pds"], k=0)["total_pds_from_plolo"]
        s["price_macro"] = pm; s["price_agg"] = pa
        er = s["audit"].get("energy_removed_frac", float("nan"))
        print(f"{arm:30s} {s['n']:2d} {100*s['rel_agg']:+9.3f} {100*s['rel_macro']:+10.3f} "
              f"{s['d_pds']:+10.6f} {s['d_qw']:+10.6f} {s['cos_pos']:>4d} {s['pds_pos']:>4d} "
              f"{er:6.3f} {pm:+12.5f} {pa:+11.5f}")


def main() -> None:
    rows, pt, meta = load()
    arms = sorted({a for _, a in rows if a != BASE})
    have_cancer = [f for f in CANCER if (f, BASE) in rows]
    print(f"folds present: {[f for f in FOLDS if (f, BASE) in rows]}")

    fam = lambda pref: [a for a in arms if a.startswith(pref)]
    families = [
        ("A: deflated SELF-residual on top of the shipped stack (+controls)", fam("A_")),
        ("B: shipped residuals deflated at the shipped fraction (+controls)", fam("B_")),
        ("B2: CD4 deflation k-sweep, re-solved dose, energy-matched + shuffled-projector "
         "controls", fam("B2_")),
        ("B3: split-half projector transfer (basis fitted on the OTHER half of CD4's rows)",
         fam("B3_")),
        ("DOSE: the CD4 dose axis on the same folds (216's lever) for the exchange rate",
         fam("DOSE_")),
        ("D: MECHANISM -- deflated self-residual on the BARE base (218 self0.1 = none)",
         fam("D_")),
        ("C: CD4 partial-condition rows on top of the stack (+controls)", fam("C")),
        ("E: composed", fam("E_")),
        ("known-answer / reference rows", [a for a in arms if a in
                                          ("f0", "self0.1_bare", "cd4only")]),
    ]
    report = {}
    for title, fa in families:
        if not fa:
            continue
        s = summarise(rows, fa, have_cancer)
        fmt_table(s, f"{title}   -- {len(have_cancer)} CANCER folds vs STACK")
        report[title] = s
        if ("kolf2", BASE) in rows:
            sk = summarise(rows, fa, ("kolf2",))
            print("   kolf2 apart: " + "; ".join(
                f"{a} agg {100*v['rel_agg']:+.2f}% pds {v['d_pds']:+.5f}"
                for a, v in sk.items() if not a.endswith("_SHUF") and "rand" not in a
                and "SHUFDEF" not in a))
            report[title + " | kolf2"] = sk

    # falsifiers for the arms that matter
    key_arms = [a for a in arms if a.startswith("B2_cd4_svd") and "SHUF" not in a] + \
               [a for a in arms if a.startswith("B_cd4_svd")] + fam("B3_") + fam("DOSE_") + \
               fam("C2_partial@") + fam("C_partial@") + fam("E_") + fam("D_svd") + fam("A_svd2")
    print("\n" + "=" * 118)
    print("CONCENTRATION (cancer folds): leave-top-k-out retention, pds and cosine on their "
          "own orderings; last col = 222's cos-on-pds-order")
    print("=" * 118)
    conc = {}
    for arm in key_arms:
        c = concentration(pt, arm, have_cancer)
        if not c:
            continue
        conc[arm] = c
        c3 = c[3]
        print(f"{arm:30s} d pds full {c3['pds_full']:+.6f} LT3O {c3['pds_loo']:+.6f} "
              f"({100*c3['pds_retained']:6.1f}% ret) | d cos full {c3['cos_full']:+.6f} "
              f"LT3O {c3['cos_loo']:+.6f} ({100*c3['cos_retained']:6.1f}% ret) | "
              f"cos LT3O(pds order) {c3['cos_loo_on_pds_order']:+.6f}")

    print("\n" + "=" * 118)
    print("NON-CARRYING ROWS must not move (rel change of mean per-target cosine on rows the "
          "residual does not touch; cancer folds)")
    print("=" * 118)
    ncr = {}
    for arm, key in [(a, "carried_cd4") for a in arms if a.startswith(("B_cd4", "B2_cd4", "B3_"))
                     and "rand" not in a and "SHUF" not in a] + \
                    [(a, "carried_self") for a in arms if a.startswith("A_") and "SHUF" not in a
                     and "rand" not in a] + \
                    [(a, "carried_self") for a in arms if a.startswith("D_") and "rand" not in a
                     and "SHUF" not in a]:
        n = noncarrying(pt, meta, arm, have_cancer, key)
        if not n:
            continue
        ncr[arm] = n
        worst = max(abs(v["rel_out"]) for v in n.values())
        print(f"{arm:30s} worst |rel_out| {100*worst:8.5f} pp   "
              + "  ".join(f"{f[:8]} out {100*v['rel_out']:+.4f}% in {100*v['rel_in']:+.3f}%"
                          for f, v in n.items()))

    # leave-one-fold-out selection for the B2 k (dose) axis
    print("\n" + "=" * 118)
    print("LEAVE-ONE-FOLD-OUT k selection for B2 (fixed and resolved), priced on the held fold")
    print("=" * 118)
    lofo = {}
    for variant in ("", "_resolved"):
        cand = [a for a in arms if a.startswith("B2_cd4_svd") and a.endswith(variant)
                and "SHUF" not in a and (variant or "resolved" not in a)]
        held_prices = []
        for held in have_cancer:
            others = [f for f in have_cancer if f != held]
            best, bestv = None, -1e9
            for a in cand:
                s = summarise(rows, [a], others).get(a)
                if not s:
                    continue
                v = price(s["rel_macro"], s["d_pds"], k=0)["total_pds_from_plolo"]
                if v > bestv:
                    best, bestv = a, v
            if best is None:
                continue
            s = summarise(rows, [best], (held,))[best]
            v = price(s["rel_macro"], s["d_pds"], k=0)["total_pds_from_plolo"]
            held_prices.append((held, best, v, s["rel_macro"], s["d_pds"]))
        if held_prices:
            m = float(np.mean([v for _, _, v, _, _ in held_prices]))
            lofo[variant or "fixed"] = dict(mean=m, per_fold=held_prices)
            print(f"B2{variant or ' fixed':10s} out-of-sample price {m:+.5f}  " + ", ".join(
                f"{h[:8]}->{b.replace('B2_cd4_','')} {v:+.5f}" for h, b, v, _, _ in held_prices))

    (WORK / "analyse246.json").write_text(json.dumps(
        dict(report=report, concentration=conc, noncarrying=ncr, lofo=lofo,
             gates=dict(shape=GATE_SHAPE, amplitude=GATE_AMPL)), indent=1, default=float))
    print(f"\nwrote {WORK/'analyse246.json'}")


if __name__ == "__main__":
    main()
