"""299 -- the six-member price, the controls and the frozen gate for every reliability-selected-row operator arm,
over `cd4svd5`, on 261's recipe and the three COHERENT scenarios.  Known-answer / retrodiction rows first.

    PYTHONDONTWRITEBYTECODE=1 uv run python price299.py

261's `price261` is imported for its `compare`, `retention`, `price`, `factors`, `scenarios` and 257's constants,
so the arithmetic is the one that retrodicts BR249; 265's `boot265` for the paired target bootstrap.  Read-only
everywhere except this experiment's `work/`.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E261 = ROOT / "experiments/261-union-deflation"
E263 = ROOT / "experiments/263-donor-selection-arbiters"
E265 = ROOT / "experiments/265-majority-permutation-test"
sys.path.insert(0, str(E261))
sys.path.insert(0, str(E265))
import price261 as PR  # noqa: E402   (imports price257 as PR.P)
import boot265 as BT  # noqa: E402

P257 = PR.P
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
BASE = "cd4svd5"
CONTROL_KINDS = ("rand", "mag", "rev")
AMP_SHAPE_ONLY = 0.005          # the mandate: "delivered amplitude within 0.5 %" -> the shape-only gate
LINES: list[str] = []


def W(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def load():
    rows, pt, idn, meta = {}, {}, {}, {}
    for held in CANCER:
        p = WORK / f"sweep299-{held}.json"
        if not p.exists():
            continue
        d = json.loads(p.read_text())
        meta[held] = d["meta"]
        for r in d["rows"]:
            rows[(held, r["arm"])] = r
        z = np.load(WORK / f"pertarget299-{held}.npz")
        for k in z.files:
            pt[(held, k)] = z[k]
        z = np.load(WORK / f"identity299-{held}.npz")
        for k in z.files:
            idn[(held, k)] = z[k]
    return rows, pt, idn, meta


def groups_and_knee(pt, arm, base, folds) -> dict:
    """262's groups on the INCUMBENT's per-target transfer cosine (below chance <= 0; marginal 0 < c <= 0.05;
    saturated > 0.05) and the 0.02-knee crossings (R262 section 4.2 / R273 section 0.6)."""
    out = {g: [] for g in ("below", "marginal", "saturated")}
    n = {g: 0 for g in out}
    up = dn = 0
    for f in folds:
        c0 = pt[(f, f"{base}|cos")]
        c1 = pt[(f, f"{arm}|cos")]
        d = pt[(f, f"{arm}|pds")] - pt[(f, f"{base}|pds")]
        m = {"below": c0 <= 0.0, "marginal": (c0 > 0.0) & (c0 <= 0.05), "saturated": c0 > 0.05}
        for g in out:
            out[g].append(float(d[m[g]].sum() / d.size))          # contribution to the fold-mean d pds
            n[g] += int(m[g].sum())
        up += int(((c0 <= 0.02) & (c1 > 0.02)).sum())
        dn += int(((c0 > 0.02) & (c1 <= 0.02)).sum())
    return dict(contrib={g: float(np.mean(v)) for g, v in out.items()}, n_targets=n,
                knee_up=up, knee_down=dn, knee_net=up - dn)


def arm_report(rows, pt, idn, arm, folds, targets, rng=None) -> dict:
    # a fresh generator per arm, seeded deterministically from the arm name, so an arm's interval does not move
    # when other arms are added to the sweep
    rng = np.random.default_rng(BT.SEED + int.from_bytes(arm.encode(), "little") % 100000)
    c = PR.compare(rows, idn, arm, BASE, folds)
    ret = PR.retention(pt, arm, BASE, folds, targets, n=3)
    diffs = [pt[(f, f"{arm}|pds")] - pt[(f, f"{BASE}|pds")] for f in folds]
    bf = BT.boot_folds(diffs, rng)
    per = {f: BT.boot(d, rng) for f, d in zip(folds, diffs)}
    amp = c["delivered_ratio"] - 1.0
    sc = PR.scenarios(c)
    gk = groups_and_knee(pt, arm, BASE, folds)
    gate = P257.GATE_SHAPE if abs(amp) <= AMP_SHAPE_ONLY else P257.GATE_AMPL
    return dict(arm=arm, cmp=c, retention=ret, boot=bf, boot_per_fold=per, scen=sc, groups=gk,
                amp_pct=100 * amp, gate=gate, shape_only=bool(abs(amp) <= AMP_SHAPE_ONLY))


def main() -> None:
    t0 = time.time()
    rows, pt, idn, meta = load()
    folds = [f for f in CANCER if f in meta]
    targets = {f: meta[f]["targets"] for f in folds}
    rng = np.random.default_rng(BT.SEED)
    arms = []
    for (f, a) in rows:
        if a not in arms and a not in (BASE, "STACK"):
            arms.append(a)
    arms = [a for a in arms if all((f, a) in rows for f in folds)]
    order = {"full": 0, "sel:q0.5": 1, "sel:q0.25": 2, "rand:q0.5": 3, "mag:q0.5": 4, "rev:q0.5": 5, "wgt:mild": 6, "wgt:mildrev": 7, "wgt:R": 8}
    arms.sort(key=lambda a: (order.get(a, 3), a))

    W("=" * 132)
    W("299 -- the whitening operator fitted on RELIABILITY-SELECTED training rows, over cd4svd5 (PLOLO 162, three cancer folds)")
    W("=" * 132)
    W(f"folds read: {folds}; arms: {arms}")
    fails = []
    for f in folds:
        m = meta[f]
        W(f"  {f:14s} {m['n_targets']} targets, wcols {m['n_wcols']}, known-answer {m['ka_pass']}/{m['ka_total']}, "
          f"failures {m['failures']}, {m['seconds']:.0f}s, peak(ru_maxrss incl. mmap page cache) {m['peak_rss_gib']:.2f} GiB")
        if m["failures"] or m["ka_pass"] != m["ka_total"]:
            fails.append(f"sweep {f}")

    # ------------------------------------------------------------------ known-answer rows
    W("\n" + "=" * 132)
    W("KNOWN-ANSWER AND RETRODICTION ROWS (nothing is priced unless they pass)")
    W("=" * 132)
    ka_all = [k for f in folds for k in meta[f]["ka"]]
    for k in ka_all:
        W(f"  {'PASS' if k['ok'] else 'FAIL'}  {k['label'][:104]:104s} |d| {k['abs_diff']:.2e} (tol {k['tol']:.0e})")
        if not k["ok"]:
            fails.append(k["label"])
    # R2: the BR249 retrodiction, cd4svd5 over STACK (h1union) on the POINT set
    c249 = PR.compare(rows, idn, BASE, "STACK", folds)
    # 261 R2: the retrodiction is priced on h1union's OWN board terms (it is the base of that comparison)
    s249 = PR.p257_sets(c249, P257.COS_H1UNION, P257.PER_H1UNION)["mine"]
    ref = json.loads((E261 / "work/price261.json").read_text()) if (E261 / "work/price261.json").exists() else None
    W(f"\nR2  BR249 retrodiction (cd4svd5 over h1union, 261's POINT set on h1union's own terms): "
      f"overall {s249['point']['total']:+.7f} (pds {s249['point']['pds']:+.5f}, mse {s249['point']['mse']:+.5f}, "
      f"DE-three {s249['point']['de3']:+.5f}, nmae {s249['point']['nmae']:+.5f})   board realised +0.0042634")
    W(f"    legs: d pds {c249['d_pds']:+.9f} (per fold {['%+.6f' % x for x in c249['per_fold_d_pds']]}), "
      f"cos_agg {100*c249['rel_agg']:+.6f}%, d nmae89 {c249['d_nmae']:+.9f}, delivered ratio {c249['delivered_ratio']:.9f}")
    r2_ok = abs(s249["point"]["total"] - 0.00472) < 5e-5
    if ref is not None:
        r261 = ref["retrodiction"]["R2"]
        dcmp = max(abs(float(r261["compare"][k]) - float(c249[k])) for k in ("d_pds", "rel_agg", "rel_macro", "d_nmae", "delivered_ratio"))
        dpri = abs(float(r261["price"]["point"]["total"]) - float(s249["point"]["total"]))
        W(f"    against 261's own R2 row: compare max|d| {dcmp:.1e}; POINT total |d| {dpri:.1e} "
          f"(261 {r261['price']['point']['total']:.9f})")
        r2_ok = r2_ok and dcmp == 0.0 and dpri <= 1e-15
    W(f"    R2 {'PASS' if r2_ok else 'CHECK'} (257/261 quote +0.00472)")
    s249 = dict(POINT=s249["point"], PESSIMISTIC=s249["conservative"], OPTIMISTIC=s249["point"])
    if not r2_ok:
        fails.append("R2")

    # ------------------------------------------------------------------ the arm table
    W("\n" + "=" * 132)
    W("ARM x FOLD -- d pds per fold, cosine, amplitude (everything over cd4svd5)")
    W("=" * 132)
    W(f"{'arm':12s} " + " ".join(f"{f[:11]:>11s}" for f in folds) +
      f" {'mean dpds':>10s} {'3/3':>4s} {'cos_agg%':>9s} {'macro%':>8s} {'d nmae89':>9s} {'amp':>8s} {'ret%':>6s}")
    reps: dict[str, dict] = {}
    for a in arms:
        rep = arm_report(rows, pt, idn, a, folds, targets, rng)
        reps[a] = rep
        c = rep["cmp"]
        W(f"{a:12s} " + " ".join(f"{x:+11.5f}" for x in c["per_fold_d_pds"]) +
          f" {c['d_pds']:+10.5f} {c['pds_pos']:>3d}/{len(folds)} {100*c['rel_agg']:+9.3f} {100*c['rel_macro']:+8.3f} "
          f"{c['d_nmae']:+9.5f} {rep['amp_pct']:+7.3f}% {100*rep['retention']['retained']:+6.1f}")

    # ------------------------------------------------------------------ per-arm detail
    W("\n" + "=" * 132)
    W("PRICE (261's recipe, three coherent scenarios), BOOTSTRAP, GROUPS, KNEE, GATE")
    W("=" * 132)
    verdicts = {}
    for a in arms:
        rep = reps[a]
        c, s, g = rep["cmp"], rep["scen"], rep["groups"]
        W(f"\n--- {a}  ({'shape-only' if rep['shape_only'] else 'amplitude-changing'}, delivered {rep['amp_pct']:+.3f}%, "
          f"gate {rep['gate']:+.4f})")
        W(f"    d pds {c['d_pds']:+.6f}  pooled bootstrap {BT.fmt_b(rep['boot'])}")
        for f in rep["boot_per_fold"]:
            W(f"      {f:14s} {BT.fmt_b(rep['boot_per_fold'][f])}")
        W(f"    cos_agg {100*c['rel_agg']:+.4f}%  macro {100*c['rel_macro']:+.4f}%  d nmae89 {c['d_nmae']:+.6f}  "
          f"d prec89 {c['d_prec89']:+.6f}  d a15_sgn_w {c['d_qw']:+.6f}")
        W(f"    262 groups (contribution to the fold-mean d pds): below-chance {g['contrib']['below']:+.6f} "
          f"(n {g['n_targets']['below']}), marginal {g['contrib']['marginal']:+.6f} (n {g['n_targets']['marginal']}), "
          f"saturated {g['contrib']['saturated']:+.6f} (n {g['n_targets']['saturated']}); "
          f"0.02-knee crossings up {g['knee_up']} down {g['knee_down']} net {g['knee_net']:+d}")
        W(f"    leave-top-3-out: full {rep['retention']['full']:+.6f} -> {rep['retention']['leave_top_out']:+.6f} "
          f"= {100*rep['retention']['retained']:+.1f}% retained (top-3 share {100*rep['retention']['share']:.1f}%)")
        for w in ("PESSIMISTIC", "POINT", "OPTIMISTIC"):
            x = s[w]
            W(f"    {w:12s} total {x['total']:+.5f}  = pds {x['pds']:+.5f} + mse {x['mse']:+.5f} + DE-three {x['de3']:+.5f} "
              f"+ nmae {x['nmae']:+.5f}   (k_pds {x['k_pds']}, k_cos {x['k_cos']}, d cos board {x['dcos_board']:+.6f})")
        # ---- the gate
        cond = {}
        cond["price"] = (s["PESSIMISTIC"]["total"] >= rep["gate"], f"{s['PESSIMISTIC']['total']:+.5f} vs {rep['gate']:+.4f}")
        cond["3/3"] = (c["pds_pos"] == len(folds) and rep["boot"]["lo"] > 0,
                       f"{c['pds_pos']}/{len(folds)} folds positive, pooled interval [{rep['boot']['lo']:+.5f}, {rep['boot']['hi']:+.5f}]")
        cond["retention"] = (rep["retention"]["retained"] >= 0.5, f"{100*rep['retention']['retained']:+.1f}% vs 50%")
        MIRROR = {"wgt:mild": {"mirror (wgt:mildrev)": "wgt:mildrev"},
                  "wgt:mildrev": {"mirror (wgt:mild)": "wgt:mild"},
                  "wgt:R": {"mirror (wgt:rev)": "wgt:rev"},
                  "wgt:steep": {"mirror (wgt:steeprev)": "wgt:steeprev"},
                  "wgt:steeprev": {"mirror (wgt:steep)": "wgt:steep"}}
        if a.startswith("sel") or a in MIRROR:
            q = a.split(":")[1]
            ctl = MIRROR[a] if a in MIRROR else {k: next((x for x in arms if x.startswith(f"{k}:{q}")), None) for k in CONTROL_KINDS}
            ctl = {k: (v if v in reps else None) for k, v in ctl.items()}
            det = []
            ok = True
            for k, nm in ctl.items():
                if nm is None:
                    det.append(f"{k} ABSENT")
                    ok = False
                    continue
                d = c["d_pds"] - reps[nm]["cmp"]["d_pds"]
                need = 0.5 * abs(c["d_pds"])
                ok = ok and (d > need)
                det.append(f"{k} {reps[nm]['cmp']['d_pds']:+.5f} (margin {d:+.5f} vs {need:+.5f})")
            cond["controls"] = (ok, "; ".join(det))
        for k, (ok, txt) in cond.items():
            W(f"    GATE {k:10s} {'PASS' if ok else 'FAIL'}  {txt}")
        v = all(ok for ok, _ in cond.values())
        verdicts[a] = dict(clears=bool(v), conditions={k: dict(pass_=bool(ok), detail=txt) for k, (ok, txt) in cond.items()})
        W(f"    -> {a} {'CLEARS' if v else 'DOES NOT CLEAR'}")

    # ------------------------------------------------------------------ headline
    prim = [a for a in arms if a.startswith("sel") or a.startswith("wgt")]
    best = max(prim, key=lambda a: reps[a]["cmp"]["d_pds"]) if prim else None
    W("\n" + "=" * 132)
    if best is not None:
        r = reps[best]
        c = r["cmp"]
        ctls = {k: (v if v in reps else None) for k, v in (MIRROR[best].items() if best in MIRROR else {k: next((x for x in arms if x.startswith(f"{k}:{best.split(chr(58))[1]}")), None) for k in CONTROL_KINDS}.items())}
        cs = ", ".join(f"{k} {reps[n]['cmp']['d_pds']:+.5f}" if n else f"{k} n/a" for k, n in ctls.items())
        W(f"RELIABLE-ROW WHITENER: best arm {best} d pds {c['d_pds']:+.5f} [{r['boot']['lo']:+.5f}, {r['boot']['hi']:+.5f}], "
          f"3/3: {'YES' if c['pds_pos'] == len(folds) else 'NO'} ({c['pds_pos']}/{len(folds)}), controls ({cs}), "
          f"priced {r['scen']['PESSIMISTIC']['total']:+.5f} pessimistic vs gate {r['gate']:+.4f} "
          f"-> {'CLEARS' if verdicts[best]['clears'] else 'DOES NOT CLEAR'}")
    W("=" * 132)
    out = dict(folds=folds, arms=arms, R2=dict(point_total=s249["POINT"]["total"], compare=c249, scen=s249, ok=bool(r2_ok)),
               ka=ka_all, failures=fails, verdicts=verdicts,
               reports={a: dict(cmp=reps[a]["cmp"], scen=reps[a]["scen"], boot=reps[a]["boot"],
                                boot_per_fold=reps[a]["boot_per_fold"], retention=reps[a]["retention"],
                                groups=reps[a]["groups"], amp_pct=reps[a]["amp_pct"], gate=reps[a]["gate"],
                                shape_only=reps[a]["shape_only"]) for a in arms},
               seconds=round(time.time() - t0, 1))
    (WORK / "price299.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "price299.txt").write_text("\n".join(LINES) + "\n")
    W(f"\nfailures: {fails}")


if __name__ == "__main__":
    main()
