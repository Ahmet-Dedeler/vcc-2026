"""300 -- the six-member price, the controls, the frozen gate and the frozen P1 / P2 / P3 rules for every arm,
over `cd4svd5`, on 261's recipe and the three COHERENT scenarios.  Known-answer / retrodiction rows first.

    PYTHONDONTWRITEBYTECODE=1 uv run python price300.py

261's `price261` supplies `compare` / `retention` / `price` / `factors` / `scenarios` and 257's constants, so the
arithmetic is the one that retrodicts BR249; 265's `boot265` supplies the paired target bootstrap.  Read-only
everywhere except this experiment's `work/`.
"""
from __future__ import annotations

import json
import re
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
E299 = ROOT / "experiments/299-reliable-training-rows"
sys.path.insert(0, str(E261))
sys.path.insert(0, str(E265))
import price261 as PR  # noqa: E402   (imports price257 as PR.P)
import boot265 as BT  # noqa: E402

P257 = PR.P
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
BASE = "cd4svd5"
AMP_SHAPE_ONLY = 0.005
A_DOSES = (0.5, 1.0, 1.25, 1.5, 1.8, 2.0, 2.5)     # 1.8 is 299's `steeprev`, carried as KA-W3 and as a sixth A point
LINES: list[str] = []


def W(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def rng_for(name: str) -> np.random.Generator:
    return np.random.default_rng(BT.SEED + int.from_bytes(name.encode(), "little") % 100000)


def load():
    rows, pt, idn, meta = {}, {}, {}, {}
    for held in CANCER:
        p = WORK / f"sweep300-{held}.json"
        if not p.exists():
            continue
        d = json.loads(p.read_text())
        meta[held] = d["meta"]
        for r in d["rows"]:
            rows[(held, r["arm"])] = r
        z = np.load(WORK / f"pertarget300-{held}.npz")
        for k in z.files:
            pt[(held, k)] = z[k]
        z = np.load(WORK / f"identity300-{held}.npz")
        for k in z.files:
            idn[(held, k)] = z[k]
    return rows, pt, idn, meta


def groups_and_knee(pt, arm, base, folds) -> dict:
    out = {g: [] for g in ("below", "marginal", "saturated")}
    n = {g: 0 for g in out}
    up = dn = 0
    for f in folds:
        c0 = pt[(f, f"{base}|cos")]
        c1 = pt[(f, f"{arm}|cos")]
        d = pt[(f, f"{arm}|pds")] - pt[(f, f"{base}|pds")]
        m = {"below": c0 <= 0.0, "marginal": (c0 > 0.0) & (c0 <= 0.05), "saturated": c0 > 0.05}
        for g in out:
            out[g].append(float(d[m[g]].sum() / d.size))
        for g in out:
            n[g] += int(m[g].sum())
        up += int(((c0 <= 0.02) & (c1 > 0.02)).sum())
        dn += int(((c0 > 0.02) & (c1 <= 0.02)).sum())
    return dict(contrib={g: float(np.mean(v)) for g, v in out.items()}, n_targets=n,
                knee_up=up, knee_down=dn, knee_net=up - dn)


ARMPAT = re.compile(r"(tilt|spread|prune|rand):d(-?[\d.]+)(?::s(\d+))?")
IDX: dict = {}


def index_arms(arms: list[str]) -> None:
    IDX.clear()
    for a in arms:
        m = ARMPAT.fullmatch(a)
        if m:
            IDX[(m.group(1), float(m.group(2)), m.group(3))] = a


def AN(kind: str, dose: float, seed: str | None = None) -> str | None:
    """the arm name as run300 wrote it (`tilt:d1.0`), looked up by (kind, dose, seed)"""
    return IDX.get((kind, float(dose), seed))


def mirror_of(arm: str) -> str | None:
    m = ARMPAT.fullmatch(arm)
    if not m or m.group(1) == "rand":
        return None
    d = float(m.group(2))
    return AN(m.group(1), -d) if d != 0 else None


def arm_report(rows, pt, idn, arm, folds, targets) -> dict:
    rng = rng_for(arm)
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
                amp_pct=100 * amp, gate=gate, shape_only=bool(abs(amp) <= AMP_SHAPE_ONLY), diffs=diffs)


def paired_diff(reps, a: str, b: str, label: str) -> dict:
    """the paired target bootstrap of (arm a - arm b), both measured against `cd4svd5` on the same targets"""
    d = [x - y for x, y in zip(reps[a]["diffs"], reps[b]["diffs"])]
    return BT.boot_folds(d, rng_for(label))


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean()
    ry -= ry.mean()
    den = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    return float((rx * ry).sum() / den) if den > 0 else float("nan")


def spearman_boot(reps, a: str, b: str, label: str, B: int = 20000) -> dict:
    xs = np.concatenate(reps[a]["diffs"])
    ys = np.concatenate(reps[b]["diffs"])
    obs = spearman(xs, ys)
    rng = rng_for(label)
    n = xs.size
    draws = np.empty(B)
    for i in range(B):
        j = rng.integers(0, n, n)
        draws[i] = spearman(xs[j], ys[j])
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return dict(rho=obs, lo=float(lo), hi=float(hi), n=int(n))


def main() -> None:
    t0 = time.time()
    rows, pt, idn, meta = load()
    folds = [f for f in CANCER if f in meta]
    targets = {f: meta[f]["targets"] for f in folds}
    # only arms measured on EVERY read fold enter the table
    arms = []
    for (f, a) in rows:
        if a not in arms and a not in (BASE, "STACK"):
            arms.append(a)
    partial = [a for a in arms if not all((f, a) in rows for f in folds)]
    arms = [a for a in arms if all((f, a) in rows for f in folds)]

    def okey(a: str):
        m = re.fullmatch(r"(tilt|spread|prune|rand):d(-?[\d.]+)(?::s(\d+))?", a)
        grp = {"tilt": 0, "prune": 1, "spread": 2, "rand": 3}
        return (grp.get(m.group(1), 9), float(m.group(2)), m.group(3) or "") if m else (9, 0.0, a)
    arms.sort(key=okey)
    index_arms(arms)

    W("=" * 140)
    W("300 -- the SPREAD mechanism behind the whitener's rank-retrieval gain, over cd4svd5 (PLOLO 162, three cancer folds)")
    W("=" * 140)
    W(f"folds read: {folds}")
    W(f"arms on every fold: {arms}")
    W(f"arms NOT on every fold (excluded from the table, listed in RESULTS.md): {partial}")
    fails = []
    for f in folds:
        m = meta[f]
        W(f"  {f:14s} {m['n_targets']} targets, wcols {m['n_wcols']}, known-answer {m['ka_pass']}/{m['ka_total']}, "
          f"failures {m['failures']}, {m['seconds']:.0f}s, peak(ru_maxrss incl. mmap page cache) {m['peak_rss_gib']:.2f} GiB")
        if m["failures"] or m["ka_pass"] != m["ka_total"]:
            fails.append(f"sweep {f}")

    # ------------------------------------------------------------------ known-answer rows
    W("\n" + "=" * 140)
    W("KNOWN-ANSWER AND RETRODICTION ROWS (nothing is priced unless they pass)")
    W("=" * 140)
    ka_all = [dict(fold=f, **k) for f in folds for k in meta[f]["ka"]]
    for k in ka_all:
        W(f"  {'PASS' if k['ok'] else 'FAIL'}  [{k['fold'][:13]:13s}] {k['label'][:96]:96s} |d| {k['abs_diff']:.2e} (tol {k['tol']:.0e})")
        if not k["ok"]:
            fails.append(k["label"])
    c249 = PR.compare(rows, idn, BASE, "STACK", folds)
    s249 = PR.p257_sets(c249, P257.COS_H1UNION, P257.PER_H1UNION)["mine"]
    ref = json.loads((E261 / "work/price261.json").read_text()) if (E261 / "work/price261.json").exists() else None
    W(f"\nR2  BR249 retrodiction (cd4svd5 over h1union, 261's POINT set on h1union's own terms): "
      f"overall {s249['point']['total']:+.7f}   board realised +0.0042634")
    r2_ok = abs(s249["point"]["total"] - 0.00472) < 5e-5
    if ref is not None:
        r261 = ref["retrodiction"]["R2"]
        dcmp = max(abs(float(r261["compare"][k]) - float(c249[k])) for k in ("d_pds", "rel_agg", "rel_macro", "d_nmae", "delivered_ratio"))
        dpri = abs(float(r261["price"]["point"]["total"]) - float(s249["point"]["total"]))
        W(f"    against 261's own R2 row: compare max|d| {dcmp:.1e}; POINT total |d| {dpri:.1e} (261 {r261['price']['point']['total']:.9f})")
        r2_ok = r2_ok and dcmp == 0.0 and dpri <= 1e-15
    W(f"    R2 {'PASS' if r2_ok else 'CHECK'} (257/261 quote +0.00472)")
    if not r2_ok:
        fails.append("R2")

    # ------------------------------------------------------------------ the arm table
    reps: dict[str, dict] = {}
    for a in arms:
        reps[a] = arm_report(rows, pt, idn, a, folds, targets)
    W("\n" + "=" * 140)
    W("ARM x FOLD -- d pds per fold, cosine, amplitude, price, gate (everything over cd4svd5)")
    W("=" * 140)
    W(f"{'arm':20s} " + " ".join(f"{f[:11]:>11s}" for f in folds) +
      f" {'mean dpds':>10s} {'95% pooled':>22s} {'3/3':>4s} {'cos_agg%':>9s} {'d nmae89':>9s} {'amp%':>7s} {'ret%':>7s} {'pessim':>9s} {'gate':>8s}")
    for a in arms:
        r, c = reps[a], reps[a]["cmp"]
        W(f"{a:20s} " + " ".join(f"{x:+11.5f}" for x in c["per_fold_d_pds"]) +
          f" {c['d_pds']:+10.5f} [{r['boot']['lo']:+.5f},{r['boot']['hi']:+.5f}] {c['pds_pos']:>3d}/{len(folds)} "
          f"{100*c['rel_agg']:+9.3f} {c['d_nmae']:+9.5f} {r['amp_pct']:+6.3f}% {100*r['retention']['retained']:+6.1f} "
          f"{r['scen']['PESSIMISTIC']['total']:+9.5f} {r['gate']:+8.4f}")

    # ------------------------------------------------------------------ per-arm detail and the gate
    W("\n" + "=" * 140)
    W("PRICE (261's recipe, three coherent scenarios), BOOTSTRAP, GROUPS, KNEE, GATE")
    W("=" * 140)
    verdicts = {}
    RAND = [a for a in arms if a.startswith("rand:")]
    for a in arms:
        r = reps[a]
        c, s, g = r["cmp"], r["scen"], r["groups"]
        W(f"\n--- {a}  ({'shape-only' if r['shape_only'] else 'amplitude-changing'}, delivered {r['amp_pct']:+.3f}%, gate {r['gate']:+.4f})")
        W(f"    d pds {c['d_pds']:+.6f}  pooled bootstrap {BT.fmt_b(r['boot'])}")
        for f in r["boot_per_fold"]:
            W(f"      {f:14s} {BT.fmt_b(r['boot_per_fold'][f])}")
        W(f"    cos_agg {100*c['rel_agg']:+.4f}%  macro {100*c['rel_macro']:+.4f}%  d nmae89 {c['d_nmae']:+.6f}  d prec89 {c['d_prec89']:+.6f}")
        W(f"    262 groups (contribution to the fold-mean d pds): below-chance {g['contrib']['below']:+.6f} (n {g['n_targets']['below']}), "
          f"marginal {g['contrib']['marginal']:+.6f} (n {g['n_targets']['marginal']}), saturated {g['contrib']['saturated']:+.6f} "
          f"(n {g['n_targets']['saturated']}); 0.02-knee up {g['knee_up']} down {g['knee_down']} net {g['knee_net']:+d}")
        W(f"    leave-top-3-out: {r['retention']['full']:+.6f} -> {r['retention']['leave_top_out']:+.6f} "
          f"= {100*r['retention']['retained']:+.1f}% retained (top-3 share {100*r['retention']['share']:.1f}%)")
        for wname in ("PESSIMISTIC", "POINT", "OPTIMISTIC"):
            x = s[wname]
            W(f"    {wname:12s} total {x['total']:+.5f}  = pds {x['pds']:+.5f} + mse {x['mse']:+.5f} + DE-three {x['de3']:+.5f} "
              f"+ nmae {x['nmae']:+.5f}   (k_pds {x['k_pds']}, k_cos {x['k_cos']})")
        cond = {}
        cond["price"] = (s["PESSIMISTIC"]["total"] >= r["gate"], f"{s['PESSIMISTIC']['total']:+.5f} vs {r['gate']:+.4f}")
        cond["3/3"] = (c["pds_pos"] == len(folds) and r["boot"]["lo"] > 0,
                       f"{c['pds_pos']}/{len(folds)} folds positive, pooled interval [{r['boot']['lo']:+.5f}, {r['boot']['hi']:+.5f}]")
        cond["retention"] = (r["retention"]["retained"] >= 0.5, f"{100*r['retention']['retained']:+.1f}% vs 50%")
        mir = mirror_of(a)
        ctl = {}
        if mir is not None:
            ctl["mirror"] = mir if mir in reps else None
        if not a.startswith("rand:"):
            for rd in RAND:
                ctl[rd] = rd
        if ctl:
            det, ok = [], True
            need = 0.5 * abs(c["d_pds"])
            for k, nm in ctl.items():
                if nm is None:
                    det.append(f"{k} ABSENT")
                    ok = False
                    continue
                dd = c["d_pds"] - reps[nm]["cmp"]["d_pds"]
                ok = ok and (dd > need)
                det.append(f"{k}={nm} {reps[nm]['cmp']['d_pds']:+.5f} (margin {dd:+.5f} vs {need:+.5f})")
            cond["controls"] = (ok, "; ".join(det))
        for k, (ok, txt) in cond.items():
            W(f"    GATE {k:10s} {'PASS' if ok else 'FAIL'}  {txt}")
        v = all(ok for ok, _ in cond.values())
        verdicts[a] = dict(clears=bool(v), conditions={k: dict(pass_=bool(ok), detail=txt) for k, (ok, txt) in cond.items()})
        W(f"    -> {a} {'CLEARS' if v else 'DOES NOT CLEAR'}")

    # ------------------------------------------------------------------ P1 / P2 / P3
    W("\n" + "=" * 140)
    W("THE THREE PREREGISTERED PREDICTIONS (rules frozen in PREREGISTRATION.md section 4)")
    W("=" * 140)
    out_p: dict = {}
    A = [AN("tilt", d) for d in A_DOSES if AN("tilt", d) in reps]
    shape_A = [a for a in A if reps[a]["shape_only"]]
    W(f"\nP1 -- the A dose curve ({len(A)} doses measured: {A}; shape-only: {shape_A})")
    W(f"    {'dose':>6s} {'d pds':>10s} {'95%':>22s} {'3/3':>5s} {'cos_agg%':>9s} {'amp%':>8s} {'shape-only':>11s}")
    curve = []
    for a in A:
        r, c = reps[a], reps[a]["cmp"]
        curve.append(c["d_pds"])
        W(f"    {a.split('d')[-1]:>6s} {c['d_pds']:+10.5f} [{r['boot']['lo']:+.5f},{r['boot']['hi']:+.5f}] "
          f"{c['pds_pos']:>3d}/{len(folds)} {100*c['rel_agg']:+9.3f} {r['amp_pct']:+7.3f}% {str(r['shape_only']):>11s}")
    mono = all(curve[i + 1] >= curve[i] for i in range(len(curve) - 1)) if len(curve) > 1 else False
    best_shape = max(shape_A, key=lambda a: reps[a]["cmp"]["d_pds"]) if shape_A else None
    p1 = dict(doses=A, curve=curve, monotone=bool(mono), best_shape_only=best_shape)
    A1 = AN("tilt", 1.0)
    if best_shape and A1 in reps:
        b = paired_diff(reps, best_shape, A1, "P1")
        p1["vs_d1"] = b
        W(f"    monotone non-decreasing across the measured doses: {mono}")
        W(f"    best shape-only dose {best_shape} minus d=1.0: {BT.fmt_b(b)}")
        p1["pass"] = bool(mono and b["lo"] > 0)
        if best_shape == A1:
            p1["pass"] = False
            W("    (the best shape-only dose IS d = 1.0, so the second clause cannot be met)")
    W(f"    P1 {'PASSES' if p1.get('pass') else 'DOES NOT PASS'} by its frozen rule")
    out_p["P1"] = p1

    dstar = float(ARMPAT.fullmatch(best_shape).group(2)) if best_shape else 1.0
    Bd = AN("prune", dstar)
    P0 = AN("prune", 0.0)
    W(f"\nP2 -- the essential-gene prune's interaction with the tilt (d* = {dstar})")
    p2: dict = dict(dstar=dstar)
    if Bd in reps and best_shape in reps:
        letter = paired_diff(reps, Bd, best_shape, "P2letter")
        p2["letter"] = letter
        W(f"    prune alone ({P0}): {reps[P0]['cmp']['d_pds']:+.6f} {BT.fmt_b(reps[P0]['boot'])}"
          if P0 in reps else "    prune alone: NOT MEASURED")
        W(f"    the letter of the rule, (B at d*) - (A at d*): {reps[Bd]['cmp']['d_pds']:+.6f} - "
          f"{reps[best_shape]['cmp']['d_pds']:+.6f} = {BT.fmt_b(letter)}")
        if P0 in reps:
            # the interaction: [B(d*) - B(0)] - [A(d*) - A(0)]; A(0) is the incumbent, so A(d*) - A(0) = A's own d pds
            d = [(pt[(f, f'{Bd}|pds')] - pt[(f, f'{P0}|pds')]) - (pt[(f, f'{best_shape}|pds')] - pt[(f, f'{BASE}|pds')])
                 for f in folds]
            inter = BT.boot_folds(d, rng_for("P2inter"))
            p2["interaction"] = inter
            W(f"    the interaction [B(d*)-B(0)] - [A(d*)-A(0)]: {BT.fmt_b(inter)}")
            p2["mechanism_confirmed"] = bool(inter["hi"] < 0)
            W(f"    P2's mechanism prediction (the tilt gains LESS on top of the prune) "
              f"{'CONFIRMED' if inter['hi'] < 0 else ('REFUTED' if inter['lo'] > 0 else 'NOT SEPARATED FROM ZERO')}")
    else:
        W(f"    NOT MEASURED ({Bd} or {best_shape} missing on some fold)")
    out_p["P2"] = p2

    W("\nP3 -- the direct spread manipulation")
    p3: dict = {}
    C1, A1 = AN("spread", 1.0), AN("tilt", 1.0)
    if C1 in reps and A1 in reps:
        sp = spearman_boot(reps, C1, A1, "P3rho")
        p3["spearman"] = sp
        W(f"    Spearman(per-target d pds of {C1}, of {A1}) pooled over {sp['n']} targets: "
          f"{sp['rho']:+.4f} [{sp['lo']:+.4f}, {sp['hi']:+.4f}]")
        Cs = [a for a in arms if a.startswith("spread:d") and float(a.split("d")[-1]) > 0]
        bestC = max(Cs, key=lambda a: reps[a]["cmp"]["d_pds"]) if Cs else None
        p3["bestC"] = bestC
        if bestC:
            dd = paired_diff(reps, bestC, A1, "P3size")
            p3["vs_A"] = dd
            W(f"    C's best dose {bestC} {reps[bestC]['cmp']['d_pds']:+.6f} vs A at the matched dose "
              f"{reps[A1]['cmp']['d_pds']:+.6f}: {BT.fmt_b(dd)}")
            p3["pass"] = bool(sp["lo"] > 0 and dd["lo"] >= 0 or (sp["lo"] > 0 and dd["hi"] > 0 and reps[bestC]["cmp"]["d_pds"] >= reps[A1]["cmp"]["d_pds"]))
            W(f"    C's own gate: {'CLEARS' if verdicts.get(bestC, {}).get('clears') else 'DOES NOT CLEAR'}")
    else:
        W("    NOT MEASURED")
    W(f"    P3 {'PASSES' if p3.get('pass') else 'DOES NOT PASS'} by its frozen rule")
    out_p["P3"] = p3

    # ------------------------------------------------------------------ headline
    prim = [a for a in arms if a.startswith(("tilt:d", "spread:d")) and float(a.split("d")[-1]) > 0]
    best = max(prim, key=lambda a: reps[a]["cmp"]["d_pds"]) if prim else None
    bestshape = max([a for a in prim if reps[a]["shape_only"]], key=lambda a: reps[a]["cmp"]["d_pds"], default=None)
    W("\n" + "=" * 140)
    if bestshape:
        r, c = reps[bestshape], reps[bestshape]["cmp"]
        mir = mirror_of(bestshape)
        cs = []
        if mir in reps:
            cs.append(f"mirror {reps[mir]['cmp']['d_pds']:+.5f}")
        for rd in RAND:
            cs.append(f"{rd} {reps[rd]['cmp']['d_pds']:+.5f}")
        W(f"SPREAD MECHANISM: P1 {'PASS' if out_p['P1'].get('pass') else 'FAIL'} / "
          f"P2 {'mechanism confirmed' if out_p['P2'].get('mechanism_confirmed') else 'mechanism not confirmed'} / "
          f"P3 {'PASS' if out_p['P3'].get('pass') else 'FAIL'}; best shape-only arm {bestshape} d pds {c['d_pds']:+.5f} "
          f"[{r['boot']['lo']:+.5f}, {r['boot']['hi']:+.5f}], 3/3 {c['pds_pos']}/{len(folds)}, retention "
          f"{100*r['retention']['retained']:+.1f}%, controls ({'; '.join(cs)}), priced {r['scen']['PESSIMISTIC']['total']:+.5f} "
          f"vs {r['gate']:+.4f} -> {'CLEARS' if verdicts[bestshape]['clears'] else 'DOES NOT CLEAR'}")
    W("=" * 140)
    out = dict(folds=folds, arms=arms, partial=partial, ka=ka_all, failures=fails, verdicts=verdicts,
               predictions=out_p, best=best, best_shape_only=bestshape,
               R2=dict(point_total=s249["point"]["total"], compare=c249, ok=bool(r2_ok)),
               reports={a: dict(cmp=reps[a]["cmp"], scen=reps[a]["scen"], boot=reps[a]["boot"],
                                boot_per_fold=reps[a]["boot_per_fold"], retention=reps[a]["retention"],
                                groups=reps[a]["groups"], amp_pct=reps[a]["amp_pct"], gate=reps[a]["gate"],
                                shape_only=reps[a]["shape_only"]) for a in arms},
               seconds=round(time.time() - t0, 1))
    (WORK / "price300.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "price300.txt").write_text("\n".join(LINES) + "\n")
    W(f"\nfailures: {fails}")


if __name__ == "__main__":
    main()
