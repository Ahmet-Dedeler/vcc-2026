"""534 Part B -- aggregation and the frozen reading rules (PREREGISTRATION-B sections 6 and 7). Light (jsons + arrays).

Inputs: every R6 panel whose three folds and full panel are priced (S1..S4, in order), efp's known-answer run (the
exact milestone designs on efp, reported beside the record), record534.json (the record's values).
Primary: 438's LEVELLED route at k = 1 (price534). Per panel: the mean of its three cancer folds. Pooled over R6: the mean
over every fold of every priced panel; paired target bootstrap, 4,000 draws, seed 534 (price534.pooled_bootstrap).
Section 7:
  1. verdict per step (pooled): TRANSFERS if D > 0 and the 95 % lower bound > 0; FRAGILE if D > 0 and lb <= 0; FAILS if D <= 0
  2. retention = R6 pooled D / A/B/C D (the record's value), per step and for the cumulative M0 -> M6 (sum of the six steps)
  3. tuned-on: real/efp/efr OVER-READ if, for >= 2 of M3->M4, M4->M5, M5->M6, R6 pooled D < 0.5 x mean(record real/efp/efr)
  4. controls: C-a and C-b must read <= 0 pooled; C-a > 0 with its interval above 0 flags 438's CD4 x1 decision
  5. consequence: FAILS -> candidate to switch off; FRAGILE -> flagged. Nothing is changed here.
Projection (section 6): each step at its calibrated transfer (shape k 0.65 / k_pds 0.42, composition k 1.3, key k 1.31),
same draws. Secondary for M5->M6: 522's own route (0.0710 k r ret + 0.1245 dpds_lit, unlevelled, k 1).

    uv run python aggregate534.py [--panels S1,S2,S3,S4]     -> work/partB/results534B.json, work/partB/tables534B.md
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard534  # noqa: E402,F401
import argparse  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import designs534 as D  # noqa: E402
import price534 as P5  # noqa: E402

HERE = P5.HERE
NB, SEED = 4000, 534


def complete(p):
    return (HERE / p / "full534.json").exists() and all(
        (HERE / p / f"comp534-{f}.json").exists() and "peak_rss_gib" in json.loads((HERE / p / f"comp534-{f}.json").read_text())
        for f in P5.CANCER)


def route522(P, f):
    """522's route for M5 -> M6 on fold f: r = cos_agg(M6)/cos_agg(M5) - 1 (unlevelled, realised 16), dpds literal."""
    a, b = P.stats(f, D.M["M6"]), P.stats(f, D.M["M5"])
    r = a["cos"] / b["cos"] - 1
    dp = a["pds"] - b["pds"]
    return dict(r=r, dpds=dp, ret=P.ret[f], k1=0.0710 * r * P.ret[f] + 0.1245 * dp, k131=0.0710 * 1.31 * r * P.ret[f] + 0.1245 * dp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", default="S1,S2,S3,S4")
    ap.add_argument("--out", default="534B", help="output stem: results<stem>.json, tables<stem>.md")
    ap.add_argument("--nb", type=int, default=NB)
    a = ap.parse_args()
    t0 = time.time()
    want = [p for p in a.panels.split(",") if p]
    priced = [p for p in want if complete(p)]
    assert priced, "no R6 panel is priced yet"
    P = {p: P5.Panel(p) for p in priced}
    rec = json.loads((HERE / "record534.json").read_text())
    txt, stored = rec["text"], rec["stored"]

    def abc_record(sid):
        return stored[sid]["abc"]["mean"] if sid in stored else txt[sid]["abc"]

    def oct_record(sid, t):
        return stored[sid][t]["mean"] if sid in stored else txt[sid][t]

    prim = P5.point_table(P, level=True)
    cal = P5.point_table(P, level=True, cal=True)
    unl = P5.point_table(P, level=False)
    bp, bpp = P5.pooled_bootstrap(P, a.nb, SEED, level=True)
    bc, bcp = P5.pooled_bootstrap(P, a.nb, SEED, level=True, cal=True)
    efp = P5.Panel("efp") if (HERE / "efp" / "full534.json").exists() else None
    efp_tab = P5.point_table({"efp": efp}, level=True)["efp"] if efp else None

    steps = {}
    for st in D.STEPS:
        sid = st[0]
        pooled = prim["pooled"][sid]
        ci = P5.summ(bp[sid])
        d, lb = pooled["mean"], ci["ci95"][0]
        verdict = "TRANSFERS" if (d > 0 and lb > 0) else ("FRAGILE" if d > 0 else "FAILS")
        abc = abc_record(sid)
        row = dict(adds=st[3], cls=st[4], arm=D.M[st[1]], ref=D.M[st[2]],
                   record=dict(abc=abc, real=oct_record(sid, "real"), efp=oct_record(sid, "efp"), efr=oct_record(sid, "efr"),
                               board=txt[sid]["board"], src=txt[sid]["src"], exact_stored=sid in stored),
                   efp_exact=efp_tab[sid] if efp_tab else None,
                   panels={p: dict(mean=prim[p][sid]["mean"], n_pos=prim[p][sid]["n_pos"], per_fold=prim[p][sid]["per_fold"],
                                   ci95=P5.summ(bpp[p][sid])["ci95"]) for p in priced},
                   pooled=dict(mean=d, n_pos=pooled["n_pos"], n_folds=pooled["n_folds"], **ci),
                   verdict=verdict if st in D.MAIN_STEPS else None,
                   unlevelled_pooled=unl["pooled"][sid]["mean"],
                   calibrated=dict(k=D.K_CAL[st[4]], pooled=cal["pooled"][sid]["mean"], **P5.summ(bc[sid]),
                                   panels={p: cal[p][sid]["mean"] for p in priced}))
        if abc:
            rr = np.asarray(bp[sid]) / abc
            row["retention"] = dict(point=d / abc, ci95=[float(np.quantile(rr, 0.025)), float(np.quantile(rr, 0.975))])
        steps[sid] = row

    # cumulative M0 -> M6 (the sum of the six milestone steps; the same draws)
    main_ids = [s[0] for s in D.MAIN_STEPS]
    abc_cum = sum(abc_record(s) for s in main_ids)
    cum = dict(pooled=float(sum(steps[s]["pooled"]["mean"] for s in main_ids)), **P5.summ(bp["cum"]),
               panels={p: float(sum(steps[s]["panels"][p]["mean"] for s in main_ids)) for p in priced},
               abc_record=abc_cum, record_sum={t: float(sum(oct_record(s, t) for s in main_ids)) for t in ("real", "efp", "efr")},
               calibrated=dict(pooled=float(sum(steps[s]["calibrated"]["pooled"] for s in main_ids)), **P5.summ(bc["cum"])),
               direct_M6_vs_M0=dict(pooled=float(np.mean([P[p].price(f, D.M["M6"], 1.0, D.M["M0"], True)["ov"]
                                                          for p in priced for f in P[p].folds]))))
    rr = np.asarray(bp["cum"]) / abc_cum
    cum["retention"] = dict(point=cum["pooled"] / abc_cum, ci95=[float(np.quantile(rr, 0.025)), float(np.quantile(rr, 0.975))])
    board_sum = {s: steps[s]["record"]["board"] for s in main_ids}
    cum["board_known"] = dict(steps=board_sum, sum_M0_M5=float(sum(v for s, v in board_sum.items() if v is not None)),
                              note="board reads of the bulk steps; M1->M2 (the H1 half alone) and M5->M6 (495, 10-03) have no board read")

    # section 7 rule 3: the tuned-on check
    tuned = {}
    for s in D.OCTOBER_SWITCHES:
        m = float(np.mean([oct_record(s, t) for t in ("real", "efp", "efr")]))
        tuned[s] = dict(r6=steps[s]["pooled"]["mean"], mean_real_efp_efr=m, half=0.5 * m,
                        below_half=bool(steps[s]["pooled"]["mean"] < 0.5 * m))
    over = sum(v["below_half"] for v in tuned.values())
    # rule 4: controls
    ca, cb = steps["C-a"]["pooled"], steps["C-b"]["pooled"]
    controls = dict(Ca=dict(mean=ca["mean"], ci95=ca["ci95"], ok=ca["mean"] <= 0,
                            flag_cd4=bool(ca["mean"] > 0 and ca["ci95"][0] > 0)),
                    Cb=dict(mean=cb["mean"], ci95=cb["ci95"], ok=cb["mean"] <= 0))
    r522 = {p: {f: route522(P[p], f) for f in P[p].folds} for p in priced}
    all522 = [v for p in priced for v in r522[p].values()]
    secondary = dict(M5_M6_route522=dict(pooled_k1=float(np.mean([v["k1"] for v in all522])),
                                         pooled_k131=float(np.mean([v["k131"] for v in all522])),
                                         pooled_r=float(np.mean([v["r"] for v in all522])),
                                         per_panel={p: dict(k1=float(np.mean([v["k1"] for v in r522[p].values()])),
                                                            r=float(np.mean([v["r"] for v in r522[p].values()])))
                                                    for p in priced}, per_fold=r522))
    # ---- secondary: levelling sensitivity. 438's route levels each design with ONE ntclevel step, which lands short of
    # rho* when the factor is large (438 section 2.3: <= 0.4 % short at factors <= 1.07, ~1.2 % at >= 1.09). The October
    # re-level is a secant loop (448) that converges. Energy-only extrapolation to rho*: per design and context,
    # t = sqrt(rho* claim / B_level); A_level, B_level x t^2; the amplitude scale x t; PLOLO statistics held fixed.
    import copy
    landed, conv_tab = {}, {}
    for p in priced:
        lv = json.loads((HERE / p / "level534.json").read_text())
        rho_t = lv["rho_target"]
        landed[p] = {n: {c: P[p].E[n]["rho_level"][c] for c in P5.SLOTS} for n in P[p].E}
        Q = copy.copy(P[p])
        Q.E = copy.deepcopy(P[p].E)
        for n, e in Q.E.items():
            for c in P5.SLOTS:
                tt = (rho_t * P[p].claims[c] / e["AB_level"][c]["B"]) ** 0.5
                e["AB_level"][c] = dict(A=e["AB_level"][c]["A"] * tt * tt, B=e["AB_level"][c]["B"] * tt * tt)
        per_fold = {f: P5.step_prices(Q, f, level=True) for f in Q.folds}
        conv_tab[p] = {st[0]: float(np.mean([per_fold[f][st[0]] for f in Q.folds])) for st in D.STEPS}
    conv_pooled = {st[0]: float(np.mean([conv_tab[p][st[0]] for p in priced])) for st in D.STEPS}
    secondary["levelling_sensitivity"] = dict(
        note="energy-only extrapolation of every design to rho* (the converged October re-level); PLOLO statistics fixed",
        pooled=conv_pooled, per_panel=conv_tab, cum_pooled=float(sum(conv_pooled[s] for s in main_ids)),
        landed_rho=landed)
    # ---- supplementary (not preregistered): the exact milestone designs on A/B/C's own panel, this instrument
    if complete("abc"):
        PA = P5.Panel("abc")
        at = P5.point_table({"abc": PA}, level=True)["abc"]
        sup = {}
        for st in D.STEPS:
            v = at[st[0]]["mean"]
            rr = np.asarray(bp[st[0]]) / v
            sup[st[0]] = dict(abc_exact=v, n_pos=at[st[0]]["n_pos"], per_fold=at[st[0]]["per_fold"],
                              retention=dict(point=steps[st[0]]["pooled"]["mean"] / v,
                                             ci95=[float(np.quantile(rr, 0.025)), float(np.quantile(rr, 0.975))]))
        ac = float(sum(sup[s]["abc_exact"] for s in main_ids))
        rr = np.asarray(bp["cum"]) / ac
        sup["cum"] = dict(abc_exact=ac, retention=dict(point=cum["pooled"] / ac,
                                                       ci95=[float(np.quantile(rr, 0.025)), float(np.quantile(rr, 0.975))]))
        secondary["abc_exact"] = sup
    surprise = [s for s in ("M2->M3", "M3->M4", "M4->M5", "M5->M6") if steps[s]["pooled"]["mean"] <= 0]
    res = dict(written=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), panels_priced=priced, n_boot=a.nb, seed=SEED,
               route="438 LEVELLED, k = 1 (price534 == price438 on efp, KA-P)", steps=steps, cumulative=cum,
               tuned_on=dict(per_switch=tuned, n_below_half=int(over), OVER_READ=bool(over >= 2)), controls=controls,
               secondary=secondary, ret={p: P[p].ret for p in priced},
               surprise=dict(sign_flips_M3_M6=[s for s in surprise if s != "M2->M3"],
                             cum_retention_below_0_3=bool(cum["retention"]["point"] < 0.3)),
               seconds=time.time() - t0)
    (HERE / f"results{a.out}.json").write_text(json.dumps(res, indent=1, default=float))
    tables(res, a.out)
    print(f"aggregate: panels {priced}; {time.time() - t0:.0f} s")


def f5(x):
    return "—" if x is None else f"{x:+.5f}"


def tables(res, out="534B"):
    L = []
    pr = res["panels_priced"]
    L.append(f"Panels priced: {', '.join(pr)}. Levelled k = 1, cancer-fold means; pooled = mean over {3 * len(pr)} folds; "
             f"95 % intervals from {res['n_boot']} paired target-bootstrap draws (seed {res['seed']}).\n")
    L.append("| step | adds | A/B/C (record) | real (record) | efp (record) | efr (record) | efp, exact design (this instrument) | "
             + " | ".join(pr) + " | **R6 pooled** [95 % CI] (n+) | verdict | retention [95 % CI] |")
    L.append("|---|---|---:|---:|---:|---:|---:|" + "---:|" * len(pr) + "---|---|---|")
    for sid, s in res["steps"].items():
        r = s["record"]
        ex = s["efp_exact"]
        cells = [f"{s['panels'][p]['mean']:+.5f} ({s['panels'][p]['n_pos']}/3)" for p in pr]
        po = s["pooled"]
        ret = s.get("retention")
        L.append(f"| {sid} | {s['adds']} | {f5(r['abc'])} | {f5(r['real'])} | {f5(r['efp'])} | {f5(r['efr'])} | "
                 + (f"{ex['mean']:+.5f} ({ex['n_pos']}/3)" if ex else "—") + " | " + " | ".join(cells)
                 + f" | **{po['mean']:+.5f}** [{po['ci95'][0]:+.5f}, {po['ci95'][1]:+.5f}] ({po['n_pos']}/{po['n_folds']}) | "
                 + f"{s['verdict'] or 'control'} | "
                 + (f"{ret['point']:.2f} [{ret['ci95'][0]:.2f}, {ret['ci95'][1]:.2f}]" if ret else "—") + " |")
    c = res["cumulative"]
    L.append(f"| **M0 → M6** | the six steps summed | {c['abc_record']:+.5f} | {c['record_sum']['real']:+.5f} | "
             f"{c['record_sum']['efp']:+.5f} | {c['record_sum']['efr']:+.5f} | — | "
             + " | ".join(f"{c['panels'][p]:+.5f}" for p in pr)
             + f" | **{c['pooled']:+.5f}** [{c['ci95'][0]:+.5f}, {c['ci95'][1]:+.5f}] | — | "
             + f"**{c['retention']['point']:.2f}** [{c['retention']['ci95'][0]:.2f}, {c['retention']['ci95'][1]:.2f}] |")
    L.append("\nCalibrated board-unit projection (shape k 0.65 / k_pds 0.42; composition k 1.3; key k 1.31), R6 pooled:\n")
    L.append("| step | k (X, pds) | " + " | ".join(pr) + " | pooled [95 % CI] | A/B/C board read |")
    L.append("|---|---|" + "---:|" * len(pr) + "---|---:|")
    for sid, s in res["steps"].items():
        k = s["calibrated"]
        L.append(f"| {sid} | {k['k'][0]:g}, {k['k'][1]:g} | " + " | ".join(f"{k['panels'][p]:+.5f}" for p in pr)
                 + f" | {k['pooled']:+.5f} [{k['ci95'][0]:+.5f}, {k['ci95'][1]:+.5f}] | {f5(s['record']['board'])} |")
    L.append("| **M0 → M6** | | " + " | ".join(f"{sum(res['steps'][s]['calibrated']['panels'][p] for s in ('M0->M1', 'M1->M2', 'M2->M3', 'M3->M4', 'M4->M5', 'M5->M6')):+.5f}" for p in pr)
             + f" | **{c['calibrated']['pooled']:+.5f}** [{c['calibrated']['ci95'][0]:+.5f}, {c['calibrated']['ci95'][1]:+.5f}] | "
             + f"{c['board_known']['sum_M0_M5']:+.5f} (the four direct reads) |")
    (HERE / f"tables{out}.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
