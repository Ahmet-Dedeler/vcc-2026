"""534 Part C -- aggregation and the frozen reading rules (PREREGISTRATION-C section 4). Light (jsons + arrays).

  1. Step value Delta_s = B0 - B_s (price438's levelled route, k = 1, arm B0 vs ref B_s), per panel as the mean of its
     three folds, and pooled over real / efp / efr (9 folds) with a paired target bootstrap: 4,000 draws, seed 534.
  2. Verdict on the October-like pool: HOLDS if Delta > 0 and the 95 % lower bound > 0; FRAGILE if Delta > 0 and the lower
     bound <= 0; FAILS if Delta <= 0.
  3. The abc row beside the pool: if its sign disagrees with the board's positive read, the instrument disagrees with the
     board for this class (reported; no conclusion drawn).
  4. Consequence: FAILS with a positive board read -> A/B/C-specific, a candidate to switch off for October; FRAGILE ->
     flagged. Nothing in production changes here.
abc's interval: its own paired bootstrap (4,000 draws, a separate Generator seeded 534). Secondary: the unlevelled route
and each step's decomposition (relative cosine, board-emulated pds, raw mse).

    uv run python aggregate534c.py        -> work/partC/results534C.json, work/partC/tables534C.md
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard534  # noqa: E402,F401
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import designs534c as D  # noqa: E402
import price534c as P5  # noqa: E402

HERE = P5.HERE
NB, SEED = 4000, 534


def summ(x):
    x = np.asarray(x)
    return dict(ci95=[float(np.quantile(x, 0.025)), float(np.quantile(x, 0.975))], p_le0=float(np.mean(x <= 0)), sd=float(np.std(x)))


def main():
    t0 = time.time()
    P = {p: P5.Panel(p) for p in D.PANELS}
    pt = {p: P5.point(P[p], level=True) for p in D.PANELS}
    unl = {p: P5.point(P[p], level=False) for p in D.PANELS}
    decomp = {p: {sid: {k: float(np.mean([P[p].price(f, arm, 1.0, ref, True)[k] for f in P[p].folds]))
                        for k in ("d_rel", "dpds_bd", "dpds_lit", "dmse")} for sid, arm, ref, *_ in D.STEPS} for p in D.PANELS}
    boot_oct = P5.bootstrap([P[p] for p in D.OCT], NB, SEED)
    boot_abc = P5.bootstrap([P["abc"]], NB, SEED)
    steps = {}
    for sid, arm, ref, what, board, bsrc in D.STEPS:
        folds = [pt[p][sid]["per_fold"][f] for p in D.OCT for f in P[p].folds]
        d = float(np.mean(folds))
        ci = summ(boot_oct[sid])
        lb = ci["ci95"][0]
        verdict = "HOLDS" if (d > 0 and lb > 0) else ("FRAGILE" if d > 0 else "FAILS")
        abc = pt["abc"][sid]["mean"]
        steps[sid] = dict(
            what=what, arm=arm, ref=ref, design=f"{D.COMP}|{D.KEY}+sh{D.Q:g}", toggles_arm=D.BASES[arm], toggles_ref=D.BASES[ref],
            pooled=dict(mean=d, n_pos=int(sum(x > 0 for x in folds)), n_folds=len(folds), **ci), verdict=verdict,
            panels={p: dict(mean=pt[p][sid]["mean"], n_pos=pt[p][sid]["n_pos"], per_fold=pt[p][sid]["per_fold"],
                            unlevelled=unl[p][sid]["mean"], decomposition=decomp[p][sid]) for p in D.PANELS},
            abc=dict(mean=abc, n_pos=pt["abc"][sid]["n_pos"], **summ(boot_abc[sid]),
                     agrees_with_board_sign=bool((abc > 0) == (board > 0))),
            board=dict(value=board, source=bsrc),
            unlevelled_pooled=float(np.mean([unl[p][sid]["per_fold"][f] for p in D.OCT for f in P[p].folds])),
            consequence=("A/B/C-specific: candidate to switch off for October (PREREGISTRATION-C 4.4)" if (verdict == "FAILS" and board > 0)
                         else "flagged (FRAGILE)" if verdict == "FRAGILE" else "none (HOLDS)"))
    # known answers, from the run jsons
    ka = {}
    for p in D.PANELS:
        if p not in D.FOLD_ROUTE:
            fj = json.loads((HERE / p / "full534c.json").read_text())
            ka[f"{p}/full"] = dict(composition_vs_meta=fj["ka"]["composition_vs_meta"]["worst"], key=fj["ka"]["key410_vs_artifact"],
                                   energies=fj["ka"]["energies"], uniform=fj["states"]["B2"]["pre_equals_uniform"],
                                   peak_rss_gib=fj["peak_rss_gib"], seconds=fj["seconds"])
        else:
            lv = json.loads((HERE / p / "level534c.json").read_text())
            ka[f"{p}/foldlevel"] = dict(claims_vs_438_efr_s1=lv["claims_vs_438_efr_s1"], factors_vs_438=lv["ka"])
        for f in D.CANCER:
            cj = json.loads((HERE / p / f"comp534c-{f}.json").read_text())
            ka[f"{p}/{f}"] = dict(n=len(cj["ka"]["rows"]), worst=cj["ka"]["worst"], fold=cj["ka"]["fold"],
                                  rows=cj["ka"]["rows"], peak_rss_gib=cj["peak_rss_gib"], seconds=cj["seconds"],
                                  states={t: dict(pre_rel_diff_vs_production=s.get("pre_rel_diff_vs_production"),
                                                  key_frac_floor=s["key_frac_floor"], pre_equals_uniform=s["pre_equals_uniform"])
                                          for t, s in cj["states"].items()})
    res = dict(written=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), n_boot=NB, seed=SEED,
               route="438 LEVELLED, k = 1; efr on the fold route (438 section 2.5); price534c == price438.Panel.price",
               steps=steps, ka=ka, ret={p: P[p].ret for p in D.PANELS}, energy_route={p: P[p].route for p in D.PANELS},
               seconds=time.time() - t0)
    (HERE / "results534C.json").write_text(json.dumps(res, indent=1, default=float))
    tables(res)
    print(f"aggregate534c: {time.time() - t0:.0f} s")


def f5(x):
    return f"{x:+.5f}"


def tables(res):
    L = ["| step | Δ pooled real/efp/efr (k = 1, board units) | 95 % CI | folds positive (of 9) | verdict | Δ on `abc` [95 % CI] | the board's A/B/C read |",
         "|---|---:|---|---|---|---|---|"]
    for sid, s in res["steps"].items():
        po, a = s["pooled"], s["abc"]
        L.append(f"| {sid}: {s['what']} (B0 − {s['ref']}) | **{f5(po['mean'])}** | [{f5(po['ci95'][0])}, {f5(po['ci95'][1])}] | "
                 f"{po['n_pos']}/9 | **{s['verdict']}** | {f5(a['mean'])} [{f5(a['ci95'][0])}, {f5(a['ci95'][1])}] ({a['n_pos']}/3) | "
                 f"{f5(s['board']['value'])} ({s['board']['source']}) |")
    L.append("\n| step | real | efp | efr | abc | pooled, unlevelled |")
    L.append("|---|---:|---:|---:|---:|---:|")
    for sid, s in res["steps"].items():
        L.append(f"| {sid} | " + " | ".join(f"{f5(s['panels'][p]['mean'])} ({s['panels'][p]['n_pos']}/3)" for p in ("real", "efp", "efr", "abc"))
                 + f" | {f5(s['unlevelled_pooled'])} |")
    (HERE / "tables534C.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
