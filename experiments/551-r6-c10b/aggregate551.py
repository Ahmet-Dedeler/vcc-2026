"""551 (C10b) -- the claim's aggregation and frozen reading rules (PREREGISTRATION.md sections 5-7). Light (jsons + arrays).

Reads the configuration F recorded BEFORE the access (work/c10b-config.json: the coordinator's message verbatim, its
timestamp, the resolved designs551.CONFIGS key), then every R6 panel priced by comp551 --repro --config F.

  1. KA-REPRO (pooled): C10a's six milestone steps and two controls, re-priced from the re-fitted operators, equal C10a's
     pooled point values in experiments/534-generalization-audit/work/partB/results534B.json to |d| <= 1e-12. Abort otherwise.
  2. The claim: F's own step vs its base (S_P2@M6, P2g@M6, S_P2@M5 or P2g@M5), 438's LEVELLED route at k = 1, the mean over
     R6's 12 folds; paired target bootstrap, 4,000 draws, seed 551. Verdict: TRANSFERS if D > 0 and the 95 % lower bound
     > 0; FRAGILE if D > 0 and lb <= 0; FAILS if D <= 0.
  3. Its control: the same per-target q permuted over targets (seed 493), priced vs the same base, same draws. Reported,
     with the paired difference (claim - control) and its interval: the value of the TARGETING.
  4. Cumulative M0 -> F = C10a's six steps + F's step (same draws); retention vs A/B/C (C10a's record values for the six
     steps; F's step on A/B/C's own panel, this instrument, the abc run), and the calibrated board-unit projection.
  5. Consequence: FAILS -> the coordinator gets `sharpen_pt` as a candidate to switch off for October; FRAGILE -> flagged;
     TRANSFERS -> nothing. Nothing is changed here.
If F is a C10a design (M5 or M6: no per-target step), C10b does not access R6: the script says so and exits before reading.
Only the configuration's own steps are priced and printed; the other per-target candidates are never computed on R6.

    uv run python aggregate551.py [--panels S1,S2,S3,S4] [--nb 4000]   -> work/results551.json, work/tables551.md
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard551  # noqa: E402,F401
import argparse  # noqa: E402
import copy  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import designs551 as D  # noqa: E402
import price551 as P5  # noqa: E402

HERE, WORK, ROOT = P5.HERE, P5.WORK, P5.ROOT
W534 = ROOT / "experiments/534-generalization-audit/work/partB"
NB, SEED = 4000, 551
TOL = 1e-12


def config():
    c = json.loads((WORK / "c10b-config.json").read_text())
    assert c["config"] in D.CONFIGS, c["config"]
    return c


def claim_steps(cfg):
    base, sid = D.CONFIGS[cfg]
    if sid is None:
        return base, None, None
    st = {s[0]: s for s in D.STEPS_C10B}[sid]
    ctrl_id = ("C-S@" if sid.startswith("S_P2") else "C-P@") + base
    return base, st, {s[0]: s for s in D.STEPS_C10B}[ctrl_id]


def complete(p):
    return (WORK / p / "full551.json").exists() and all(
        (WORK / p / f"comp551-{f}.json").exists() and "peak_rss_gib" in json.loads((WORK / p / f"comp551-{f}.json").read_text())
        for f in P5.CANCER)


def verdict(d, lb):
    return "TRANSFERS" if (d > 0 and lb > 0) else ("FRAGILE" if d > 0 else "FAILS")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", default="S1,S2,S3,S4")
    ap.add_argument("--nb", type=int, default=NB)
    a = ap.parse_args()
    t0 = time.time()
    cfg = config()
    base, st, ctrl = claim_steps(cfg["config"])
    if st is None:
        print(f"C10b: F = {cfg['config']} is a C10a design; its R6 value is C10a's. R6 is NOT accessed. Nothing read.")
        return
    want = [p for p in a.panels.split(",") if p]
    priced = [p for p in want if complete(p)]
    assert priced == want, f"not every R6 panel is priced: {priced} of {want}"
    for p in priced:                                          # every fold ran with the recorded configuration and KA-REPRO
        for f in P5.CANCER:
            j = json.loads((WORK / p / f"comp551-{f}.json").read_text())
            assert j.get("config") == cfg["config"], (p, f, j.get("config"))
            assert j["ka"].get("repro", {}).get("worst", 1.0) <= TOL, (p, f, j["ka"].get("repro"))
    P = {p: P5.Panel(p) for p in priced}

    # 1. KA-REPRO (pooled): C10a's steps reproduce C10a's pooled point values exactly
    c10a = json.loads((W534 / "results534B.json").read_text())
    prim_a = P5.point_table(P, level=True, steps=D.STEPS_C10A)
    worst = max(abs(prim_a["pooled"][s[0]]["mean"] - c10a["steps"][s[0]]["pooled"]["mean"]) for s in D.STEPS_C10A)
    print(f"KA-REPRO pooled: C10a's 8 steps re-priced from the re-fitted operators, worst |d| {worst:.1e}")
    assert worst <= TOL, f"KA-REPRO pooled FAILED: {worst}"

    # 2-4. the claim, its control, the cumulative (one set of draws)
    steps = D.STEPS_C10A[:6] + (st, ctrl)
    cum_ids = tuple(s[0] for s in D.STEPS_C10A[:6]) + (st[0],)
    prim = P5.point_table(P, level=True, steps=steps)
    cal = P5.point_table(P, level=True, cal=True, steps=steps)
    bp, bpp = P5.pooled_bootstrap(P, a.nb, SEED, level=True, steps=steps, cum_ids=cum_ids)
    bc, _ = P5.pooled_bootstrap(P, a.nb, SEED, level=True, cal=True, steps=steps, cum_ids=cum_ids)
    rec = json.loads((W534 / "record534.json").read_text())
    txt, stored = rec["text"], rec["stored"]
    abc_rec = {s: (stored[s]["abc"]["mean"] if s in stored else txt[s]["abc"]) for s in cum_ids[:6]}
    abc_new = None
    if complete("abc"):                                       # F's step on A/B/C's own panel (tuned-on; not an R6 access)
        abc_new = P5.point_table({"abc": P5.Panel("abc")}, level=True, steps=(st,))["abc"][st[0]]["mean"]

    def row(s):
        sid = s[0]
        po = prim["pooled"][sid]
        ci = P5.summ(bp[sid])
        return dict(adds=s[3], cls=s[4], arm=D.M[s[1]], ref=D.M[s[2]],
                    panels={p: dict(mean=prim[p][sid]["mean"], n_pos=prim[p][sid]["n_pos"], per_fold=prim[p][sid]["per_fold"],
                                    ci95=P5.summ(bpp[p][sid])["ci95"]) for p in priced},
                    pooled=dict(mean=po["mean"], n_pos=po["n_pos"], n_folds=po["n_folds"], **ci),
                    calibrated=dict(k=D.K_CAL[s[4]], pooled=cal["pooled"][sid]["mean"], **P5.summ(bc[sid])))

    claim = row(st)
    claim["verdict"] = verdict(claim["pooled"]["mean"], claim["pooled"]["ci95"][0])
    if abc_new:
        rr = np.asarray(bp[st[0]]) / abc_new
        claim["abc_this_instrument"] = abc_new
        claim["retention"] = dict(point=claim["pooled"]["mean"] / abc_new,
                                  ci95=[float(np.quantile(rr, 0.025)), float(np.quantile(rr, 0.975))])
    control = row(ctrl)
    diff = np.asarray(bp[st[0]]) - np.asarray(bp[ctrl[0]])
    targeting = dict(point=claim["pooled"]["mean"] - control["pooled"]["mean"], **P5.summ(diff))
    abc_cum = sum(abc_rec.values()) + (abc_new or 0.0)
    cum = dict(pooled=float(sum(prim["pooled"][s]["mean"] for s in cum_ids)), **P5.summ(bp["cum"]),
               panels={p: float(sum(prim[p][s]["mean"] for s in cum_ids)) for p in priced},
               abc=abc_cum if abc_new else None, c10a_M0_M6=c10a["cumulative"]["pooled"],
               calibrated=dict(pooled=float(sum(cal["pooled"][s]["mean"] for s in cum_ids)), **P5.summ(bc["cum"])))
    if abc_new:
        rr = np.asarray(bp["cum"]) / abc_cum
        cum["retention"] = dict(point=cum["pooled"] / abc_cum, ci95=[float(np.quantile(rr, 0.025)), float(np.quantile(rr, 0.975))])

    # secondary: levelling sensitivity for the claim step (534's energy-only extrapolation to rho*)
    conv = {}
    for p in priced:
        lv = json.loads((WORK / p / "level551.json").read_text())
        Q = copy.copy(P[p])
        Q.E = copy.deepcopy(P[p].E)
        for n, e in Q.E.items():
            for c in P5.SLOTS:
                tt = (lv["rho_target"] * P[p].claims[c] / e["AB_level"][c]["B"]) ** 0.5
                e["AB_level"][c] = dict(A=e["AB_level"][c]["A"] * tt * tt, B=e["AB_level"][c]["B"] * tt * tt)
        conv[p] = float(np.mean([P5.step_prices(Q, f, level=True, steps=(st,))[st[0]] for f in Q.folds]))
    res = dict(written=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), config=cfg, panels_priced=priced, n_boot=a.nb,
               seed=SEED, route="438 LEVELLED, k = 1 (price551 == price534 == price438)",
               ka_repro_pooled_worst=worst, claim=claim, control=control, targeting=targeting, cumulative=cum,
               secondary=dict(levelling_sensitivity=dict(per_panel=conv, pooled=float(np.mean(list(conv.values()))))),
               consequence={"TRANSFERS": "nothing changes", "FRAGILE": "flag sharpen_pt to the coordinator",
                            "FAILS": "sharpen_pt goes to the coordinator as a candidate to switch off"}[claim["verdict"]],
               seconds=time.time() - t0)
    (WORK / "results551.json").write_text(json.dumps(res, indent=1, default=float))
    L = [f"C10b, F = {cfg['config']} (recorded {cfg['recorded_at']}); R6 panels {', '.join(priced)}; levelled k 1; pooled over "
         f"{3 * len(priced)} folds; {a.nb} paired target-bootstrap draws, seed {SEED}. KA-REPRO pooled worst |d| {worst:.1e}.\n",
         "| step | adds | " + " | ".join(priced) + " | **R6 pooled** [95 % CI] (n+) | verdict | calibrated [95 % CI] |",
         "|---|---|" + "---:|" * len(priced) + "---|---|---|"]
    for s, r in ((st, claim), (ctrl, control)):
        po, k = r["pooled"], r["calibrated"]
        L.append(f"| {s[0]} | {s[3]} | " + " | ".join(f"{r['panels'][p]['mean']:+.5f} ({r['panels'][p]['n_pos']}/3)" for p in priced)
                 + f" | **{po['mean']:+.5f}** [{po['ci95'][0]:+.5f}, {po['ci95'][1]:+.5f}] ({po['n_pos']}/{po['n_folds']}) | "
                 + f"{r.get('verdict', 'control')} | {k['pooled']:+.5f} [{k['ci95'][0]:+.5f}, {k['ci95'][1]:+.5f}] |")
    L.append(f"\nTargeting (claim − control, paired): {targeting['point']:+.5f} [{targeting['ci95'][0]:+.5f}, {targeting['ci95'][1]:+.5f}].")
    L.append(f"Cumulative M0 → F: {cum['pooled']:+.5f} [{cum['ci95'][0]:+.5f}, {cum['ci95'][1]:+.5f}] (C10a's M0 → M6 {cum['c10a_M0_M6']:+.5f})"
             + (f"; retention vs A/B/C {cum['retention']['point']:.2f} [{cum['retention']['ci95'][0]:.2f}, {cum['retention']['ci95'][1]:.2f}]" if abc_new else "")
             + f"; calibrated {cum['calibrated']['pooled']:+.5f} [{cum['calibrated']['ci95'][0]:+.5f}, {cum['calibrated']['ci95'][1]:+.5f}].")
    L.append(f"Consequence: {res['consequence']}.")
    (WORK / "tables551.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
