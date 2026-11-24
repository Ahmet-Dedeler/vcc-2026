"""R7a's claim aggregation (552's reserve, panels T1-T4; the claim is preregistered in 552 and recorded by r7a551.py).
Light (jsons + per-target arrays). 551's pricer unchanged (price551 == price534 == price438: LEVELLED route, k = 1).

Reads work/r7a/claim.json (the bulk configuration F, recorded before the access) and the four T panels priced by
comp551 --config F (no --repro: R7a has no earlier values to reproduce; the instrument's identity with the efp-KA'd code is
proven by r7a551.py plan / record / access instead).

  PRIMARY    F against M5, the bulk recipe before the switches decided after C10a: the sum of F's post-C10a steps
             (M5 -> M6 if F's base is M6; F's per-target step on its base, if any), mean over R7a's 12 folds; paired target
             bootstrap, 4,000 draws, seed 552. Verdict: TRANSFERS if D > 0 and the 95 % lower bound > 0; FRAGILE if D > 0
             and lb <= 0; FAILS if D <= 0. If F = M5 (no post-C10a step), the primary row is the cumulative M0 -> M5.
  STEPS      each post-C10a step on its own, and the per-target step's control (its q permuted over targets, seed 493):
             the targeting value (step - control, paired).
  SECONDARY  the cumulative October bulk recipe M0 -> F (C10a's milestones up to F's base + F's per-target step), with its
             retention against R6 (C10a's pooled steps; C10b's step if it has been read) and against A/B/C (534's record;
             F's step on A/B/C's own panel if 551's abc run exists), and the calibrated board-unit projection (K_CAL).

    uv run python experiments/551-r6-c10b/aggregate_r7a551.py                         -> work/r7a/results_r7a.json, tables_r7a.md
    uv run python experiments/551-r6-c10b/aggregate_r7a551.py --test-config M6+P2g --panels efp --nb 200   (smoke test on
        efp's known-answer outputs; writes work/r7a/test/ only, never a T panel)
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard551  # noqa: E402,F401
import argparse  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import designs551 as D  # noqa: E402
import price551 as P5  # noqa: E402

HERE, WORK, ROOT = P5.HERE, P5.WORK, P5.ROOT
R7A = WORK / "r7a"
W534 = ROOT / "experiments/534-generalization-audit/work/partB"
NB, SEED = 4000, 552
C10A = {s[0]: s for s in D.STEPS_C10A}
C10B = {s[0]: s for s in D.STEPS_C10B}
MILESTONES = tuple(s[0] for s in D.STEPS_C10A[:6])          # M0->M1 ... M5->M6


def plan_steps(F: str):
    """-> (post-C10a step ids, the per-target step or None, its control or None, the cumulative M0 -> F ids)."""
    base, sid = D.CONFIGS[F]
    post = ["M5->M6"] if base == "M6" else []
    st = ctrl = None
    if sid is not None:
        st = C10B[sid]
        ctrl = C10B[("C-S@" if sid.startswith("S_P2") else "C-P@") + base]
        post.append(sid)
    upto = MILESTONES[:6] if base == "M6" else MILESTONES[:5]
    cum = tuple(upto) + ((sid,) if sid else ())
    return tuple(post), st, ctrl, cum


def complete(p):
    return (WORK / p / "full551.json").exists() and all(
        (WORK / p / f"comp551-{f}.json").exists() and "peak_rss_gib" in json.loads((WORK / p / f"comp551-{f}.json").read_text())
        for f in P5.CANCER)


def verdict(d, lb):
    return "TRANSFERS" if (d > 0 and lb > 0) else ("FRAGILE" if d > 0 else "FAILS")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", default="T1,T2,T3,T4")
    ap.add_argument("--nb", type=int, default=NB)
    ap.add_argument("--test-config", default=None, help="smoke test only: a designs551.CONFIGS key, with --panels efp")
    a = ap.parse_args()
    t0 = time.time()
    if a.test_config:
        assert a.panels == "efp", "the smoke test runs on efp's known-answer outputs only"
        F, claim_rec, out_dir = a.test_config, None, R7A / "test"
    else:
        # 555 finding 7: the reading is made once, with the preregistered options only
        assert a.nb == NB and a.panels == "T1,T2,T3,T4", "the R7a reading uses the preregistered --nb 4000 and T1-T4 only"
        assert not (R7A / "results_r7a.json").exists(), "results_r7a.json exists: the R7a reading is made once"
        claim_rec = json.loads((R7A / "claim.json").read_text())
        F, out_dir = claim_rec["config"], R7A
    assert F in D.CONFIGS, F
    post, st, ctrl, cum_ids = plan_steps(F)
    primary_ids = post if post else MILESTONES[:5]
    primary_what = "F vs M5 (the post-C10a switches)" if post else "cumulative M0 -> M5 (F = M5: no post-C10a switch)"
    want = [p for p in a.panels.split(",") if p]
    priced = [p for p in want if complete(p)]
    assert priced == want, f"not every panel is priced: {priced} of {want}"
    if not a.test_config:
        for p in priced:                                      # every fold ran with the recorded configuration, no repro
            for f in P5.CANCER:
                j = json.loads((WORK / p / f"comp551-{f}.json").read_text())
                assert j.get("config") == F and j.get("op_axis_checked"), (p, f, j.get("config"))
    P = {p: P5.Panel(p) for p in priced}

    steps = tuple(C10A[s] for s in MILESTONES) + ((st, ctrl) if st else ())
    prim = P5.point_table(P, level=True, steps=steps)
    cal = P5.point_table(P, level=True, cal=True, steps=steps)
    bp, bpp = P5.pooled_bootstrap(P, a.nb, SEED, level=True, steps=steps, cum_ids=tuple(primary_ids))
    bc, _ = P5.pooled_bootstrap(P, a.nb, SEED, level=True, cal=True, steps=steps, cum_ids=tuple(primary_ids))

    def draws_sum(B, ids):
        return np.sum([np.asarray(B[s]) for s in ids], axis=0)

    def row(s):
        sid = s[0]
        po = prim["pooled"][sid]
        return dict(adds=s[3], cls=s[4], arm=D.M[s[1]], ref=D.M[s[2]],
                    panels={p: dict(mean=prim[p][sid]["mean"], n_pos=prim[p][sid]["n_pos"], per_fold=prim[p][sid]["per_fold"],
                                    ci95=P5.summ(bpp[p][sid])["ci95"]) for p in priced},
                    pooled=dict(mean=po["mean"], n_pos=po["n_pos"], n_folds=po["n_folds"], **P5.summ(bp[sid])),
                    calibrated=dict(k=D.K_CAL[s[4]], pooled=cal["pooled"][sid]["mean"], **P5.summ(bc[sid])))

    def cum_row(ids):
        return dict(steps=list(ids), pooled=float(sum(prim["pooled"][s]["mean"] for s in ids)), **P5.summ(draws_sum(bp, ids)),
                    panels={p: float(sum(prim[p][s]["mean"] for s in ids)) for p in priced},
                    n_pos_folds=int(sum(v > 0 for v in np.sum([list(prim["pooled"][s]["per_fold"]) for s in ids], axis=0))),
                    calibrated=dict(pooled=float(sum(cal["pooled"][s]["mean"] for s in ids)), **P5.summ(draws_sum(bc, ids))))

    primary = cum_row(primary_ids)
    primary.update(what=primary_what, verdict=verdict(primary["pooled"], primary["ci95"][0]) if post else "NO-POST-C10A-SWITCH")
    rows = {sid: row(C10A[sid] if sid in C10A else C10B[sid]) for sid in post}
    targeting = None
    if st:
        rows[ctrl[0]] = row(ctrl)
        diff = np.asarray(bp[st[0]]) - np.asarray(bp[ctrl[0]])
        targeting = dict(point=rows[st[0]]["pooled"]["mean"] - rows[ctrl[0]]["pooled"]["mean"], **P5.summ(diff))
    secondary = cum_row(cum_ids)

    # retention against R6 (C10a's pooled milestone values; C10b's step once read) and against A/B/C (534's record; F's
    # step on A/B/C's own panel from 551's abc run, if it exists)
    c10a = json.loads((W534 / "results534B.json").read_text())
    r6 = {s: c10a["steps"][s]["pooled"]["mean"] for s in MILESTONES}
    c10b_path = WORK / "results551.json"
    if st and c10b_path.exists():
        c10b = json.loads(c10b_path.read_text())
        if c10b.get("config", {}).get("config") == F:
            r6[st[0]] = c10b["claim"]["pooled"]["mean"]
    rec = json.loads((W534 / "record534.json").read_text())
    txt, stored = rec["text"], rec["stored"]
    abc = {s: (stored[s]["abc"]["mean"] if s in stored else txt[s]["abc"]) for s in MILESTONES}
    abc = {s: v for s, v in abc.items() if v is not None}
    abc_src = {s: ("this instrument, levelled (record534 stored)" if s in stored else
                   "record534 text, NOT this instrument: " + str(txt[s].get("src", ""))[:120]) for s in abc}
    if st and complete("abc"):
        abc[st[0]] = P5.point_table({"abc": P5.Panel("abc")}, level=True, steps=(st,))["abc"][st[0]]["mean"]
        abc_src[st[0]] = "this instrument, levelled (551's abc run)"

    def retention(r, ids, ref):
        if not all(s in ref for s in ids):
            return None
        den = float(sum(ref[s] for s in ids))
        rr = draws_sum(bp, ids) / den
        return dict(reference=den, point=r["pooled"] / den, ci95=[float(np.quantile(rr, 0.025)), float(np.quantile(rr, 0.975))],
                    sources={s: (abc_src.get(s, "") if ref is abc else "R6, this instrument, levelled") for s in ids})
    for r, ids in ((primary, primary_ids), (secondary, cum_ids)):
        r["retention_vs_R6"] = retention(r, ids, r6)
        r["retention_vs_ABC"] = retention(r, ids, abc)

    res = dict(written=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), config=F, claim=claim_rec and
               {k: claim_rec[k] for k in ("claim", "claim_sha256", "recorded_at")}, test=bool(a.test_config),
               panels_priced=priced, n_boot=a.nb, seed=SEED, route="438 LEVELLED, k = 1 (price551 == price534 == price438)",
               primary=primary, steps=rows, targeting=targeting, secondary=secondary,
               consequence={"TRANSFERS": "nothing changes", "FRAGILE": "flag the post-C10a switches to the coordinator",
                            "FAILS": "the post-C10a switches go to the coordinator as candidates to switch off",
                            "NO-POST-C10A-SWITCH": "nothing to confirm (F = M5)"}[primary["verdict"]],
               seconds=time.time() - t0)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "results_r7a.json").write_text(json.dumps(res, indent=1, default=float))
    L = [f"R7a{' SMOKE TEST (not a reading)' if a.test_config else ''}, F = {F}; panels {', '.join(priced)}; levelled k 1; pooled "
         f"over {3 * len(priced)} folds; {a.nb} paired target-bootstrap draws, seed {SEED}.\n",
         f"**PRIMARY ({primary_what}):** {primary['pooled']:+.5f} [{primary['ci95'][0]:+.5f}, {primary['ci95'][1]:+.5f}] "
         f"({primary['n_pos_folds']}/{3 * len(priced)} folds > 0): **{primary['verdict']}**; calibrated {primary['calibrated']['pooled']:+.5f}.",
         "", "| step | adds | " + " | ".join(priced) + " | pooled [95 % CI] (n+) | calibrated |", "|---|---|" + "---:|" * len(priced) + "---|---|"]
    for sid, r in rows.items():
        po = r["pooled"]
        L.append(f"| {sid} | {r['adds']} | " + " | ".join(f"{r['panels'][p]['mean']:+.5f} ({r['panels'][p]['n_pos']}/3)" for p in priced)
                 + f" | {po['mean']:+.5f} [{po['ci95'][0]:+.5f}, {po['ci95'][1]:+.5f}] ({po['n_pos']}/{po['n_folds']}) | "
                 + f"{r['calibrated']['pooled']:+.5f} |")
    if targeting:
        L.append(f"\nTargeting ({st[0]} − its control, paired): {targeting['point']:+.5f} [{targeting['ci95'][0]:+.5f}, {targeting['ci95'][1]:+.5f}].")
    for name, r in (("Primary", primary), ("Cumulative M0 → F", secondary)):
        for lab, k in (("R6", "retention_vs_R6"), ("A/B/C", "retention_vs_ABC")):
            if r.get(k):
                L.append(f"{name}: retention vs {lab} {r[k]['point']:.2f} [{r[k]['ci95'][0]:.2f}, {r[k]['ci95'][1]:.2f}] "
                         f"({lab} {r[k]['reference']:+.5f}).")
    L.append(f"Cumulative M0 → F: {secondary['pooled']:+.5f} [{secondary['ci95'][0]:+.5f}, {secondary['ci95'][1]:+.5f}]; "
             f"calibrated {secondary['calibrated']['pooled']:+.5f}.")
    L.append(f"Consequence: {res['consequence']}.")
    (out_dir / "tables_r7a.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
