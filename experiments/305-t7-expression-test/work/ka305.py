"""305 step 3 -- the known-answer rows (PREREGISTRATION sec 10).  Run BEFORE any candidate conclusion.

    PYTHONDONTWRITEBYTECODE=1 uv run python ka305.py

KA-1 construct identities + 13/21 labels + cell counts   KA-5 chunked sums vs a direct row slice
KA-2 gene order and missingness                          KA-6 same-NTC-halves null
KA-3 301's 4/13, 4/21, 4/31 reproduced                   KA-7 34 own-target effects vs 293's kd_full
KA-4 synthetic known fold change

A failure aborts the candidate conclusions (`expr305.py` refuses to run).
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib305 as X  # noqa: E402

E293 = X.ROOT / "experiments/293-challenge-construct-ceiling/work"
E298 = X.ROOT / "experiments/298-reviewer-checks/work"
SEED305 = 20260907
KA7_TOL = 1e-9


def load(split):
    z = np.load(X.HERE / f"sums305_{split}.npz", allow_pickle=False)
    meta = json.loads((X.HERE / f"sums305_{split}.json").read_text())
    S, names, halves = z["sums"], z["names"].astype(str), z["halves"]
    genes = z["genes"].astype(str)
    tot = {}
    hf = {}
    for i, (c, h) in enumerate(zip(names, halves)):
        hf[(c, int(h))] = S[i]
        tot[c] = tot.get(c, 0) + S[i]
    return meta, genes, tot, hf, z["cells"], names, halves


def main() -> None:
    P = X.L.Logger(X.HERE / "ka305.log")
    ka, fails = [], []

    def row(label, ok, detail):
        ka.append(dict(label=label, ok=bool(ok), detail=detail))
        P(f"  {'PASS' if ok else 'FAIL'}  {label}: {detail}")
        if not ok:
            fails.append(label)

    D = {s: load(s) for s in X.SPLITS}
    mapping = json.loads((X.HERE / "mapping305.json").read_text())
    arms = list(csv.DictReader(open(E298 / "seed298.csv")))

    # ---------------------------------------------------------------- KA-1
    ident = sum(1 for a in arms if a["construct"] in D[a["split"]][0]["cells_per_half"].keys() or
                a["construct"] in D[a["split"]][2])
    nfail = sum(1 for a in arms if a["failing"] == "True")
    cellok = sum(1 for a in arms
                 if int(a["cells"]) == int(D[a["split"]][0]["cells_per_half"][f"{a['construct']}|0"])
                 + int(D[a["split"]][0]["cells_per_half"][f"{a['construct']}|1"]))
    row("KA-1 34 construct identities, 13/21 labels, 34/34 cell counts",
        ident == 34 and nfail == 13 and len(arms) - nfail == 21 and cellok == 34,
        f"identities {ident}/34, failing {nfail}, working {len(arms) - nfail}, cell counts {cellok}/34")

    # ---------------------------------------------------------------- KA-2
    gsets = [D[s][1] for s in X.SPLITS]
    same = all(np.array_equal(gsets[0], g) for g in gsets[1:])
    cands = sorted({r["candidate_gene"] for r in csv.DictReader(open(X.HERE / "candidates305.csv"))})
    present = [c for c in cands if c in set(gsets[0].tolist())]
    row("KA-2 gene order identical across splits; candidate genes present",
        same and len(present) == len(cands),
        f"order identical {same} ({len(gsets[0])} genes); candidates present {len(present)}/{len(cands)}")

    # ---------------------------------------------------------------- KA-3
    cm = json.loads((X.HERE / "candidates305_meta.json").read_text())
    ok3 = (cm["counts"]["failing"] == [4, 13] and cm["counts"]["working"] == [4, 21]
           and cm["counts"]["ntc"] == [4, 31])
    row("KA-3 301's primary contrast reproduced", ok3,
        f"failing {cm['counts']['failing']}, working {cm['counts']['working']}, ntc {cm['counts']['ntc']}")

    # ---------------------------------------------------------------- KA-4 synthetic known fold change
    meta, genes, tot, hf, cells, names, halves = D["train"]
    gpos = {g: i for i, g in enumerate(genes)}
    ntc_par = meta["ntc_parity"]
    ctrl = sum(tot[c] for c in meta["ntc_constructs"])
    c0, j = "CAST_P1", gpos["UTP18"]
    f = 0.5
    e_before = float(X.effect(tot[c0], ctrl)[j])
    Sp = tot[c0].astype(np.float64).copy()
    Sp[j] *= f
    e_after = float(X.effect(Sp, ctrl)[j])
    # independent closed form, pure-python math, no numpy and no lib285
    T = float(tot[c0].sum()); s = float(tot[c0][j]); Tc = float(ctrl.sum()); sc = float(ctrl[j])
    an_before = math.log1p(X.L.TS * s / T) - math.log1p(X.L.TS * sc / Tc)
    an_after = math.log1p(X.L.TS * (f * s) / (T - s + f * s)) - math.log1p(X.L.TS * sc / Tc)
    d1, d2 = abs(e_before - an_before), abs(e_after - an_after)
    row("KA-4 synthetic known fold change (f=0.5 on one column) recovered",
        d1 <= 1e-9 and d2 <= 1e-9,
        f"|recovered-analytic| before {d1:.2e}, after {d2:.2e}; effect moved "
        f"{e_before:+.4f} -> {e_after:+.4f} (analytic delta {an_after - an_before:+.6f}, "
        f"log(f) = {math.log(f):+.6f})")

    # ---------------------------------------------------------------- KA-6 same-NTC-halves null
    ka6 = {}
    for s in X.SPLITS:
        meta_s, genes_s, tot_s, hf_s, *_ = D[s]
        par = meta_s["ntc_parity"]
        cA = sum(tot_s[c] for c in meta_s["ntc_constructs"] if par[c] == 0)
        cB = sum(tot_s[c] for c in meta_s["ntc_constructs"] if par[c] == 1)
        e = X.effect(cA, cB)
        ka6[s] = dict(median=float(np.median(e)), mean=float(e.mean()), sd=float(e.std()),
                      p2_5=float(np.percentile(e, 2.5)), p97_5=float(np.percentile(e, 97.5)),
                      candidates={g: float(e[gpos[g]]) for g in cands})
    worst = max(abs(v["median"]) for v in ka6.values())
    row("KA-6 same-NTC-halves null centred on zero", worst <= 0.01,
        "; ".join(f"{s} median {v['median']:+.5f} sd {v['sd']:.4f} "
                  f"[2.5,97.5]=[{v['p2_5']:+.3f},{v['p97_5']:+.3f}]" for s, v in ka6.items()))

    # ---------------------------------------------------------------- KA-7 own-target effects vs 293
    r293 = {s: json.loads((E293 / f"run293_{s}.json").read_text()) for s in X.SPLITS}
    diffs, detail = [], []
    for a in arms:
        s, c, t = a["split"], a["construct"], a["target"]
        meta_s, genes_s, tot_s, *_ = D[s]
        gp = {g: i for i, g in enumerate(genes_s)}
        ctrl_s = sum(tot_s[x] for x in meta_s["ntc_constructs"])
        mine = float(X.effect(tot_s[c], ctrl_s)[gp[t]])
        theirs = float(r293[s]["kd_full"][c]["lfc"])
        diffs.append(abs(mine - theirs))
        detail.append(dict(construct=c, split=s, mine=mine, r293=theirs, absdiff=abs(mine - theirs)))
    row(f"KA-7 34 own-target effects reproduce 293's kd_full (tol {KA7_TOL:.0e} abs, full depth, "
        "all cells, pooled-NTC reference -- the exactly comparable configuration)",
        max(diffs) <= KA7_TOL, f"max |diff| {max(diffs):.3e} over 34 arms "
        f"(median |diff| {float(np.median(diffs)):.3e})")

    # ---------------------------------------------------------------- KA-5 (from the stream)
    ka5 = {s: D[s][0]["ka5"] for s in X.SPLITS}
    row("KA-5 chunked group sums == direct per-row slice", all(v["ok"] for v in ka5.values()),
        "; ".join(f"{s} {v['construct']} ({v['n_cells']} cells) max|diff| {v['max_abs_diff']}"
                  for s, v in ka5.items()))

    out = dict(ka=ka, failures=fails, ka6=ka6, ka7=detail, seed=SEED305)
    X.dump(out, X.HERE / "ka305.json")
    P(f"\n{len(ka) - len(fails)}/{len(ka)} known-answer rows pass"
      + ("" if not fails else f"; FAILURES: {fails}"))
    P.save()
    if fails:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
