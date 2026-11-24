"""534 Part C -- SECONDARY (not preregistered): the same three steps when the toggled base keeps the PRODUCTION base's
LODO key (B1pk..B3pk). The primary designs re-fit the key on each design's own base (418's convention); an October knob
that switched a step off in the builder alone would still load the key artifact that def_lodo_key.py computes on the
production base (it hard-codes normrel cap -1.5 and PARTIAL_G 0.5), so this is what such a knob would deliver.
Same route, same rule applied descriptively; the decision table (results534C.json, RESULTS-C.md section 1) is unchanged.

KNOWN ANSWER: every primary step value, recomputed from the pk run's files, equals results534C.json's exactly.

    uv run python aggregate534cpk.py        -> work/partC/results534Cpk.json, work/partC/tables534Cpk.md
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
from aggregate534c import NB, SEED, summ, f5  # noqa: E402

HERE = P5.HERE


def main():
    t0 = time.time()
    P = {p: P5.Panel(p, sfx="pk") for p in D.PANELS}
    prim = json.loads((HERE / "results534C.json").read_text())["steps"]
    worst = 0.0
    for p in D.PANELS:                                      # KA-PRIMARY on the priced steps
        pt = P5.point(P[p], level=True, steps=D.STEPS)
        for sid, *_ in D.STEPS:
            for f in D.CANCER:
                worst = max(worst, abs(pt[sid]["per_fold"][f] - prim[sid]["panels"][p]["per_fold"][f]))
    assert worst == 0.0, f"KA-PRIMARY step values: {worst}"
    print(f"KA-PRIMARY: the primary step values reproduce from the pk run's files exactly ({len(D.PANELS) * 9} values)")
    steps = P5.STEPS_PK
    pt = {p: P5.point(P[p], level=True, steps=steps) for p in D.PANELS}
    boot_oct = P5.bootstrap([P[p] for p in D.OCT], NB, SEED, steps=steps)
    boot_abc = P5.bootstrap([P["abc"]], NB, SEED, steps=steps)
    out = {}
    for sid, arm, ref, what, board, bsrc in steps:
        folds = [pt[p][sid]["per_fold"][f] for p in D.OCT for f in D.CANCER]
        d = float(np.mean(folds))
        ci = summ(boot_oct[sid])
        rule = "HOLDS" if (d > 0 and ci["ci95"][0] > 0) else ("FRAGILE" if d > 0 else "FAILS")
        out[sid] = dict(what=what + " -- the toggled base with the production base's LODO key", arm=arm, ref=ref,
                        pooled=dict(mean=d, n_pos=int(sum(x > 0 for x in folds)), n_folds=9, **ci), rule_reading=rule,
                        primary=dict(mean=prim[sid.replace("_pk", "")]["pooled"]["mean"], verdict=prim[sid.replace("_pk", "")]["verdict"]),
                        panels={p: dict(mean=pt[p][sid]["mean"], n_pos=pt[p][sid]["n_pos"], per_fold=pt[p][sid]["per_fold"]) for p in D.PANELS},
                        abc=dict(mean=pt["abc"][sid]["mean"], n_pos=pt["abc"][sid]["n_pos"], **summ(boot_abc[sid])),
                        board=dict(value=board, source=bsrc))
    res = dict(written=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), note=__doc__.split("\n\n")[0], n_boot=NB, seed=SEED,
               ka_primary_worst=worst, steps=out, seconds=time.time() - t0)
    (HERE / "results534Cpk.json").write_text(json.dumps(res, indent=1, default=float))
    L = ["| step (production key kept) | Δ pooled real/efp/efr | 95 % CI | folds positive (of 9) | the rule would read | primary (key re-fitted) | Δ on abc [95 % CI] | real | efp | efr |",
         "|---|---:|---|---|---|---|---|---:|---:|---:|"]
    for sid, s in out.items():
        po, a = s["pooled"], s["abc"]
        L.append(f"| {sid} | {f5(po['mean'])} | [{f5(po['ci95'][0])}, {f5(po['ci95'][1])}] | {po['n_pos']}/9 | {s['rule_reading']} | "
                 f"{f5(s['primary']['mean'])} {s['primary']['verdict']} | {f5(a['mean'])} [{f5(a['ci95'][0])}, {f5(a['ci95'][1])}] | "
                 + " | ".join(f"{f5(s['panels'][p]['mean'])} ({s['panels'][p]['n_pos']}/3)" for p in D.OCT) + " |")
    (HERE / "tables534Cpk.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
