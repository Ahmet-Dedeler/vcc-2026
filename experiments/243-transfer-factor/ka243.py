"""243 Phase 0 — the three mandated known-answer rows.  ABORT on any failure.

  KA-1  216 §4.4's PLOLO +17.98% (`cos_agg`, f 0.10 -> 0.34) and its k_cos 1.384, from
        experiments/216-cd4-hard/work/{sweep216,board216}.json, on 216's own convention
        (four folds, a23 tag, ratio of fold-means, board floor on the TRANSCRIBED raws).
  KA-2  222's +10.3065% (`cos_agg`) and +7.6803% (macro), three cancer folds, mean of
        per-fold ratios, from experiments/222-h1-union-build/work/repro222.json.
  KA-3  197's inversion: the two constants from (Kc, sumD, K); rc_hybrid's MEASURED cosine
        0.079248 from a15's own terms; the FLOOR on cd4f034 (0.83635 -> 0.09027) and
        h1union (0.81726 -> 0.09947) with 197's constants.

    uv run python ka243.py
"""
from __future__ import annotations

import json

import numpy as np

from lib243 import (CANCER, COEFF, CTX, FLOOR0, FOLDS, INV197, K_CAP, KC, P197, ROOT, SUMD,
                    WORK, board_row, floor_cos, rel_gain, repro222, sweep216)

rows = []


def check(name, got, want, tol, note=""):
    ok = abs(got - want) <= tol
    rows.append(dict(name=name, got=got, want=want, tol=tol, ok=bool(ok), note=note))
    print(f"{'PASS' if ok else 'FAIL'}  {name:58s} got {got:+.9f}  want {want:+.9f}  "
          f"|d| {abs(got-want):.2e}  {note}")
    return ok


def main() -> None:
    # ---------------------------------------------------------------- KA-1  (216 §4.4)
    s216 = sweep216()
    b216 = json.loads((ROOT / "experiments/216-cd4-hard/work/board216.json").read_text())
    before = {f: s216[(f, "f0.1")] for f in FOLDS}
    after = {f: s216[(f, "f0.34")] for f in FOLDS}
    rel = rel_gain(before, after, "cos_agg", FOLDS, "pooled")
    check("KA-1a PLOLO cos_agg f0.10->0.34, 4 folds, pooled", rel, b216["rel_plolo"], 1e-12,
          "(216 §4.4: +17.98%)")
    check("KA-1a' ...quoted to 2 dp of percent", round(100 * rel, 2), 17.98, 1e-9)
    m01 = float(np.mean([before[f]["cos_agg"] for f in FOLDS]))
    m34 = float(np.mean([after[f]["cos_agg"] for f in FOLDS]))
    check("KA-1b 4-fold mean cos_agg at f=0.10", m01, 0.020086, 5e-7)
    check("KA-1c 4-fold mean cos_agg at f=0.34", m34, 0.023698, 5e-7)
    # 216 used the TRANSCRIBED raws 0.87367 / 0.83635 (board216.py BOARD dict) and the
    # ROUNDED constants FLOOR_A, FLOOR_B = 1.023611, 2.074495 (board216.py line 47); the
    # full-precision constants differ from those at the 6e-8 level, so the row is
    # reproduced on 216's own constants and the full-precision variant is printed after.
    fl216 = lambda raw: (1.023611 - raw) / 2.074495  # noqa: E731
    c_ship, c_inc = fl216(0.87367), fl216(0.83635)
    print(f"      (full-precision constants would give {floor_cos(0.87367):.9f} / "
          f"{floor_cos(0.83635):.9f}, k_cos {(floor_cos(0.83635)/floor_cos(0.87367)-1)/rel:.6f})")
    check("KA-1d board floor at raw 0.87367", c_ship, b216["c_ship"], 1e-12)
    check("KA-1e board floor at raw 0.83635", c_inc, b216["c_inc"], 1e-12)
    relb = c_inc / c_ship - 1.0
    check("KA-1f board relative cosine (transcribed raws)", relb, b216["rel_board"], 1e-12,
          "(216 §4.4: +24.89%)")
    check("KA-1g k_cos = board/PLOLO", relb / rel, b216["k_cos"], 1e-12, "(216 §4.4: 1.384)")

    # ---------------------------------------------------------------- KA-2  (222)
    r222 = repro222()
    arms = {(r["fold"], r["family"], r["g"]): r for r in r222["arms"]}
    sel = {f: arms[(f, "h1PAIRCARRY", 0.28)] for f in FOLDS}
    base = {f: dict(cos_agg=sel[f]["base_cos_agg"], a23_cos=sel[f]["base_a23_cos"],
                    pds=sel[f]["base_pds"]) for f in FOLDS}
    agg3 = rel_gain(base, sel, "cos_agg", CANCER, "perfold")
    mac3 = rel_gain(base, sel, "a23_cos", CANCER, "perfold")
    summ = r222["summary"]["h1PAIRCARRY@0.28"]
    check("KA-2a 222 rel cos_agg, 3 cancer folds, per-fold mean", agg3,
          summ["rel_agg_nokolf"], 1e-12, "(222: +10.3065%)")
    check("KA-2a' ...to 4 dp of percent", round(100 * agg3, 4), 10.3065, 1e-9)
    check("KA-2b 222 rel macro, 3 cancer folds, per-fold mean", mac3,
          summ["rel_macro_nokolf"], 1e-12, "(222: +7.6803%)")
    check("KA-2b' ...to 4 dp of percent", round(100 * mac3, 4), 7.6803, 1e-9)
    dp3 = float(np.mean([sel[f]["pds"] - base[f]["pds"] for f in CANCER]))
    check("KA-2c 222 d pds, 3 cancer folds", dp3, 0.004886737506532087, 1e-12)

    # ---------------------------------------------------------------- KA-3  (197)
    kc_bar = float(np.mean([KC[c] for c in CTX]))
    check("KA-3a floor0 = 1 + mean Kc / sumD", 1.0 + kc_bar / SUMD, FLOOR0, 1e-12)
    check("KA-3b coeff = 2 sqrt(K / sumD)", 2.0 * np.sqrt(K_CAP / SUMD), COEFF, 1e-12)
    # rc_hybrid: a15's own per-context terms (196's streamed totals, 191's claim), only
    # `claim` moved onto the cap (rchybrid197.py), board raw from the live row.
    per = {c: dict(INV197["ours"]["a15"]["per"][c]) for c in CTX}
    c191 = json.loads((ROOT / "experiments/191-rc-hybrid-build/work/claim191.json").read_text())
    for c in CTX:
        per[c]["claim"] = per[c]["n"] * c191["rc_hybrid"][c]["claim_per_row"] / 18532.0
    rc = board_row("rc_hybrid")
    num = float(np.mean([per[c]["A"] + per[c]["Kc"] - min(per[c]["claim"], per[c]["B"])
                         for c in CTX]))
    X = (num - (rc["raw"]["mse"] - 1.0) * SUMD) / 2.0
    cos_rc = X / (np.sqrt(SUMD) * float(np.mean([np.sqrt(per[c]["A"]) for c in CTX])))
    rch = json.loads((P197 / "rchybrid197.json").read_text())
    check("KA-3c rc_hybrid board raw mse (live row)", rc["raw"]["mse"], rch["rc_raw"], 1e-12,
          "(mandate: 0.8738045)")
    check("KA-3d rc_hybrid MEASURED cosine, 197's inversion", cos_rc, rch["rc_cos"], 1e-9,
          "(mandate: 0.079248)")
    check("KA-3d' ...to 6 dp", round(cos_rc, 6), 0.079249, 1.5e-6, "(197 §1: 0.079249)")
    check("KA-3e FLOOR at raw 0.83635 (cd4f034, transcribed)", floor_cos(0.83635), 0.09027, 5e-6,
          "(mandate: 0.09027)")
    check("KA-3f FLOOR at raw 0.81726 (h1union, transcribed)", floor_cos(0.81726), 0.09947, 5e-6,
          "(mandate: 0.09947)")
    # the exact raws, for the record
    for arm in ("rc_callset", "cd4f034", "h1union"):
        b = board_row(arm)
        rows.append(dict(name=f"exact raw mse {arm}", got=b["raw"]["mse"], want=None,
                         floor=floor_cos(b["raw"]["mse"]), ok=True))
        print(f"      {arm:11s} exact raw mse {b['raw']['mse']:.10f} -> floor "
              f"{floor_cos(b['raw']['mse']):.7f}   raw pds {b['raw']['pds']:.10f}   "
              f"[{b['snapshot']}, {b['id']}]")

    ok = all(r["ok"] for r in rows)
    (WORK / "ka243.json").write_text(json.dumps(rows, indent=1))
    print(f"\n{'ALL PASS' if ok else 'FAILURES — STOP'}  ({sum(1 for r in rows if not r['ok'])} fail)")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
