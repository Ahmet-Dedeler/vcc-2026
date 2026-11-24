#!/usr/bin/env python3
"""711 reader: the preregistered readings and the October K rule (PREREGISTRATION.md) for ku38 (K 4,200).

Written by the coordinator [2d6f46] on 2026-10-16, before 711's build or read exists; the rule is copied from the
preregistration, which fixed it before 711's read and before the H1 size curve.
Base 446bj2z (be90smZUDcIV7j5MgQ03, K 2,600); ku30 (HAbelEdgYYziqwfjTgdI, K 3,400) is reported beside it.
Primary: Δreach raw vs 446bj2z. Also: Δ(fid + jac) in overall units (701's statistic), Δnmae raw, Δoverall.
Falsifiers (STOP before reading): |Δpds| or |Δmse| raw > 1e-6 vs 446bj2z (a bulk leak); |Δnmae| raw > 0.002. The two
build-side falsifiers (config diff; realised calls grow by < +400 over ku30, via 701's deliver701.py) are checked separately.
October K = the size with the largest board Δoverall among {2,600: 0, 3,400: ku30, 4,200: 711}, only if
  (i)  experiments/701-kup-callset/work/reach701.json says verdict_i_at_f087 = true (POST HOC H1 size curve). AMENDMENT 2: if
       ku38's board Δreach ≤ +0.0019 raw, the f 0.61 dial is validated out of sample and (i') replaces (i), per size, from
       reach701_f061.json (paired Δreach vs 2,600 > 0 on 3/3 seeds, mean > 2 se); AND
  (ii) AMENDMENT 1 (step-down, written before 711's read): start at the best size; if its Δoverall exceeds the next-smaller
       size's by ≥ +0.0005 adopt it, else step down and test that size against ITS next-smaller one, until one passes or 2,600;
otherwise K stays 2,600 (rc_callset_cpm off). If K 3,400 → DEF_DAY_RC_CALLSET_CPM=3000; if K 4,200 → 3800.
usage: uv run python experiments/711-kup-dose/read711.py <711 entry json>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BOARD = ROOT / "data/board-entries"
BASE, KU30 = "be90smZUDcIV7j5MgQ03", "HAbelEdgYYziqwfjTgdI"
REACH701 = ROOT / "experiments/701-kup-callset/work/reach701.json"
# Amendment 2: f 0.61 per-size paired verdicts, written by science BEFORE ku38's read, format
# {"regime": 0.61, "per_size": {"3400": {"dreach_per_seed": [..], "mean": x, "se": y, "pass": bool}, "4200": {...}}}
REACH061 = ROOT / "experiments/701-kup-callset/work/reach701_f061.json"
DIAL_THRESHOLD = 0.0019   # Amendment 2: ku38 board Δreach raw ≤ this validates the f 0.61 dial out of sample
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "jac": "de_wilcoxon_sig_jaccard"}
CPM = {2600: None, 3400: 3000, 4200: 3800}


def entry(p):
    d = json.loads(Path(p).read_text())
    return d.get("entry", d)


def deltas(arm, base):
    raw = {k: arm[v] - base[v] for k, v in RAW.items()}
    fj = ((arm["score_fid"] - base["score_fid"]) + (arm["score_jac"] - base["score_jac"])) / 6.0
    return raw, fj, arm["score_avg"] - base["score_avg"]


def main():
    arm = entry(sys.argv[1])
    base, ku30 = entry(BOARD / f"{BASE}.json"), entry(BOARD / f"{KU30}.json")
    raw, fj, dov = deltas(arm, base)
    raw30, fj30, dov30 = deltas(ku30, base)
    stop = []
    if abs(raw["pds"]) > 1e-6 or abs(raw["mse"]) > 1e-6:
        stop.append(f"bulk leak: Δpds {raw['pds']:+.2e}, Δmse {raw['mse']:+.2e} raw (> 1e-6)")
    if abs(raw["nmae"]) > 0.002:
        stop.append(f"Δnmae {raw['nmae']:+.5f} raw (> ±0.002)")
    cond_i = None
    if REACH701.exists():
        cond_i = bool(json.loads(REACH701.read_text()).get("verdict_i_at_f087"))
    sizes = {2600: 0.0, 3400: dov30, 4200: dov}
    best = max(sizes, key=lambda k: sizes[k])
    order = sorted(sizes)
    # Amendment 2: ku38's unseen board reach picks the H1 dial. At or below the threshold, f 0.61 is validated out of sample and
    # condition (i') replaces (i), per size; above it, f 0.87's global condition (i) stands.
    dial061 = raw["reach"] <= DIAL_THRESHOLD + 1e-12   # tolerance: the board reports ~16 digits; a value AT the threshold counts
    f061 = json.loads(REACH061.read_text()).get("per_size", {}) if REACH061.exists() else None

    def cond(size):
        if dial061:
            return None if f061 is None else bool(f061.get(str(size), {}).get("pass"))
        return cond_i

    # Amendments 1 + 2, the step-down: walking down from the best size, adopt the first size that clears +0.0005 over its
    # next-smaller size AND passes the dial's condition at that size (2,600 needs neither).
    cand, steps, pending = best, [], False
    while cand != order[0]:
        prev = order[order.index(cand) - 1]
        step = sizes[cand] - sizes[prev]
        c = cond(cand)
        steps.append({"size": cand, "vs": prev, "step": step, "condition": c})
        if c is None:
            pending = True
            break
        if step >= 0.0005 and c:
            break
        cand = prev
    dial = "f 0.61 (validated by ku38: Δreach ≤ +0.0019)" if dial061 else "f 0.87 (preregistered; ku38 Δreach > +0.0019)"
    if stop:
        k_verdict = "STOP: " + "; ".join(stop) + " (diagnose before reading)"
    elif pending:
        need = "reach701_f061.json" if dial061 else "reach701.json"
        k_verdict = f"PENDING: dial {dial}; {need} does not exist yet; K stays 2,600 until it does"
    elif cand == 2600:
        k_verdict = f"K stays 2,600 (rc_callset_cpm off): dial {dial}; no larger size clears its step and its condition"
    else:
        last = steps[-1]
        k_verdict = (f"ADOPT K {cand:,}: DEF_DAY_RC_CALLSET_CPM={CPM[cand]} (dial {dial}; {cand:,} beats {last['vs']:,} by "
                     f"{last['step']:+.7f}{'; stepped down from ' + format(best, ',') if cand != best else ''}); "
                     "enters October only after R-ABC3 rehearses it")
    out = {"arm_entry": arm.get("entry_id"), "arm_overall": arm["score_avg"],
           "primary_delta_reach_raw_vs_446bj2z": raw["reach"], "ku30_delta_reach_raw": raw30["reach"],
           "delta_fid_plus_jac_overall": fj, "ku30_delta_fid_plus_jac_overall": fj30,
           "delta_raw_vs_446bj2z": raw, "delta_overall": {"2600": 0.0, "3400": dov30, "4200": dov},
           "condition_i_reach701_f087": cond_i, "dial": dial, "f061_present": f061 is not None, "best_size": best, "adopted_candidate": cand, "step_down": steps,
           "prediction": "Δreach ≥ ku30's +0.0054 raw if monotone; fid/jac near flat, jac may slip", "october_k": k_verdict}
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
