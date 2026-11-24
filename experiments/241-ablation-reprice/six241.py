"""241 — the four DE members BEFORE (200's renderer) and AFTER (the shipped renderer).

Every pair is matched-seed and every base is the base of ITS OWN renderer.  The `art_*` base
for the shipped renderer is 200's own `de200.art_nomf.*` — this experiment re-emitted seed 1
of it and verified the table is BIT-IDENTICAL, so 200's five seeds are the same object.

Also reports 188's `nmae` identity split (`M`/`U`) and `e = dM/dU` per row: §8.1d's test of
whether a component's contribution is SIGNS (`e` well above 0.5) or MAGNITUDES (`e` ~ 0.07).
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
E200 = os.path.abspath(os.path.join(HERE, "..", "200-production-ablation"))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, HERE)
from lib241 import (Surface, efficiency, nmae_parts, paired_de,           # noqa: E402
                    row_for, seed_names)

os.makedirs(WORK, exist_ok=True)

# (row label, BEFORE arm, BEFORE base, AFTER arm, AFTER base, n seeds, note)
ROWS = [
    ("union whitening operator (γ 0.65)", "reb_ship", "reb_noWhiten",
     "reb_ship_nomf", "reb_noWhiten_nomf", 3, ""),
    ("equal donor weights (1,1,1)", "reb_ship", "reb_w41",
     "reb_ship_nomf", "reb_w41_nomf", 3, ""),
    ("--overdispersion 4.0", "art_mf", "art_phi0_mf",
     "art_nomf", "art_phi0_nomf", 5, ""),
    ("tiered spread policy (8000/1.0)", "art_mf", "art_flat_mf",
     "art_nomf", "art_flat_nomf", 5, ""),
    ("PerturbMatch K562 replacement (0.5)", "reb_ship", "reb_noPM",
     "reb_ship_nomf", "reb_noPM_nomf", 3, ""),
    # the same row on the SHIPPED EMITTER: BEFORE = shipped renderer + `quasi`,
    # AFTER = shipped renderer + `rc_callset` (202's CallSet(400, 2200, 1.90))
    ("PerturbMatch on rc_callset", "reb_ship_nomf", "reb_noPM_nomf",
     "rc_ship_nomf", "rc_noPM_nomf", 3, "emitter quasi -> rc_callset"),
    ("CPM gate (--min-cpm 10)", "reb_ship_gate10", "reb_ship_gate0",
     "reb_ship_nomf_gate10", "reb_ship_nomf_gate0", 3, ""),
    ("post-whitening CD4 residual", "reb_ship", "reb_noCD4",
     "reb_ship_nomf", "reb_noCD4_nomf", 1, "coverage 4/85 -> solved fraction 0"),
    ("a23 mass fix (minnorm)", "art_mf", "art_nomf",
     "art_mf", "art_nomf", 5, "OFF in the shipped stack: read as an ADD-IN"),
]
PREFIX = {True: os.path.join(E200, "de200."), False: os.path.join(HERE, "de241.")}


def path_of(arm: str) -> str:
    p241 = os.path.join(HERE, f"de241.{arm}.npz")
    p200 = os.path.join(E200, f"de200.{arm}.npz")
    if os.path.exists(p241):
        return p241
    if os.path.exists(p200):
        return p200
    return None


def main() -> None:
    t0 = time.time()
    need = set()
    for _, ba, bb, aa, ab, n, _ in ROWS:
        for stem in (ba, bb, aa, ab):
            need |= set(seed_names(stem, n))
    cache = os.path.join(WORK, "six241_rows.json")
    rows = json.loads(open(cache).read()) if os.path.exists(cache) else {}
    for arm in sorted(need):
        if arm in rows:
            continue
        p = path_of(arm)
        if p is None:
            print(f"  (missing {arm} — row will be reported as not emitted)")
            continue
        r = row_for(p)
        r.update(nmae_parts(Surface(p)))
        r["src"] = "200" if "de200." in p else "241"
        rows[arm] = r
        json.dump(rows, open(cache, "w"), indent=1, default=float)
        print(f"  scored {arm:<28} src {r['src']}  [{time.time()-t0:.0f}s]", flush=True)

    print("\n" + "=" * 130)
    print("FOUR DE MEMBERS PER COMPONENT — BEFORE (VCC_MASS_FIX=minnorm) vs AFTER (shipped, unset)")
    print("positive = what removing the component would cost us, in overall board units")
    print("=" * 130)
    out = {}
    hdr = (f"  {'component':<34}{'when':>7}{'reach':>10}{'fid':>10}{'jac':>10}{'nmae':>10}"
           f"{'SUM4':>11}{'jac_np':>10}{'draw reach':>11}{'t':>7}{'e=dM/dU':>9}")
    print(hdr)
    for label, ba, bb, aa, ab, n, note in ROWS:
        rec = {}
        for when, arm, base in (("BEFORE", ba, bb), ("AFTER", aa, ab)):
            an, bn = seed_names(arm, n), seed_names(base, n)
            if any(x not in rows for x in an + bn):
                print(f"  {label:<34}{when:>7}   (not emitted)")
                continue
            d = paired_de(rows, an, bn)
            eff = efficiency(rows[an[0]], rows[bn[0]])
            d["e"] = eff["e"]
            d["dM"] = eff["dM"]
            d["dU"] = eff["dU"]
            d["d_nmae_identity"] = eff["d_nmae"]
            s4 = sum(d[m] for m in ("reach", "fid", "jac", "nmae"))
            d["sum4"] = s4
            rec[when] = d
            print(f"  {label:<34}{when:>7}{d['reach']:>+10.6f}{d['fid']:>+10.6f}"
                  f"{d['jac']:>+10.6f}{d['nmae']:>+10.6f}{s4:>+11.6f}{d['jac_np']:>+10.6f}"
                  f"{d['d_reach_raw']:>+11.5f}{d['reach_t']:>+7.2f}{d['e']:>9.3f}")
        if "BEFORE" in rec and "AFTER" in rec:
            print(f"  {'':<34}{'Δ':>7}"
                  + "".join(f"{rec['AFTER'][m] - rec['BEFORE'][m]:>+10.6f}"
                            for m in ("reach", "fid", "jac", "nmae"))
                  + f"{rec['AFTER']['sum4'] - rec['BEFORE']['sum4']:>+11.6f}"
                  + (f"   {note}" if note else ""))
        out[label] = rec

    json.dump(dict(rows=rows, table=out), open(os.path.join(WORK, "six241.json"), "w"),
              indent=1, default=float)
    print(f"\nwrote work/six241.json [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
