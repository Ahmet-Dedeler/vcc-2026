"""236 — tables.  Reads `work/stack236.json` (+ `work/six236.json` if present); writes
`work/analyse236.json` and prints every table that goes into RESULTS.md.

CONVENTIONS, fixed to reproduce 233's published table exactly (verified: all eight rows):
  rel cos = mean over folds of 100 * (arm - BASE) / BASE  on `cos_agg`, at a stated dose
  d pds   = mean over folds of (arm - BASE)

FOLD SETS.  `k562_essential` is banned on the `k562_gwps` fold (same cell line, same lab), so
subsets containing it are defined on two folds only.  Every number is therefore reported on
BOTH:
  CLEAN2  orion_hct116 + orion_hek293t -- defined for all 15 subsets, so the k-curve is
          apples-to-apples.  **This is the headline set.**
  DEFINED 233's convention: whatever folds the arm is defined on (3 where possible).
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
OCT = ("nadig_hepg2", "nadig_jurkat", "rpe1", "k562_essential")
CLEAN2 = ("orion_hct116", "orion_hek293t")
ALL3 = ("k562_gwps", "orion_hct116", "orion_hek293t")
DOSES = (0.05, 0.10, 0.20, 0.34, 0.45, 0.55)


def load():
    rows = json.loads((WORK / "stack236.json").read_text())
    by = {(r["panel"], r["fold"], r["arm"]): r for r in rows if r["arm"] != "__overlap__"}
    ov = {(r["panel"], r["fold"]): r for r in rows if r["arm"] == "__overlap__"}
    return rows, by, ov


def rel(by, panel, arm, folds, key="cos_agg"):
    """mean over `folds` of the relative (cos) or absolute (pds) gain vs BASE."""
    vals = []
    for f in folds:
        b, a = by.get((panel, f, "BASE")), by.get((panel, f, arm))
        if b is None or a is None:
            return None, 0
        vals.append(100 * (a[key] - b[key]) / b[key] if key == "cos_agg"
                    else a[key] - b[key])
    return (float(np.mean(vals)), len(vals)) if vals else (None, 0)


def main() -> None:
    rows, by, ov = load()
    out: dict = {}
    panels = sorted({r["panel"] for r in rows})

    # ---------------------------------------------------------------- Phase 0
    # recomputed here from the two raw JSONs so it does not depend on the run having finished
    o233 = {(r["panel"], r["fold"], r["arm"]): r for r in json.loads(
        (HERE.parents[1] / "experiments/233-october-retention/work/oct233.json").read_text())}
    ka = []
    for k, a in by.items():
        b = o233.get(k)
        if b is None:
            continue
        dc, dp = abs(b["cos_agg"] - a["cos_agg"]), abs(b["pds"] - a["pds"])
        ka.append(dict(panel=k[0], fold=k[1], arm=k[2], d_cos=dc, d_pds=dp,
                       ok=bool(dc < 1e-12 and dp < 1e-12)))
    out["phase0"] = dict(n_rows=len(ka), n_ok=sum(1 for r in ka if r["ok"]),
                         max_dcos=max([r["d_cos"] for r in ka], default=0.0),
                         max_dpds=max([r["d_pds"] for r in ka], default=0.0))
    print(f"PHASE 0 — {out['phase0']['n_ok']}/{out['phase0']['n_rows']} rows reproduce 233; "
          f"max |dcos| {out['phase0']['max_dcos']:.1e}, max |dpds| {out['phase0']['max_dpds']:.1e}")

    print("\nCOVERAGE (fold target counts, and each source's carrying rows)")
    cov = {}
    for p in panels:
        for f in ALL3:
            k = (p, f)
            if k not in ov:
                continue
            cov[f"{p}/{f}"] = dict(n_targets=ov[k]["n_targets"], carrying=ov[k]["carrying"])
            print(f"  {p:5s} {f:14s} n_targets {ov[k]['n_targets']:3d}   "
                  + "  ".join(f"{n}:{v}" for n, v in ov[k]["carrying"].items()))
    out["coverage"] = cov

    print("\nCARRYING-ROW OVERLAP (|A n B|; the mechanism behind the saturation)")
    for p in panels:
        for f in ALL3:
            if (p, f) not in ov:
                continue
            o = ov[(p, f)]["overlap"]
            names = list(ov[(p, f)]["carrying"])
            short = {"nadig_hepg2": "hepg2", "nadig_jurkat": "jurkat",
                     "rpe1": "rpe1", "k562_essential": "k562ess"}
            print(f"  {p}/{f}: " + "; ".join(
                f"{short.get(a, a)}^{short.get(b, b)}={o[f'{a}|{b}']}"
                for i, a in enumerate(names) for b in names[i + 1:]))
            full_arm = by.get((p, f, "+".join(n for n in OCT if n in names) + "@0.34"))
            print(f"      union of all {len(names)}: "
                  f"{full_arm['carrying_rows'] if full_arm else '?'} rows "
                  f"(singles {', '.join(str(v) for v in ov[(p, f)]['carrying'].values())}"
                  f"; fold has {ov[(p, f)]['n_targets']} targets)")
    out["overlap"] = {f"{p}/{f}": ov[(p, f)]["overlap"] for p in panels for f in ALL3
                      if (p, f) in ov}

    # ---------------------------------------------------------------- singles
    print("\nSINGLE SOURCES @ g=0.34 (233's table, reproduced; DEFINED fold set)")
    out["singles"] = {}
    for p in panels:
        for s in OCT:
            r, n = rel(by, p, f"{s}@0.34", ALL3)
            if r is None:
                r, n = rel(by, p, f"{s}@0.34", CLEAN2)
            dp, _ = rel(by, p, f"{s}@0.34", ALL3 if n == 3 else CLEAN2, "pds")
            sh, _ = rel(by, p, f"{s}_shuf@0.34", ALL3 if n == 3 else CLEAN2)
            c2, _ = rel(by, p, f"{s}@0.34", CLEAN2)
            if r is None:
                print(f"  {p:5s} {s:16s} (incomplete)")
                continue
            out["singles"][f"{p}/{s}"] = dict(rel=r, n_folds=n, d_pds=dp, shuf=sh, rel_clean2=c2)
            print(f"  {p:5s} {s:16s} n={n} rel {r:+7.2f}%  dpds {dp:+.5f}  "
                  f"shuf {sh:+6.2f}%  | CLEAN2 {c2:+7.2f}%")

    # ---------------------------------------------------------------- all subsets
    subsets = [s for k in range(1, 5) for s in itertools.combinations(OCT, k)]
    tab: dict = {}
    for p in panels:
        for names in subsets:
            t = "+".join(names)
            for g in DOSES:
                for foldset, lab in ((CLEAN2, "clean2"), (ALL3, "all3")):
                    r, n = rel(by, p, f"{t}@{g:g}", foldset)
                    dp, _ = rel(by, p, f"{t}@{g:g}", foldset, "pds")
                    sh, _ = rel(by, p, f"{t}_shuf@{g:g}", foldset)
                    if r is not None:
                        tab[(p, t, g, lab)] = dict(rel=r, d_pds=dp, shuf=sh, n=n, k=len(names))
    out["subsets"] = {f"{p}|{t}|{g}|{lab}": v for (p, t, g, lab), v in tab.items()}

    # ---------------------------------------------------------------- the stack curve
    print("\n" + "=" * 92)
    print("PHASE 1 — THE STACK CURVE.  CLEAN2 fold set (hct116 + hek293t), all 15 subsets.")
    print("=" * 92)
    curve: dict = {}
    for p in panels:
        print(f"\n--- panel {p} ---")
        # best single, at 0.34 and at its own argmax
        singles = {n: tab[(p, n, 0.34, "clean2")]["rel"] for n in OCT
                   if (p, n, 0.34, "clean2") in tab}
        best_single_034 = max(singles.values())
        best_single_name = max(singles, key=singles.get)
        arg = {}
        for names in subsets:
            t = "+".join(names)
            cand = [(tab[(p, t, g, "clean2")]["rel"], g) for g in DOSES
                    if (p, t, g, "clean2") in tab]
            if cand:
                arg[t] = max(cand)
        best_single_arg = max(arg[n] for n in OCT if n in arg)
        print(f"  best single @0.34 : {best_single_name} {best_single_034:+.2f}%")
        print(f"  best single @argmax: {best_single_arg[0]:+.2f}% (g={best_single_arg[1]})")

        # GREEDY cumulative order at fixed g = 0.34
        chosen: list[str] = []
        greedy = []
        for k in range(1, 5):
            cands = [n for n in OCT if n not in chosen]
            scored = []
            for n in cands:
                t = "+".join([x for x in OCT if x in chosen + [n]])
                if (p, t, 0.34, "clean2") in tab:
                    scored.append((tab[(p, t, 0.34, "clean2")]["rel"], n, t))
            if not scored:
                break
            v, n, t = max(scored)
            chosen.append(n)
            greedy.append(dict(k=k, added=n, stack=t, rel=v,
                               d_pds=tab[(p, t, 0.34, "clean2")]["d_pds"],
                               shuf=tab[(p, t, 0.34, "clean2")]["shuf"],
                               sat=v / best_single_034,
                               rel_argmax=arg[t][0], g_argmax=arg[t][1],
                               sat_argmax=arg[t][0] / best_single_arg[0]))
        print("\n  GREEDY cumulative, fixed g=0.34:")
        print("   k  added                stack rel%   sat x   d_pds     shuf%  |  argmax rel%  g*   sat x")
        for r in greedy:
            print(f"   {r['k']}  {r['added']:16s} {r['rel']:+8.2f}%  {r['sat']:5.2f}x "
                  f"{r['d_pds']:+.5f}  {r['shuf']:+6.2f}%  |  {r['rel_argmax']:+8.2f}%  "
                  f"{r['g_argmax']:.2f}  {r['sat_argmax']:5.2f}x")

        # BEST subset of each size, and whether ORDER matters (all 24 orders)
        bestk = {}
        for k in range(1, 5):
            cand = [(tab[(p, "+".join(s), 0.34, "clean2")]["rel"], "+".join(s))
                    for s in subsets if len(s) == k and (p, "+".join(s), 0.34, "clean2") in tab]
            if cand:
                bestk[k] = max(cand)
        print("\n  BEST subset of each size @0.34:  " + "  ".join(
            f"k={k}: {v[0]:+.2f}% ({v[1]})" for k, v in bestk.items()))

        orders = []
        for perm in itertools.permutations(OCT):
            path = []
            ok = True
            for k in range(1, 5):
                t = "+".join([x for x in OCT if x in perm[:k]])
                if (p, t, 0.34, "clean2") not in tab:
                    ok = False
                    break
                path.append(tab[(p, t, 0.34, "clean2")]["rel"])
            if ok:
                orders.append((perm, path))
        if orders:
            arr = np.array([p_ for _, p_ in orders])
            print(f"\n  ORDER SWEEP — all {len(orders)} orders of the four:")
            for k in range(4):
                print(f"    after {k+1} source(s): min {arr[:, k].min():+.2f}%  "
                      f"median {np.median(arr[:, k]):+.2f}%  max {arr[:, k].max():+.2f}%  "
                      f"spread {arr[:, k].max()-arr[:, k].min():.2f} pts")
            out.setdefault("orders", {})[p] = dict(
                n=len(orders),
                per_k={k + 1: dict(min=float(arr[:, k].min()),
                                   median=float(np.median(arr[:, k])),
                                   max=float(arr[:, k].max())) for k in range(4)})

        # leave-one-out of the full stack
        full = "+".join(OCT)
        if (p, full, 0.34, "clean2") in tab:
            fv = tab[(p, full, 0.34, "clean2")]["rel"]
            print("\n  LEAVE-ONE-OUT of the 4-stack @0.34:")
            for n in OCT:
                t = "+".join([x for x in OCT if x != n])
                if (p, t, 0.34, "clean2") in tab:
                    print(f"    -{n:16s} {tab[(p, t, 0.34, 'clean2')]['rel']:+8.2f}%  "
                          f"(costs {fv - tab[(p, t, 0.34, 'clean2')]['rel']:+.2f} pts)")

        curve[p] = dict(best_single_034=best_single_034, best_single_name=best_single_name,
                        best_single_argmax=best_single_arg[0], g_single=best_single_arg[1],
                        greedy=greedy, best_of_size={str(k): v for k, v in bestk.items()})
    out["curve"] = curve

    # ---------------------------------------------------------------- dose monotonicity
    print("\nDOSE MONOTONICITY of the 4-stack (CLEAN2):")
    for p in panels:
        full = "+".join(OCT)
        line = []
        for g in DOSES:
            if (p, full, g, "clean2") in tab:
                line.append(f"{g:g}:{tab[(p, full, g, 'clean2')]['rel']:+.2f}%")
        print(f"  {p:5s} " + "  ".join(line))
        sh = [f"{g:g}:{tab[(p, full, g, 'clean2')]['shuf']:+.2f}%" for g in DOSES
              if (p, full, g, "clean2") in tab]
        print(f"        shuf " + "  ".join(sh))

    (WORK / "analyse236.json").write_text(json.dumps(out, indent=2, default=float))
    print(f"\n-> {WORK / 'analyse236.json'}")


if __name__ == "__main__":
    main()
