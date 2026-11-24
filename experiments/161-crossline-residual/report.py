"""161 — aggregate the fold-level rows into the tables the verdict is read off."""
from __future__ import annotations

import argparse
import json
import os
from collections import defaultdict

import numpy as np

from common import WORK


def load(path):
    with open(path) as fh:
        return json.load(fh)


def table(rows, key, cols, lines, title, fmt="{:+.5f}"):
    """worst-line / mean over lines, for each level of `key`."""
    by = defaultdict(dict)
    for r in rows:
        by[key(r)][r["line"]] = r
    out = []
    for k, d in by.items():
        if len(d) < len(lines):
            continue
        v = {c: np.array([d[l][c] for l in lines]) for c in cols}
        out.append((float(v[cols[0]].min()), k, v, d))
    out.sort(reverse=True)
    print(f"\n{title}")
    hdr = f"  {'arm':<28} {'worst':>9} {'mean':>9}  " + "  ".join(f"{l:>9}" for l in lines)
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for worst, k, v, d in out:
        a = v[cols[0]]
        print(f"  {str(k):<28} {worst:+9.5f} {a.mean():+9.5f}  "
              + "  ".join(f"{x:+9.5f}" for x in a))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    a = ap.parse_args()
    for path in a.files:
        d = load(path)
        rows = d["rows"]
        lines = sorted({r["line"] for r in rows})
        print(f"\n{'=' * 118}\n{os.path.basename(path)}  —  surface {d['surface']}, "
              f"lines {', '.join(lines)}")

        base = {r["line"]: r for r in rows if r["regime"] == "baseline"}
        print(f"\nCONSENSUS BASELINE (the bar), delivered:")
        print(f"  {'line':<12} {'pds':>9} {'cos':>9} {'row-row':>9} {'sse*':>9}")
        for l in lines:
            r = base[l]
            print(f"  {l:<12} {r['pds']:9.5f} {r['cos']:9.5f} {r['rr']:+9.5f} {r['sse']:9.6f}")
        print(f"  {'MEAN':<12} "
              f"{np.mean([base[l]['pds'] for l in lines]):9.5f} "
              f"{np.mean([base[l]['cos'] for l in lines]):9.5f}")

        kn = [r for r in rows if r["regime"] == "known"]
        if kn:
            print(f"\nKNOWN-ANSWER ROWS (mean over lines):")
            for nm in dict.fromkeys(r["fset"] for r in kn):
                v = [r for r in kn if r["fset"] == nm]
                print(f"  {nm:<44} pds {np.mean([x['pds'] for x in v]):.5f}  "
                      f"cos {np.mean([x['cos'] for x in v]):.5f}")

        for reg in dict.fromkeys(r["regime"] for r in rows):
            if reg in ("baseline", "known"):
                continue
            sub = [r for r in rows if r["regime"] == reg and r.get("arm") == "orth"
                   and r.get("tgroup") is None]
            if sub:
                table(sub, lambda r: f"{r['fset']} w={r['w']}", ["d_pds"], lines,
                      f"{reg}: Δpds vs the consensus, by feature set and dose")
                table(sub, lambda r: f"{r['fset']} w={r['w']}", ["d_cos"], lines,
                      f"{reg}: Δaggregate-cosine vs the consensus")
            # held-target rows are averaged over the target groups
            grp = [r for r in rows if r["regime"] == reg and r.get("arm") == "orth"
                   and r.get("tgroup") is not None]
            if grp:
                agg = defaultdict(lambda: defaultdict(list))
                for r in grp:
                    agg[(r["fset"], r["w"])][r["line"]].append(r["d_pds"])
                out = []
                for k, dd in agg.items():
                    if len(dd) < len(lines):
                        continue
                    v = np.array([np.mean(dd[l]) for l in lines])
                    out.append((v.min(), k, v))
                out.sort(reverse=True)
                print(f"\n{reg}: Δpds on HELD TARGETS (mean over target groups)")
                print(f"  {'arm':<28} {'worst':>9} {'mean':>9}  "
                      + "  ".join(f"{l:>9}" for l in lines))
                for worst, k, v in out:
                    print(f"  {f'{k[0]} w={k[1]}':<28} {worst:+9.5f} {v.mean():+9.5f}  "
                          + "  ".join(f"{x:+9.5f}" for x in v))

        # nulls and bootstraps for the rows that carry them
        nb = [r for r in rows if "null_d_pds" in r and r.get("tgroup") is None]
        if nb:
            print(f"\nPERMUTATION NULL and TARGET-CLUSTER BOOTSTRAP "
                  f"(rows where they were computed)")
            print(f"  {'arm':<22} {'line':<10} {'d_pds':>9} {'boot 95%':>22} {'NULL d_pds':>11}")
            for r in sorted(nb, key=lambda r: (r["fset"], r["w"], r["line"])):
                lbl = f"{r['fset']} w={r['w']}"
                print(f"  {lbl:<22} {r['line']:<10} {r['d_pds']:+9.5f} "
                      f"[{r['boot_lo']:+8.5f},{r['boot_hi']:+8.5f}] {r['null_d_pds']:+11.5f}")

        # model diagnostics
        diag = {}
        for r in rows:
            if r.get("arm") == "orth" and r.get("tgroup") is None:
                diag.setdefault((r["regime"], r["fset"]), {})[r["line"]] = r
        print(f"\nMODEL DIAGNOSTICS (per feature set, mean over lines)")
        print(f"  {'regime/fset':<22} {'R^2 on y':>10} {'cos(yhat,C)':>12} "
              f"{'|orth|/|C|':>11} {'cos(orth, TRUE orth)':>21}")
        for (reg, fs), d in diag.items():
            v = list(d.values())
            print(f"  {f'{reg}/{fs}':<22} {np.mean([x['eval_r2'] for x in v]):10.4f} "
                  f"{np.mean([x['cos_yhat_C'] for x in v]):12.3f} "
                  f"{np.mean([x['relnorm'] for x in v]):11.3f} "
                  f"{np.mean([x['cos_orth_true'] for x in v]):21.4f}")


if __name__ == "__main__":
    main()
