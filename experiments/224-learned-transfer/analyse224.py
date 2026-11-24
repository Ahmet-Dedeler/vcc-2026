"""224 — read the ladder: rung x k, the slope in k, per-held-line breakdown, and the
paired target bootstrap that the pre-registered gate requires.

    uv run python analyse224.py --tag S7
"""
from __future__ import annotations

import argparse
import collections
import json

import numpy as np

from core224 import WORK

GATE = 0.029                     # AGENTS.md §8.4 corrected shape-only bar, conservative end
BOOT = 4000
SEED = 20260904


def load(tag):
    d = json.load(open(WORK / f"ladder224-{tag}.json"))
    pt = np.load(WORK / f"ptcos-{tag}.npz", allow_pickle=True)
    return d, pt


def rungs_of(d):
    out = []
    for r in d["rows"]:
        for m in r:
            if m not in ("cid", "held", "k", "donors") and m not in out:
                out.append(m)
    return out


def boot_rel(pt, rows, tag, rng, n=BOOT):
    """Paired bootstrap over TARGETS of the relative cosine gain, pooled over configs.
    Targets are the resampling unit and are resampled jointly across configs, so the
    same-target correlation between configs is preserved (a config-level bootstrap would
    understate the variance)."""
    a, b = [], []
    for r in rows:
        ka, kb = f"{r['cid']}|{tag}", f"{r['cid']}|R0"
        if ka not in pt or kb not in pt:
            return None
        a.append(pt[ka].astype(np.float64))
        b.append(pt[kb].astype(np.float64))
    A, B = np.stack(a), np.stack(b)              # (configs, targets)
    nt = A.shape[1]
    idx = rng.integers(0, nt, size=(n, nt))
    ga = A[:, idx].mean(axis=(0, 2))
    gb = B[:, idx].mean(axis=(0, 2))
    rel = 100.0 * (ga - gb) / gb
    obs = 100.0 * (A.mean() - B.mean()) / B.mean()
    return dict(rel=float(obs), lo=float(np.percentile(rel, 2.5)),
                hi=float(np.percentile(rel, 97.5)),
                p_pos=float((rel > 0).mean()), p_gate=float((rel > 100 * GATE).mean()))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="S7")
    args = ap.parse_args()
    rng = np.random.default_rng(SEED)
    d, pt = load(args.tag)
    rungs = rungs_of(d)
    print(f"surface {d['surface']}  lines {d['members']}\n"
          f"{d['n_targets']} targets x {d['n_genes']} genes; T_eval {d['t_eval']}; "
          f"permute_labels={d.get('permute_labels')}\n")

    by = collections.defaultdict(list)
    for r in d["rows"]:
        for m in rungs:
            if m in r:
                by[(m, r["k"])].append(r)

    print(f"{'rung':24s}{'k':>3}{'n':>4}{'cos':>9}{'rel%':>8}{'decoy':>9}{'exc.rel%':>10}"
          f"{'dpds':>9}{'pert.rel%':>10}")
    table = {}
    for m in rungs:
        for k in sorted({kk for (mm, kk) in by if mm == m}):
            R = by[(m, k)]
            v = np.array([[r[m]["line"], r["R0"]["line"], r[m]["pds"], r["R0"]["pds"],
                           r[m]["decoy"], r[m]["excess"], r["R0"]["excess"],
                           r[m]["pert"], r["R0"]["pert"]] for r in R])
            rel = 100 * (v[:, 0] - v[:, 1]).mean() / v[:, 1].mean()
            exr = 100 * (v[:, 5] - v[:, 6]).mean() / v[:, 6].mean()
            pr = 100 * (v[:, 7] - v[:, 8]).mean() / v[:, 8].mean()
            table[(m, k)] = dict(n=len(R), cos=v[:, 0].mean(), rel=rel,
                                 decoy=v[:, 4].mean(), exc_rel=exr,
                                 dpds=(v[:, 2] - v[:, 3]).mean(), pert_rel=pr)
            print(f"{m:24s}{k:3d}{len(R):4d}{v[:, 0].mean():9.5f}{rel:+8.2f}"
                  f"{v[:, 4].mean():9.5f}{exr:+10.2f}{(v[:, 2] - v[:, 3]).mean():+9.4f}"
                  f"{pr:+10.2f}")

    print(f"\nSLOPE IN k  (relative cosine gain, %, regressed on k) — the decision variable")
    print(f"{'rung':24s}{'slope/line':>12s}{'k range':>10s}{'rel@kmax':>10s}")
    slopes = {}
    for m in rungs:
        ks = sorted({kk for (mm, kk) in table if mm == m})
        if m == "R0" or len(ks) < 2:
            continue
        y = np.array([table[(m, k)]["rel"] for k in ks])
        s = np.polyfit(np.asarray(ks, float), y, 1)[0]
        slopes[m] = float(s)
        print(f"{m:24s}{s:+12.3f}{f'{ks[0]}-{ks[-1]}':>10s}{y[-1]:+10.2f}")

    print(f"\nPAIRED TARGET BOOTSTRAP at k = kmax  ({BOOT} draws; gate = +{100 * GATE}%)")
    print(f"{'rung':24s}{'rel%':>9}{'95% CI':>20}{'P(>0)':>8}{'P(>gate)':>10}  verdict")
    kmax = max(r["k"] for r in d["rows"])
    boots = {}
    for m in rungs:
        if m == "R0":
            continue
        rows = by.get((m, kmax), [])
        if not rows:
            continue
        b = boot_rel(pt, rows, m, rng)
        if b is None:
            continue
        boots[m] = b
        ok = b["lo"] > 0 and b["rel"] >= 100 * GATE
        print(f"{m:24s}{b['rel']:+9.2f}{f'[{b[chr(108)+chr(111)]:+.2f}, {b[chr(104)+chr(105)]:+.2f}]':>20}"
              f"{b['p_pos']:8.3f}{b['p_gate']:10.3f}  {'CLEARS' if ok else 'no'}")

    print(f"\nPER-HELD-LINE, at k = kmax  (relative cosine gain %)")
    lines = d["members"]
    hdr = f"{'rung':24s}" + "".join(f"{ln[:11]:>13s}" for ln in lines)
    print(hdr)
    perline = {}
    for m in rungs:
        if m == "R0":
            continue
        row = f"{m:24s}"
        perline[m] = {}
        for ln in lines:
            R = [r for r in by.get((m, kmax), []) if r["held"] == ln]
            if not R:
                row += f"{'-':>13s}"
                continue
            a = np.mean([r[m]["line"] for r in R])
            b = np.mean([r["R0"]["line"] for r in R])
            perline[m][ln] = float(100 * (a - b) / b)
            row += f"{100 * (a - b) / b:+13.2f}"
        print(row)

    (WORK / f"analyse224-{args.tag}.json").write_text(json.dumps(dict(
        table={f"{m}|{k}": v for (m, k), v in table.items()}, slopes=slopes,
        bootstrap=boots, per_line=perline, kmax=int(kmax), gate=GATE), indent=1))


if __name__ == "__main__":
    main()
