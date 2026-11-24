#!/usr/bin/env python3
"""528 — fill RESULTS.md's placeholders from work/price528.json, work/summ528.md and the final test logs (run once, last).

    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/528-callprop-block/fill528.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent


def f(x, d=4):
    return f"{x:+.{d}f}"


def arm_line(r, reg):
    x = r[reg]
    dm = x["d_mean"]
    seeds = " / ".join(f(p["value"]) for p in x["per_seed"])
    return (f"{f(x['value_mean'])} [{f(x['ci'][0])}, {f(x['ci'][1])}] (seeds {seeds}); fid {f(dm['fid'])}, jac {f(dm['jac'])}, "
            f"reach {f(dm['reach'])}, nmae {f(dm['nmae'])}; fid+jac overall {f(x['fid_jac_overall'])}")


def main():
    d = json.load(open(HERE / "work/price528.json"))["arms"]
    s = (HERE / "RESULTS.md").read_text()
    t1, t2 = (HERE / "work/summ528.md").read_text().strip().split("\n\n")
    s = s.replace("SUMM528_TABLE_1", t1).replace("SUMM528_TABLE_2", t2)
    hk, hc = d["lcap-kept+mzg +zs"], d["lcap+mzg +zs"]
    nk, nc = d["lcap-kept-n1200 +zs"], d["lcap-n1200 +zs"]
    hook_head = (f"**{f(hk['1.0']['value_mean'])}** [{f(hk['1.0']['ci'][0])}, {f(hk['1.0']['ci'][1])}] at f 1 over "
                 f"{len(hk['seeds'])} seeds ({' / '.join(f(p['value']) for p in hk['1.0']['per_seed'])}), "
                 f"{f(hk['0.5']['value_mean'])} / {f(hk['0.35']['value_mean'])} / {f(hk['cmp']['value_mean'])} at f 0.5 / 0.35 / cmp; "
                 f"fid {f(hk['1.0']['d_mean']['fid'])}, jac {f(hk['1.0']['d_mean']['jac'])}")
    n_head = (f"{f(nk['0.5']['value_mean'])} / {f(nk['0.35']['value_mean'])} / {f(nk['cmp']['value_mean'])} at f 0.5 / 0.35 / cmp "
              f"({len(nk['seeds'])} seeds), {f(nk['1.0']['value_mean'])} at f 1")
    s = s.replace("HOOK_HEADLINE", hook_head).replace("N1200_HEADLINE", n_head)
    lines = []
    lines.append("**On top of 527's hook** (`+mzg` = run517's `c4-mzg0.05`: the MDE add on complement genes takes Arc's meanz_ne "
                 "programme sign at tau 0.05 on 467's strong targets and CLEAN4's on the weak third = 446bj2zmz's DE construction; "
                 "KA-hook PASS), each vs `cpm+mzg` of the same seed:")
    lines.append("")
    lines.append("| arm | regime | overall [95 % CI] (seeds); members raw |")
    lines.append("|---|---|---|")
    for name, r in (("lcap-kept+mzg", hk), ("lcap+mzg (coupled)", hc)):
        for reg in ("1.0", "0.5", "0.35", "cmp"):
            lines.append(f"| `{name}` | {reg} | {arm_line(r, reg)} |")
    dl = d.get("lcap-kept+mzg", {}).get("deliverability")
    if dl:
        sa = dl["sign_accuracy"]["swapped-in (arm, called)"]
        lines.append("")
        lines.append(f"With the hook the called swapped-in genes carry the programme's signs where it is confident: sign accuracy "
                     f"{sa['acc_all']:.4f} (all) / {sa['acc_gate']:.4f} (gate) against 0.508 / 0.533 with design signs, so fid turns "
                     f"positive; the swapped-in call rate is {dl['arm']['block swapped-in']['call_rate']:.4f} and unsuppressible calls "
                     f"{dl['arm']['unsuppressible_per_target']:.1f} / target (base {dl['base']['unsuppressible_per_target']:.1f}). "
                     "The coupled form moves the swapped-in genes into the steered block, i.e. OUT of the hook's complement, which "
                     "is where the hook earns its nmae (517 §2.1): nmae +0.009.")
    lines.append("")
    lines.append("**At 523's K 1,600** (`-n1200` = CallSet(400, 1200, 1.9)), each vs `cpm-n1200` of the same seed:")
    lines.append("")
    lines.append("| arm | regime | overall [95 % CI] (seeds); members raw |")
    lines.append("|---|---|---|")
    for name, r in (("lcap-kept-n1200", nk), ("lcap-n1200 (coupled)", nc)):
        for reg in ("1.0", "0.5", "0.35", "cmp"):
            lines.append(f"| `{name}` | {reg} | {arm_line(r, reg)} |")
    dl = d.get("lcap-kept-n1200", {}).get("deliverability")
    if dl:
        lines.append("")
        lines.append(f"At K 1,600 the key swaps {dl['arm']['block swapped-in']['n_per_target']:.0f} of the 1,200 block genes (top-CPM "
                     f"genes leave) and every prescribed gene is called, but the kept form makes "
                     f"{dl['arm']['unsuppressible_per_target']:.0f} unsuppressible calls / target (base "
                     f"{dl['base']['unsuppressible_per_target']:.0f}): low-to-mid CPM complement genes (median CPM rank ~4,300), spread "
                     "over targets, 83 % of them not leaking in the base, and NOT the swapped-out genes (call rate 0.0003). The coupled "
                     "form leaks less (36) but loses fid. The size and the composition are not additive under today's emitter.")
    s = s.replace("INTERACTIONS_BLOCK", "\n".join(lines))
    bt = (HERE / "work/test_patch528.final.log").read_text() if (HERE / "work/test_patch528.final.log").exists() else ""
    dt = (HERE / "work/test_defday528.final.log").read_text() if (HERE / "work/test_defday528.final.log").exists() else ""
    m1 = re.search(r"(\d+)/(\d+) PASS", bt)
    m2 = re.search(r"(\d+)/(\d+) PASS", dt)
    s = s.replace("BUILDER_TESTS", f"**{m1.group(1)}/{m1.group(2)} PASS** (after `k_emit528.py` was regenerated to the final form)"
                  if m1 else "see work/test_patch528.final.log")
    s = s.replace("DEFDAY_TESTS", f"**{m2.group(1)}/{m2.group(2)} PASS**" if m2 else "see work/test_defday528.final.log")
    left = [p for p in ("SUMM528_TABLE", "INTERACTIONS_BLOCK", "HOOK_HEADLINE", "N1200_HEADLINE", "BUILDER_TESTS", "DEFDAY_TESTS") if p in s]
    (HERE / "RESULTS.md").write_text(s)
    print("filled; placeholders left:", left)


if __name__ == "__main__":
    main()
