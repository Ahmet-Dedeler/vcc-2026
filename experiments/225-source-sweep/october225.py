"""225 — INBOX-OCTOBER: record BREADTH, not just today's overlap.

The 300-target panel is replaced on 2026-10-22 and the replacement is **disjoint** from the
current one, so a source selected because it happens to overlap today's 300 is worth **zero** in
the only submission that counts.  For every census member this records:

  1. total distinct perturbed targets in the source (the breadth number);
  2. overlap with the UNION of all donor coverage — the pool a new panel is actually drawn from;
  3. projected coverage under **DEF-RUNBOOK §5's two synthetic draws**, reusing 199's own bundles
     (`def-bundle-r` = realistic, drawn from K562 ∩ Orion ∩ CD4; `def-bundle` = pessimistic, drawn
     from any donor's coverage) rather than inventing a sampling model, exactly as the INBOX asks.

The ranking that matters in October is column (3), not today's overlap.

    uv run python october225.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from harness225 import ROOT, WORK, sources218

E199 = ROOT / "experiments/199-def-runbook/work"


def panel_of(p: Path) -> set[str]:
    return {line.split(",")[0].strip()
            for line in p.read_text().splitlines()[1:] if line.strip()}


TODAY = panel_of(ROOT / "data/challenge/pert_counts.csv")
REALISTIC = panel_of(E199 / "def-bundle-r/pert_counts.csv")     # K562 ∩ Orion ∩ CD4
PESSIMISTIC = panel_of(E199 / "def-bundle/pert_counts.csv")     # any donor
assert len(TODAY) == len(REALISTIC) == len(PESSIMISTIC) == 300
# DEF-RUNBOOK §5's stated property, re-asserted here rather than trusted
assert not (TODAY & REALISTIC), "the realistic draw is supposed to be disjoint from A/B/C"
assert not (TODAY & PESSIMISTIC), "the pessimistic draw is supposed to be disjoint from A/B/C"

# every source 225 priced, plus the production donors and 223's tables, for the reference rows
MEMBERS = ("dld1_lowmoi", "igvf_cardiomyocyte", "jiang_5line_noK562", "jiang_6line",
           "myllia_a549", "sunshine_calu3", "tian_ipsc_neuron", "xucao_hek293",
           "tian_ipsc_neuron_crispra", "norman_k562", "frangieh_melanoma",
           "gasperini_k562", "wessels_thp1", "tian2019_ipsc",
           "cmap_sh_cgs", "cmap_xpr", "kolf2", "telohaec", "altair_a549",
           "k562_gwps", "orion_hct116", "orion_hek293t", "h1_2025", "gse295214_h1",
           "feng_ipsc", "nadig_hepg2", "nadig_jurkat", "rpe1", "k562_essential")

DONOR_UNION_SOURCES = ("k562_gwps", "orion_hct116", "orion_hek293t")


def targets_of(name: str) -> set[str] | None:
    meta = sources218.SOURCES.get(name)
    if meta is None:
        return None
    p = ROOT / meta["path"] / "perts.txt"
    if not p.exists():
        return None
    return {x.strip() for x in p.read_text().split("\n") if x.strip()}


def main() -> None:
    donor_union: set[str] = set()
    for n in DONOR_UNION_SOURCES:
        t = targets_of(n)
        if t:
            donor_union |= t
    # CD4's banks are read through the production path; its target list is on disk
    cd4p = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"
    for cand in (cd4p / "perts.txt", cd4p / "targets.txt"):
        if cand.exists():
            donor_union |= {x.strip() for x in cand.read_text().split("\n") if x.strip()}
            break
    print(f"donor-coverage union: {len(donor_union):,} distinct targets "
          f"(the pool a new panel is drawn from)\n")

    rows = []
    for name in MEMBERS:
        t = targets_of(name)
        if t is None:
            print(f"{name:26s} — no table on disk, skipped")
            continue
        row = dict(
            source=name, breadth=len(t),
            genome_wide=bool(len(t) >= 5000),
            today=len(t & TODAY),
            donor_union_overlap=len(t & donor_union),
            donor_union_frac=len(t & donor_union) / max(len(donor_union), 1),
            realistic=len(t & REALISTIC),
            pessimistic=len(t & PESSIMISTIC),
        )
        # the decision-relevant disagreement: does today's ranking survive October?
        row["oct_vs_today"] = row["pessimistic"] - row["today"]
        rows.append(row)

    rows.sort(key=lambda r: -r["pessimistic"])
    print(f"{'source':26s}{'breadth':>9s}{'today':>7s}{'∩donorU':>9s}"
          f"{'realistic':>11s}{'PESSIM':>8s}{'Δ vs today':>12s}")
    for r in rows:
        print(f"{r['source']:26s}{r['breadth']:>9,d}{r['today']:>7d}"
              f"{r['donor_union_overlap']:>9,d}{r['realistic']:>11d}"
              f"{r['pessimistic']:>8d}{r['oct_vs_today']:>+12d}")

    # the specific question the INBOX asks: does the ranking change?
    by_today = [r["source"] for r in sorted(rows, key=lambda r: -r["today"])]
    by_pess = [r["source"] for r in sorted(rows, key=lambda r: -r["pessimistic"])]
    moved = [(s, by_today.index(s), by_pess.index(s)) for s in by_today
             if abs(by_today.index(s) - by_pess.index(s)) >= 3]
    print("\nrank changes of 3+ places between today's panel and the pessimistic draw:")
    for s, a, b in moved:
        print(f"  {s:26s} #{a + 1} today -> #{b + 1} in October ({b - a:+d})")

    out = dict(donor_union_size=len(donor_union), rows=rows,
               rank_by_today=by_today, rank_by_pessimistic=by_pess,
               moved_3_or_more=moved,
               draws=dict(realistic="199 def-bundle-r (K562 ∩ Orion ∩ CD4)",
                          pessimistic="199 def-bundle (any donor)"))
    (WORK / "october225.json").write_text(json.dumps(out, indent=2))
    print(f"\n-> {WORK / 'october225.json'}")


if __name__ == "__main__":
    main()
