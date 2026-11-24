"""592b — Draft v16.1: the eight corrections from the 593 cloud audit (experiments/593-cloud-v16-audit/RESULTS.md),
each re-checked locally before writing (see 592 RESULTS.md §5). Anchored replacements on the v16 copy in paper/
(sha b0a6a916… / a36dbb46…, identical to the published experiments/294-writeup v16).
Usage: uv run python experiments/592-manuscript-v16/apply592b.py
"""

from __future__ import annotations

import hashlib
from pathlib import Path

P = Path(__file__).resolve().parent / "paper"


def edit(path: Path, pairs: list[tuple[str, str]], name: str) -> None:
    s = path.read_text()
    for old, new in pairs:
        n = s.count(old)
        if n != 1:
            raise SystemExit(f"{name}: anchor found {n} times:\n{old[:200]}")
        s = s.replace(old, new)
    path.write_text(s)
    print(f"{name}: {len(pairs)} edits, sha256 {hashlib.sha256(s.encode()).hexdigest()[:16]}")


NOTE = """**Draft v16.1 (2026-10-01, after an adversarial audit of v16 by a second model, experiment 593).** The board
arithmetic and the 198-block bound were confirmed. Eight statements were corrected, each re-checked against its file:
H1's programme doses (83 of 85 targets at full dose, not 85), the identity census (five of six wrong calls are documented
duplicates, not all six), the identity floor's quantiles (correct unanimous calls) and its cost (stated per read), the
H1 re-roll formula behind 0.00123, the per-target correlations of a bulk re-roll, 458g's dose wording in §4.1, and the
consensus builds' dose medians (392 is the exception).

"""

PAPER_EDITS = [
    ("**Round 15 revision of 2026-10-01 (Draft v16;", NOTE + "**Round 15 revision of 2026-10-01 (Draft v16;"),
    # revision note
    ("  458d's 900 blocks, against 1.0 for every consensus-prior entry and 1.0 for the same priors on H1. At least 198 of the 648\n",
     "  458d's 900 blocks, against 1.0 for the non-programme builds on the current emitter and 1.0 for the same priors on 83 of\n"
     "  85 H1 targets. At least 198 of the 648\n"),
    ("validated on 402 labelled public profiles from several chemistries. It is\n  right whenever the line is in CCLE; a line",
     "validated on 402 labelled public profiles from several chemistries. It is\n  right on 98.5 % of in-CCLE profiles (99.5 % when the eight methods agree); a line"),
    # §3(h): the limit of the 446b reading (593 attack b)
    ("  channel's re-roll for this emission law (≈ 0.00025 overall; §3(k), C-96). With the day-of key on both sides the October\n",
     "  channel's re-roll for this emission law (≈ 0.00025 overall; §3(k), C-96). The reading rests on three bulk departures\n"
     "  from one file, not on a seed pair, and 446b itself sits 0.00023 below its own parent (433k): its draw is high against its\n"
     "  descendants, not against the whole chain. With the day-of key on both sides the October\n"),
    # §3(j) amendment
    ("dose without any error. The build records give λ's pooled median over the 900 blocks: 1.0 for 446b and every\n"
     "consensus-prior entry, **0.433 for 458d**, 0.243 for 458 (the earlier bank), 0.953 for 458g. On H1 the same programme\n"
     "priors were delivered at λ 1.0 on 85 of 85 targets (462, 467, 517).",
     "dose without any error. The build records give λ's pooled median over the 900 blocks: 1.0 for 446b and for every other\n"
     "non-programme build on the current emitter (one earlier consensus build, 392 on an older emitter, reads 0.901; and a\n"
     "median of 1.0 still leaves single blocks clipped: 446b's minimum is 0), **0.433 for 458d**, 0.243 for 458 (the earlier\n"
     "bank), 0.953 for 458g. On H1 the same programme priors were delivered at λ 1.0 on 83 of 85 targets (462, 467, 517; MED23\n"
     "at 0.85–0.88 and PHF12 at 0.97–0.99 are the exceptions)."),
    # §3(k)
    ("    (on H1, 83–84 of 85 targets, with per-target reach uncorrelated with the base's). Bulk-only board pairs are re-rolls.\n",
     "    (on H1, 83–84 of 85 targets). Per-target reach correlates −0.21 with the base's and nmae +0.15, while fid and jac stay\n"
     "    correlated (+0.42, +0.60) because the call set is prescribed, not drawn. Bulk-only board pairs are re-rolls.\n"),
    ("    by 0.000013 raw: the four DE members together 0.00027 overall, including a known ballast term in nmae. H1's own\n"
     "    re-roll, scaled from 85 to the board's 900 target rows, predicts 0.00123 raw for reach; an eighth pair, preregistered\n",
     "    by 0.000013 raw: the four DE members together 0.00027 overall, including a known ballast term in nmae. H1's\n"
     "    per-target re-roll, taken as a paired difference over the board's 900 target rows (√2 × 0.0261 / √900, from the 82\n"
     "    targets with a finite reach), predicts 0.00123 raw for reach; an eighth pair, preregistered\n"),
    # §4.1
    ("    Both entries delivered the programme well under full dose on these lines (C-97), so what is established is the harm\n"
     "    of the programme as deliverable here, not of a fully delivered one.\n",
     "    458d's pooled median dose was 0.43. 458g's pooled median (0.95) is lifted by its 252 consensus blocks, and at least\n"
     "    198 of the 648 programme blocks it shares with 458d sat at λ ≤ 0.433 (C-97). What is established is the harm of the\n"
     "    programme as deliverable here, not of a fully delivered one.\n"),
    # §5.9
    ("- Every wrong call on an in-CCLE line names a documented duplicate or a line from the same patient: PC-14 = PC-9,\n"
     "  KCI-MOH1 = HPAC, KPL-1 = MCF-7, GP2d / GP5d [CELLOSAURUS]. 146 profiles",
     "- Of the six wrong top-1 calls on an in-CCLE line, five name a documented duplicate or a line from the same patient:\n"
     "  PC-14 = PC-9, KCI-MOH1 = HPAC, KPL-1 = MCF-7, GP2d / GP5d (two profiles) [CELLOSAURUS]. The sixth is a non-unanimous\n"
     "  K562 → COR-L23 call on a 66-gene targeted (TAP-seq) profile, at a lowest top-1 correlation of 0.334. 146 profiles"),
    ("  ≥ 0.68 for 99 % of correct calls and 0.719–0.772 for A/B/C (Flex, 18,400 cells, measured exactly as the final lines\n"
     "  will be).\n",
     "  ≥ 0.68 for 99 % of correct unanimous calls (1st percentile 0.681 on the plain read, 0.678 after the Flex-minus-3′\n"
     "  bias; across every correct call, non-unanimous included, 97 % are ≥ 0.68), and 0.719–0.772 for A/B/C (Flex, 18,400\n"
     "  cells, measured exactly as the final lines will be).\n"),
    ("- A floor of 0.65 on that quantity removes every false switch in this set, at the cost of 2 of 19 correct switches, both\n"
     "  from profiles unlike a challenge context (a partial gene axis; sci-RNA-seq3 or ECCITE chemistry). It does not protect\n",
     "- A floor of 0.65 on that quantity removes both false switches in this set. On the plain read it costs 2 of 19 correct\n"
     "  registry switches (Nadig HepG2, a partial axis, 0.540; sci-RNA-seq3 A549, 0.636); on the Flex-simulated read 2 of 18\n"
     "  (the same HepG2 profile, 0.396; ECCITE THP-1, 0.636). None is shaped like a challenge context. It does not protect\n"),
]

CLAIM_EDITS = [
    ("DE members inside the ordering noise", "DE members inside the re-roll noise (≈ 0.00025 overall, C-96)"),
    ("H1's re-roll scaled by √(85/900) predicts 0.00123 raw",
     "H1's per-target re-roll as a paired difference over 900 rows (√2 × 0.0261 / √900, 82 finite targets) predicts 0.00123 raw"),
    ("pooled median λ over 900 blocks 1.0 for every consensus-prior entry, **0.433** (458d)",
     "pooled median λ over 900 blocks 1.0 for every non-programme build on the current emitter (392, an earlier consensus build on an older emitter, 0.901; medians, not every block: 446b's minimum is 0), **0.433** (458d)"),
    ("the same programme priors at λ 1.0 on 85 / 85 H1 targets",
     "the same programme priors at λ 1.0 on 83 / 85 H1 targets (MED23 0.85–0.88, PHF12 0.97–0.99)"),
    ("is right whenever the line is in CCLE: top-1 **0.985**", "is right on 98.5 % of in-CCLE profiles: top-1 **0.985**"),
    ("every wrong in-CCLE call is a documented duplicate or same-patient line;",
     "five of the six wrong in-CCLE calls are documented duplicates or same-patient lines (the sixth: a non-unanimous call on a 66-gene TAP-seq profile, 0.334);"),
    ("against ≥ 0.68 for 99 % of correct calls and 0.719–0.772 for A/B/C; a floor of 0.65 removes every false switch at 2 of 19 correct ones",
     "against ≥ 0.68 for 99 % of correct unanimous calls (97 % of all correct calls) and 0.719–0.772 for A/B/C; a floor of 0.65 removes both false switches at 2 of 19 correct registry switches (plain read; 2 of 18 Flex-simulated)"),
]


def main() -> None:
    edit(P / "PAPER.md", PAPER_EDITS, "PAPER.md")
    edit(P / "CLAIMS.md", CLAIM_EDITS, "CLAIMS.md")


if __name__ == "__main__":
    main()
