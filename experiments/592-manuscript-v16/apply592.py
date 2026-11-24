"""592 — Draft v16 edits to the manuscript COPY in experiments/592-manuscript-v16/paper (anchored, each exactly once).

Run once on a fresh copy of experiments/294-writeup's PAPER.md / CLAIMS.md / REFERENCES.md (v15.1, sha 9aeefb62… /
83d3c7e7…). Every number below is traced in RESULTS.md §2. Usage: uv run python experiments/592-manuscript-v16/apply592.py
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


REV = """**Round 15 revision of 2026-10-01 (Draft v16; prepared in experiment 592 from the record of 2026-09-30 … 10-01).** Two
scored entries, three corrections of the record and one methods validation enter. Every number was traced to its primary
file (`experiments/592-manuscript-v16/RESULTS.md` §2).
- §3(h): the October path's day-of leave-one-donor-out key transfers (460dr, +0.0028694 overall against 453dr; mse at
  0.90× the A/B/C key's own board effect; C-94). With the key on both sides, the October path and the tuned path agree in
  the four DE members to 0.0000139 overall (460dr − 446bj2). The reach the previous draft charged to the refit path was
  most likely 446b's own realisation (C-82 amended).
- §3(j): C-89 amended. The producer's programme reached the board at under half its dose: a pooled median λ 0.433 over
  458d's 900 blocks, against 1.0 for every consensus-prior entry and 1.0 for the same priors on H1. At least 198 of the 648
  programme blocks that 458g shares with 458d sat at λ ≤ 0.433, so the gated entry does not test full dose either. The
  read establishes that the programme, as this emitter can deliver it on these lines, hurt; whether a fully delivered
  programme would hurt is open (C-97).
- §3(k): the DE channel's read noise corrected (C-96). The 0.0045 overall quoted since September is the spread of six
  different emission laws (230), not the re-roll of one design. For the law in use since 410b a re-roll moves the four DE
  members by about 0.00025 overall, measured on seven bulk-only board pairs and matched by H1's re-roll scaled to the
  board. A bulk-only change is itself a re-roll: the ballast count shifts the file's single emission stream.
- §3(l): the call-propensity block alone read NULL on the board (446bj2zcp: Δ(fid + jac) −0.00018 overall against an
  H1-derived +0.0020 … +0.0025). Its swapped-in genes were called (jac +0.0040 raw) but carried no sign steering (fid
  −0.0036 raw; C-95; T12 in §4.5).
- §5.9 is new: the day-of line identification, validated on 402 labelled public profiles from several chemistries. It is
  right whenever the line is in CCLE; a line absent from CCLE is still called unanimously, at a lower match quality that a
  similarity floor separates (C-98).
- The abstract, the plain-language summary, §4.1 and §4.7 carry the dose caveat on C-89.
- Pending, with markers: the consensus uncalled-gene hook (496a, 10-05), the K pair (523, from 10-06), the block + hook arm
  (446bj2zcph, 10-07).

"""

PAPER_EDITS = [
    ("**Round 14 revision of 2026-09-29, with the 09-30 read integrated", REV + "**Round 14 revision of 2026-09-29, with the 09-30 read integrated"),
    # abstract
    ("of the same line from another laboratory. The benchmark producer's own screen of a different line hurt (fid −0.0055 raw,\nthe preregistered harm bin).",
     "of the same line from another laboratory. The benchmark producer's own screen of a different line hurt (fid −0.0055 raw,\nthe preregistered harm bin), as delivered: on these lines the pipeline's emitter carried it at under half its dose."),
    # plain-language summary
    ("no such public screen. The organisers' own screen of a different cell line made it worse.",
     "no such public screen. The organisers' own screen of a different cell line made it worse, at least at the reduced\nstrength our pipeline could deliver it on these lines (under half)."),
    # §3(h)
    ("- Net, the refit at the tuned composition sits −0.0005332 from the tuned entry: bulk +0.0000141, and the rest reach, inside\n  the 0.0045 re-roll sd of the ordering channel.\n",
     "- Net, the refit at the tuned composition sits −0.0005332 from the tuned entry: bulk +0.0000141, and the rest (−0.0005473)\n"
     "  reach. That reach is most likely 446b's own realisation, not a cost of the path (corrected 2026-10-01). Every later entry\n"
     "  that kept 446b's DE design and changed only the bulk or the path reads 0.0004–0.0006 overall below it in the four DE\n"
     "  members (446bj2 −0.00039, 446bd2 −0.00048, 453dr −0.00060, 457dr −0.00055), about two standard deviations of the DE\n"
     "  channel's re-roll for this emission law (≈ 0.00025 overall; §3(k), C-96). With the day-of key on both sides the October\n"
     "  path and the tuned path agree in the four DE members to 0.0000139 overall (460dr − 446bj2), and their whole difference is\n"
     "  the CD4 composition (bulk −0.0030321).\n"
     "- The day-of leave-one-donor-out key itself transfers on the path (460dr − 453dr, 2026-10-01; C-94): +0.0028694 overall,\n"
     "  mse −0.013060 raw (0.90× the A/B/C key's own board effect, 446bj2 − 446b, −0.014473) and pds +0.001044 raw, with the DE\n"
     "  members inside their preregistered 3σ noise bands (|Δfid| ≤ 0.0001, |Δjac| ≤ 0.00005, |Δreach| ≤ 0.0040 raw).\n"),
    # §3(j) amendment
    ("(`experiments/458-arc-lab-programme/BOARD-RESULT-458d-458g.md`).\n\n*A benchmark fact recovered from public submissions.*",
     "(`experiments/458-arc-lab-programme/BOARD-RESULT-458d-458g.md`).\n\n"
     "*Amendment (2026-10-01): the programme reached the board at under half its dose (C-97).* The emitter bisects a block\n"
     "prior's log2 move (λ) down until every library class keeps a fixed 6,000-count margin, so a prior can arrive below full\n"
     "dose without any error. The build records give λ's pooled median over the 900 blocks: 1.0 for 446b and every\n"
     "consensus-prior entry, **0.433 for 458d**, 0.243 for 458 (the earlier bank), 0.953 for 458g. On H1 the same programme\n"
     "priors were delivered at λ 1.0 on 85 of 85 targets (462, 467, 517). The likely cause is the programme's balanced signs\n"
     "(about 56 % down, against the consensus's 84 %), which need more class capacity than these lines' sub-gate mass gives.\n"
     "- So 458d compared a programme at roughly 0.43 of its dose with the consensus at full dose; part of its −0.0085 overall\n"
     "  may be lost steering rather than wrong signs.\n"
     "- 458g does not settle it. Its programme blocks are byte-identical to 458d's on the same targets (the exact decomposition\n"
     "  above), and its pooled median is lifted by the 252 consensus blocks of its weak third. Because 458d's median is 0.433\n"
     "  and its weak third holds only 252 of its 900 blocks, at least 198 of the 648 strong-target programme blocks sat at\n"
     "  λ ≤ 0.433 in both entries.\n"
     "- What the read establishes: the producer's H1 programme, as this emitter can deliver it on these three lines, hurt fid,\n"
     "  reach and nmae. Whether a fully delivered programme would carry is open. A dose-only rebuild with the emitter's\n"
     "  per-block audit (the scored file unchanged) would measure the doses exactly; no October decision depends on it, since\n"
     "  the programme is off.\n"
     "(`experiments/458-arc-lab-programme/BOARD-RESULT-458d-458g.md`, addendum; `experiments/554-dose-guard/RESULTS.md`; the\n"
     "build metas' `output.rc_callset.prior_lambda`.)\n\n"
     "*A benchmark fact recovered from public submissions.*"),
    # §3(k)
    ("- The four bulk-only changes moved nmae by at most 0.0002 raw, fid by at most 0.00004 and jac by at most 0.000016. Reach\n"
     "  moved by up to 0.0034 raw. Reach is the ordering channel, whose re-roll sd between digest-identical entries is about 0.027\n"
     "  in scaled reach (0.0045 overall), and a bulk change re-draws which genes head the |z| order.\n",
     "- The four bulk-only changes moved nmae by at most 0.0002 raw, fid by at most 0.00004 and jac by at most 0.000016. Reach\n"
     "  moved by up to 0.0034 raw (453dr against 446b; §3(h) says why most of that sits on 446b's own draw).\n"
     "- *How large is the DE channel's read noise? (corrected 2026-10-01; 537 §2, 547, 548; C-96).* The 0.0045 overall this\n"
     "  record quoted since September as the re-roll sd \"between digest-identical entries\" is the sample sd of 230's six pairs,\n"
     "  and every one of those pairs changed the emission law (overdispersion, resampling, the |z| assignment, the call set). It\n"
     "  is the prior spread of a change of law, not the re-roll of one design. For the law in use since 410b the re-roll is\n"
     "  measured directly.\n"
     "  - The builder draws one emission stream per file, and each block's ballast count sets how many carriers that block takes\n"
     "    from it. A bulk change that moves any block's count therefore re-draws the first emission pass for every later block\n"
     "    (on H1, 83–84 of 85 targets, with per-target reach uncorrelated with the base's). Bulk-only board pairs are re-rolls.\n"
     "  - Over seven of them reach moved by rms 0.00130 raw [0.00086, 0.00264] (0.00024 overall), fid by 0.000024 raw and jac\n"
     "    by 0.000013 raw: the four DE members together 0.00027 overall, including a known ballast term in nmae. H1's own\n"
     "    re-roll, scaled from 85 to the board's 900 target rows, predicts 0.00123 raw for reach; an eighth pair, preregistered\n"
     "    before it was scored (460dr − 453dr), fell inside its 3σ bands.\n"
     "  - DE knobs that leave the ballast counts alone keep the first pass and read at 0.3–0.6× a re-roll; prior-only and\n"
     "    post-pass knobs share every draw and carry no re-roll at all (C-62, 468).\n"),
    # §3(l)
    ("*Why it matters for benchmark design.* Two of the four DE members turn on quantities a model cannot learn from public data.",
     "*The call set's composition on the board (2026-10-01; 528, 535; C-95).* The call set's size is fixed by design; which\n"
     "genes fill it is a choice. On H1, choosing the 2,200-gene block by call propensity (control expression plus excess\n"
     "dropout, from each context's own controls), with the consensus sign steering kept only on retained high-expression genes,\n"
     "priced +0.0051 [+0.0033, +0.0070] overall-equivalent (3 seeds; jac +0.0091, reach +0.0096, fid −0.0013 raw), and the\n"
     "preregistered board prediction was Δ(fid + jac) +0.0020 … +0.0025 overall. The board read NULL (446bj2zcp − 446bj2z):\n"
     "- Δ(fid + jac) −0.00018 overall. jac +0.003986 raw (0.44× H1's), fid −0.003580 raw (2.8× H1's loss), reach −0.003353 raw\n"
     "  (the opposite sign); pds identical and mse within 2.3 × 10⁻⁸ raw (a few whole counts moved with the block's emission):\n"
     "  a DE-only change for every practical purpose.\n"
     "- The mechanism is visible: the ~550–680 genes swapped into each context's block are called (jac up) but carry no sign\n"
     "  steering (about 0.51 sign accuracy on H1), so the extra calls dilute fid's directional precision. A version that gives\n"
     "  the swapped-in genes the consensus sign through the uncalled-gene hook is built and scheduled\n"
     "  [[PENDING 10-07: 446bj2zcph, the block + hook]].\n"
     "- H1 over-read jac by 1/0.44 and missed the fid loss: the size of error a single-line surface made on a change of which\n"
     "  genes fill a fixed-size call set (T12).\n\n"
     "*Why it matters for benchmark design.* Two of the four DE members turn on quantities a model cannot learn from public data."),
    # §4.1
    ("every context's block prior hurt (fid −0.0055 raw); its design-strength gate recovered part of the loss (+0.0036), §3(j).\n",
     "every context's block prior hurt (fid −0.0055 raw); its design-strength gate recovered part of the loss (+0.0036), §3(j).\n"
     "    Both entries delivered the programme well under full dose on these lines (C-97), so what is established is the harm\n"
     "    of the programme as deliverable here, not of a fully delivered one.\n"),
    # §4.5 T12
    ("  the board (417c: −0.0068 predicted on H1, +0.0089 raw on the board). This is not a contradiction inside 516: its\n"
     "  cross-study cancer folds read the same ≈ 0 that the board implies.\n",
     "  the board (417c: −0.0068 predicted on H1, +0.0089 raw on the board). This is not a contradiction inside 516: its\n"
     "  cross-study cancer folds read the same ≈ 0 that the board implies.\n"
     "- **T12 — two more reads against their local prices (2026-10-01).**\n\n"
     "  | read | class | local price | board | verdict |\n"
     "  |---|---|---|---|---|\n"
     "  | 460dr, day-of LODO key on the October path | bulk magnitude key | the A/B/C key's board effect: mse −0.014473 raw | mse −0.013060; pds +0.001044 | transferred, mse 0.90× |\n"
     "  | 446bj2zcp, call-propensity block | DE call-set composition | Δ(fid + jac) +0.0020 … +0.0025 overall (H1 × 0.72) | −0.00018 | NULL: jac 0.44×, fid 2.8× H1's loss, reach of opposite sign |\n\n"
     "  The block's H1 measurement is not contradicted (three seeds, interval clear of zero). What failed is H1 as an instrument\n"
     "  for which genes fill a fixed-size call set: the swapped-in genes are steered by nothing, and how much that costs fid\n"
     "  depends on how many of them the reference calls, which H1 does not share with the cancer lines.\n"),
    # §4.7 item 3
    ("   It did not (458d, 09-30: the producer's H1 programme hurt on the challenge's cancer lines): a line's programme is\n",
     "   It did not (458d, 09-30: the producer's H1 programme hurt on the challenge's cancer lines, as delivered at under half\n"
     "   its dose; C-97): a line's programme is\n"),
]

SEC59 = """

**5.9 Line identification on the day, and its validation (590; C-98).** A new context's line is identified from its
control cells alone. Their pseudobulk (log2 CPM on the official ∩ CCLE axis, 17,900 genes) is matched to the 1,600 CCLE
bulk profiles (Broad 2025 release, DepMap Public 25Q3) by four feature sets × Pearson / Spearman. The eight rankings are
combined by a Borda count weighted by how well each method ranks three calibration anchors (K562, HCT116 and Jurkat
controls from public screens). A call switches same-line donor rows on only when it names a registered line and all eight
methods agree.
- Validation on 402 labelled public profiles, 390 of them of lines in CCLE (200 MIX-seq lines [MIXSEQ20]; sci-Plex A549,
  MCF7 and K562; Calu-3, THP-1, Jurkat, K562 and A549 screens; donor control profiles; 10x 3′ and 5′, sci-RNA-seq3,
  CROP-seq, scifi-RNA-seq, TAP-seq and Flex): top-1 accuracy 0.985, and 0.995 when the eight methods agree.
- Every wrong call on an in-CCLE line names a documented duplicate or a line from the same patient: PC-14 = PC-9,
  KCI-MOH1 = HPAC, KPL-1 = MCF-7, GP2d / GP5d [CELLOSAURUS]. 146 profiles of unregistered lines that share a lineage with a
  registered one produced no false switch.
- A line absent from CCLE is still called unanimously: both HEK293 profiles were matched to HAP1, a registered line. Such
  calls match less well. The lowest top-1 correlation over the eight methods is 0.51–0.61 for them, also after a per-gene
  Flex-minus-3′ bias (fitted on the challenge's own Flex controls, cross-fitted over the MIX-seq lines) is added. It is
  ≥ 0.68 for 99 % of correct calls and 0.719–0.772 for A/B/C (Flex, 18,400 cells, measured exactly as the final lines
  will be).
- A floor of 0.65 on that quantity removes every false switch in this set, at the cost of 2 of 19 correct switches, both
  from profiles unlike a challenge context (a partial gene axis; sci-RNA-seq3 or ECCITE chemistry). It does not protect
  against a line absent from CCLE whose nearest CCLE line is a close relative: primary T cells match a T-cell lymphoma line
  at 0.719."""

CLAIM_ROWS = r"""| **C-94** | The October path's day-of leave-one-donor-out key transfers: 460dr − 453dr = **+0.0028694** overall (mse **−0.013060** raw, 0.90× the A/B/C key's board effect −0.014473; pds +0.001044 raw); DE members inside their preregistered 3σ noise bands; at matched key the October and tuned paths agree in the four DE members to **0.0000139** overall (460dr − 446bj2; bulk −0.0030321, the CD4 composition) | `E460/BOARD-RESULT-460dr.md`; `E460/PREREGISTRATION.md`; `E537/RESULTS.md` §2.5; `data/board-entries/{LrkhsIXWMjEsvRXSGToa,q3VuTImh4FRGKXH2waCp,YvVru3DAF89Jiyr34olR}.json`; `E592/RESULTS.md` §2 (recomputed) | board entries · 3 | exact (board); DE re-roll sd ≈ 0.00025 overall (C-96) | not graded (engineering validation) | **preregistered**: bin "same sign, outside the band: keep it ON" (mse 0.90× inside 0.7–1.3; pds 5.0× outside) | a re-read of the day-of key whose mse effect falls below 0.7× the A/B/C key's |
| **C-95** | The call set's composition (a call-propensity block from each context's own controls, consensus steering kept on retained high-expression genes) read **NULL** on the board: 446bj2zcp − 446bj2z Δ(fid + jac) = **−0.00018** overall against a preregistered +0.0020 … +0.0025; jac +0.003986 raw (0.44× H1), fid −0.003580 (2.8× H1's loss), reach −0.003353 (opposite sign to H1), pds identical, mse within 2.3e-8 raw; mechanism: swapped-in genes called but unsteered (~0.51 sign accuracy on H1) | `E535/BOARD-RESULT-446bj2zcp.md`; `E535/PREREGISTRATION.md`; `E528/RESULTS.md`; `data/board-entries/{K1Q0y2OFYjTNRqh9jKfk,be90smZUDcIV7j5MgQ03}.json` | board entries · 2; H1 85 targets × 3 seeds | exact (board); H1 +0.0051 [+0.0033, +0.0070] | NEW (negative on the board, positive on H1) | **preregistered**: NULL bin (−0.0004 … +0.0004) | the block + hook arm (446bj2zcph) reading Δ(fid + jac) ≥ +0.0012 would place the loss on the missing steering, not the composition [[PENDING 10-07]] |
| **C-96** | The DE channel's board read noise for the K-class ballast emission law (410b onward) is ≈ **0.00025 overall** per re-roll. Bulk-only board pairs are re-rolls (a moved ballast count shifts the file's single emission stream; on H1 it re-draws pass 1 on 83–84 of 85 targets). Over seven pairs: reach rms **0.00130** raw [0.00086, 0.00264] (0.00024 overall), fid 0.000024, jac 0.000013 (DE total 0.00027 overall incl. a ballast nmae term); H1's re-roll scaled by √(85/900) predicts 0.00123 raw; a preregistered eighth pair fell inside its 3σ bands. The 0.0045 quoted before is the sd of six emission-**law** changes (230) | `E537/RESULTS.md` §2.1–§2.6 (`noise537b.py`, `redraw537.py`, `work/noise537b.json`, `work/redraw537.json`); `E547/RESULTS.md` (with its correction); `E548/RESULTS.md`; `E230/work/fit230.json` | board pairs · 7 (+1 preregistered); H1 85 targets × 3 seeds | χ² 95 % interval of the rms | NEW (a property of this pipeline and scorer) | 230's 0.0045 reproduced to 1e-12 as the known answer; the preregistered check (537 §2.5): 1 of 7 pairs read, inside all bands | a builder-seed re-roll pair on the board moving the DE members by > 0.001 overall |
| **C-97** | A block prior can reach the board below full dose without error (the emitter bisects its log2 move until every library class keeps a 6,000-count margin): pooled median λ over 900 blocks 1.0 for every consensus-prior entry, **0.433** (458d), **0.243** (458), 0.953 (458g, lifted by its 252 consensus blocks); the same programme priors at λ 1.0 on 85 / 85 H1 targets; at least **198 of the 648** strong-target programme blocks shared by 458d and 458g at λ ≤ 0.433 (from 458d's median and the exact decomposition of C-89) | build metas `submissions/{446b,458,458d,458g}-…h5ad.meta.json` (`output.rc_callset.prior_lambda`); `E458/BOARD-RESULT-458d-458g.md` (addendum); `E554/RESULTS.md` (calibration); `E468/RESULTS.md` (compare458g); `E592/RESULTS.md` §2 | build records · 4 (+ the H1 diagnostics of 462 / 467 / 517) | pooled min / median / max only (no per-block dump exists for these builds) | NEW (a delivery property of this emitter) | not preregistered (found by a calibration) | a per-block dose dump of 458d showing its strong-target blocks at λ ≥ 0.9 |
| **C-98** | Line identification from control cells (041's CCLE matcher with Oct 22's weights) is right whenever the line is in CCLE: top-1 **0.985** over 390 labelled public profiles (**0.995** when the 8 methods agree); every wrong in-CCLE call is a documented duplicate or same-patient line; 0 false same-line switches over 146 same-lineage unregistered lines; a line absent from CCLE is still called unanimously (both HEK293 profiles → HAP1, registered) at a lowest top-1 correlation of **0.51–0.61** (plain and Flex-adjusted), against ≥ 0.68 for 99 % of correct calls and 0.719–0.772 for A/B/C; a floor of 0.65 removes every false switch at 2 of 19 correct ones | `E590/RESULTS.md`; `E590/PREREGISTRATION.md`; `E590/work/{match590.csv,match590_flexbias_fold0.csv,match590_flexbias_fold1.csv,summary590.json}` | profiles · 402 (390 in CCLE; 200 MIX-seq lines) | exact counts (a census; no sampling interval) | KNOWN in part (expression-based cell-line authentication); the unanimous call of an absent line and its floor are candidate NEW for this use | **preregistered** (590); known answer exact (A/B/C reproduce `def_identity.match`) | a profile of a line absent from CCLE matched unanimously to a registered line at ≥ 0.65 |
"""

CLAIM_EDITS = [
    ("net 457dr − 446b = **−0.0005332** (bulk +0.0000141) |",
     "net 457dr − 446b = **−0.0005332** (bulk +0.0000141); amended 2026-10-01: that DE share (−0.0005473) is most likely 446b's own realisation (every later entry keeping 446b's DE design reads 0.0004–0.0006 below it in DE; at matched key the October and tuned paths agree in DE to 0.0000139, C-94) |"),
    ("| board entries · 3 | exact; ordering-channel re-roll sd 0.0045 overall |",
     "| board entries · 3 | exact; DE re-roll sd ≈ 0.00025 overall (C-96; the 0.0045 once quoted here is the spread of six emission-law changes) |"),
    ("| a ballast-carried bulk arm whose DE members move beyond the ordering channel's re-roll noise |",
     "| a ballast-carried bulk arm whose DE members move beyond the DE re-roll noise (≈ 0.00025 overall, C-96) |"),
    ("**Read 2026-09-30: HARM** (fid −0.0055 raw; the gate +0.0036 vs 458d; strong two thirds −0.0049, weak third −0.0036) |",
     "**Read 2026-09-30: HARM** (fid −0.0055 raw; the gate +0.0036 vs 458d; strong two thirds −0.0049, weak third −0.0036); **amended 2026-10-01**: delivered at a pooled median λ 0.433 (458d) against 1.0 for the consensus, with at least 198 of the 648 strong-target programme blocks shared with 458g at λ ≤ 0.433, so established for the programme as this emitter delivers it on these lines, not at full dose (C-97) |"),
    ("| a second producer screen of a cancer line carrying where H1's did not; an emission re-realisation moving fid by > 0.004 |",
     "| a second producer screen of a cancer line carrying where H1's did not; an emission re-realisation moving fid by > 0.004; the same programme delivered at full dose (per-block λ ≥ 0.9) reading null or positive |"),
    ("| board entries · 2 | exact; ordering sd 0.0045 overall |",
     "| board entries · 2 | exact; DE re-roll sd ≈ 0.00025 overall (C-96) |"),
]

REFS = """

## Added for Draft v16 (592, 2026-10-01)

| key | reference |
|---|---|
| **CELLOSAURUS** | Bairoch A. *The Cellosaurus, a cell-line knowledge resource.* J Biomol Tech 29:25–38 (2018). Records read via api.cellosaurus.org on 2026-10-01 [F]: CVCL_1640 PC-14 ("Problematic cell line: Contaminated. Shown to be a PC-9 derivative"), CVCL_2090 KCI-MOH1 (an HPAC derivative), CVCL_2094 KPL-1 (an MCF-7 derivative), CVCL_2450 GP2d (same individual as CVCL_1235 GP5d). §5.9. |
| **CCLE-2025** | Cancer Cell Line Encyclopedia, Broad 2025 release as mirrored by cBioPortal (`ccle_broad_2025`, DepMap Public 25Q3): `data_mrna_seq_tpm.txt` (1,600 lines with expression) and `data_clinical_sample.txt`; on disk under `data/external/ccle_broad_2025` (cache sha in `ccle_cache.meta.json`). §5.9's reference panel. |
"""

MANIFEST_EDITS = [
    ("with bulk +0.00001: the DE remainder is reach, inside the ordering channel's re-roll sd of 0.0045 overall.",
     "with bulk +0.00001: the DE remainder is reach, most likely 446b's own realisation (every later entry keeping its DE design reads 0.0004–0.0006 below it in DE; with the day-of key on both sides the October and tuned paths agree in DE to 0.00001, §3(h))."),
    ("Exact board values, no intervals; the ordering-channel noise between digest-identical arms is sd 0.0045 overall.",
     "Exact board values, no intervals; the DE channel's re-roll between two builds of this emission law is ≈ 0.00025 overall (C-96; the 0.0045 quoted in earlier drafts was the spread of six emission-law changes)."),
]


def main() -> None:
    edit(P / "PAPER.md", PAPER_EDITS, "PAPER.md")
    s = (P / "PAPER.md").read_text()
    i = s.index("**5.8 Data and code availability.**")
    j = s.index("\n\n---", i)
    s = s[:j] + SEC59 + s[j:]
    (P / "PAPER.md").write_text(s)
    print(f"PAPER.md: §5.9 inserted, sha256 {hashlib.sha256(s.encode()).hexdigest()[:16]}")
    edit(P / "CLAIMS.md", CLAIM_EDITS, "CLAIMS.md")
    c = (P / "CLAIMS.md").read_text()
    k = c.index("| **C-93** |")
    e = c.index("\n", k) + 1
    c = c[:e] + CLAIM_ROWS + c[e:]
    (P / "CLAIMS.md").write_text(c)
    print(f"CLAIMS.md: C-94 … C-98 inserted, sha256 {hashlib.sha256(c.encode()).hexdigest()[:16]}")
    r = (P / "REFERENCES.md").read_text().rstrip("\n") + REFS
    (P / "REFERENCES.md").write_text(r)
    print(f"REFERENCES.md: 2 entries added, sha256 {hashlib.sha256(r.encode()).hexdigest()[:16]}")
    edit(P / "figs/MANIFEST.md", MANIFEST_EDITS, "figs/MANIFEST.md (F-12 caption)")


if __name__ == "__main__":
    main()
