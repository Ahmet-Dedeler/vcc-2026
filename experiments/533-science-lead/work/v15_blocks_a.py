"""533: Draft v15 text blocks, part A (revision note, abstract, plain-language summary, §1.4, §3(h)).

Each block is (anchor, mode, text): mode 'after' inserts text right after the unique anchor, 'before' right before it,
'replace' swaps the anchor for text. apply_v15.py asserts every anchor occurs exactly once in the v14 file.
Every number below is traced in experiments/533-science-lead/RESULTS.md §2.
"""

REVISION_NOTE = """

**Round 14 revision of 2026-09-30 (Draft v15; prepared in experiment 533 from the record of 2026-09-27 … 09-29).** Six
scored entries and experiments 513–527 and 540–546 enter the manuscript. Every number was traced to its primary file, and
the board arithmetic was recomputed from the captured entries (`experiments/533-science-lead/RESULTS.md` §2,
`work/verify_board.py`).
- §3(h): a zero-shot pipeline re-fitted from the panel alone reproduces its tuned self on the scored board, within 0.0005
  overall, once its one composition choice matches (C-82; Figure F-12). A fifth donor line in the bulk reads flat (C-92).
  Dependency does not predict cross-line response size among non-essential targets (C-81). Two field facts: no public source
  names or screens the final lines, and no public emission design explains the leaders' fold-change members (C-90, C-91).
- §3(i): a second pure fold-change-magnitude change transferred at 1.19× its H1 price with the other five members
  bit-identical (C-74). Called genes' magnitudes sit at H1's L1 optimum (C-75). The uncalled-gene add sourced from a
  same-lab programme lowers H1's fold-change error 12× more than the consensus-sourced add, and its bar is the design's own
  signs rather than the consensus (C-76, C-77; Figure F-13). Its board test is preregistered and unread.
- §3(j): target-specific public sign sources are negative on every cancer fold (C-79). Recipient-side, library-size and
  dropout features carry no transferable sign (C-80). A same-lab programme from another line helps where the design is
  weakest (C-78).
- §3(k): the decoupling holds over six more entries (C-83). At its cliff the profile-error member reads the aggregate
  cosine, and the per-target magnitude pot is sized (C-84, C-85, C-93).
- §3(l) is new: the rank member is sign-limited, and per-target call-set sizing is closed by its payoff (F20; C-86 … C-88).
- §4.4 and §4.5 absorb the negatives and the board's verdicts on local prices (T9 … T11). §4.7 (open problems) is new.
- **[[PENDING 09-30]]** The 09-30 entries test F18's open question on the challenge's own lines: the benchmark producer's
  H1 programme as every context's block prior (458d) and its design-strength gate (458g). Their result belongs in §3(j),
  the abstract, the plain-language summary, §4.1, §4.7 and claim C-89; every place is marked.
- Corrections carried: F-9's caption said the DE members stayed "within 0.0006" (0.0007: 503's C-61 correction had not
  reached the figure); F-10's source line named `labline450.json` while the figure plots `labline450_v3.json`; the figure
  index said "all eight"; §1.4 said board numbers are confined to §3(h), which §3(i)–(l) no longer are."""

ABSTRACT_ADD = """ On the benchmark's own scorer, the fold-change members reward sign information on weak and uncalled responses, and none of
the public sources we could assemble supplies it across cancer lines: target-specific combinations of six same-target rows
from four public screens gain +0.016 [+0.007, +0.024] in pooled sign accuracy on the stem-cell reference and −0.0002 …
+0.0033 on three held-out cancer lines. On the stem-cell reference the one source that closed most of that gap was a screen
of the same line from another laboratory. **[[PENDING 09-30:
whether the benchmark producer's own screen of a different line supplies a challenge line's programme on the scored
member.]]**"""

PLAIN_ADD = """

Since the last revision we asked a narrower question about the competition itself: in a cell line nobody has perturbed,
which way do the many weakly affected genes move? Every public source we could combine either fails outright or helps on
the stem-cell reference and not on the cancer lines we could test. That includes other laboratories' screens of the same
genes, the cell line's own unperturbed cells, and technical features of the measurement. The one thing that carried the
direction, on the stem-cell line, was another laboratory's screen of that same line, and the competition's cell lines have
no such public screen. **[[PENDING 09-30: whether the organisers' own screen of a different cell line carries it.]]**"""

S14_OLD = ("numbers from a private competition scorer are engineering calibration, not science, and are confined to §3(h) "
           "and the Discussion.")
S14_NEW = ("numbers from the competition's scorer are engineering calibration, not biology: they are confined to "
           "§3(h)–(l) and the Discussion, and §3(i)–(l) use them only as evidence about what the scorer's members measure.")

S3H_ADD = """

**Later scored entries (2026-09-27 … 09-29): a zero-shot pipeline re-fitted from the panel alone reproduces its tuned self
(C-82; Figure F-12a).** The final test panel arrives with its own targets. Every panel-dependent component must then be
re-fitted from that panel alone: the recipient's identity match, the CD4 donor rows, the whitening operator (a function of
the panel that fails silently on a new one), the per-target bulk key, and the level at which each context sits on the
profile-error member's cliff, set from control cells because no truth exists. That whole path was run by one script on A/B/C
relabelled as new contexts (453dr) and scored.
- At the composition the day-of rules choose (CD4 donor rows at ×1), it cost −0.0036816 overall against the tuned entry
  446b. The bulk members' share was −0.0030866, inside the preregistered band (−0.0045 … −0.0005); the DE share −0.000595 was
  reach, the ordering channel.
- The same path at the tuned composition (CD4 ×1.5; 457dr) recovered +0.0031484, all of it bulk (pds +0.001713, mse
  −0.014313 raw). That is 1.83× its leave-one-library-out (PLOLO) price at k 1 (+0.00172), above the k 1.3 point
  (+0.0022) and inside its preregistered band (+0.0005 … +0.0035).
- Net, the refit at the tuned composition sits −0.0005332 from the tuned entry: bulk +0.0000141, and the rest reach, inside
  the 0.0045 re-roll sd of the ordering channel.
- So on this panel the operator refit, the key refit, the control-proxy re-level and the relabelling cost nothing
  measurable. The only cost was a deliberate composition choice. The day-of rules keep CD4 ×1 because on synthetic panels
  drawn like the final one CD4 ×1.5 loses (438), and the A/B/C gain did not reach the preregistered re-open threshold
  (+0.0040).
- Independent code reviews of the same path found defects that would have failed or silently degraded the final build
  (§5.6).

**A fifth donor line in the bulk reads flat (446bd2; C-92).** DLD-1 (a public colorectal CRISPRi screen covering 67 of the
300 panel targets) was added as a second union residual. The entry read +0.0004200 overall: pds +0.002148 raw, at the top of
its preregistered band, and mse −0.000471 raw against a preregistered −0.001 … −0.004. Flat by its bins; the bulk keeps its
sources. The board thereby confirms a closure the programme had reached on local folds: donor sources saturate.

**Response size does not follow dependency among non-essential targets (514; C-81).** Across lines, the log ratio of a
target's response energy in HCT116 and in K562 does not follow its DepMap dependency difference [DEPMAP-CHRONOS]: r +0.001
over the 7,051
targets non-essential in both (−0.041 over all 8,781 shared). Within a line, dependency tracks response size only through
essential genes (K562 r −0.354 over all targets, −0.062 over non-essential ones). A dependency-based per-(target, line)
magnitude correction is closed at its construction for an essential-free panel.

**Three weeks after the 306 sweep (520, 525, 532; 2026-09-28/29; C-90, C-91).**
- No public document names the final three lines. The organisers' 2026 commentary [VCC26-CELL] and news post [VCC26-NEWS]
  say only that six lines
  from different tissues were screened with CRISPRi on 10x Flex. Every public genome-scale, non-essential cancer CRISPRi
  bank that covers a challenge-style panel is already registered in the pipeline's same-line mode (K562 272, HCT116 290 and
  DLD-1 67 of the A/B/C panel's 300 targets). The lines judged most plausible (A549, MCF7, HepG2, HT-29, BxPC3) have public
  banks covering 0–14 of 300. The organisers' 2025 design note [VCC-BTD] chose its line because K562 and A375 dominate public
  Perturb-seq. No qualifying screen was released between 2026-09-15 and 09-29 among the records opened (520, re-checked by
  532). The final entry therefore rests on general transfer.
- No public emission design explains the leaders' fold-change members. The most reproduced public emitter [ATLASSHIFT] and its
  bit-identical board copies lose all four DE members to this pipeline (raw fid 0.5076, nmae 0.9183, reach 0.1493, jac
  0.0300, against 0.5144, 0.8949, 0.2852 and 0.0408). So does every other opened emitter whose scored numbers are known,
  and one team whose generator is not public. One opened decoder has no mapped board row.
- A cluster of six teams runs a private pipeline: two carry one model string, and four more agree with it at seed scale in
  every member (fid 0.528, nmae 0.881, jac 0.033). The string is absent from every repository opened, from PyPI and from the
  Hugging Face index. Board descriptions are self-reports and are treated as such."""

ANCHORS_A = [
    # the title line: insert the newest revision note right after it
    ("# Three reproducibilities: what a single-cell perturbation reference measures, which responses are conserved across "
     "cell contexts, and what bounds zero-shot transfer", "after", REVISION_NOTE),
    ("Finally we show why a rank metric rewards exactly the directions in which a fitted whitener has learnt its own "
     "training noise.", "after", ABSTRACT_ADD),
    ("including eight datasets we had deliberately never opened until the test was written down.", "after", PLAIN_ADD),
    (S14_OLD, "replace", S14_NEW),
    ("and the public constructions that can be read use our four sources and land 0.05–0.08 below us.\"", "after", S3H_ADD),
]
