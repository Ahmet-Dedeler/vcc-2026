# 310 — RESULTS: Round 6 folded into the manuscript (F16, 308, 306); page rebuilt, nothing published

**2026-09-07.** Editorial only. No measurement was made, no `vcc`, no `cell-eval2`, no build, no
submission, no download, nothing under `submissions/` or `data/`, **nothing published**. `uv run python`
throughout (one invocation of `experiments/294-writeup/build_html.py`, one 6-line in-place reorder of a
markdown table). Files written: the three authorised manuscript files and this directory. Verified by
`git status`: `M experiments/294-writeup/{PAPER.md,CLAIMS.md,REFERENCES.md}` plus
`?? experiments/310-writeup-round6/work/` and nothing else attributable to this agent.

---

## THE LINE

The draft now carries Round 6. **§3(g) is a new results subsection** for F16 — the aggregate
seed-match repression signal with the **construct** as the unit — carrying C8 (`UNRESOLVABLE`, aborted
by 2 of its own 7 known-answer rows) and C9 (**PASS**, `S` = **+0.71293 [+0.32419, +1.07351]**,
p = **1.99998e-05**, 9/9 rows) with the bias repair's **−9.5 %** stated as declared-in-advance, the
partner-construct test as the load-bearing design element, and "two tests of one hypothesis on one data
set — **not an independent replication**" in the body. **§4.5's T7 is rewritten**: 301's negative stands
(4/13 vs 4/21, reproduced exactly by 307's KA-2), 305 produced no accepted discovery (0 of 11 pass Holm,
p = 0.10476, and its design could not attain its own FWER — attainable minimum 0.0506 > 0.05), and F16
supplies the aggregate measurement, with one sentence saying **why 301's negative and F16 are not in
conflict**. **§3(h)** (the renamed old §3(g)) carries 308's closure *with its mechanism* and 306's
reconnaissance negative about the field. `CLAIMS.md` gains rows **C-55, C-56, C-57** and a Round 6
addendum; `REFERENCES.md` gains **[HORLBECK16]** and the Excel-escaping defect. The page rebuilds at
1,612 KB, 26 TOC entries, all 8 figure captions found.

**140 provenance entries were logged in `work/verify310.log`, each naming the source file, the line and a
verbatim anchor from that line; 138 of them carry measured values and 2 are structural counts. Seven
candidate readings were considered and rejected — they are listed in §5.**

---

## 1. What changed, file by file

### 1.1 `experiments/294-writeup/PAPER.md`

| # | where | change |
|---|---|---|
| 1 | masthead, after the Draft v8 note | **New Round 6 revision note** — names §3(g) as new, says F16 does not establish off-target binding, and points at §3(h) for 308 and 306. Marked "local; not republished". |
| 2 | Abstract | **One new sentence** carrying F16 bounded: construct as the unit (n = 34), `S` = +0.713 [+0.324, +1.074], p = 2.0e-5 against a 1.0e-5 floor, the partner-construct contrast, and — in the same sentence — "does not establish off-target binding" and "two tests of one hypothesis on one data set". |
| 3 | Plain-language summary | **Four new sentences** in the paper's plain register: the seed match, the partner-guide contrast, "not proof that the guide sticks to those genes", and the sequence-mutant experiment as what would settle it. |
| 4 | §1.4 | fourteen findings → **sixteen (F1–F16)**; five ledger claims → **nine (C1–C9)** (both counted from the source files, §G of the log); the `§3(g)` cross-reference for scorer numbers → **§3(h)**. |
| 5 | §2 Data and instruments | **New fifth measurement surface**: the 2025 challenge library's own cells — 65 constructs (34 targeting over 17 paired intended targets, 31 non-targeting), three splits sharing one pool of 38,176 control cells, read only through saved integer sufficient sums. |
| 6 | **§3(g) (new)** | **"What a construct does besides its intended gene"** — F16 in eight paragraphs: the question and the unit; the null and its attainable floor against 305's; C8 and its self-abort; C9's two repairs and the −9.5 %; the partner-construct test and the real-spacer null; three awkward facts (the seed-length ordering, the absent failing-vs-working contrast, the nominally positive non-targeting constructs); the multiplicity, stated plainly; the limits **in the body**; the decisive next experiments with their costs; and the Excel-escaping defect. |
| 7 | §3(g)→**§3(h)** | old §3(g) renamed. Its internal "§3(a)–(f)" updated to "§3(a)–(g)"; its T7 sentence now says the per-construct question is still open while §3(g) measures the aggregate. |
| 8 | §3(h), "what has been measured not to move it" | **308 added to the closed list, then a full paragraph** stating it as a closure with a mechanism: the premise holds (|kd| build-time observable, defined for 226/290/294 targets, absent from the base, ρ +0.03/+0.03/−0.005 with the norm and +0.04…+0.19 with reliability, acting set 106/300 — so not a power problem) and the arm still fails (0 of 16 arms, 51 cells, best +0.00077 [−0.00114, +0.00283], priced −0.00145 pessimistic vs +0.0089, indistinguishable from entropy-matched random re-weighting, 64/64 KA rows), with the measured reason: the higher-\|kd\| donor's row is **not** the one closer to the held truth (down > up in all six cut × fold cells, 1.5–2.4×; 25 up vs 53 down on rel/HEK293T, mean Δ`pds` −0.031). |
| 9 | §3(h), leaderboard paragraph | the live 15:48:55Z board read added (leader 0.2393600, ours 0.2160944 rank 3, gap +0.0232656), then **two new paragraphs from 306**: no qualifying dataset became downloadable between 2026-08-30 and 2026-09-07 (`usedStorage` exactly 0 on the one atlas that would matter; six GEO queries, ArrayExpress, Zenodo, scPerturb all empty; zero acquisitions clear), and the frontier is silent (769 entries, 292 descriptions, **0 of the 2 entries ≥ 0.220 describe a method**, best describing entry rank 42 at 0.17543 and it is an emission recipe; the four source-naming entries score 0.13085–0.14954 and name only sources already on disk; the "21/300" assertion is not reproducible, the defensible range being 20–25). **Rival descriptions are labelled unverified self-reports.** |
| 10 | §4.1 | **New paragraph**: F16 postdates both audits and carries no provenance verdict; the phenomenon is published ([SEED24], [HARTMAN26]) and what is ours is the construct-as-unit test with this null and this control — "until that has been audited … that is all we claim". |
| 11 | §4.4 | **New "not claimed" bullet** for F16: not off-target binding, n = 34, no biological replicate, the partner control cannot separate seed binding from a spacer-correlated artefact, and the failing-vs-working contrast stays unsupported. |
| 12 | **§4.5 T7 — rewritten in full** | three numbered parts: (i) 301's negative **stands** and the aggregate re-reading agrees; (ii) 305 produced **no accepted discovery** — 0 of 11 at Holm 0.10476 and an unattainable FWER (0.0506 > 0.05 on the 8 testable pairs at a 0.0063 rank-p floor); (iii) F16 supplies what 305 could not. Then the required sentence: **301's negative and F16 answer different questions and are not in conflict.** Closes with what remains open (no occupancy, no replicate, no separation of artefact from binding, no 2026 rate) and the sequence-mutant control as the next step. |
| 13 | §4.6 Limitations | the small-n list now ends with **"34 constructs — one experiment's cells, one shared control pool, no biological replicate — behind F16, where the interval is set by the number of constructs"**. |
| 14 | §5.6 | C8 (16:00:01Z) and C9 (16:27:26Z) added to the preregistration-hash list; **a new sentence on the known-answer discipline overriding a positive** — C8 met both pass conditions and was still reported unresolvable — including that the unchanged no-signal row governs in C9 by reproduction, not by level. |
| 15 | §5.8 | the two public sequence tables added (Replogle Table S1; [HORLBECK16] supp. 3 and 2, CC BY 4.0) as the source of the 209,070-protospacer rescan **and of the Excel defect**; the 307/309 derived tables added to the derived-object list; the reserve rule now records that **C8 and C9 declared and spent no reserve**. |

### 1.2 `experiments/294-writeup/CLAIMS.md`

- **C-55 (claim C8)** — the bounded estimate, its unit and n, the null's construction and floor against
  305's 0.0696, the two failed known-answer rows, and `UNRESOLVABLE-AT-THIS-RESOLUTION` as the
  confirmation status. Provenance column says plainly: **postdates both audits, no verdict**, the
  phenomenon is published, the construct-as-unit test on this library is not.
- **C-56 (claim C9)** — PASS by its frozen rule, all five headline numbers, the −9.5 % repair declared in
  advance, the partner test, the sequence-resampling null, and the multiplicity in the confirmation
  column: **not an independent replication**, Bonferroni satisfied arithmetically, the real multiplicity
  analytic, and the wrong-recipient permutation lost as an independent check.
- **C-57** — a *not claimed* row in the file's existing convention: no off-target binding, no
  seed-length-graded mechanism (L\* = 11 at +0.963 above L\* = 12 at +0.685), no failing-construct
  explanation, no population off-target rate, no independent replication; plus the non-targeting reading
  (+0.43522, p 0.00666) explicitly marked "not a negative control and not a verdict cell".
- **C-50 amended** (touched by 305): T7 now "partly answered at the aggregate level by C-55/C-56 and
  still open for any individual construct", with the Round 5 addendum's restriction reaffirmed — C-55/C-56
  do **not** supply a population off-target rate.
- **C-53 amended** (touched by 301): 301's 4/13 vs 4/21 scan now recorded as **reproduced exactly by
  307's KA-2**, with both aggregate re-readings (−0.077 and +0.1949, opposite signs) and the statement
  that C-55/C-56 add an aggregate association, **not** a failing-construct explanation.
- **Preamble amended**: rows C-55…C-57 postdate the audit and its extension and carry no verdict.
- **Round 6 addendum** added after the Round 5 one: T7's new state, why 301 and F16 do not conflict, why
  C-56 is not a replication, and one bullet each for 308 and 306 (which get no claim row — 308 is
  engineering, 306 is reconnaissance).
- **Summary table**: C-56 added to "confirmed by its own frozen rule" *with the caveat in the cell*; a
  **new status line** "preregistered, aborted by its own known-answer rows" for C-55; C-57 added to
  "explicitly not claimed".

### 1.3 `experiments/294-writeup/REFERENCES.md`

- **[HORLBECK16]** added to §A: eLife 5:e19760 (2016), doi 10.7554/eLife.19760, CC BY 4.0, with what
  Supplementary files 3 and 2 supply, both cdn.elifesciences.org URLs, both sha256 prefixes and byte
  counts, **and the Excel-escaping defect (`="SYMBOL"`, 0.107 → 0.597) named in the entry itself**.
- **§E** gains three rows: the hCRISPRi-v2 composition + CAGE TSS annotation, Replogle Table S1 as the
  2025 library's protospacer source, and the Ensembl GRCh37 REST promoter fetch.
- **§F** gains six artefact rows (301, 305 + its coordinator review, 307, 309, 308, 306).
- The two stale pointers corrected: FINDINGS "(v1.9) F1–F14" → **v3.0, F1–F16** (with a note that every
  F15/F16-era number in the manuscript was read from the experiment's own `RESULTS.md`, not from the
  summary), and the ledger "claims C1–C4" → **C1–C9**.
- **Provenance discipline preserved**: a Round 6 note states that nothing was fetched from the web for
  this revision and that [HORLBECK16] is assembled entirely from on-disk records — doi/locus/licence/
  digests from `E301/work/sources301.md`, author list and title from the cached PMC page in the same
  directory. The entry's "accessed" column says so too.

### 1.4 This directory

- `work/verify310.log` — the number-by-number provenance log (sections A–H).
- `work/manuscript-round6.html` — 1,612 KB, built by
  `uv run python experiments/294-writeup/build_html.py experiments/310-writeup-round6/work/manuscript-round6.html`.
  26 TOC entries (was 25 — the new §3(g)), 8/8 figure captions, `3g-what-a-construct-does-besides-its-intended-gene`
  and `3h-what-bounds-zero-shot-prediction-in-this-challenge` both present, C-55/C-56/C-57 rendered.

---

## 2. Every number verified — the summary; the log is `work/verify310.log`

| block | source read | provenance entries |
|---|---|---|
| A | `experiments/307-offtarget-enrichment/RESULTS.md` | 41 |
| B | `experiments/309-c9-calibration-repair/RESULTS.md` | 48 |
| C | `experiments/305-t7-expression-test/RESULTS.md` + its `COORDINATOR-REVIEW-2026-09-07.md`; `experiments/301-metadata-acquisitions/RESULTS.md` | 6 |
| D | `experiments/308-kd-weighted-donors/RESULTS.md` | 16 |
| E | `experiments/306-recon-oct/RESULTS.md` | 19 |
| F | `experiments/301-metadata-acquisitions/work/sources301.md` and `work/pmc9380471.html` | 8 |
| G | `docs/FINDINGS-2026-09.md` and `docs/CONFIRMATION-LEDGER.md` — **structure only** (the F1–F16 range and the C1–C9 range), no measured value | 2 |
| | **total** | **140 entries, of which 138 carry measured values** |

Some entries carry more than one number (a by-level row, a three-fold triple), so "140 entries" is a
count of provenance lines, not of digits; every distinct value written into the manuscript is on one of
them. Section H lists what was read and deliberately **not** written.

Every line of the log carries the path, the line number and a verbatim anchor from that line.
`docs/FINDINGS-2026-09.md` v3.0 was read **for framing only**; no number in the manuscript was taken
from it, and the two counts taken from `docs/` are structural (how many findings, how many claims).

### 2.1 Discrepancies found while verifying, and how each was resolved

1. **The mandate said "of 770 entries"; the source says 769.** `306/RESULTS.md:18` records
   `work/lb_live_306.json` with **769 entries** and `:59` repeats "Of 769 entries, 292 write a
   `description`". **769 is what the manuscript says.**
2. **The mandate said rivals "score 0.136–0.150"; two ranges exist in 306.** Its verdict line (`:10`)
   says 0.136–0.175, which mixes the source-naming entries with the rank-42 emission-recipe entry.
   The four measured rows (`:74, :78, :80, :82`) are 0.14008, 0.14200, 0.13085, 0.14954. **The
   manuscript quotes the four measured values as 0.13085–0.14954 and reports rank 42 / 0.17543
   separately**, which is both narrower and exactly what the file supports.
3. **Two different floors exist for 305's design.** 307 (`:132`) characterises it as a Holm floor of
   11/158 ≈ 0.0696 on the declared family of eleven; 305 itself (`:134–136`) says the realised family
   was 8 testable pairs, giving an attainable minimum of **0.0506 > 0.05** at a rank-p floor of 0.0063.
   **Both are now stated together in §3(g)** so the two figures elsewhere in the paper cannot read as a
   contradiction.
4. **The exact binomial p for 26/34 appears twice in 309.** The abstract (`:53`) rounds to 0.0015; the
   primary table (`:195`) gives **0.00147**. The table's value is used.
5. **307's own write-up misstates its split sizes** (21/7/6 in prose against the 16/10/8 its code read),
   found by C9's audit (`309:124–130`). C8's numbers are unaffected. **Disclosed in §3(g)** in one
   clause, since a reader comparing the two write-ups would otherwise hit it.
6. **The parent library has 209,078 sgRNAs; the rescan used 209,070 protospacers.** The reference entry
   carries the file's 209,078 (`sources301.md:37`); the body carries the analysis's 209,070
   (`307:117`, `309:93`). Both are correct for what they describe and are not conflated.

---

## 3. Judgement calls, stated

- **Where F16 went.** As a new **§3(g)**, before the competition-facing section, so the results run
  biology-then-benchmark; the old §3(g) became §3(h) and all four cross-references (§1.4, §4.4, C-48's
  row, and the internal "§3(a)–(f)") were updated. The alternative — appending it as §3(h) after the
  benchmark section — would have put a biological finding after the paper's closing summary of ceilings.
- **No provenance verdict was invented for F16.** §1.4 already commits the paper to saying so when a
  claim postdates the audit, and that is what it now does. The prior art that *does* exist ([SEED24],
  [HARTMAN26]) is named in the same sentence as the claim, not buried.
- **No new figure.** A figure for F16 would require new computation, which this mandate does not
  authorise. The figure gallery is unchanged at F-1…F-8 and `build_html.py` needed no change.
- **308 got no claims-ledger row.** It is a competition-track engineering closure, and `CLAIMS.md`'s
  own convention (C-48) keeps scorer engineering out of the claims table. It is in the Round 6 addendum
  and in §3(h) in full, mechanism included.
- **306 got no claims-ledger row either.** It is reconnaissance about the world, not a claim about
  perturbation biology. It is in §3(h) where the paper discusses the field.
- **The non-targeting result was included, framed as 309 frames it** — a within-group consistency check
  on constructs where intended-target biology cannot be the explanation, explicitly **not** a negative
  control and **not** a verdict cell. Omitting it would have hidden the reading a reviewer would ask for;
  presenting it as a control would have been wrong.
- **Bounded language throughout.** Every F16 sentence that could be read as mechanism carries its limit
  in the same sentence or the next: "does not establish off-target binding" appears in the abstract, the
  plain-language summary, §3(g)'s heading paragraph, §3(g)'s limits paragraph, §4.4 and §4.5.

---

## 4. One thing this experiment could not fix

**The built page's masthead still reads "Round 5".** `build_html.py:86–88` hard-codes the eyebrow and
the two `.meta` strings — `Sources: FINDINGS v2.9 · confirmation ledger · audits 287, 297, 304` and
`Round 5 reviewed: … republished as Draft v8 (13:48Z, 2026-09-07)` — and that file is outside this
mandate's write scope, so the rebuilt page's header understates the body it now wraps. The body itself
is correct: the Round 6 revision note is the first thing after the title. **The coordinator should
update those three strings in `build_html.py` before republishing** (suggested: `FINDINGS v3.0 ·
confirmation ledger · audits 287, 297, 304` and `Round 6: F16 — C8 aborted, C9 PASS; 308 closes donor
quality weighting; 306 finds no acquirable data`).

Two smaller pre-existing items, noted and **not** changed because they are outside the brief:
`PAPER.md` §1.2 has an unclosed parenthesis ("(§3(a), §3(b) — three quantities"), §2's PLOLO bullet has
another ("(§3(a)."), and the TOC renders "Figures" twice because `PAPER.md` carries a `## Figures`
heading and `build_html.py` also injects one.

---

## 5. Considered and left out as too weakly supported

| candidate | source | why not |
|---|---|---|
| 307's P6 two-half agreement (Pearson +0.718, same sign 27/34) as corroboration | `307:196–204` | It is a **cell-sampling** split, not a biological replicate; 307 says so itself. Quoting it beside the primary would read as replication, which is precisely what the paper must not imply. |
| 309's Sec-2 wrong-recipient permutation p = 2.0e-4 as an independent check | `309:217, 269–273` | Under the per-construct centring repair it is the **same object** as the repaired no-signal row. The manuscript instead says the check was *lost* in the repair. |
| A "one truth in five is an off-target" style rate | — | Explicitly refused: 83 candidates over 34 constructs and 17 targets cannot support a population rate (C-57), and the Round 5 addendum already blocks the 2025→2026 extrapolation. |
| 306's ASD cortical-neuron preprint as a coming data source | `306:42, :173–179` | ANNOUNCED with **no accession verified**; 306's own rule is "do not download on the preprint alone". Not used as a claim about the field. |
| 306's IGC preprint as support for the incumbent's axis | `306:45, :197–201` | A reading of another group's abstract with no data behind it. Corroboration, not evidence; the manuscript does not lean on it. |
| 305's post-hoc Poisson-binomial rescue | `305:5–7` | Rejected by its own coordinator review as not an accepted confirmatory test. §4.5 says so. |
| "The variance cost of averaging two rows instead of three exceeds the bias benefit of removing a failed construct" | `docs/READ-FIRST.md` §5 | A coordinator gloss on 308, **not present in `308/RESULTS.md`**. Replaced with 308's own sentence: "the construct term is a ceiling on the truths; it is not a usable selector among donors" (`308:419–420`). |

---

## 6. Files

| file | what |
|---|---|
| `work/verify310.log` | every number written, with source path, line and a verbatim anchor (sections A–H) |
| `work/manuscript-round6.html` | the rebuilt page, 1,612 KB, 26 TOC entries, 8/8 captions |
| `RESULTS.md` | this file |

Edited outside this directory, under explicit authorisation and nowhere else:
`experiments/294-writeup/PAPER.md`, `experiments/294-writeup/CLAIMS.md`,
`experiments/294-writeup/REFERENCES.md`.

**Nothing was published.** The artifact at the URL recorded in `docs/FINDINGS-2026-09.md` is untouched;
republication is the coordinator's call, and §4 above lists what to fix in `build_html.py` first.
