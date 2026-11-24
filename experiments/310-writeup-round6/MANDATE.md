# 310 — Fold Round 6 into the manuscript: F16 (the seed-match repression signal, C8 aborted / C9 PASS), 308's closure, and 306's reconnaissance. Then rebuild the HTML. Do NOT publish.

You are a research subagent. Work in `experiments/310-writeup-round6/` **and** you are authorised to edit
`experiments/294-writeup/PAPER.md`, `experiments/294-writeup/CLAIMS.md` and `experiments/294-writeup/REFERENCES.md`
— those three files only, plus your own directory. **Edit nothing else.** No `vcc`, no `cell-eval2`, no build of a
submission, nothing under `submissions/` or `data/`, no download, no publish (you have no publishing tool; the
coordinator republishes). `uv run python`. **Time box: 2 hours.**

## Why
The draft is at v8 (Codex's Round 5 revision). Round 6 has landed and the manuscript does not carry it:
- **F16** (`docs/FINDINGS-2026-09.md` v3.0) — the new finding. Sources: `experiments/307-offtarget-enrichment/RESULTS.md`
  (claim C8, `UNRESOLVABLE-AT-THIS-RESOLUTION`) and `experiments/309-c9-calibration-repair/RESULTS.md` (claim C9,
  **PASS**, S = +0.71293 [+0.32419, +1.07351], p = 1.99998e-05, 9/9 known-answer rows).
- **308** (`experiments/308-kd-weighted-donors/RESULTS.md`) — the competition-facing closure of donor quality
  weighting, and the mechanistic reason (the higher-|kd| donor is not the one closer to the held truth).
- **306** (`experiments/306-recon-oct/RESULTS.md`) — the reconnaissance negative, which is itself a finding about
  the field: no acquirable data exists, and no entry above 0.20 publishes a method.

## What to write
1. **§ T7 must be rewritten.** It currently records T7 as unresolved (305). The honest new state: 301's
   failing-vs-working negative **stands** (307's KA-2 reproduced its selection exactly: failing 4/13, working 4/21,
   NTC 4/31), 305 produced **no accepted discovery** (0/11 pass Holm; its design could not attain its own FWER),
   and F16 now supplies a measured aggregate association with the **construct** as the unit. **301's negative and
   F16 are not in conflict** — different questions — and the paper must say why in one sentence.
2. **A new results subsection for F16**, in the paper's existing voice. It must carry: the unit (construct, n = 34);
   the joint null's construction and its attainable floor (1.0e-5 vs 305's 0.0696); both C8 and C9 with the reason
   there are two (C8 aborted itself on 2 of 7 known-answer rows; C9 repaired exactly those two and re-tested at the
   same α); **the fact that the bias repair moved the estimate DOWN 9.5 % and was declared to do so in advance**;
   the partner-construct specificity test as the key design element; and the **two-tests-one-dataset** multiplicity
   stated plainly, with "this is not an independent replication".
3. **The limits, in the body and not only in a limitations paragraph:** it does NOT establish dCas9-KRAB off-target
   binding; expression is observational; no construct has a biological replicate; all three splits share one
   control pool; the partner test narrows the cause to something carried by the protospacer but cannot separate
   seed-mediated binding from a spacer-correlated construct artefact; with the construct as the unit, more cells do
   not narrow the interval, more constructs would. Name the decisive next experiments and their cost (a
   sequence-mutant control, ~$3–6 k; an occupancy assay, ~$5–10 k).
4. **308 into the competition-track narrative** as a closure with a mechanism, not a null result.
5. **306 where the paper discusses the field**: no acquirable data since the freeze; of 770 entries none above
   0.220 describes a method; every rival who does describe one names only sources already on our disk and scores
   0.136–0.150. Treat rival descriptions as unverified self-reports and say so.
6. **`CLAIMS.md`**: add rows for C8 and C9 with their provenance verdicts, and update any row 305 or 301 touched.
7. **`REFERENCES.md`**: add Horlbeck et al. (the hCRISPRi-v2 library) if not present, and note in the text the
   **Excel-escaping defect in its supplementary gene labels (`="SYMBOL"`)** that inflated a per-spacer hit rate
   from 0.107 to 0.597 and was caught by a known-answer row — a small but genuine contribution.

## Rules
- **Every number you write must be read from the RESULTS.md it came from.** Do not restate a number from memory or
  from another summary. Keep a log `work/verify310.log` recording each number and the file and line it came from.
- Do not overstate. Where the finding is bounded, the sentence must be bounded. Prefer the weaker true sentence.
- Do not touch `docs/`, the ledger, or any experiment directory but your own.

## Deliverable
The three edited files; `work/verify310.log`; a rebuilt page at
`experiments/310-writeup-round6/work/manuscript-round6.html` via
`uv run python experiments/294-writeup/build_html.py experiments/310-writeup-round6/work/manuscript-round6.html`
(read that script first; it may take arguments); and `RESULTS.md` listing every change made, every number verified,
and anything you judged too weakly supported to include. **Do not publish anything.** Then stop.
