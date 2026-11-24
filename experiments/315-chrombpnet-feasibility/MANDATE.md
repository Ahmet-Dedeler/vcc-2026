# 315 — Feasibility: can PER-CELL-TYPE promoter chromatin features say anything about the held line that basal expression could not? (ChromBPNet / ENCODE4; the axis 289 closed for expression)

You are a research subagent. Work in `experiments/315-chrombpnet-feasibility/`. **Edit nothing outside it.**
**NEVER run `vcc submit`; never build a submission; nothing under `submissions/`; no `cell-eval2`.** You MAY use web
search/fetch and you MAY download into `experiments/315-chrombpnet-feasibility/work/` up to **8 GB total** — but
**nothing into `data/`**, and you must justify every download against a question below before making it.
`uv run python`. Peak RSS under **8 GB** (nothing else heavy is running; keep headroom). **Time box: 4 hours.**
A well-supported NO-GO is a valid and expected result and must be delivered as clearly as a GO.

## Why
**Established here, do not re-derive:** F5 (268) — on Broad-16 a reproducible line-specific response component
exists, is high-rank, and is **not predictable from basal EXPRESSION** by a linear model at 11 training lines.
289 then tested the source paper's own state features (p53 status, HPV-E6, lineage, proliferation programmes, arm
CN) and got **0/25 cells transcriptome-wide**; its conclusion was that a stage 2 needs **more LINES, not more
features**. 306 (today) established there is **no new perturbation dataset to acquire anywhere**. 314 (today) found
the new board leader's public label names an external promoter model, that PromoterAI itself is gated and
**cell-type-blind by construction** (its released score means over the track axis), and that **ChromBPNet is the
better instrument for the same idea: MIT, downloadable, 1,512 ENCODE4 models over 408 biosamples, base-resolution
accessibility PER CELL TYPE**. Borzoi (Apache-2.0), scooby (MIT) and AlphaGenome (click-through since 2026-01-28)
are the alternates.

**The question this experiment exists to answer is NOT "can we run ChromBPNet".** It is: **does a per-cell-type
promoter-accessibility feature carry information about the held line that the basal expression 289 tried did
not?** If it cannot, say so and stop — that is the deliverable.

## Gate 1 — identity (do this FIRST; if it fails, most of the rest is unanswerable)
Can the held contexts be tied to a named cell line, or at least to a nearest ENCODE/CCLE biosample, **from data we
already hold**? The challenge's held contexts are A/B/C (`data/challenge/`), and the six 2026 lines' identities may
not be disclosed. Use the contexts' own **control** cells and `data/external/ccle_broad_2025` / CCLE expression to
rank candidate lines by correlation, and state the margin between the best and second-best match. **Report the
identification as a ranked list with a confidence statement, never as a fact.** If A/B/C cannot be matched better
than chance, say so plainly — a cell-type-specific feature is then unusable for the competition track and the rest
of this mandate is a scientific-track question only.

## Gate 2 — coverage
For whichever biosamples Gate 1 nominates (and for the donor lines K562, HCT116, HEK293T, plus Broad-16's 12–16
lines, which are the surfaces where the claim could actually be TESTED): does ChromBPNet publish a model? Give the
ENCODE accessions, the assay (DNase vs ATAC), the file sizes, and the licence you observed. **Do not download 1,512
models.** If a handful of the relevant biosamples are covered, download at most a few and say why those.

## Gate 3 — the only test that matters
On **Broad-16** (`data/external/broad16line`, 82 Cas9 knockouts × 12 working lines, the object 268/289 used), ask
whether a per-cell-type promoter feature predicts the **line-specific residual** where 289's state features scored
0/25. Use 289's frozen residual object and its exact LOLO protocol so the comparison is like-for-like — read
`experiments/289-broad16-state-features/RESULTS.md` and reuse its splits. Controls, all required: the wrong
recipient's features, shuffled line labels, and 289's own feature set as the baseline to beat. Depth-condition per
programme rule 9. **If ChromBPNet models do not exist for enough Broad-16 lines to run this, say so and report the
count — that is the answer, and it is worth knowing.**

**A cheaper substitute if Gate 2 is thin:** ENCODE's *measured* DNase/ATAC signal for those biosamples is public
and needs no model at all. Measured accessibility at the 300 panel targets' promoters would test the same
hypothesis with less machinery. Prefer it if it is available for more lines.

## Rules
Programme rules 1–9 bind: name the experimental unit (the LINE is the unit that limits every interval here);
strong simple baselines (289's features, basal expression, magnitude); a failed feasibility test closes its
construction at that scope, not the field. Carry a known-answer row for anything you reproduce from 268/289.
**Nothing competition-facing is selected, tuned or built.**

## Deliverable
`RESULTS.md` with: Gate 1's ranked identification and its confidence; Gate 2's coverage table with accessions,
sizes and licence; Gate 3's result with its controls, or an explicit statement of why it could not run; and one
line — **315 (per-cell-type promoter accessibility for the held line): GO / NO-GO (identity …; coverage …/12
Broad-16 lines; residual prediction … vs 289's 0/25)** — plus what a stage 2 would cost if GO. Then stop.
