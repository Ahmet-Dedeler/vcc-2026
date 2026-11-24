# 533 — Science lead: bring the manuscript up to date as an audited scientific record (worker mandate)

Written 2026-09-29 13:00Z by the coordinator, at Ahmet's request for a second long-running Claude Code (Opus 5.5, max
effort). You are the programme's science lead for this run. Read `AGENTS.md` first (above all §4 scoring, §11 method rules,
§14: nothing leaves this machine), then `docs/FINDINGS.md` §0–§2, then the manuscript in `experiments/294-writeup/`
(`PAPER.md`, `CLAIMS.md`, `REFERENCES.md`, `build_html.py`, `make_figs.py`) and the audit that last corrected it,
`experiments/503-cloud-claims-audit/`.

## Why this matters
Ahmet funds this project to win the Virtual Cell Challenge (final D/E/F entry, Nov 5) AND to produce scientifically useful
results, not only a rank. The manuscript "Three reproducibilities: what a single-cell perturbation reference measures, which
responses are conserved across cell contexts, and what bounds zero-shot transfer" (Draft v14, private artifact version 31,
last revised 2026-09-27; claims ledger C-01 … C-73) is that deliverable. Two days of results have not reached it, several of
them negative results that close whole families of ideas, and at least one (tonight's 09-30 read) answers the paper's central
open question in §3(j): does a same-lab, same-chemistry perturbation programme carry across cell lines?

## Who else is working right now (context; do not interfere with any of it)
- **The coordinator** (an interactive Claude Code session, Opus 5.5). It decides and arms leaderboard slots, reads board
  results, applies and commits patches, and integrates workers' outputs. It reads your `RESULTS.md`; it cannot message you
  while you run (you are headless), so write anything it must act on into `RESULTS.md` §2.
- **Worker 528** (Claude Code, Opus 5.5, `experiments/528-callprop-block/`): the DE call set's composition (a call-propensity
  block from per-gene dropout). It runs real-emitter jobs on the H1 harness (~7.5 GB each).
- **Cloud agents** (Cursor, Grok 4.7 extra-high; they work on GitHub branches `cloud/530-…`, `cloud/531-…`, `cloud/532-…`,
  with no data on their VMs): 530 implements the October pass-through bundle as patch files; 531 builds final-phase tooling
  and a tabletop run of the Oct 22 → Nov 5 runbook; 532 scans literature and data released since 09-15. Their results reach
  master later as `experiments/530…532` directories; if one lands while you work and is relevant, you may cite it.
- **Automated chains** (detached bash): leaderboard slot waiters (`arm0930.sh` submits at 00:01Z each day), the day-rule
  chains (`chain1001b.sh` … `chain1006.sh`), `recipe_dr524.sh` (the 10-05 October dress rehearsal), `fix526.sh` (a
  known-answer build running now), and `commit_after_chain.sh` committers. They detect builders by PROCESS NAMES: never start
  a long-lived command whose arguments contain their script names or production script names (AGENTS §6.2), never touch
  their files, never run builds.
- Finished in the last two days (your raw material, all on master): 513–527 and the board reads listed below.

## Your territory
- Work directory: `experiments/533-science-lead/`. At the start, copy the manuscript's sources
  (`PAPER.md`, `CLAIMS.md`, `REFERENCES.md`, `build_html.py`, `make_figs.py`, `figs/`) from `experiments/294-writeup/` into
  `experiments/533-science-lead/paper/` and edit ONLY that copy (in your copy of `build_html.py`, point its `W` at your
  `paper/` directory). The coordinator reviews your copy and moves it into `294-writeup` and republishes the private artifact.
- Write nowhere else: no git, no builds, no `vcc`, no `cell-eval2`, nothing under `data/`, `submissions/`, `docs/`,
  `scripts/`, or another experiment's directory.

## Objective
1. **Inventory** (RESULTS.md §1). Every result since the last revision (commit `f31ad98a`, 2026-09-27): the board reads
   453dr (the October path VALIDATED), 446bj2z (the zs rider transfers, pure nmae), 457dr (the CD4 half transfers on the
   October path, k ≈ 1.43), 446bd2 (DLD-1 as a second union residual: flat), and the 09-30 read (item 4); experiments 513
   (reach anatomy: sign-limited), 514, 515, 516 (NEGATIVE: target-specific public sign sources), 517 (the lab/Flex programme,
   the meanz_ne estimator, the programme-sourced uncalled-gene hook, Flex technical components, dropout predicts calls not
   signs), 518 (the bulk mse identity and the magnitude pot), 519/526 (engineering reviews: at most a Methods note), 520
   (NEGATIVE: D/E/F lines and same-line data), 521 (n_real predictability; per-target sizing closed by payoff; the global
   call-set size is regime-sensitive), 522 (October's bulk key; the October-like tails are real), 523 (the K pair, pending
   the board), 525 (NEGATIVE: public competitor DE designs), 527 (the hook arm built; pending the board), and 540–546 if they
   post-date the last revision. One line each: the finding, where it lives, its status (board-validated / H1-only / analytic /
   negative / pending).
2. **Verify before you write** (RESULTS.md §2; the 503 standard). Trace every number you will state to its primary artifact
   (a RESULTS table → the json / csv / log that produced it; a board number → `data/board-entries/<id>.json`). Recompute the
   cheap arithmetic yourself (deltas; scaled members s = (u − b)/(r − b); the overall as the mean of six) and show the
   recomputation. List every discrepancy with file:line and the correct value, for the coordinator. Re-plot saved tables only;
   compute no new statistic on raw matrices.
3. **Write** (in your `paper/` copy). Integrate each verified result where it belongs (§3(h)–(k), a new subsection only if
   none fits), add claims C-74 onward with the ledger's columns (evidence file, unit · n, interval, provenance NEW / EXTENDED
   / KNOWN / CONTRADICTED, confirmation status, what would falsify it), and update the abstract, plain-language summary and
   Discussion only where the evidence changes them. Negative results get their mechanism and their exact scope ("a negative
   row closes the construction it names, never a field"). §4.4 (explicitly not claimed) and §4.5 (contradictions) must
   absorb 516, 521 and 525, and anything the board contradicted (e.g. H1 prices that did or did not transfer). Add a revision
   note at the top (Draft v15, dated).
4. **The 09-30 read.** 458d (Arc's H1 programme as every context's DE block prior, vs 446b) and 458g (+ the design-strength
   gate, vs 458d) submit from 00:01Z 09-30; the coordinator reads them with `experiments/468-prior-gate/read468.py` and writes
   BOARD-RESULT files. Do NOT wait for it: finish the draft first. If the captures (`experiments/459-h1-axis-def/work/submit0930s{1,2}.log`
   "CAPTURED own entry …" lines and `data/board-entries/<id>.json`) or BOARD-RESULT files already exist when you are done,
   integrate them; otherwise leave explicit `[[PENDING 09-30]]` markers at every place the result belongs, and say in
   RESULTS.md exactly what to insert for each outcome (positive / partial / null / harm, using 458's bins).
5. **Figures.** Add a figure only if it carries a result better than a table (candidates: the October path's A/B/C
   decomposition 446b → 453dr → 457dr; the hook's H1 price against attenuation from 517's saved tables). Keep every figure
   reproducible from saved tables by `make_figs.py`.
6. **Open problems** (a short Discussion subsection): what would move the frontier, grounded only in the record. The
   rivals' nmae implies +0.06 … +0.08 gate sign accuracy over ours (434f / 480); which public information was shown
   insufficient (516, 525, 174, 184, 450h), and what kind of data would supply it.
7. **Build.** `uv run python experiments/533-science-lead/paper/build_html.py experiments/533-science-lead/paper.html` (with
   your copy's `W` repointed). The paper is PRIVATE until after Nov 5: never post, upload, email or share it anywhere.

## Rules
- Memory: builds run on this machine now and later. Keep every process ≤ 4 GB RSS; check `vm_stat` free+inactive ≥ 33 GiB
  before anything above 2 GB. No emitter runs.
- Every sentence traceable to a file; no number without its source; state uncertainty and scope plainly; British or American
  spelling as the paper already uses.
- Deliverables: `experiments/533-science-lead/paper/` (the updated manuscript, claims, references, figures, build scripts),
  `experiments/533-science-lead/paper.html`, and `RESULTS.md` with §1 inventory, §2 verification and discrepancies, §3 what
  changed (sections, claims C-74+, figures), §4 pending items (with the 09-30 insertion text per outcome). Finish with a final
  message summarising: the new draft number, claims added, discrepancies found, and pending markers.
- Take the time the work needs (hours); the first complete draft should exist by 2026-09-30 12:00Z.
