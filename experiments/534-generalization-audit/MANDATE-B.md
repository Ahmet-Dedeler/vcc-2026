# 534 Part B — worker mandate: the sealed mock final (bulk channel)

Written 2026-09-29 13:55Z by the session "Public data scope [2e83f3]", which owns experiment 534. You are a headless
Opus worker. Nobody can message you mid-run. If something blocks you, write `experiments/534-generalization-audit/work/partB/BLOCKED.md`
with the exact reason and stop.

## Objective

Ahmet asked whether our score climb solves the task or only fits A/B/C. Every October bulk switch was chosen on A/B/C's
board and then re-priced on three October-like panels (`real`, `efp`, `efr`). Those panels have informed dozens of
decisions since, so they are tuned-on.

Your job is to execute `experiments/534-generalization-audit/PREREGISTRATION-B.md` exactly:
1. draw four fresh Arc-like panels (reserve R6) that no one has ever used;
2. price the October bulk recipe's cumulative milestones M0 → M6, plus two controls, on them with 438's instrument;
3. read the result with the frozen rules in its §7.

The preregistration is binding. Do not change a seed, a design, a rule or the order. If something in it cannot be
executed as written, stop and write BLOCKED.md.

The number you produce is the honest projection of how much of the September bulk climb survives a new panel. The
incumbent is 446bj2z = 0.2728260 overall on the board. The per-step A/B/C and tuned-on values are in 438's RESULTS
(§0, §3.3–3.5) and 522's RESULTS; quote them from there, do not recompute them.

## Read first

- `AGENTS.md`: all of it, especially 6.2 (process-name traps), 6.5, 9 (memory gate), 11.1 and the new 11.1c.
- `docs/FINDINGS.md` §0–§2.
- `experiments/534-generalization-audit/PREREGISTRATION-B.md`.
- 438: `experiments/438-october-stack-repricing/RESULTS.md` (§0–§5), `lib438.py`, `full438.py`, `dayof438.py`,
  `price438.py`, `comp438.py`.
- 522: `experiments/522-october-bulk-key/RESULTS.md`, `dump522.py`, `arms522.py`, `price522.py`, `boot522.py`
  (the percentile-cap route and its bootstrap).
- 238: `experiments/238-essential-free-panel/bundle238.py` (the draw and bundle rules) and whatever 238 and 233 used to
  fit the per-fold held-out operators and extract the CD4 banks (`cd4_238.py`, `experiments/233-october-retention/`
  `buildfolds233.py`, `whiten233.py`, `panels233.py`). Also 438's §2.1 on what follows the panel.

## Steps

1. **Isolate the instrument.** Use a copy of `lib438.py` (and of any 522 / 233 / 238 module that writes files), or
   monkeypatch it, so that `WORK`, `bind()`'s work dir, `bankstats/` and every other cache point under
   `experiments/534-generalization-audit/work/partB/`. Never write into another experiment's directory. The coordinator
   warned that `lib438.bind()` writes `bctrl-{ctx}.npy` and `bankstats/{name}.npz` under 438's WORK.
2. **Known answers on `efp` (preregistration §5), before any draw.**
   - M0–M5 and C-a must reproduce 438's stored values, and M6 522's `pctPROD` row, to |d| ≤ 1e-12.
   - On failure, stop. Debug only on `efp`, never on R6.
   - Record the KA table in RESULTS-B.md.
3. **Draw R6** with a script named `draw534.py`, per preregistration §2.
   - At draw time, exclude every target in every `pert_counts.csv` under `experiments/` and `data/`, and each earlier R6
     draw.
   - Assert 300 unique targets, 0 essential, 0 overlap with any excluded panel, all on the gene axis.
   - Write `work/partB/draws.json` with each target list, pool size and sha256.
   - Build each synthetic bundle the way `bundle238.py` does (relabelled A/B/C controls, byte-identical to 199's).
4. **Per panel, in the order S1, S2, S3, S4:**
   - fit the three per-fold held-out operators and whatever else 438's route needs for that panel: CD4 extraction,
     H1 union masking, and the 410 / LODO keys re-fitted per fold. The panel's own targets are never in an operator or
     key fit;
   - assert each operator's axis for its own panel (preregistration §5);
   - price M0–M6, C-a and C-b, levelled at k = 1;
   - dump the per-target arrays the bootstrap needs.
5. **Aggregate.**
   - Pooled paired target bootstrap: 4,000 draws, seed 534.
   - The step verdicts, retention, tuned-on check and controls, exactly as in §7.
   - The calibrated board-unit projection per §6.
6. **Write** `experiments/534-generalization-audit/RESULTS-B.md`:
   - Lead with one sentence: how much of the M0 → M6 bulk gain survives on never-used panels (retention, with an
     interval) and which switches TRANSFER / are FRAGILE / FAIL.
   - Then the table: step × (A/B/C, real, efp, efr from the record; S1–S4 and pooled from you), the KA table, the
     controls, memory and timing.
   - Also write `work/partB/results534B.json`.
   - Also write a draft of the ledger outcome paragraph in `work/partB/ledger-C10a.md`. The owner appends it to
     `docs/CONFIRMATION-LEDGER.md`; you do not.

If time runs short, finish S1 and S2 fully (one panel per rule) before starting S3. A clean two-panel answer beats four
half-done ones. Say in RESULTS-B.md which panels were priced.

## Rules

- **Write scope.** Write only under `experiments/534-generalization-audit/`; your files go in `work/partB/` and
  `RESULTS-B.md`. Do not edit `PREREGISTRATION-B.md`, `MANDATE-B.md`, `RESULTS.md` (Part A, the owner's) or anything
  outside 534. No git, no `vcc`, no `cell-eval2`, no builder edits, nothing in `submissions/`.
- **Memory gate** (preregistration §9, binding):
  - start a heavy process (operator fit or scoring) only when `memory_pressure` ≥ 40 % free, `vm_stat` free+inactive
    ≥ 33 GiB, and no builder or chain runs;
  - check for running builders with `uv run python -c "import sys; sys.path.insert(0,'scripts'); import def_pipeline as d; print(d.other_builders())"`,
    or the equivalent `ps -axo comm,args` filter on `make_submission_pm_cd4.py` / `build_union_whitener.py`. Also
    check for bash processes running `run_arm_generic.sh`, `run_ka_*`, `def_day.sh`, `chain4NN`/`chain5NN*.sh`,
    `relevel_sp2.sh`, `run_dump4*`, `build527.sh`, `fix526.sh`;
  - if the gate is closed, wait by polling every 60 s (max 2 h per wait), then re-check;
  - strictly serial; run each heavy process under `/usr/bin/time -l` and report peak RSS;
  - operator fits at `--max-rss-gib 14` (~12–13 GiB `ru_maxrss`, as 238 measured) are AUTHORISED for this mandate,
    above the launcher's ~10 GB default. Say so in RESULTS-B.md.
- **Coordinator's windows:**
  - no heavy process starts after 2026-09-29 23:45Z or during 2026-09-30 00:30–06:00Z;
  - if unfinished at 23:45Z, write `work/partB/RESUME.md` (exactly what is done and what is next) and stop.
- **Heavy log.** Append every heavy process's start and end (UTC, command, peak RSS) to `work/partB/heavy.log`, so the
  owner can tell the coordinator when a whitener fit runs. The coordinator's build gates count any
  `build_union_whitener.py` process as a builder and wait for it; that is expected.
- **Names.** Name your scripts `*534*.py` / `*534*.sh`. Never name a script `chain5NN*`, `run_ka_*` or `def_day*`,
  and never put those strings in a long-lived process's argv (AGENTS 6.2: gates match on names).
- **Donor hygiene.** Never use h1_2025 on an h1 surface (none is in scope). Exclude each panel's own targets from any
  operator or key fit.
- **Verify, don't narrate.** Every number in RESULTS-B.md comes from a log or json in `work/partB/`. If a result
  surprises you, re-run its known answer before believing it.

## What others are doing (do not duplicate)

- The coordinator: slots, builds and the October recipe. Its next heavy window starts 2026-09-30 00:30Z.
- Worker 528: the DE call-set composition on H1. It runs intermittent emitter jobs of ~7.5 GB each; your gate covers
  them.
- Worker 533: the manuscript.
- The owner of 534 (Part A): an audit, from the written record, of board Δ vs held-out Δ for every step and switch.
  You do not need to do that part.
