# 553 — R7b's instrument: the real-emitter DE harness on a sealed Calu-3 surface (worker mandate)

Written 2026-09-30 ~22:00Z by the session "Explore contribution areas" (a peer of the coordinator; it owns reserve R7). You
are a headless research worker (Claude Opus 5.5, max effort). Read `AGENTS.md` first (§6.2 process traps, §6.4–6.5 harness
rules, §9 memory, §11 method rules, §14), then `docs/FINDINGS.md` §0–§2, then **`experiments/552-r7-reserve/RESERVE.md` in
full (its §3, §5 and §6 are binding on you)** and the ledger's "Reserve R7" section in `docs/CONFIRMATION-LEDGER.md`.

## Why this matters
The prize is the last D/E/F entry on Nov 5 (three unseen cancer lines, a dark leaderboard). The October recipe's DE switches
(528's call-propensity block, 496's CLEAN4 uncalled-gene hook, 523's call-set size K) have been priced off the board only on
H1 (h1-test), a pluripotent line that has chosen many settings (tuned-on). R7b is a sealed, never-used cancer-line surface
(Calu-3, Sunshine et al. 2023, CRISPRi; its cells are SARS-CoV-2-infected, a stated caveat) on which the final DE recipe
gets ONE confirmation read after the switches settle (not before the 10-07 board reads). Reading rule (RESERVE §5.2): R7b
can only VETO a switch on a preregistered HARM beyond its own re-roll band; null or positive reads never change a switch.
Your job is to build the instrument that read will use, and prove it on h1-test, WITHOUT reading the reserve.

## What exists
- **The H1 instrument** (the thing to generalise): 323 (`experiments/323-current-family-validation/`: `run323.py` builds the
  production-builder effects for h1-test's 85 targets; `evaluate323.py` builds `work/reference.npz` from the truth — note
  its hard checks: 18,400 controls, 400 truth cells per target — and scores emitted blocks with cell-eval2's own
  functions) → 400 (`run400.py`, the per-target emit/score chain) → 410b's key-scaled bulk (`--bulk ps`) → 416 (carriers) →
  426 (MDE) → 429 (zs) → 435 (coarse) → 496 (`run496.py`, the hook) → 528 (`run528.py`, the call-propensity key/block) →
  537 (`experiments/537-arm-pipeline/run537.py`: the MODE grammar `cpm|lcap[-kept][-n<cpm>][+c4h][...]`, its known
  answers in `ka537.py` / `drive537.sh`, prices in its RESULTS §3). Read-only for you.
- **The surface** (built by 552's `build_r7b.py`; if `experiments/552-r7-reserve/work/r7b/build_r7b.json` does not exist yet
  when you start, wait for it, polling every 10 min, and work on the code meanwhile):
  - `work/r7b/bundle/context_R.h5ad`: the 3,536 non-targeting control cells on the official 18,533-gene axis (obs
    `target_gene` = non-targeting, `ntc_id`). **Allowed**: any control-only statistic (RESERVE §6).
  - `work/r7b/bundle/pert_counts.csv`: the ~110 eligible targets; `build_r7b.json`: cells per target (metadata).
  - `work/r7b/truth_r7b.h5ad`: **SEALED. Do not open it, load it, stat its matrix, or compute anything from it.** You may
    check its sha256 against `build_r7b.json`. Only the access script (below) opens it, and only when the owner runs it.

## Objective (deliver in `experiments/553-r7b-harness/`, updating RESULTS.md as you go)
1. **A surface-parametric harness** that prices 537's arms on a SURFACE = (control bundle, target list, gene axis, truth
   path) instead of H1's frozen inputs. Generalise, do not fork logic: same emitter module (`src/vcc2026/rc_emit_k.py`, sha
   0f56209d asserted), same CallSet semantics, same scoring functions. It must handle any number of controls (3,536 here)
   and any number of truth cells per target (100–300 here); the emitted prediction stays 400 cells per target, as a
   submission is.
2. **KA-H1 (abort-on-failure, before anything else counts):** pointed at h1-test's inputs, the generalised harness
   reproduces 537's recorded rows bit for bit: at least `cpm` base seed 0 (pre-zs and zs) == 528's `h1_cpm`; `lcap-kept`
   == 528's `h1_lcap-kept`; `cpm+c4h` == 517's `h1_c4-c4t0.05`; one `-n<cpm>` row from 537 §3.2. Use their stored stats
   files; report |d| per member.
3. **The Calu-3 prediction side, built from controls + public donors only:** per-target donor effects for the R7b targets by
   the SAME construction 323/410b used for h1-test's (the three production donors `k562_gwps`, `orion_hct116`,
   `orion_hek293t`; the whitening operator for this panel refitted the way the October day does, `scripts/build_union_whitener.py`
   with the panel's genes off the axis — check what 323 used for H1 and mirror it; record path + sha of every operator);
   the carriers' CLEAN4 q80 prior on the official axis; the call-propensity key from R7b's controls (528's formula); the
   emitter's control block from R7b's controls. **`sunshine_calu3` / `SunshineHein2023` must never enter any donor mix,
   bank, prior, registry or key.** Audit donor ↔ surface target overlap and write it down. Emitting the arms' predicted
   counts is allowed (it reads no truth).
4. **`access553.py` (written and tested, NOT run on R7b):** computes the reference statistics from a truth file the way
   `evaluate323.reference()` does (generalised to any control count and truth cells per target), scores each arm × seed,
   and builds the noise band the reading rule needs: (a) emission seeds (≥ 3 per arm) and (b) split-half re-draws of the
   REFERENCE (truth cells split in halves per target, each half scored as the truth; ≥ 5 splits, seeded). Output: per arm
   and member, the paired effect vs the incumbent (`cpm` base with zs), its seed sd and split-half sd, and the overall-
   equivalent (mean of six scaled members, scaled by the surface's own baseline/replicate anchors computed inside the access).
   **Test it end to end on h1-test as the stand-in truth** (h1-test is not a reserve) and show the output format. Never run
   it on R7b.
5. **Price the arms' EMISSIONS on R7b** only as far as needed to prove the prediction side runs (e.g. counts per arm and
   seed written, call-set sizes, K per target: emission-side diagnostics only; nothing compared to the truth).
6. RESULTS.md: what was generalised (file:line), KA-H1 |d| table, the R7b prediction-side audit (targets, donor coverage,
   operator shas, control-block stats), how to run the access (exact commands, time and memory), and anything that makes the
   surface unfit (say so plainly; a well-argued "R7b cannot be priced faithfully because …" is a valid result).

## Hard rules
- Write ONLY in `experiments/553-r7b-harness/`. No git, no `vcc`, no `cell-eval2` CLI, no builds of submission packages, no
  edits to `scripts/`, `src/`, `docs/`, `data/`, `submissions/` or any other experiment's directory (including 552: its
  files are inputs). Reading anything else is fine.
- **Sealing (RESERVE §3/§6):** nothing that involves a perturbed Calu-3 cell's expression, ever. If any code path you adapt
  would read the truth file or a pseudobulk of target cells, stub it for R7b and run it only on h1-test.
- Memory and concurrency: ≤ 8 GB RSS per process (a real-emitter H1 run is ~7.5 GB). Before any job above 2 GB require
  `vm_stat` free+inactive ≥ 24 GiB AND `uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P;
  print(len(P.other_builders()))'` prints 0 AND no bash process runs `run_arm_generic.sh`, `def_day.sh`, `build5NN*.sh` or
  `rebuild5NN*.sh` (match interpreter + script path, AGENTS §6.2). Operator fits (~13 GiB) need ≥ 27 GiB and the same checks.
  Heavy jobs strictly serial, under `/usr/bin/time -l`. **No heavy job 2026-10-04 05:00Z → 2026-10-05 18:00Z** (R6's C10b
  access and the coordinator's build window). Worker 537 (`experiments/537-arm-pipeline/`) also runs real-emitter jobs:
  never pause or signal it; alternate through the gate.
- Process-name traps: name your drivers `*553*`; never put a production or chain script's name in a long-lived argv.
- Known answers first (AGENTS rule 10); report both KA and new numbers with seeds; negative results written with their
  mechanism.
- Deadline: 2026-10-06 12:00Z for items 1–4 (the access happens after the 10-07 reads, on the owner's go). Keep RESULTS.md
  current; the owner reads it (you cannot receive messages; if blocked, write `BLOCKED.md` with the question and keep
  working on what you can).
