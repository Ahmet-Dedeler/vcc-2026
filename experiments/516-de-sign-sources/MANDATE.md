# 516 — Target-specific DE sign sources: find one that beats production on held-out surfaces (worker mandate)

Written 2026-09-28 13:15Z by the coordinator. You are a research worker (Claude Opus 5.5, max effort). Read `AGENTS.md`,
then `docs/FINDINGS.md` §0–§2, then the experiments named below, before writing code. Your directory:
`experiments/516-de-sign-sources/` (write ONLY there).

## Objective
The DE channel (nmae, fid, reach, jac) is where we lose to the leaders, and it is **sign-limited** (480, 513): on H1 our
emitted sign is right on 0.561 of the reference's DE-pool genes (0.600 on genes we call, 0.537 on genes we do not); the
board calibration implies ~0.60–0.62 on A/B/C; the nmae leaders (0.83–0.85 raw vs our 0.897) imply ~0.65. Find a
**target-specific** per-(target, gene) sign source (plus a calibrated confidence) that beats our production signs by a
margin that would move nmae and reach on the board, on surfaces that transfer.

## Territory (disjoint from the other workers)
- YOURS: target-specific information — the same target's rows in other lines' screens, related targets' rows,
  learned combinations/stacking of donors, denoising (Wiener/low-rank), reliability weighting, conditioning on the
  recipient's own control expression and co-expression as FEATURES.
- NOT yours: target-independent programmes / lab or chemistry signatures (worker 517), the bulk channel pds/mse (518),
  any production code, builds or slots (the coordinator).

## Instruments and known answers (reproduce before reading anything new)
- `experiments/513-reach-walk-anatomy/walk513.py` reproduces the four H1 DE members of the production emitter exactly
  (fid 0.5089291277653674, jac 0.0746985342, nmae 0.9427040518, reach 0.2692569776) from
  `experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz` and
  `experiments/323-current-family-validation/work/reference.npz`. Reuse its loaders; your first log line must show
  these four numbers reproduced.
- `experiments/513-reach-walk-anatomy/purity513.py` (donor agreement vs sign accuracy) and
  `experiments/515-lowrank-de-signs/lowrank515.py` (Wiener: 0.572 → 0.581 on H1 pool genes) are today's baselines.
- `experiments/480-nmae-instrument/common480.py` (nmae on a fixed gate; `load_donor_rows` converts donor banks to
  log2 ratios) and 480's cancer-fold gates (`experiments/480-nmae-instrument/work/gate-*.npz`, `fold480.py`).
- Surfaces: **H1** (h1-test truth; NEVER use `data/donors/h1_2025` or `h1_2025_full` as a donor when H1 is the
  recipient — they cover h1-test's targets: oracle leakage) and **cancer folds** (Orion HCT116 / HEK293T, K562-GWPS held
  out, the others as donors; read 480 §4 for the Orion technical blob and the power caveats — relative direction only).
  A source must win on H1 AND at least one cancer fold to count.
- Price every candidate by (i) sign accuracy on the recipient's DE pool by pool band (< 20, 20–100, ≥ 100), with a
  target-bootstrap 95 % CI of the difference vs production; (ii) analytic H1 nmae with the candidate's sign at our
  delivered magnitudes (and at L1-optimal magnitude ∝ (2P − 1) if you produce calibrated P); (iii) H1 reach with the
  candidate's signs substituted in our order (walk513's machinery). Report raw deltas and their H1 → board calibration
  (`docs/H1-BOARD-CALIBRATION.md`: magnitude levers transfer ~1.2–2.3×; sign-prior reach is not priceable on H1).

## Prior art — read, do not repeat under a new name
434 (sign-confidence GBM, AUC 0.589; no learned magnitude beats the MDE), 477/480 (consensus rider: −0.003 nmae),
486 (a learned reliability gate does not transfer), 379 (similarity-weighted donors), 190 (network signs < |H|),
399/403 (neighbourhood tilts and denoised rows: failed on the board), 430 (unwhitened donor means), 513 (reach is
sign-limited; the block reorder loses), 515 (Wiener +0.009, hard low-rank hurts small pools). FINDINGS §2's closed list.

## Directions worth trying (yours to choose and extend)
1. Supervised stacking across donors, per gene or per gene class, with leave-one-line-out learned weights and shrinkage;
   include the CD4 atlas (3 conditions × 4 human donors: an internal replicate structure for reliability), RPE1, Jiang
   2025's lines and the Broad-16 where they cover targets (`data/external/…`, `data/donors/…`).
2. Denoise each donor row (Wiener with a per-donor noise model from its own cell counts) before combining.
3. Recipient-aware gating: a gene barely expressed in a donor carries no response there; weight donors per gene by the
   donor's detection power for that gene relative to the recipient's control expression.
4. Related-target borrowing: targets whose donor responses correlate (same complex/pathway, learned from the donor
   screens' own correlation structure) to impute or shrink a noisy row (403 failed as a bulk neighbourhood tilt; this
   is a per-gene SIGN source for the DE channel, which is a different use — price it as such).
5. Calibrated P(sign) with magnitude ∝ (2P − 1), which is nmae's L1 optimum under uncertainty.

## Deliverables
- `RESULTS.md`: the known answers, every source tried with its numbers (negatives with their mechanism), and a verdict.
- If a source beats production by ≥ +0.02 pool sign accuracy on H1 AND on a cancer fold (CIs excluding 0): write, in
  `work/`, the per-(target, gene) sign and confidence for the A/B/C panel's 300 targets on the challenge gene axis
  (`data/challenge/gene_names.csv` order; context-agnostic unless your method uses the recipient's controls), with a
  one-page build proposal (which emitter input it replaces, expected config diff, predicted board deltas with ranges).

## Rules
- Memory: a production build queue runs on this machine until ~02:00Z 2026-09-29. Until then keep EVERY process ≤ 4 GB
  RSS and check `vm_stat` free+inactive ≥ 33 GiB before any step you expect above 2 GB. After 02:00Z, ≤ 10 GB per
  process, still one heavy step at a time. Measure peak RSS with `/usr/bin/time -l` on heavy steps and report it.
- No builds, no `vcc`, no `cell-eval2`, no git, no edits outside your directory. Never start long-lived shell commands
  whose arguments name production scripts (the build chains watch process lists).
- Known-answer row in every screen; target-bootstrap CIs on every H1 claim; negative results are results.
- Finish with RESULTS.md and a final message summarising the numbers.
