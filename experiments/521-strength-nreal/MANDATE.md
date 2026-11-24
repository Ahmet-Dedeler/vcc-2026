# 521 — Per-target DE call-set sizing, reopened by a better n_real predictor (worker mandate)

Written 2026-09-28 14:30Z by the coordinator. You are a research worker (Claude Opus 5.5, max effort). Read `AGENTS.md`,
then `docs/FINDINGS.md` §0–§2, then `experiments/192-predicted-callset/RESULTS.md` IN FULL (it is the prior art and its
machinery is yours to reuse), `experiments/385-per-target-callset-size/RESULTS.md`, `experiments/513-reach-walk-anatomy/
RESULTS.md` and `experiments/518-bulk-mse-gap/RESULTS.md`. Your directory: `experiments/521-strength-nreal/` (write ONLY
there; `nreal521.py` in it is the coordinator's first measurement — keep it).

## Why this is open again
192 showed that sizing each perturbation's DE call set to its own n_real is worth **+0.036 overall (oracle, m ≈ 2 ×
n_real, fid + jac)**, and closed it "by its input, not by its payoff": its n_real predictor (features of our design H)
had LOTO Pearson +0.24 (6 % of variance, 11.8× multiplicative error), so the buildable policy lost −0.097 against the
one-sided penalty (under-calling a strong target collapses fid's yield linearly). Today (`nreal521.py`) 518's donor
response-strength key — a noise-normalised count of strongly responding genes in each donor's own same-target row
(`experiments/518-bulk-mse-gap/h1skey518.py: strength_for`) — predicts H1's log n_real at **Pearson +0.60 (K562-GWPS
alone) / +0.55 (mean rank of K562, Orion HCT116, Orion HEK293T)**, AUC 0.75 for negligible (n_real < 20) vs not; the
weakest 20 % by strength are 82 % negligible (base rate 41 %). The key needs only donor rows, so it can be fitted for any
panel, including D/E/F on Oct 22.

## Objective
Price, honestly and with the asymmetric risk in full, whether a per-target call-set policy driven by this predictor
beats the shipped global call set (192's `top 400 by |H| ∪ top 2,200 by control CPM`, `rc_callset`), on all four DE
members, at the board's n_real regime.

## Required, in order
1. **Known answers.** Reproduce 192's sizing table exactly with 192's own code (ORACLE m 2.0 → +0.0361 overall;
   PREDICTED m 2.0 → −0.0972) and 513's four H1 DE members (walk513.py).
2. **The predictor, done properly.** Out-of-sample (leave-one-target-out on H1; and on the cancer folds where a
   same-target truth exists — Orion HCT116 / HEK293T and K562 as held-out lines with the others as strength sources;
   480's fold gates), with the full predictive distribution (not just a point): residual sd in log space, calibration.
   Try the strength variants 518 built (per-donor, robust, CD4, VIPerturb, KOLF2) and simple combinations; report r,
   AUC(negligible), and the multiplicative error. Never use H1 data (h1_2025, h1_2025_full) as a strength source when
   H1 is the truth.
3. **Policies.** The Bayes action under the predictive distribution (choose n_pred per target to maximise expected
   fid + jac + reach + nmae given the one-sided penalty), and simpler robust forms: shrink only the bottom q % by
   strength to a floor K_min (never below it), leave everything else at the global set; or a two-level set. The call set
   CONTENT when shrunk must be deliverable by the emitter (which genes are dropped: block genes by CPM rank? the |H| head
   kept?) — read 192 §5.3 ("the protected design") and 395's q80 protection first.
4. **All four DE members.** 192 priced fid + jac; 513 shows our reach lives in small pools (n_real < 20: reach 0.64 on
   H1), exactly the targets a shrink touches — price reach through the walk (called genes first, then p order) with the
   dropped genes moved out of the called set, and nmae on the scorer's gate (the delivered lfc of dropped genes changes if
   the emitter stops shifting them). pds and mse must be pinned by the ballast construction; say whether that holds.
5. **The board's regime.** H1's n_real is 2–3× A/B/C's (192 §4.2: A/B/C median 31–45). Price at f = 1.0 AND at the
   A/B/C-like f 0.35–0.5 with 192's PowerSurface (proportional scaling, never a uniform n_real target).
6. **Target bootstrap** CIs on every headline, and the fraction of the gain carried by the top 3 targets.

## Deliverables
- `RESULTS.md`: known answers, predictor table, policy table (per member and overall, per f, with CIs), the verdict.
- If a policy clears **+0.003 overall at f 0.35–0.5 with a CI excluding 0 AND is not worse at f 1.0**: a build spec —
  the emitter/builder change it needs (as a patch file in your dir, NOT applied; `src/vcc2026/rc_emit_k.py` and
  `scripts/make_submission_pm_cd4.py` are production files you may read but not edit), the strength artifact for the
  A/B/C panel on the challenge gene axis (`data/challenge/pert_counts.csv` order), the expected config diff vs 446bj2,
  and predicted board deltas per member with ranges and the H1→board calibration used (`docs/H1-BOARD-CALIBRATION.md`:
  call-set levers retrodicted at 0.72× on fid/jac and inverted on reach in 192 §3.2 — carry that into your prediction).

## Rules
- Memory: a production build queue runs until ~02:00Z 2026-09-29. Until then keep EVERY process ≤ 4 GB RSS and check
  `vm_stat` free+inactive ≥ 33 GiB before any step expected above 2 GB; after that ≤ 10 GB per process, one heavy step
  at a time. Report peak RSS (`/usr/bin/time -l`). If 192's frontier needs more, subsample targets and say so.
- No builds, no `vcc`, no `cell-eval2` CLI, no git, no edits outside your directory. Never start long-lived shell
  commands whose arguments name production scripts (the build chains watch process lists).
- Known-answer row in every screen; negative results are results. Finish with RESULTS.md and a final message
  summarising the numbers.
