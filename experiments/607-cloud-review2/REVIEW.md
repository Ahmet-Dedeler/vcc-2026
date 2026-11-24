# 607 — second adversarial review of the 10-06 chain

Reviewed commit: `7dbb1007fc0b182c487402711bb113a0b4d40c72` (`git rev-parse HEAD` at the start of this review).
Branch `cloud/607-review2`. No `vcc`. No challenge h5ad, donor matrices, or harness. Numpy only, plus the tracked npz files.

605 (`experiments/605-cloud-council-review/REVIEW.md`) read `2788ccd9`. Its 630 donor-gain blocker and its 612 tie-at-cut finding in `patch612.py` are fixed on this commit (69972aa1, 740f5bcf). They are not re-opened. Two 605 findings are still present in non-shipping harnesses; those are marked below.

Apply order checked: 641 → 598 → 612 → 612dd → 630 → 642 → 604 → 644, built by `test_skip644.chain_tree` (which is `test598.canon_tree` plus `compose606.PATCHES`). Final shas match 606's PINNED 13:35Z table. The live repo was not modified by that build.

## Findings

### MINOR — a pinned reference that shares a context's name is scored as a sister of itself and dropped

`experiments/604-identity-prior/def_identity_prior.py:111` and `:116–120`.

The pool and the sister similarities are dicts keyed by profile name. The context is inserted first; a later reference with the same name overwrites it. `sims[name]` then correlates that reference with itself and reads 1.0, and the drop deletes every reference with that name.

Executed on the composed `def_sameline.py`, with `ctrl_cpm` fed from `flex_refs604.npz` (no h5ad) and context B called `K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE` while the pin named `B` is still HELA:

- provenance for B: `sister_similarity["B"] = 1.0`, `sister_dropped = ["B"]`, references left `["A", "C", "h1"]`. HELA was removed. The call guard would have kept it (K562 ≠ HELA).
- provenance for A (bundle B is K562, pin B is HELA): `references = ["B", "C", "B", "h1"]`. Both profiles stay in the score, under one name. The similarity dict holds only the second.

This does not touch the checked-in A/B/C artifacts: there the context calls equal the pin calls, the second profile is removed by the call guard, and the reference lists are `B,C,h1` / `A,C,h1` / `A,B,h1` with nothing dropped. It does not touch D/E/F either: those context ids are not `A`/`B`/`C`, so they cannot collide with the pin names. It does touch any rehearsal that reuses ids A/B/C with different calls.

Fix: key the pool by an id the context cannot share (`ctx` vs `ref:{i}`), and drop that one entry.

### MINOR — `DEF_DAY_SAMELINE_SKIP=E, F` passes the bash check and then the day dies

`experiments/644-runbook-amendments/patch_skip644.py:77–78` (bash) and `:55` (Python).

Bash word-splits on the comma-to-space translation, so the tokens are `E` and `F` and both match a manifest `D,E,F`. The same string is passed through as `--skip-contexts`. Python splits on the comma and does not strip, so the second token is `" F"`.

Executed: bash printed `token=[E]`, `token=[F]`, `bash_ok=1`. The composed `def_sameline.py` then raised `SystemExit: --skip-contexts names [' B'], which are not contexts of the call ['A', 'B', 'C']` for `"A, B"`. Fail-closed: the build stops, it does not skip the wrong context. `E,F` with no space works.

Fix: `x.strip()` in the Python set, and the same strip in the bash loop.

### MINOR — `VCC_DE_UPRN=b_k100_q20r0p10` is accepted and moves nothing

`experiments/630-rank-neutral-up/src_rank_up.py:48–49` and the loop at `:140`. The builder regex in `patch_builder630.py:59–60` accepts the same string (`rounds` is only `\d+`).

`parse_spec` returns `rounds=0`. `range(0)` does not run. Executed on an all-deep block that the same function moves when `rounds≥1`: the block is unchanged and `n_moved=0`. A meta record of "0 genes moved" then looks like an empty selection. `quantum=0` is different: the promotion still happens (`counts=1` on that block), because the cap is `max(k_promo, quantum)` (`src_rank_up.py:200–202`).

Fix: reject `rounds < 1` in `parse_spec` and in the builder's range check. Reject `quantum < 1` the same way if a round is not allowed to move only the promotion.

### NOTE — still present: `run612.py` keeps a tie at the cut on the strong side

605's MAJOR on `patch612.py` is fixed (`patch612.py:79–82`, commit 740f5bcf). The pricing harness was not.

`experiments/612-key-gated-hook/run612.py:102–103` is still `weak = (~nz) | (a < cut)`.

Executed at q=0.5: key `[3,3,3,3]` → patch weak count 2 (`floor(0.5·4)`), run612 weak count 0, cut 3.0. Zeros stay weak and the fill is stable index order: `[0,0,5,5,5,5]` → weak on indices 0–3.

The shipped gate matches run612 wherever no nonzero key equals the cut (0 mismatches on the off-tie draws below). A price from `run612.py` and a build from the patch diverge exactly on targets that tie at the cut. The design-key matrices are not in this checkout, so the count of such targets is unchecked. What to count: on the keys `run612.py` priced, how many targets have some nonzero `|key|` equal to `np.quantile(|key|[nonzero], q)`.

Fix: call `_kg_mask_sign` from the harness, or copy the three tie lines.

### NOTE — still present: the 630 pricing harness fails open; v1 still clamps the gain

605's MAJOR guard was the only implementation it could see. The production entry point is now `apply()`, and that one fails safe (confirmed below). The harness main loop was not updated.

`experiments/630-rank-neutral-up/rnup630.py:449`: `if not flips or guard_rounds >= 5: break`, then the flipped block is scored and a WARNING is printed (`:466–467`). v1 `up_pass` still clamps the denominator at `rnup630.py:144` (`max(gain, 1e-300)`). v2 `up_pass2` in that file has the eligibility filter (`:287`).

Fix: on exhaustion, break only when `not flips`, and keep the pre-pass block when the round cap is hit, as `apply()` does.

### NOTE — the chosen carrier is the only one tried

`experiments/630-rank-neutral-up/src_rank_up.py:171–176`.

If the carrier nearest the control max is at least as deep as every ballast row, the gene is skipped. There is no second candidate. Executed: 30 carriers at library 10,000, one carrier at 500,000 sitting closest to the control max, ballast at 100,000. The block is unchanged, `n_no_pos_donor=1`. The same 30 shallow carriers with the deep carrier removed move the gene (`n_moved=1`, 1,845 counts). Nothing with a lower mean is accepted. This is a missed raise.

It cannot bind while every carrier library is below every ballast library. The owner's "85/85 H1 blocks, ballast ≥ 994,978, carriers ≤ 24,917" is not rechecked here: those blocks are not in git. What to check: per target, `min(ballast library) > max(carrier library)` on the H1 seed-0 dumps and on the A/B/C capacity blocks.

Fix, if wanted: on `not elig.any()`, try the next carrier in the candidate list before giving up.

### NOTE — `identity@` and `identity:` with an empty dose raise ValueError

`experiments/604-identity-prior/patch604.py:43`. `float("")` is not caught. Executed for both spellings. A negative dose is accepted (`identity@-1` → −1). The day then dies with a traceback rather than `SystemExit`.

Fix: parse in a `try`, require `lam >= 0`, and `raise SystemExit` with the source string.

### NOTE — the skip log says the context is unchanged while a programme source then writes one

`experiments/644-runbook-amendments/patch_skip644.py:65`.

Executed on the composed tree with `--programme-source identity:0.5 --skip-contexts B`. The log line is `B: … skipped (--skip-contexts); context unchanged (644)`. The artifact then has B's identity programme, `w=0`, and no rows. The arrays are the unmatched-source contract (confirmed below). The sentence in the log is false whenever a programme source is set.

Fix: log `unmatched (644)`.

### NOTE — score 0 copies the donor's side; P = 0.5 leaves the destination gene's CLEAN4 sign

`experiments/604-identity-prior/def_identity_prior.py:140–144`.

`side_dst == 0` keeps `Psrc`, so a covered gene whose score is exactly 0 inherits a donor side whenever that donor's P is not 0.5. `P == 0.5` does not update `psign`, which stays the destination gene's own CLEAN4 sign. Neither edge occurs on the three v3 artifacts: 0 genes with score exactly 0, 0 covered genes with P exactly 0.5, at λ 0, 0.5 and 1.

Fix: when `side_dst == 0`, write P = 0.5 so `psign` stays that gene's CLEAN4 sign, which is what the docstring already says for P = 0.5.

### NOTE — still true, already recorded by the owner: fixed totals do not pin mse

605's MAJOR. `mse` is a jackknife over cells. Equal per-gene totals and equal libraries leave it free to move. 630's RESULTS treats the H1 mse column as a measurement. Not re-measured here (no blocks, no scorer).

## Confirmed

### 630 — no accepted gene has a lower mean or a flipped sign

Module `experiments/630-rank-neutral-up/src_rank_up.py`, sha `d8c39c3f4e1a0fc17be3105e609fcea0bb34f6fdaff6df7400330d269d334671`. Pre-fix copy `work/src_rank_up_prefix605.py` still shows the 605 bug. `test_neggain630.py`: 11/11.

| case | result |
|---|---|
| N1, 8 carriers at 22k, 32 ballast at 18k, both `cross` values | pre-fix lfc +0.00144 → −0.28614, 12,736 counts. Fixed module: block unchanged, `n_no_pos_donor=1` |
| N2, half 18k / half 995k | only the deep rows give; lfc → +1.00146; libraries and column sums exact |
| N3, all ballast 995k, 396 carriers | fixed == pre-fix bit for bit; lfc → +1.00152 |
| equal libraries (16 ballast at the carrier library) | fixed moves nothing (`n_no_pos_donor=1`); pre-fix still moves 7,777 counts and the mean does not rise (gain 0) |
| 4 ballast tied with the carrier, 8 at 400k, `n_car=40`, both `cross` | the tied rows are bitwise unchanged; lfc +0.00144 → +1.00152; one carrier row gains the gene; every donor row is deeper than that carrier |
| shallow rows hold 800 of the gene and sort first (`argsort(-xb)`), deep rows hold 40, both `cross` | shallow rows bitwise unchanged; deep rows give (8 and 30 counts); mean rises; sign holds |
| carrier `log1p` set equal to the control max | the promoted row lands strictly above it (`log1p` 12.15 > cmax 4.52); lfc → +1.00149 |
| `quantum=3, rounds=4, rec_pool=6` | exactly one carrier row gains the gene (12 counts) |
| negative base lfc, deep donors | the move that would cross 0 is reverted (`n_mean_revert=1`); block unchanged |
| 80 random mixed-library blocks × both `cross`, with `quantum=5, rounds=3` | 160 moved blocks, 0 mean drops, 0 sign flips, 0 negative counts, libraries and column sums exact, one recipient row, `n_mean_revert=0` |
| 0 ballast rows, and `n_car=0` | no exception; nothing moved |

Eligibility is `cw[bal] < cw[a]` (`src_rank_up.py:171`), with `cw = 1/(n · library/1e6)`, so a donor is used only when its library is strictly larger. The mean check is `a_new < a_old or sign_new != sign(l_base)` (`:218`), and a hit restores the column and the sub-gate budget (`:221–225`). `ZS.fill` never takes from a non-positive cap (`rank_mag.py:175`).

`own_col`. `apply` clears `sel[own_col]` before the pass (`src_rank_up.py:266–267`). On a make_case block whose real Wilcoxon q is 4e-31, forcing that gene's q to 0.4 (so it sits in the uncalled-UP set the pass is defined to edit): `own_col=-1` selects it and moves it (mean 5198 → 5275, sign held, sums exact); `own_col=0` leaves `sel[0]` false and the block identical. The builder computes that index only for a target gene on `rc_scorer_control.keep` (`patch_builder630.py:85–87`).

Guard. `GUARD_ROUNDS = 5` (`src_rank_up.py:37`). With `de_stats`/`up_pass2` stubbed so every attempt flips: 6 `up_pass2` calls (the attempt at `rounds==5` still runs), then `X0.copy()`, `ur_guard_failed=1`, `ur_guard_rounds=5`, `ur_n_moved=0`, no raise. A stub that flips once and is then clean returns the modified block, `ur_guard_failed=0`, `ur_guard_rounds=1`. There is no path that returns a block whose BH calls still differ. Official N4 matches.

Knob. `parse_spec` and the builder regex agree on `b_k100`, `s_k050`, `b_k300`, `b_k100_q20r3p10` (accept), and on `b_k049`, `b_k301`, `x_k100`, `b_k100_q20`, `""` (refuse). Unset: the composed builder (sha `c2ee0e618876`) sets `DE_UPRN = None` and puts both the import and `_rank_up.apply` in `if DE_UPRN is not None`. The meta key is the ternary `{} if DE_UPRN is None`. No count is written when the variable stays None. A byte-identical unset submission was not rebuilt: that needs the challenge h5ad and a full emit. What to check before the first production build: an all-unset known-answer rebuild against the pre-630 builder, same recipe, same prediction hash.

`ballast_count` returns k in `[1, BALLAST_KMAX]` (`scripts/make_submission_pm_cd4.py:1023`), so the production call is not the K=0 case. The module still survives K=0.

"Exposure none structurally" overstates the filter. The filter is correct for any block; it is idle on a block only when every ballast library is above every carrier library. That inequality is an empirical claim about H1, not a property of the code.

### 612 / 612dd — the tie fill, the knob, the no-op

`_kg_mask_sign` (`patch612.py:79–82`): weak = key 0, or `|key| <` the linear quantile of the nonzero keys; if that count is short of `floor(q · n_nonzero)`, genes with `|key| == cut` fill the gap in index order, and only that many.

Executed: 16,160 draws (n = 1..8, q in {0.1, 0.25, 0.5, 0.75, 0.9}, constants and random multisets). Nonzero-weak count was never below `floor(q·n)`. Off a tie, the mask equals the run612 rule (0 mismatches). On a tie it reaches the floor. The largest excess over the floor, which is the old linear rule when the cut falls in a gap, was 1. All-zero key returns the hook sign unchanged (`:76–77`).

`_kg_apply` with the env var unset returns the same object (`patch612.py:88–89`). q of `0`, `1`, `1.5`, `nope`, `-0.2` raises `SystemExit`. The knob set with `DE_MDE_SIGN is None` raises `SystemExit`. The per-target call sits in the `de_mde_sign is None` else-branch (`patch612.py:96–98`), so an unset hook does not evaluate it. The meta key is absent when the knob is None (`:105–107`).

612dd on the composed `def_day.sh`: the quantile check is inside `[ -n "${DEF_DAY_MDE_SIGN_KEYGATE:-}" ]`, it requires `DEF_DAY_MDE_SIGN`, and that check precedes `export VCC_DE_MDE_SIGN_KEYGATE` (`patch612dd.py:45–54`).

### 604 — the v3 artifacts match the stated map

Recomputed with `identity_prior` and, separately, an assignment written from the claim (zc centred at 1/2, zi mean/sd over every gene where I is finite, side = sign(score), confidence = CLEAN4's `|P−1/2|` in order of `|score|`, ties by own `|P−1/2|` then index, mirror when the sides differ). Both match the npz arrays on programme and pdown, `equal_nan=True`.

Context CPM is the pinned `cpm_A/B/C` in `flex_refs604.npz`. `build604.py:23–24` documents those as the challenge-control CPMs, which this checkout cannot recompute (no h5ad). Given that identification, the artifacts are exactly this function.

| λ | sha256 | flips vs CLEAN4, A/B/C | confidence max abs diff vs CLEAN4 `|P−1/2|` |
|---|---|---|---|
| 0 | `f65d12fc3a7d8a7557e071ea6f86dfc1e89847cd0eea859b6ac14b26ef6fc6e6` | 0 / 0 / 0 | 0 |
| 0.5 | `6cd84dc258cb609d08d93ec4c6c6ec7c7134ff8e3093eccde5ac15bbd2ac618f` | 597 / 601 / 648 | 5.6e-17 |
| 1 | `3b1e22f037aff26b7a13e06328d011757423a51a8cf4acaef7d26cd206fd6a07` | 1,102 / 1,125 / 1,199 | 5.6e-17 |

These are the shas and the flip counts in `experiments/604-identity-prior/RESULTS.md`. λ 0 is CLEAN4 on covered P and on every gene's psign; the function raises if it is not (`def_identity_prior.py:145–147`). zc's mean on covered genes is +0.152, so the centre is 1/2 and not the mean (`:128`). Uncovered genes keep CLEAN4's psign and a NaN P. `identity_only_scored` is 0 (6,346 covered). References on this panel: A→`B,C,h1`, B→`A,C,h1`, C→`A,B,h1`, nothing dropped. zi's sd is taken over every finite I, including genes CLEAN4 does not cover; restricting it to covered genes would not reproduce the artifacts. The side is applied only on covered genes (`:132`).

`sign(P−1/2) == sign(score)` on every covered gene with score ≠ 0 and P ≠ 0.5, at all three λ (0 mismatches).

Parse, executed: bare `identity` → 0.5; `identity:0.5` and `identity@0.5` → 0.5; `identity@1` → 1; `identity@0` → 0. `:` and `@` are both one character, so `src[len("identity:"):]` is the dose for both (`patch604.py:43`). Hookbank split, the 530 rule still in the composed day: `hookbank:identity:0` → bank `identity`, τ 0, accepted; `hookbank:identity@1:0` and `hookbank:identity@0.5:0` accepted; `hookbank:identity:0.5` refused because 0.5 is read as τ and τ must be `< 0.5`. The `@` spelling is the one that carries a dose through that split.

Anchor: `patch604.py` replaces `if a.programme_source not in BANK_DIRS:` once. Matched contexts keep rows, w, and targets; unmatched contexts get w 0 and empty rows (`patch604.py:54–59`). The inserted comment still says "quantile-mapped". RESULTS line 39 says that wording was left so the post-644 `def_sameline.py` sha would not move. The map that runs is the v3 function, not that sentence.

### 644 — a skipped context takes the unmatched path, including the identity source

Composed tree, fake `load_bank` (one non-panel row, one panel row; no `delta.npy`), `ctrl_cpm` from `flex_refs604.npz`. Calls: A Jurkat (matched), B K562 (matched), C CAL33 (unmatched). All unanimous, min top-1 0.75.

- `identity:0.5`, `identity`, `identity:0`, `identity@0`, `identity@1` each match `identity_prior` at 0.5, 0.5, 0, 0, 1 on those calls.
- No skip: A and B keep a row and w 2; C has w 0, no rows, and C's identity programme.
- `--skip-contexts B`: B has w 0, 0 targets, 0 rows, 0 signs, `programme_ok`, mode contains `644`, and a `programme_source` record. A's rows and w are unchanged. Programmes and pdown of A, B, and C are identical to the no-skip run, so B is still in the identity reference pool. "Unmatched" here is the row contract, not "drop this context's CPM and call."
- `--skip-contexts A` with `identity@0.5`: A is w 0 with no rows; B's rows survive; programmes unchanged.

`test_skip644.py` covers this shape for `orion_hct116` and does not import the identity source. The identity case above is the one that was missing. It holds.

### The composed chain

Shas of the built tree, 12 hex, against 606 RESULTS "PINNED 2026-10-02 ~13:35Z". All twelve match.

| file | sha |
|---|---|
| `scripts/def_day.sh` | `3918e0537ee5` |
| `scripts/def_sameline.py` | `d964aa03f6e2` |
| `scripts/def_identity.py` | `7073406f6aad` |
| `scripts/def_day_helpers.py` | `fb715e3833f3` |
| `scripts/def_pipeline.py` | `76314af6b249` |
| `scripts/october_recipe.py` | `9dbcd9149020` |
| `scripts/make_submission_pm_cd4.py` | `c2ee0e618876` |
| `scripts/def_preflight.py` | `17f0836fd1e1` |
| `scripts/def_postbuild.py` | `557c2b522b56` |
| `experiments/524-october-recipe-dr/gate524.py` | `b16e60c08e48` |
| `scripts/def_identity_prior.py` | `a2d48da00f01` |
| `src/vcc2026/rank_up.py` | `d8c39c3f4e1a` |

`bash -n` / compile is what `compose606.py` already gates; this build applied every patch in order, including 644 last, and the file hashes are the pin. `compose606.py` itself deletes its tree. The checks here used a tree kept long enough to hash and to run `def_sameline.py`.

## Not checked, and what would settle it

- H1 "ballast library > every carrier library" on 85/85. Needs the seed-0 block dumps. The module is safe without that inequality; the inequality is only the claim that the new branch never binds in production.
- An unset `VCC_DE_UPRN` submission is byte-identical. Needs an all-unset rebuild against the pre-630 builder. The source-level guard is in place.
- How many priced targets tie at the 612 cut. Needs the design-key matrix `run612.py` used. The two gates differ on that set only.
- `test612.py` T4b (85 real h1-test keys). The harness npz is not in git; T4b skips. The synthetic tie rule was re-run here.
