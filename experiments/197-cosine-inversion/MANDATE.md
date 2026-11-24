# 197 — Invert the board: how much better is the leader's DIRECTION, in cosine, and what would it take to match it?

You are a research subagent. Work ONLY in `experiments/197-cosine-inversion/`.
No `vcc`, no `cell-eval2`, no submission build, no writes outside your directory.
Peak RSS ≤ 8 GB. `uv run python`. This is mostly algebra on public board data plus our own frozen
artifacts; it should not need heavy compute.

## READ FIRST
`AGENTS.md` §4, §4.1, §4.1b, §8, §8.1, §11.1 rules 1–13. Then:
- `experiments/096-mse-floor/` **in full** — the `mse` decomposition and its four board-validated
  known-answer gates. This is your instrument.
- `experiments/187-joint-amplitude/RESULTS.md` — the ρ-cliff identity, and the fact that our global
  amplitude is already at the joint optimum, which is what makes the inversion well-posed.
- `docs/BOARD-LEDGER.md` — every scored arm of ours, all six members, raw and scaled.
- `experiments/159-transfer-ceiling/RESULTS.md` (κ² ≈ 0.34, the 1/k law),
  `184-tahoe-transfer/RESULTS.md` (line adaptation closed), `scripts/analyse_submission.py`
  (an existing, partial attempt at exactly this inversion — read it and say what it got right).
- `data/leaderboard/2026-09-03T0030Z.json` — 612 scored teams, raw and scaled for all six members.

## The question

Every strategic decision now rests on an unmeasured quantity: **how much better is the leading
direction than ours?** `AGENTS.md` §8 records an estimate of *"at least 1.28×"*, flagged as a late
correction that reversed a more flattering reading. That number is doing a lot of work — it is the
basis for "calibration buys rank 3–5, first place needs the direction" — and it has never been
re-derived on the current board.

**It is now invertible.** `expr_mse_unbiased_capped_norm` is, per 096,

```
mse − 1 = ( s²A + Kc − min(claim, s²B) − 2 s X ) / ΣD
```

and 187 established that our global amplitude sits at the joint optimum, so `s` is pinned and the
branch is known. For our own rows every term is measured. For other teams only `mse` is observed —
but `pds`, `nmae`, `fid`, `jac` and `reach` are observed too, and several of them constrain the
same underlying direction. The inversion is therefore over-determined and should be done as a fit
with stated assumptions, not as a single algebraic step.

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
Reproduce 096's four gates: `claim(0.25, φ0)` **0.11860**, `claim(1.00, φ0)` **0.19216**,
`claim(0.25, φ4)` **0.16478**, `claim(a15 as shipped)` **0.24581**; and a15's raw `mse` / `score_mse`
**0.88228 / 0.11112**. Then reproduce the out-of-sample record that makes the decomposition
credible: it predicted a10's raw `mse` at **1.3567** (actual 1.3524) and a11's at **0.9021**
(actual 0.9100) *before* either was built.

### Phase 1 — invert our own rows first
Using `docs/BOARD-LEDGER.md`, solve for our implied transfer cosine at each of our board points
(a8, a11, a12, a14, a15, a23, abc-knee, and `rc_hybrid` when its row lands). **The trajectory is
the check**: arms that changed only the emitter must invert to the *same* cosine, and arms that
changed the operator or the profile must move it. If the inversion assigns different cosines to
two arms with a digest-identical pseudobulk, the inversion is wrong and you must say so.

### Phase 2 — invert the field
Apply it to every board team with raw `mse < 1.2` (the shrunk cohort, 49 teams), and separately
bound the unshrunk ones. State the assumptions each inversion needs — the largest are that `ΣD`
and `Kc` are properties of the *reference*, hence shared across teams, and that `s` is near each
team's own optimum. **Test the second assumption rather than asserting it**: a team far off its
own amplitude optimum will invert to an implausibly low cosine, and 187's ρ-cliff gives you the
shape to check against.

Cross-check every inverted cosine against two independent observables that also read the direction:
- `pds_cosine` (a retrieval statistic — monotone in the direction but not a cosine);
- `nmae`, through 187 §5's exact identity `nmae = 1 − 2 q_w W_hi − (2 q_lo − 1) K_lo + K_hi`, which
  gives the |t|-weighted **sign accuracy** `q_w`.
**Three routes to one latent quantity. Report the disagreement, do not average it away.**

### Phase 3 — the target
Convert the answer into what it demands of us:
- What cosine (and what `q_w`) does the leader hold, and what is the ratio to ours?
- 159 says we sit at **31.1% / 28.5%** of the infinite-donor context-independent ceiling with 3.1
  effective donors, approached as `1/k`. **How many effective donors would the leader's cosine
  correspond to under that law** — and is that number reachable at all, given 166 found no further
  panel-covering line and 184 closed line adaptation on 50 lines?
- If the leader's cosine is *not* reachable by donor count, what else can produce it? Say what the
  arithmetic permits, and be explicit about what it rules out.

### Phase 4 — falsifiers
1. **A team whose `mse` is clamped carries no information** — exclude them explicitly rather than
   letting them enter a fit at the clamp value.
2. **Our own byte-identical resubmission** (`5QfX14axLAG8FF2rSDJe`) must invert to the identical
   cosine. Free known-answer row.
3. **The `1.28×` claim in `AGENTS.md` §8 must be adjudicated** — confirmed, corrected, or shown to
   be unidentifiable from board data. Whichever it is, say it plainly; that sentence currently
   governs the project's allocation.

## Deliverable
`RESULTS.md`: the known-answer gates, our own trajectory of implied cosines with the
digest-identical consistency check, the field inversion with its assumptions tested, the three-route
cross-check, the effective-donor translation, and an adjudication of §8's 1.28×. **A clean "the
inversion is not identifiable from board data" is a full result** — it would retire a number the
project has been leaning on.

## Territory
192 owns the call set. 193 owns the operator's eigenvalue map. 194 owns the donor tables. 195 owns
fitting the transform. 196 owns emitted depth. You own **the inversion and the target-setting**.
Touch no production code and build nothing.
