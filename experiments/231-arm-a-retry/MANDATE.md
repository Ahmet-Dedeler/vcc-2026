# 231 — 214's Arm A, re-measured against a gate that does not double-charge. The largest unharvested number we have.

You are a research subagent. Work in `experiments/231-arm-a-retry/`. **No `vcc`, no submission
build, no production change. No `cell-eval2` while another agent runs one — check first.** Peak RSS
≤ 10 GB. `uv run python`. Five other agents are live (221, 226, 228, 229, 230) — **229 is editing
`scripts/make_submission_pm_cd4.py` right now**, so if you read it, take a byte copy into your own
`work/` and pin its sha. Write only inside your own directory.

## The lead

**214's Arm A — drop PerturbMatch — prices ALL SIX members at +0.0088 native**, +0.0054 after 214's
own measured **×1.533** h1-test deflator, **+0.0046** if `reach` transfers at 164's **0.4277**.

**It was vetoed by a gate defect, and 219 §2.5 says so in terms.** 214 was veto'd by comparing rule
13's **two-member** composite (−0.000830) against the **+0.0089** floor. But:

> *"Rule 13's floors are, in `AGENTS.md`'s own words, 'the realised RMS of the **unpriced** channels
> for each arm class.' An arm that has priced all six has **no unpriced channels**. Comparing a
> six-member estimate against the RMS of the channels it just priced **double-charges the same
> uncertainty** — and comparing a *two*-member composite against it while a six-member estimate
> exists throws away the measurement. **This is a defect in the gate, not in 214**, which reported
> both numbers and chose the conservative one."*

219 ranked it **row 5: +0.0046…+0.0088, effort MED**, needing exactly two things — *"its `reach`
transfer factor, which is the whole arm; and a `pds` re-read on `cd4f034`"* — and was then
**explicitly forbidden** from measuring it (*"any re-measurement of a Class-A arm — out of
mandate"*). **Nobody went back. You are the follow-up.**

For scale: **`h1union` scored +0.0069114 and took us to rank 3.** The gap to first is **+0.01248**.

## READ FIRST
- `experiments/214-cheap-positives/RESULTS.md` — **Arm A in full**, its native price, its ×1.533
  deflator and how it was measured, and its own two reported numbers.
- `experiments/219-reprice-closures/RESULTS.md` **§2.5 and §3 row 5** — the gate defect and the
  re-price.
- `experiments/220-ordering-keys-retry/RESULTS.md` — **the channel bound you must reconcile with.**
- `AGENTS.md` §8.1 (the two-channel split), §11.1 (rule 13 and the digest test), §4.1c.

## 🔴 The reconciliation that decides whether this is real

**220 bounded the ENTIRE ordering channel from the outside**: board raw `reach` ours **0.24983**,
field best (`luka`/PertMind) **0.2908**, so *"matching the best team in the world on `reach` is
worth 0.0410 × 1.12927 / 6 = **+0.0077 overall**."*

**219 says every point of Arm A is in `reach`.** So Arm A at +0.0046…+0.0088 claims **60–114% of the
entire world-best-`reach` prize** from one component removal. **That is either the best arm on our
board or an artefact, and the two are distinguishable.**

**Do this reconciliation FIRST, before any re-measurement**, and report it as Phase 0:
1. Take our **current** board `reach` (the shipped `h1union` row, raw and scaled — note it has moved
   since 220 wrote its bound) and compute what raw `reach` Arm A's claim implies.
2. **If that implied raw `reach` exceeds the field best, the arm is refuted on arithmetic alone** —
   report that and stop. It is a complete result and it costs one afternoon rather than a slot.
3. If it does not, say precisely where it would place us and proceed.

## Then, the two things 219 named
1. **The `reach` transfer factor — this is the whole arm.** 214 used 164's **0.4277**. Re-derive it
   on today's instruments and say what it is, with an error bar. **A number this arm rests on
   entirely deserves its own known-answer row**, not an inherited constant.
2. **A `pds` re-read on `cd4f034`.** 219 §3 row 2 flags that 214's `+0.002069` `pds` was reproduced
   against `cd4f034` rather than `rc_callset`. Establish which baseline Arm A's numbers are against
   and re-read on the **currently shipped** stack (`h1union`, board 0.2118310), not on a superseded
   one — 208's rider rule.

## Falsifiers, all mandatory
1. **Price all six members, on the shipped baseline, in the challenge frame.** Never carry a member
   at zero — §8.4's whole lesson, and the reason 214's own two-member number was misleading.
2. **Rule 13's digest test**, and report the digest. Dropping PerturbMatch may change the emission
   class; **230 is concurrently testing whether that toll is a step or scales — read its `RESULTS.md`
   if it has landed, and if not, report your result BOTH ways** (toll charged and not).
3. **Leave-top-3-out retention on `reach` and cosine separately.**
4. **What does dropping PerturbMatch COST?** It is a retrieval component we ship deliberately.
   Report what it was buying — on `pds` and cosine, not only on `reach` — because an arm that trades
   direction for ordering is the trade `AGENTS.md` §4.1c explicitly warns against (MT Zeng holds the
   board's best direction and ranks 8th for exactly this reason).

## Reporting
**Do not convert PLOLO to board with a scalar** — the factor is **not identified** (1.384 on
`cd4f034`, 1.000 on `h1union`, matched `cos_agg`). **Only invert a board row whose `score_mse` is
strictly positive** — 15 of the 55 teams above 0.15 are clamped and 197's inversion returns cosines
for them that were never measured. Both traps cost this project real errors today.

## Deliverable
`RESULTS.md`: the Phase-0 reconciliation against 220's bound, the re-derived `reach` transfer factor
with its error bar, the six-member price on the shipped baseline, the digest, the falsifiers, and a
one-line verdict. **"Refuted on arithmetic against 220's bound" is a complete and valuable result.**
**Propose no slot; build nothing.**
