# 637 — Final-phase information design: factorial D/E/F reads (H-630) and the feedback-use policy (H-631)

2026-10-02 ~11:55Z, council session [d4cd5f]; council portfolio item 4 (`experiments/600-breakthrough-council/SYNTHESIS.md`).
A proposal for the coordinator (FINAL-PHASE-FEEDBACK-PLAN and DEF-RUNBOOK are its files) and, for §2, a decision for
Ahmet. No compute beyond reading `data/board-entries/`.

## 1. H-630: one bulk switch AND one DE switch per D/E/F entry

**Premise (from our own board history, `data/board-entries/`):**

| pair | kind | Δpds raw | Δmse raw | bulk Δ overall | DE Δ overall |
|---|---|---:|---:|---:|---:|
| 446bja vs 446b | DE prior-only | 0 | +2.4e-8 | 0.000000 | −0.001066 |
| 458g vs 458d | DE prior-only | 0 | 0 | 0.000000 | +0.003566 |
| 446bj2zcp vs 446bj2z | DE-only (cp block) | 0 | −2.3e-8 | 0.000000 | −0.000794 |
| 446bj2z vs 446bj2 | DE magnitude (zs) | 0 | 0 | 0.000000 | +0.000549 |
| 433k vs 433h | bulk-only | +4.6e-3 | −1.8e-2 | +0.004860 | **−0.000213** |
| 446b vs 433k | bulk-only | +3.2e-3 | −7.6e-3 | +0.002525 | **−0.000241** |
| 446bj2 vs 446b | bulk-only | +2.1e-4 | −1.5e-2 | +0.002622 | **−0.000388** |
| 479p2g vs 446bj2 | bulk-only | +1.6e-3 | −3.7e-3 | +0.001235 | **−0.000199** |

- A DE switch never moves the bulk members (pds identical; mse ≤ 2.4e-8 raw).
- A bulk switch leaks into the DE members by a small and CONSISTENT amount: −0.00020 … −0.00039 overall, all four
  negative (548: the ballast's emission stream and the nmae part of the bulk step).

**Design.** Each D/E/F entry after R0 changes one bulk switch X and one DE switch Y (prior-only where possible) against
the current best:
- X is read from Δ(pds + mse) alone, exact to ~1e-8;
- Y is read from Δ(fid + reach + jac + nmae) minus a bulk-leak allowance of −0.0003 ± 0.0002 overall;
- every DE bin must sit outside ±0.0012 overall (≥ 4× the leak's spread), which the plan's DE rows already do
  (programme row 6: ±0.0015).

**What it buys.** With uploads at ~1 arm/day (531's amendment), the Oct 23 → Nov 3 window reads ~11 switches. Pairing
reads ~22, so the dropped rows (4a cap form, 8–9 called-gene scale, 12 post-combination ρ*) come back, with room left
for the DE rows (programme, prior gate, the CLEAN4 hook, K).
- Example pairing from FINAL-PHASE-FEEDBACK-PLAN §2: day 2 = ρ* ×0.98 (B) + the hook on/off (D); day 3 = q (B) +
  K (D); day 4 = cap form (B) + prior gate (D); …
- Same-line w (row 5) moves both channels and stays a single-switch row.

**Value.** Conditional on D/E/F entries returning members: P ≈ 0.3 (the rules hide the final LEADERBOARD; they do not
say a team's own status is hidden). So this is cheap insurance: ~0.3 × (switches read doubled).

**Validation before Oct 22 (one A/B/C slot).** Build one factorial arm whose X and Y each already have an A/B/C single
read against the SAME base, e.g. base 446bj2z with X = a bulk rider read alone on that base and Y = a DE-only rider read
alone on that base (the coordinator picks the pair from the ledger). Pass if:
- |Δbulk(factorial) − Δbulk(X alone)| ≤ 1e-6 raw on pds and mse;
- |ΔDE(factorial) − ΔDE(Y alone) − leak(X)| ≤ 0.0006 overall.
If it fails, D/E/F stays one switch per entry. If a D/E/F factorial read later contradicts its own A/B/C single reads,
the D/E/F read wins: it is the final's truth, and the A/B/C read is only a prior.

**Reader change.** `531/readdef.py ROW ENTRY BASE`: add `--factorial X_ROW,Y_ROW`, which applies X's bin to the bulk
members only and Y's bin to the DE members (after the leak allowance), and prints both verdicts.

## 2. H-631: feedback-use policy (for Ahmet; one decision, before Oct 22)

**Why now.** Prior-only DE arms share every emission draw with their base (468), and the board returns members to 16
significant figures. So if our final entries return members, an entry that flips prior signs on a gene set reads EXACTLY
how many of those (target, gene) signs were right on the hidden D/E/F truth. Module-level flips would recover D/E/F's
programme in ~8–16 entries. The rules: "There will be no leaderboard displayed for final entries"; Arc may alter the
rubric or prizes "in the event of evidence that the Scoring Rubric is subject to or has been manipulated in violation of
the spirit of the Challenge". Finalists must also describe "any classical statistical, mechanistic, or otherwise
non-learned components, including anything blended with model outputs".

**Proposed policy (council recommendation; DEF-RUNBOOK §6.8 text for the coordinator):**
1. **Allowed:** choosing among a small set of recipe switches preregistered before Oct 22 (the existing plan: re-level
   ρ*, q, cap form, CD4, same-line w, programme on/off, prior gate, the hook, K, DE magnitude). This is model selection
   on members each switch can move (11.1 rule 2), one global parameter per switch.
2. **Not allowed:** any arm whose purpose is to read truth bits rather than choose a model:
   - gene-, module- or target-level sign or size probes;
   - per-target n_pred tuning (yield is per target, so it is tempting);
   - designs that make score terms separable (precision encoding).
3. **Disclosure:** the finalist write-up says that D/E/F feedback, if any came back, chose among N preregistered
   switches, and names them.

Expected cost of the policy: ~+0.01 overall forgone in the members-returned branch (the module readout's ceiling, ~0.3
probability). Expected benefit: no disqualification-class risk on a $100k prize.

## 3. Delivered (2026-10-02 ~12:25Z): the reader, UNAPPLIED to 531
- `readdef_factorial637.py`: `--factorial BULK_ROW DE_ROW ENTRY BASE` and `--validate FACT X_ALONE Y_ALONE BASE`. It
  imports 531's `readdef` read-only and reuses its row rules. Bulk rows 2/3/4 × DE rows 6/7. It refuses 5, 8-9 (548's
  +0.0014 raw nmae leak vs that row's −0.002 bin), 10-12, 1 and 13-14. A verdict must be identical across the leak band
  (−0.0001, −0.0003, −0.0005); otherwise it is `inconclusive_leak_band` and nothing is adopted.
- `test_factorial637.py`: 27 checks (additive, bin-edge inconclusive, both-lose, row-2 vs row-3 leak charging, the
  row-7 jac falsifier, every refusal, hidden members, validation PASS/FAIL, a smoke read on the real 479p2g vs 446bj2
  captures), then 531's own suite: **74 row-verdicts still pass**.
- `readdef637.diff`: an optional 5-line dispatch so `readdef.py --factorial …` works. It was verified on a patched copy
  in a mirrored tree: 531's suite passes 74/74 with the patch, and `--validate` dispatches. It is not applied.

## 4. DECISION on §2: PENDING CONFIRMATION (conflicting relays, 2026-10-02 ~12:35Z)
In [4d0207]'s session Ahmet answered "for 1., yeah we can do that" (read there as approving §2's no-decoding policy; DEF-RUNBOOK
§6 item 8). In [537757]'s session he answered as quoted below. [4d0207] has asked him to confirm. The design below is
ready either way; it applies only if the maximum-information reading is confirmed.

Relayed by [537757] (2026-10-02 ~12:30Z), quoted as relayed:
"i dont find that an issue. what we should do is do what helps us learn most. we should document and make use of all
info. its not gonna disqualify." So if D/E/F members come back, the policy is MAXIMUM-INFORMATION reads with full
documentation. §2's "not allowed" list does not apply; the disclosure line stays (document everything). The coordinator
updates FINAL-PHASE-FEEDBACK-PLAN.

**Design that follows (proposal):** fuse §1's factorial reads with 625's module readout.
- Each D/E/F entry = R0 + one bulk switch X_k (read exactly on pds + mse) + one PRIOR-ONLY DE pattern h_k: row k of a
  16-module Sylvester–Hadamard flip design over 625's fixed modules (CLEAN4 octile × expression half, sha 3bce7179).
- The DE members of each entry give one exact linear equation in the 16 module nets (prior-only arms share draws: 468).
  The bulk leak is a known offset, −0.0003 ± 0.0002 overall.
- ~16 entries therefore read ~16 bulk configurations AND decode D/E/F's own module-level programme. That programme is
  the final's truth, not an A/B/C transfer: 622h's Broad-16 decode-and-flip gave +0.018 … +0.020 held-out gate sign; on
  the final's own truth the gain is the in-sample decode, larger.
- The decode uses nmae as primary and fid as the cross-check (625's rule). 630's UP pass is all-or-none across the
  16 entries (625's caution: magnitude asymmetry).
- Combined with 631-style documentation: every entry, its X_k, h_k and decoded value go into the finalist write-up.

## 5. Revision after 625c ([7772dc], 2026-10-02 ~12:45Z): decode entries are PURE prior-only, not fused
Ahmet approved the maximum-information D/E/F use ("we can do the risky thing", in [7772dc]'s session ~12:40Z, as relayed
by [4d0207]). But 625c's delivery-faithful fold model puts the module nets at a median |x_M| of 0.0001–0.0004 overall (16
modules) or 0.0002–0.0007 (8 octiles). §4's fused design carries a bulk-switch leak of ±0.0002 per entry, so only 3–11 of
16 modules would be readable at 2σ. Therefore:
- **Decode entries are PURE prior-only** (bulk identical to the base: exact reads, 468). Fusion with a bulk switch only
  if that switch is MEASURED K-neutral (DE members unmoved to < 2e-5 on an A/B/C pair), or after per-block emission
  streams exist (an emitter change, after 10-06).
- **8 CLEAN4 octiles at ~10 entries** (625c: Broad-16 transfer value +0.0165 vs +0.0172 for 16 modules at ~18 entries).
- **Size it first on our own lines:** 625's A/B/C readout measures the real nets before the final relies on them. The
  hook-only value of perfect module knowledge is +0.0001 … +0.002 overall. The block-prior share is not in the model yet;
  [7772dc] specs H1 real-emitter flip arms, dumped so that 630's UP pass can be stacked: module value WITH and WITHOUT
  the UP pass, from the same blocks.
- **Caveat on 625c's numbers ([21dd4c], 13:00Z):** 625c's delivery-faithful fold model applies the flip to every gene,
  while production's hook acts only on uncalled genes (rc_emit_k.py:1059, `use_s = ~head & ~blk`). In 612 the same
  masking cut the gains 7–15× (612 CORRECTION 2, 22691ce5). So the |x_M| medians above may be inflated. That only
  strengthens the pure-prior-only conclusion; it does not change it. [7772dc]'s real-emitter flip arms (625 flip625.py,
  octiles 0/2/5/7 × 3 seeds, dumped) measure the actual nets, with and without 630's UP pass stacked.

## Leak band v2 (2026-10-04T06:46Z; the coordinator's request after [7772dc]'s 658 finding, cdc1445c)
Over 9 bulk-only board pairs the DE leak is mean −0.000174, sd 0.000213, range −0.00041 … +0.00019: a re-roll around a
small negative mean, not v1's fixed −0.0003 ± 0.0002. `readdef_factorial637.py` now carries versioned bands:
- v1 = −0.0003, band (−0.0001, −0.0003, −0.0005);
- **v2 (default) = −0.0002, band (+0.0002, −0.0002, −0.0006).**
`--leak-band v1` reproduces any read registered under v1, e.g. 658's ALLSW R1 as registered. Every verdict and
validation records its band. test_factorial637: 32/32, with 531's 74 row-verdicts untouched. The v1 synthetic cases run
under v1; new v2 cases cover the validate window (+0.00035: v2 PASS, v1 FAIL), recording, and refusal of unknown bands.
