# Contradictions v2 — the head documents vs the primary sources (2026-09-04 17:05Z)

Scope: `AGENTS.md` §8.x and the live header of `docs/CURRENT_STATE.md` (lines 1–127), re-checked against the primary
sources after today's edits, plus the primary-vs-primary rows 245 listed that the head documents inherit. Nothing is
resolved here; both sides are shown. Key: A=`AGENTS.md`, CS=`docs/CURRENT_STATE.md`, BL=`docs/BOARD-LEDGER.md`,
RB=`docs/DEF-RUNBOOK.md`, E<n>=`experiments/<n>-*/`, R<n>=E<n>`/RESULTS.md`, BR<n>=E<n>`/BOARD-RESULT*.md`,
LB=`experiments/253-read-first-v2/work/ranks.json` (board fetched 16:53Z), RF1=`docs/READ-FIRST.md` (245's file as moved).

## Status of 245's fifteen rows

| 245 row | status now | where |
|---|---|---|
| 1 next-slot candidate | RESOLVED: the armed 229-stack6 slot is recorded in the header | CS:87-90 |
| 2, 3 pricing rule 1.80× / `pds` cost 2.0× | RESOLVED by banner: §8.4 now carries "SUPERSEDED BY 243 … kept for the record" above the rule | A:1222-1239, 1244-1246 |
| 4 §8.2 cosine 0.079249 / +31% | NOT resolved in place; §8.4 corrects the bar, 243 measures the cosine | row 5 below |
| 5 §8.1 table "attacks" | NOT resolved | row 2 below |
| 6 "exactly two levers" | NOT resolved | row 3 below |
| 7 "two live hypotheses" | banner only (CS:3-4 marks the file as history); text unchanged | row 10 below |
| 8 "221 owns the reconciliation" | banner only; `AGENTS.md` now carries 243's answer at A:1222-1239 | row 11 below |
| 9 fifth donor line +0.0140–0.0163 | NOT resolved | row 4 below |
| 10 world ranks | banner only; the board moved again today | row 12 below |
| 11 `cd4f034` 0.2049216 typo | RESOLVED in the primaries: BL:203 and RB:251 read 0.2049196; RF1:21 still repeats the old comparison | row 15 below |
| 12 BL:94 `cd4f034` "pending" | NOT resolved (primary) | row 6 below |
| 13 `reach` factor 1.290 | NOT resolved in BL/BR185; A:980-988 records 221's 0.769× | row 7 below |
| 14 PerturbMatch "214 is testing" | RESOLVED by the §3.1 update that re-decides the row; the old row text remains | RB:232-238 vs RB:254 |
| 15 §5.-2 source ranking on essential-contaminated bundles | PARTLY: RB:416-421 carries 235's caveat, A:905-954 the rule; the table stands | row 9 below |
| note: `SUBMITTED_LEDGER.json` 5 entries | RESOLVED: 10 entries, five backfilled by 247 | SL:31-77; R247:261-275 |

## Rows, ranked by decision relevance

| # | claim (path:line) | what the primary source says (path:line) | what depends on it |
|---|---|---|---|
| 1 | CS:89-90 — "DLD-1 alone (+0.00228, 225) is the alternative if the gate is validated" | R248:1,417-421: the DLD-1 rider through the builder prices +0.00216, 1.14× the shape-only gate against a required 1.5×, leave-top-3-out retention 25.2% against a 50% bar — "RIDER NOT BUILDABLE"; R229:22-24: ARM3 (H1 pair + `dld1_lowmoi`) +0.00067 | what the second 2026-09-05 slot or the next one carries if the 237 reading validates the gate |
| 2 | A:714-715 (§8.1 table) — direction "attacks: Tahoe-100M (184), per-gene magnitude shape (188)"; ordering "−1.92σ … attacks: 185, 189, 190" | R184:1 (the operator "does not exist in our CRISPRi data"); R228:1 (per-gene gain killed by `nmae`); BL:125 (−1.92σ → +0.70σ after `rc_callset`); R220:1 and A:1001-1025 (ordering closed); R231:19-22 (`reach` prize +0.002689 vs +0.0089 floor) | where a new agent reading §8.1 is sent |
| 3 | A:1099-1103 — "the remaining direction levers are exactly two": the estimator (194) and the eigenvalue floor (193) | R193:1 "The geometry is CLOSED"; R194:1 "the plain pseudobulk difference wins"; RB:257-260 lists both as closed; the arms that did move the board were CD4 dose (BR210:5) and the H1 union (BR222:6), neither on that list | allocation of the direction programme |
| 4 | A:1208-1210 — a fifth independent donor line is worth "+0.0140–0.0163, 1.6–1.8× above the corrected bar" | R218:1 independence does not predict contribution; R225:1,21 DLD-1 is the only source clearing alone at +0.00228, saturation ~1.8×; R229:22-24 ARM6 −0.00219, ARM4 −0.00589 | donor-acquisition spend |
| 5 | A:1073-1077 (§8.2, dated 2026-09-03) — "our transfer cosine 0.079249"; "a single direction arm must deliver +31% cosine" | R243:28-30: measured board cosines `rc_callset` 0.079249 → `cd4f034` 0.097850 → `h1union` 0.107332; A:1205-1206 (§8.4) says the +31% bar is 2.0–2.3× too high | sizing of any direction arm read from §8.2 |
| 6 | BL:94 (§3b, the table pricing models are told to fit against) — `cd4f034` "*pending*, pred +0.0044 / +0.0057" | BR210:5,20: realised +0.0130557, 3.0× / 2.3× the prediction; SL:35 carries the scored value | any model fitted to BL §3b |
| 7 | BL:134-137 and BR185:50 — h1-test→board `reach` factor 1.290, "h1-test understates the board by 3×" | R221:444: 0.769× on the real scorer, h1-test overstates; A:980-988 records it as a "correction candidate, not a settled overturn" | every `reach` prediction for a designed-order arm |
| 8 | BL:195 — `h1union` "submitted 2026-09-04 03:52 UTC, awaiting score"; BL:213-214 — "my board estimate ≈ +0.0035 overall → ≈ 0.2084" | BR222:1-6: scored 0.2118310, Δ +0.0069114, rank 3; BR222:35-38: the reissued band was wrong, the original +0.0055–0.0105 was right | the ledger's own record of the current incumbent |
| 9 | RB:315-360 (§5.-2) — the October source ranking, computed on 199's pessimistic bundle (30 essential targets, RB:434) | RB:416-421: on an essential-free draw `k562_essential` 30 → 0, the ranks change; R238:3-10: the October four are 70.7–90.4% essential-carried, robust value 16%, `pds` sign flips, "the strategy claim is retracted" | which October sources are wired and how they are priced |
| 10 | CS:122-126 — "The two live hypotheses, both running: 224 … 225" | R228:1 kills 224's gain on `nmae`; R225:1 "The space is swept"; today's running set is 252 and 254 (`E252/MANDATE.md:1`, `E254/MANDATE.md:1`) | what a reader of the header believes is running |
| 11 | CS:78 — "221 has both points and owns the reconciliation" of the 1.384 / 1.000 factor | R221:776: a different instrument, not resolved there; R243:1,44-47: resolved as a class constant (1.36 re-dose, 0.94 new source), now at A:1222-1239 | who to cite for the transfer factor |
| 12 | CS:9-12 — rank 3 of 656, leader 0.22431, gap +0.01248; CS:23-24 — `mse` 10th, `reach` 10th, `pds` 7th, `nmae` 59th, `fid` 129th; A:739-740 — "`reach` world 7th" | LB:3-7,53-57 (673 entries, 16:53Z): leader `cqawesome` `m0904b` 0.22454757, gap +0.01272; LB:9-49: `pds` 8, `mse` 7, `nmae` 58, `fid` 107, `reach` 8, `jac` 1 | per-member headroom statements |
| 13 | A:891-893 — the shared PLOLO silent zero: "Status 2026-09-04 11:05: NOT FIXED … because 241 is running against it" | R241 has landed (R241:502-508); the fix is now 252's mandate, running (`E252/MANDATE.md:1`); until it lands the status line is right about the code and wrong about the reason | whether the next agent hand-sizes partial-coverage sources (A:901-903) |
| 14 | RB:92 (§1 table) — builder `scripts/make_submission_pm_cd4.py` at HEAD `ab8aa3e` sha `6cc9150…` | RB:41-42 says three commits stale ("249 moves it again"); R251:475: the builder is `10d0177e…` after commit `c321f76` (working tree matches: `shasum` 2026-09-04 16:5xZ, `10d0177e33e5…`) | the §1 monthly check, and the meta's `hashes.implementation` on the day |
| 15 | RF1:21 — "`cd4f034` overall: 0.2049216 (BL:203; RB:200) vs 0.2049195830 (BR210:5)" | BL:203 and RB:251 now read 0.2049196 (commit `6985807` "fix the 0.2049216 typo"); the disagreement no longer exists in the cited files | nothing, once v2 replaces RF1 |

Notes.
- Rows 6–9 and 14 are primary-vs-primary; the head documents inherit them.
- `AGENTS.md` §8.1b/§8.1c and CS:3-4 now carry banners; a banner leaves the contradicted sentence in place, so a reader
  who lands on the sentence (A:714, A:1099, A:1208, CS:122) still reads the superseded claim without the correction.
