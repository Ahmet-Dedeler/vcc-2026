# Contradictions — the summary layer vs the primary sources (2026-09-04)

Scope: claims in `AGENTS.md` §8.x and in the live header of `docs/CURRENT_STATE.md` (lines 1–121) that a
primary source contradicts or does not support. Ranked by decision relevance. Both sides are shown; nothing
is resolved here. Key: A=`AGENTS.md`, CS=`docs/CURRENT_STATE.md`, BL=`docs/BOARD-LEDGER.md`,
RB=`docs/DEF-RUNBOOK.md`, E<n>=`experiments/<n>-*/`, R<n>=E<n>`/RESULTS.md`, BR<n>=E<n>`/BOARD-RESULT*.md`,
LB=`experiments/245-read-first/work/ranks.json`.

| # | claim (path:line) | what the primary source says (path:line) | what depends on it |
|---|---|---|---|
| 1 | CS:83 — "Leading candidate for the next one: 223's DLD-1 rider (+0.00228)" | `E237/submit-armed.log:1`: `ARMED pkg=submissions/229-stack6-g035.vcc when=202609050000Z` (pid 6813, alive at 11:57Z); `E237/PREREGISTRATION.md` §5 names that package; R229:17-27 measured it at −0.00219 and marked it NOT RECOMMENDED; R225:21 has DLD-1 at +0.00228 | the 2026-09-05 00:00Z slot: which package fires, and what the reading means |
| 2 | A:1216 — pricing rule `board_cosine_gain ≈ 1.80 × PLOLO_cosine_gain` | BR222:61-64: the factor is 1.384 on `cd4f034` and 1.000 on `h1union`, "not a validated instrument"; R219:22: 1.38×, range 1.17–1.39 | every direction-arm price and band |
| 3 | A:1218 — `pds cost = 2.0 × the PLOLO prediction` | BR222 §2: `h1union`'s `pds` realised 70% of its 3-cancer prediction (+0.00344 vs +0.004887), i.e. under-, not over-realised | the `pds` leg of every direction-arm price |
| 4 | A:1064,1068 — "our transfer cosine 0.079249"; "a single direction arm must deliver +31% cosine" | BR222:19: implied cosine 0.09947 after two arms; A:1196-1197 (§8.4 itself) says the +31% bar is 2.0–2.3× too high | sizing of any direction arm read from §8.2 |
| 5 | A:705 — direction "attacks: Tahoe-100M (184), per-gene magnitude shape (188)"; A:706 — ordering "−1.92σ … attacks: 185, 189, 190" | R184:1 (Tahoe: the operator "does not exist in our CRISPRi data"); BL:125: −1.92σ → +0.70σ after `rc_callset`; R220:1: ordering-key family CLOSED; R231:19-22: whole `reach` prize +0.002689 vs +0.0089 floor | where a new agent is sent |
| 6 | A:1090 — "the remaining direction levers are exactly two": estimator (194) and eigenvalue floor (193) | RB:206-209 lists 193 and 194 as closed; R193:1 "The geometry is CLOSED"; R194:1 "the plain pseudobulk difference wins" | allocation of the direction programme |
| 7 | CS:116-120 — "The two live hypotheses, both running: 224 … 225" | R224 exists (verdict R224:12-17); R228:1 kills 224's gain on `nmae`; R225:1 "The space is swept" | what is running; the closed list |
| 8 | CS:71-74 — "221 has both points and owns the reconciliation" of the 1.384/1.000 factor | R221:776: the PLOLO cosine→board question "is a different instrument" and was not resolved; `E243/MANDATE.md:1` now owns it, running | who resolves the factor before any arm is sized |
| 9 | A:1200 — a fifth donor line is worth "+0.0140–0.0163, 1.6–1.8× above the corrected bar" | R218:1 "Independence does not predict contribution"; R225:1 DLD-1 is the only source clearing alone, sources saturate at ~1.8×; R229:22-23 ARM6 −0.00219, ARM4 −0.00589 | donor-acquisition spend |
| 10 | CS:19-20 — world rank per member "jac 1st, mse 10th, reach 10th, pds 7th, nmae 59th, fid 129th" (04:20Z, 656 entries); A:730-731 (§8.1b) "`reach` world 7th" | LB:8-48 (11:49Z, 664 entries): `pds` 8, `mse` 6, `nmae` 56, `fid` 106, `reach` 7, `jac` 1 | per-member headroom statements |
| 11 | CS:16 (also BL:203, RB:200) — `cd4f034` overall 0.2049216 | BR210:5: 0.2049195830; BR222:3: baseline 0.2049196 | `h1union`'s Δ (+0.0069114 uses 0.2049196) |
| 12 | BL:94 (the §3b table BL says pricing models must be fit against) — `cd4f034` "pending, pred +0.0044 / +0.0057" | BR210:5,20: realised +0.0130557, 3.0×/2.3× the prediction | any model fitted to BL §3b |
| 13 | BL:134-135 and BR185:50 — h1-test→board `reach` factor 1.290, "h1-test understates the board by 3×" | R221:23,444: 0.769× on the real scorer, h1-test overstates | every `reach` prediction for order arms |
| 14 | RB:203 (§3.1 PerturbMatch row) — "Experiment 214 is testing whether it turned negative because 207 tripled the CD4 dose" | R214:5-8: refuted; PerturbMatch's `pds` value rises with the dose, weak complements | the PerturbMatch row of the blind decision table |
| 15 | RB:264-309 (§5.-2, source ranking on 199's pessimistic bundle with 30 essential targets) and R236:3-11 (+0.0066…+0.0104) | RB:358-369 (235: essential-free draw changes the ranks, `k562_essential` 30→0); R238:3-10: robust value falls to 16%, `pds` changes sign, "the strategy claim is retracted" | which October sources are wired and how they are priced |

Notes.
- Rows 11–15 are primary-vs-primary disagreements that the summary layer inherits; they are listed because the
  summary layer repeats one side of each.
- `submissions/SUBMITTED_LEDGER.json` holds 5 entries ending 2026-09-02 (lines 5-33); the five entries scored since
  (BL:92-94; BR180:3; BR222:3) are absent, and R242:45 demonstrates the novelty guard passes an already-scored
  package as a consequence. Not a summary-layer claim, but it affects row 1.
