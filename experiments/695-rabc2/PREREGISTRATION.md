# 695 — R-ABC2: the October recipe after the 10-17/18 reads, on def_day, on A/B/C (preregistered 2026-10-16 ~09:50Z, readiness [4592ee])

For 10-19 s1 (coordinator [2d6f46]'s slot plan, docs/PLAN-2026-10-15.md 09:33Z, 4a984060). Written before any of its inputs
exists: the formula, the bins and the revert order are fixed now. The inputs (which switches were adopted, and their board Δs)
are filled in "Resolved inputs" below from the coordinator's recorded verdicts, and that section is committed before the build starts.

## What is built
`october_recipe.py command` (plain) as resolved after the reads, run by def_day.sh on 180's def-bundle with 453's operator 66610be7
pre-placed (same-panel branch), exactly as R-ABC (694). It is R-ABC's env plus, for each switch that its read adopted:

| switch (read, slot) | adopted iff (its own prereg) | export | acts on |
|---|---|---|---|
| lodo_strength θ 0.25 (705, 10-17 s1) | 705 primary > −0.0005 (not HARM) | `DEF_DAY_LODO_STRENGTH=0.25` | the LODO key → dump weights (bulk) |
| strong-tier swap hook (707, 10-17 s2) | 707 Δ(pds+mse) ≥ +0.0010 (TRANSFERS) | `DEF_DAY_DUMP_HOOK=experiments/707-strong-tier-direction/dumphook707.sh` | the dump's strong-tier direction (bulk) |
| call-set K (711 / ku38, 10-18 s2) | 711's fixed K rule picks K ≠ 2,600 | `DEF_DAY_RC_CALLSET_CPM=<K − 400>` (3,400 → 3000; 4,200 → 3800) | the DE call set (re-draws pass 1) |

H1 stays x2 (`h1_stack` DECIDED x2): R-ABC2 is R0's recipe, and an H1 move can only come from 708's day-of rule as a later entry.
710's read is reported only. All three exports are env-only in def_day.sh (530, 690g), so 690k's R-KA proof (def_day.sh 11ff102e)
stands. If NO switch is adopted, R-ABC2 equals R-ABC (already scored): the driver stops "nothing new" and the slot is freed.

## Prediction
pred = R-ABC 0.2743318 (entry PooEbFClgTLMQoXiOCQP) + Σ over adopted switches of that switch's own board Δoverall vs its base
(705 and 707: vs ALLSW's entry; K: the chosen size's Δoverall vs 446bj2z, i.e. vs K 2,600). Additivity is the null hypothesis:
R-ABC composed UPRN + DNMAG at 95 %; ALLSW's DE switches added exactly on the cap6+P2g bulk.

## Bins on D = R-ABC2 − pred (overall, board)
Noise: every adopted switch re-draws pass 1 (a bulk change moves ballast K; a K change re-draws the call set), so the read carries
the rule-13 re-roll sd ≈ 0.00025 overall. ±0.0005 = 2 sd.
| bin | D | consequence for October (R0 builds from docs/october-decisions.json) |
|---|---|---|
| ABOVE | ≥ +0.0005 | ships as decided; record the interaction |
| **COMPOSES** | (−0.0005, +0.0005) | ships as decided |
| BELOW | ≤ −0.0005 | the coordinator's revert rule: revert the switches R-ABC2 added, **latest adopted first** (K, then 707's hook, then 705's θ), one at a time, until 10-21's R-ABC3 read |

## Falsifiers / checks (pre-slot unless marked)
- The recipe check: the resolved env is R-ABC's env plus ONLY the adopted exports above, with those values, and the same args.
  Anything else (e.g. a DEF_DAY_H1_STACK export) stops the build.
- The content-aware config diff vs R-ABC's final meta (`scripts/def_metadiff.py`, allow/need per adopted switch):
  θ → bulk_load (renamed `_s0.25`), bulk_load_sha256, pert_scale(+sha), bulk_ctx_scale, bulk_sharpen_pt; hook → bulk_load_sha256,
  bulk_ctx_scale; K → rc_callset_cpm. Nothing else.
- Per-switch evidence: θ: def_day logs "LODO strength rider: theta 0.25", sidecar theta 0.25, dump name `_s0.25`. Hook: def_day logs
  "690g: dump hook …dumphook707.sh" and "dump replaced by the hook's" (def_dumphook_check PASS). K: meta rc_callset_cpm == K − 400.
- If K is the ONLY adopted switch: bulk identity vs R-ABC's package (bulkid658: per-(context, target, gene) sums identical), and
  on the board pds and mse raw equal R-ABC's (exact members). With θ or the hook, the bulk moves by design.
- CLIFF OK, postbuild all PASS, zero-gene audit 0, novelty rc 0, package sha logged.

## Class and tries (AGENTS 11.1c)
A composition read and the set-value rehearsal on A/B/C of switches whose own reads are single board arms: 705 (transfer, safety
read; 9/9 October-like folds), 707 (transfer; Arc-panel proxy), K (scorer mechanics; 711's rule picks among 3 board sizes with a
+0.0005 step margin and an H1 mechanism gate). The new-panel rehearsal of θ and the hook is rider690 (690), judged against R-NEW.

## Resolved inputs (filled from the coordinator's verdicts; committed before the build)
- 705: **CLEARED** (10-17 s1, lstr705 entry o7ObvpGu05zvmCwsFtiC 0.2757977; primary Δ(pds+mse) −0.0002397 > −0.0005). Adopted:
  `DEF_DAY_LODO_STRENGTH=0.25` (october-decisions 00:31:58Z). Its own Δoverall vs ALLSW = **−0.0002874**.
- 707: **WEAK, not adopted** (10-17 s2, sw707 entry AYbAiHRjuHdil2zCjlsk 0.2771181; primary Δ(pds+mse) +0.0008417 in +0.0005…+0.0010).
  dump_hook stays off for R0; 708's AMENDMENT 2 swap route stays live (not HARM).
- 711 / K: stays 2,600 (Amendment 1: condition (i) failed); a K 3,400 adoption can only feed 696 R-ABC3.
- 710 (reported only): reads 10-18 s1, after R-ABC2's build; it does not enter R-ABC2.
- → **adopted exports: `DEF_DAY_LODO_STRENGTH=0.25` only. pred = 0.2743318 − 0.0002874 = 0.2740445**; bins ±0.0005;
  BELOW (≤ 0.2735445) reverts θ. Frozen 10-17 ~01:30Z, before GO-rabc2.

## Amendment 1 (10-16 ~12:40Z; before any of the reads, before any build)
- Slot moved to **10-18 s1** (coordinator 298522ca; waiter armready 1018s1 rabc2, floor 10-18 00:01Z, prior submit-1017s2.log).
- **K resolved: stays 2,600** (`rc_callset_cpm` off). 711's condition (i) failed at f 0.87. K 3,400 can return only through 711's
  Amendment 2 dial test (ku38, 10-18 s2) and then into a later R-ABC3 (gate695's K-only path covers it). So R-ABC2 = R-ABC + θ 0.25
  (if 705 clears) + 707's hook (if TRANSFERS). The revert order for BELOW becomes: hook, then θ.
- GO-rabc2 is touched right after the 707 verdict is applied (10-17 ~02:00Z), not after ku38.
- If neither 705 nor 707 adopts, the driver stops "nothing new" and the coordinator writes NOGO-1018s1.
