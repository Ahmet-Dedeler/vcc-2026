# 702 follow-ups: centred (702c) and per-context (702x). DRAFT preregistration, conditional on 702's read

Written 2026-10-14 ~23:58Z by the science session [bf003b]. At this time no 702 entry exists, and neither follow-up has a
package. Exactly one of them is built, chosen by the rule below. Both are the PIE-bulk family's tries 2 and 3; tries are
counted.

## Which one (fixed now)
| 702's bin (read702.py) | next arm |
|---|---|
| TRANSFERS | **702x** (per-context), only if ctxspec700 shows PIE's context-specific part is informative: own-context wins > 0.55 AND interaction cosine > its permutation null by ≥ 0.02. Otherwise none: 702 enters October as is. |
| WEAK | **702c** (centred): the variant that removes the pds cost the instruments saw. |
| NULL / HARM | **702c** only if the read's pds and mse moved in opposite directions (the uncentred blend's mse-up / pds-down trade); otherwise none, and the PIE bulk family closes. |

## 702c — centred PIE rows
- = ALLSW + VCC_BULK_LOAD = `work/dump_pie702c_w025_seen.npz` (sha256 3d32944b11c8b622938d548c5902f919ccdc759dd2ba9653e4b9a54547267950):
  `make_pie_dump700.py --center` (each gene's mean over the panel's targets removed from PIE's rows before the blend), same
  w 0.25, seen rule, 5 CPM rule, and norm preservation as 702. Re-levelled r1 → ntclevel → r2 (build702.sh with PIE702_TAG
  `pie702c` and this dump).
- Held-out (found AFTER 702's reads, so it is a post hoc variant):
  - h1-test all ew 0.1031 / mini-pds 0.8867 (base 0.0961 / 0.8847; uncentred 0.1046 / 0.8812);
  - K562 all ew 0.0377 / mini-pds 0.7387 (base 0.0377 / 0.7370; uncentred 0.0391 / 0.7357).
- Prediction (518's identity, k 0.5; pds by span): Δmse +0.0000 … +0.0026, Δpds +0.0006 … +0.0007, **primary +0.0006 …
  +0.0033**.
- Bins: as 702's, on Δ(pds + mse) vs ALLSW (TRANSFERS ≥ +0.0015; WEAK +0.0005 … +0.0015; NULL ±0.0005; HARM ≤ −0.0005).

## 702x — per-context PIE rows (672's VCC_BULK_LOAD_CTX)
- = 702's env + `VCC_BULK_LOAD_CTX=work/ctx_pie702u_w025_seen.npz` (sha256 0e0050eaf90ee33aa02eb94c2efb0fb6231ff9e55ead7debf55533649aa4e4b5). Each context's ballast
  render takes its OWN context's PIE row in the same blend; carriers, tercile keys and audits keep the shared dump.
- It needs 672's patch applied with its KA0/KA1 (the coordinator's offer), then r1 → ntclevel → r2.
- Its gate is ctxspec700 (Replogle K562/RPE1/Jurkat/HepG2, PIE-clean): see the table above.
- Bins: vs 702's own entry (one variable: shared → per-context), on Δ(pds + mse): ≥ +0.0010 adopt per-context; ±0.0010
  keep the shared row (simpler on the day); ≤ −0.0010 keep shared and record.
