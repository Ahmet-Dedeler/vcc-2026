# 674 — W2: W1 with CD4 at ×1.5 (R0's recipe + H1 ×2.5 + CD4 ×1.5, def_day / A/B/C), the joint weights guard

Written 2026-10-20 04:31Z (by the clock) by the coordinator [3c7a7e], before any build, any D/E/F data, and any W1 or W2 number.

## Why this slot
- 708's 10-16 rehearsal (r708abcdef, production operator) recommended **H1 ×2.5, CD4 ×1.5**. W1 (673) reads H1 ×2.5 alone on R0's recipe, and the
  CD4 ×1.5 part has only 457dr (+0.0031484 vs 453dr), on the September October path, before the LODO percentile key, tercilec and θ. No board
  read has the joint cell on the recipe that ships (SIXA-DECISION Amendment 2, "coverage").
- W2 differs from W1 by one export (rule 5), so W2 − W1 reads the CD4 step on R0's recipe at H1 ×2.5.
- It uses 10-21 s2, the last A/B/C slot, which nothing else needs. Guard only: it can block a 708 move, never create one.

## Arm
- W2 = R-ABC2's resolved recipe (frozen decisions 04be2af3) + `DEF_DAY_H1_STACK=x2.5` + `DEF_DAY_CD4_STACK=x1.5`. No hook. def_day on 180's A/B/C bundle.
- Config diff vs W1 within {bulk_load, bulk_load_sha256, bulk_ctx_scale} (the dump is renamed again); the only changed export vs W1 is
  DEF_DAY_CD4_STACK. The named dump's env must record VCC_H1_UNION_G = 0.7 and the CD4 ×1.5 stack (as readiness's gate673 checks for G).
  Frozen tree shas before and after. Build after W1's driver exits; no start after 10-21 12:00Z.

## Prediction and bins (fixed now)
- pred = W1's own board read + 457dr's Δoverall vs 453dr (+0.0031484).
- D = W2 − pred. **COMPOSES** |D| < 0.0005; **ABOVE** D ≥ +0.0005; **BELOW** D ≤ −0.0005. A second BELOW condition: W2 < W1 − 0.0005.
- **BELOW blocks:** no 708 move may raise CD4 above ×1 in the final (any H1 dose). The route is the same as W1's: readiness re-arms wait708
  with the block before the release (a new `WAIT708_*` switch if decide708 has none for CD4; else the coordinator applies it at SIXA time).
  COMPOSES / ABOVE: no change. CD4 ×1.5 still needs AMENDMENT 1's fold confirmation on D/E/F.
- If W1 itself reads BELOW, W2 is still read (the block on H1 > ×2 already covers the joint cell; W2 then reads the CD4 step for the record).
- Noise: same emission law, pass 1 redrawn (sd ≈ 0.00025 overall), so ±0.0005 is ~2 sd.

## Class and tries
- Class: transfer (composition weights), used only as a guard. Known to be A/B/C-positive: 438 found CD4 ×1.5 wins on A/B/C and loses on every
  October-like panel, so a COMPOSES/ABOVE read says nothing about D/E/F. Only a BELOW is informative for the final.
- Tries for the CD4 dose on the board: ×1 (many), ×1.5 (417e, 457dr, the ALLSW family), and this one.
- Reader: `read695.py --w2 <w2 entry> <w1 entry>` (to be added and dry-run before the read).

## Amendment 1 (2026-10-20 04:32Z by the clock, before any build or number): the block's route
- decide708 has no CD4 switch (only R710, for H1). Adding one two days before the release would change a script on the day path for a
  low-probability block, so no new switch. **On a W2 BELOW the coordinator applies the block at SIXA decision time:** X1 is only ever
  built through checklist 6a, which reads SIXA-DECISION first. Recorded there as Amendment 3.

## Not built (2026-10-20 11:34Z, coordinator)
- W1's gate673 false-FAILed at 11:15Z, so w2_674 stopped by its own start rule ("W2 has no reference run"). The coordinator dropped W2 instead of
  restarting it. It was the lowest-value guard (CD4 ×1.5 is already ruled out on October-like panels by 438), and gate674 carries the same allowed-set
  error. 10-21 s2 stays unspent. No W2 number exists; SIXA Amendment 3's block never applies.
