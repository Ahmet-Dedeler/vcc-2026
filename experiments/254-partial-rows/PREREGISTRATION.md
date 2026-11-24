# PREREGISTRATION — `cd4svd5+partial`: CD4's partial-condition rows as their own residual on top of `cd4svd5`. Drafted 2026-09-04 by 254; FROZEN 2026-09-09 23:10Z by the lead, before the entry exists.

**Status: FROZEN 2026-09-09 23:10Z.** Submitted as the second 2026-09-09 UTC slot. The draft below is unchanged except this header and §7.

## 1. The arm — one variable on top of `cd4svd5`
The package `submissions/254-cd4svd5-partial-f034.vcc` is `cd4svd5` (`VCC_CD4_DEFLATE_K=5`, 249) with ONE more variable: `VCC_CD4_PARTIAL_G=0.34`. The 15 A/B/C targets carried by one or two of the three frozen CD4 conditions (4 one-condition, 11 two-condition: ANKRD52, EAPP, GMEB2, KMT2A, PDK1, PIN1, PLAGL2, PRMT7, SAP30L, STK3, TAF4, TIMM50, TMF1, UPF3B, ZKSCAN2) — the rows production's `complete mean3 only` rule zeroes — enter as their own residual: the mean over the conditions that carry them, `gene_norm`, the banks' own rescale, no re-normalisation (246 `lib246.cd4_partial(equalise=False)`), sized by the carrying-row rule (218 §7.1) at g = 0.34 against the same whitened base (solved fraction 0.11244119388619626 on A/B/C; the shipped all-rows rule would solve exactly 0.0), accumulated LAST into the one merged add (deflated CD4 at its re-solved fraction 0.34104860341941284, the H1 union at 0.19207161964444322, then the partial rows), one global preserve. The complete-row CD4 object, its k = 5 deflation (`U_5` `f06e0262…`, energy removed 0.19866046506100588) and its solved fraction are byte-identical to `cd4svd5`'s; the union's frame is byte-identical; the only CD4-block scalars that move are the merged preserve pair (0.9115882187626339 → 0.9101933222360135). α 0.982, operator `e2f65ecf…`, emitter `rc_callset`, `VCC_H1_UNION_G` 0.28, `VCC_CD4_TARGET_F` 0.34 — everything else the shipped configuration.

## 2. Why, and why it is NOT a cleared arm
246 §3.7: on 162's PLOLO (three cancer folds) the partial rows on top of the shipped stack read macro cosine +0.587%, aggregate +1.156%, `pds` +0.00178 (3/3 folds each); target-shuffled ≈ 0; the row-matched self control reproduces the cosine (+0.57%) but not the `pds` (−0.0001) — the `pds` is CD4-specific, the cosine is "more residual on rows that had none". 254 reproduced that row through the builder to 0.0e+00 and measured the composition the package actually carries (the RESOLVED deflation + the partial rows): vs `h1union` macro −0.222%, aggregate +0.702%, `pds` +0.01405 (3/3, leave-top-3-out 61.7%); **the increment over `cd4svd5` is exactly the partial-row row: macro +0.586%, aggregate +1.162%, `pds` +0.00178 (3/3: +0.00111 / +0.00051 / +0.00373), Δq_w −0.00017, leave-top-3-out of the increment −6.3%** (top-1 carries half; `kolf2` apart +0.018% / +0.00033). Priced on 246's rate card the increment is **+0.00089 / +0.00109 / +0.00175 overall at `pds` transfer 0.7× / 1.0× / 2.0× — 0.47× / 0.57× / 0.92× the +0.0019 shape-only gate, and inside the 0.0045 ordering noise on every reading. It is not a cleared arm (246: "NOT build-worthy alone").** What a slot on it buys is the second board point on the deflation class with one more variable whose PLOLO sign is 3/3 — and the answer to whether 15 rows of CD4 information on rows no residual reaches read on the board at all (the shape `h1union` had on 77 rows, 222). It should be spent only if the `cd4svd5` point lands AND the lead has no candidate with a margin over the gate.

## 3. Predictions vs the `cd4svd5` board point (246 §4's rate card on the increment; `pds` from PLOLO at three transfer factors)
| member | predicted Δ overall vs `cd4svd5` |
|---|---:|
| `pds` | **+0.00047** (0.7×, `h1union`'s realised factor) · +0.00066 (1.0×) · +0.00133 (2.0×) |
| `mse` | +0.00019 (macro +0.586% through `cd4f034`'s card; the package stream says raw `mse` moves by -0.0030 raw (+0.00052 overall) at k_cos 1 — see §5) |
| `reach` | +0.00013 |
| `nmae` | +0.00005 (rate card) / +0.00005 (exact identity, `e` 0.75 — a magnitude bound, not a sign) |
| `fid` | +0.00004 |
| `jac` | +0.00001 |
| **overall** | **+0.00089 / +0.00109 / +0.00175** |

Board raw `pds` predicted **+0.0012** raw at 0.7× (range +0.0012 … +0.0036) on top of `cd4svd5`'s realised raw `pds`. Delivered amplitude vs `cd4svd5` off the package: +0.227% (bound ±0.5%).

## 4. Falsifiers — pre-committed (each read against the `cd4svd5` board row, the deterministic scorer)
1. **Board raw `pds` rises by less than +0.0006 over `cd4svd5`'s** (half the 0.7× prediction): the partial rows' `pds` does not transfer; close the partial-row family (and note the increment was 44% top-1-concentrated on PLOLO).
2. **Raw `mse` worsens by more than +0.0010 against `cd4svd5`'s** (the arm ADDS +1.16% aggregate cosine, so raw `mse` should not rise; the identity at the package's own terms predicts 0.81610 vs `cd4svd5`'s 0.81908 at k_cos 1).
3. **Overall Δ vs `cd4svd5` outside [−0.0010, +0.0040].**
4. **Config diff against `249-cd4svd5-f034.h5ad.meta.json` is not exactly `cd4_partial_g` plus its consequences** — the `effects.cd4_partial` block and the merged preserve pair; `cd4_deflate_k` 5, `effects.cd4.fraction`, the CD4 frame, the deflation block, `operator.sha256`, α, `cd4_target_f`, `h1_union_g` and every `effects.h1_union.*` field except the preserve pair byte-identical. (Measured at build time by 254 `metadiff254.py`; recorded in RESULTS.md.)
5. **Delivered amplitude vs `cd4svd5` moves by more than 0.5%** (243's stream off the package): the arm is claimed shape-only; if amplitude moved it is not.

## 5. What the slot buys even if flat
Whether the direction channel pays on rows that carry NO other residual (the `h1union` shape) when the information is CD4's rather than H1's; and a second point on the `cd4svd5` class at one more variable. A flat result closes the partial rows (they were never priced above the gate); a positive one says the carrying-row rule on small row sets transfers and re-opens the dose g and the D/E/F composition (`VCC_CD4_DEFLATE_BASIS=V5_abc` + `VCC_CD4_PARTIAL_G`, which the builder accepts together).

## 6. Discipline
Joke model name, **empty description** (§14). Never in flight together with another entry; only after the `cd4svd5` point is read, and only if the lead decides the information is worth a slot given the margin in §2.

## 7. Reading bins, fixed 2026-09-09 23:10Z before the entry exists
**Reference for Δ:** `cd4svd5`, entry `Xv0ayvuoHpsKYWv2iu42`, **0.2160944022**, same panel `vcc2026-val-1` and same
anchor set `vcc2026-valA-r4+valB-r4+valC-r4` as every entry since (verified 23:09Z). One variable:
`cd4_partial_g` absent → 0.34 (verified on the two packages' meta). Package SHA-256
`eea5ed1802a2337dc3d69b38bec4e7962835ba15b194c650015049b7ff035bf0`, novelty guard OK.

**Why now.** The 09-05 decision tree (`experiments/249-cd4-deflate/DECISION-TREE.md`) said `+partial` "needs its own
transfer reading"; the base then moved to normrel/cap (312–319) and the reading was never bought. The mechanism is
untested on the board and unclosed by later evidence. This is a transfer reading of the mechanism on its own base,
**not** a candidate package: the displayed row will fall to ≈0.216–0.218 (below 319 and 319d3) and is not to be
corrected by re-submission.

| bin (Δ overall vs 0.2160944) | reading |
|---|---|
| **Δ ≥ +0.0010** | the partial rows transfer at ≥ the drafted 1.0× card. Register `VCC_CD4_PARTIAL_G=0.34` as a one-variable rider on the 319 composition; price it on 323 before any slot. |
| **−0.0005 … +0.0010** | inside the drafted band's floor or noise: the increment does not earn a slot on the new base. Closed unless a cheaper instrument resolves it. |
| **Δ ≤ −0.0005** | the partial rows hurt on the board; the flag stays inert for good. |

**Secondaries, no verdict:** `pds` (drafted +0.00047 / +0.00066 / +0.00133 overall at 0.7× / 1.0× / 2.0×), `mse`
(drafted +0.0002 … +0.0005), the four DE members expected flat. Read with `scripts/read_board_entry.py` against
the captured reference record. Capture the new record immediately after it publishes.
