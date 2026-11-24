# 619 — 611's H1 → board multipliers, re-derived on 618's validated h1-trainval surface

[21dd4c], 2026-10-03 ~00:00Z (coordinator's 619 item). Read-only scoring of runs that already existed: the gate pairs
618 ran (`618/work/g1`, `g2nomde`, `g2mde`, `g3cp`, 3 seeds each). `score619.py`, `work/score619.json`.

## Method
- The same member functions 611 used through 537's price537: 528's `rows_at` (fid / jac / reach / nmae per target on
  521's truth pools, realized order), averaged over targets and over 3 seeds.
- Only the state changes: lib521's REF / EFF / ST435 point to 618's reference, faithful key and g2mde base (as
  `618/price618.py` patches them).
- **KA:** the f 1 column equals 618's gate (`gate618.json`, cell-eval2's own member code through evaluate323.score) on
  every member of every pair: zs −0.00163; MDE nmae −0.01473 and reach +0.00050; cp fid −0.00408 and jac +0.00872.
- Seed SDs are tiny: ≤ 2e-4 on nmae, ≤ 4e-5 on fid / jac. Reach for the MDE pair has SD 4e-4 … 2e-3. So every
  multiplier below is precise for this surface; its uncertainty is between panels.

## Board vs 618 (raw Δ; board / 618) beside 611's board / h1-test at f 0.35
| pair | member | board Δ | 618 f 1 | f 0.5 | **f 0.35** | cmp | h1-test at f 0.35 (611) |
|---|---|---|---|---|---|---|---|
| zs (446bj2z − 446bj2) | nmae | −0.00218 | −0.00163 (1.34×) | −0.00110 (1.98×) | **−0.00186 (1.18×)** | −0.00242 (0.90×) | 1.18× |
| coarse MDE (433a − 417e) | nmae | −0.01679 | −0.01473 (1.14×) | −0.01705 (0.98×) | **−0.01610 (1.04×)** | −0.01721 (0.98×) | 1.59× |
| coarse MDE | reach | +0.00240 | +0.00050 (noisy) | −0.00332 (sign ✗) | **−0.00163 (sign ✗, 3 / 3 seeds)** | −0.00630 (✗) | +0.00207 (1.16×, sign ✓) |
| cp call set (WORN on 618) | fid | −0.00358 | −0.00408 (0.88×) | −0.00451 (0.79×) | **−0.00468 (0.77×)** | −0.00468 (0.77×) | 1.96× |
| cp call set (WORN on 618) | jac | +0.00399 | +0.00872 (0.46×) | +0.00778 (0.51×) | **+0.00702 (0.57×)** | +0.00795 (0.50×) | 0.48× |

## Reading
1. **nmae, delivery / magnitude levers:**
   - **618 at f 0.35 is ≈ board-scale.** zs reads 1.18× (the same as h1-test) and the coarse MDE 1.04×, against h1-test's
     1.59×.
   - So h1-test under-reads add-based nmae levers. 618 reads them at about the board's size.
   - **For October prices of add-based / magnitude levers (630's UP pass, the α / DN stacks):** use 618's f 0.35 with
     nmae × ~1.0–1.2, not h1-test's × 1.2–1.6.
   - f 0.5 on 618 is noisier across levers (zs 1.98×, MDE 0.98×). f 0.35 stays the regime to read, as in 611.
2. **Reach does NOT transfer for delivery pairs on 618.** The coarse MDE's reach is mis-signed at f 0.5 and f 0.35 on 3 / 3
   seeds, where the board moved +0.0024.
   - h1-test's single sign-right match (611) now looks like luck rather than a calibration.
   - **Reach stays unpriced for delivery levers on either surface.** Price their overall-equivalent on nmae (and fid / jac
     where calls move), and state reach as unmeasured.
3. **Call-set pair (worn on 618, because 528's β was selected on trainval):**
   - jac × ~0.5 holds on both surfaces (618 0.46–0.57×, h1-test 0.48×).
   - fid differs: 618 over-reads the board's fid cost (0.77–0.88×), while h1-test under-reads it (1.96×).
   - So fid's call-set multiplier is panel-specific, between 0.8× and 2×. Price call-set fid costs with that range, not a
     point.
4. **The evidence is still three pairs.** Each multiplier rests on one board pair, and the panel-to-panel spread (h1-test
   vs 618) is now the honest uncertainty: nmae 1.04–1.59× for the coarse MDE.

## For the coordinator's October numbers
- **630 and the α / DN stacks** (priced on 618 next, PRED-618-stacks.md): board nmae ≈ 618's f 0.35 nmae × 1.0–1.2.
  - Reach: unpriced.
  - fid / jac: ≈ 0 for these rank-neutral / delivery arms (630 is rank-neutral by construction; DN / ×1.25 move few calls,
    which the stack runs count).
- **The 611 rule stays for h1-test-priced levers, with two amendments:**
  - (i) nmae multipliers are panel-specific (1.0–1.6×), and 618's are closer to 1;
  - (ii) the "reach × ~1 for delivery" line is withdrawn.

## AMENDMENT (2026-10-04, after K20's board read; [21dd4c], coordinator's ask): call-set changes read at f ≈ 0.61, not f 0.35
- **Withdrawn:** "fid's call-set multiplier is panel-specific, between 0.8× and 2×, price call-set fid costs with that range"
  (Reading 3) and its use at f 0.35.
  - K20 (K 2,000 vs 446bj2z, 10-04) moved fid −0.00169 raw on the board.
  - H1 at f 0.35 said +0.0035, the wrong sign, and so does every level ≤ 0.5.
  - Pricing call-set changes at f 0.35 is void.
- **The level the board matches for call-set changes, from 446bj2zcp + K20: f ≈ 0.61 (0.60–0.62) on h1-test. It is NOT f 1.**
  - f 1 gets K20's fid sign right but over-reads it 6× (−0.0100 vs −0.0017) and mis-signs its jac (−0.0029 vs +0.0002).
  - At f 1, cp needs fid × 2.83 and K20 × 0.17: there is no single call-set multiplier there.
  - At f 0.61 one fid multiplier fits both pairs (cp 2.70, K20 2.69), their jac ratios agree (0.43 / 0.49), and every priced
    sign of the four board DE pairs is right (cp, K20, zs, coarse MDE).
  - Source: `experiments/662-callset-regime/RESULTS.md`, `fit662.py`.
- **Call-set multipliers at f 0.61:**

| member | h1-test | 618 (cp only, WORN) | so |
|---|---|---|---|
| fid | × 2.7 (cp and K20) | × 0.81 | **× 0.8–2.7, panel-specific** (replaces × 0.8–2) |
| jac | × 0.43 (cp) / 0.49 (K20) | × 0.49 | **× 0.43–0.49** (≈ the old × 0.5) |

- **Delivery-lever nmae at f 0.61:** × 0.92–1.37 across both panels (h1-test zs 0.92 / MDE 1.37; 618 zs 1.30 / MDE 1.05).
  - At f 0.35 the range was × 1.04–1.59, so f 0.61 is no worse for delivery levers.
  - Reach for delivery pairs stays mis-signed at f 0.61: unpriced, as before.
- **Not settled:** which level reads nmae for strongly level-dependent magnitude levers (DN). 661's f 0.35 predictions stand
  as registered; an f 0.61 alternative may be added before their reads (662 RESULTS, Consequence 4).
- **Registered out-of-sample tests:**
  - K16 (10-06), `662/PRED-662-K16.md`: I-662 −0.0067 overall vs the f 1 reading −0.0030;
  - CPH (10-07), `662/PRED-662-CPH.md`: the level does not move CPH's price (≤ 0.0002 overall).
