# 625 — A/B/C module readout: preregistered, approved, decode committed. No read yet.

> **AUDIT 616 (7dac7a98, [21dd4c], 10-02): the thinned Orion gate tables' DE truth was ~98 % selection artefact.** Here the hct116 / hek293t rows of the module-net table are VOID. The decisions stand: Amendment 3, sized by H1's real emitter (625 FINAL), and the Broad-16 8-octile read.

Status, 2026-10-02 ~12:40Z (council session [7772dc]):
- **PREREGISTRATION.md** with Amendment 1, sha256 `bd8fdd61…` at f235cca5. Modules: `work/modules625.npz`, sha
  `3bce7179…`.
- **Coordinator APPROVED** with four conditions (via [4d0207]):
  - (a) one base B, named by package sha and frozen before the first arm;
  - (b) the decode committed with synthetic known answers BEFORE any read. This file, below.
  - (c) the arms enter 537's queue behind every decision arm, obey the 35 GiB build floor, and each package is swept
    after capture + ledger;
  - (d) every `-m` from the MODEL-NAMES queue via the claim tool, `-d` empty.
  Ahmet was told about the ~18 public probe rows and has not vetoed.

## Decode (condition b): `decode625.py`
`--selftest` prints **SELFTEST PASS**:
- 200 random additive 16-module effect vectors (nmae, fid, reach) are recovered from simulated 18-arm board reads to
  1e-12;
- the §6 transfer rule is exact on a constructed case: flip where both members agree, a fid veto, the 0.25×median bar,
  and NULL when CLEAN4 agrees everywhere;
- a planted module-0 × module-15 interaction trips the §5 additivity check, giving the fallback mask and the ABORT
  reading.

## Partition check before any read (`partition625.py`, `work/partition625.json`): the preregistered modules stay
Broad-16 (12 working lines; reserved 4 never read), sim625's procedure at bar 0.10, 40 triples per held-out line:

| 16-module partition | held-out gate-sign Δ vs CLEAN4 | worst line | lines > 0 | triples < 0 |
|---|---:|---:|---:|---:|
| **CLEAN4 octile × expression half (preregistered)** | **+0.0166** | +0.0025 | 12/12 | 0.04 |
| 16 co-response clusters (K562 GWPS programmes) | −0.0012 | −0.0033 | 1/12 | 0.16 |
| 8 co-response clusters × expression half | +0.0000 | −0.0007 | 7/12 | 0.18 |
| CLEAN4 sign × 8 co-response clusters | +0.0093 | −0.0084 | 11/12 | 0.18 |

Modules defined by CLEAN4's own value recover the most. The readout corrects CLEAN4 where CLEAN4 is systematically wrong
for the lab, so grouping genes by how CLEAN4 calls them is the right axis. Biological co-response modules mix CLEAN4-right
and CLEAN4-wrong genes. **No amendment.**

## Final-phase use (proposal 637 §4, d4cd5f)
637 §4 proposes fusing 637's factorial reader with these modules on D/E/F (a Hadamard row per final entry), which would
decode D/E/F's own module programme. That depends on Ahmet's H-631 decision (decoding D/E/F feedback). The A/B/C readout
here would double as its end-to-end rehearsal.

## 625c — Final-phase feasibility after Ahmet's decision (2026-10-02 ~12:40Z: "we can do the risky thing" = 637 §4 approved)
**Can a single module flip be read through final-entry noise?** (`snr625.py`, `work/snr625.json`)
- Model: the delivery-faithful fold model (480/612b `small_base+mde`) with the October hook (key-gated identity
  hybrid). x_M is the change in mean-target nmae from flipping module M's CLEAN4 component.
- Units: overall-equivalent, board multiplier 1.
- Noise per decoded module in 637's fused design (a bulk switch per entry, leak ±0.0002 per entry): σ_x ≈ 0.00010 with
  16 rows, 0.00014 with 8.

| surface | design | median \|x_M\| | max \|x_M\| | detectable at 2σ_x (fused) | value of flipping every helpful module |
|---|---|---:|---:|---:|---:|
| h1 | 16 mod | 0.00013 | 0.00057 | 6/16 | +0.0004 |
| h1 | 8 oct | 0.00021 | 0.00093 | 4/8 | +0.0001 |
| hct116 | 16 mod | 0.00037 | 0.0038 | 11/16 | +0.0002 |
| hct116 | 8 oct | 0.00072 | 0.0063 | 6/8 | +0.0001 |
| hek293t | 16 mod | 0.00028 | 0.0046 | 9/16 | +0.0021 |
| hek293t | 8 oct | 0.00073 | 0.0063 | 5/8 | +0.0016 |
| k562 | 16 mod | 0.00012 | 0.00091 | 3/16 | +0.0005 |
| k562 | 8 oct | 0.00023 | 0.0012 | 2/8 | +0.0005 |

**Reading**
1. **Through the HOOK alone, perfect module knowledge is worth only +0.0001 … +0.002 overall.** The hook moves signs
   only through a small MDE add on uncalled genes. The block prior's share (called block genes: fid filler precision,
   reach order, gate genes inside the block) is not in this model and is probably the larger part. It needs a
   real-emitter estimate, a few flip arms on H1 through 537's harness, before anyone sizes the final-phase decode.
2. **The fused design reads only 3–11 of 16 modules at 2σ.** If final entries do return members, run the module decode
   as PURE prior-only entries (shared draws, exact reads) and read the bulk switches in other entries. Fusing costs
   most of the small modules.
3. **8 CLEAN4 octiles** carry larger per-module signals (they aggregate two modules) at half the entries (~10). On the
   Broad-16 transfer simulation they keep nearly all the value: +0.0165 vs +0.0172 for the preregistered 16, and 16
   CLEAN4 hexadeciles give +0.0190.
4. **If a UP pass (630) delivers larger magnitudes on uncalled genes,** sign knowledge is worth more there.
   Re-price after 630 is in the recipe.

## CORRECTION to 625c (2026-10-02 ~13:00Z, from [21dd4c]'s 612 CORRECTION 2)
625c's per-module nets came from the same mis-specified fold shape. It applied the flipped hook sign to EVERY gene; in
production the hook acts only on uncalled genes, and the head and CPM block follow the key and the block prior. **625c's
|x_M| values and readability counts are therefore not reliable** (likely inflated on the Orion surfaces). The real
instrument is 625d (`flip625.py`), production-emitter H1 runs that flip an octile in BOTH the block prior and the hook
(KA: `none` == 612's c4kg bit for bit). 637's fused-vs-pure conclusion (pure prior-only entries) still holds on design
grounds: a fused bulk switch is not K-neutral.

## 625d — Real-emitter H1 octile flips (block prior + hook), 2026-10-02 (`flip625.py`, `score625d.py`, `delivery625.py`)
Base = 612's c4kg: production rc_emit_k 0f56209d, CallSet(400, 2200, 1.9), CLEAN4 q80 block prior, coarse MDE, zs, and the
CLEAN4 hook (τ 0.05) on each target's weak-key half. Arm k flips octile k (modules625_v8.npz `octile`) in BOTH the block
prior and the hook. Prior-only arm, so the bulk sums equal `none` (asserted).
**KA:** `none` s0 == 612 c4kg s0 bit for bit (q, p, lfc, sums): PASS.
Chain: chain625d.sh ran none/0/2/5 s0 (07 s0 queued). chain625d2.sh adds seeds 1–2 for none/0/5, because the gate602
queue was long.

### Seed 0 (21dd4c's 521 members at 4 regimes; overall-equivalent by 192's INV; "w/o reach" drops reach, see below)

| octile (P range) | CLEAN4 vs H1 truth, gate genes (\|l\|-w) | f 1: Δfid / Δnmae (raw) | f 1 overall w/o reach | f 0.35 overall w/o reach | Δreach f 1 / 0.35 |
|---|---:|---|---:|---:|---|
| 0 (P < 0.440, strongest UP) | 0.634 | −0.0019 / +0.0032 | −0.0019 | −0.0018 | +0.001 / +0.016 |
| 2 (0.465–0.485) | 0.575 | +0.0001 / +0.0014 | −0.0004 | −0.0003 | +0.012 / +0.013 |
| 5 (0.525–0.548) | 0.537 | −0.0001 / −0.0009 | −0.0005 | −0.0037 | (see json) |

H1's CLEAN4 agreement by octile (f 1): 0.634, 0.596, 0.575, 0.527, 0.500, 0.537, 0.567, 0.629. **On H1 no octile is
wrong**, so a decode on H1 would read NULL. That is the expected answer for a line whose programme agrees with CLEAN4.

**Reach is excluded from every reading.** The Δreach of an octile flip on H1 comes from a few tiny-n_real targets' discrete
jumps. For octile 2 s0, 5 targets with n_real 6–25 carry all of it (ZNF32, n_real 6: +0.67), and 72 of 85 targets are
unchanged. These jumps are tied to specific genes in specific targets and do not transfer across panels. The
preregistered transfer rule already reads nmae + fid only.

### 625e — delivery calibration: how much of a sign source reaches the score (`delivery625.py`, seed 0)
| octile | source Δ gate-sign (if fully delivered) | delivered Δ gate-sign | delivery ratio | delivered sign flips, called / uncalled | delivered \|lfc\| vs real \|l\| on flips |
|---|---:|---:|---:|---|---|
| 0 (block + hook) | −0.0151 | −0.0067 | 0.44 | 1,089 / 949 | 0.054 / 0.239 |
| 2 (block only: \|P−½\| < τ) | −0.0076 | −0.0010 | 0.14 | 1,105 / 12 | 0.050 / 0.187 |
| 5 (block only) | −0.0062 | −0.0038 | 0.61 | 1,704 / 22 | 0.080 / 0.187 |

- **Exchange rate:** +0.01 of SOURCE gate-sign accuracy ≈ **+0.0011 overall** without reach (f 1; f 0.35 +0.0012), i.e.
  ≈ +0.0025–0.004 per +0.01 of DELIVERED gate-sign.
- **So:** the "responding-gene programme" prize (+0.045 gate sign, 635) is ≈ **+0.005 overall through today's routes**, not
  +0.008–0.010. Every sign source is discounted the same way.
- **Delivery coverage and magnitude multiply sign quality.** The head follows the key. The hook acts only on uncalled
  weak-key genes with |P − ½| ≥ τ. Delivered magnitudes on these routes are a quarter or less of the real |l|. This is
  the mechanism behind 630 × KGID compounding.
- **Side effect:** flipping block-prior genes DOWN → UP EXPANDS the call set; UP → DOWN does not shrink it.
  - Octile 5 (793 genes DOWN → UP): +281 calls per target (median +189, max +1,449), of which only +32 are octile-5 genes.
    It is a BH / UP-design cascade.
  - Octiles 0 and 2: +6 and +5.
  - Any prior change with many DOWN → UP flips (ID05: 149–302 per context) carries a call-count side effect that moves fid
    and jac independently of sign quality (sent to [d4cd5f] for ID05's reading).
- Off-octile delivered-sign spillover is small: 500–950 sign changes against 11k–31k on the octile.
- **Mechanism of the call expansion (14:15Z): CPM coupling from the block's UP mass.** Block genes are called at 98.5 %
  either way. The new calls are NON-block, non-head genes, whose call rate nearly doubles under the octile-5 flip:

  | non-block genes, CLEAN4 prior | call rate, none | call rate, octile 5 flipped |
  |---|---:|---:|
  | DOWN | 0.068 | 0.131 |
  | UP | 0.040 | 0.064 |
  | octile 5 | 0.056 | 0.116 |

  Making ~250 block genes UP adds carrier mass. With every library fixed, that compresses every other gene's CPM in the
  carriers, and 400 cells against 18,400 controls detect the shift. This is 610's CPM coupling, here driven by a
  PRIOR-ONLY change. Any prior that moves block mass UP or DOWN changes the call set and the null genes' signs through
  this route.

### 625e (cont.) — ROUTE-MASK sizing: price a sign source offline (14:20Z)
Delivery is explained by WHERE the source's changed genes sit. Define the routes on the base run's realised calls:
- block = called genes outside the head (head = 400 by |key|);
- hook = uncalled, weak-key half, |P − ½| ≥ τ.

Then compare the source's gate-sign change restricted to those routes with what is delivered:

| octile flip | source Δ (all gate) | route-restricted Δ (block / hook) | measured delivered Δ | overall w/o reach |
|---|---:|---|---:|---:|
| 0, f 1 | −0.0151 | −0.0090 (−0.0029 / −0.0061) | −0.0067 | −0.0019 |
| 2, f 1 | −0.0076 | −0.0014 (−0.0014 / 0) | −0.0010 | −0.0004 |
| 5, f 1 | −0.0062 | −0.0040 (−0.0040 / 0) | −0.0038 | −0.0005 |
| 5, f 0.35 | −0.0032 | −0.0040 | −0.0043 | −0.0037 (call-set side effect) |

- Delivered ≈ 0.7–1.0 × the route-restricted change.
- Overall ≈ **0.13–0.26 × the route-restricted gate-sign change** (≈ 0.2), excluding the call-expansion outlier.
- **Rule of thumb for sizing a sign source:** compute its gate-sign change on route genes only, then multiply by ~0.2.
  Seeds and the oracle (improvement-direction) arms of 625f will tighten the factor.

**Applied to the identity rider on H1** (ID05-like: rider signs in the BLOCK prior only, I vs A/B/C Flex):
- λ 0.5: source +0.0048 (f 1) / +0.0100 (f 0.35) on all gate genes, but only +0.0015 / +0.0030 on block genes. Predicted
  ≈ **+0.0003 … +0.0006 overall**.
- λ 1: +0.0005 … +0.0010.
- Both are below ID05's TRANSFERS bar (+0.0015).
- Most of identity's source gain sits on uncalled genes the current hook does not reach (strong-key or weak-prior genes,
  hook-route share +0.0001) and on head genes. That is the coverage trade-off 622o measured for idall3.

### 625d/e update, all four octiles at seed 0 (14:45Z; octile 7 landed 14:39Z)
Octile 7 is P > 0.583, the strongest DOWN prior, 9 % of gate pairs.
- Flipping it: f 1: fid −0.0051, jac −0.0033, nmae +0.0042 raw. Overall w/o reach −0.0055 (f 1) / −0.0117 (f 0.35).
- Delivery ratio 0.66–0.68.
- **Calls per target 2,571 → 3,857 (+1,287): the CPM-coupling call expansion at full strength.**

Decomposition into overall-equivalent parts (INV: fid 3.33, jac 2.72, nmae 1.65; ÷6):

| octile flip | Δ calls/target | nmae part | fid part | jac part | source Δ gate-sign | nmae-only per +0.01 | nmae+fid per +0.01 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0 (UP→DOWN), f 1 | +5 | −0.0009 | −0.0010 | 0 | −0.0151 | +0.0006 | +0.0013 |
| 2 (UP→DOWN), f 1 | +4 | −0.0004 | +0.0001 | 0 | −0.0076 | +0.0005 | +0.0005 |
| 5 (DOWN→UP), f 1 | +281 | +0.0003 | 0 | −0.0008 | −0.0062 | −0.0004 | −0.0004 |
| 7 (DOWN→UP), f 1 | **+1,287** | −0.0012 | −0.0028 | −0.0015 | −0.0243 | +0.0005 | +0.0016 |
| 7, f 0.35 | +1,287 | −0.0012 | −0.0074 | −0.0032 | −0.0282 | +0.0004 | +0.0031 |

**The exchange rate splits in two:**
1. **The sign-quality term.** With the call set stable (UP→DOWN flips), +0.01 of source gate-sign is worth ≈ +0.0005
   overall through nmae alone, and ≈ +0.0005 … +0.0013 through nmae + fid. The nmae-only term is the same on octiles
   0, 2 and 7 (+0.0004 … +0.0006), so it is the robust number. FINDINGS §0's "+0.0011 per +0.01" is this term at the
   upper end.
2. **The call-set term.** A prior change that moves block mass DOWN → UP expands the call set through CPM coupling
   (+281 calls for octile 5's weak DOWN priors, +1,287 for octile 7's strong ones). That costs jac and fid on top of any
   sign effect: −0.004 overall at f 1 for octile 7. UP → DOWN changes leave the call set unchanged. Size this term
   separately for any prior change; it scales with the strength of the flipped DOWN priors.

### 625f — the exchange rate is SYMMETRIC: an oracle IMPROVEMENT pays at the same nmae rate (14:47Z; `oracle625.py`, `delivery625.py`)
Arm `oracle0`: flip only the 152 octile-0 genes where H1's own in-sample programme (P_down over the 521 f 1 pools,
≥ 3 pools) disagrees with CLEAN4. That flip IMPROVES gate signs. It is leaky by design and calibrates delivery only.
The flip set was checked identical to the runner's (152 genes), and the bulk sums are bit-identical to `none`.

| arm, seed 0 | source Δ gate-sign | delivered Δ | ratio | Δ calls/target | overall parts nmae / fid / jac | nmae-only per +0.01 | nmae+fid per +0.01 |
|---|---:|---:|---:|---:|---|---:|---:|
| oracle0, f 1 | **+0.0068** | +0.0040 | 0.59 | +1 | +0.0003 / +0.0007 / 0 | **+0.0005** | +0.0016 |
| oracle0, f 0.35 | +0.0049 | +0.0028 | 0.57 | +1 | +0.0003 / +0.0008 / 0 | +0.0006 | +0.0022 |
| oct0 (degrade), f 1 | −0.0151 | −0.0067 | 0.44 | +5 | −0.0009 / −0.0010 / 0 | +0.0006 | +0.0013 |

- **Symmetric on nmae:** +0.0005 … +0.0006 overall per +0.01 of source gate-sign, whether signs get better or worse.
  The sign-quality exchange is a property of the delivery routes, not of the direction of the change.
- With fid included, the improvement direction pays a little more (+0.0016 … +0.0022 against +0.0012 … +0.0013).
  Oracle flips concentrate on called genes, whose sign fid reads directly.
- **Sizing rule, final form:** a sign source that leaves the call set unchanged is worth ≈ +0.0005 (nmae) to
  ≈ +0.0015 (nmae + fid) overall per +0.01 of SOURCE gate-sign accuracy, through today's block prior + hook. Add the
  call-set term separately for DOWN → UP block-mass shifts. Octile 5's oracle arm (298 genes DOWN → UP) prices that
  case next.

### Readout arm artifacts ready (2026-10-02 ~15:00Z)
- [4d0207] built the 16 single-octile arms, λ 0 and λ 0.5 × octiles 0–7. Manifest:
  `experiments/604-identity-prior/work/flip625_artifacts.json`, at 85d4be68 (module 428f4a32).
- Two contamination routes at λ 0.5 were found and fixed before any use:
  1. zc's scale had been recomputed from the flipped prior (4d0207's catch).
  2. v3's rank-based confidence reshuffled pdown magnitudes on ~5k non-octile genes per context (7772dc's catch). The fix
     takes the base artifact and re-decides only the masked octiles' sides; every gene keeps the base's |P − ½|.
- Independent check (7772dc): all 16 shas match the manifest; outside the octile, pdown and psign are byte-identical to
  the base; the kept |P − ½| matches to ≤ 5.6e-17 (1 ULP); octile sides agree with pdown.
- The arms are therefore valid for any hook τ the 10-08 base carries. The D/E/F day-of path is 604f (8465eb14):
  `DEF_DAY_DE_PRIOR_FLIP`.

### 625f, octile 5 oracle (14:57Z): an IMPROVING DOWN → UP flip of 298 genes leaves the call set unchanged
| arm, seed 0 | source Δ | delivered Δ | ratio | Δ calls/target | parts nmae / fid / jac | nmae-only per +0.01 | nmae+fid per +0.01 |
|---|---:|---:|---:|---:|---|---:|---:|
| oracle5, f 1 | **+0.0104** | +0.0030 | 0.29 | **+2** | +0.0008 / +0.0012 / −0.0001 | +0.0007 | +0.0019 |
| oracle5, f 0.35 | +0.0095 | +0.0029 | 0.30 | +2 | +0.0008 / +0.0012 / −0.0001 | +0.0008 | +0.0021 |
| oct5 (whole octile, degrade), f 1 | −0.0062 | −0.0038 | 0.61 | +281 | +0.0003 / 0 / −0.0008 | −0.0004 | −0.0004 |

- **The sign-quality exchange holds in the improvement direction for a DOWN-prior octile too:** nmae-only +0.0007 …
  +0.0008 per +0.01 of source gate-sign; nmae + fid +0.0019 … +0.0021. Overall gain of this oracle arm: **+0.0019 / +0.0020**.
- **The call-set term is NONLINEAR.**
  - The 298 oracle genes sit in the realised block about as often as the whole octile does (35 % vs 32 % of target-gene
    pairs). Yet they add +2 calls per target, against +281 for the whole octile's 793 genes and +1,287 for octile 7's 794
    strong-prior genes.
  - So the CPM-coupling call expansion behaves like a threshold or cascade (BH's step-up is self-reinforcing once enough
    genes cross) rather than scaling with the number of DOWN → UP flips.
  - For identity's ID05 (149–302 DOWN → UP flips per context, offset by 314–448 UP → DOWN), expect it to stay below the
    cascade, but measure the realised calls (R2).
- **Pooled sizing rule** (5 arms × 2 regimes): nmae-only ≈ **+0.0004 … +0.0008** overall per +0.01 of source
  gate-sign; nmae + fid ≈ **+0.0012 … +0.0022**. The rate is symmetric in direction. The call-set term is added
  separately and nonlinearly.

## 625d/e/f FINAL (2026-10-02 ~16:05Z): 3 seeds for octiles 0 and 5, all four oracle arms, KA s0 + s1
**KA:** `none` s0 AND s1 == 612's c4kg s0 / s1 bit for bit (q, p, lfc, sums): PASS.

**Realisation noise is negligible.**
- Octile 0's Δnmae over seeds 0/1/2: +0.00322 / +0.00327 / +0.00328 raw (sd 3.1e-5).
- Octile 5's: −0.00091 / −0.00090 / −0.00083 (sd 3.9e-5).
- Overall w/o reach: sd 0.00002 (octile 0) and 0.0001 (octile 5).
- On the board's scale, σ_k = sd × √(85/900) ≈ **1e-5 raw nmae**, two orders below the effects. Amendment 3's 2σ term is
  therefore ≈ 2e-5. The 0.25 × median bar binds, and the A/B/C read's realisation noise is not a concern.
- Reach does vary by seed (octile 0, f 0.35: overall sd 0.0004 from reach). Reach stays excluded.

**The exchange rate, final: CORRECTION to the "+0.0011 per +0.01" quoted at 13:50–14:45Z.**
That figure came from three degradation flips. The four oracle arms show that IMPROVEMENTS pay more, because the genes
where CLEAN4 is wrong concentrate on called and strongly responding genes, which fid reads.

| arms (seed 0, f 1 / f 0.35) | overall w/o reach per +0.01 of SOURCE gate-sign | per +0.01 of ROUTE-restricted change |
|---|---|---|
| oracle0 / 2 / 5 / 7 (improve) | nmae+fid: +0.0016 / +0.0018 / +0.0019 / +0.0027 (f 1) | **median +0.0035** (range +0.0014 … +0.0050) |
| oct0 / 2 / 7 (degrade, call set ~stable or excluded) | nmae+fid: +0.0013 / +0.0005 / +0.0016 | +0.0010 … +0.0027 |
| pooled, all 8 arms | +0.0016 (f 1) / +0.0027 (f 0.35) | |

- **Sizing rule for a sign source (an improvement): ≈ +0.002 overall per +0.01 of source gate-sign** (nmae + fid, w/o
  reach; range +0.0015 … +0.0027). Equivalently ≈ 0.35 × its gate-sign change on route genes (realised block +
  hookable weak-key uncalled).
- The nmae-only part is ≈ +0.0003 … +0.0014.
- **So the programme prize (+0.045 clean-oracle gate sign) is ≈ +0.007 … +0.012 overall through today's routes.** The
  council's original +0.008 … +0.010 stands. The "+0.005" given at 13:50Z is withdrawn.
- In-sample ceiling, for scale: each octile's oracle arm gains +0.0010 … +0.0022 overall on H1. Summed over 8 octiles,
  if additive, that is ≈ +0.01. This is a leaky in-sample programme, so it is an upper bound.
- **Revised predictions:**
  - 625's A/B/C → D/E/F mask, if GO: Broad-16's +0.0165 gate sign × 0.2 ≈ **+0.003** (Amendment 3 said +0.002).
  - ID05-like identity block prior on H1, route-restricted (block) change × 0.35: **λ 0.5 +0.0005 / +0.0011 (f 1 /
    f 0.35); λ 1 +0.0009 / +0.0017.** λ 0.5 stays below the +0.0015 TRANSFERS bar. λ 1 reaches it only at f 0.35.
- **The call-set term** (DOWN → UP block-mass shifts, through CPM coupling) is nonlinear and separate:
  - oracle7: +47 calls;
  - oracle5: +2;
  - whole octile 5: +281;
  - whole octile 7: +1,287.
