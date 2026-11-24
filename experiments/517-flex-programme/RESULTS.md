# 517 — Target-independent DE sign information: the lab/Flex programme, ready for the 09-30 read

Worker 517 (Claude Opus 5.5), mandate `MANDATE.md`. Analysis only: no build, no `vcc`, no `cell-eval2`, no git, nothing
written outside this directory. **Status (2026-09-28 15:40Z): §1–§3 done on the analytic instruments (validated against
real emitter runs, §1) and public-data instruments. Real-emitter confirmation runs of the hook arms (`drive517.sh`, armed,
pid in `work/drive517.pid`) start after the build queue (≥ 02:00Z 09-29); their numbers go in §9 when they land.**

## 0. Read this before 09-30 00:30Z (pre-read note, written before any 09-30 entry exists)

1. **458d's jac falsifier (±0.0010 vs 446b) may fire from a known emitter leak, not a bug.** On H1, Arc's programme as the
   block prior makes the K-class emitter leak calls it cannot suppress: `call_set_unsuppressible` median **837 per target**
   (ungated, 462l) and 422 (gated33, 467) against **0** for CLEAN4 (435). That is H1's jac −0.0027 / −0.0022 raw. On A/B/C
   the 458d build leaks much less but not nothing: its meta records unsuppressible **median 1, max 1,655** (446b: median 0,
   max 86). If a tenth of A/B/C targets leak ~1,000 calls, jac falls ~0.002. So a jac read of −0.001 … −0.003 is the leak's
   expected size; diagnose it as the leak (compare `call_set_unsuppressible` of 458d/458g vs 446b), then read fid as
   preregistered. On H1 the leak also costs fid ~−0.0016 (leaked calls are coin flips; analytic fid without the leak
   +0.0213 vs the real run's +0.0198).
2. **Mapping the 458d read onto H1's attenuation curve** (§2.4): the preregistered POS / PARTIAL / NULL / HARM fid bins
   correspond to Arc's programme having lost roughly r ≤ 0.23 / 0.23–0.25 / 0.25–0.27 / ≥ 0.28 of its H1 sign information
   on A/B/C (random-flip emulation, board fid ≈ 0.67 × H1 fid as for 392). The block prior's bar is narrow; the hook's is not.
3. **The most valuable NEW arm in every branch except HARM is the complement hook sourced from Arc's programme** (§2.1): H1 nmae
   −0.0148 raw gated (10× 496's CLEAN4 hook, −0.0015), and it beats the CLEAN4 hook down to r ≈ 0.38, far past the point
   where the block prior stops beating CLEAN4 (r ≈ 0.24–0.26), because it replaces the design's ~0.53-accurate signs on uncalled
   genes, not CLEAN4's. A public within-lab analogue shows exactly this (§2.5: HEK293T's programme hurts HCT116's block
   and helps its uncalled genes). Patch ready (`patch_builder517.py`, tested, not applied), env knobs in §6.
4. **Use a better estimator for the hook's source: mean z over Arc's non-essential rows (`meanz_ne`, §2.6).** It wins on H1
   (nmae −0.0183 vs −0.0148, reach 0, gain spread over 47 of 58 targets) AND on the cross-line Xaira analogue with Arc-sized
   programmes (+0.02 … +0.03 over the P_down rule; the essential rows carry nothing for another line's non-essential genes).
   The A/B/C artifact is built (`work/arch1meanz_ne_hook_abc.npz`, `165c6b78…`); `patch_sameline517.py` gives October the same.
5. **Null / harm branch:** nothing target-independent beats CLEAN4 robustly. The recipient-side axes (identity, own-control
   PCs) win on H1 but REVERSE in 3 of 6 public cancer lines (§3); lab consensus and quality-weighted consensus do not beat
   CLEAN4. Keep CLEAN4 + 496a; the Arc (meanz_ne) hook stays a candidate because its bar is the design, not CLEAN4.
6. **q and dose (§2.2–2.3):** keep q80 with the gate; q70 is the ungated programme's reach contingency (H1 reach −0.041 →
   −0.010). w 1 wins at r < 0.2, w ≤ 0.5 insures r ≥ 0.25: 458e's read vs 458d's sets it.
7. **Flex technical components (§7):** library-size / composition effects are real (eps_g reliable, strong perturbations
   shrink libraries) but two orders of magnitude too small to carry signs (gate accuracy 0.52–0.53, and they hurt every
   programme they are added to). Per-gene excess dropout is a shared probe/state property that predicts DE CALLS, not signs
   (a call-set lead for another mandate).

## 1. Known answers (every screen carries one; all pass)

| check | reproduced | recorded |
|---|---|---|
| 490's recon, every base run × seed (kf, mde_fu/fd, n_up/dn, flips, n_block == diag), 435 c4 / 467 gated33 / 462l arcfull | exact (worst diff 0.0) | — |
| delivered tables → metrics json (score_de, cell-eval2's own DE functions) | ≤ 1e-9 | 435 / 467 / 462l metrics |
| MDE with the shipped signs == recon's A_int; "no change" variant Δ | ≤ 1e-9; Δ ≤ 1e-17 | — |
| **496's hook, CLEAN4 τ 0.05, seed 0: nmae / reach / flips / fd / ov [CI]** | **−0.00147 / +0.0035 / 810 / 0.771 / +0.0011 [−0.0001, +0.0027]** | 496 real run −0.00147; 478 analytic +0.0011 [−0.0001, +0.0027], 810, 0.77 |
| 462f fid-like: CLEAN4 / Arc-full (`acc517.py`) | **0.5244 / 0.5545** | 0.5244 / 0.5545 |
| 450j gate-block CLEAN4 | 0.5799 | 0.5799 |
| halves: Arc programme vs all controls == the h1_2025_full bank's (215 rows) | max \|ΔP\| 0, signs 1.000000 | — |
| halves: h1-test truth recomputed vs all controls == the reference lfc | max 2.6e-15 | — |

**Validation of the analytic block-prior path against real emitter runs** (`blockval517.py`: from the c4 base, swap the
block prior; compare with the real runs that did exactly that; 3 seeds each):

| arm | nmae pred / real | fid pred / real | reach pred / real | jac real (not seen analytically) |
|---|---|---|---|---|
| Arc ungated (462l) | −0.01484 / −0.01479 | +0.0213 / +0.0198 | −0.0405 / −0.0388 | −0.0026 |
| gated33 (467) | −0.01026 / −0.01020 | +0.0167 / +0.0161 | +0.0018 / +0.0015 | −0.0022 |
| gated50 | −0.00765 / −0.00760 | +0.0119 / +0.0126 | +0.0031 / +0.0029 | −0.0012 |
| gated67 | −0.00416 / −0.00410 | +0.0069 / +0.0067 | +0.0018 / +0.0017 | −0.0002 |
| wdown (470) | −0.01065 / −0.01059 | +0.0163 / +0.0157 | +0.0007 / −0.0003 | −0.0022 |

nmae reproduces to 5e-5, fid to ~10 % (the gap is the leak), reach to ~0.002; jac is the leak and needs real runs.

## 2. §1 — the positive / partial branch: exploitation arms (H1, analytic, 3 seeds; Δ raw vs the base; overall-equiv with
467's spans; [ ] = seed-0 target-bootstrap 95 %; "ne" = the 67 non-essential targets, the panel-like subset)

### 2.1 The complement hook sourced from Arc's programme (`price517.py`, `hookgrid517.py`)
Arc = P_down over `h1_2025_full` minus h1-test's 85 targets' rows (467's rule). "gated" = the hook follows the block gate:
467's strong 57 targets take Arc's sign, the weak 28 (strength ≤ 1, = `le1` on H1) take CLEAN4 τ 0.05.

| base (board analogue) | hook | nmae | fid | reach | overall-equiv [CI] | ne ov | flips/target |
|---|---|---:|---:|---:|---|---:|---:|
| c4 (446bj2z) | CLEAN4 τ .05 (496a) | −0.00146 | 0 | +0.0029 | +0.0010 [−0.0001, +0.0027] | +0.0008 | 810 |
| c4 | Arc τ .05, ungated | −0.0167 | 0 | −0.0175 | +0.0014 [−0.0044, +0.0054] | −0.0004 | 2,586 |
| c4 | **Arc τ .05, gated** | **−0.0148** | 0 | −0.0031 | **+0.0036 [+0.0014, +0.0054]** | +0.0025 | 2,562 |
| c4 | Arc τ .10, gated | −0.0106 | 0 | 0 | +0.0030 [+0.0018, +0.0040] | +0.0023 | 1,323 |
| c4 | Arc τ .15, ungated | −0.0065 | 0 | 0 | +0.0018 [+0.0012, +0.0024] | +0.0016 | 568 |
| g33 (458g) | CLEAN4 τ .05 | −0.0015 | +0.0007 | +0.0035 | +0.0014 [+0.0002, +0.0032] | +0.0012 | 810 |
| g33 (458g) | **Arc τ .05, gated** | **−0.0157** | +0.0044* | −0.0031 | **+0.0062 [+0.0034, +0.0092]** | +0.0043 | 2,562 |
| g33 | Arc τ .10, gated | −0.0112 | +0.0034* | 0 | +0.0050 [+0.0033, +0.0068] | +0.0037 | 1,323 |
| arc (458d) | **Arc τ .05, ungated** | **−0.0177** | +0.0057* | −0.0091 | **+0.0064 [+0.0024, +0.0100]** | +0.0044 | 2,586 |
| arc (458d) | Arc τ .05, gated | −0.0156 | +0.0047* | −0.0031 | +0.0064 [+0.0035, +0.0093] | +0.0045 | 2,562 |
| c4 | **meanz_ne (§2.6), gated** | **−0.0183** | 0 | 0 | **+0.0051 [+0.0035, +0.0067]** | — | ~2,560 |
| g33 (458g) | **meanz_ne, gated** | **−0.0182** | +0.0038* | 0 | **+0.0072 [+0.0048, +0.0101]** | — | ~2,560 |
| arc (458d) | meanz_ne, ungated / gated | −0.0201 / −0.0182 | +0.0050* / +0.0042* | −0.0061 / 0 | +0.0072 / +0.0074 | — | — |

\* fid on the g33/arc bases comes from flipping LEAKED complement calls (§0.1: 422–837 per target on H1); A/B/C's builds leak
~1 per target, so expect ~0 fid there. The nmae (the hook's real effect) is the same on every base.
- Source accuracy on uncalled gate genes (the hook's territory): Arc τ .05 **0.685** (25,471 pairs) vs the delivered design
  signs 0.532 on the same genes; CLEAN4 τ .05 0.604 on 6,352 pairs.
- **Not inflated by shared controls** (`halves517.py`, 462's rule): Arc's accuracy on h1-test's gate genes with the programme
  and the truth taken against DISJOINT control halves equals the shared-half value (gate-block 0.6556 vs 0.6556; gate
  non-block 0.6237 vs 0.6238; non-block τ .05 0.6694 vs 0.6695). Only non-significant genes inflate (fid-like block +0.0071,
  i.e. about a quarter of Arc's H1 fid-like edge over CLEAN4 is control noise, which matters for the fid mapping in §2.4).
- Seeds: gated τ .05 on c4: nmae −0.01475 / −0.01481 / −0.01494; ov +0.00353 / +0.00355 / +0.00359.

### 2.2 Where q80 should sit when the prior is Arc's (`grid517.py`, block only, w 1)
| protect_q (share of the block NOT protected) | gated: fid / nmae / reach / ov | ungated: fid / nmae / reach / ov |
|---|---|---|
| 0.5 | +0.0068 / −0.0062 / −0.0024 / +0.0051 | +0.0095 / −0.0076 / −0.0143 / +0.0048 |
| 0.6 | +0.0103 / −0.0069 / −0.0024 / +0.0072 | +0.0136 / −0.0088 / −0.0143 / +0.0074 |
| 0.7 | +0.0136 / −0.0085 / +0.0020 / +0.0103 | +0.0173 / −0.0112 / **−0.0098** / **+0.0109** |
| **0.8 (shipped)** | +0.0167 / −0.0103 / +0.0018 / **+0.0125** | +0.0213 / −0.0148 / −0.0405 / +0.0085 |
| 0.9 | +0.0182 / −0.0119 / +0.0014 / +0.0137 | +0.0236 / −0.0168 / −0.0401 / +0.0104 |
(at full H1 quality; vs CLEAN4 q0.8). At attenuation r 0.1 / 0.2 / 0.3 the gated q0.9 − q0.8 difference is +0.0006 /
+0.0004 / +0.0001, and q0.7 is below q0.8 for every r ≤ 0.3 (§2.4's grid).
- **Gated (458g's construction): keep q80.** q90 buys +0.0012 at full same-line quality and ≤ +0.0006 at any plausible
  cross-line quality: inside the board's noise. Not worth a slot.
- **Ungated (458d's construction): q70 is the reach contingency's setting.** It moves the four small-gate targets' reach-
  critical genes under protection: reach −0.0405 → −0.0098, fid −0.0040, nmae +0.0036, ov +0.0024 on H1. If 458d reads
  fid-positive with reach ≤ −0.005 and the gate (458g) reads null or harm, `VCC_DE_SIGN_PROTECT_Q=0.7` on 458d is the
  REACH-CONTINGENCY arm (a board number, since sign-prior reach does not price on H1).

### 2.3 Programme-strength dose (block only, q0.8; w = weight of z(Arc) in sign(w z(Arc) + (1−w) z(CLEAN4)); 458e is w 0.5)
| r (Arc sign info lost) | gated w 1 / 0.75 / 0.5 / 0.25 | ungated w 1 / 0.75 / 0.5 / 0.25 |
|---|---|---|
| 0.0 | +0.0125 / +0.0113 / +0.0098 / +0.0065 | +0.0085 / +0.0082 / +0.0060 / +0.0017 |
| 0.1 | +0.0076 / +0.0069 / +0.0062 / +0.0048 | +0.0036 / +0.0038 / +0.0023 / −0.0000 |
| 0.2 | +0.0015 / +0.0011 / +0.0016 / +0.0017 | −0.0015 / −0.0007 / −0.0009 / −0.0011 |
| 0.3 | −0.0026 / −0.0028 / −0.0010 / +0.0000 | −0.0063 / −0.0051 / −0.0037 / −0.0027 |
| 0.4 | −0.0085 / −0.0087 / −0.0049 / −0.0021 | −0.0114 / −0.0101 / −0.0057 / −0.0050 |
- w 1 wins while the programme keeps ≥ 80 % of its H1 information; around r 0.2 the dose curve is flat; beyond, the blend
  insures the downside (w 0.25 costs nothing at r 0.3 where w 1 costs −0.0026 gated). So 458e's read against 458d's
  sets w: if Δfid(458e) ≥ Δfid(458d), A/B/C sits at r ≥ 0.2 and October should use w ≤ 0.5.

### 2.4 Attenuation: what each 458d read implies (block curves at q0.8, 3 seeds; Arc with a fixed random fraction r of genes' signs flipped)
| r | ungated block: fid / nmae / reach | gated block: fid / nmae / ov | gated block + Arc hook τ .05: ov (hook increment) | hook alone on c4, gated τ .05: nmae / ov |
|---|---|---|---|---|
| 0.0 | +0.0213 / −0.0148 / −0.0405 | +0.0167 / −0.0103 / +0.0125 | +0.0164 (+0.0039) | −0.0148 / +0.0036 |
| 0.1 | +0.0129 / −0.0105 / −0.0350 | +0.0097 / −0.0070 / +0.0076 | +0.0106 (+0.0029) | −0.0114 / +0.0026 |
| 0.2 | +0.0046 / −0.0063 / −0.0312 | +0.0033 / −0.0035 / +0.0015 | +0.0032 (+0.0017) | −0.0072 / +0.0015 |
| 0.3 | −0.0027 / −0.0031 / −0.0305 | −0.0027 / −0.0011 / −0.0026 | −0.0017 (+0.0008) | −0.0045 / +0.0007 |
| 0.4 | −0.0097 / +0.0001 / −0.0325 | −0.0085 / +0.0019 / −0.0085 | −0.0088 (−0.0004) | −0.0007 / −0.0004 |
- Board 458d Δfid ≈ 0.67 × H1 Δfid (392's same-quality transfer) gives the r column for each preregistered bin (§0.2). The
  H1 fid at r 0 is ~25 % control-noise inflated (§2.1), which shifts every r by ≈ −0.03: the mapping is coarse.
- The block's nmae stays helpful to r ≈ 0.35, beyond its fid (part of it is the prior's less down-skewed kf balance, a
  magnitude effect), so a 458d read with fid null and nmae ≤ −0.0015 is expected, not contradictory.
- **The hook's bar:** vs 496's CLEAN4 hook (nmae −0.0015), the Arc hook wins on nmae until r ≈ 0.38; the gated block stops beating CLEAN4 at r ≈ 0.24 (ungated fid: 0.26).

### 2.5 A cross-line, same-lab analogue on public cancer banks (`labhook517.py`; no H1 in it)
The A/B/C situation is "design from other labs' donor rows, hook from a same-lab OTHER-line programme". Xaira (HCT116,
HEK293T) and Weissman (K562-GWPS; the RPE1/Jurkat/HepG2 essential screens) reproduce it: hold out line L, design_t = the
other lab's donor rows for t, calls = top 4 % by |design| + top 22 % by L's control expression (production proportions),
truth = L's delta signs on |delta| ≥ 0.2 pairs, non-essential targets for the genome-wide screens (3,000 per line).
Sign accuracy; "@63 %" = the source's most confident 63 % of genes (Arc's τ .05 coverage on H1):

| held-out L (same-lab partner) | complement (hook's territory): design / same-lab @100 % / @63 % / other-lab consensus @63 % | block: design / same-lab @100 % / other-lab @100 % |
|---|---|---|
| HCT116 (HEK293T) | 0.506 / **0.542 / 0.561** / 0.534 | 0.523 / **0.488** / 0.557 |
| HEK293T (HCT116) | 0.510 / **0.591 / 0.622** / 0.504 | 0.518 / **0.619** / 0.488 |
| K562 (RPE1, an essential screen) | 0.519 / 0.537 / 0.557 / 0.567 | 0.554 / 0.573 / 0.555 |
| RPE1, Jurkat, HepG2 (K562-GWPS) | 0.587 / 0.633 / 0.686; 0.583 / 0.573 / 0.603; 0.588 / 0.592 / 0.629 | (essential targets) |

- In the hook's territory a same-lab other-line programme beats the other-lab design in every genome-wide case: +0.036 /
  +0.081 / +0.019 with every gene, **+0.055 / +0.112 / +0.038 at 63 % coverage** (H1, same line: +0.153). It beats the
  other-lab consensus (CLEAN4's role) in 2 of 3.
- **HEK293T → HCT116 is the null-block, positive-hook case, observed:** HEK293T's programme HURTS HCT116's block (0.488 vs
  the design's 0.523 and the other lab's 0.557) and still HELPS HCT116's uncalled genes (0.561 vs 0.506 at 63 %).
- No expression-similarity filter helps (`simexpr517.py`): stratified by |Δ ln control expression| between the two lines,
  the partner programme's complement accuracy is flat (HCT116 0.554 / 0.556 / 0.552; HEK293T 0.598 / 0.598 / 0.614 from
  similar to dissimilar genes). The programme transfers as a whole, not on "shared-expression" genes.
- On H1 the hook's nmae is linear in this gain: −0.097 per unit accuracy (0.153 → −0.0148; 0.046 → −0.0045). The
  analogues' gains map to H1-equivalent nmae −0.0037 … −0.011 (−0.0053 for the HCT116-like case), against the CLEAN4
  hook's −0.0015; in grid terms r ≈ 0.1–0.3.

### 2.6 A better programme ESTIMATOR for the hook: mean z over the non-essential rows (`progest517.py`, `variant517.py`, `variantblock517.py`)
**Cross-line, within a lab (Xaira, Arc-sized programmes: 275 rows, 36 % essential = Arc's H1 mix; 5 draws; non-essential
targets' uncalled |delta| ≥ 0.2 genes at 63 % coverage; the other-lab design ≈ 0.51):**

| programme from the partner line | → HCT116 | → HEK293T |
|---|---:|---:|
| P_down, all ~17k rows (upper reference) | 0.551 | 0.600 |
| P_down, 275 random rows (~7 % essential) | 0.541 | 0.532 (sd 0.028) |
| **P_down, Arc-like mix (99 essential + 176 non-essential): the production rule on Arc's composition** | **0.528** | **0.544** |
| … its 99 essential rows only | 0.500 | 0.480 |
| … its 176 non-essential rows only | 0.536 | 0.546 |
| … mean z over the mix | 0.539 | 0.520 |
| **… mean z over the 176 non-essential rows (meanz_ne)** | **0.548** | **0.573** |

Essential rows carry nothing (or the opposite) for another line's non-essential responses; mean z (magnitude-weighted)
beats the sign count. **On H1 (same line) the same estimator also wins** (hook on the CLEAN4 block, gated, Arc τ .05
coverage, 3 seeds): meanz_ne **nmae −0.0183, reach 0, overall-equiv +0.0051 [+0.0035, +0.0067]** vs the production P_down
−0.0148 / −0.0031 / +0.0036; complement-gate accuracy 0.702 vs 0.686. (P_down over the non-essential rows alone is worse on
H1, −0.0136, so it is the mean z that matters there.) As an UNGATED block prior, meanz_ne on H1 gives fid +0.0178, nmae
−0.0160 and **reach +0.0171 instead of 458d's −0.0405** (overall-equiv +0.0175 vs +0.0085; gated +0.0155 vs 458g's
+0.0125). Reach is a small-gate lottery for sign priors on H1 and the estimator was chosen after seeing it, so the block use
is a candidate, not a claim.
- **Broad, not a few targets (`pertarget517.py`, seed 0, the 58 nmae-scored targets):** meanz_ne improves 47 of 58 (median
  −0.0149, the top five carry 27 % of the gain; non-essential mean −0.0160, essential −0.0238); Arc's P_down 44 of 58 (top
  five 35 %); 496's CLEAN4 hook 41 of 58 with 81 % of its gain in five targets (non-essential mean −0.0005).
- **Attenuated (`mzgrid517.py`, the same random-flip emulation, 3 seeds):** the meanz_ne hook alone (gated, CLEAN4 block)
  gives nmae −0.0183 / −0.0133 / −0.0091 / −0.0054 / −0.0020 at r 0 / 0.1 / 0.2 / 0.3 / 0.4, ~1.2× the P_down hook at every r;
  it beats the CLEAN4 hook (−0.0015) to r ≈ 0.42. Its reach becomes −0.009 once r ≥ 0.2 (a small-gate term). Block + hook
  from the one meanz_ne artifact (gated): overall-equiv +0.0202 / +0.0139 / +0.0044 / −0.0009 / −0.0077; ungated, the block's
  full-quality reach gain (+0.017) turns into −0.039 by r 0.2, so **the gate stays** for any meanz_ne block use.
- **No CLEAN4 blend needed (`mzblend517.py`):** a hook source w·z(meanz_ne) + (1−w)·z(CLEAN4) at the same coverage gives
  nmae −0.0150 / −0.0083 / −0.0052 / −0.0022 (w 0.5) against −0.0183 / −0.0091 / −0.0054 / −0.0020 (w 1) at r 0 / 0.2 / 0.3 /
  0.4: the pure source is as good or better on nmae at every quality; the blend only moves the small-gate reach term.
- **Artifact (hook source for A/B/C):** `work/arch1meanz_ne_hook_abc.npz`, sha `165c6b78a16de2990ed87ca0e3436ab8c4023f7dddcdb30cf00b8f5eebda9f6b`
  (`artifact517.py`; 275 rows after the A/B/C panel, 197 non-essential; pdown encoded so that τ .05 covers exactly the
  production programme's 11,313 official-axis genes; on each context's gated genes 6,271–6,562 signs; sign agreement with
  458d's programme 0.81). KA: `artifact517.py --estimator pdown_all` rebuilds 458d's pdown_A/B/C exactly.
- **October:** `patch_sameline517.py` (NOT applied) adds `def_sameline.py --programme-estimator {pdown,meanz_ne}`; default
  pdown rebuilds 458d's programme exactly, meanz_ne rebuilds the artifact above exactly (both KAs on a patched copy, A/B/C
  bundle, `--programme-source h1_2025_full`). `patch_defday517.py` passes it as `DEF_DAY_HOOK_ESTIMATOR=meanz_ne`.

## 3. §2 — the null / harm branch: other target-independent sources, ranked (H1, analytic, 3 seeds, vs the incumbent's c4
base; `rank517.py`, `acc517.py`; the H1 identity axis is H1's controls vs A/B/C's in the same chemistry; "own" is 491c's
z(CLEAN4) + z(proj_k2 of CLEAN4 on H1's control PCs) + z(−identity), PCs by a streamed randomized SVD)

| source | fid-like / gate-block / gate-nonblock acc | block prior ungated: ov [CI] (nmae, fid) | gated (strong 2/3): ov [CI] | hook, CLEAN4-matched coverage: ov (nmae) |
|---|---|---|---|---|
| CLEAN4 (incumbent) | 0.5244 / 0.5799 / 0.5674 | 0 | 0 | +0.0009 (−0.0015) |
| **own (491c)** | **0.5303 / 0.5926 / 0.5656** | +0.0017 [−0.0034, +0.0053] (−0.0020, +0.0054) | **+0.0026 [+0.0011, +0.0042]** (−0.0017, +0.0044) | +0.0007 (**−0.0042**) |
| CLEAN4 + 2·identity | 0.5273 / 0.5905 / 0.5639 | **+0.0033 [+0.0006, +0.0076]** (−0.0015, +0.0026) | +0.0022 [−0.0006, +0.0059] | +0.0004 (−0.0011) |
| identity alone (−I) | 0.5247 / 0.5848 / 0.5588 | +0.0026 [−0.0015, +0.0076] (−0.0019, +0.0007) | +0.0009 | −0.0007 (−0.0002) |
| CLEAN4 + identity | 0.5258 / 0.5872 / 0.5635 | +0.0009 | +0.0001 | −0.0005 |
| Weissman + Xaira (lab units, 450j) | 0.5257 / 0.5860 / 0.5491 | −0.0009 | +0.0007 | −0.0020 (+0.0026) |
| + Satija (VIPerturb, Flex-style) | 0.5248 / 0.5818 / 0.5444 | −0.0021 | −0.0004 | −0.0017 |
| six public lines, equal weights (`wcons517.py`) | 0.5261 / 0.5858 / 0.5514 | −0.0013 (−0.0012, +0.0020; reach −0.015) | — | −0.0020 (+0.0010) |
| six lines weighted by measured cross-line quality (xline517: HCT116 .60, HepG2 .59, K562 .58, RPE1 .58, Jurkat .54, HEK293T .52 → w = acc − .5) | 0.5272 / 0.5889 / 0.5519 | −0.0014 (−0.0012, +0.0029; reach −0.018) | — | −0.0019 (+0.0009) |
| all-down | 0.5175 / 0.5587 / 0.5080 | −0.0032 [−0.0053, −0.0012] | −0.0029 | −0.0042 (+0.0095) |
| (Arc, same line: reference) | 0.5545 / 0.6581 / 0.6308 | +0.0084 | +0.0125 | +0.0023 (−0.0090) |

- Only the recipient-side axes beat CLEAN4 on H1; every other lab's programme loses (as 450j/462j found at accuracy level).
  Weighting the six public lines by their measured cross-line quality helps over equal weights (gate-block 0.5889 vs
  0.5858; fid +0.0029 vs +0.0020 as a block prior) but still does not beat CLEAN4 overall (a small-gate reach term) and is
  worse than no hook as a hook source.
  The ranking by the robust members (nmae, fid): **own (gated) > CLEAN4 + 2·identity > CLEAN4**, each worth ≤ +0.003
  overall-equiv on H1, i.e. ≤ +0.001–0.002 on the board at the 0.07–0.67× sign-prior transfer. `own`'s k = 2 was chosen on
  H1 and H1's PC2 is spread over several PCs in B/C (491b); the identity axis has no free parameter.
- As hook sources, only `own` (nmae −0.0042) beats CLEAN4's (−0.0015) on nmae; its reach is a lottery term (−0.0025).
- **H1's recipient-side gain is honest on H1** (`halves517b.py`: the identity axis from one control half vs the truth from
  the other: gate-block 0.5841 disjoint = 0.5841 shared; c4id2 0.5909 = 0.5909) **but does not generalise across cancer
  lines** (`xline517.py`, six public lines, sign accuracy on (target, gene) pairs with |delta| ≥ 0.2, each line held out):

  | held-out line | consensus of the other lines | + the line's own −I (best λ) | −I alone |
  |---|---:|---:|---:|
  | K562 (non-essential targets) | 0.5549 | 0.5421 (λ 0.5, **−0.013**) | 0.4730 |
  | RPE1 | 0.7556 | 0.7445 (**−0.011**) | 0.4133 |
  | Jurkat (context A's line) | 0.5745 | 0.5720 (−0.003) | 0.5181 |
  | HepG2 | 0.6649 | 0.6578 (−0.007) | 0.4527 |
  | HCT116 (non-essential) | 0.5419 | 0.5448 (+0.003) | 0.5215 |
  | HEK293T (non-essential) | 0.5289 | 0.6334 (λ 2, **+0.105**) | 0.6360 |

  "Identity genes go down" holds for H1, HEK293T and weakly HCT116, and REVERSES in RPE1, HepG2 and K562. It is a property
  of some lines, not a law, and nothing on the controls says which kind a new line is. **So none of the recipient-side
  sources is a safe A/B/C or D/E/F prior**; under a null 458d read the ranking is: CLEAN4 block (board-validated, 392) +
  496a's CLEAN4 hook (10-05), then the Arc-sourced hook (§2.4: its bar is the design's uncalled-gene accuracy), and no
  recipient-side axis. `own`'s k = 2 PCs are H1-fitted as well (491b).
- CLEAN4-hook variants (`c4var517.py`, H1 analytic, 3 seeds; for the harm branch): τ .05 (496a) nmae −0.00146; the four
  banks' unanimous genes −0.00131; unanimous & τ .05 −0.00126 (accuracy 0.618); τ .03 −0.00214; τ .02 −0.00167. All within
  ±0.001 of 496a: nothing to change there.
- The same instrument re-measures 462's lab effect in the panel's regime (non-essential targets): a same-lab other line
  predicts HEK293T at 0.595 (HCT116) against 0.541 for the best other-lab line (RPE1) and 0.508 for the other-lab
  consensus; for HCT116, HEK293T (0.538) ties the best other-lab line (RPE1 0.542). Lines can even anti-predict each other
  (Jurkat → HEK293T 0.448). A same-lab programme helps, but by line-pair-dependent amounts: the 09-30 read is the only way
  to learn Arc's.

## 4. The October artifact (confirmed on disk, no new artifact needed)
- `scripts/def_sameline.py --programme-source h1_2025_full` writes `programme_<c>` AND `pdown_<c>` for every context
  (lines 251–253), with the D/E/F panel's own rows excluded (396's rule) and 481 F7's guard (≥ 150 rows, coverage ≥ 0.5).
  461dr/469dr rehearsed it on the relabelled bundle: `experiments/461-arcprog-october-path/work/run/sameline.npz` carries
  `pdown_D/E/F` (coverage 0.975, 56.8 % down, 275 rows); `def_day.sh` passes `DEF_DAY_PROGRAMME_SOURCE` and
  `DEF_DAY_PRIOR_GATE` through. Arc's programme uses no D/E/F cells; the gene axis is the official 18,533 (fixed across
  rounds), so it exists for any panel. Only risk: a D/E/F panel overlapping > 150 of H1's 300 perturbations (A/B/C: 25).
- The hook needs the builder patch (§5) and a def_day pass-through: **`patch_defday517.py`** (NOT applied) goes on top of
  502's `DEF_DAY_MDE_SIGN` (502 = 497's `e29e793a` / `60b2a777` + `patch_defday502.py` → `610db782` / `20142142`; 517 →
  `ff54380c` / `3bdee4ce`, copies in `work/patchtest/defday/`). New forms: `DEF_DAY_MDE_SIGN=programme:<tau>[:<tau_c4>]`
  (needs `DEF_DAY_PROGRAMME_SOURCE`; reads the day's `sameline.npz`) and `DEF_DAY_MDE_SIGN=hookbank:<bank>:<tau>[:<tau_c4>[:le<N>]]`
  (the hook alone on a CLEAN4 block: a second `def_sameline.py --programme-source <bank>` artifact, `hookprog.npz`,
  exported as `file:<npz>:<sha>:…`). `riders_match` checks every form against the final meta. Tests: `bash -n` OK;
  `test_defday517.py` 13/13; 502's own `test_riders502.py` 9/9 on the 517 copy (backward compatible).
- On A/B/C, Arc τ .05 signs 5,750–6,050 of each context's ~9,600–10,100 gated genes (CLEAN4 τ .05: ~2,500); where both
  are confident they agree 78–79 % (93 % at τ .1).

## 5. The patch (`patch_builder517.py`, NOT applied; `test_patch517.py` 26/26 PASS on copies)
On top of 496's builder patch (HEAD `b4975024` + `patch_builder496.py` = `24b95dec`; emitter `1658ae2c`, which already
reads a per-block `mde_add_sign`). Patched copy sha `affde905…` (`work/patchtest/bld517.py`). `VCC_DE_MDE_SIGN` grammar:
- `clean4:<tau>` — 496's, parsed identically (same config dict; tested against 496's patched builder);
- `programme:<tau>[:<tau_c4>]` — the VCC_SAMELINE artifact's `pdown_<c>` where the programme replaced CLEAN4; follows
  `VCC_DE_PRIOR_GATE` (weak targets → CLEAN4 at tau_c4, default 0.05); a programme-less context → CLEAN4 at tau_c4;
- `file:<npz>:<sha256>:<tau>[:<tau_c4>[:le<N>]]` — `pdown_<c>` from any same-line-format artifact WITHOUT switching the
  block prior (the hook alone on a CLEAN4 block); optional own strength gate `le<N>` (468's measure), refused together
  with `VCC_DE_PRIOR_GATE`.
- `sameline_load` additionally carries `pdown` (nothing else reads it); config `de_mde_sign` gains a per-context audit
  (source used, n_sign, n_down, weak targets) for the programme/file forms only.
- Known answers the coordinator must run after applying: KA-1 unset = the base's output sha; KA-496 `clean4:0.05` = 496a's
  output sha; KA-2 the base's VCC_SAMELINE(+gate) env without the hook = 458d/458g's output sha.

## 6. Ranked next arms per 09-30 outcome (bins: `458/PREREGISTRATION.md`, `468/PREREGISTRATION.md`)
Base names: B0 = the incumbent at the time (446bj2z or later); every arm below is a DE-only change on a bulk-identical base
(pds/mse must be bit-identical), so there is no re-level. `ARC=experiments/458-arc-lab-programme/work/arch1d_programme_abc.npz`
(sha `6951dbc7…`, 458d's block-prior artifact); `MZ=experiments/517-flex-programme/work/arch1meanz_ne_hook_abc.npz`
(sha `165c6b78…`, the hook source, §2.6). Every hook below uses `$MZ` (better on both instruments); `$ARC`'s P_down in the
same slot is the fallback (H1 nmae −0.0148 instead of −0.0183).

| 458d (fid) | 458g−458d | arm 1 (next slot) | arm 2 | October |
|---|---|---|---|---|
| POS (≥ +0.0015) | POS / NULL (le1 kept) | **B0 + `VCC_SAMELINE=$ARC:6951dbc7…:0` + `VCC_DE_PRIOR_GATE=le1` + `VCC_DE_MDE_SIGN=file:$MZ:165c6b78…:0.05`** (the hook follows the gate) vs the same without the hook | 458e (w 0.5) as planned; q90 not worth a slot | `DEF_DAY_PROGRAMME_SOURCE=h1_2025_full`, `DEF_DAY_PRIOR_GATE=le1`, hook `programme:0.05` if arm 1 transfers |
| POS | HARM | B0 + `VCC_SAMELINE=$ARC:6951dbc7…:0` + `file:$MZ:165c6b78…:0.05` (ungated) | if reach ≤ −0.005: `VCC_DE_SIGN_PROTECT_Q=0.7` on 458d (§2.2) | programme, no gate |
| PARTIAL | any | 458e (planned 10-01) | **B0 + `VCC_DE_MDE_SIGN=file:$MZ:165c6b78…:0.05:0.05:le1`** (hook alone, CLEAN4 block) | decided by 458e + the hook's read |
| NULL | any | **B0 + `file:$MZ:165c6b78…:0.05:0.05:le1`** (the hook's bar is lower than the block's: §2.4) | none: the recipient-side axes do not generalise across cancer lines (§3) | CLEAN4 block; hook source per the arm-1 read |
| HARM | — | 496a (CLEAN4 hook, 10-05, as planned) | `file:$MZ:165c6b78…:0.05:0.05:le1` only if 496a transfers | CLEAN4 block and hook |

**At 00:30Z 09-30, after `read468.py`:** `uv run python experiments/517-flex-programme/predict517.py --fid-458d <Δfid>
--fid-458g-vs-458d <Δfid> --nmae-458d <Δnmae>` prints the implied r and each hook arm's H1-equivalent and board-band nmae
(e.g. a POS-edge read, Δfid +0.0015 → r 0.23 → gated τ .05 hook −0.0064 H1-equivalent, board −0.0022 … −0.0113; a null
read, +0.0002 → r 0.26 → −0.0056; a harm read, −0.0015 → r 0.30 → −0.0046; the CLEAN4 hook: −0.0015).
Caveat on every hook price: the bar is the DESIGN's accuracy on the recipient's uncalled gate genes (H1 0.532; the public
analogues' single-donor designs 0.506–0.519). If A/B/C's design is better there (513f infers A/B/C's gate sign accuracy
above H1's from the MDE's 1.76× transfer), both hooks, Arc's and CLEAN4's, shrink by the same gap; their ratio (~4×) holds.

**October knobs per branch:** POS: `DEF_DAY_PROGRAMME_SOURCE=h1_2025_full DEF_DAY_PRIOR_GATE=le1` + the hook as
`DEF_DAY_MDE_SIGN=hookbank:h1_2025_full:0.05 DEF_DAY_HOOK_ESTIMATOR=meanz_ne` (follows the gate) or `programme:0.05` (P_down);
PARTIAL / NULL (if the hook arm transfers): `DEF_DAY_MDE_SIGN=hookbank:h1_2025_full:0.05:0.05:le1 DEF_DAY_HOOK_ESTIMATOR=meanz_ne` (block stays CLEAN4);
HARM: `DEF_DAY_MDE_SIGN=clean4:0.05` only if 496a transfers. Each also needs its rehearsal (a 497dr-like build) first.

**Arm-1 details (the hook arms), for the preregistration:**
- Config diff vs its no-hook base: exactly `de_mde_sign` + `de_sign_emitter_sha256` (the 496 emitter pin), plus nothing else.
- H1 price at full same-line quality, meanz_ne source: nmae −0.0183 (on a CLEAN4 block, gated) / −0.0182 (on 458g) raw,
  reach 0; with 458d's P_down as the source −0.0148 / −0.0157, reach −0.003. pds/mse exact. At the r implied by 458d's
  read (§2.4, P_down curves; meanz_ne ~1.2× these): r 0.2 → −0.0072, r 0.25 → ~−0.006, r 0.3 → −0.0045.
- Board transfer: the hook's own H1 → board nmae factor is unknown until 496a reads (10-05; its H1 price is −0.00146).
  Until then use the MDE family's 1.76× as the high end and 446bja's programme 0.35× as the low end.
- **Draft bins (Δ vs the no-hook base, raw):** falsifiers pds/mse bit-identical, |Δfid| and |Δjac| ≤ 0.001 (a call leak
  otherwise); primary **nmae ≤ −0.0030 TRANSFERS** (October adds the hook), −0.0030 … −0.0008 PARTIAL (keep only if 496a
  also transferred), > −0.0008 FAILS; reach reported, not binned (±0.01 is small-gate noise).

## 7. §3 — Flex technical components: library size / composition (`tech517.py`, `stream517.py`)
The scorer's DE reads per-cell CPM, so a perturbation that shifts its cells' libraries moves every gene whose CPM depends on
the library (per-gene elasticity eps_g = d ln(mean CPM) / d ln L) by eps_g × delta_t, with no biology of its own for g.
Both factors are real; together they carry almost no sign information.

| measurement | value |
|---|---|
| KAs | per-perturbation cell counts == the official pert_counts CSVs (300); 38,176 unique controls; h1-test truth recomputed vs the reference lfc, max 2.3e-15 |
| library shift delta_t = E[ln L \| t] − E[ln L \| NTC], 300 public H1 perturbations | mean −0.020, median −0.005, sd 0.061 (sampling se at 1,071 cells 0.011), 58 % negative; h1-test's 85: mean −0.014, 67 % negative, range −0.29 … +0.17 |
| delta_t vs response strength | Spearman −0.45 with the bank row's norm (300), −0.44 with h1-test n_real: strong perturbations shrink libraries |
| eps_g reliability (control split halves) | H1 0.969; A / B / C 0.85 / 0.84 / 0.87 (a real per-gene property; sd 0.12 on expressed genes) |
| eps_g agreement across Flex contexts (genes ≥ 5 CPM in both) | H1–A 0.41, H1–B 0.22, H1–C 0.25, A–B 0.36, A–C 0.48, B–C 0.33: partly shared |
| size of the effect | at delta −0.05, \|eps × delta\| is 0.0045 log2 (median; q95 0.016) against a median gate-gene \|lfc\| of 0.27 |
| h1-test gate genes: sign accuracy of eps_g × delta_t (per-target ORACLE delta) | 0.527 (65,045 pairs); by \|delta\| tercile 0.494 / 0.521 / 0.537 |
| target-independent form eps_g × median delta | 0.523 (CLEAN4 on its genes 0.573; Arc 0.640) |
| added to a programme (z-sum, λ 0.25 … 1) | Arc 0.640 → 0.639 … 0.618; CLEAN4 0.573 → 0.571 … 0.563: it only hurts |
| correlation with the programmes | r(eps·sign, Arc 0.5−P) +0.06; with CLEAN4 +0.14; with the targets' mean truth lfc −0.09 |
| share of a target's gate-lfc variance along eps (least squares) | median R² 0.015 (q90 0.09); CYFIP1 0.44, LBR 0.33, ZNF92 0.22 — but with fitted coefficients 20–280× the library shift and of either sign: those responses lie ALONG the eps axis (a cell-size / composition programme), they are not caused by the shift |
| NTC guides (46 × 400 per context): between-guide variance of mean ln L vs sampling | F 0.92 / 0.90 / 0.89 (A / B / C): no guide-level library effect |
| NTC guides: guide pseudobulk deviations projected on eps vs the guide's mean ln L | r +0.50 / +0.60 / +0.42: the mechanism is real even in cells with no perturbation (a group's sampling library drift moves its CPM along eps) |

**Reading.** Flex data do carry a composition component (reliable eps_g, systematic library shrinkage in strong perturbations),
but it is two orders of magnitude below gate-gene effects, predicts gate signs at 0.52–0.53, and subtracts from every
programme it is added to. It is not the missing part of Arc's programme (r 0.06 with it) and not a transferable sign
source. Testh's "library-size tilt" is an emission device that re-shapes the per-cell CPM distribution at fixed column totals
(our K-class A-channel already does this by design); their own board decomposition gave it −0.0026 net, which fits a device
with no signal behind it. The one real, unmodelled feature is that strong perturbations' cells are smaller (delta down to
−0.29): our carriers take control libraries. Through the per-cell CPM that matters only via eps (above: negligible) and via
the Wilcoxon's zero fraction on sparse genes (calls, not signs). Not pursued: no sign information, and the call-set side
is 192/385's closed territory.

**Per-gene detection (dropout), `dropout517.py`.** Excess dropout x_g = Poisson-expected detection rate (given the gene's
pooled fraction and each control cell's library) minus the observed one, in each context's controls:
- a shared per-gene property: r(H1, A/B/C) 0.34–0.41, among A/B/C 0.51–0.69 (probe / chemistry, or state-dependent genes);
  A/B/C have more of it than H1 (median 0.015–0.020 vs 0.009; 8–9 % of expressed genes > 0.05 vs 1.5 %);
- **no sign information:** h1-test gate genes are 0.43–0.58 down across expression × dropout cells, with no trend;
- **but it predicts CALLS:** at matched expression, the most dropout-prone third of genes is in h1-test's DE sets 2–3× as
  often (lowest expression quintile: call rate 0.033 / 0.067 / 0.109 by dropout tercile). A 2,200-gene block chosen by a
  call propensity (log CPM + excess dropout, fitted on half the targets, `callprop517.py`) captures 36.7 % of the held-out
  targets' gate genes against 33.9 % for the shipped top-CPM block (per-target mean 0.446 vs 0.390), with no loss in the
  block's median |reference lfc|. **A lead outside this mandate's territory (calls, jac; not signs):** jac ≈ overlap / n_pred
  at our call-set size, so +8 % overlap would be ~+0.006 raw jac on H1, IF the emitter can make those lower-count genes
  significant (441 PS3 lost jac re-choosing the block by prior confidence, from power). It needs an emitter block-key hook
  and a real H1 run; D/E/F's own controls supply the feature.

## 7b. Files (all in this directory; outputs in `work/`)
| script | what | peak RSS |
|---|---|---|
| `price517.py` | the analytic pricer (490's recon, exact KAs on 3 bases) and the hook grid on c4 / g33 / arc | 1.3 GB |
| `blockval517.py` | validates the analytic block-prior swap against 462l / 467 / 470's real runs | 1.4 GB |
| `grid517.py`, `hookgrid517.py`, `mzgrid517.py`, `mzblend517.py` | q / w / attenuation grids (block, hook, block + hook; P_down and meanz_ne) | 1.3 GB |
| `rank517.py`, `acc517.py` | §2's sources (lab consensus, identity, own-control PCs), emitter prices and 462f-exact accuracies | 3.5 GB |
| `halves517.py`, `halves517b.py` | 462's control-halves rule for Arc's programme and for the identity axis | 5.0 / 1.8 GB |
| `xline517.py`, `labhook517.py`, `progest517.py`, `simexpr517.py` | public cross-line instruments: identity axis, same-lab hook analogue, programme estimators, expression-similarity filter | 0.8–1.9 GB |
| `variant517.py`, `variantblock517.py`, `mzg33_517.py`, `pertarget517.py`, `c4var517.py`, `wcons517.py` | estimator variants as hook / block, on each base; per-target breadth; CLEAN4-hook variants | ≤ 1.4 GB |
| `tech517.py`, `dropout517.py`, `callprop517.py`, `stream517.py` | §3: library size / composition, dropout, call propensity; the memory-bounded CSR streamer | ≤ 3.0 GB |
| `artifact517.py` | the meanz_ne hook-source artifact for A/B/C (`work/arch1meanz_ne_hook_abc.npz`) + its KA | 0.2 GB |
| `patch_builder517.py` + `test_patch517.py` | the builder patch (on 496's) and its 26 parse/load tests | — |
| `patch_defday517.py` + `test_defday517.py`, `patch_sameline517.py` | October pass-throughs (on 502's def_day; def_sameline's estimator) and tests | — |
| `run517.py`, `drive517.sh`, `summ_real517.py`, `append_real517.py` | the real-emitter H1 runs (after the build queue) and their summary (§9) | ~7.5 GB each |
| `predict517.py` | 09-30 read → implied r → each hook arm's expected nmae | — |
| `PREREGISTRATION-draft-hook.md` | the hook arms' draft bins | — |

## 8. Deviations
- Memory: `halves517.py`'s first version peaked at **6.7 GB RSS for ~16 s** (13:45Z, under the ≤ 4 GB rule; free+inactive
  was ~29 GiB, nothing else was affected). macOS did not reuse freed per-chunk allocations; the rerun with preallocated
  buffers still peaked at 5.0 GB (the h1-test half). `stream517.py` (fully preallocated, `np.take(mode="clip")`) now
  streams h1-test at ≤ 1.8 GB; every later script uses it. `rank517.py` / `acc517.py` peaked at 3.5 GB (under 4 GB,
  over the 2 GB step size, while free+inactive was 29.6 GiB < 33 GiB). `tech517.py` peaked at 3.0 GB (same caveat).
  `xline517.py`'s first run peaked at **6.6 GB for 30 s** (np.load mmap of the two 2.7 GB Orion banks: file-backed pages
  count in RSS; free+inactive was 30.6 GiB); the rerun reads row chunks with np.fromfile and peaks at 0.84 GB. All other
  runs ≤ 1.4 GB. `wcons517.py` peaked at **4.7 GB** (run twice, ~40 s each; `def_sameline.load_bank` mmaps the Orion banks,
  file-backed pages). No step was run while a build was starting; nothing else was affected.

## 9. Real emitter runs on H1 (auto-written by `append_real517.py` from `drive517.sh`)

Known answers: PASS h1_c4-none == ~/Code/virtual-cell-challenge-2026/experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz; PASS h1_g33-none == ~/Code/virtual-cell-challenge-2026/experiments/467-strength-gated-programme/work/h1_gated33/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz; PASS h1_c4-c4t0.05 == ~/Code/virtual-cell-challenge-2026/experiments/496-uncalled-sign-hook/work/h1/stats-base_s0_bps_hookem_mde0.1cap0.1min0.3coarseas0.05.npz. Peak RSS per run: 7.0–7.6 GB.

Δ raw vs the no-hook base run of the same seed (c4 = 435, g33 = 467's gated33, arc = 462l; mz/mzg block swaps vs 435).
unsupp = mean `call_set_unsuppressible` per target (base → arm): the leak the analytic path cannot see.

| mode | seeds | pds / mse | nmae | fid | reach | jac | overall-equiv (seed range) | unsupp base → arm | analytic (H1) |
|---|---:|---|---:|---:|---:|---:|---|---|---|
| arc-none-q0.7 | 1 | exact | +0.00367 | -0.00245 | +0.03261 | +0.00054 | +0.00390 (+0.00390 … +0.00390) | 866 → 1314 | vs 458d's construction: fid −0.0040, nmae +0.0036, reach +0.0308 |
| c4-arcg0.05 | 3 | exact | -0.01489 | -0.00002 | -0.00696 | -0.00004 | +0.00282 (+0.00216 … +0.00326) | 5 → 3 | nmae −0.0148, reach −0.0031 |
| c4-arcg0.05-r0.25 | 3 | exact | -0.00627 | -0.00001 | -0.00519 | -0.00004 | +0.00076 (+0.00007 … +0.00122) | 5 → 3 | nmae ≈ −0.0059 (interp.) |
| c4-arcg0.1 | 1 | exact | -0.01064 | +0.00009 | -0.00143 | -0.00001 | +0.00274 (+0.00274 … +0.00274) | 5 → 5 | nmae −0.0106, reach 0 |
| c4-arct0.05 | 1 | exact | -0.01657 | -0.00003 | -0.02440 | -0.00003 | +0.00006 (+0.00006 … +0.00006) | 5 → 3 | nmae −0.0167, reach −0.0175 |
| c4-c4t0.05 | 1 | exact | -0.00147 | +0.00007 | +0.00193 | -0.00005 | +0.00078 (+0.00078 … +0.00078) | 5 → 3 | nmae −0.00146 (= 496's real run) |
| c4-mzg0.05 | 3 | exact | -0.01832 | +0.00022 | -0.00387 | -0.00012 | +0.00444 (+0.00374 … +0.00488) | 5 → 12 | nmae −0.0183, reach 0 |
| c4-none | 1 | exact | +0.00000 | +0.00000 | +0.00000 | +0.00000 | +0.00000 (+0.00000 … +0.00000) | 5 → 5 | — |
| g33-arcg0.05 | 3 | exact | -0.01566 | +0.00475 | -0.00681 | +0.00025 | +0.00584 (+0.00522 … +0.00632) | 591 → 393 | nmae −0.0157, reach −0.0031 |
| g33-arcg0.1 | 1 | exact | -0.01124 | +0.00391 | -0.00151 | -0.00009 | +0.00497 (+0.00497 … +0.00497) | 588 → 576 | nmae −0.0112, reach 0 |
| g33-mzg0.05 | 1 | exact | -0.01810 | +0.00712 | -0.00161 | +0.00017 | +0.00876 (+0.00876 … +0.00876) | 588 → 793 | nmae −0.0182, reach 0 |
| g33-none | 1 | exact | +0.00000 | +0.00000 | +0.00000 | +0.00000 | +0.00000 (+0.00000 … +0.00000) | 588 → 588 | — |
| mz-none | 1 | exact | -0.01593 | +0.01652 | +0.01962 | -0.00078 | +0.01689 (+0.01689 … +0.01689) | 5 → 148 | fid +0.0178, nmae −0.0160, reach +0.0171 (vs CLEAN4) |
| mzg-none | 1 | exact | -0.01224 | +0.01327 | +0.02188 | -0.00075 | +0.01448 (+0.01448 … +0.01448) | 5 → 124 | fid +0.0146, nmae −0.0124, reach +0.0213 |

With the incumbent's zs rider (hook + zs vs 435's coarse + zs, seed 0; additive if nmae ≈ the plain run's):

| mode | nmae | fid | reach | jac | pds | mse |
|---|---:|---:|---:|---:|---:|---:|
| c4-arcg0.05_zs (seed 0) | -0.01481 | +0.00001 | -0.00458 | -0.00003 | +0.00000 | +0.00000 |
| c4-mzg0.05_zs (seed 0) | -0.01827 | +0.00022 | -0.00146 | -0.00012 | +0.00000 | +0.00000 |

## Correction (2026-09-29, the 533 claims audit)
D11 (§7: Spearman −0.4448 with the bank row's norm and −0.4348 with h1-test n_real per work/tech517.json, i.e. −0.44 / −0.43). Source: `experiments/533-science-lead/RESULTS.md` §2.3.
