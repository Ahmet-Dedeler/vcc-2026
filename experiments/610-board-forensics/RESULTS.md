# 610 — Board forensics: the whole gap to #2–#5 is the DE channel (council 600, [21dd4c], 2026-10-02)

Analysis only: no build, no `vcc`, no `cell-eval2`. Every process < 2 GB, seconds. The narrative and the hypotheses are in
`experiments/600-breakthrough-council/21dd4c.md`. The H1 → board calibration (611) is in `RESULTS-611.md`. This file holds
the measured numbers and the negatives.

## Config
- Board: the live API at 10:43Z (`work/lb_2026-10-02T1043Z.json`) plus 14,548 distinct entries from our 211 snapshots
  (`hist610.py`).
- Truth surface: h1-test (85 targets × 400 cells vs 18,400 controls, the board's shape), with 213's real DE table and
  521's regimes.

## Results
1. **Channel split, best entry per team, overall units.**
   - Our bulk (pds + mse)/6 = 0.1993, ~3rd: Illumina 0.264, cqawesome 0.214, Testh 0.206.
   - Our DE (fid + reach + jac + nmae)/6 = 0.0736, ~24th of the top 35: cqawesome 0.114, Testh 0.100, chill 0.098,
     CellWM 0.097.
   - Against #2–#5 the DE gap is +0.025 … +0.041 overall.
   - Our DE members have been pinned since 09-25 (fid 0.5144, jac 0.0408).
2. **The fid attractor.**
   - Ten teams sit at fid 0.5275–0.5295. **Correction (4d0207, 608, 13:27Z):** five of them (Ava Brown, HY H, Anna Xie,
     SQ W, Vivian Xu) are one clone cluster with identical DE members (fid 0.528, jac 0.0337). That leaves about six
     independent designs. It is a cascade: the cellpred multi-account family moved in ONE step on 09-20/22,
     then cqawesome, Testh, The_Way, sunny and Pythia followed (`work/attractor_timeline610.txt`).
   - No public code exists (gh search; `ArcInstitute/gpudge` is only a DE library). The VCC Discord was unreachable from
     the built-in browser.
   - The pack's jac (0.032–0.037) sits BELOW ours (0.041) and has been flat for a week (608), so its fid edge is not
  more real-DE genes among its calls. Illumina's 09-28 pure-DE step is the opposite route: fid +0.016 and jac +0.008
  together, i.e. a better call set.
- Pythia reached 0.5287 within 7 entries, so the level is reachable from public information. chill's day-1 described
     design gives 0.5225 / nmae 0.835.
3. **fid is not call-set size on this board.**
   - Clean size-only curves are flat above ~2,500 calls: Chiu 0.5097 / 0.5164 / 0.5148 at ~1,075 / 2,470 / 4,300; ours
     on 09-03 0.5037 → 0.5037 at 4,480 → 2,599.
   - Disclosed big model-signed call sets land at 0.505–0.517 (`work/callset_hints610.txt`).
4. **Rivals' DE-only steps** (201 cataloged, `destep610.py`):
   - pure-nmae steps −0.015 … −0.070 raw;
   - reach +0.014 … +0.074;
   - fid +0.006 … +0.016.
   - The pure-fid steps move only FILLER calls: nmae, reach and jac stay flat.
5. **fid decomposition on our real H1 build** (`phidecomp610.py`; 435's stats, 521's regimes; f 0.35):
   - φ (the DE share of calls) is 0.051, p_DE 0.675, p_null 0.521, fid 0.525;
   - the CPM block carries MORE truly-DE genes than a target-specific top-K (0.052 vs 0.030), so the composition is right;
   - our board fid implies p_null ≈ 0.507–0.510 on A/B/C.
6. **Same gate genes, sign sources** (`gatesign610.py`, f 0.35):
   - delivered signs 0.615; head 0.838;
   - uncalled genes: the key beats CLEAN4 in the top |key| quintile (0.63–0.72 vs 0.51–0.56), and CLEAN4 beats the key in
     the bottom half (0.572 vs 0.527). This led to 612.

## Negatives (closed at this construction)
- **Null sign skew from controls** (`nullskew610.py`): random 400-cell control groups predict the reference's non-DE
  signs at 0.500–0.503, worse than all-down. The context "programme" is not a sampling artifact of the group mean.
  635 later showed the LOO programme is mostly shared CONTROL-POOL noise, which is unknowable.
- **Filler chosen by CLEAN4 sign confidence instead of call propensity** (`fillerchoice610.py`, K 2,200): fid +0.008,
  jac −0.033, a net loss.
- **The yield explanation of the pack's fid** (d4cd5f's 4b) is refuted by the flat board size curves (item 3).
- **Rank-neutral UP delivery as a large lever:** the analytic dose curve (`dosecurve610.py`) gave +0.001 … +0.005
  overall. On the real emitter (630, d4cd5f) no carrier sits above the control max, so delivery is 10–26 % of the ask:
  b_k100 +0.0030 overall at f 0.35, 3 seeds.

- **CPM coupling as the pack's null-sign edge** (`cpmcouple610.py`, `work/cpmcouple610.out`; H1 truth + 435's stats).
  - The mechanism is real. Per-cell CPM couples every null gene to its perturbation's net library-mass change.
    - The share of a target's null genes with positive real lfc spans 0.38–0.75 (sd 0.045).
    - It tracks the real-DE genes' net mass change at Spearman −0.86 (p 2e-26).
    - The mass is carried by highly expressed identity genes: CD24, TUBB2B, ACTB, DNMT3B, POU5F1, EEF2. MED19, EP300,
      the Mediator subunits and CTNNB1 lose 0.5–1.5 % of the library, and their null genes go 55–75 % up.
  - It is not a lever for us:
    - our design key's implied mass change does not predict the real one (Spearman −0.22), nor does our emitted one
      (−0.14; the emitter mass-balances);
    - even the ORACLE rule "null calls take −sign(real mass change)" gives p_null 0.5183 on our null CALLS against our
      0.5169 today (+0.0014). Our expression-ranked block carries per-gene responses that dominate the global offset.
    - On ALL null genes the oracle reaches 0.5219, but fid only counts calls.
  - So it does not explain the pack's +0.014. Per gene this is the same identity axis 639 found, and 639d's block-prior
    route is the reachable form.

- **Larger call sets that extend down the expression ranking** (`decile610.py`, `work/decile610.out`; 636's
  surfaces). CLEAN4's sign precision on NULL genes is highest on the top-expressed 2,200, which is our block, and falls
  with expression rank on every surface:
  - h1 0.518 → 0.511;
  - HCT116 0.543 → 0.501;
  - HEK293T 0.550 → 0.494;
  - K562 0.514 → 0.500.
  - So the expression-ranked block is already the precision-optimal filler, and the pack's larger call sets (lower jac)
    cannot buy fid this way. After size, filler choice, identity, CPM coupling, block sign sources (613c) and expression
    range, the pack's +0.014 fid stays unexplained by public-data levers we can test.

## One honest sentence
Our bulk is level with the #2–#4 public-data teams, and the whole remaining gap is our DE channel. Its sign information is
at the public limit. The levers left are delivery and source combinations (611, 612, 630, 622), each worth +0.001 … +0.009,
plus the unexplained ~+0.014 fid edge of the pack.
