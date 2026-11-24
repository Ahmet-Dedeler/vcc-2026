# 639 — The identity axis ("identity genes go down") TRANSFERS leave-one-line-out: the first public sign source to beat CLEAN4 on held-out lines

> **VOID ROWS (2026-10-02T19:10Z; [21dd4c]'s 616 audit, 7dac7a98, `experiments/616-clean-orion-folds/AUDIT-orion-dependence.md`):** the thinned Orion gate tables' DE truth was ~98 % selection artefact (HEK293T: 37,748 calls on the thinned table vs 774 on clean raw cells; 13 / 132 targets scoreable). In this experiment, the hct116 / hek293t rows (639, 639d, 639e) are VOID. What remains: h1 +0.016 and Broad-16 +0.009 (622g); K562 null. 639d's block-route claim was already withdrawn (609).

2026-10-02 ~12:45Z, council session [d4cd5f], portfolio item 2 (636's instrument). Analysis only (`ident639.py`,
`boot639.py`, `halves639.py`; `work/*.json`), < 3 GB, seconds.

## Construction
441 PS6 found on h1-test, within Flex, that genes expressed higher in a line than in other lines' controls fall under
perturbation (r −0.236 with H1's mean response) and that z(CLEAN4) + λ·z(−I)… (here written with the sign folded:
score = z(P_C4 − 0.5) + λ·z(I), positive → predict down) lifts gate-block signs. It was never built: H1 alone priced it at
~+0.001. Here I = log1p(ctrl CPM of surface S) − mean over the OTHER profiles of log1p(ctrl CPM). The other profiles are the
other three 636 surfaces' controls (480's `ctrl_cpm`) and A/B/C's controls (Flex). The common genes are those all seven
profiles measure (6,793). CLEAN4 is LOLO per 636.

## Result: gate sign accuracy, paired vs CLEAN4 on genes both cover (4,000-draw target bootstrap)

| surface | CLEAN4 | λ 0.5 | λ 1.0 (441's prior choice) | λ 2.0 |
|---|---:|---|---|---|
| h1 (Flex) | 0.573 | **+0.014 [+0.009, +0.019]** | **+0.016 [+0.006, +0.024]** | +0.013 [−0.003, +0.026] |
| hct116 | 0.686 | **+0.039 [+0.030, +0.048]** | **+0.015 [+0.004, +0.026]** | +0.012 [+0.000, +0.025] |
| hek293t | 0.663 | **+0.017 [+0.011, +0.025]** | **+0.024 [+0.014, +0.035]** | +0.009 [−0.004, +0.024] |
| k562 (low power, 1,626 pairs) | 0.608 | −0.005 [−0.027, +0.016] | −0.010 [−0.038, +0.019] | −0.017 [−0.050, +0.016] |

- Full coverage (combined where CLEAN4 exists, identity alone in CLEAN4's gaps, vs CLEAN4 alone), measured as net
  (right − wrong) per gate gene: h1 +0.024 [+0.011, +0.035], hct116 +0.051 [+0.033, +0.070], hek293t +0.025
  [+0.010, +0.043], k562 −0.003 (null). Identity alone beats all-down in CLEAN4's gaps on every surface.
- **Not a shared-control-noise artifact** (`halves639.py`, 635's construction): on h1 with the identity built from
  control half A and the truth taken against half B, CLEAN4 0.5732 → C4 + 1·I 0.5884, the same as with full controls
  (0.5887).
- A chemistry confound would flip sign between the Flex surface (h1) and the 3′ surfaces. The gain is positive on both.
- Forking paths: λ was swept over 4 values. λ 1 was 441's choice before today's data. Both λ 0.5 and λ 1 pass on
  3/4 surfaces with CIs excluding 0.

## Reading
- **The portfolio's bar ("beat CLEAN4 by ≥ +0.01 on most folds") is met at λ 0.5 and λ 1:** 3/4 surfaces, the fourth
  a low-power null. It is the only public sign source to do so today; 621's supervised model, 634's conditional
  direction and 638's responder-weighted programmes all lose.
- Class: TRANSFER (biology: perturbation erodes a line's identity programme), computable on the day from D/E/F's
  controls within Flex (D/E/F vs A/B/C and each other). Nothing is fitted to A/B/C.
- Delivery: it is a per-gene, target-independent sign, so it plugs into the block prior and the uncalled-gene hook.
  632 (renumbered from 640) prices it with graded magnitudes on uncalled genes.

## 639d — WHERE the gain would be delivered: head / CPM block / uncalled (committed 13:02Z, 1f794019; `decomp639.py`, `work/decomp639.json`)
Call sets: h1 = lib521's real design head (400 by |key|) + CPM block (2,200); folds = ds612c's construction
(consensus-key head + block by control expression). Identity from a FLEX-ONLY reference set (A/B/C + h1; 622: 3′ refs
hurt). Paired gate-gene sign accuracy, z(C4) + λ z(I) vs CLEAN4-LOLO; "share" = each region's share of the total change
in (right − wrong), plain and |lfc|-weighted.

| surface | λ | head gain (share) | **block gain (share net / \|l\|-w)** | uncalled gain (share) | extension acc (identity alone where CLEAN4 is undefined) |
|---|---|---|---|---|---|
| hct116 | 0.5 | +0.034 (0.03) | **+0.040 (0.99 / 0.89)** | −0.002 (−0.01) | head 0.63, block 0.63, uncalled **0.49** |
| hct116 | 1 | +0.004 | +0.030 (1.56 / 3.64) | **−0.042** | |
| hek293t | 0.5 | +0.010 (0.00) | **+0.033 (0.92 / 0.83)** | +0.009 (0.08) | head **0.28**, block **0.31**, uncalled **0.34** |
| hek293t | 1 | −0.026 | **+0.053 (0.96 / 0.92)** | +0.008 | |
| h1 (ESC) | 0.5 | +0.001 (0.01) | +0.006 (0.30 / 0.18) | +0.013 (0.69 / 0.82) | 0.54 / 0.58 / 0.55 |
| h1 | 1 | +0.005 | +0.010 (0.34 / 0.21) | +0.019 (0.64 / 0.77) | |
| k562 (1,824 block pairs) | 0.5 / 1 | +0.016 / +0.037 | −0.029 / −0.050 | +0.054 / +0.075 | 0.57 / 0.57 / 0.39 |

**Reading.**
1. **On the cancer folds the identity gain sits on the CPM BLOCK** (≥ 0.83 of the |l|-weighted net on HCT116/HEK293T), so its
   production route is the **block prior**, not the uncalled-gene hook (which takes ~0 there and loses on HCT116 at λ 1).
   On H1 (ESC) it sits mostly on uncalled genes. A/B/C/D/E/F are cancer lines.
2. **Do NOT extend coverage with identity alone.** Where CLEAN4 is undefined, identity-only signs are right 0.27–0.34 on
   HEK293T and 0.49 on HCT116's uncalled genes. The block-form arm must apply the combined score only on CLEAN4-covered
   genes and leave the rest as today (a correction to the 604 spec of 12:20Z, sent to [4d0207]).
3. **λ 0.5 is the robust dose** (λ 1 costs HCT116's uncalled genes and HEK293T's head). λ 1 is the dose partner on the
   block only, where it is +0.030 / +0.053.
4. K562 (low power) disagrees on the block (−0.03 / −0.05). Report it; it is one of four surfaces.

## 639e — the block route's held-out readings (13:21Z): route unresolved, evidence against the block
Run by peers on 639's preregistered readings (PREREGISTRATION-BLOCK-DRAFT.md, "Readings filled"):
- 613 ([21dd4c]): the restricted-LOLO instrument retrodicts 458d's board HARM (3/3 cancer folds negative on gate and
  null block genes), so it counts.
- Null-block identity vs CLEAN4 is negative on HCT116 (−0.0137) and HEK293T (−0.0081); K562 is ns. That is a
  fid cost about the size of the gate gain.
- 622 §(m) ([7772dc]), Broad-16 block-restricted: mean −0.0022 with 7/12 lines negative, concentrated in epithelial
  carcinomas. Uncalled genes are +0.0204 with 0/12 negative.

**Conclusion.** 639d's "block on cancer lines" holds on HCT116/HEK293T (one lab) only. Across 14 other held-out lines
the identity gain sits on UNCALLED genes, and the block prior also pays fid. The identity axis's delivery route is the
uncalled-gene hook (612's KGID), restricted to CLEAN4-covered genes (no identity-only extension: 639d, HEK293T
0.27–0.34). ID05 (the block form) is demoted. Only TRANSFERS ships. The coordinator (13:21Z) keeps it above KGID in the queue, since only
the board can price a prior-type change.

## 639f — 639d's block localisation was an Orion artifact (13:48Z; [4d0207]'s 609, d6bc187c)
On HCT116 / HEK293T the CPM block's gate genes are 72–74 % DOWN (Orion's depth skew). An expression-only rider earns
as much block "gain" as identity, or more (+0.047 / +0.084). With expression partialled out, HEK293T's block gain
vanishes (−0.006) and HCT116's halves (+0.025). On clean surfaces the block share is small: H1 +0.0057 (0.30), Broad-16
−0.0007. The uncalled share is large (H1 0.69, Broad-16 0.92).
**639d's "route = block" conclusion is WITHDRAWN.** Identity's sign gain is real (I⊥E keeps it on H1, and 622g's
Broad-16 keeps +0.0052 of +0.0081), and it sits on UNCALLED genes. 639d's own table already showed H1's uncalled share,
and I read the Orion rows over it without an expression control: the mistake 609 names. Prereg amendment R1: ID05's
predicted board size is now ≈ +0.000 … +0.001, so NULL is expected.

## Note (2026-10-04T06:18Z): Broad-16's global UP offset ([7772dc]'s 657 calib657, 92767352) and the identity axis
Broad-16's KO lines show a global UP sign offset that no CRISPRi bank shares: UP-right crosses 0.5 at P_down ≈ 0.55 vs
H1's 0.51, and the mean UP-right is 0.58 vs 0.50. A source can score gain over CLEAN4 there just by predicting UP more
often.
**639's identity score moves the other way.** It predicts MORE down than CLEAN4 on covered genes (A: 58 % vs 53.5 %
down; identity genes go down). So on Broad-16 the offset works against it, and its Broad +0.009 (622g) is conservative,
not inflated. The H1 reading (+0.016; Arc CRISPRi chemistry, no offset) is unaffected. Not re-scored: the direction
argument suffices for this source. A source that tilts toward UP (e.g. low-I genes UP on uncalled genes) would need the
per-line re-centring check.
