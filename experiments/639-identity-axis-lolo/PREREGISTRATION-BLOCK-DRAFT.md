# 639 — preregistration DRAFT: the identity axis as the DE BLOCK prior (block-form arm)

Drafted 2026-10-02, first committed 13:02Z (1f794019) by [d4cd5f] (639's owner), on the coordinator's request routed by [21dd4c]. It becomes a
preregistration when the coordinator fixes base, day and names and hashes it into the ledger. Bins fixed now.

## Arm (ONE variable vs the base): `ID05`, and `ID10` as the dose partner
- Base and day: the coordinator's. Base = the best scored entry when the slot is armed. ~~ID05 first; ID10 only if ID05
  reads ≥ WEAK~~ (superseded by R3: ID10 first; ID05 only if ID10 reads TRANSFERS). Queue: behind the 10-07 tree and 630's dose pair, ahead of the KGID read (coordinator, 13:09Z).
- ID05 = base + `VCC_SAMELINE=experiments/604-identity-prior/work/ident604v3_flex_lam0.5_ABC.npz:6cd84dc258cb609d08d93ec4c6c6ec7c7134ff8e3093eccde5ac15bbd2ac618f:0`
  (w 0: no same-line rows; the artifact's `programme_<c>` is each context's block prior).
  ID10 = the same with `experiments/604-identity-prior/work/ident604v3_flex_lam1_ABC.npz:3b1e22f037aff26b7a13e06328d011757423a51a8cf4acaef7d26cd206fd6a07:0`.
  Shas verified with `shasum -a 256` at 13:10Z; 604's v3 KAs 39/39 PASS (`604/work/ka604v3.log`, ALL PASS). NOT v1 (`ident604_flex_lam0.5_ABC.npz`: its
  pdown is finite on 10,962 genes, i.e. it carries the identity-only extension 639d found harmful). v2 (6,346 genes =
  CLEAN4's coverage) fixed the extension; v3 also fixes the side rule (iii). ([21dd4c]'s review, applied in 187ee529 at
  13:07Z.) v3 at λ 0.5 vs CLEAN4: flips 597 / 601 / 648 (block 164 / 182 / 187) in A / B / C; block down-fraction
  0.739 / 0.757 / 0.727 (CLEAN4 0.677 / 0.678 / 0.659).
- **KA for the arm:** the v3 λ 0 artifact (`ident604v3_flex_lam0_ABC.npz`, sha
  f65d12fc3a7d8a7557e071ea6f86dfc1e89847cd0eea859b6ac14b26ef6fc6e6, byte-identical to v2's λ 0) through the builder
  rebuilds the base bit for bit. `programme_<c>` is ±1 on all
  18,533 genes in every artifact, so this KA is what proves the out-of-coverage default is unchanged.
- **Artifact requirements (from 639d):**
  - (i) Flex-only references (A/B/C + h1), per 622;
  - (ii) the combined score is applied ONLY on CLEAN4-covered genes; identity-only genes get NO prior. 639d: the
    extension is right 0.27–0.34 on HEK293T;
  - (iii) the SIDE is sign(z(P_C4 − ½) + λ·z(I)), the rule 639 / 639d validated (centred at P = ½, not rank-mapped:
    a rank map keeps CLEAN4's down-count and would undo part of the identity down-shift in A, where it is largest:
    A 58 % down by the sign rule vs 53.5 % rank-mapped). The CONFIDENCE |P − ½| is CLEAN4's own |P − ½| distribution over
    covered genes, assigned by the rank of |score|, so the builder's confidence calibration (q80, shrink, τ) is kept
    ([4d0207]'s v3 form, decided and committed 13:06Z, 99ce3c38);
  - KAs: λ 0 reproduces CLEAN4's P exactly (stable ranks); v3 equals v2 wherever the two side rules agree;
    uncovered genes keep CLEAN4's behaviour (P NaN).
- Expected config diff vs base: exactly `sameline` (the artifact path + sha, w 0).

## Mechanism and class
"Identity genes go down": a perturbation erodes the line's identity programme. Seen in 639 (LOLO), 622g (Broad-16,
another lab, Cas9 KO) and 627 (Tahoe's 50 lines under drugs, Spearman +0.147). I = log1p ctrl_c − mean log1p
ctrl_others, from the context's own controls, within Flex. **Class: TRANSFER** (a better prediction of the responses),
computable on the day for D/E/F (vs A/B/C and the other two D/E/F contexts).

## Gates (before the slot)
- **554 dose gate:** each context's median prior_lambda ≥ 0.9 in the build's audit dump. Below it, the arm is not
  submitted (458d's 0.43 lesson).
- Block sign changes vs CLEAN4 and the block down-fraction reported per context (604's audit).

## Instrument validation FIRST (the A1 lesson; [21dd4c], run by [21dd4c])
Before ID05's numbers are read, the same restricted-LOLO instrument (block genes only: gate genes AND null genes,
q ≥ 0.05) must retrodict the one board pair of a block-prior source change: **458d** (Arc's H1 programme as every block's
prior, dose ~0.4; board fid −0.0055 raw, reach −0.0267 raw, i.e. HARM). **Pass rule (fixed now):** 458's programme vs
CLEAN4 is negative on BOTH the gate-gene and the null-gene block readings on ≥ 2 of the 3 cancer folds
(HCT116 / HEK293T / K562). **If the instrument fails, it cannot vouch for ID05:** the WEAK bin below becomes "October stays
CLEAN4", and only TRANSFERS ships (mirroring KGID's A1 FAIL branch).

## The readings (11.1c.2)
0. **Null-block (fid filler) reading ([21dd4c]):** identity (v3, λ 0.5) vs CLEAN4 sign precision on NULL block genes
   (q ≥ 0.05) on HCT116 / HEK293T / K562. The block is ~2,200 of 2,600 calls and ~95 % null, so this IS fid. H1 so far:
   0.5186 → 0.5201 (612/fillerid612), but the block is not identity's route on H1. Caveat: identity is computed from a
   surface's own control pool, the same pool its truth lfc is taken against, so a full-pool reading can carry a small
   shared-control-noise agreement (635). It is real on the board too, since A/B/C's identity comes from the reference
   controls, but it should be reported separately where control halves exist (H1).
   **Rule:** if identity's null-block precision is below CLEAN4's on ≥ 2 of 3 cancer folds, fid is expected to eat the gate
   gain, and the arm ships only on TRANSFERS.
1. **Held-out gate genes (fixed now, from 639d):** LOLO gate-sign accuracy restricted to BLOCK genes, λ 0.5: HCT116 +0.040, HEK293T
   +0.033, H1 +0.006, K562 −0.029 (low power). 2 of 3 cancer folds positive, robustly.
   **Tail ([7772dc], 622g's Broad-16):** the same block-restricted gain per held-out line (block = top 2,200 by control
   expression after a 400-gene head by the other lines' consensus |lfc|): the count of negative lines out of 12 is the
   estimate for a hostile D/E/F line, **broken down by lineage** (coordinator, 13:09Z). K562's −0.029 matters because
   context A is Jurkat, also haematopoietic, and the board pools A with B/C. If the losses cluster in blood lines, the
   record says so. A day-of lineage gate would be a separate preregistered switch, not a tweak inside this arm.
2. **Board:** the prior-only arm shares every emission draw with its base (468), so the read is exact. The primary
   statistic is Δ(fid + reach + nmae) overall, as 458's bins.

## Readings filled (13:21Z, before any entry with this knob exists)
- **458d retrodiction: PASS** ([21dd4c], 613/RESULTS.md, `block613.py` committed 2ac94af2 before its first run; KA =
  decomp639's identity gate-block gains, exact on all 4 surfaces). The H1 programme minus CLEAN4 on cancer-fold block
  genes, gate / null: HCT116 −0.170 / −0.026, HEK293T −0.194 / −0.038, K562 −0.085 / −0.008. Negative on both readings,
  3/3 folds. Size: the fid estimate on applied cells at full dose is −0.004 … −0.034 raw, ×0.43 dose → brackets the
  board's −0.0055. On H1 itself the same programme reads +0.087 / +0.028 (own-line programmes do not transfer). The
  instrument counts.
- **Reading 0, null block: identity BELOW CLEAN4 on 2/3 cancer folds** ([21dd4c], 613): HCT116 −0.0137
  [−0.0149, −0.0125], HEK293T −0.0081 [−0.0097, −0.0066], K562 +0.0001 (ns). H1 control halves: +0.0012 vs +0.0013
  full pool, so the shared-control effect is negligible. **By the rule above, ID05 ships only on TRANSFERS.** The
  descriptive fid exposure on applied cells (Δfid ≈ Δcorrect / 2,600) is HCT116 −0.0060 raw (−0.0033 overall), HEK293T
  −0.0004, K562 −0.0002, H1 +0.0009. That is the size of the predicted gate gain.
- **Reading 1 tail, Broad-16 block-restricted: NEGATIVE** ([7772dc], 622 RESULTS §(m), `block622.py`, commit 4cd40535;
  approximate production mask; λ 0.5, same-lab refs, no identity-only extension). Block: mean −0.0022, 7/12 lines
  negative, min −0.0265 (T3M5, head & neck, CI [−0.044, −0.0095]), CI < 0 on 2 lines, CI > 0 on 1. The same lines'
  head +0.0133 (3/12 negative) and **uncalled +0.0204 (0/12 negative, CI > 0 on 3)**. Pooled share of the net gain:
  block −0.11, head +0.28, uncalled +0.83. Lineage: no blood lines in the set. **Negatives cluster in epithelial
  carcinomas** (6/9 negative, block mean −0.0069; cervix SKGII −0.014, head & neck T3M5 −0.027: the B/C look-alikes);
  other lineages 1/3 negative, mean +0.012. Raw-lfc head as a check: block −0.0024, 7/12 negative. The gain is identity,
  not expression level (C4 + 0.5·z(log1p CPM) +0.0011; I with expression regressed out keeps +0.0052 of +0.0081).
- **Route: UNRESOLVED, and the weight of evidence is AGAINST the block.** The block carries the gain on 2 lines from one
  lab (HCT116, HEK293T). The uncalled genes carry it on 13 lines (Broad-16's 12 at 0/12 negative, K562 +0.054) and on H1
  (0.69 of the net). With the null-block cost on top, ID05's expected board value is ≈ 0 with HARM plausible.
  **Owner's recommendation (d4cd5f): demote ID05 below KGID.** The identity axis's robust route is the UNCALLED-gene
  hook with CLEAN4-covered genes only (no extension), which is 612's KGID with a v2/v3 artifact.
  **Coordinator's queue (13:21Z):** ID05 is demoted but kept ABOVE KGID: 630 → 524dr2 → factorial validation → 625 →
  ID05 → KGID. KGID is a delivery-class change that 611's H1 multipliers already price, while ID05 is a prior-type
  change only the board can price. Only TRANSFERS ships; the bins below stand.

## Amendment R1 (13:48Z, before any entry with this knob exists): reading 1 is ARTIFACT-DOMINATED on Orion
[4d0207]'s 609 (`609/RESULTS.md`, d6bc187c; its `rider` row reproduces decomp639.json exactly, and block609 reproduces
block622 exactly):
- On HCT116 / HEK293T the block's gate genes are 72–74 % DOWN (the Orion depth skew; AGENTS 6.5,
  HARNESS-DEPTH-CONFOUND). An expression-only rider (zc + 0.5·z(log1p own CPM)) earns block "gain" +0.047 / +0.084,
  as much as identity's +0.040 / +0.033 or more.
- With expression partialled out, HEK293T's block gain is −0.006 and HCT116's halves to +0.025.
- Clean surfaces: H1 block +0.0057 (+0.0080 with I⊥E; expression alone −0.0078); Broad-16 under production's uncentred
  zc block −0.0007, uncalled +0.0231 (share 0.92).

**Reading 1 is replaced by: H1 +0.006, Broad-16 −0.001 (block-restricted).** The predicted board size drops ~5–6×, to
**≈ +0.000 … +0.001**, below the TRANSFERS bin (+0.0015). Since only TRANSFERS ships, **the expected decision is NULL →
October stays CLEAN4 on the block**. The arm remains an exact, cheap prior-only read of a small effect. The owner's
recommendation (d4cd5f): keep it only as an empty-slot arm; the queue position is the coordinator's call.
Identity's gain lives on UNCALLED genes, and its delivery route is hook signs × 630's magnitudes (630's stack on 612's
idkgd: nmae −0.0105 vs −0.0080 raw at f 1).
Lesson: any Orion-fold result about highly expressed genes needs an expression-only control first.
**Independent sizing (14:09Z; [7772dc]'s 625e, 982b00d4):** route-mask sizing calibrated on 625d's real-emitter octile flips
(delivered ≈ 0.7–1.0× the source's gate-sign change on route genes; overall ≈ 0.2× that) predicts an ID05-like block
prior on H1 at **≈ +0.0003 / +0.0006 overall** (f 1 / f 0.35, λ 0.5; λ 1: +0.0005 / +0.0010). Identity gains +0.0048 /
+0.0100 on all gate genes, but only +0.0015 / +0.0030 on block genes; most of the rest sits on uncalled genes the
current hook cannot reach (hook route +0.0001). Two independent estimates agree: **NULL is the modal outcome.**
**Correction (15:59Z; 625 FINAL, 73c0b648):** for IMPROVEMENTS the route-mask factor is ≈ 0.35 (four oracle arms,
median; range 0.14–0.50), not 0.2 (which came from degradation flips). Revised ID05-like prediction on H1:
λ 0.5 +0.0005 / +0.0011 (f 1 / f 0.35); λ 1 +0.0009 / +0.0017. λ 0.5 stays below TRANSFERS (+0.0015), so NULL is still
modal for ID05. λ 1 reaches the bar only at f 0.35. The registered order is unchanged (ID05 first; ID10 only if ID05 ≥
WEAK). An ID10-first switch would need its own amendment before any entry, and the coordinator's call. Call-set
term: improving flips of 150–300 genes add only +1 … +47 calls per target, below the whole-octile cascade.

## Amendment R2 (13:55Z, before any entry exists): the block prior moves the CALL SET ([7772dc]'s 625d)
On H1's real emitter, flipping block-prior genes DOWN→UP EXPANDS the realized call set, while UP→DOWN does not shrink it.
- Octile 5's flip (793 genes) adds +281 calls per target (median +189, up to +1,449). Only 32 of them are octile-5
  genes; the rest are a BH / UP-design cascade.
- fid moved −0.0049 at f 0.35 from octile 5 alone.
ID05's artifact flips DOWN→UP on 149 / 287 / 302 genes in A / B / C (UP→DOWN 448 / 314 / 346), so expect roughly +50–100
calls per target.
**Changes:**
- (a) The falsifier "jac moving > 0.0010 raw" is WITHDRAWN: its premise ("the block prior should not change the call
  set") is false.
- (b) The read records each context's realized calls per target, ID05 vs base (rc diagnostics). A jac / fid move that
  tracks the call-count change is call-set mechanics, not sign quality.
- (c) **nmae on the reference's DE genes is the primary sign read.** The bins stay on Δ(fid + reach + nmae) as
  registered; the NULL expectation from R1 stands.

## Amendment R3 (2026-10-02T16:01Z, the coordinator's, before any read): ID10 FIRST
If an empty slot opens for the identity block prior, **ID10 (λ 1.0; `ident604v3_flex_lam1_ABC.npz`, sha 3b1e22f0…) runs
FIRST**. It is the larger dose, so it has more power to show whether the block route exists on Arc lines at all; a NULL at
λ 1 closes the route more firmly than a NULL at λ 0.5. **ID05 runs only if ID10 reads TRANSFERS**, to pick the dose.
- The bins below apply to whichever arm is read, with ID10 in ID05's place.
- The dose rule is unchanged (both ≥ +0.0015 → adopt the larger).
- The 554 dose gate, the KA (v3 λ 0 builds the base bit for bit), the realized-calls record (R2) and the
  TRANSFERS-only rule (R1) all apply to ID10.
- Predicted for ID10 (625 FINAL, H1): +0.0009 / +0.0017 at f 1 / f 0.35.

## Note (2026-10-02T19:11Z): the 616 audit voids the Orion-based readings above
[21dd4c]'s 616 audit (7dac7a98) found the thinned Orion gate tables' DE truth ~98 % selection artefact. Void:
- 613's 458d retrodiction PASS, now on K562 only (−0.085 / −0.008, still negative on both readings);
- the null-block rule's HCT116 / HEK293T rows (K562 +0.0001 ns remains);
- 639d's Orion block gains (already withdrawn).
The arm's decision rules are unchanged (TRANSFERS-only, ID10 first, bins as written): they are now conservative rather
than evidence-driven. The coordinator is re-ruling the open calls from the audit, and this prereg follows that ruling.

## Bins (primary: ID05's board Δ(fid + reach + nmae) overall vs the base)
- **≥ +0.0015: TRANSFERS.** Adopt for A/B/C; October `programme_source = identity` (def_day computes I from D/E/F's
  controls; 604's patch).
- **+0.0005 … +0.0015: WEAK.** Adopt only if ID10 is ≥ ID05 − 0.0005 (no dose collapse); October adopts, flagged.
- **−0.0005 … +0.0005: NULL.** October stays CLEAN4.
- **≤ −0.0005: HARM.** Close the block route; record the 458d analogy.
- **Falsifiers:** |Δpds| or |Δmse| > 1e-6 raw (a prior-only arm moved the bulk); jac moving > 0.0010 raw (call-set
  composition changed, which the block prior should not do).
- **Dose rule:** if both are ≥ +0.0015, adopt the larger. ID10 is read only on the block (its uncalled-gene loss does not
  apply: the hook stays as the base has it).

## Try count (11.1c.4)
λ: 4 values on 639 (0.5/1/2 + 0), 2 on 622 / 639d. Reference sets: 2 (mixed, Flex-only). Coverage forms: 2 (with and
without extension). Surfaces: 636's four + Broad-16 (622g) + Tahoe (627). Board: 0.

## Predicted size (board)
Block gate-gene signs +0.033 … +0.040 on ~50 % of the folds' gate |l| mass. Through the block's delivered magnitudes
(~0.1–0.15 vs truth ~0.25) and reach's head: **+0.001 … +0.004** from the gate genes. **fid is the main risk, not
negligible:** 458d, the only board pair of a block-prior source change, moved fid −0.0055 raw and reach −0.0267 raw
(≈ −0.003 and −0.005 overall), larger than this gain. The null-block reading (0) and the 458d retrodiction decide whether
the predicted sign can be trusted.
