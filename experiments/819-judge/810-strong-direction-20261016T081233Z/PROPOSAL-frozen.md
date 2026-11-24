# 810 PROPOSAL-2: strength-gated unwhitening of the donor base, one arm (sgu1), priced beside the closed global form

Written 2026-10-16, started 08:10:53Z by the clock, by the 810 worker (lane B1). This is the resubmission after
`VERDICT-1.md` (REVISE, judge run 819; sha256 b5d87daf607f19603e7d2956561fcbe44fc7f67c4ce794995250cbac538e425b).
**FROZEN once written: nobody edits this file.** Its sha256 is recorded in `work/SHA256SUMS.txt` and goes into the run's
candidate.json (`prereg_sha256`). No candidate has been priced. The measurements it cites are in `RESULTS.md` §1–§2, and
every script and work file it names carries a sha256 in `RESULTS.md` §7.

## 1. Mechanism (corrected per VERDICT-1 item B)
- **The measured facts** (composition level, 22 folds; `work/table810.md`, `work/explore810-*.txt`):
  - On every z_src tier, the whitened donor base points at the held line's truth worse than the unwhitened, normrel'd
    donor mean `pre` it is built from. The per-target loss on BULK-W is +0.0096 / +0.0103 / +0.0193 (weak / mid / strong).
  - In the strong tier, the loss rises with source strength (Spearman +0.22 W / +0.16 ABC cancer) and with z_src
    (+0.20 / +0.23).
  - H1 and CD4 help on every tier: dropping either part loses tier-X on 18/18 W and 3/3 ABC folds.
- **The candidate swap** replaces the base with `pre` at the base's norm and keeps cd4 / h1 / part and the row norm. It
  gains at composition level on EVERY tier, not only the strong one:
  - BULK-W tier-X +0.0055 (14/18) / +0.0074 (15/18) / +0.0203 (18/18);
  - BULK-ABC +0.0066 (2/3) / +0.0090 (3/3) / +0.0092 (3/3).
  - Applied to all targets, this is the closed whitening relaxation (097, 099 §5.3, 170, 212): cosine rises and pds pays.
- **A correction to v1.** v1's "≈ 0 on weak/mid" and the strong-agree vs strong-disagree split (+0.0055 vs +0.0005
  fold-X) are reads of a DIFFERENT construction, the pre-only swap (the whole row replaced by `pre`). They are not evidence
  for this candidate.
- **The claim.** The gate is justified ONLY by pds exposure.
  - mse reads the aggregate cosine, which is ENERGY-weighted. The strong tertile holds 0.59 (A/B/C) to 0.73 (W) of the
    realised energy a.
  - pds is COUNT-weighted: each target counts 1/n.
  - Swapping only the strong tertile keeps 47 % (A/B/C cancer: +8.6 % of +18.4 % fold-level composition r) to 80 % (W:
    +33.3 % of +41.9 %) of the global swap's cosine gain. It exposes one third of the targets to pds loss, and these are
    the targets with the largest pds margins (literal pds 0.944 A/B/C, 0.863 W).
  - So the gated form should have a better pds-for-mse exchange than the global one. That is the hypothesis. It is
    decided by comparing it with the global form in the same run.
- **Members.**
  - mse moves through the aggregate cosine. Rows keep their norms, so the energy A, the key, P2g's terciles and the cliff
    do not move (170's mechanism does not apply).
  - pds loses some specificity on the moved rows.
  - fid / jac / reach / nmae are pinned by the unchanged carriers. Ballast K moves only through the render (dK reported).

## 2. Class
**Transfer.** A better bulk direction for strong targets from the same public donors and operator. The gate (z_src
tertile) is computed from the day's donor banks. Nothing is fitted to A/B/C.

## 3. The run: one candidate, three comparators, honest tries (item C)
Code: `price810.py` (da001560…), `cloud/job810p.sh` (ea9e743f…), jobs `cloud/jobs-p.txt` (22 folds), data
`cloud/sync810p.txt`. Part rows are read from `s3://<your-s3-bucket>/output/lab/810parts/` (pinned in
`work/s3parts810.tsv`).

| arm | role | construction |
|---|---|---|
| **sgu1** | **THE candidate** (family `810-sgu`) | top z_src tertile (707's tiering: over targets with a donor row); base' = \|base\| · pre/\|pre\|; fin' = base' + cd4 + h1 + part, rescaled to \|fin\| |
| inc | base of every comparison | L folds: pct+P2g (October). ABC folds: cap6+P2g(abc) (ALLSW's bulk), as ka804 |
| 707 | comparator (a separate family, board read pending) | ka804's 707 row (U blend β 0.25 of the final row, top tertile) |
| **sgu1all** | **known-answer comparator: the CLOSED global form** | the same swap on EVERY target with a donor row |
| sgu0 | known answer, not saved | β 0 on the top tertile, == inc |

- **Known answers** (abort unless |d| < 1e-12, per fold): the parts sum to F.fin; sgu0 == inc (x, a, pds, pdsl, K); inc and
  707 == 804's stored arrays (`experiments/804-known-answers/work/run/bulk/<fold>/`). The engine's own KAs run on load.
- **Tries.** sgu1 is charged **N = 6**, e(6) = 1.28 sd, because β 1 and the tertile gate were picked after
  composition-level reads of six forms on these same 22 folds (`work/foldlevel810.txt`), next to the ablation grid in
  `work/table810.md`. β 0.5 and the agreement gate (P2) are withdrawn and will not be priced. sgu1all, inc and 707 are not
  tries of this family.
- **Surfaces.**
  - BULK-ABC cancer (ABC/k562_gwps, orion_hct116, orion_hek293t) is the deciding surface.
  - BULK-W (18 L folds) is the transfer check.
  - ABC/kolf2 is reported only.

## 4. Expected answers before data

### 4a. sgu1all, the closed global form (item A)
- From 212 (`experiments/212-gamma-cosine/RESULTS.md:26-30`): along γ, pds −0.3517 per unit cos against mse +0.3655, so
  pds eats ≈ 96 % of the mse gain.
- From 097 §9h (`experiments/097-direction/RESULTS.md:1272-1277`): γ 0 nets −0.0033 (h1-test) and −0.0080 (orion);
  gain/spend 0.77 / 0.37, i.e. pds eats 130–270 %.
- **Expected on BULK-ABC cancer:**
  - mse leg ≈ +0.0106 overall. Realised r ≈ 0.81 × 18.4 % ≈ +15 %, at stats800's +0.000715 overall per % r. This
    projection is NOT blind (composition read × 707's realised/composition ratio).
  - The pds leg eats ≥ 95 % of it: E(sgu1all) = −(pds leg)/(mse leg) ≥ 0.95.
  - Net point ≤ +0.0005; range −0.016 … +0.0005.
- **STOP rule:** if the checker reads sgu1all clearly net-positive on BULK-ABC cancer (95 % lower bound > 0, literal
  pds), the checker does not see the pds mechanism that closed γ on this class. Then no read of sgu1 is made. The result
  is reported to the lab lead as "unpriceable on the bulk checker", and the family stops.
- If sgu1all reads net-positive but its lower bound ≤ 0, that is flagged as unexpected and the rule below still applies.
- BULK-W: reported only; no closure was measured on these panels.

### 4b. sgu1
Not blind: the mse legs (composition reads × 707's realised/composition ratio, 0.81 ABC / 0.71 W). Blind: the pds legs,
the render and the ballast realisation.

| surface | mse leg (not blind) | pds leg (blind) | net (literal pds) | E = −pds/mse |
|---|---|---|---|---|
| BULK-ABC cancer | +0.0050 (r ≈ +7.0 %) | −0.0015 … −0.0050 (at 707's E 0.69: −0.0035) | +0.0015 (0 … +0.0035) | 0.3 … 1.0, point 0.7 |
| BULK-W | +0.017 (r ≈ +24 %) | −0.002 … −0.008 | +0.009 … +0.015, ≥ 16/18 folds positive | — |
| kolf2 | +0.0018 (r ≈ +2.5 %) | — | reported only | — |

- Paired on BULK-ABC: sgu1 − sgu1all > 0 (+0.001 … +0.015); sgu1 − 707 +0.0005 … +0.0020.
- With a per-fold bootstrap sd ≈ 0.0009 (707's was ≈ 0.0007), the lower bound after e(6) is expected ≤ 0. The most likely
  outcome under §6 is **HOLD-707**, not DOSSIER.

## 5. Falsifiers (any one CLOSEs the family for October), in the checker's units (BULK-ABC cancer)
- **C1:** sgu1's net point ≤ 0.
- **C2 (falsifier 2, rewritten against sgu1all):** E(sgu1) ≥ E(sgu1all) at the points, with E = −(overall pds leg)/(overall
  mse leg). The gate does not concentrate mse relative to pds. P(E_sgu1 < E_sgu1all) over the paired draws is reported
  beside it.
- **C3:** sgu1 − sgu1all paired point ≤ 0: the gate does not beat the closed global form.
- **C4:** sgu1 − 707 paired point ≤ 0: the base-only construction adds nothing over 707.
- **C5:** sgu1's BULK-W net point ≤ 0.

## 6. Decision rule (item D), applied by `decide810.py` (724e8059…) in this order
- **Units.** Overall = mse leg (518's identity, stats800) + literal pds leg, in board units (stats800.INV). Literal pds
  DECIDES. pdsl (λ 0.5) is reported beside it and is required positive for a DOSSIER (D2).
- **Bootstrap.** Paired target bootstrap with 4,000 draws, seed 810. Every fold is resampled independently, the same
  targets are used for all arms within a draw, and folds carry equal weight.
- `check800.py aggregate --cand sgu1 --base inc` also runs, family `810-sgu`, with candidate.json carrying this file's
  sha256, so the wear ledger records the try. check800 computes N from the ledger (it will say 1); decide810's N = 6 rules.

The rule:
0. **KA gate:** every KA in §3 passes on all 22 folds, or there is no read.
1. **STOP** if sgu1all's BULK-ABC 95 % lower bound > 0 (§4a).
2. **DE-H1 FIRST** if |dK_rel| > 5 % for any arm on any fold. The bulk engine cannot price it.
3. **CLOSE** if any of C1–C5 holds. The family is recorded closed with its mechanism, and October does not use it.
4. **DOSSIER** if none of the above holds and all of these do:
   - D1: sgu1's BULK-ABC lower bound after the e(6) shift is > 0;
   - D2: the pdsl net on BULK-ABC is > 0;
   - D3: BULK-W has ≥ 14/18 folds positive.

   The dossier states the board projection as a range including 0 (702 rule 4) until 707's board read gives this class
   a ρ.
5. **HOLD-707** otherwise: positive and through every falsifier, but not D1–D3. It is resolved by 707's board read
   (10-17/18, primary statistic Δ(pds + mse) overall vs ALLSW, per 707's prereg), call it b707, against the checker's
   prospective 707 read c707 = +0.00044 on BULK-ABC (804 grade):
   - b707 ≤ −0.0005 (HARM): **CLOSE**. On the board, this class's pds cost exceeds what the checker sees, and sgu1 moves
     the same rows further.
   - |b707| < 0.0003, or sign(b707) ≠ sign(c707): no ρ. The projection includes 0, so it is **parked**: WEAK is
     recorded, there is no DOSSIER, and there is no A/B/C slot request.
   - b707 ≥ +0.0003 with the same sign: ρ = c707 / b707. Projected sgu1 = (sgu1's BULK-ABC point after e(6)) / ρ.
     **DOSSIER** if the projected point ≥ +0.0005 and C4 held; otherwise parked.

## 7. Cost
- One EC2 run of 22 folds × 5 arms through engine801, ≈ 3 min per fold and ≤ 6 GB.
- r7i.8xlarge, ~30 min cap, ≈ $1.5; ≤ $5 with one re-run.
- decide810 and check800 aggregate run locally on MBs, in seconds.
- Lane spend so far: $0.73 of $50.

## 8. Day-of path (unchanged from v1; the judge accepted it)
- The parts are recomposed on the day with lib801's calls (decomp810.compose's code, 20–40 s, ≤ 6 GB) from the day's
  bundle, the day's 3-donor operator and the day's CD4 banks. KA: the parts sum to the builder's dump to ≤ 1e-12
  (3.6e-14 measured here).
- z_src comes from 522's def_strength_key522. The swap is price810.sgu_rows at β 1 on the top tertile.
- Nothing is fitted to A/B/C or D/E/F truth. It stays UNAPPLIED until a dossier and the coordinator's decision.

## 9. Novelty against the record (unchanged in substance; see VERDICT-1 check 1)
- **Closed:**
  - a GLOBAL γ/λ change (170, 193, 212: "do not reopen");
  - global whitening relaxation (097, 099 §5.3);
  - learned per-target routing (066, 263, 269).
- **Not closed:** a fixed, strength-gated swap of the base part only, with H1/CD4 and row norms kept. Its case rests
  entirely on beating the closed global form in the same run (C2, C3).
