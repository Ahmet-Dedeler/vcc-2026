# 446bja BOARD RESULT — 0.2689767 (−0.0010661 vs 446b): context A's same-line programme (Nadig Jurkat) does NOT transfer through fid

Entry `AsHteaDAXCWQo57fLYti`, submitted 2026-09-26 00:38Z by `446-q175-level/arm0926_s2b.sh` ("the vortex mixer wrote a
sonnet", description empty) after 446b's capture. Scored and captured 01:03Z. Package `446bja-…-mde010coarse-zr5.vcc`
sha256 `512a7e41…`. Reference 446b `UDgjOKL9gI1udAexq6ue` 0.2700428. Bins: `PREREGISTRATION-446bja.md`. Reader:
`446-q175-level/read0926.py`.

| member | Δ raw vs 446b | Δ scaled | Δ overall | predicted (A only) |
|---|---:|---:|---:|---|
| pds | 0.000000 | 0 | 0 | pinned (falsifier: pass) |
| mse | 0.000000 | 0 | 0 | pinned (falsifier: pass) |
| nmae | **−0.001026** | +0.001411 | +0.000235 | −0.0015 … −0.0060 |
| fid | **+0.000219** | +0.000669 | +0.000111 | **+0.0015 … +0.0040** |
| reach | −0.007891 | −0.008470 | −0.001412 | −0.010 … +0.010 (unpriced) |
| jac | −0.000002 | −0.000006 | −0.000001 | ≤ 0.0005 (falsifier: pass) |
| **overall** | | | **−0.0010661** | +0.0005 … +0.0030 |

## Against the bins
- **fid (primary): no transfer** (|Δ| < 0.0005). The H1 price, taken to one context at 392's transfer rate, was about
  +0.002. The board gave +0.0002, about a tenth.
- nmae: neutral bin (|Δ| < 0.0015), in the predicted direction (−0.0010).
- reach −0.0079 raw, the one member that moved, in the direction H1 showed. On H1 the programme cost reach on four
  small-gate targets, where gate sign flips are coin flips.
- **Overall: harm bin (≤ −0.0010), by 0.00007.** Nearly all of it is reach, which is the noisiest member (ordering
  channel sd ≈ 0.0045 overall). The fair reading is **no gain**, and the prereg's consequence is the same either way.

## Why the programme did not carry here, when it did on H1 (candidates, none tested yet)
1. **CLEAN4 already carries a quarter Jurkat.** Nadig Jurkat is one of CLEAN4's four banks, and the two priors differ on
   only ~25 % of A's block genes. The marginal information is smaller than Pan's over CLEAN4 on H1, which is not in
   CLEAN4 (a 27 % difference).
2. **Within-Jurkat agreement was weak even internally.** 370's within-Jurkat split read 0.577 against the cross-line
   ~0.55 on A's block, a smaller same-line edge than H1's 0.625 vs 0.581.
3. **Screen and state.** Nadig's Jurkat is an essential-gene screen (2,270 perturbations, 10x 3′). On H1, Pan's
   essential-only programme still read 0.627, so essential-only is not by itself the cause. Culture state and lab
   differences may matter more for a T-cell line.

## Consequences
- **October: the same-line mode keeps `noprog`** (programme off; CLEAN4 stays the default block prior). The flags stay
  `pt:noprog`, and `def_day.sh`'s default is unchanged. Only the per-target rows (and their per-target prior) switch
  on for a matched context.
- The rows part, the large H1 lever (+0.093 at 35 % coverage), is not tested by this entry. A/B/C has no bank that
  covers its panel. Its transfer rests on H1 alone, and this read is a caution: the H1-priced programme fid gain did
  not reach the board. The rows' H1 evidence is much larger and more direct (per-target gate signs 0.764 vs 0.551, 26/26
  targets; centred direction 0.260 vs 0.069, 30/30). The programme's H1 net was +0.0001, i.e. already nothing.
- The 09-27 rule is **Branch 1b** (446bja − 446b < +0.0010). The base stays 446b: s1 = 446bj2 (LODO key), s2 = 446br2
  (floor −2). Both were built on 09-25, and `chain0926.sh` armed them at 01:04Z.

## Addendum 01:20Z — candidate 1 (dilution in CLEAN4) does NOT explain the null (`weightcurve450.py`, H1)
Setup: on H1, the sign of (1 − w)·z(CLEAN4 P_down) + w·z(Pan P_down) over the 6,007 genes both priors cover, scored on
h1-test.

| w (same-line weight) | 0 | 0.1 | 0.2 | **0.25** | 0.33 | 0.5 | 0.67 | 0.8 | **1** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| gate-block sign accuracy | 0.5811 | 0.5874 | 0.5948 | **0.5972** | 0.6037 | 0.6094 | 0.6185 | 0.6234 | **0.6239** |

- At w 0.25, about Jurkat's weight inside CLEAN4 (its row share is ~14 %, so if anything less), the consensus holds
  only **37 %** of the w 0 → 1 gain on gate-block genes (49 % on all gate genes).
- The H1 analogue therefore predicts that going from CLEAN4 to Nadig alone would still gain about +0.027 gate-block
  accuracy. That is exactly the move 446bja made on A, and the board saw fid +0.0002.
- **So the null is about the bank, not dilution.** Nadig 2025's Jurkat programme does not carry Arc's context-A
  programme the way Pan's H1 programme carries VCC H1's. Candidates 2 and 3 remain: lab, culture or screen differences,
  which matter more for a suspension T-cell line than for an ESC line.
- For October, the programme part stays off. Whether the rows part carries a line from another lab is the open
  question for any matched D/E/F context. H1 (Pan: rows +0.093) is still the only scorer-level evidence.
