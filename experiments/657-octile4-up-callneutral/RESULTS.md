# 657 — CLEAN4 octile 4 delivered UP, call-neutrally: every route passes the no-harm gates on 618, and B is the board route by the rule. A post-hoc check says the cancer-line case was a KO-data offset, so B goes last

> **COORDINATOR DECISION (10-04, FINDINGS §2): 657 is CLOSED.**
> - Every route is negative on 618, and calib657's Broad-offset diagnosis is the reading taken.
> - A rider whose expected value is negative gets no slot, not even an empty one.
> - So B's eligibility below is moot. No 657 arm will be built.

Runs 2026-10-03 21:51Z – 10-04 02:16Z (from the chain logs), council session [7772dc], assigned by the coordinator after
625 Amendment 4. PREREGISTRATION ed36f45f and the harness (5496c2a7, 999ec545) were committed before any result.
Written 10-04 ~06:30Z.

**Class:** transfer, as attempted. On 618 this measures no-harm only; the value is unmeasured. Nothing here is fitted to
A/B/C.

**Tries:** octile 4 is the best of 8 octiles from 625's scan. Here: 3 routes × 2 stacks (630, 630+DN) × 3 seeds on 618,
and no board arm.

## Result (618: h1-trainval, faithful key, zero repair −5, zs; overall-eq = 192's INV ÷ 6 on the DE members; paired per seed)
| route vs base | f 0.35 mean (seeds 0 / 1 / 2) | f 0.5 mean (seeds) | Δcalls / target | members at f 0.35 [fid jac reach nmae] | verdict |
|---|---:|---:|---:|---|---|
| **A** hook UP, + 630 | −0.00044 (−0.00059 / −0.00036 / −0.00036) | −0.00028 (−0.00036 / −0.00029 / −0.00019) | +3.6 (3.4 / 3.3 / 4.2) | −0.00022 −0.00015 −0.00065 +0.00044 | ELIGIBLE |
| **B** 630 selects octile 4 | **−0.00011** (−0.00013 / −0.00013 / −0.00007) | −0.00013 (−0.00017 / −0.00015 / −0.00006) | **0** (0 / 0 / 0) | 0 0 0 +0.00040 | **ELIGIBLE** |
| **AB** both | −0.00050 (−0.00065 / −0.00035 / −0.00049) | −0.00032 | +3.6 | −0.00022 −0.00015 −0.00065 +0.00066 | ELIGIBLE (0.000003 inside the bound) |
| A0 hook only, no 630 (component) | −0.00030 | −0.00009 | +3.6 | −0.00022 −0.00015 −0.00065 −0.00005 | (component) |
| A + DN vs 630 + DN | −0.00044 | −0.00029 | +3.6 | −0.00022 −0.00015 −0.00064 +0.00046 | ELIGIBLE |
| B + DN vs 630 + DN | −0.00011 | −0.00013 | 0 | 0 0 0 +0.00039 | ELIGIBLE |
| AB + DN vs 630 + DN | −0.00050 | −0.00032 | +3.6 | −0.00022 −0.00015 −0.00064 +0.00068 | ELIGIBLE (0.000001 inside) |

- **Gates (prereg), all routes:** |Δcalls| ≤ 10 per target on every seed (max 4.2); |Δfid|, |Δjac| ≤ 0.0005 at f 0.35 and
  0.5 (max 0.00022); overall-eq ≥ −0.0005 at f 0.35 and 0.5. Every route is ELIGIBLE. AB passes only by 0.000003 (0.000001
  with DN).
- **The bulk is unchanged in every route.**
  - Every 630 stack keeps each per-gene sum equal to its emission's, on every seed.
  - Route A's emission equals the base emission's sums except 4 gene × target cells on seed 2, by 1 count each (the
    zero-repair floor).
  - So pds and mse do not move, and the DE-only overall-eq above is the whole price.
- **Board route, by the prereg rule:** "if both A and B are eligible, the larger f 0.35 mean goes to the board". That is
  **B**: −0.00011, vs A −0.00044 and AB −0.00050.
- **B's board price if A/B/C behave like 618:**
  - B is nmae-only: calls 0, fid / jac / reach exactly 0, sums identical.
  - With 619's nmae ×1.0–1.2, that is −0.00011 … −0.00013 overall.
  - The read would be exact if B rides the same build as the 630 arm, differing only in `_o4`.
- **DN changes nothing.** Each route's Δ with DN equals its Δ without DN to within 0.00002. DN works on called DOWN genes,
  which are disjoint from both routes.

## What the routes did to octile 4 (uncalled, non-own cells; seed means; `analyze657`, `zero657`)
| arm | sign agreement with the truth (f 1 pool) | mean delivered lfc | mean \|lfc\| (all / in the reference's call set) |
|---|---:|---:|---|
| base (c4h + 630) | 0.580 | +0.020 | 0.157 / 0.162 |
| A | 0.538 | +0.108 | 0.139 / 0.132 |
| B | **0.580** (unchanged) | +0.039 | 0.137 / 0.142 |
| AB | 0.538 | +0.111 | 0.136 / 0.129 |
| emission c4h → A0 (no 630) | 0.580 → 0.538 | +0.001 → +0.073 | 0.132 → 0.096 |

- **What A does:** it overrides the model's own delivered sign, which is right 58 % of the time on these cells, with UP.
  CLEAN4's DOWN sign is right 52.4 % on 618 (so UP is right 47.6 %), and agreement falls to 0.538.
- **What B does:** it flips no sign. It shrinks DOWN-delivered octile-4 genes toward 0: mean |lfc| 0.157 → 0.137, mean
  signed +0.020 → +0.039, agreement unchanged.
  - This confirms the prereg's mechanism: 630 cannot raise a DOWN-delivered gene past 0 without recipients above the
    control max.
  - So **even where UP is right, B's gain is bounded.** It can only remove DOWN error, never deliver UP.
- **Zero totals:** 0 non-own and 0 own zero-total cells in every arm on every seed (`zero657`).

## Prediction check (prereg)
- A, |Δcalls| ≤ 10: **holds**, +3.6.
- A, Δnmae in 0 … +0.0010: **holds** at f 0.35 / 0.5 (+0.00044 / +0.00037). It edges over at f 1 (+0.00105).
- A, overall-eq in −0.0003 … 0: **misses at f 0.35** (−0.00044). The +3.6 calls cost fid / jac / reach, which I had put at
  ~0. It holds at f 0.5 (−0.00028).
- A, "630 amplifies the harm": **holds**. Without 630, nmae is −0.00005; with it, +0.00044.
- B, "smaller than A, Δcalls 0": **holds**.

## Post-hoc (NOT preregistered; `calib657.py`, `work/calib657.{out,json}`): the octile-4 case is a Broad KO offset, not cancer biology
The prereg rested 657's value on Broad-16 (12 KO cancer lines: CLEAN4's DOWN sign wrong on 12/12 at octile 4). Reading the
same evidence across all octiles and banks:

| CLEAN4 octile | 0 | 1 | 2 | 3* | 4 | 5 | 6 | 7 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| UP right, Broad-16 KO (mean of 12, gate cells) | 0.847 | 0.768 | 0.743 | 0.583 | **0.601** | 0.517 | 0.403 | 0.199 |
| UP right, H1 (CRISPRi Flex = Arc's chemistry) | 0.623 | 0.584 | 0.564 | 0.515 | **0.498** | 0.455 | 0.406 | 0.367 |

*Octile 3 straddles P_down 0.5; its entry is CLEAN4-sign agreement.

- **The Broad curve is shifted.** It crosses 0.5 between octiles 5 and 6 (P_down ≈ 0.55); H1's crosses at octile 4
  (≈ 0.51). The whole Broad curve sits about 1.5 octiles to the UP side.
- **Octile 4 is not special.** Its 0.601 sits beside octile 3's 0.583. It is simply the first DOWN-signed octile below
  Broad's shifted crossing, so a global offset flips exactly this octile on every line: 12/12.
- **Every CRISPRi bank is neutral or leans DOWN at octile 4.** Gene level (371's P_banks; each bank read against a prior
  that excludes it), fraction of genes leaning UP in the octile-4 band:

  | bank | prior it is read against | fraction leaning UP |
  |---|---|---:|
  | K562-GWPS | LOLO-3 | 0.41 |
  | K562-ess | LOLO-3 | 0.50 |
  | Jurkat (context A's line) | LOLO-3 | 0.51 |
  | HepG2 | LOLO-3 | 0.26 |
  | H1 2025 | CLEAN4 | 0.48 |
  | Orion HCT116 | CLEAN4 | 0.46 |
  | Orion HEK293T | CLEAN4 | 0.45 |

  - None shows Broad's ~0.60.
  - K562, Jurkat, HepG2 and HCT116 are cancer lines, so "cancer line" does not carry the offset. The Broad lab or the KO
    modality does.
- **Caveats:**
  - The gene-level P_down and gate-cell agreement are different statistics.
  - LOLO bands regress toward each bank's own mean.
  - HEK293T's curve is flat and inverted (616's Orion audit), and HepG2 carries a global DOWN offset (mean 0.534).
  - None of this changes the direction: no CRISPRi bank leans UP at octile 4.
- **Consequence for 657:** the best available prior for Arc's CRISPRi Flex lines is H1 and the CRISPRi banks, i.e.
  octile 4 ≈ coin flip or slightly DOWN. B's expected A/B/C Δ is then about 618's, ≈ −0.0001. Its read would mainly test
  whether Arc's lines carry the Broad offset. That is a real question, but a small one.
- **Wider caution (for the owners of Broad-16-based evidence):** Broad-16 gate cells carry a global UP offset of about
  +0.08 (mean UP-right over octiles: 0.58 vs H1 0.50).
  - Sign evidence learned or validated on Broad-16 near P_down ≈ 0.5 can be the offset rather than transferable biology.
  - Examples: 625's 8-octile transfer simulation (+0.0165, already dropped), and the 622 / 639 identity-track gains where
    they lean on Broad lines.
  - The check: re-read with each line's global sign offset removed (re-centre each line's UP-right curve at P_down 0.5)
    and see what survives.

## Decision / recommendation
- **Verdict (prereg):** A, B and AB are ELIGIBLE. The board route is **B** (CLEAN4 octile 4 UP via 630's selection).
- **Recommendation:** B ranks LAST, in an otherwise empty slot only.
  - Build it as a paired rider on the 630 arm, the same build plus `_o4`, so its nmae-only read is exact.
  - Expected Δ ≈ −0.0001 by the post-hoc calibration evidence.
  - Nothing from 657 enters the October recipe without that board read AND a held-out read on a CRISPRi cancer surface.
- **Route A / AB:** not recommended. They cost calls (+3.6) and reach, and pass the downside bound only narrowly (AB by
  0.000003).

## Known-answer checks
- **KA-E** (seeds 0, 1, 2): 657's `c4h` emission == 618's `hk_c4h` stats bit for bit (q, p, lfc, sums; pre-zs and zs).
- **KA-R** (seed 0): the original rnup630 (7ab0447e, `rnup657orig.py`) == the copy `rnup630_657.py` without `_o4`, on
  657's own c4h dumps (q, p, lfc, sums). 618's own stk_base dumps had been cleaned, so the KA was reworked onto these dumps
  (999ec545). The first attempt failed on the missing dumps: `work/rnup-kar-s0.log`, `work/seq657.v1-failed.log`.
- Route A's override is recorded in every emission log (`o4_n_gated` 775, `o4_was_hooked` 0): today's hook never touched
  octile 4, as the prereg said.

## Files
- **Code:**
  - `run657e.py` (618's run618e + mode `tvc4h_o4`; diff `work/run657e.diff`);
  - `rnup630_657.py` (rnup630 + the `_o4` suffix; diff `work/rnup630_657.diff`), `rnup657.py`, `rnup657orig.py`;
  - `ka657.py`, `seq657.sh`, `seq657b.sh`, `analyze657.py`, `zero657.py`, `calib657.py`.
- **Outputs:**
  - `work/analyze657.{out,json}`, `work/zero657.{out,json}`, `work/calib657.{out,json}`;
  - chain logs `work/seq657{,b}.log`, per-job logs `work/{emit,rnup,rnupdn,rnuporig}-*.log`;
  - diagnostics `work/{c4h,c4h_o4,stacks,stacks_orig}/*.json`.
  - Stats npz files (1.1 GB) and dumps (3.6 GB) stay local.
- **Peak RSS** (`/usr/bin/time -l`): emissions 7.8–8.3 GiB, 630 stacks 3.8–4.9 GiB. All through gate602, before the
  10-04 05:00Z blackout. The post-hoc reads are < 0.2 GB.
