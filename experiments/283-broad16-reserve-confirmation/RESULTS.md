# 283 — Claim C2: F5 confirmed on the four RESERVED Broad-16 lines (one access, preregistered)

Date: 2026-09-06 (access 08:48–08:58Z). Analysis only: no build, no `vcc`, no `cell-eval2`, no submission, nothing written outside this
directory except the claim and its result in `docs/CONFIRMATION-LEDGER.md`. The four reserved lines' single-cell
members were fetched by CRC-verified HTTP range request into `work/sc/` (2.42 GB, never under `data/`). Peak RSS 4.83 GB
(the UMRC3 build); every analysis step ≤ 4.2 GB. Wall clock ≈ 30 min end to end (prereg 08:46Z, access 08:48Z, all results 08:58Z).

## 0. Headline

**C2a FAIL, C2b PASS, C2c PASS.**

C2a fails **the letter of its frozen rule** (1 of 4 reserved lines clears a median split-half residual cosine > 0.028
with the bootstrap lower bound also > 0.028; the rule required 3). It does **not** fail on the substance: every
preregistered secondary says the line-specific component is there. Against each reserved line's **own** permutation
null the residual is reproducible at +0.0347 [+0.0229, +0.0468], P(≤0) < 0.0001; all four lines sit above their own
OR-control (target-independent) null, which reads the same on the reserve as on the working lines (0.0148 vs 0.0147);
three of the four medians lie inside the working-line range; and the pooled reserved median is 0.0345 [0.0265, 0.0451]
against the working lines' 0.0455 [0.0414, 0.0497] — 76 % of the working-line size. The frozen rule's own calibration,
disclosed before the access, was that it passes on only 8 of 12 **working** lines; 1 of 4 has probability 0.111 under
that rate. The post-hoc reason is visible and measured (§5): **sequencing depth**, on which the reserve drew badly
(median 17,428 UMI/cell vs the working lines' 30,944) — the five shallowest of all 16 lines all fail the rule, three of
them reserved.

## 1. The preregistration, and what was frozen before the access

`PREREGISTRATION.md`, sha256 `65f74d562ee6d4a55d8b93e467338b1fd8622e7723d4f37c6f21cf9ce3808fc0`, written
2026-09-06T08:46:12Z and hashed at 08:47:21Z. It was appended to `docs/CONFIRMATION-LEDGER.md` as claim C2 **before the
first byte of the reserved lines' data was fetched** (first fetch 08:48Z, log `work/sc/fetch283.log`). The reserve
itself was declared by 268 on 2026-09-05T16:20Z (`reserve.json` sha256 `1ce639b5…`) by a name-sort rule, before any
response vector existed. `docs/RESEARCH-PROGRAMME-2026-09.md` rule 2 is satisfied: one claim, one access, logged.

Frozen before looking: the pass rules; the threshold 0.028; the working-line LOLO transfer range [0.0400, 0.0937];
the calibration of the C2a rule on the working lines (8/12); the four secondaries (working-range membership, the
within-line permutation null, the within-line OR-control null, the guide-half version); the bootstrap (4,000 paired
target draws, seed 283); the target set (268's frozen 82); the gene axis rule; and the two known-answer rows.

## 2. Known-answer rows (both before any reserved number was computed)

| KA | what | result |
|---|---|---|
| **KA-1** | 268's `build_halves2.py` (copied byte-identical, sha256 `cd81abca…`) re-run on **working** line `T3M5` from the already-fetched h5, `EXCLUDE_ARM=1`, seed 268 → compare all 9 arrays of `halves/T3M5.npz` with 268's | **0.0e+00** on `full`, `A`, `B`, `v1`, `v2`, `ncell`, `libmean`; `targets`/`genes` identical; NaN patterns identical. 8 s, peak RSS 3.13 GB |
| **KA-2** | `work/lolo283.py` run in KA mode (12 folds over the working lines, train = the other 11) → compare with 268's `lolo268_K2.csv` | **0.0e+00** over **45 numeric columns × 984 rows**, including the `SHUFFLED_labels` arms — i.e. the RNG stream is identical too. Peak RSS 3.17 GB (`work/ka2-diff.txt`) |
| **KA-3** (on the reserve, part of the build) | each reserved line's FULL per-KO pseudobulk must equal the on-disk `pseudobulk_sum_by_ko` column exactly | **100/100 columns exact on each of the four lines**, median relative deviation 0.0 (`work/halves/<line>.json`) |

## 3. The data (the single access)

Figshare 33273600 `single_cell_data.zip` (`https://ndownloader.figshare.com/files/67643568`, 3,171,208,158 B), the four
`*_GE_filtered_feature_bc_matrix` members only — 817 MB compressed at 27–40 MB/s, 2.42 GB on disk, every member
CRC32/size-verified. Cell assignment from the authors' `all_cells_table.csv` already on disk from 268.

| line | batch | lineage | h5 cells | cells kept | KOs | cells/KO median | min | max | control cells | UMI/cell | targets < 6 cells | parse KA |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| JVE127 | apollo34 | Bowel | 44,095 | 11,366 | 100 | 117 | 17 | 195 | 1,171 | **9,981** | 0 | 100/100 |
| MG63 | apollo34 | Bone | 35,429 | 8,780 | 100 | 87 | 4 | 179 | 945 | 19,457 | 1 | 100/100 |
| SLR23 | apollo32 | Kidney | 6,154 | 2,736 | 100 | **27** | 7 | 58 | 311 | 65,265 | 0 | 100/100 |
| UMRC3 | apollo35 | Kidney | 41,388 | 8,884 | 100 | 92 | 11 | 172 | 959 | 15,399 | 0 | 100/100 |

31,766 cells kept. All 82 of 268's frozen targets have ≥ 6 cells on all four lines, so **n = 82 on every line**, no
target dropped. The gene axis came out at **12,199 genes on all four folds** — identical to 268's working-line axis, so
no gene passing the working-line filter is missing from any reserved line. Every fit (gene filter, CPM filter,
centring, basal PCA, ridge, residual programme basis) used the 12 working lines only.

## 4. The three sub-claims

### C2a — reproducibility of the line-specific residual: **FAIL by the frozen rule, confirmed by every secondary**

Statistic: per line, the median over the 82 targets of `cos(r_A, r_B)`, `r = (reserved line's half response) − (12
working lines' mean of the same half)`; paired target bootstrap, 4,000 draws, seed 283.

| line | median | 95 % CI | rule (med > 0.028 **and** lo > 0.028) | in working range [0.0211, 0.0652] | frac. of targets above the line's own permutation null p95 | the line's own OR-control median | guide-half median |
|---|---:|---|---|---|---:|---:|---:|
| JVE127 | 0.0173 | [+0.0131, +0.0265] | no | **NO** (below) | **0.378** | +0.0051 | 0.0120 |
| MG63 | 0.0502 | [+0.0293, +0.0579] | **YES** | yes | 0.232 | +0.0129 | 0.0161 |
| SLR23 | 0.0563 | [+0.0190, +0.0794] | no (CI) | yes | 0.146 | +0.0315 | 0.0305 |
| UMRC3 | 0.0362 | [+0.0273, +0.0495] | no (CI, by 0.0007) | yes | 0.122 | +0.0168 | 0.0141 |
| *working lines (12), for comparison* | *0.0211–0.0652* | — | *8 of 12 pass* | — | *0.146–0.512* | *median +0.0147* | *median 0.0341 pooled* |

**1 of 4 → C2a FAIL** (the rule needed ≥ 3). Under the disclosed working-line rate 8/12, P(≤ 1 of 4) = 0.111.

The preregistered secondaries, all pointing the other way:

- **S2, the within-line permutation null** (each target's half-A residual paired with a permuted target's half-B
  residual on the same line): paired over targets, the reserved lines' residual cosine exceeds their own null by
  **+0.0347 [+0.0229, +0.0468], P(≤ 0) = 0.0000**. Per line, the fraction of targets above that null's 95th percentile
  is 0.378 / 0.232 / 0.146 / 0.122 against 0.05 expected — 2.4× to 7.6× the null, three of the four inside the working
  lines' 0.146–0.512.
- **S3, the within-line OR-control null** (each of the 10 OR knockouts against the other nine, same pipeline — the
  target-independent artefacts a KO-vs-controls comparison carries): reserved medians +0.0051 / +0.0129 / +0.0315 /
  +0.0168, median **0.0148**, against the working lines' **0.0147**. The null is the same on the reserve. Every
  reserved line's target residual is above its own OR-control value (3.4× / 3.9× / 1.8× / 2.2×); pooled, 0.0345 vs
  0.0148 = **2.3×**, against 3.1× on the working lines (268's headline 4× was an energy-fraction ratio, 0.115/0.028).
- **S1, size**: pooled median **0.0345 [0.0265, 0.0451]** (n = 328) vs the working lines' **0.0455 [0.0414, 0.0497]**
  (n = 984) — 76 %, intervals overlapping. Three of four line medians inside the working range; JVE127 below it.
- **S4, the guide unit**: reserved pooled guide-half median 0.0157 [0.0097, 0.0233], **0.45×** the cell-half value,
  against 0.75× on the working lines. F5's "79 % survives across guides" (an energy fraction) is **weaker** on the
  reserve in cosine terms; the deficit runs in the same direction as, and is explained by, the depth deficit (§5:
  Spearman(UMI/cell, guide-half cosine) = +0.644 over 16 lines).

### C2b — the 12-line shared response transfers: **PASS**

Statistic: per line, the median over targets of `cos(shared response, reserved line's full response)`. Frozen range
from 268's working-line LOLO folds, stated before looking: **[0.0400, 0.0937]**.

| line | median | 95 % CI | in range | mean | retrieval mean |
|---|---:|---|---|---:|---:|
| JVE127 | 0.0313 | [+0.0228, +0.0365] | **NO** | 0.0485 | 0.684 |
| MG63 | 0.0595 | [+0.0398, +0.0817] | yes | 0.0999 | 0.768 |
| SLR23 | 0.0543 | [+0.0354, +0.0972] | yes | 0.0914 | 0.701 |
| UMRC3 | 0.0493 | [+0.0359, +0.0705] | yes | 0.0903 | 0.743 |
| *working lines* | *0.0400–0.0937* | — | — | *0.0741–0.1450* | *0.700–0.819* |

**3 of 4 → C2b PASS.** The one line outside is JVE127, the shallowest line in the whole 16 (§5), and it is *below* the
range, not above: the shared response transfers to the reserved lines at the same size it transfers between working
lines. Retrieval on the reserve (0.684–0.768) sits at the bottom of the working lines' band (0.700–0.819), the same
ordering. Declared before the access and not corrected: the reserved lines' shared response is a mean of **12** lines
where the working folds used 11, a bias in the reserve's favour.

### C2c — basal state does not predict the residual: **PASS** (268's negative reproduces)

268's model, fitted on the 12 working lines only (K = 2 basal programmes, per-gene ridge λ = 1 for c1; J = 5 residual
programme basis for c2), predicting the reserved lines' residuals; baseline (a) = the shared response alone. Paired
target bootstrap pooled over the four lines (each target carries its four folds).

| arm | metric | reserved lines | 268's 12 working lines |
|---|---|---|---|
| c1 target × basal, ridge | Δcos(full response) | **−0.0050 [−0.0077, −0.0022]** | −0.0071 [−0.0093, −0.0050] |
| c1 | Δretrieval | **−0.0033 [−0.0168, +0.0103]** | +0.0003 [−0.0097, +0.0110] |
| c2 target × basal, programmes | Δcos(full response) | **+0.0000 [−0.0026, +0.0029]** | +0.0000 [−0.0016, +0.0017] |
| c2 | Δretrieval | **−0.0047 [−0.0178, +0.0086]** | −0.0010 [−0.0109, +0.0099] |

**No interval lies above 0 → C2c PASS**: on data never used for fitting or selection, conditioning a target's response
on the recipient's basal programmes buys nothing over the shared response — c1 is harmful, c2 is exactly zero, both
retrieval readings are negative with intervals across zero. Every number reproduces 268's working-line value inside its
interval.

The preregistered secondary — the "detectable but not usable" half of F5 — also reproduces, at about half the size:

| arm | cos(predicted residual, reserved residual) | WRONG recipient | SHUFFLED labels | minus WRONG |
|---|---|---|---|---|
| c1 | +0.0084 [+0.0006, +0.0164] | +0.0003 [−0.0004, +0.0010] | +0.0036 [+0.0010, +0.0061] | +0.0081 [−0.0003, +0.0167], P = 0.030 |
| c2 | +0.0101 [+0.0000, +0.0205] | +0.0003 [−0.0005, +0.0011] | −0.0007 [−0.0044, +0.0028] | +0.0098 [−0.0010, +0.0207], P = 0.041 |

(268: c1 +0.0181, c2 +0.0174, wrong +0.0012/+0.0007.) The direction is right and the wrong-recipient control is at
zero, but on 4 folds instead of 12 the "minus WRONG" intervals brush zero — the effect is confirmed in sign, not
sharply in size.

For the record, the other arms on the reserve: the recipient-only offset Δcos −0.0010 [−0.0021, +0.0000], Δretrieval
−0.0083 [−0.0163, −0.0005]; the nearest-basal-line rule Δcos −0.0379, Δretrieval −0.0643 (net harmful, as in 268); the
**oracle** line offset (uses the held line's own data) Δcos **+0.1935 [+0.1804, +0.2058]** but Δretrieval **−0.0345
[−0.0634, −0.0071]** — 268's signature exactly: a target-independent offset is real and large, and is the one component
a `pds`-type metric cannot reward.

## 5. Post-hoc (NOT preregistered): why the frozen C2a rule failed — depth

Computed from per-line summaries already produced under the preregistration; no further access to the reserve. Over all
16 lines, **Spearman(UMI per cell, median residual cosine) = +0.650, p = 0.0064**; cells per KO is uninformative
(−0.182, p = 0.50); the same depth correlation holds for the guide-half cosine (+0.644, p = 0.0071) and is weaker for
the transfer cosine (+0.453, p = 0.078). `work/line_table283.csv`:

- The **five shallowest lines of the sixteen — JVE127 (9,981), SNU8 (11,723), KMRC20 (12,611), UMRC3 (15,399), C4I
  (15,507) — all fail the 0.028 rule**, and three of those five are working lines that failed it in 268's own data.
- Below 21,000 UMI/cell: 2 of 7 lines pass, median residual cosine 0.0315. At or above: 7 of 9 pass, median 0.0555.
- The name-sort reserve rule drew a **depth-poor quarter**: reserved median 17,428 UMI/cell vs working median 30,944;
  JVE127 is the shallowest line in the release and UMRC3 the fourth shallowest.
- The fourth failure, **SLR23**, fails on interval width, not on level: its median 0.0563 is the fourth highest of the
  sixteen, but it has the **fewest cells per KO of any line (27 median, 6,154 cells total)**, so its bootstrap lower
  bound is 0.0190.

So of the three reserved-line failures, two are the depth tail that 268's own working lines already contained and one
is a cell-count tail. This is a reading of the *rule*, not a rescue of the claim: C2a's frozen letter failed and is
recorded as failed.

## 6. What this does and does not establish

**Establishes**, on four cell lines never used for fitting or selection, chosen by a rule fixed before any residual
existed:

1. The line-specific response component of F5 is **present and reproducible** on the reserve — above each line's own
   permutation null (+0.0347, P < 0.0001) and above each line's own OR-control null (4/4 lines) — at **76 %** of the
   working-line size, with the null itself reading identically on both sets (0.0148 vs 0.0147).
2. The 12-line shared response **transfers to unseen lines at the size measured between seen lines** (3/4 inside
   [0.0400, 0.0937]) — the strongest of the three results and the one with the cleanest pass.
3. F5's negative — **basal state does not usefully predict the line-specific residual** — reproduces on held data with
   the same signs and sizes (c1 −0.0050, c2 +0.0000 in Δcos; both Δretrieval negative with intervals across 0). This is
   the confirmation the programme owed: the closure in `docs/FINDINGS-2026-09.md` F5 is not an artefact of the 12 lines
   it was fitted on.
4. **How much of it is measurable is set by depth**, not by cells per condition (ρ +0.650 vs −0.182 over 16 lines) —
   consistent with F10's finding that response magnitude, not cell count, predicts reproducibility, and adding
   UMI depth as the second axis.

**Does not establish**: that the component is biology rather than batch (268 §7 limit 1 stands — every reserved line's
batch also contains a working line, so the *comparison* is not confounded with an unseen batch, but a line measured in
two experiments is still what would separate biology from batch); anything about the 300-target VCC panel (overlap with
the 90 Broad targets is still 0); or any competition-track quantity. Nothing here is buildable and nothing was built.

## 7. Files

- `PREREGISTRATION.md` (sha256 `65f74d562ee6d4a55d8b93e467338b1fd8622e7723d4f37c6f21cf9ce3808fc0`), the ledger's claim
  C2 and its result.
- `work/prereg-calibration.txt` — the rules' working-line calibration, computed before the access.
- `work/lib283.py` (268's `load268.py` with the reserved-line assertion lifted inside this experiment, plus the
  OR-control response builder), `work/lolo283.py` (268's LOLO in KA / RESERVE modes), `work/stats283.py`,
  `work/build_halves2.py` (268's, byte-identical, sha256 `cd81abca…`), `work/fetchzip.py`, `work/zipdir.py`.
- Data: `work/sc/` (four h5 members, `fetch283.log`), `work/halves/{JVE127,MG63,SLR23,UMRC3,T3M5}.npz`
  (`231deee7…`, `f65bf72e…`, `64eb45c4…`, `c4f706d9…`, `7c2afcfe…`).
- Results: `work/lolo283_RESERVE.csv` (`4a368bbe…`, 328 rows), `work/lolo283_KA.csv` (`89d6abe1…`, 984 rows),
  `work/lolo283_{KA,RESERVE}_folds.json`, `work/c2_results.json` (`03de57ff…`), `work/line_table283.csv`,
  logs `work/ka2.log`, `work/ka2-diff.txt`, `work/reserve.log`, `work/stats283.log`, `work/posthoc-depth.log`.

**The reserve is now spent.** The four lines have been accessed once, for claim C2, and are not to be used again.

---

**C2a FAIL, C2b PASS, C2c PASS.**
