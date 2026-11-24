# 374 — Jurkat pays no retrieval cost up to β≈0.5; recommend β=0.35 (zero pds price, 64% of sign gain)

Status: **COMPLETE**. Numerical, six-line source frame only. No emission, no
323, no board forecast beyond the stated pds translation. DE-member
translation unknown.

## Known-answer reproduction (gate PASSED, own code, saved 341 predictions)

Recomputed from `experiments/341-known-line-calibration/work/<line>/predictions.npz`
with 341/334's own imported functions (`core341.sign_endpoints`,
`core334.endpoints`), all 1,024 targets, vs `work/results.json`:

| check | mine | record | diff |
|---|---|---:|---:|
| same_line_full / weighted_sign (6-line macro) | 0.621947489881444 | 0.621947489881444 | +1.1e-16 |
| same_line_full / retrieval | 0.746316014224910 | 0.746316014224910 | 0.0 |
| same_line_full / cosine | 0.268227504773924 | 0.268227504773924 | 0.0 |
| same_line_full / top-64 | 0.682760716529909 | 0.682760716529909 | 0.0 |
| crossline_diag / weighted_sign | 0.589397759427262 | 0.589397759427262 | 0.0 |
| crossline_diag / retrieval | 0.794613855286738 | 0.794613855286738 | 0.0 |
| crossline_diag / cosine | 0.198734908480320 | 0.198734908480320 | +2.8e-17 |
| crossline_diag / top-64 | 0.642516695930019 | 0.642516695930019 | 0.0 |
| Jurkat Δ wsign vs diag | +0.03380218 | +0.03380218 | <1e-8 |
| Jurkat Δ wsign vs full | +0.02800431 | +0.02800431 | <1e-8 |
| pooled Δ retrieval vs diag | −0.04829784 | −0.04829784 | <1e-8 |

Blend identities: β=0 vs `crossline_diag` max abs diff **1.78e-15**
(PLAN said bit-exact; in fact one ulp from the divide–multiply round trip —
tolerance <1e-12 used, stated here). Assignment's alternate wording
(scale-ridge-first then blend) agrees with 346's unit-blend form to
**2.7e-15**. β=1 row norms equal the diagonal's to 1.8e-15.

## What was measured (inputs, estimator, controls)

- Frame: 341's `selection.json` — 6 lines, 1,024 targets, 512 coordinates.
  Per-line saved `truth`, `crossline_diag` (5-donor mean × diagonal gain),
  `same_line_full` (ridge on the held-out line's own training half),
  both target halves assembled.
- Blend (346's fixed-norm rule, per target row): `m = (1−β)·d/‖d‖ +
  β·p/‖p‖`, `pred = m·(‖d‖/‖m‖)`, β ∈ {0, 0.1, 0.25, 0.35, 0.5, 0.75, 1}.
- Endpoints per β: weighted sign + top-64 (`core341.sign_endpoints`),
  retrieval + own cosine (`core334.endpoints`, midrank/(n−1) over 1,024).
- Intervals: `core341.interval` (seed 3410909, 10,000 draws) on Δ vs β=0 —
  target 95% CI reported (crossed CI also computed for pooled5; Jurkat-only
  crossed == target by construction with one line).
- Translation only: expected pds raw Δ ≈ retrieval Δ × **1.009** (346's
  H1-measured slope). No board conversion, no DE claims.
- Control: saved `labels_shuffled_full` on Jurkat reads wsign 0.5653 /
  top-64 0.6055 / retrieval 0.6372 / cosine 0.1361 — far below every blended
  β, so the curve is not a generic ridge artefact.

## Results — Jurkat (nadig_jurkat, 1,024 targets)

Means; deltas vs β=0 with paired-target 95% CI. pds price = retrieval
Δ × 1.009 (raw pds units, translation, not a score claim).

| β | wsign (Δ [CI]) | top-64 (Δ [CI]) | retrieval (Δ [CI]) | own cosine (Δ [CI]) | pds raw price |
|---:|---|---|---|---|---:|
| 0 | 0.591160 | 0.643790 | 0.819913 | 0.202936 | 0 |
| 0.1 | 0.597851 (+0.00669 [+0.00599,+0.00741]) | 0.653031 (+0.00924 [+0.00778,+0.01075]) | 0.823836 (+0.00392 [+0.00248,+0.00541]) | 0.216698 (+0.01376 [+0.01287,+0.01468]) | **+0.00396** |
| 0.25 | 0.607287 (+0.01613 [+0.01473,+0.01752]) | 0.664986 (+0.02120 [+0.01864,+0.02374]) | 0.825757 (+0.00584 [+0.00242,+0.00933]) | 0.236243 (+0.03331 [+0.03103,+0.03567]) | **+0.00590** |
| 0.35 | 0.612719 (+0.02156 [+0.01979,+0.02342]) | 0.671213 (+0.02742 [+0.02428,+0.03063]) | 0.824980 (+0.00507 [+0.00038,+0.00985]) | 0.247776 (+0.04484 [+0.04164,+0.04819]) | **+0.00511** |
| 0.5 | 0.619462 (+0.02830 [+0.02590,+0.03079]) | 0.678466 (+0.03468 [+0.03078,+0.03866]) | 0.821598 (+0.00169 [−0.00490,+0.00839]) | 0.261530 (+0.05859 [+0.05408,+0.06332]) | +0.00170 |
| 0.75 | 0.625118 (+0.03396 [+0.03081,+0.03726]) | 0.685495 (+0.04170 [+0.03666,+0.04672]) | 0.810248 (−0.00967 [−0.01885,−0.00049]) | 0.272417 (+0.06948 [+0.06324,+0.07610]) | −0.00975 |
| 1.0 | 0.624963 (+0.03380 [+0.03014,+0.03762]) | 0.684792 (+0.04100 [+0.03535,+0.04675]) | 0.795811 (−0.02410 [−0.03532,−0.01290]) | 0.269207 (+0.06627 [+0.05889,+0.07398]) | −0.02432 |

Share of Jurkat's β=1 sign gain captured: 0.1 → 20%, 0.25 → 48%,
0.35 → 64%, 0.5 → 84%, 0.75 → 100% (sign peaks at 0.75, not 1.0).
Retrieval peaks at β=0.25.

## Results — other five lines pooled (5×1,024 targets)

| β | wsign (Δ [CI]) | top-64 (Δ) | retrieval (Δ [CI]) | own cosine (Δ) | pds raw price |
|---:|---|---|---|---|---:|
| 0 | 0.589045 | 0.642262 | 0.789554 | 0.197895 | 0 |
| 0.1 | +0.00674 [+0.00630,+0.00719] | +0.01066 | −0.00163 [−0.00258,−0.00069] | +0.01423 | −0.00164 |
| 0.25 | +0.01597 [+0.01499,+0.01695] | +0.02432 | −0.00599 [−0.00827,−0.00377] | +0.03448 | −0.00604 |
| 0.35 | +0.02106 [+0.01977,+0.02238] | +0.03054 | −0.01067 [−0.01374,−0.00768] | +0.04646 | −0.01077 |
| 0.5 | +0.02713 [+0.02543,+0.02890] | +0.03736 | −0.01939 [−0.02358,−0.01536] | +0.06085 | −0.01956 |
| 0.75 | +0.03272 [+0.03046,+0.03502] | +0.04206 | −0.03583 [−0.04161,−0.03034] | +0.07263 | −0.03616 |
| 1.0 | +0.03230 [+0.02970,+0.03490] | +0.04009 | −0.05314 [−0.06012,−0.04645] | +0.07014 | −0.05362 |

Per-line appendix (Δ vs own β=0; β = 0.1/0.25/0.35/0.5/0.75/1.0):

| line | Δ wsign | Δ retrieval |
|---|---|---|
| k562_gwps | +0.0067/+0.0167/+0.0223/+0.0292/+0.0347/+0.0337 | +0.0013/+0.0007/−0.0013/−0.0060/−0.0164/−0.0288 |
| orion_hct116 | +0.0032/+0.0069/+0.0088/+0.0105/+0.0099/+0.0057 | +0.0014/+0.0021/+0.0006/−0.0033/−0.0133/−0.0264 |
| orion_hek293t | +0.0035/+0.0086/+0.0112/+0.0139/+0.0151/+0.0120 | −0.0012/−0.0057/−0.0118/−0.0237/−0.0454/−0.0668 |
| rpe1 | +0.0141/+0.0328/+0.0436/+0.0573/+0.0740/+0.0805 | −0.0078/−0.0210/−0.0309/−0.0477/−0.0752/−0.1008 |
| **nadig_jurkat** | +0.0067/+0.0161/+0.0216/+0.0283/+0.0340/+0.0338 | **+0.0039/+0.0058/+0.0051/+0.0017/−0.0097/−0.0241** |
| nadig_hepg2 | +0.0062/+0.0148/+0.0193/+0.0248/+0.0298/+0.0296 | −0.0019/−0.0061/−0.0100/−0.0162/−0.0289/−0.0429 |

Jurkat is the most forgiving line in the frame: the only one with
retrieval gains through β=0.5. The pooled curve falls monotonically
because the other five lines dominate it.

## Verdict against the PLAN.md kill criterion

- Gate passed (max |diff| 1.1e-16); curve read.
- Trade-off exists (sign up everywhere; retrieval down at high β) — not killed.
- Knee premise **confirmed with a Jurkat-specific twist**: efficiency
  Δwsign/|Δretr| on pooled5 falls monotonically (4.14 → 2.67 → 1.97 →
  1.40 → 0.91 → 0.61), maximum at the smallest β. **On Jurkat the ratio
  is undefined below β=0.75 because there is no loss** — retrieval
  improves at 0.1/0.25/0.35 (resolved) and is neutral at 0.5. Among
  loss-carrying Jurkat β, efficiency is maximal at **0.75 (3.51)** vs 1.0
  (1.40) — but that is the wrong optimum: it pays −0.0098 pds raw for
  the last 16% of sign.

**Recommended β = 0.35.** Largest β with *resolved* retrieval non-loss
(CI [+0.00038,+0.00985]): sign +0.0216 (64% of the β=1 gain), top-64
+0.0274, own cosine +0.0448, expected pds raw **+0.0051 (a gain, not a
cost)**. β=0.5 buys 84% of the sign at a point-neutral price (+0.0017,
CI straddles zero) — acceptable fallback if the PI values sign over
certainty. The PI's planned **β=0.25 is safe but leaves 16 points of
sign gain on the table** (+0.0161, 48%) at the same zero price; prefer
0.35. Do **not** ship β≥0.75 on this evidence (resolved pds cost
−0.0098 to −0.0243 raw).

## Limits

- Source-frame development evidence only (unit-RMS raw directions,
  512 genes, 1,024 non-panel targets). Not H1, not A/B/C, not emitted
  counts, not the six competition members.
- The pds translation borrows 346's H1 cross-study slope (1.009,
  Pearson 0.9991 over five arms). Applied here to a same-frame,
  same-study blend; the slope's transport to Jurkat/context-A is
  assumed, not measured. Treat prices as preregistration inputs with
  that caveat.
- Jurkat is atypical in this frame (only line with a retrieval gain
  region). The pooled curve — and every other single line — is
  strictly worse. If context A is not Jurkat, this recommendation
  does not transfer.
- Fixed-norm blending only; no amplitude claim (row norms pinned to
  the diagonal by construction). No DE-member (nmae/fid/reach/jac)
  translation — unknown by mandate.
- `core341.interval` resamples targets (and lines for pooled5);
  correlated genes/targets and study confounding remain as in 341.

## Exact commands to reproduce

```
cd ~/Code/virtual-cell-challenge-2026
memory_pressure | tail -1   # proceed only if free >= 50%
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/374-jurkat-knee-on-frame/run374.py
# writes work/curve.json + work/curve.csv; gate prints PASS/FAIL first
```

Top-64/own-cosine paired intervals were produced by the inline
heredoc in this session (same imports, `core341.interval` on Δ vs
β=0); values are quoted in the tables above. `run374.py` stores
wsign/retrieval intervals in `work/curve.json` permanently.

## Peak memory and wall time

- `run374.py` (gate + 7 β × 6 lines + 10k-draw intervals):
  wall **12.0 s**, peak RSS **362,741,760 bytes (0.34 GiB)**
  (`/usr/bin/time -l`), single process, 76–77% system memory free
  throughout. Follow-up top-64/cosine interval pass: seconds, <1 GiB.
- Wrote only `experiments/374-jurkat-knee-on-frame/`
  (`PLAN.md`, `run374.py`, `work/curve.json`, `work/curve.csv`,
  `RESULTS.md`). No git, no builder, no submission, no reserve read,
  no download. Budget: ≪ 8 GiB / 2 h.
