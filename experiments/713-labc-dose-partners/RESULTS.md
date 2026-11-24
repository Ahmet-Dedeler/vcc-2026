# 713 — LABC dose partners of 705 (θ) and 707 (β): no partner qualifies; the reference doses stand

Run 2026-10-16 11:53–11:54Z through gate602 (`work/price713.txt`, `work/price713.json`, `work/boot713.npz`; peak RSS 2.41 GB).
Registered before any number (REGISTRATION.json, 032a91ad; logged by the lab 2923077b). The driver gained only the paired
sd output (ed6a0df5, sha256 bfbcdb8b…) before the run. Written 2026-10-16T11:57Z by the clock.

**Instrument.** LABC (A/B/C's 300 targets, def rules, day-recipe operators; the lab's Arc-drawn surface), score708's exact
arithmetic (522 units, k 1, ret 1; percentile-cap incumbent, q 1.75). **KA0 PASS:** the incumbent's cos_agg and pds equal the
lab's stored `pct` arm on all 3 folds (|d| 0.0). **KA1 PASS:** prices recomputed from per-target arrays equal score708's (15
rows). Paired target bootstrap within folds, 4,000 draws, seed 713.

| config | K562 | HCT116 | HEK293T | mean (n+) | 95 % CI |
|---|---:|---:|---:|---|---|
| s_b0.125 vs incumbent | +0.00080 | +0.00151 | +0.00086 | +0.00106 (3/3) | [+0.00055, +0.00158] |
| s_b0.25 vs incumbent | +0.00123 | +0.00254 | +0.00099 | +0.00159 (3/3) | [+0.00061, +0.00261] |
| s_b0.5 vs incumbent | +0.00080 | +0.00272 | -0.00129 | +0.00074 (2/3) | [-0.00113, +0.00275] |
| r_t0.25 vs incumbent | +0.00041 | +0.00004 | -0.00029 | +0.00005 (2/3) | [-0.00024, +0.00036] |
| r_t0.5 vs incumbent | +0.00057 | -0.00009 | -0.00069 | -0.00007 (1/3) | [-0.00065, +0.00053] |

| paired (partner − reference) | K562 | HCT116 | HEK293T | mean (n+) | 95 % CI | sd | e(N) | best-of-N lower bound |
|---|---:|---:|---:|---|---|---|---|---|
| s_b0.125-s_b0.25 | -0.00044 | -0.00104 | -0.00012 | -0.00053 (0/3) | [-0.00105, -0.00005] | 0.00026 | 1.163 | -0.00135 |
| s_b0.5-s_b0.25 | -0.00044 | +0.00017 | -0.00228 | -0.00085 (1/3) | [-0.00183, +0.00015] | 0.00051 | 1.163 | -0.00243 |
| r_t0.5-r_t0.25 | +0.00016 | -0.00013 | -0.00041 | -0.00013 (1/3) | [-0.00042, +0.00018] | 0.00015 | 1.029 | -0.00057 |

**Decision (the registered rule).** No partner is > 0 against its reference on 3/3 folds. Every best-of-N lower bound is
negative (707 family N 5, e 1.163; 705 family N 4, e 1.029; the lab's counts). So neither family sends a 10-19 / 10-20 slot
candidate, and the reference doses (β 0.25, θ 0.25) stand.

**What it says (class: transfer; one Arc-drawn panel, so direction and size only).**
- **707's swap has an interior optimum near β 0.25** on an Arc-drawn panel. β 0.25 reads +0.00159 [+0.00061, +0.00261], 3/3
  (the lab's LABC read of 707 vs pct+P2g: +0.00115). Halving it costs −0.00053 (0/3 folds) and doubling it −0.00085 (1/3).
  706's October-like panels favoured MORE β. On Arc's panel more β is worse, which is the lab's fresh-panel over-read again.
- **705's rider is null on an Arc-drawn panel at both doses**: θ 0.25 +0.00005 [−0.00024, +0.00036], θ 0.5 −0.00007. 522's
  October-like panels preferred θ 1–2. That preference does not survive on Arc's panel, consistent with the lab's LABC read
  (≈ 0).
- Both families' October use goes through 708's day-of folds on D/E/F's own panel (the swap) or 705's board read (the
  rider). Nothing here changes a decision.

## Part 2 — the H1 dose on LABC (the coordinator's fallback slot (2)): both points ≤ 0, so slot (3) takes the slot (recorded 2026-10-16T13:11Z)
- Run 13:10:11–13:10:24Z through gate602 (`work/price713h.txt`, `work/price713h.json`; peak RSS 2.62 GB).
- Registered before any number: REGISTRATION-H1.json (ee791df7), with the LABCd shas added before reading (88d988fe).
- Dumps: the lab's LABCd rebuild. Science's own check: every original LABC array is byte-identical, and the only new arrays
  are fin_h2.5c1 and fin_h3c1.
- **KA0 PASS:** the incumbent equals the lab's LABC pct arm (|d| 0.0) on 3 folds. **KA1 PASS.**
- All prices in 522 units, k 1, ret 1, against h2c1 (H1 ×2, CD4 ×1) under the incumbent percentile-cap key. Best-of-N
  lower bound uses the lab's N 4 (e 1.029).

| design | K562 | HCT116 | HEK293T | mean (n+) | 95 % CI | best-of-N lb | HARM on LABC |
|---|---:|---:|---:|---|---|---|---|
| A:h2.5c1 (H1 ×2.5, CD4 ×1) | +0.00009 | +0.00036 | −0.00081 | −0.00012 (2/3) | [−0.00048, +0.00028] | −0.00068 | no |
| A:h3c1 (H1 ×3, CD4 ×1) | −0.00022 | +0.00034 | −0.00203 | −0.00064 (1/3) | [−0.00133, +0.00014] | −0.00172 | **yes** |
| step ×2.5 → ×3 (paired) | −0.00031 | −0.00003 | −0.00122 | −0.00052 (0/3) | [−0.00085, −0.00014] | −0.00103 | — |

- **Decision (the registered rule): both points ≤ 0 on the mean, so science reports so and slot (3), R-ABC2's seed-only
  re-roll, takes the fallback slot.** A:h3c1 is also HARM on LABC (mean ≤ −0.0005, 2/3 folds negative).
- **What it says (transfer; one Arc-drawn panel):**
  - At CD4 ×1 and under the production key, raising H1 above ×2 does not pay on an Arc-drawn panel. ×2.5 is null and ×3 hurts.
  - This sits beside 438's levelled K3 table (062 operators, CD4 ×1.5), which recommended H1 ×2.5 / CD4 ×1.5 (+0.0022). That
    gap is what 708's AMENDMENT 1 exists for: a Part A H1 move must be confirmed under the production key.
  - It agrees with the lab's prospective NULL for 710 (H1 ×2.5 on the board, at the A/B/C CD4 dose ×1.5).
