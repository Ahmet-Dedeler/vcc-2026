# 524dr BOARD RESULT: 0.2724054 (+0.0060441 vs 453dr; Σ +0.0054536). The October recipe composes as predicted: PASS

Entry `KTzIDR0Occ2jYzklkARE`, rank 23. 524dr is the October recipe as resolved on 10-05: def_day on 180's relabelled bundle with LODO 1,
lodo_hi_pct 96.647, sharpen_pt tercile:2.5,2.0,1.75 and de_mag dn_k250_zs. mde_sign was still PENDING, so OFF. Build: 664's re-run,
gate524 PASS; package cbeadd65.

**Delivery.**
- First entry VqBgraiBxlaK0WP5efOF (00:01Z 10-06) sat queued with no machine for 5.5 h and was cancelled (667).
- Resubmitted with 10-06's second slot at 18:32Z. The 18:53Z reboot cut its upload, which was resumed at 22:11Z (same entry).
- Scored and captured 23:04:05Z.

**Read** (`experiments/658-prereg-1007-builds/read658.py dr KTzIDR0Occ2jYzklkARE --recipe experiments/524-october-recipe-dr/work/recipe524.json --hook HoQY4xrzYer5y232oh85`, rc 0):

| | overall | bulk | DE | raw pds | raw mse | raw nmae |
|---|---:|---:|---:|---:|---:|---:|
| Σ of the switch reads (de_mag +0.00055, lodo +0.00223, lodo_hi_pct +0.00163, sharpen_pt +0.00104) | +0.0054536 | +0.0054876 | −0.0000341 | +0.000505 | −0.030202 | −0.001384 |
| **524dr − 453dr** | **+0.0060441** | +0.0053541 | +0.0006901 | −0.000609 | −0.031812 | −0.001599 |

- **SECONDARY (vs 453dr, 524's bins): d − Σ = +0.0005906, inside ±0.0020 → PASS.** No falsifier.
- The four October switches compose on the board the way their single reads predicted. The bulk lands at Σ (−0.00003); the DE part
  is +0.0007 above Σ, within the same-law re-roll band of the reach member (+0.0015 raw vs Σ's −0.0021).
- **R6 hand step:** PASS, so no HOLD; 10-07 s1 (CPH) goes ahead. The read landed at 23:4xZ, after chain1007's 20:00Z (re-armed 22:10Z
  after the reboot), and before CPH's 00:01Z submit.
- It is the base of 524dr2's primary read (branch N3: tercilec + mde_sign new vs 524dr).
