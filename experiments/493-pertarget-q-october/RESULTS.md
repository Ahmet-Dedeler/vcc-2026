# 493 — Per-target sharpening on October-like panels: the direction gain holds everywhere. The knob needs 492's floor guard, and the strength-based assignment is the most robust

2026-09-26 23:30Z, third session. Analysis only.
- Instrument: 438's panel-agnostic PLOLO (`lib438`) with the October stack h2c1 (H1 ×2, CD4 ×1), the K3 key and
  ballast realisation at KMAX 16, on three cancer folds of each panel.
- Script `comp493.py` (runner `run493.sh`, gated on vm_stat free+inactive ≥ 33 GiB; each run 15 s, peak RSS 5.6 GB).
- v1 (no guard) is kept in `work/v1/`; the final run adds the guarded arms.

## Known answer
The per-target path at q 1.5 for every row reproduces 438's stored `h2c1|K3+sh1.5` row (design and realised 16) with
|d| 0 on all 12 panel/fold pairs.

## Result (direction members vs the q 1.75 base on the same panel; mean over the 3 cancer folds, folds positive in brackets)
| panel | K3 floor share (3 folds) | P2 rel cos / Δpds_bd | **P2g** (floor guard) | L30 → L30g | **S_P2** (strength tercile) | CTRL (P2 permuted) |
|---|---|---|---|---|---|---|
| abc | 0.26–0.28 | +1.12 % (3/3) / +0.0087 | +1.12 % / +0.0087 (identical) | +0.51 % → +0.51 % | +1.67 % (3/3) / +0.0077 | −0.44 % (0/3) / +0.0078 |
| efp | 0.18–0.35 | +0.28 % (3/3) / +0.0034 | **+0.35 % (3/3) / +0.0050** | +0.07 % → +0.21 % | **+1.85 % (3/3) / +0.0063** | −1.85 % (0/3) |
| efr | **0.35–0.41** | +0.62 % (3/3) / +0.0022 | **+0.78 % (3/3) / +0.0019** | **+0.00 % (no-op) → +0.36 %** | **+2.01 % (3/3) / +0.0052** | −0.51 % (1/3) |
| real | 0.32–0.39 | +0.54 % (3/3) / +0.0039 | **+0.62 % (3/3) / +0.0068** | +0.03 % → +0.23 % | **+1.06 % (3/3) / +0.0073** | −2.48 % (0/3) |

- **The mechanism ports to October-like panels.** P2's direction gain is positive on all 9 October folds, and its
  board-emulated pds is positive on 9/9.
- The permuted control loses direction on every panel. It still gains pds, because it too sharpens a third of the rows
  harder.
- **The floor guard is necessary, not cosmetic.**
  - On October-like panels 18–41 % of targets sit at the K3 floor.
  - When more than a third do, `np.quantile`/`digitize` puts the whole floor in MID and leaves LOW empty. L30 is then an
    exact no-op on `efr`, and P2 loses part of its gain.
  - With `low = (w < q1) | (w <= 0.3)`, P2g ≥ P2 on every October panel and equals P2 on A/B/C.
- **The strength-based assignment (S_P2: the unkeyed row norm's tercile) is the strongest and the most panel-robust:**
  +1.1 … +2.0 % direction and +0.005 … +0.007 pds on every October panel.
  - Strength is continuous, so it has no floor ties, and it does not depend on the key's clip.
  - It moves more A/B/C energy (×0.93 in 492d), which the mandatory day-of re-level absorbs.

## Recommendation (the coordinator's decision)
- 479's knob should implement the floor guard: every target at the floor goes to LOW. The A/B/C build is unchanged by it.
- For October, the assignment variable is worth a board comparison.
  - Key tercile with the guard (P2g, the preregistered family) against response-strength tercile (S_P2).
  - Both are computable on the day. S_P2 is the better October default on this evidence (3 panels × 3 folds, 9/9
    positive, largest gains).
  - It should be preregistered on its own after P2 reads.

## Addendum 493c (00:15Z, EXPLORATORY): floor 0.5 stacked with the per-target q on the October panels (`comp493.py <panel> <fold> --stack`)
The K3 key is re-clipped at 0.5 and the tercile stays the K3 one (guarded). Known answer: h2c1|K3+sh1.5 == 438, |d| 0 on
every run. 8 of 9 panel/folds ran; real/orion_hek293t was skipped when the memory gate closed at the coordinator's
rehearsal. Direction members vs the q 1.75 base (K3):

| panel | F5 alone (floor 0.5, q 1.75) | P2g | **F5 + P2g** | S_P2 | **F5 + S_P2** |
|---|---|---|---|---|---|
| efp | +0.30 % (3/3) | +0.35 % (3/3) | **+0.74 % (3/3)**, pds_bd +0.0065 | +1.85 % (3/3) | **+2.21 % (3/3)**, pds_bd +0.0073 |
| efr | +1.13 % (3/3) | +0.78 % (3/3) | **+2.06 % (3/3)**, +0.0052 | +2.01 % (3/3) | **+3.14 % (3/3)**, +0.0083 |
| real (2 folds) | +0.20 % (2/2) | +0.82 % (2/2) | **+1.15 % (2/2)**, +0.0101 | +0.81 % (2/2) | **+1.03 % (2/2)**, +0.0085 |

- The floor and the per-target sharpening stack on October-like panels, as on A/B/C (492c).
- As before, the floor alone moves direction through magnitude reallocation but not literal pds (+0.0000 … +0.0002).
- Render energy on frame A: ×0.99–1.00 for the P2g stacks and ×0.91–0.95 for the S_P2 stacks. The mandatory re-level
  refunds it.
- For the October bulk default, after the P2g / S_P2 board pair reads: the winning assignment, then the floor as its own
  step.

## Addendum 493g / 493h (04:35–04:50Z): the same arms under the OCTOBER key (the LODO key at clip 0.5, 3; `DEF_DAY_LODO=1`)
`comp493g.py` / `comp493h.py` (lib438, October stack h2c1, the fold's LODO-V1 at (0.5, 3) computed as comp438 does).
Known answers: the per-target path at q 1.5 == 438's stored `h2c1|V1_5_3+sh1.5`, |d| 0 on 9/9 panel/folds. 493h's P2g
== 493g's, |d| 0 on 9/9. Peak RSS 5.65 GiB, ~20 s per run, after 04:30Z at ≥ 32 GiB free+inactive.

| panel | LODO floor share | non-carried | P2g (guard) | P2 (no guard) | **S_P2** | L30g | CTRL (P2g permuted) |
|---|---|---|---|---|---|---|---|
| efp | 0.37 / 0.22 / 0.20 | 6 / 153 / 156 | +0.48 % (2/3) | +0.41 % (2/3) | **+1.32 % (3/3)** | +0.18 % (2/3) | −1.84 % (0/3) |
| efr | 0.41 / 0.45 / 0.39 | 0 | +1.18 % (3/3) | +0.74 % (3/3) | **+2.36 % (3/3)** | +1.00 % (3/3) | −1.01 % (1/3) |
| real | 0.38 / 0.41 / 0.36 | 0 | +0.91 % (3/3) | +0.70 % (3/3) | +0.89 % (3/3) | +0.47 % (3/3) | −5.39 % (0/3) |

(relative cos vs the q 1.75 base, mean of the three cancer folds; literal Δpds is positive for every arm on every panel,
+0.004 … +0.011.)

- **The floor guard is necessary under the LODO key.** The LODO floor holds 36–45 % of targets on efr/real, and
  without the guard LOW is empty. P2g ≥ P2 on every panel.
- **S_P2 (strength tercile) is the most robust October assignment:** 9/9 folds positive, the largest gains on efp and
  efr, and tied on real. It has no ties to handle.
- The permuted control loses direction on every panel, so the targeting is real under the October key too.
- **A second tie (493h):** targets carried by fewer than 2 donors sit at w = 1 exactly. On efp (140 of ~300 carried)
  that block straddles a cut point, and the key tercile's MID is empty on 2 folds.
  - Carried-only cut points with non-carried → MID repair it: efp +0.84 % (3/3), vs +0.48 % (2/3) as is, +0.63 % with
    non-carried → MID on the plain cut points, and +0.69 % (2/3, most ballast) with non-carried → LOW.
  - The rule and the knob guidance are under 492's knob definition.
