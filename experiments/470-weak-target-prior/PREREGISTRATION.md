# 470 — which block prior do design-weak targets want? (preregistration, 2026-09-26 17:13Z, before any 470 number)

H1 only, production emitter, 467's harness and q33 strength gate, 3 seeds. The known answer is `wc4` seed 0, which must
equal 467's gated33 seed 0 bit for bit. Deltas are raw vs 435's CLEAN4 run at the same seed; the reference row is
gated33: fid +0.0161, nmae −0.0102, reach +0.0015, overall-equiv +0.0110.

## Arms
- `wdown`: Arc on the strong two thirds, every block gene DOWN on the weak third.
- `warcw` / `warcs`: the weak third gets Arc's programme from the bank's weaker / stronger half of perturbations.
  - Response size = count of |delta| > 0.2 over the 4,000 most expressed genes.
  - The median is 17. The halves are 108 / 107 of the 215 non-h1-test perturbations.
- `sarcs`: the strong two thirds get Arc's strong-half programme; the weak third keeps CLEAN4.
- Priors on the full 18,080-gene axis:
  - down fraction: CLEAN4 0.837, arc 0.549, arcw 0.521, arcs 0.556;
  - agreement with CLEAN4: down 0.837, arc 0.595, arcw 0.556, arcs 0.602.

## Predictions (vs gated33)
- `wdown`: reach within ±0.005 of gated33 and fid ≥ gated33. CLEAN4 is mostly down, and 361 found all-down buys fid on
  H1. Overall-equiv +0.000 … +0.002.
- `warcw`, `warcs`: reach 0.02–0.04 below gated33, like Arc's full programme on weak targets (−0.0388 ungated). A
  per-gene programme from this line misleads its weak targets, whichever half it is estimated from.
- `sarcs`: fid within ±0.002 of gated33. The strong half dominates the programme's non-noise content.

## Rule
- A weak-third prior replaces CLEAN4 in the gate only if all three hold:
  - overall-equiv ≥ gated33 + 0.001 (5 × the paired seed sd);
  - reach ≥ gated33 − 0.003;
  - pds and mse deltas stay exactly 0.
- Even then it is not a board arm before 458g's read on 09-30. Its board test would be 458g vs 458g + the new weak
  prior, on the same day.
- `sarcs` replaces `arc` on the strong side under the same rule.

## Side analysis 18:39Z 09-26 (before any 470 number; `downfrac470.py` → `work/downfrac470.json`)
Is a weak CRISPRi response down-dominated in other lines too?
- Proxy DE set per bank perturbation: |delta| > 0.3 on the bank's 2,000 most expressed genes, target excluded.
- Mean down fraction per perturbation:

  | bank | weak (1–20 genes) | mid (21–100) | strong (> 100) |
  |---|---:|---:|---:|
  | K562 GWPS | **0.652** (5,192) | 0.548 | 0.519 |
  | RPE1 | **0.633** (508) | 0.651 | 0.730 |
  | Jurkat (Nadig) | **0.666** (793) | 0.618 | 0.628 |
  | HCT116 (Orion) | **0.580** (8,728) | 0.584 | 0.607 |
  | HEK293T (Orion) | **0.624** (8,102) | 0.742 | 0.695 |
  | H1 (all cells) | **0.605** (147) | 0.555 | 0.666 |

- Weak responses lean down in every line (0.58–0.67), so a down-leaning prior for design-weak targets rests on a
  cross-line regularity. It is not weak-specific: strong responses also lean down in several lines (K562 is the
  exception at 0.52).
- This does not change 470's predictions or rule. It says `wdown`, if it wins on H1, has a reason to transfer.

## Addendum 20:00Z 09-26 (before any 470b number): 470b, the strong side by essentiality
The 2026 panel is essential-free, and Arc's bank mixes essential and non-essential perturbations: 60 and 155 of the 215
non-h1-test ones.
- `sarcne` gives the strong two thirds Arc's programme over the 155 non-essential perturbations. `sarcess` gives them
  the programme over the 60 essential ones. The weak third keeps CLEAN4. `run470b.py` is a copy of `run470.py` with
  the two priors added; `chain470b.sh` runs after chain471; 3 seeds each.
- Down fraction on the full axis: arcne 0.523, arcess 0.576 (arc 0.549).
- **Prediction** (vs gated33, on the 67 non-essential targets):
  - `sarcne` within ±0.002 overall-equiv;
  - `sarcess` fid lower by 0.002–0.006, since essential biology fits an essential-free panel's strong targets less
    well. On the 18 essential targets it gains.
- **Rule:** `sarcne` becomes the candidate programme for an essential-free panel only if all hold on the non-essential
  subset:
  - overall-equiv ≥ gated33 + 0.001;
  - reach ≥ gated33 − 0.003;
  - pds and mse exactly 0.

  A board test comes after 09-30's read and uses the same exact-split design.
