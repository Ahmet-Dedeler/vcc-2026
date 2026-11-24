# 551 — C10b, R6's second and last access: does the final October bulk configuration's last step hold on the sealed panels?

Written 2026-09-30 by the session "Explore contribution areas", before any C10b access and before the 10-02 … 10-04
reads that fix the configuration. Coordinator agreement: its message of 2026-09-30 ~21:10Z (R6's C10b is this session's;
the ledger section for C10b is this session's to write). Timing amendment: `docs/CONFIRMATION-LEDGER.md`, "Claim C10b —
timing amendment" (da543a7c).

## 1. Question

C10a (534 Part B, 2026-09-29) priced the October bulk recipe's milestones M0 → M6 on R6 once: every step TRANSFERS, and
the recipe keeps 0.70 [0.62, 0.78] of its A/B/C gain. After C10a, two October bulk switches remain open:

- `lodo_hi_pct` (10-03, 495): the LODO key's cap at 3 (M5) or at the 96.647th percentile (M6, priced by C10a);
- `sharpen_pt` (10-02, 479's pair): per-target sharpening on top of the base, S_P2 or P2g, or none.

The October path has no union-source knob, so 471c (10-03) cannot enter the October bulk without a new pass-through.
C10b asks, once: **does the one bulk step the final configuration adds after C10a (per-target sharpening) hold on R6?**

## 2. The configuration F, and how it is fixed

- F is one of six: `M5`, `M6`, `M5+S_P2`, `M5+P2g`, `M6+S_P2`, `M6+P2g` (`designs551.CONFIGS`).
- The coordinator names F after the 10-04 read (~04:30Z) by message. That message is copied verbatim, with its time,
  into `work/c10b-config.json`, whose sha256 is appended to the ledger's C10b section BEFORE the first R6 byte is read.
- F is chosen by the A/B/C board reads alone (479's, 495's and 498⟨P⟩ / 497dr's preregistered rules). No R6 number enters
  that choice.
- **If F is `M5` or `M6`** (no per-target step), F is a C10a design and its R6 value is C10a's. C10b then does not access R6;
  the ledger records F and cites C10a. R6 stays sealed.
- **If the coordinator names a bulk switch outside these six** (for example the strength rider of 522, or the union
  source), C10b prices only the part of F these designs express and records what is not covered. The uncovered part goes
  to R7.

## 3. Designs (lib438 names; `designs551.py`)

| id | design | what it is |
|---|---|---|
| M5 | `h2c1\|V1_5_3+sh1.75` | the LODO key, clip (0.5, 3), q 1.75 (C10a) |
| M6 | `h2c1\|V1pct+sh1.75` | the LODO key, clip (0.5, the P 96.64697522769474-th percentile of carried key/median) (C10a) |
| b+S_P2 | `h2c1\|<key>+ptS_P2` | q (2.5, 2.0, 1.75) by tercile of response strength (the unkeyed h2c1 row's non-panel norm) |
| b+P2g | `h2c1\|<key>+ptP2g` | q (2.5, 2.0, 1.75) by the carried-only tercile of the base key's weights; non-carried MID; floor LOW |
| control | `CTRL:h2c1\|<key>+pt<form>perm` | F's own q vector permuted over targets (`np.random.default_rng(493)`), at its own level factors |

The per-target forms are 493h's `P2g_carried` / `S_P2` and 546o's `p96.6 + P2g` / `p96.6 + S_P2`: the forms priced on the
October-like panels. On each PLOLO fold the terciles are cut on the fold's targets; for the level factors, on the full
panel (as 534 did for the percentile cap).

## 4. Known answers

Before any R6 access (on `efp`, and `abc` for the supplementary column), abort on any failure, |d| ≤ 1e-12:
- KA-FOLD, KA-438, KA-493, KA-546 (M5), KA-522 (M6): 534's rows, unchanged (`comp551.py efp … --ka`).
- KA-534: every 534 design and its levelled twin == 534's own efp `comp534-<fold>.json`; the full panel's energies, rho,
  level factors and keys == 534's efp `full534.json` / `level534.json` / `keys534-full.npz`.
- KA-PT0: the per-target path at a constant q 1.75 == the scalar path, bit for bit (M5 and M6 energies).
- KA-546o: P 96.647 (543's rounded P) at q 1.75, + P2g, + S_P2 == 546o's `p96.6`, `p96.6 + P2g`, `p96.6 + S_P2` on efp
  (design and realised 16; K_mean).
- KA-493h: M5 + P2g / + S_P2 == 493h's `P2g_carried` / `S_P2` on efp (design and realised 16).

Inside the access, before any new design's number is read (abort on failure):
- every R6 operator is re-fitted with 534's exact command (`state534.json`) and `--expect-sha256` = C10a's sha256;
- every rebuilt fold cache == C10a's `fold-<held>-held.npz`, bit for bit;
- KA-REPRO: every 534 design and its levelled twin == C10a's stored `comp534-<fold>.json` on every fold of every panel,
  and C10a's eight pooled step values == `results534B.json`, |d| ≤ 1e-12.

## 5. Units and statistics

- Primary: 438's LEVELLED route at k = 1, Δoverall in board units (`price551` == `price534` == `price438`).
- Per panel: the mean over its three cancer folds. Pooled: the mean over R6's 12 folds.
- Paired target bootstrap, 4,000 draws, one generator seeded **551**, consumed panel by panel (S1 … S4), fold by fold
  (k562, hct116, hek293t), draw by draw.
- Also reported: the calibrated board-unit projection (shape: k 0.65 for X, 0.42 for pds), the retention vs A/B/C
  (C10a's record values for its six steps; F's step on A/B/C's own panel with this instrument), and the levelling
  sensitivity (534's energy-only extrapolation to ρ*).

## 6. What is priced and printed

Only F's own step (b → b + form), its control (b → the permuted control), and the cumulative M0 → F (C10a's six steps +
F's step, the same draws). The other per-target candidates are never computed on R6 (`comp551.py --config F`), so no one
can select among them on R6 scores.

## 7. Reading rules (frozen now)

1. **Verdict** on F's step, R6 pooled: TRANSFERS if Δ > 0 and the 95 % lower bound > 0; FRAGILE if Δ > 0 and the lower
   bound ≤ 0; FAILS if Δ ≤ 0.
2. **Targeting:** claim − control, paired, with its interval. If its interval covers 0, the report says the step's R6
   value does not depend on WHICH targets get which q (only on the q mix).
3. **Consequence:** FAILS → `sharpen_pt` goes to the coordinator as a candidate to switch off for October; FRAGILE → flagged;
   TRANSFERS → nothing. This experiment changes no production file, recipe value or slot.
4. The outcome is appended to the ledger's C10b section whatever it is.

## 8. Prior expectation (written before the access; from 493h / 546o on real, efp, efr)

- On the percentile-cap base, per-target sharpening adds +0.24 … +2.56 % relative cosine on all 9 October-like folds
  (546o), with literal Δpds +0.004 … +0.011. In board units at k 1 on R6: expect > 0, TRANSFERS for S_P2 more likely than
  for P2g (S_P2 is the more panel-robust form, 493g/493h/546o).
- The permuted control near 0 (479's PLOLO note: "+0.0000, direction −0.5 %").
- A surprise is F's step ≤ 0.

## 9. Memory, disk and time

- Heavy processes: the 16 operator re-fits (~12–13 GiB RSS, ~5.3 min each; `--max-rss-gib 14`) and the pricing runs
  (≤ 8 GiB). Strictly serial, each under `/usr/bin/time -l`.
- Gate before each heavy process (amended 2026-09-30 ~22:40Z, before any access, with the coordinator's approval):
  `memory_pressure` ≥ 55 % free, `vm_stat` free+inactive ≥ 22 GiB, swap free ≥ 400 MiB, no builder
  (`def_pipeline.other_builders()` empty; no bash process running `run_arm_generic.sh`, `def_day.sh`, `run_ka_*`, chain or
  prep builds; no `vcc prep`), data volume ≥ 25 GiB free. (The original 27 GiB floor cannot be met while desktop apps hold
  ~12 GB in the compressor, which free+inactive does not count; the pressure floor is the real guard.) The four operators of a panel are deleted after that panel is priced (sha256 and
  bytes recorded), so the footprint stays ≤ ~6 GiB.
- Window (the coordinator's, message 2026-09-30 ~21:50Z): heavy steps start only in 10-04 05:00–12:00Z, after its
  configuration message (~04:30Z); if its read slips, the window slips with it but ends by 10-04 16:00Z (`gate551`,
  PHASE `access`). 10-04 12:00Z – 10-05 18:00Z are the coordinator's build windows. Not after 10-06 12:00Z without a new
  ledger amendment.
