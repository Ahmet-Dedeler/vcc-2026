# 534 Part C — do the three small base steps hold off the board? RESULTS

**Normrel −1.5 holds on the October-like panels. Its cap form and the CD4 partial rows fail there:** their pooled value
is at or below zero. Both still read positive on A/B/C, in this instrument and on the board. By PREREGISTRATION-C §4.4
they are A/B/C-specific, and both go to the coordinator as candidates to switch off for October.

Worker for 534 Part C (the resumed Part B session, headless Opus), 2026-09-29, 17:40Z to 18:35Z. It ran
`PREREGISTRATION-C.md` (sha256 `274a5512…`, unchanged) under `MANDATE-C.md`.
- Panels priced: `abc`, `real`, `efp`, `efr`, each on its three cancer folds with its existing held-out operators, read
  only. **Reserve R6 was not touched.**
- No operator fits; no `vcc`, no `cell-eval2`, no git, no builder edit.
- Nothing written outside `experiments/534-generalization-audit/`: the write guard was active in every process.
- Two owner deviations (§6): the gate floor, and an exception for the armed build535.sh waiter.

## 1. Decision table (PREREGISTRATION-C §4; cite this)

Δ_s = B0 − B_s: 438's levelled route, k = 1, overall board units. B0 is the October bulk configuration (h2c1, the LODO
key clip (0.5, 3), q 1.75). B_s is the same configuration with one base step removed. Pooled means the nine folds of
real, efp and efr, with a paired target bootstrap (4,000 draws, seed 534).

| step | Δ pooled real/efp/efr (k = 1, board units) | 95 % CI (paired target bootstrap, 4,000 draws, seed 534) | folds positive (of 9) | verdict | Δ on `abc` [95 % CI] (folds +) | the board's A/B/C read |
|---|---:|---|---|---|---|---|
| **normrel −1.5 with its cap** (B0 − B2: per-target donor weights vs equal weights) | **+0.00112** | [+0.00041, +0.00182] | 8/9 | **HOLDS** | +0.00116 [+0.00046, +0.00193] (3/3) | +0.0026: 312 +0.0016 (off → plain −1.5), then 319 +0.0009 (plain → cap) |
| **the cap form** (B0 − B1: cap vs plain at −1.5) | **−0.00011** | [−0.00049, +0.00028] | 3/9 | **FAILS** | +0.00028 [−0.00017, +0.00074] (3/3) | +0.0009 (319) |
| **the CD4 partial rows at g 0.5** (B0 − B3: vs none) | **−0.00026** | [−0.00083, +0.00026] | 4/9 | **FAILS** | +0.00094 [+0.00017, +0.00173] (3/3) | +0.0015: 319p +0.0012 (g 0.34), then 319p2 +0.0003 (→ 0.5) |

- **§4.2, the verdicts, applied exactly.** The cap and the partial rows have pooled Δ ≤ 0, so they FAIL. Their intervals
  include 0: P(Δ ≤ 0) is 0.71 and 0.82 in the bootstrap. The rule does not separate "≈ 0" from "clearly negative", and no
  rounding was applied.
- **§4.3.** The `abc` row agrees in sign with the board's positive read for all three steps. So the instrument does see
  these steps on A/B/C's own panel. For the cap and the partial rows, it does not see them on the October-like panels.
- **§4.4, the consequence.**
  - **The cap and the CD4 partial rows are A/B/C-specific.** Both FAIL on the October-like pool, and both have positive
    board reads. They go to the coordinator as candidates to switch off for October. That is a `def_day` change: the
    coordinator's call, needing a pass-through and a rehearsal.
  - Normrel HOLDS: nothing to do.
  - Nothing in production changed here.
- **§5 prior, as expected.**
  - normrel was expected to HOLD or be FRAGILE: it HOLDS.
  - The cap was expected FRAGILE or FAILS: it FAILS.
  - The partial rows were expected FRAGILE, small, either sign: they FAIL (small, negative).

**What each step means for October:**
- **normrel −1.5 (HOLDS).** Per-target donor weighting is worth +0.0011 [+0.0004, +0.0018] off A/B/C, positive on 8/9
  folds and on every panel's mean. On A/B/C's own panel it is +0.0012 in the same instrument. It is transfer, and it stays
  on for October.
- **The cap (FAILS).** The cap on top of the weighting is worth −0.0001 [−0.0005, +0.0003] off A/B/C (3/9 folds), against
  +0.0003 on A/B/C's panel and +0.0009 on the board.
  - It was chosen from the board's own dose curve and never helped off the board.
  - Switching it off for October (plain −1.5) is expected to be neutral, about +0.0001, inside noise.
  - Keeping it costs about as little.
- **The CD4 partial rows (FAILS).** They are worth −0.0003 [−0.0008, +0.0003] off A/B/C (4/9 folds), against +0.0009
  [+0.0002, +0.0017] on A/B/C's own panel and +0.0015 on the board.
  - This is the same A/B/C-only pattern as 438's CD4 ×1.5. The CD4 residual rows point at the held-out truth worse on
    October-like panels than on A/B/C (438 §4.5).
  - Switching them off for October is expected to be neutral to slightly positive (about +0.0003), inside noise.

## 2. Per-panel values (levelled, k = 1; mean of the three cancer folds, folds positive in brackets)

| step | real | efp | efr | pooled real/efp/efr | abc | pooled, unlevelled route |
|---|---:|---:|---:|---:|---:|---:|
| normrel −1.5 with its cap | +0.00038 (3/3) | +0.00154 (3/3) | +0.00145 (2/3) | +0.00112 (8/9) | +0.00116 (3/3) | +0.00112 |
| the cap form | −0.00007 (1/3) | +0.00023 (2/3) | −0.00049 (0/3) | −0.00011 (3/9) | +0.00028 (3/3) | −0.00005 |
| the CD4 partial rows | −0.00033 (1/3) | +0.00013 (2/3) | −0.00059 (1/3) | −0.00026 (4/9) | +0.00094 (3/3) | −0.00021 |

Per fold (k562 / hct116 / hek293t), from `results534C.json`:

| step | real | efp | efr | abc |
|---|---|---|---|---|
| normrel | +0.00025 / +0.00086 / +0.00003 | +0.00142 / +0.00208 / +0.00112 | −0.00005 / +0.00256 / +0.00183 | +0.00066 / +0.00029 / +0.00251 |
| cap | +0.00016 / −0.00028 / −0.00009 | −0.00053 / +0.00085 / +0.00036 | −0.00022 / −0.00069 / −0.00056 | +0.00006 / +0.00059 / +0.00019 |
| partial | −0.00059 / +0.00152 / −0.00190 | −0.00051 / +0.00061 / +0.00031 | +0.00038 / −0.00061 / −0.00153 | +0.00080 / +0.00012 / +0.00191 |

Decomposition (panel means: relative aggregate cosine d_rel, board-emulated Δpds, raw Δmse):
- **normrel** gains through direction on every panel: d_rel +0.65 % (abc) to +2.6 % (efp), and Δmse −0.0022 to −0.0063.
- **The cap** moves direction by ≤ 0.7 %. On `efr` it costs board-emulated pds (−0.0016).
- **The partial rows** lose direction on `real` (d_rel −0.42 %) and pds on `efr` (−0.0017). On `abc` both are positive.

## 3. The instrument

- **The copy.** `lib534c.py` is Part B's `lib534.py` (lib438's redirected copy) with its caches under `work/partC/` and
  three base toggles read at call time:
  - `NORMREL_MODE` and `NORMREL_EXP` enter `lib418.pre_normrel` through `pre_base()` / `lodo_bases_base()`.
  - `PARTIAL_G` enters `compose()`, now passed to `compose_arm` on the abc-rules path too.
- **Designs.**
  - **B0 (production):** cap −1.5, partial g 0.5.
  - **B1:** plain −1.5.
  - **B2:** exponent 0, i.e. equal weights, bit-identical to the uniform mean.
  - **B3:** partial g 0.
  - Every B is h2c1 with its own LODO key clip (0.5, 3) and q 1.75.
- **What follows the toggle, and what does not.**
  - Each design re-derives its own donor mean, whitened rows, compositions, LODO key and carriers from its own base, as
    418 priced its normrel and partial arms.
  - The fold's calibration constants stay the production base's: the board-emulated pds beta, the pricer's X anchors
    h1c1|K3 and h1c1|none, and ret_f.
- **Energies and levels.**
  - `abc`, `real` and `efp` use the full-panel route: each design's energy on the panel's own 3-donor operator
    (`full534c.py`).
  - `efr` has no 3-donor operator, so it uses 438's fold route: fold-averaged energies and fold-route factors
    (`foldlevel534c.py`, two passes).
  - The pricer is price438.Panel.price, verbatim (`price534c.py`).

## 4. Known answers (PREREGISTRATION-C §3, all before any B1–B3 number; every script aborts on failure)

On every fold, the production checks run and pass before any toggled base is built. Each toggled base passes its own
check before its B row is scored.

| known answer | what | result |
|---|---|---|
| KA-FOLD | my fold caches, built from scratch, vs 438's `fold-<held>-held.npz` | bit-identical, 12/12 folds |
| KA-438 (toggles at production == the unpatched copy) | h1c1\|none, h1c1\|K3 and h2c1\|V1_5_3+sh1.5 (the nearest design 438 stored; §3's fallback), design, realised 16 and levelled twins, vs 438's stored rows | \|d\| = 0 on every fold of all four panels (efr levelled only what 438 levelled) |
| KA-546 (the exact B0 design) | B0 vs 546o's `base q1.75` (unpatched lib438), design and realised 16 | \|d\| = 0, 9/9 folds (real, efp, efr) |
| KA-534B | B0 and B0:L vs Part B's h2c1\|V1_5_3+sh1.75 rows; full-panel energies and level factor | \|d\| = 0 (efp, abc: 6 folds + 2 full panels) |
| KA-UNIFORM (§3: B2 = the uniform-weight composition) | B2's base at exponent 0 vs lib418.pre_normrel's uniform mean (the builder with `VCC_NORMREL_EXP` unset) | bit-identical on 12/12 folds and 3/3 full panels |
| KA-418 (stored toggled compositions, abc) | the h1c1\|K3 design rows of B1 / B2 / B3 vs 418's `f2:plain-1.5\|K3` / `e:exp0\|K3` / `d:g0\|K3` | \|d\| = 0, 3/3 folds each |
| KA-COMP | production PP constants vs the builder's meta | abc 8.87e-06, real 2.57e-05, efp 8.93e-06 (tolerance 5e-5; 438's values) |
| KA-KEY | full-panel 410 key vs the builder's artifact | identical (abc `d7f76a22`, real `4f1f0209`, efp `c3e429dd`) |
| KA-E438 | full-panel energies, anchored AB / AB_level, ρ and level factors of the production designs vs 438's `full438.json` / `level438.json` | relative \|d\| = 0 (abc, real, efp) |
| efr fold route | claims vs 438's `efr_s1`; fold-route factors of h1c1\|K3 and h2c1\|V1_5_3+sh1.5 vs 438's efr `level438.json` | \|d\| = 0 |
| KA-RESUME (efr) | the levelled second pass vs the first pass's rows | 7 rows × 3 folds identical |
| KA-BOOT | the bootstrap's identity draw vs the json statistics | < 1e-12, every fold |

Per fold: 14 (real), 18 (efp), 19 (abc, efr) checks. The worst |d| is 0 everywhere, apart from KA-COMP's stated
tolerance.

## 5. Secondary (not preregistered): the same steps with the production LODO key kept

An October knob that switched a step off in the builder alone would still load the key artifact that `def_lodo_key.py`
computes on the production base: it hard-codes normrel cap −1.5, and its composition uses PARTIAL_G 0.5. So I also priced
each toggled base with B0's LODO key (`--pk` runs).
- Every primary row reproduces bit for bit in those runs, and so does every primary step value.
- These results are in `results534Cpk.json`. The same rule, applied descriptively:

| step (production key kept) | Δ pooled | 95 % CI | folds + | the rule would read | Δ on abc | real / efp / efr |
|---|---:|---|---|---|---|---|
| normrel | +0.00098 | [+0.00028, +0.00168] | 7/9 | HOLDS | +0.00114 | +0.00011 / +0.00129 / +0.00153 |
| the cap | −0.00024 | [−0.00062, +0.00015] | 3/9 | FAILS | +0.00022 | −0.00028 / +0.00010 / −0.00055 |
| the partial rows | −0.00021 | [−0.00075, +0.00029] | 4/9 | FAILS | +0.00097 | −0.00006 / −0.00009 / −0.00048 |

**The verdicts do not depend on whether the key follows the base.** A builder-only knob and a knob passed through to
`def_lodo_key.py` would read the same.

## 6. Deviations (the owner's), memory and timing

1. **Gate floor, the owner's amendment of 2026-09-29T15:43Z.**
   - The floor is `vm_stat` free+inactive ≥ 27 GiB, read at every poll from `work/partB/GATE-OVERRIDE.json`.
   - PREREGISTRATION-B §9 said 33 GiB, which was unreachable with the desktop apps open.
   - The owner's `swapguard534.sh` returns the floor to 33 if swap free drops below 400 MiB. It never did during Part C:
     its log shows 27 throughout.
2. **The armed build535.sh waiter, the owner's exception of 2026-09-29T17:50Z** (`work/partC/CHAIN-EXCEPTION.json`).
   - MANDATE-C counts any bash running `build5NN*.sh` as a chain.
   - The coordinator's build535.sh (pid 73717) had been armed since 16:27Z, sleeping until 09-30 06:00Z. It would have
     blocked Part C's gate until ~09:30Z tomorrow.
   - At my request (message 17:48Z), the owner allowed my gate to ignore it only while its log's last line still reads
     "armed" and only before 09-30 06:00Z.
   - All Part C heavy work ran 17:55–18:25Z. Its log still read "armed" throughout, and no build535 step had started.
3. **Unchanged:**
   - memory_pressure ≥ 40 % free;
   - no builder and no other chain;
   - strictly serial;
   - each heavy process under `/usr/bin/time -l`, logged to `work/partC/heavy.log`;
   - MANDATE-C's windows (none reached).
   - All 36 heavy starts opened at 27.7–30.3 GiB free+inactive (`gate.log`).

| runs | heavy processes | when (UTC) | wall | peak RSS (ru_maxrss) |
|---|---|---|---|---|
| primary (§1–§4) | 18: 3 full-panel + 15 fold (efr: 2 passes) | 17:55:42–18:07:25 | 700 s | full 7.79 / 7.97 / 7.85 GiB; folds 4.61–7.82 GiB |
| secondary `--pk` (§5) | 18 | 18:11:44–18:24:41 | 775 s | full 5.51–5.73 GiB; folds 4.79–6.08 GiB |
| aggregation (4,000-draw bootstraps) | light | 18:08, 18:26 | 5 s each | 0.12 GiB |

Scoring stayed ≤ 8 GiB, under the 9 GB cap. There were no operator fits. Disk used: 0.43 GB under `work/partC/`.

## 7. Files (`work/partC/`)

- **Instrument:**
  - `guard534.py` (the write guard)
  - `lib534c.py`
  - `designs534c.py`
  - `comp534c.py`
  - `full534c.py`
  - `foldlevel534c.py`
  - `price534c.py`
  - `gate534c.py`
  - `drive534c.py`
  - `aggregate534c.py`
  - `aggregate534cpk.py`
- **Outputs:**
  - `results534C.json`: the decision table, per panel and per fold, the KA summary, ret_f and energy routes.
  - `tables534C.md`
  - `results534Cpk.json` and `tables534Cpk.md`: the secondary.
  - `state534c.json`, `heavy.log`, `gate.log`, `logs/`
  - `<panel>/{comp534c*-<fold>.json, arrays534c*-<fold>.npz, full534c*.json, level534c*.json, fold-*-held.npz}`
- **Reproduce:**
  - `uv run python experiments/534-generalization-audit/work/partC/drive534c.py`, then `aggregate534c.py`. It is gated
    and skips done steps.
  - Add `--pk` to the driver, then run `aggregate534cpk.py`, for §5.
