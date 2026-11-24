# 362 — The oracle common vector buys +0.056 sign accuracy and pays −0.37 pds for it; the route is dead for every line

**Status: NEGATIVE (complete measurement; all six arms emitted, scored, bootstrapped).**

## Known-answer reproduction

Unmodified cap15, run through this directory's own emission replica, reproduces
`experiments/323-current-family-validation/work/metrics-cap15.json` **bit-exactly**
(`work/known-answer.json`, `pass: true`): emitted-count digest `547cf0ab…`,
panel digest `f699b3d5…`, all six raw members with absolute difference **0.0**,
mse_legs exact, every stats array (`targets/sums/jk/p/q/lfc`) identical, every
per-target row identical. Design-stage check first: s_oracle recomputed as **2.7549**,
matching 345 to 4 decimals, and the design wacc grid reproduces 345 exactly
(cap15 0.6002, α=0.5 0.6527, α=1 0.6564). The construction is right; what follows is
the idea failing, not the plumbing.

## What was measured (inputs, estimator, controls)

- **Inputs (all read-only):** `stage345.npz` (`Cmean` 85×10,779 LOO same-line common
  vectors, `Dfull` cap15 design on the universe, `genes`, `targets`, `elig` 58 indices,
  `gate`, `R`); `effects.npz::cap15` (85×18,533, 319 composition) with universe genes
  mapped **by name** (10,777/10,779 land; `TBCE-1`/`HSPA14-1` have no slot and are
  dropped — their Dfull entries are 0, sign 0 = miss, the 343/345 convention; verified
  `max|mapped cap15 − Dfull| = 0.000e+00`); `jurkat_common.npz` (350; mean log2-CPM-ratio
  over 300 seeded Jurkat targets, 8,284/18,533 challenge-mapped → **7,560/10,779** on
  the universe, 0 where unmapped).
- **Scales (345's rule, pooled-median |D| / pooled-median |C| over the eligible gate,
  zeros included):** s_oracle **2.7549** (= 345 exactly; my first attempt filtered zeros
  and got 2.6935 — corrected before any emission), s_generic **1.1085**.
- **Arms** (universe genes only; all other effects-axis genes bit-identical to cap15):
  `cap15`; `oracle_a05`/`oracle_a1` = D + α·s·Cmean, α = 0.5/1.0; `oracle_a1_normpreserve`
  = oracle_a1 rows rescaled to the original row norms (factors n0/n1 median 0.57 —
  the added component carries ~1.75× the design's energy); `generic_a1` = D +
  s_generic·J (same Jurkat vector for all targets); `shuffled_a1` = oracle_a1 with one
  global gene permutation π (seed 362) applied identically to every Cmean row.
- **Emission + scoring:** own-directory replica of `evaluate323.emit_and_score`
  (production builder, H1 union OFF; 400 cells, min_cpm 10, tail_spread 1.0, keep_top
  8000, keep_spread 0.25, block map [0, 0.60], dispersion 0.1, overdispersion 4.0,
  CallSet(400, 2200, 1.9) first, jitter 1000+i, seed 0), `E.block_stats` + `E.score`
  per arm, sequential in one process. Paired 10,000-target bootstrap vs cap15
  (seed 20260909, MSE panel-ratio recomputed per resample, report342 pattern).

## Results

### Design stage (343's |lfc|-weighted sign accuracy on the gate, mean over 58 targets)

| arm | wacc | Δ vs design |
|---|---:|---:|
| cap15 (design) | 0.6002 | — |
| oracle_a05 | 0.6527 | +0.0525 |
| oracle_a1 | 0.6564 | +0.0562 |
| oracle_a1_normpreserve | 0.6564 | +0.0562 (signs identical — rescaling is sign-preserving, as intended) |
| generic_a1 | 0.5842 | −0.0160 (below design; PI's 0.534 re-derivation confirmed in combination) |
| shuffled_a1 | 0.5880 | −0.0122 |

### Six raw members after emission (units: raw metric values; mse/nmae lower-better)

| arm | pds ↑ | mse ↓ | nmae ↓ | fid ↑ | reach ↑ | jac ↑ |
|---|---:|---:|---:|---:|---:|---:|
| cap15 | 0.88964 | 0.97495 | 0.95084 | 0.49755 | 0.25435 | 0.07476 |
| oracle_a05 | 0.59622 | 1.07804 | 0.90918 | 0.51348 | 0.28166 | 0.07847 |
| oracle_a1 | 0.52283 | 1.59107 | 0.91454 | 0.52085 | 0.25884 | 0.08069 |
| oracle_a1_normpreserve | 0.52171 | 1.07559 | 0.92223 | 0.51882 | 0.25814 | 0.08070 |
| generic_a1 | 0.64860 | 5.19962 | 1.19418 | 0.50514 | 0.13505 | 0.07173 |
| shuffled_a1 | 0.84062 | 1.90675 | 0.96961 | 0.49298 | 0.22650 | 0.07354 |

Raw Δ vs cap15: oracle_a05 pds **−0.293**, mse +0.103, nmae −0.042, fid +0.016,
reach +0.027, jac +0.004; oracle_a1 pds **−0.367**, mse +0.616, nmae −0.036,
fid +0.023, reach +0.004, jac +0.006; normpreserve pds **−0.368**, mse +0.101,
nmae −0.029, fid +0.021, reach +0.004, jac +0.006.

### H1-local scaled overall + paired bootstrap (10,000 target resamples)

| contrast vs cap15 | Δ overall (point) | 95% CI | verdict |
|---|---:|---|---|
| oracle_a05 − cap15 | −0.08824 | [−0.11514, −0.06552] | entirely negative |
| oracle_a1 − cap15 | −0.11942 | [−0.15249, −0.09033] | entirely negative |
| oracle_a1_normpreserve − cap15 | −0.12260 | [−0.15393, −0.09567] | entirely negative |
| generic_a1 − cap15 | −0.18327 | [−0.22623, −0.15034] | entirely negative |
| shuffled_a1 − cap15 | −0.03064 | [−0.04738, −0.01812] | entirely negative |

0 non-finite resamples in every contrast; all MSE/bootstrap identity checks pass
(identity abs ≤ 1e-10).

### Kill-criterion accounting (scaled-unit point deltas)

| arm | pds loss (scaled) | nmae+fid+reach gains (scaled) | ratio |
|---|---:|---:|---:|
| oracle_a05 | −0.661 | +0.072 +0.019 +0.032 = +0.123 | **5.4×** |
| oracle_a1 | −0.827 | +0.062 +0.029 +0.005 = +0.096 | **8.6×** |
| oracle_a1_normpreserve | −0.829 | +0.049 +0.026 +0.004 = +0.080 | **10.4×** |

## Verdict against the PLAN.md kill criterion

**The route is dead for every line — record it so.** Neither requirement is met, and it
is not close: (1) no oracle arm improves the H1-local overall — all three are
−0.09…−0.12 with paired intervals entirely below zero; (2) pds pays **5–10×** what the
DE members gain in scaled units *even for the oracle*, which uses other targets' H1
truth and is therefore the ceiling of anything buildable. A lower dose does not save it
(α=0.5 still pays 5.4×), and removing the added energy does not save it (normpreserve
pays 10.4× — the damage is directional, not energetic). `generic_a1` moves plenty and
all of it harmfully (mse +4.22 raw, nmae +0.24 raw, reach −0.12 raw, pds −0.24 raw; it
improves **neither** mse **nor** pds, with or without sign gain — design wacc 0.5842,
below the design alone). Do not build a same-line common-vector arm for A or any line,
in additive form at any dose, with or without norm preservation.

**Mechanism (one line):** adding any shared vector compresses across-prediction
differences, and pds is pure retrieval among 300 — the oracle rows are dominated by the
common component (median energy ratio 1.75×), so all 85 predictions point nearly the
same way and discrimination collapses to near-chance (0.52 vs 0.51 baseline). This is
§4.2's whitening lesson in reverse: whitening strips shared directions and pds rises;
here we inject one and pds dies. The +0.056 sign gain is real but pds-incompatible in
additive form. corroboration from the controls: shuffled_a1 (energy, no information)
costs pds −0.049 raw and mse +0.93 raw while gaining nothing — so the oracle's small DE
gains *are* the information (oracle_a05 nmae −0.042 vs shuffled +0.019), swamped 5–10×
by the retrieval cost of carrying them this way.

Two further notes for the record: (a) the H1-local scaled overall **hides** a raw mse
cost via the clip (all arms clip to scaled 0.0 including cap15) — oracle_a1 adds +0.616
raw mse and generic_a1 +4.22 raw; on A/B/C, where our mse member is positive, the true
price is even worse than the −0.09…−0.12 printed here. (b) The scorer issued an
across-perturbation-spread warning on normpreserve only (spread/claim ratio 0.48 < 0.7;
bounded correction, score not inflated) — another symptom that this construction fights
the emission geometry, not a scoring artefact.

## Limits

H1 development panel, one line, 85 targets / 58 nmae-eligible; H1-local anchors, not an
A/B/C forecast (and the forecast would be worse per note (a) above). The oracle uses
same-panel H1 truth — a ceiling by construction; the transferable version (generic_a1)
was priced and is worse on every member it moves. Only additive combinations at
α ∈ {0.5, 1.0} were tested — the α=0.5 arm already pays 5.4× and the cost grows with
dose while the sign gain saturates (345: 0.6527 → 0.6564 → 0.6573), so smaller α can
only shrink both sides of an already-lost trade; a non-additive use of same-line
information (not tested here) would need its own mandate and its own pds accounting,
because the retrieval tax is on shared *direction*, not on amplitude.

## Exact commands to reproduce

```bash
# design matrices + design-stage wacc (fast; asserts 345's grid to 0.002)
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/362-common-vector-oracle-sixmember/build_designs.py
# emission + scoring, 6 arms sequentially (~15 min; asserts bit-exact cap15 first)
memory_pressure | tail -1  # proceed only if free >= 50%
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/362-common-vector-oracle-sixmember/run362.py
# paired bootstrap report
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/362-common-vector-oracle-sixmember/report362.py
```

All scripts read outside directories read-only and write only into
`experiments/362-common-vector-oracle-sixmember/work/` (`designs362.npz`,
`design_report.json`, `stats-<arm>.npz`, `per-target-<arm>.csv`,
`metrics-<arm>.json`, `known-answer.json`, `emission.json`, `results.json`).

## Peak memory and wall time

Emission run: **910 s wall, 5.22 GiB peak RSS** (6 arms, one sequential process;
per-arm 128–178 s). Design build < 60 s, negligible memory. Bootstrap report ~2 min.
Memory ≥ 70% free at every launch (`memory_pressure`); cap was 12 GiB. No scratch
outside `work/`; no downloads, no cloud, no reserve paths.
