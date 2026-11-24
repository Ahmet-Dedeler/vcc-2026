# 164 — HANDOFF (transferred mid-run, 2026-09-02 ~04:09 UTC+3)

No submission, no build, no `cell-eval2`, no existing file modified. Everything is in this dir.

## 1. Scripts written

| file | what it does | status |
| --- | --- | --- |
| `blocks.py` | Block census on the REAL challenge controls: the scorer's DE universe (mean-per-cell CPM > 5) vs production's spread gate (**pseudobulk** CPM > 10, `--min-cpm 10.0`), and how `keep_top 8000` splits them. | **ran cleanly**, 11 s → `blocks.log` |
| `realize164.py` | h1-test realization of 7 arms (a15 ×2 seeds, **armA** per-gene phi=4 keep / 0 tail, **armB** tail-spread 0.5, **armC** = a8 flat ×2 seeds), saving reference + predicted `p/padj/lfc` **and the per-target spread vector** so the call set can be split into keep/tail blocks exactly. Derived from 157's `realize_h1test.py`, same code path, same gate convention (so 157's board-validated `reach` calibration still applies). Donors are 095's h1-test artifact = `k562_gwps` + `orion_hek293t`; `h1_2025` never enters. | **RUNNING in background, PID 94461**, log `realize164.log`, expect ~10–12 min (157 did 3 arms/85 perts in 276 s; this is 7 arms). Leave it. |

## 2. Numbers measured so far

**`blocks.py`, real challenge controls (`data/challenge/context_{A,B,C}.h5ad`):**

| ctx | scorer DE universe (cell-CPM>5) | production spread gate (pb-CPM>10) | DE-visible band BELOW the production gate | effect-ranked tail inside the gate (`gate − 8000`) | total scorer-visible tail |
| --- | ---: | ---: | ---: | ---: | ---: |
| A | 9929 | 9135 | **794** | 1135 | 1929 (19.4%) |
| B | 9626 | 8668 | **958** | 668 | 1626 (16.9%) |
| C | 10124 | 9127 | **997** | 1127 | 2124 (21.0%) |

`prod-gated & not DE` is **0** in all three contexts, so gate(10) ⊂ DE(5) exactly.

Two findings from this, both new relative to 157 §5b:
1. 157's tail SIZE (1929/1626/2124) is **arithmetically identical** to mine, but only by coincidence (`gate−8000 + band = DE−8000`). Its **composition** was wrong: **41 / 59 / 47 %** of the tail block is the sub-gate 5–10 CPM band, where the CPM gate sets the predicted effect to **exactly zero**. For that fraction the zero-effect calculation is exact, not an approximation.
2. Production gates the spread vector at pseudobulk CPM 10; 157's `zero_effect.py` used mean-per-cell CPM 5 for both. h1-test realization (157 and mine) also uses 5 — kept deliberately, because that convention is what reproduced two board `reach` points.

**Read from the spec (`docs/reference/vcc2026-metrics-BRIEF.md` §3–4), and it reframes the thesis:**
`fid`'s direction comes from the **lfc** (arithmetic mean of per-cell normalized values), *not* from the Wilcoxon z. The quasi emitter repairs each gene's panel total exactly and gives every cell the same library, so a gene's `m_pred = 1e6·total_g/(n·L)` is (near-)invariant to `(spread, phi)`. **Therefore spread/phi arms change WHICH genes are called, never the predicted DIRECTION of any gene.** 157's "100 % DOWN" is the z sign, not the scored direction.
Corollary worth measuring: at zero effect `lfc_0(g) ≈ log2(pb_cpm_g / cell_cpm_g)`, and the sub-gate band is *by construction* the set where mean-of-ratios most exceeds ratio-of-sums — so the band's scored direction should be systematically DOWN as a pure **normalisation** artefact, with no model content. This is exactly computable from the controls alone and was **not yet run**.

**Board anchors carried in (not measured here):** raw `fid` a15 0.5022 / a8 0.5021, board max 0.5175 (`The_Way`), leader cqawesome 0.4877; `fid` span 0.29725 (1/span = 3.3642); a8→a15 raw `reach` 0.1890 → 0.1811; `mse` scale `score = clip((0.987862879 − raw)/0.950162566, 0, 1)`; solved `X = 0.019240`, `ΣD = 0.229141`, `Kc = 0.005410`, `K = ΣwC_p = 0.24653`, `claim = 0.016188 · C̄_hat`.

## 3. What I was about to do next

1. `zero_lfc.py` — exact, challenge-only: emit the zero-effect profile at `(0.25,φ4)` / `(1.0,φ4)` / `(1.0,φ0)` / `(0.5,φ4)`, and record per gene **both** the BH call and the resulting **lfc sign**, split by the three real blocks above. Gives the artefact call rate *and* its scored direction.
2. Precision decomposition on h1-test from `de164.*.npz`: split the call set by the saved `spread` vector into keep/tail, and measure `k/n_pred` on each. Then reconcile with board 0.5022 through 096's conserved-headroom closed form `p_after − 0.5 = (p0 − 0.5)·n_before/n_after` (the artefact-at-chance hypothesis), **not** by transferring harness precision levels (h1-test headroom is 0.0437 vs the board's 0.0022 — 20×).
3. **`mse` pricing, mandatory (rule 12 / a19).** Reuse `experiments/096-mse-floor/mech_claim_grid2.py` verbatim — it already computes `C_hat`/`claim` per `(spread, phi)` on the challenge panel from cached `work/{a12,a13,a15?}-stats.npz` + `cp-proxy-ntc.npz` + `board_solve.npz`, via `cell_eval2.moments.jackknife_correction`. It needs a per-gene-phi and per-gene-spread variant for arms A/B. **Suspected binding constraint: a15 sits ON the ρ cliff (`VCC_ALPHA 0.982` is the per-context argmin, where `ρ = min(1, B/claim) = 1`), so ANY reduction in `claim` raises raw `mse` by ≈ `Δclaim/0.229141` with no amplitude compensation available (amplitude is out of my mandate).** Arm A removes φ from the tail block, which is the block a12 added *to reach the #247 cap* — price this before recommending anything.
4. Known-answer gates: `a15` must reproduce 157's median `n_pred` 6035 / `reach` 0.1828; `armC` must reproduce 6035-analogue and `reach` 0.2029; `mech_claim_grid2` must reproduce `claim` 0.1186 (0.25/φ0), 0.1668 (0.25/φ4), 0.1922 (1.0/φ0).

## 4. Traps hit / to avoid

- **`min_cpm` is 5.0 in the module and 10.0 on a15's command line** (`MIN_CPM = 5.0` at `make_submission_pm_cd4.py:431`, overridden by `--min-cpm 10.0` and passed into `emission_spread_vector`). Getting this wrong is what made 157's block composition wrong.
- **`emission_spread_vector` initialises the whole vector to `tail_spread`**, so ALL non-gated genes get spread 1.0, not just the low-`|effect|` gated ones.
- Per-gene `overdispersion` is real and needs no builder change (`donor_transfer.py:403-411`, and `_quasi_block` 290-294 for `resample_mix`) — it must be length `G_full` or it raises.
- Arm D (drop artefact genes post hoc) is **not expressible**: we submit counts, not calls. It is a bound only.
- `pred_lfc` is *not* exactly arm-invariant — per-cell library sums differ slightly between arms even with totals repaired — so assert it rather than assume it (095 saw `nmae` invariant only to 4 decimals).
