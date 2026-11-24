# 393 — mse shared component: NEGATIVE (shared energy 0.85% of A; best case 0.0087 raw < 0.01 kill line)

Status: **NEGATIVE** (complete except the H1 harness, which was deliberately not run —
the Task 3 falsifier answered first; see Verdict). Zero harness runs used of the two allowed.

## Known-answer reproduction (own code, `ka393.py` → `work/ka393.{json,log}`)

| # | row | reproduced | record | verdict |
|---|---|---:|---|---|
| KA-i | C_p-free LS solve a8+a10+a11 → X, ΣD | 0.019240 / 0.229141 (bit-identical to 211's `mse187.json` rerun, \|d\| = 0.00e+00 both) | 0.019240 / 0.229141 | **PASS** |
| KA-i residuals | +5.99e-05 / −6.12e-05 / +1.69e-06 | 211 prints +6.0e-5 / −6.1e-5 / +0.2e-5 | **PASS** (equal at the record's own 2-s.f. reporting; full-precision gate used 6.5e-5) |
| KA-ii | q_c(a15) recomputed from stored per-row Ĉ (`chat211.json`, MAXT24 sel) | 1.000000 / 0.999904 / 1.000000 | 1.000000 / 0.999904 / 1.000000 | **PASS** ×3 to 6 decimals |

Two honest footnotes, both resolved in the log: (1) the middle LS residual is
6.12e-05 in full precision vs the mandate's 6.1e-05 — 211 itself rounds it to
6.1e-05, and X/ΣD are bit-identical, so this is print-rounding, not a mismatch.
(2) q must be recomputed on the MAXT24 24-row subset (`sel = linspace(0, 299,
24)`); the stored dict is the 8+24 union (30 rows), whose mean is 0.999774 on B.

## What was measured (inputs, estimator, controls)

- **Task 1** (`decomp393.py`): streamed the shipped 319p2 artifact
  (`experiments/319-cap-variant/work/build_p2/the-photocopier-demands-a-raise.h5ad`;
  360,000 rows = 900 blocks × 400, 18,533 genes) per (context, target) block
  via h5py bincount — never densified. Per block: pseudobulk total,
  μ = log1p(5e4·tot/Σtot), d = μ − b_ctrl(context), b_ctrl from the pooled
  real controls (`data/challenge/context_{A,B,C}.h5ad`, gene axis verified
  **identical order** to the artifact on all three contexts). Per context:
  w = scorer's own `_row_weights`, A_c = Σw(‖d‖² − target-coord²),
  B_c = scorer's own `_across_pert_budget`, d̄ = weight-mean effect,
  S_c = (Σw)·‖d̄‖², cos(d̄) across lines.
- **Control**: A_c/B_c land within +0.46%/−0.19%/−0.24% (A) and
  +0.53%/−0.14%/−0.19% (B) of 211's rc_callset values on the same instrument;
  B/A sits at 0.9913–0.9918 (211's 0.990–0.999 band); wsum = 0.01618821 =
  300/18532 exactly (300/300 targets resolve); A−B matches S_c to 1.3%
  (target-exclusion residue, as expected).
- **Task 2**: our 12 scored board entries tabulated (mse_raw 0.8094–0.8218);
  closest effect-changing pair 319p→319p2 (partial G 0.34→0.5) moves the whole
  profile (A, B, X jointly) — one equation in two unknowns per context, and no
  historical pair varies only the shared scale. X_shared is **not identified**
  from board pairs; no restores performed (nothing to identify with them).
  Bracket adopted per mandate: X_shared ∈ [0, S_c] (c* ∈ [0, 1]).
- **Task 3** (`price393.py`): closed-form Δraw(c) = mean_c[(c²−1)S_c −
  2(c−1)Xs_c]/ΣD, ΣD = 0.229141, c ∈ {0,…,1.25} × Xs/S ∈ {0, 0.5, 1}.
- **Task 4**: not run (0 runs). Justification in Verdict.

## Results tables (board units)

Per-context decomposition of 319p2 (scorer bulk_lognorm space, TS = 5e4):

| ctx | A_c | B_c | B/A | S_c = P‖d̄‖² | (A−B)/A (shared share) | S_c/A (direct) |
|---|---:|---:|---:|---:|---:|---:|
| A | 0.237515 | 0.235566 | 0.99179 | 0.001974 | 0.8205% | 0.8309% |
| B | 0.225984 | 0.224023 | 0.99132 | 0.001988 | 0.8679% | 0.8796% |
| C | 0.236052 | 0.234041 | 0.99148 | 0.002035 | 0.8517% | 0.8623% |

Shared-component commonality across lines (mandate's parenthetical — yes, but moot):

| pair | cos(d̄_x, d̄_y) |
|---|---:|
| A–B | 0.8151 |
| A–C | 0.8113 |
| B–C | 0.8381 |

Δ raw mse (panel mean over contexts; negative = improvement):

| Xs/S (c*) | c=0 | c=0.25 | c=0.5 | c=0.75 | c=1 | c=1.25 |
|---|---:|---:|---:|---:|---:|---:|
| 0.0 (c*=0) | **−0.008724** | −0.008178 | −0.006543 | −0.003817 | 0 | +0.004907 |
| 0.5 (c*=0.5) | 0 | −0.001636 | −0.002181 | −0.001636 | 0 | +0.002726 |
| 1.0 (c*=1) | +0.008724 | +0.004907 | +0.002181 | +0.000545 | 0 | +0.000545 |

H1 six-member table: **not measured** — no harness run was started (see Verdict).

## Verdict against the PLAN.md kill criterion

**KILL, clause 1.** The bracket's best case (X_shared = 0, c = 0 — i.e. the
truth's mean effect orthogonal to ours, the most favourable alignment the
bracket allows) gains **0.0087 raw mse < 0.01**. In score units that is
+0.0092 scaled = **+0.0015 overall**, and it requires deleting the shared
component entirely. If c* = 1 (shared as aligned as its energy allows), any
c ≠ 1 strictly hurts. The mechanism reason is one number: the shared
component carries **< 0.9% of predicted energy** on all three contexts, so no
shared-only rescaling can move the numerator by more than ~0.002 against
ΣD = 0.229141. The cross-line cosines (0.81–0.84) confirm the shared component
is largely one common direction — but a common direction with 0.85% of the
energy cannot price cqawesome's 0.711→0.658 move (worth 0.027 overall) or any
material part of it.

Task 4 was therefore not run: per the protocol ("stop when the falsifier
answers") the idea is dead on mse arithmetic before any pds price is relevant,
and the two spared harness runs (~142 s / 5 GiB each) were left to the
production build running on this Mac. The H1 pds cost at c = 0.5 is unmeasured
and moot. For the record: 380's finding (deflating shared principal directions
drops pds at k ≥ 3) is consistent with this — the shared direction matters for
retrieval rank, but its mse lever is capped at 0.0087 raw by construction.

## Limits

- ΣD = 0.229141 is carried as common across contexts (187's assumption (i),
  from a pooled solve; cannot be separated from pooled raws).
- The pricing uses S_c = (Σw)‖d̄‖² (direct); the identity A−B differs by 1.3%
  (target-gene exclusion residue). Using A−B instead moves the best case from
  0.008724 to 0.008611 — same verdict.
- The bracket Xs ∈ [0, S] is the mandate's, not a measurement; a negative
  X_shared (mean predicted effect anti-aligned with the truth's mean effect)
  would raise the ceiling, but nothing in the record supports anti-alignment,
  and the mandate's bracket is the stated decision frame.
- All energetics are on the 319p2 profile; other profiles' shared shares could
  differ, but 211 shows A/B/C all sit at B/A ≈ 0.991–0.992 across emitters and
  profiles, leaving no room for a large shared share anywhere in this family.
- H1 pds cost unmeasured (moot after clause-1 kill).

## Exact commands to reproduce

```
cd experiments/393-mse-shared-component
uv run python ka393.py      # KA-i + KA-ii  → work/ka393.{json,log}
uv run python decomp393.py  # Task 1        → work/decomp393.{json,log}, work/dbars393.npz
uv run python price393.py   # Task 3        → work/price393.{json,log}
```

All with `NUMBA_CACHE_DIR=experiments/393-mse-shared-component/work/numba-cache`
(set inside each script), threads capped at 4. Nothing outside
`experiments/393-mse-shared-component/` written; no restores; no harness runs.

## Peak memory and wall time

| script | peak RSS | wall |
|---|---:|---|
| ka393.py | < 0.5 GiB | ~60 s (096 `measure` reads) |
| decomp393.py | **2.92 GiB** (getrusage) | ~12 s stream + controls (free ≥ 79% at start) |
| price393.py | negligible | < 2 s |
