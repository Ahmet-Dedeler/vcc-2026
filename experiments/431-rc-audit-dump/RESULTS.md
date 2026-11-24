# 431 — per-block rc diagnostics dump (`VCC_RC_AUDIT_DUMP`)

Status: patch prepared 2026-09-23 07:20Z. Applied with 429's zs patch after the 428k build; first used on 429a's build.

**Why.** A build's meta keeps only the median, min and max of each rc diagnostic over the 900 blocks. The open question
it cannot answer is where the A/B/C MDE dose binds: 426a's `mde_lambda` has median 0.84 and min 0, against 1.0 on
every H1 target median. This option writes one JSON line per block (context, target, every numeric rc diagnostic).

**What it does.** Diagnostics only. The emission is untouched, so a knob-on build's X equals the knob-off build's. It
writes no config key.

**Files.** `patch_rcdump431.py` (anchored, refuses a double patch), `analyse431.py` (per-context distributions and
Spearman of `mde_lambda` against block features).

## Result: 429a's build (2026-09-23 09:10Z; 428h + zs; MDE 0.1,0.1,0.3 without `coarse`)
`work/rc_audit_429a.jsonl` in `experiments/429-zs-production/`, 900 lines (one block, A/CHKB, has no MDE keys).

| context | `mde_lambda` median | q10 | share at 0 | `prior_lambda` < 1 | sub-gate counts/cell | clip_frac | forced |
|---|---:|---:|---:|---:|---:|---:|---:|
| A | **0.000** | 0.000 | 0.92 (max 0.742) | 0.32 | 75.1 | 0.029 | 3,745 |
| B | 0.936 | 0.866 | 0.02 | 0.017 | 85.5 | 0.133 | 7,651 |
| C | 0.855 | 0.762 | 0.01 | 0.010 | 88.7 | 0.099 | 7,898 |

- **Context A gets no MDE.** It is the context with the least sub-gate slack per cell, and even the prior already binds
  on a third of its blocks.
- Within B and C, the dose follows the slack: Spearman(`mde_lambda`, sub-gate counts per cell) +0.62 (B) and **+0.83**
  (C); with `per_cell_total` +0.59 / +0.81; with `clip_frac` +0.69 / +0.61. The flip count and the up/down gene counts
  explain almost nothing (|rho| ≤ 0.23).
- So the binding constraint is class capacity, as the H1 diagnosis said: `split_for_A`'s polish drains the middle
  library class. The fix is 426h's `coarse` delivery (count-balanced by construction), in production since 5568ee84.
  433a (417e + coarse MDE) carries the dump and gives the same table for the fixed emitter.
- Across contexts the pooled Spearman (mde_n_dn −0.70, forced +0.61) is mostly the A-vs-B/C split. Use the per-context
  numbers.
