# 641 — the day-of depth row (UNAPPLIED patch + tests, on top of 554)

Date 2026-10-02, session [537757], at the coordinator's request. Class (AGENTS §11.1c): **scorer mechanics / delivery
guard**. It changes no prediction. It stops the day-of path when D/E/F's controls arrive at a depth the recipe was not
tuned at. No build, no `vcc`, no `cell-eval2`. Live files are untouched; the tests hash them before and after.

## Why

- Arc's 2026-10-01 post ([Behind the Data 2026](https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge-2026),
  Fig. 2) gives the full data's median depth as **29,482 UMI/cell**.
- The released A/B/C controls sit at **20,109 / 19,946 / 20,034**. Means 21,134 / 19,990 / 21,157; q10 9.0–10.6k; q90
  29.6–34.3k.
- The three medians agree to 0.4 %, so Arc downsampled every released context to one ~20k median (its stated target
  was "median 20,000").
- The recipe's depth-tuned constants were all chosen on A/B/C at that depth, above all the DE call-set size. 521: the
  global call-set size is regime-sensitive; K 1,600 is +0.0048 at H1 f 0.35 and −0.015 at f 1.
- If D/E/F come at another depth (an un-downsampled release would be ×1.47), the build still runs and every per-gene
  quantity is re-fitted from the new controls. But the call-set size would no longer sit in the regime it was tuned in.
  This patch turns that into an explicit decision.

## The patch (`patch_depth641.py`)

It applies only AFTER 554: its anchors are 554's inserted lines. It refuses a missing 554, a second apply, a partial tree
or an anchor miss, and it writes atomically.

| file | change | default |
|---|---|---|
| `scripts/def_preflight.py` | per context: `library_median_ratio` = control median / 20,030 (the mean of A/B/C's three medians), always recorded with `library_median_ref` and `library_median_tol` | a FAIL row only when \|ratio − 1\| > **0.15**. `DEF_DAY_ALLOW_DEPTH=1` makes it a passing "ALLOWED" row (JSON `library_median_allowed`) |
| `scripts/def_day.sh` (538 + 554) | the override is logged at start and written to `def_day.final.env`, next to 554's two | logs only when set |

**Why 15 %.**
- A correctly downsampled D/E/F context should read within ~1 % (A/B/C: 0.9958–1.0039).
- An un-downsampled one reads ×1.47. 15 % is ~37× A/B/C's spread and a third of the raw-data gap.
- Inside 15 %, the MDE (∝ depth^−1/2) moves by ≤ 7.5 %. That is inside 554's 10 % dose tolerance.

**What to do on a FAIL.** Read the three medians.

1. If all three are off by the same factor, Arc released at another depth.
2. Then decide (coordinator) whether the board-tuned call-set size still applies. Price the shift on H1's regime
   instrument (553's h1r7 how-to) at the new depth ratio.
3. Then set `DEF_DAY_ALLOW_DEPTH=1` and record the decision.
4. One context off alone is a different problem (a damaged or mislabelled file). Stop.

## Tests (`test_depth641.py -v`): 4/4 OK, 47 s, peak RSS 0.49 GB

- **A, mechanics** (on a 538 + 554 temp tree):
  - `--check` passes and writes nothing; the patch then applies;
  - a second apply, a tree without 554, a partial tree and an anchor miss are all refused, with nothing written;
  - `bash -n` and `py_compile` pass after.
- **B, the row through the real CLI** (synthetic bundles: data/challenge's gene axis and panel, 240 cells per context):

  | median | ratio | rows | rc (554 tree → 641 tree) |
  |---|---|---|---|
  | 20.0k | 1.000 | none, text identical | 0 → 0 |
  | 17.5k | 0.874 | none, text identical | 0 → 0 |
  | 23.0k | 1.148 | none, text identical | 0 → 0 |
  | 29.5k | 1.472 | 3 FAIL rows | 0 → 1; `DEF_DAY_ALLOW_DEPTH=1` → 3 ALLOWED rows, rc 0 |
  | 16.0k | 0.799 | 3 FAIL rows | 0 → 1; ALLOWED → rc 0 |

- **C, the two def_day.sh lines** (bash with a stub `log`): logged and recorded only when the knob is set.
- **D, real data** (`def_preflight`, base 538 + 554 vs patched 538 + 554 + 641, knobs unset):
  - **identical text and exit code on all six bundles 554 used**: data/challenge, 180 def-bundle, efp, efr, R6 S1, R7a T1;
  - ratios A/D 1.0039, B/E 0.9958, C/F 1.0002, so the row is inert on every A/B/C-derived bundle;
  - R7b's as-built Calu-3 controls as every context (median 11,418): 3 FAIL depth rows (×0.570), beside 554's spread rows.
    With `DEF_DAY_ALLOW_DEPTH=1` they become ALLOWED, while 554's spread rows still fail without their own override.

## Order on the day

Apply after 538's applier and 554's patch, any day from 10-05:

```bash
uv run python experiments/641-depth-row/patch_depth641.py --root . --check
uv run python experiments/641-depth-row/patch_depth641.py --root .
```

Then:

- `def_preflight.py --verify` must still read its 66/66 rows: the row is silent at ratio 1.0039.
- 590's patch touches `def_identity` / `def_sameline` only, so it is disjoint.

## One sentence

Arc's blog says the released files were downsampled to ~20k. A day-of check now makes a different D/E/F depth a recorded
decision instead of a silent regime change.
