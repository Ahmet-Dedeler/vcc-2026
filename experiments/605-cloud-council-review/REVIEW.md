# 605 — adversarial review of the three board-bound patches

Reviewed at `2788ccd9` (branch `cloud/605-council-review`, same commit as master). Code and logic only: this checkout has no challenge h5ad, no donor matrices, and no harness npz. Nothing here was applied to `scripts/`.

**Decision: do not apply these three as specified on 2026-10-06.**

- **630** is an H1 harness (`rnup630.py`), not an unapplied builder patch. The files named in the review brief (`src_rank_up.py`, `patch_builder630.py`, `test630.py`, `test_rank_up630.py`, `PREREGISTRATION-DRAFT.md`) are not in the tree. The pass can flip an UP log-fold change to a large DOWN while every stated invariant (per-gene totals, per-cell libraries, `|ΔAUC| ≤ 1e-4`) stays green. v1 already moved fid and jac on H1 seed 1.
- **612** is a real unapplied patch and is safe to apply after 598 as a sign change. It is not a rank-half split: a point mass at the quantile keeps the key on the whole tied mass. The knob is `VCC_DE_MDE_SIGN=clean4kg:<tau>[:<q>]`, not `VCC_DE_MDE_SIGN_KEYGATE`.
- **604** is not in this commit. `experiments/604-identity-prior/` does not exist, and there is no `identity_prior` symbol under `scripts/` or `src/`. The λ = 0 claim cannot be checked, and nothing can be applied from this tree.

Composition **598 → 612 → 630** cannot be built here. `scripts/make_submission_pm_cd4.py` carries neither `598: tercilec` nor `612: key-gated`.

---

## 630 — rank-neutral UP pass

What exists: `experiments/630-rank-neutral-up/rnup630.py` (`up_pass` v1, `up_pass2` v2), `run630.py`, `price630.py`, `RESULTS.md`, and H1 JSON for v1 only (`work/h1/rnup630_s{0,1}.json`). RESULTS says v2 is queued. There is no env knob and no unset path in the production builder, so "unset is byte-identical" is not a property of this tree.

### BLOCKER — a non-positive CPM gain is clamped away, and the gene is kept

`experiments/630-rank-neutral-up/rnup630.py:283` (v2) and `:143–144` (v1).

v2 asks for

```text
k_want = max(k_promo, ceil(rem / max(cw[a] - min(cw[bal]), 1e-300)))
```

`cw = 1/(n · library/1e6)`, so a carrier with a higher library than every ballast cell has `cw[a] < min(cw[bal])`. The denominator collapses to `1e-300`, `k_want` becomes unbounded, and `fill` donates every in-budget ballast count. Counts already above the control max do not increase `ranks_of`, so the `|ΔAUC|` cap still passes. The revert at `:297–302` checks only `abs(du) > capU` and "carrier strictly above cmax". It does not check the sign of the new log-fold change. `got = dA - rem` (`:306`) is then negative and the gene is still `n_moved`.

v1 has the same clamp against the **mean** ballast weight (`:143–144`), so one high-library ballast cell makes the ask unbounded even when other donors would have raised the mean.

**Inputs → output** (run this session, both `up_pass2` cross True/False and `up_pass` promote=True). 8 carriers at library 22,000, 32 ballast at 18,000, 400 controls, one control just above the carrier, ballast holding a large pile of the gene already above that max, base log-fold change `+0.00144`:

| | counts moved | max \|ΔAUC\| | totals exact | libraries exact | log-fold change |
|---|---:|---:|---|---|---|
| v2, cross False | 10,688 | 6.25e-5 | yes | yes | `+0.00144` → `−0.226` |
| v2, cross True | 10,692 | 6.25e-5 | yes | yes | `+0.00144` → `−0.226` |
| v1 | 10,688 | 6.25e-5 | yes | yes | `+0.00144` → `−0.226` |

`capU = floor(1e-4 · 2 · 40 · 400) = 3` (x2 units). The promotion crosses one control (`ΔU_x2 = 2`) and the extra ~10k counts add nothing to U. The asserts at `:313–322` (and v1 `:172–186`) all pass.

On a production block (`n = 400`, `n_ctrl = 18,400`) the same cap is `floor(1e-4 · 2 · 400 · 18400) = 1472` x2-units, about 736 controls. RESULTS records a median promotion of ~116 controls (`ΔU_x2 ≈ 232`), which fits the cap with room to spare. The bug fires whenever the chosen carrier's library is at least the highest ballast library: selection (`:259`) picks the carrier nearest the control max, so a single high-library cell that happens to hold the top value is enough.

`rank_mag.postpass` already refuses this for the down pass (`src/vcc2026/rank_mag.py:354–356`). 630 does not.

**Fix:** donate only from ballast cells with a strictly smaller library than the carrier, and revert the gene if the mean falls or the sign changes.

```python
pos = cw[bal] < cw[a]
dcap = np.where(pos, dcap, 0)
if not pos.any() or k_promo > int(dcap.sum()):
    # revert / skip
    ...
gain = cw[a] - float(np.min(cw[bal][pos]))
k_want = max(k_promo, int(np.ceil(rem / gain)))
# after the move, before accepting the gene:
if np.sign(l_new) != np.sign(l_base[j]) or A_new < A[j]:
    # revert snap_col / snap_sub / sub_moves
    ...
```

Apply the same donor filter in v1's second loop (`:143`), which should use per-donor `cw[bal][bi]` and never `mean(cw[bal])`.

### MAJOR — the call guard is fail-open, and v1 has none

`rnup630.py:375–397`.

The v2 loop excludes genes whose BH call flipped and re-runs, then

```python
if not flips or guard_rounds >= 5:
    break
```

If flips remain at round 5 it breaks anyway. `:396–397` only `print`s `WARNING` and the block is still scored (`:398–400`). A dirty block is a successful return. The v1 regex path (`:384–392`) never enters the loop.

H1 v1 JSON (this is the code that has actually been scored; v2 has not):

| | seed 0 `up_k100_p` | seed 1 `up_k100_p` |
|---|---|---|
| calls changed (mean over 85 targets) | 0 | 0.01176 (= 1/85, MED19 in RESULTS) |
| fid | 0.5089291277653674 → same | 0.5090793609138005 → 0.5090815061920431 (`+2.15e-6`) |
| jac | 0.07469853420125332 → same | 0.07472291848297988 → 0.07472243235401653 (`−4.86e-7`) |
| nmae | 0.9408535485045936 → 0.937020262814194 | 0.9409206615088422 → 0.937192053234356 |
| pds, mse, reach | identical at printed precision | identical at printed precision |
| jk_delta (mean) | +3.338 | +3.380 |
| max \|ΔAUC\| | 5.82e-5 | 5.85e-5 |
| max \|Δlog10 p\| | 0.00447 | 0.00447 |

RESULTS.md:38 prints seed-1 jac as `−0` and line 48 says "Only nmae moves." The JSON refutes both. `up_k100` (no promotion) is an exact no-op on both seeds: `frac_got = 0`, `jk_delta = 0`, six members identical. `zsKA_matches_run630_metrics` is true on both seeds.

I did not find a BH spillover in which an untouched gene flips and the touched gene does not. Step-up BH rejects a prefix, so decreasing an uncalled p and crossing α rejects that gene itself; increasing an already-non-rejected p does not shrink the set. The fail-open cap is still a hole: a chain longer than 5 rounds returns the last dirty block. **What to check on the dumped blocks** (not in this checkout): whether any v2 run exits with `calls_changed > 0`.

**Fix:** if flips remain after the loop, return the pre-pass block (the zs block), or `assert dg["calls_changed"] == 0` before scoring. Do not ship a warning.

### MAJOR — mse is not an invariant of fixed totals and libraries

`pds` is. `mse` is not.

`pds_cosine` is the cosine of log-pseudobulk of summed counts, panel genes dropped, scale-invariant. `up_pass` / `up_pass2` assert `X.sum(0)` unchanged (`:176`, `:315`). That pins pds. Both H1 seeds agree.

`expr_mse_unbiased_capped_norm` credits `rho * min(C_hat, C_real)`, and `C_hat` is the delete-one jackknife over **cells**. Moving counts between a library-22k cell and a library-18k cell changes the jackknife even when the column sum does not. The H1 summaries record `jk_delta ≈ +3.3` with mse repr-identical, which means on **this** panel `min(C_hat, C)` stayed equal to `C` for every perturbation whose rho is nonzero. The break is any perturbation whose `C_hat` crosses `C`, or whose rho denominator moves. **What to check:** per-target `C_hat` versus `C` before and after, on the dumped blocks. Do not write "mse exactly unchanged" into a preregistration off the column-sum assert.

### NOTE — reach is not pinned by "0 call flips"

Reach sorts the reference's DE genes by (significant first, ascending `p_adj`, ascending `p`, descending `|lfc|`, gene name). A promotion changes `p` (the cap is on pairwise U, and the normal-approximation z also has a tie term; below) and the pass exists to change `|lfc|`. Both H1 seeds left reach numerically equal, including the seed-1 call flip. That is a measurement on those 85 targets, not a structural pin. The board read RESULTS already asks for is the right check.

### NOTE — `|ΔAUC|` does not pin the scorer's p

`ranks_of` (`src/vcc2026/rank_mag.py:161–166`) returns `2·#(controls below) + #(controls equal)`. Sum/2 is Mann–Whitney U against controls only. `dAUC = d(ranks_of sum) / (2·n·n_ctrl)` and `capU = floor(max_dauc · 2 · n · n_ctrl)` (`rnup630.py:231`, `:189`) match that definition. Pred–pred ties do not change this U if no control is crossed.

The scorer's z does change. `wilcoxon_z` (`src/vcc2026/rc_emit.py:106–117`) adds `(tot³ − tot) − (cc³ − cc)` for every tied pred run, including ties that involve no control. Two blocks with the same U can have different z. The measured max `|Δlog10 p|` of 0.0045 on H1 is this gap (plus the promotion's real U change). A gene sitting just above q = 0.05 can enter the call set inside the AUC cap. That is why a guard is required, and why the guard has to be closed (the MAJOR above).

### NOTE — unknown variant names crash; `s_k100` is K = 1.0

`rnup630.py:369–386`. A name that matches neither `(s|b)_k(\d+)` nor the v1 regex leaves `m` as `None`, and `m.group(1)` raises `AttributeError`. There is no K-range check (`rank_mag` restricts k). `s_k100` / `b_k100` parse as `K = 100/100 = 1.0` and a cross flag. Fine if that is the dose; say so in the preregistration, because the variant string looks like a percentile.

### What 630 gets right

- **Libraries and per-gene totals**, on every block the function returns. Compensation is `ZS._take` on the sub-gate (`~keep`, keep = control mean CPM > 5) for the same number of counts (`:167–170`, `:308–311`), then `X.sum(1)` and `X.sum(0)` are asserted. The sign-flip examples above preserved both. This holds only if the `keep` passed in is the scorer's keep. The builder already logs disagreements between `rc_scorer_control` and `rc_control`; a future patch must call this on `rc_scorer_control`, after `rank_mag.apply`.
- **Placement relative to zs.** `postpass` selects `dn` only (`rank_mag.py:203–204`), never touches ballast rows (`:372–373`), and has its own sign guard. 630's `main` re-applies zs to the pre-zs dumps and then runs the UP pass (`rnup630.py:355–376`). That order is the right one: zs and 630 share the carrier sub-gate budget, so 630-before-zs would be a different pass. The re-applied zs matches run630's sums, p, q, lfc, and jackknife on both seeds.
- **Strict `up_k100`** (no promotion) is an exact no-op, as RESULTS says. Zero capacity is real.
- **The AUC formula's x2 units** match `ranks_of`. A one-control crossing on the synthetic block gave `|ΔAUC| = 2/(2·40·400) = 6.25e-5`.

---

## 612 — key-gated hook

The patch is `experiments/612-key-gated-hook/patch612.py`. It edits `scripts/make_submission_pm_cd4.py` in six anchors, refuses a builder without `598: tercilec` (`patch612.py:117–118`), and refuses a second apply (`:115–116`). The emitter is not touched. The live builder has neither mark, which is the intended pre-apply state.

The knob is **`VCC_DE_MDE_SIGN=clean4kg:<tau>[:<q>]`** with `0 ≤ tau < 0.5` and `0 < q < 1`, `q` default 0.5 (`patch612.py:53–63`). There is no `VCC_DE_MDE_SIGN_KEYGATE`.

### MAJOR — "top half" is a linear quantile, and ties at the cut stay with the key

`patch612.py:89–93` (same rule in `run612.py:99–103` for legal `q < 1`):

```python
cut = float(np.quantile(k[nz], q))          # numpy default: linear
return np.where((~nz) | (k < cut), s, 0.0)  # k == cut is strong
```

`k < cut` puts every gene whose `|key|` equals the cut on the strong side (hook dropped, delivered sign kept). A point mass at the median assigns **none** of the nonzero keys to the hook.

Run this session, `q = 0.5`, fraction of **nonzero** keys that receive the hook:

| `|key|` | cut | nonzero-weak fraction |
|---|---:|---:|
| `[1, 1, 1, 2]` | 1.0 | 0 |
| `[3, 3, 3, 3]` | 3.0 | 0 |
| `[1, 1, 1, 1, 5]` | 1.0 | 0 |
| `[1, 2, 3, 4]` | 2.5 | 0.5 |
| `[1, 1, 2, 2]` | 1.5 | 0.5 |

`[0, 0, 1, 1, 1, 1]` keeps the two zeros weak (correct) and the four tied nonzeros strong (nonzero-weak fraction 0). D/E/F design keys can repeat: the gate is a pooled profile, and many genes share a rounded or sparse magnitude. `test612.py` T4a (`:144–157`) locks the helper to `run612`'s gate, including rounded ties, so the tests pass while the fraction is 0. The h1 "median 0.496 of covered genes" number is in RESULTS only; T4b, which would recompute it, did not run here (MINOR below).

**Fix:** a rank split with a stable order, and an explicit tie rule. Zeros stay weak.

```python
nz_idx = np.flatnonzero(nz)
order = nz_idx[np.argsort(k[nz_idx], kind="stable")]
n_weak = int(np.floor(q * order.size))
weak = np.zeros(k.shape, dtype=bool)
weak[~nz] = True
weak[order[:n_weak]] = True          # boundary ties split by stable index, weak side first
return np.ascontiguousarray(np.where(weak, s, 0.0))
```

### MINOR — malformed `clean4kg` dies in the pre-existing outer `try`, with the old error

`scripts/make_submission_pm_cd4.py:789–793` parses `parts[2]` as `tau_c4` for every 3-field non-file source **before** the clean4kg branch the patch inserts after `:796`.

`clean4kg:x`, `clean4kg:0.05:y`, and `clean4kg:0.05:` raise `ValueError` inside that `try` and exit with `must be 'clean4:<tau>' or 'programme:<tau>[:<tau_c4>]'`, not the clean4kg message. They are still refused (`SystemExit`). Legal forms are unaffected: `clean4kg:0.05` has two fields, so the outer float of `parts[2]` is skipped and the new branch sets `q = 0.5`.

**Fix:** in the patch, replace the outer `try` anchor as well, and parse `clean4kg`'s third field only inside the clean4kg branch.

### MINOR — `test612.py` aborts the process when the harness npz is absent

`test612.py:158–161`. T4b calls `lib521.load_state(0)` and then opens `experiments/323-current-family-validation/work/reference.npz`. Neither file is in git. The exception is uncaught, so T5 and T6 never record. The module docstring describes the test as light and data-free; T4b is not.

**Fix:** catch `FileNotFoundError`, record SKIP, and continue with T5/T6.

### NOTE — the builder's weak-half source is pinned CLEAN4, not the programme and not the identity score

Per-context (`patch612.py:98–99`) `clean4kg` takes the same branch as `clean4`: the full CLEAN4 sign on `rc_control.keep`. It does not enter the programme/file arm, so it does not read `pdown` and does not set `de_mde_sign_c4`. The per-target expression (`:105–108`) then masks that CLEAN4 vector.

`run612.py` has a separate `idkg` mode that the builder patch does not implement. Stacking 604/639 on this patch does **not** put the identity score on the hook. If the October recipe wants "key's sign on the strong half, identity's sign on the weak half", that is a different patch. 612 as written is "key's sign on the strong half, CLEAN4's sign on the weak half".

### NOTE — key zeros are weak, and the key is not the scorer's 5 CPM filter

The key at the per-target call is `where(emission_rank_axis(row, b_ctrl, min_cpm=args.min_cpm), row, 0)` (builder around `:5969–5976`; `--min-cpm` default 5.0 at `:1941`). That is a pooled-profile CPM gate, not the scorer's mean per-cell CPM > 5. Every zero on the keep axis is weak (`~nz`), so the hook covers the zeros **plus** the bottom quantile of the nonzero keys: more than half the DE universe whenever many gated genes are zero. **What to check on the recipe:** the shipped `--min-cpm` (some recipes pass 10) and the fraction of 5-CPM-universe genes with key 0.

If `VCC_SAMELINE` blends rows (`w > 0`), the gate sees the blended `|row|`, not the pure design key. The October plan uses `w = 0`, which avoids this. A nonzero `w` changes which genes are "strong" without any change to `q`.

### What 612 gets right

- **Unset, `clean4:<tau>`, and `file:...` stay on their old paths.** The dict ternary (`C_NEW`, `patch612.py:72–73`) and the per-target ternary (`F_NEW`, `:105–108`) only read `_ms_kg_q` / call `_kg_mask_sign` when `source == "clean4kg"`. Python does not evaluate the unused branch. A temp-tree check of the patched per-target expression returned the same array object for `clean4` and `programme`/`file`, the CLEAN4 bypass object for a file source on a weak target, and `[1, 0, 0, 0]` for `clean4kg` on key `[0.1, 2, 0, 3]` (the two zeros and the below-cut `0.1` take the hook; `2` and `3` stay 0). Mark count after apply is 6, one per edit, and a second apply is refused.
- **Parse of legal values.** `clean4kg:0.05` → `{source, tau: 0.05, q: 0.5}`. `clean4kg:0.05:0.25` and the uppercase form → `q = 0.25`. `q = 0`, `q = 1`, `tau = 0.6`, and a 4-field value are refused inside the new branch (`:61–63`).
- **All-zero key returns the full hook** (`:90–91`). There is no strong half, so keeping the key would be keeping nothing. One nonzero gene is strong (`k < cut` is false for the only nonzero value, since the quantile of a single point is that point); the zeros still take the hook.
- **Order.** The patch will not apply before 598. 612 itself changes signs, so it moves calls; it does not claim rank neutrality. It runs at emission (`mde_add_sign` on the complement only, 0 = keep the delivered sign), which is before ballast, zero repair, and zs. 630, if it is later inserted after zs, selects `lz > 0` on the **post-zs** log-fold changes, so a gene 612 flipped to UP becomes a 630 target. Combined with the 630 blocker, that gene can then be delivered as a large DOWN while totals, libraries, and the U cap still pass. The two patches do not commute, and 630 does not know which sign 612 wrote.

---

## 604 — identity-prior plumbing

### BLOCKER — the patch is not in this commit

`experiments/604-identity-prior/` is absent. So are `def_identity_prior.py`, `patch604.py`, `test604.py`, `build604.py`, and `make_refs604.py`. Grep of `scripts/` and `src/` finds no `identity_prior`.

`docs/NEXT-SLOT-DECISION.md:1616` says artifacts were built (λ 0 equals CLEAN4 exactly; λ 1 flips 1,606–1,843 genes per context) and were being rebuilt on 622's Flex-only references. Those artifacts are not in git. None of λ = 0 exactness, Flex-only references, the quantile map, the de-duplication guard, or "every other path byte for byte" can be confirmed, and nothing can be applied on 10-06 from `2788ccd9`.

The notes below are traps in the code that **is** committed, which a 604 writer will otherwise copy. They are not bugs in a patch this review could see.

### The science code a 604 port would copy

**References are not de-duplicated by identity call.**

- `experiments/639-identity-axis-lolo/ident639.py:60`: `others = [k for k in prof if k != s]`. Exclusion is by profile name. The set includes the public 3′ profiles. 622 measured those as harmful (hek293t about −0.05 to −0.085 on the identity axis).
- `experiments/622-identity-track/refs622.py:52–62`: the Flex set is the name list `["ctxA", "ctxB", "ctxC", "h1"]` excluding `s`. Same name test, no second profile of the same line removed.

"Never the context's own" is true only as "drop this profile's key". Two profiles of one identity under different names both enter. Flex-only is the named set in `refs622.py:52`, not `ident639`'s default.

**The committed score is not a textbook z-score, and there is no quantile map.**

`ident639.py:64–65` and `refs622.py:58`:

```python
zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)    # no centering: sign matches sign(P - 0.5)
zi = (Iv - np.nanmean(Iv)) / np.nanstd(Iv)  # centered
```

Centering `(P − 0.5)` would flip genes relative to CLEAN4 wherever `mean(P) ≠ 0.5`. The uncentered form does not. Copying `sign(zc + λ zi)` is not the delivery plan in `experiments/622-identity-track/RESULTS.md:57–66`, which says quantile-map onto CLEAN4's own P.

A rank map of a strictly monotone function of P, onto those same genes' own P, recovers P at λ = 0, ties included (tied P values are equal). `np.interp` **clips** outside the observed score range, so an identity-only score past the end of the CLEAN4 range is assigned an endpoint P and will move. Filling missing `z(I)` with 0 and mapping every gene assigns order statistics to genes the identity never saw.

**Safe rule, since no implementation exists to review:** per context, quantile-map only genes where both `z(P)` and `z(I)` are finite, onto those genes' own CLEAN4 P; leave every other gene at CLEAN4's P. z-scores are computed on that context's reference intersection, so a single λ is not comparable across contexts and the map has to be per context.

**psign convention, and the silent all-down.**

`scripts/def_sameline.py:140`:

```python
psign = np.where(np.isfinite(pdown) & (pdown <= 0.5), 1.0, -1.0)
```

`P_down > 0.5` → −1 (down). `P_down ≤ 0.5` → +1 (up). **NaN → −1.** An identity writer that leaves NaN for a gene it meant to keep at CLEAN4 delivers down. The hook's own rule (`|P − 0.5| ≥ tau`, else 0 = keep delivered) is a different convention; do not reuse `programme()`'s NaN default for the identity artifact.

**`--programme-source` will reject `identity:<lam>` and would break the matched-context contract if it were forced through.**

`scripts/def_sameline.py:239–259`. The value must already be a `BANK_DIRS` key (`:240–241`) before any special case. `identity:0.5` is a `SystemExit` today. The branch then calls `programme()` (`:122–141`), which is P_down of a perturbation **bank** with panel rows removed, not a control-CPM identity. Routing the identity score through `programme()` answers the wrong question.

What the branch does preserve, and what a new branch must preserve:

- For every context, `programme_<c>` and `pdown_<c>` are overwritten (`:251–252`).
- `rows_<c>`, `signs_<c>`, `targets_<c>`, and `w_<c>` are left alone **only if** `c` was already in `matched` (`:254`).
- A context the same-line loop skipped — no bank, or a registered line whose call is not unanimous (`:180–186`) — gets `w = 0`, empty rows, empty signs, empty targets (`:255–258`).

An identity branch has to be parsed **before** line 240, must not call `programme()`, and must keep that matched/unmatched split. Overwriting `rows`/`signs`/`w` on a matched context changes the design key the builder blends, which changes 612's strong half even at `w` small.

**Board form versus a relabelled day-of bundle.**

`experiments/622-identity-track/RESULTS.md:60–63`:

- board (A/B/C): each context versus the other two A/B/C contexts **plus H1's Flex controls**;
- day-of (D/E/F): each context versus every other Flex context (A, B, C, the other two of D/E/F, plus H1).

A relabelled 3-context bundle that builds "other contexts in this bundle" drops H1 and will not hash-equal a board artifact that included H1. **What to check on the day,** once the code exists: same reference list, same gene mask, same CLEAN4 P, same map, on a relabelled copy of the A/B/C controls, hash-equal to the board npz. λ = 0 equality is about the mapped P, not about `sign(z)`.

---

## Claim checklist

| Claim | Verdict |
|---|---|
| 630 raises one carrier per uncalled UP gene just above the control max, paid by ballast | True of `up_pass2` when a candidate exists (`rnup630.py:249–294`). v1 can raise `max_rec` carriers. |
| Per-gene totals and per-cell libraries exact | True on every returned block (asserted). Preserved in the sign-flip counterexample. |
| `\|ΔAUC\| ≤ 1e-4` per gene | True of v2's accept path, in pairwise-U units. Does not pin Wilcoxon z or BH. |
| 0 call flips, so pds/mse/fid/jac exactly unchanged and only nmae moves | **False.** v1 seed 1: 1 call, fid `+2.15e-6`, jac `−4.86e-7`. Guard is fail-open and v1 has none. pds is pinned by the column sum. mse is not pinned; it happened not to move on this H1 panel (`jk_delta ≈ +3.3`). reach happened not to move. |
| Unset 630 path byte-identical | **Not implemented.** No builder patch, no knob. |
| 612 keeps the key's sign on the top half of nonzero `\|key\|`, hook elsewhere | **False for ties.** Linear `k < quantile`. A point mass at the cut keeps the key on the entire tied mass. |
| 612 no-op when unset | True for unset / `clean4` / `file` / `programme` (same objects, dead branch not evaluated). |
| All-zero key, and the composition order | All-zero key → full hook (correct). 612 applies only after 598. 630 has no builder stage, so 598 → 612 → 630 is not a patch stack. If 630 is later inserted after zs, 612's new UP genes become 630 targets and can be sign-flipped by the blocker. |
| `identity_prior`: `z(P−0.5)+λ z(I)`, Flex-only, never own context, quantile-map onto CLEAN4, λ 0 exact | **Absent.** Cannot confirm. The committed score (`ident639.py:64`) is uncentered on P and centered on I, with no quantile map and no Flex filter. |
| `patch604` adds `--programme-source identity:<lam>` and leaves every other path byte for byte | **Absent.** Current `def_sameline.py:240` rejects any source not in `BANK_DIRS`. |
