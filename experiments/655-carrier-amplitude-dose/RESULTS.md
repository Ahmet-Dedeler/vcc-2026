# 655 — carrier amplitude dose curve under zero repair: interior peak at ×1.5 for f 0.35 and ×1.25 for f 0.5; WEAK CANDIDATE by the preregistered bin

2026-10-02 20:24–21:07Z (runs; commit times), council session [7772dc]. PREREGISTRATION 18f4d2cb was committed before any 655
run. Its run condition, 654's ×1.25 − ×1 > 0 at f 0.35 on ≥ 2/3 seeds, was met 3/3 (cond655.py). Class: **scorer
mechanics** (the DE carriers' design amplitude; knob `VCC_ALPHA`, 654 §3). The dose is read on H1, so a held-out read is
owed before October.

## Result (H1 real emitter, `--zero-repair -5`, paired per seed vs ×1; overall-eq by lib521 / 192's INV)
| dose | f 1 | f 0.5 | f 0.35 | cmp | calls/target | nmae raw (f 0.35) | fid / jac raw (f 0.35) |
|---|---:|---:|---:|---:|---:|---:|---|
| ×1.25 (654) | +0.0001 | **+0.0011** (3/3) | +0.0019 (3/3) | +0.0012 (3/3) | 2,597 (+26) | −0.0042 | −0.00002 / −0.00006 |
| **×1.5** | −0.0005 | +0.0007 (3/3) | **+0.0031** (+0.0024 / +0.0036 / +0.0032) | **+0.0019** (3/3) | 2,627 (+57) | −0.0066 | +0.0007 / −0.0003 |
| ×2.0 (2 seeds) | −0.0063 | −0.0022 | +0.0005 | −0.0004 | 2,672 (+101) | −0.0019 | −0.0007 / −0.0007 |

- **D1 HOLDS.** At f 0.35 the curve has an interior peak above ×1.25: ×1.25 +0.0019 < ×1.5 +0.0031 > ×2.0 +0.0005.
- **The peak moves with the regime.** At f 0.5 ×1.25 is best (+0.0011 vs ×1.5's +0.0007). At f 1 nothing beats ×1.
- **Preregistered bin.**
  - s* = argmax at f 0.35 = ×1.5. Its f 0.5 read is +0.0007, so min(f 0.35, f 0.5) = +0.0007 and seeds agree.
  - **WEAK CANDIDATE (+0.0005 … +0.0015)**: record it; the next steps apply only if the held-out read agrees.
- **×2.0 is past the peak and also infeasible on some draws.** Seed 2 stopped in the production emitter with "rc_hybrid_k
  did not meet the per-cell totals" (`work/run655-x200-s2.log`): the design exceeds the K-class emitter's per-cell
  capacity. A production build at α ×2 could hit the same assert.
- **Regions (f 0.35):**
  - Head-DOWN delivered |lfc| rises 0.413 → 0.505 → 0.590 → 0.750 along the dose, and its nmae share is lowest at ×1.5
    (0.0666 → 0.0640 → 0.0633 → 0.0645).
  - Block DOWN improves a little at ×1.5 (0.2208 → 0.2192); uncalled UP keeps improving (0.2385 → 0.2356 at ×2).
  - Past ×1.5 the head-DOWN overshoot returns (×2.0), the same shape 654 saw at ×1.25 without zero repair.

## What it would be on the board (611's call-set rule: nmae × 1.2–1.6, fid × ~2, jac × ~0.5, reach not priced)
| dose | from f 0.35 | from f 0.5 |
|---|---|---|
| ×1.25 | +0.0014 … +0.0018 | +0.0007 … +0.0010 (fid ×2 costs ≈ −0.0002) |
| ×1.5 | +0.0028 … +0.0035 (fid ×2 adds ≈ +0.0007) | +0.0018 … +0.0022 |

The board's regime is f 0.35–0.5 (611), so both doses price positive on A/B/C. ×1.5 prices higher, and leans on fid's call-set
term (+57 calls / target). ×1.25 is the regime-robust choice (max–min over f 0.35 / 0.5: +0.0011 vs +0.0007) with half the
call-set change.

## Recommendation (to the coordinator; the decision and the slot are yours)
1. Carry **×1.25 and ×1.5** to the held-out read (618's fresh-H1 surface, per the coordinator), with zero repair. 618
   chooses between them by the same max–min rule.
2. **Builder KA for the α path**:
   - with `VCC_BULK_LOAD` set, the build at `VCC_ALPHA = 0.982 × s` must differ from today's only in the carriers (and K);
   - the config diff must be exactly `alpha`;
   - on H1-like targets the carriers must match `eff_de × s`.
   `def_postbuild` pins alpha at 0.982 (`scripts/def_postbuild.py:116`), so an October value needs that pin changed too.
3. **An A/B/C board arm** (one dose) with its own prereg: DE + carrier class, priced by 611's call-set rule, behind ALLSW.
4. [d4cd5f]'s UPDN-on-×1.25 stacking test (running on 654's dumps) decides whether the DOWN pass adds on top.
- **Tries for this switch:** 5 amplitudes on H1 (×0.8, ×1, ×1.25, ×1.5, ×2.0); 2 constructions for ×1.25; 0 on the board.
  The best of 3 doses at f 0.35 is mildly selection-biased, which is one more reason to read the held-out surface before
  the board.

## Records
- **KA:** `run655.py cpm x125 0 --zero-repair -5` == 654's x125 zr s0 bit for bit (q, p, lfc, sums).
- **Harness:** `run655_496.py` = 654's copy + `x150` / `x200` branches (diff in `work/`).
- **Files:** `analyze655.py` (654's functions, imported verbatim), `work/analyze655.{out,json}`, `work/seq655.log`. The stats
  npz files stay local.
- **Peak RSS:** 7.2–7.6 GB per run, through gate602.

## α-path KA plan and the pin patch (UNAPPLIED; requested by the coordinator, 10-02 ~21:12Z)
**How α is wired today.**
- `scripts/def_pipeline.py:84`: `SHIPPED_ENV["VCC_ALPHA"] = "0.982"`, the day-of default unless exported.
- `def_pipeline.py:276–284`: the built meta's `config.alpha` must equal the resolved `VCC_ALPHA`.
- `scripts/def_postbuild.py:116`: `SHIPPED_ARMS["alpha"] = (0.982, …)`, asserted by value. Line 452 also requires the
  `--incumbent` reference to carry 0.982.

**A/B/C board arm (needs no code change).**
- Export `VCC_ALPHA=1.2275` (×1.25; 1.473 for ×1.5) on the 446bj2z-style arm env.
- EXPECT config diff exactly `alpha` in `run_arm_generic.sh`.
- If the arm passes `def_postbuild`, name the drift: `--allow-arm-drift VCC_ALPHA`.
- The pin stays 0.982 until October adopts the value.

**Builder KA for the α path** (one serial build on A/B/C, after "builder after restore", gate as usual):
1. **KA-0:** the arm env with `VCC_ALPHA` unset reproduces the incumbent package hash. This is the existing known answer.
2. **KA-α, staging only, no emission:** build the effect stage at α 0.982 and at 1.2275 (251's staging path,
   `effects_array_sha256`).
   - The DE design rows (`context_rows`) satisfy rows(1.2275) = 1.25 × rows(0.982), max |Δ| ≤ 1e-12 relative.
   - The loaded bulk (`bulk_load` sha and the bulk rows) is IDENTICAL. That proves α reaches only the carriers.
3. **KA-α, full build with the rc audit dump on (`VCC_RC_AUDIT_DUMP`):**
   - the config diff vs the incumbent is exactly `{alpha}`;
   - the per-context composition of block totals is unchanged: mean |Δ log1p CP10k| ≤ 1e-3; library totals and ballast K
     may move (654 §2);
   - **554's dose guard PASSES in every context:** median `prior_lambda` ≥ 0.9, `mde_lambda` ≥ 0.9, `call_set_missed`
     ≤ 100. ×2.0 hit the emitter's per-cell capacity on H1, and A/B/C's capacity is ~0.70× H1's, so this check matters most;
   - zero-total cells per context are reported (`zerocount654.py` logic on the package's gated genes).
4. **D/E/F:** 554's dose guard already runs on the day. A failed guard at the adopted α falls back to 0.982 (written rule),
   never to an override.

**The pin patch, for October adoption only.** It applies after the board read and the held-out read agree, in one commit
with the `SHIPPED_ENV` change and a new incumbent reference that carries the new α:
```diff
--- scripts/def_postbuild.py:116
-    "alpha":            (0.982, "VCC_ALPHA", "the a11->a15 amplitude series, closed by 211", ""),
+    "alpha":            (1.2275, "VCC_ALPHA", "654/655: DE carrier amplitude x1.25 under zero repair (alpha now scales only the carriers; the bulk is loaded)", "<board entry>"),
--- scripts/def_pipeline.py:84
-SHIPPED_ENV = {"VCC_ALPHA": "0.982", ...
+SHIPPED_ENV = {"VCC_ALPHA": "1.2275", ...
```
- The "one-line pin" is two lines plus the reference: `--incumbent` must point at a meta built with the new α, or line 452
  reports the reference as stale.
- Also update `def_pipeline.ARM_BOARD_DELTA["VCC_ALPHA"]`'s note.
- `test_make_submission_pm_cd4.py`'s α tests read 0.982 as a parse example only, so they need no change.

## Stacking ([d4cd5f], 10-02 ~21:20Z, on 654's ×1.25-zr dumps, 3 seeds)
- 630's UP pass stays additive on ×1.25: +0.0033 at f 0.35, so ×1.25 + 630 ≈ +0.0052 there.
- 631's DOWN pass (DNMAG) on top of ×1.25 overshoots: increment −0.0005 / −0.0003 / +0.0009 (f 1 / f 0.35 / cmp). Its ×2
  target is relative to the already deeper delivery. **×1.25 and DNMAG are substitutes; one ships, not both.**
- The alternative is ×1 + 630 + DN: +0.0037 at f 0.35, as pure delivery.
- **Note for the pick:**
  - At ×1.25 the call-set term is nominal: fid −0.00002 and jac −0.00006 raw, so 611's provisional fid ×2 / jac ×0.5 barely
    enter its price, and its gain is nmae, priced like a delivery lever.
  - Its real extra risks are carrier capacity (554's dose guard per context; ×2.0 tripped the emitter on H1) and reach,
    which is unpriced for a +26 calls/target change.
  - Pricing both stacks on 618's fresh-H1 surface under the same max–min rule over f 0.35 / 0.5 settles it with one
    instrument.

## Held-out read and the October pick (10-03, [21dd4c]'s 618, 08ebeed1)
- **618 (h1-trainval, 126 unworn targets, faithful key, zero repair):**
  - ×1.25 alone reads +0.0003 … +0.0009 overall-eq (f 0.5 / f 0.35), with nmae −0.0004 … −0.0040 and 1 repaired zero-total
    cell per run.
  - That retains ~0.3–0.5 of H1's +0.0011 … +0.0019. The direction transfers; the size shrinks.
- **The coordinator's pre-fixed rule picks (b) ×1 + 630 + DN:** MIN over {f 0.35, f 0.5} +0.00404 vs (a) ×1.25 + 630
  +0.00360. **×1.25 does not enter October.** It stays recorded as the fallback if DN's builder path is not ready.
- **The bulk term the overall-eq above omits.**
  - Every "overall-eq" in 654/655 is DE-only (lib521.overall_delta: fid, jac, reach, nmae with 192's INV).
  - The carriers' amplitude moves the ballast count K, so per-gene block totals are not pinned. 654 §2 shows
    `sums` ≠ base on every arm; 618 shows the same, ~1.35M cells, up to ~5e4 counts.
  - Priced with 435's H1→board slopes on the harness's raw pds / mse (f 1, paired per seed vs base_zr):

| arm | Δ pds raw | Δ mse raw | bulk term (board units) |
|---|---:|---:|---:|
| ×0.8 | +0.0000 | −0.0006 | +0.0001 (−0.0005 / +0.0004 / +0.0004) |
| ×1.25 zr | +0.0015 | +0.0007 | **+0.0005** (+0.0002 / +0.0008 / +0.0004) |
| ×1.5 zr | +0.0010 | +0.0014 | +0.0001 (+0.0000 / −0.0006 / +0.0010) |
| ×2.0 zr | −0.0031 | +0.0025 | −0.0016 (2 seeds) |

- With the bulk term, ×1.25's H1 total is ≈ +0.0024 at f 0.35, and ×1.5's ≈ +0.0032. No conclusion changes. In
  production the loaded bulk is re-levelled per context, so this harness term is indicative only.
