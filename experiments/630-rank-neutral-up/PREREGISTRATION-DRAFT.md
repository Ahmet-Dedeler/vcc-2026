# 630 — preregistration DRAFT: the rank-neutral UP magnitude pass as a 2-point board dose (VCC_DE_UPRN)

Drafted 2026-10-02 ~12:10Z by council session [d4cd5f], before any entry with this knob exists. **DRAFT:** it becomes a
preregistration when the coordinator fixes the base, the slot day and the arm names, and hashes this file into the
ledger. The bins below are fixed now and are not edited after the first entry.

## Arms (each ONE variable vs the base)
- Base: the best scored entry at slot time, which must carry the zs pass (`VCC_DE_MAG=dn_k250_zs`); expected 479p2g
  or its 10-07 successor.
- **U100** = base + `VCC_DE_UPRN=b_k100`; **U075** = base + `VCC_DE_UPRN=b_k075`.
- Builder: `patch_builder630.py` applied, module `src/vcc2026/rank_up.py` = `630/src_rank_up.py` (sha pinned by the
  patch). KA0 (knob unset rebuilds the base byte for byte) and KA1 (`test_rank_up630.py`: the module reproduces the
  H1 pricing run bit for bit) must pass before the build.
- Expected config diff vs the base: exactly `de_uprn` (variant + module sha). Rule 10 is asserted by `def_postbuild.py`.

## Mechanism and class
After zs, for each uncalled gene the delivered design signs UP:
- one carrier is raised to just above the control maximum;
- it receives the gene's counts from the ballast rows;
- each count is paid by one of that carrier's sub-gate counts moved to the same ballast row, so per-gene totals and
  libraries are exact;
- the gene's net rank-sum change stays within |ΔAUC| ≤ 1e-4;
- a call-flip guard re-runs the block without any gene whose call flipped, so 0 flips.

Only the arithmetic mean of those genes rises, toward t = max(l, K·0.24/√count). **Class: scorer mechanics** (11.1c): it
is computed from each context's own controls and the block, nothing is fitted to A/B/C, K is one global dose, and signs
are today's delivered signs. It does not make the model better at biology.

## H1 prices (real production emitter, 446bj2z's construction; `630/RESULTS.md`, `work/h1/price630.json`)
| arm | H1 f 1 Δnmae raw | overall-equiv at f 1 / 0.5 / 0.35 / cmp | other five members |
|---|---|---|---|
| b_k100 (3 seeds) | −0.0058 | +0.0016 / +0.0027 / **+0.0030 [+0.0014, +0.0044]** / +0.0030 | exactly 0 |
| b_k075 (3 seeds) | −0.0068 | +0.0019 / — / **+0.0025 [+0.0012, +0.0037]** / +0.0028 | exactly 0 |
| (s_k100, reference) | −0.0060 | +0.0016 / +0.0025 / +0.0025 / +0.0029 | exactly 0 |
| (s_k075, reference) | −0.0063 | +0.0017 / — / +0.0022 / +0.0025 | exactly 0 |

The pair spans H1's optimum: at f 1 K 0.75 ≥ K 1.0, at f 0.35 K 1.0 > K 0.75 (the board's regime is f ≈ 0.35–0.5).

**Capacity gate (pre-committed, [4d0207]).** The H1 price rests on H1's delivered fraction of the ask (b_k100: 0.254). The
module writes `ur_frac_got`, `ur_n_moved` and `ur_n_sel` into every block's rc_diag (flat numeric, aggregated per
context like zs's `zs_*`). **Before the slot, each A/B/C context's mean `ur_frac_got` must be ≥ 0.127 (50 % of H1's), and
`ur_n_moved` ≥ 500 per target.** A context below the floor means the lever was never delivered there. The build is
not submitted as this arm; it is re-diagnosed, so a board NULL can be told apart from non-delivery (cf. 554's λ gate
for the prior; 431's context-A capacity history).

**Board prediction (U100 vs base).** Δnmae raw ≈ −0.005 … −0.011. That is H1's f 1 price × the transfer ratios of the two
magnitude levers already read on the board: zs 0.94× (board −0.00218 vs H1 −0.00232) and the coarse MDE 1.76×. It is
also H1's f 0.35 price, −0.0103, at ×1. In overall units: **+0.0014 … +0.0030**. Every other member pinned.
**Amendment (13:34Z, before any entry with this knob exists):** the shipping form is now v3 (`b_k100_q20r3p10`, module
7e6b3e52): H1 f 1 −0.0084, f 0.35 −0.0106 raw. 613e2 puts the cancer folds' uncalled-UP share of nmae's |T| mass at
~0.5–0.8× H1's. Board range: **Δnmae raw −0.004 … −0.015** (f 1 price × 0.94–1.76, × 0.5–1 for the footprint). Center
≈ −0.009 (≈ +0.0025 overall). The center stays inside TRANSFERS, so **the bins are unchanged** (the coordinator's rule:
amend only if the center leaves its band); the low end reaches PARTIAL.

## Noise
The pass consumes no RNG, so every other block of the build, and every non-touched gene, is bit-identical to the base.
The read is exact for this emission realisation, like zs (446bj2z) and the prior-only arms (468). The ~0.0003 re-roll of
bulk steps does not apply.

## Bins (primary statistic: U100's board Δnmae raw vs the base)
Overall equivalents use 192's slope (scaled nmae per raw unit 1.64787, / 6): Δnmae raw −0.0040 ≈ +0.0011 overall;
−0.0015 ≈ +0.0004; +0.0010 ≈ −0.0003.
- **≤ −0.0040 (≥ +0.0011 overall): TRANSFERS.** Adopt for A/B/C arms; October adopts `VCC_DE_UPRN` through a `def_day` pass-through
  (UNAPPLIED patch after the 10-05 bundle, inert KA).
- **−0.0040 … −0.0015 (+0.0004 … +0.0011 overall): PARTIAL.** Adopt only if U075 has the same sign and is not better by > 0.0010. October adopts
  (mechanism-only, nothing fitted), flagged.
- **> −0.0015 (< +0.0004 overall): NULL.** Do not adopt. Record the H1 → board ratio for rank-neutral magnitude levers.
- **Dose rule (both arms scored):** adopt the arm with the more negative Δnmae if both are ≤ −0.0015. If they are within
  ±0.0005 of each other, adopt K 1.0 (H1's choice before any read).
- **Falsifiers (any one: no adoption, diagnose):**
  - |Δpds| > 1e-6 raw (per-gene block totals are exact, so pds cannot move);
  - |Δmse| > 1e-4 raw. mse is NOT structurally invariant: its jackknife runs over cells, and moving counts between cells
    changes the leave-one-out profiles. It measured unchanged on H1, but that is a measurement, not a guarantee
    ([4d0207]'s 605 triage). Report the exact Δmse;
  - |Δfid| or |Δjac| > 1e-6 raw (a call changed; the guard should prevent it);
  - Δnmae > +0.0010 (wrong sign);
  - reach outside 611's band (reach is REPORTED, not binned: on H1 it moved exactly 0, measured, not structural:
    reach walks the reference's DE genes ranked by our lfc, which this pass changes).

## Try count (11.1c rule 4)
- H1: v1 3 variants (up_k100 strict, up_k100_p, up_k125_p) + 2 rejected designs (relaxed donors; R > 1 recipients,
  smoke) + v2 4 variants (s_k100, b_k100, b_k125, b_k150) + K 0.75 (b_k075, s_k075) = **11 emitter configurations**, on one
  surface (h1-test); plus the analytic pricing explorations on the same surface: 603's K × cap grid (15 cells, [4d0207]),
  610's dose curve ([21dd4c]) and 632's graded-magnitude cross-fit (mine), none with the emitter.
- Board: 0.
- The dose and the mode were chosen on H1, so the preregistered pair spans the H1 optimum.

## Closures it must escape (FINDINGS §2)
- **426/435:** MDE doses above c 0.1 leak calls. Here the call set is unchanged by construction (guard) and
  |ΔAUC| ≤ 1e-4.
- **426b:** prior-sign MDE over-fills the classes. No class fill here: the pass acts after emission, repair and zs, and
  flips no sign.
- **434/480:** no class-based magnitude rule beats the MDE; the next points are signs. 480 priced both-direction rules
  under leak-capped delivery. This is UP-only above the leak cap, and DOWN is measured harmful (603: +0.024; 632:
  +0.0035). The larger prize still needs better signs: the S-based hook (604/622) is a separate prior-only switch that
  stacks.
- **Up/down asymmetry:** measured; DOWN is untouched.

## Interaction and October
- **The hook (604/622):** separate switch, prior-only, applied in the emitter before this pass. With better uncalled
  signs, this pass amplifies more right UP genes. They stack, and each is read on its own.
- **D/E/F:** the pass needs nothing but the context's controls (control max, intervals, ref_lin, median library) and the
  block. `def_day` gets `DEF_DAY_DE_UPRN` (pass-through, UNAPPLIED until the coordinator's post-10-05 bundle).
