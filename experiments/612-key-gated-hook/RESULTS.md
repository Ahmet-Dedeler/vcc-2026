# 612 — A key-gated sign rule for uncalled genes (two CORRECTIONS below; the second voids the fold evidence)

2026-10-02 ~11:50Z, council session [21dd4c] (experiment 600's breakthrough council). Analysis only: no build, no `vcc`, no
`cell-eval2`, peak < 2 GB, 1.4 s. Script `kgh612.py`, output `work/kgh612.out`, table `work/kgh612.json`.
H1 background: `experiments/610-board-forensics/gatesign610.py`, `dosecurve610.py`.


> **CORRECTION 2 (2026-10-02 12:53Z): the "delivery-faithful" fold evidence is VOID for hook decisions. No held-out
> instrument is validated for this lever. The board decision follows the KGID prereg's A1 FAIL branch: a TIE goes to the
> plain hook.**
> - Cause: kgh612b/c/d, qscan612, and 622's hyb622b, which reproduced them, apply the hook to EVERY gene. Production adds
>   the hook's sign only on UNCALLED genes (`src/vcc2026/rc_emit_k.py`: `use_s = ~head & ~blk & (asg != 0)`), and
>   leaves alone the 400 head genes and the 2,200-gene CPM block (which the block prior steers).
> - ds612c masks that call set. kgID vs no hook, f 1:
>   - h1 +0.0057 → +0.0059 (h1 keeps 55–69 % of its gated |lfc| mass on uncalled genes);
>   - HCT116 +0.0210 → +0.0013, HEK293T +0.0203 → +0.0022 (Orion keeps 20–42 %);
>   - K562 +0.0003 → −0.0021.
> - The plain hook's "K562 −0.010 tail" becomes −0.0013 … −0.0019.
> - Most of the old Orion gain sat on weak-key genes of the CPM block. In production the block prior steers those, not the
>   hook.
> - ds612c also fails the retrodiction that amendment A1 fixed in advance (3efebfd7). On h1 it gives kgID − plain −0.0050
>   / −0.0019 / −0.0027 at f 1 / 0.5 / 0.35. The real emitter gives −0.0032 / −0.0002 / **+0.0007** (3 seeds). The f 0.35
>   sign is wrong, and the shape is not monotone in f.
> - What stands is the H1 real emitter (3 seeds; the only instrument 611 validated):
>   - idkg − c4h nmae: −0.0032 / −0.0002 / +0.0007 / −0.0037 (f 1 / 0.5 / 0.35 / cmp);
>   - reach: −0.012 / +0.004 / +0.004 / +0.001.
>   - So KGID ≈ the plain hook at board-like regimes.
> - The orthogonal knob and its tests (patch612 / patch612dd, 22/22) are still correct engineering. Whether to SHIP the
>   knob is now a board question with a plain-hook default.
> - Also closed (ds612, ds612b): the coordinator's "disagreement shrink". A never-flip clamp loses to kgID on every fold
>   (the hook's value IS its flips), and magnitude shrinks lose −0.004 … −0.037. These were measured in the call-set-blind
>   model, but every variant changes only uncalled-gene signs or magnitudes in the same direction as forms that already
>   lost, so masking cannot rescue them. Not re-run.

> **CORRECTION (2026-10-02 ~12:15Z, after the real emitter): the gain above was over-stated; the gate is a ROBUSTNESS lever,
> not a large expected-value gain.**
> - kgh612 and dosecurve610 modelled the hook as SETTING each uncalled gene's delivered sign from a zero base.
> - The real emitter (rc_emit_k, 496's `mde_add_sign`) only sets the direction of a small MDE add on top of a key-signed base
>   lfc. Measured on seed 0 (diag612): the hook flips only 2–6 % of signs; mostly it SHRINKS |lfc| where CLEAN4 disagrees.
> - On strong-key uncalled genes (sign accuracy only 0.57–0.60) that shrinkage helps the plain hook per target, and the gate
>   removes it.
> - H1, real emitter, seed 0, key-gated vs plain hook (raw):
>   - f 1: nmae +0.0004, reach −0.0002, overall −0.0001;
>   - f 0.35: nmae +0.0011, reach +0.0000, overall −0.0003.
>   - Seeds 1–2 are pending.
> - Delivery-faithful fold check (kgh612b: 480's `base + s·add` shape, base = BETA·consensus), Δnmae, positive = better:
>
>   | surface | plain hook vs no hook | key-gated vs no hook | key-gated vs plain |
>   |---|---|---|---|
>   | h1 | −0.0006 (ns) | −0.0001 (ns) | +0.0006 (ns) |
>   | hct116 | **+0.0149** | **+0.0142** | −0.0007 (ns) |
>   | hek293t | **+0.0148** | **+0.0128** | **−0.0020** [−0.0034, −0.0005] |
>   | k562 (CLEAN4 without K562's banks) | **−0.0102** | +0.0013 (ns) | **+0.0115** [+0.0075, +0.0164] |
>
> - **Reading:**
>   - The plain hook pays where CLEAN4 fits the line (Orion: +0.015) and HURTS where it does not (K562: −0.010).
>   - The key gate keeps most of the gain (+0.013 … +0.014) and removes the downside (K562 ≈ 0).
>   - Worst-case regret: key-gated 0.0020 vs the plain hook's 0.0115.
>   - So for an UNSEEN D/E/F line, the gate is insurance against a CLEAN4 mismatch. Its expected value against the plain
>     hook is ≈ 0 to slightly negative on CLEAN4-friendly lines and strongly positive on CLEAN4-hostile ones.
>   - The earlier "≈ 3× the plain hook" and "+0.016–0.018 over production" were artefacts of the zero-base model.
>     Superseded.

## Question
On H1's uncalled gate genes (610/gatesign610), the transfer key's sign beats CLEAN4's where |key| is high (top 5 %: 0.72 vs
0.51) and loses where |key| is low (bottom half: 0.527 vs 0.572). A HYBRID sign rule (CLEAN4 where |key| is in the
target's bottom half of nonzero |key|, key elsewhere) priced about 3× 496's plain CLEAN4 hook on H1 nmae (610/dosecurve610).
The half-split was chosen on H1 after seeing the table, which is a forking path. Does the rule hold on lines it never saw?

## Instrument
480's cross-study cancer folds, the same gates and code path as `fold480.py`.
- **Surfaces.** Held-out orion_hct116, orion_hek293t and k562_gwps, with h1 as the reference. The Orion truths carry a
  target-independent technical blob, so read RELATIVE comparisons only.
- **Key proxy.** The same-target donor consensus from the OTHER donors (log2 ratios, pc 1 CPM).
- **CLEAN4.** 371's banks minus the held-out line's own: the K562 fold uses Nadig Jurkat + HepG2 only, as in
  `foldsign537.py`.
- **Constructions.** Every one uses 480's shipped MDE-law magnitude (min(0.1/√cnt, 0.1), cnt ≥ 0.3), so only the SIGN
  differs: `mde_key`, `mde_clean4`, `mde_hybrid`.
- **Statistics.** nmae on each surface's gate; 4,000-sample target bootstrap.

## Result: sign accuracy by within-target |key| quantile (gate genes where both sources carry a sign)

| surface | [0, .5) key / CLEAN4 | [.5, .8) | [.8, .95) | [.95, 1] |
|---|---|---|---|---|
| h1 | 0.545 / 0.564 | 0.602 / 0.580 | 0.682 / 0.605 | 0.771 / 0.604 |
| orion_hct116 | 0.536 / **0.704** | 0.599 / 0.669 | 0.683 / 0.678 | **0.861** / 0.562 |
| orion_hek293t | 0.521 / **0.677** | 0.560 / 0.630 | 0.558 / 0.587 | **0.667** / 0.619 |
| k562_gwps | 0.558 / **0.615** | **0.714** / 0.597 | **0.836** / 0.599 | **0.938** / 0.549 |

On every surface, CLEAN4 wins the bottom half and the key wins the top 5 %. The crossover point moves by line: the median
on H1 and K562, near the 95th percentile on the Orion lines.

## Result: nmae, Δ = nmae(reference) − nmae(row); positive = the row is better

| surface | hybrid vs key | hybrid vs CLEAN4 | key vs CLEAN4 |
|---|---|---|---|
| h1 | +0.0009 [−0.0033, +0.0049] | **+0.0116** [+0.0068, +0.0166] | +0.0106 |
| orion_hct116 | **+0.0179** [+0.0141, +0.0218] | +0.0006 [−0.0038, +0.0050] | −0.0173 |
| orion_hek293t | **+0.0163** [+0.0121, +0.0205] | −0.0032 [−0.0068, +0.0003] | −0.0195 |
| k562_gwps | +0.0016 [−0.0024, +0.0056] | **+0.0197** [+0.0135, +0.0264] | +0.0182 |

- **Each single source fails badly on some line.** The key loses 0.017–0.020 to CLEAN4 on the Orion lines. CLEAN4 loses
  0.011–0.020 to the key on H1 and K562.
- **The hybrid's maximum regret against the better single source is 0.0032** (HEK293T vs CLEAN4, CI touching 0). It is
  significantly better than the key on 2 of 4 surfaces and than CLEAN4 on 2 of 4, and never significantly worse than
  either.
- So for an unseen D/E/F line, whose regime we cannot know in advance, the hybrid is the min-max-regret sign rule.

## What it means for production
- **Production.** Its uncalled-gene sign follows the key: on H1's complement, delivered accuracy equals the key's in
  every |key| bin (610/gatesign610).
- **The 10-07 candidate.** 496/537's CLEAN4 hook (τ 0.05) is the CLEAN4 row on the hook's territory.
- **What the hybrid buys:**
  - over production: ≥ 0 on every surface, +0.016 … +0.018 on the CLEAN4-favoured lines, ~+0.001 on the others;
  - over the plain hook: +0.012 … +0.020 on the key-favoured lines, −0.003 … +0.001 on the others.
- **Board expectation** (611 multiplier for nmae delivery levers × 1.2–1.6; the class is a sign rule at fixed delivery
  dose):
  - against production: +0.001 … +0.029 raw nmae, about +0.0003 … +0.008 overall, depending on where A/B/C sit;
  - the dose ×1.5 variant on top (610/dosecurve610 on H1): about +0.004 … +0.005 overall if delivery is rank-neutral.
- **Class: Transfer.** A shrinkage rule that falls back to the gene-level prior where target-specific evidence is weak.
  Nothing in it is fitted to A/B/C. Both inputs (the day-of key, CLEAN4) exist for D/E/F on the day.
- **Caveats.**
  - The key here is 480's donor consensus, not the production whitened key. The production gate should use the
    builder's own |key| rank (521's `load_state()['key']`: 323's cap15 effect key, gated cpm > 10 and b_ctrl > 0).
  - The half split stays a forking-path choice. A crossover fitted per line is NOT allowed: it is per-context tuning.
    Keep the single global split at the median, which the folds just validated as min-max-regret.
  - Orion truths are relative-only.

## Next
d4cd5f folds the hybrid into the rank-neutral UP pass (#1). It prices the hybrid with the identity axis added to the
CLEAN4 half on 435's stats, then the real-emitter run. Preregister the hybrid against the plain hook on the board as one
2-point read. Price it at f 0.35 under 611's rules (nmae × 1.2–1.6; reach read on the board).

## The builder patch (coordinator's spec, 10-02 ~11:40Z: apply on 10-06 after 598, deliver by 10-06 12:00Z)

`patch612.py` is UNAPPLIED. One file (`scripts/make_submission_pm_cd4.py`), three edits; the emitter stays at 0f56209d,
and the `VCC_DE_MDE_SIGN` parse is NOT touched (byte-identical).
- It adds ONE ORTHOGONAL knob, `VCC_DE_MDE_SIGN_KEYGATE=<q>` (0 < q < 1; it needs `VCC_DE_MDE_SIGN`). The knob masks the
  hook's per-target sign, whatever its source, to the target's weak-key genes: |design key| on the gated axis zero or
  below the target's q-quantile of nonzero |key|. The strong half gets 0, so the emitter keeps its delivered (key) sign.
- So the same knob gates:
  - the CLEAN4 hook (`clean4:0.05` + KEYGATE 0.5, the c4kg run);
  - the identity-augmented hook of 622 (`file:` or `programme:` with 604's per-context identity pdown_<c> + KEYGATE 0.5,
    the idkg run). No new source form is needed.
- The meta gains `de_mde_sign_keygate` only when the knob is set.
- It refuses a builder without 598's mark (order 538 → 554 → 590 → 598 → 612) and a second apply.
- (A first draft added a `clean4kg` source form instead. It was replaced before delivery by the orthogonal knob, which
  composes with 604's identity artifact.)

`test612.py` passes **23/23** on temp trees only (T0–T7 for 612; T8–T11 for 612dd below). Since 13:40Z (605's cloud review):
genes TIED at the gate's cut fill the weak count up to floor(q·n) in stable order. It changes nothing without a tie at the
cut: T4a covers 800 synthetic cases (50 with ties), T4b shows the 85 real h1-test keys are identical to run612's gate, and
T4c checks the review's toy keys. T4b SKIPs, without aborting, when the harness file is absent.
- composition 598 → 612 on 598's own canonical post-bundle tree;
- refusal pre-598; second apply refused; compiles;
- the diff confined to the edit blocks;
- `_kg_mask_sign` == run612's gate bit for bit on 600 synthetic cases and on all 85 real h1-test design keys (CLEAN4 kept
  on a median 0.496 of each target's covered genes);
- the `VCC_DE_MDE_SIGN` parse block byte-identical;
- the knob's parse cases (unset → None; 0.5; 0.25; 0 / 1 / 1.5 / non-numeric / without `VCC_DE_MDE_SIGN` refused);
- `_kg_apply` returns the SAME object for every source and the weak-target path when unset, and the masked sign when set;
- the meta key only when set.

**KA plan** (after applying 612 to the live tree, before committing it; three ~40-min builds, serial, 527's runner):

```bash
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
uv run python experiments/612-key-gated-hook/test612.py                         # 23/23 PASS (light; first)
uv run python experiments/612-key-gated-hook/patch612.py --root . --check        # 598 must already be applied
uv run python experiments/612-key-gated-hook/patch612.py --root .
NEWB=$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)
DL=$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz; DLS=00ccfcba1b67fa24ff743ee0c64067e01158147ffe49c5b8ca37cef67600b667
F0="A=1.0962,B=1.1461,C=1.1681"
REFZ=$REPO/submissions/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.h5ad.meta.json
REF4=$REPO/submissions/496a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-as05.h5ad.meta.json
clean(){ unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_RC_AUDIT_DUMP VCC_DE_MDE_SIGN VCC_DE_MDE_SIGN_KEYGATE VCC_RC_CALLSET_CPM VCC_RC_CALLSET_BLOCK VCC_DE_MAG VCC_BULK_CTX_SCALE; }
# KA612a: VCC_DE_MDE_SIGN unset == 446bj2z (ec3a5ae6), as KA598a
( clean; export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F0" VCC_DE_MAG=dn_k250_zs
  KA_Q=1.75 bash experiments/527-arc-hook-arm/run_ka_527.sh KA612a "$REFZ" "$NEWB" )
# KA612b: clean4:0.05 on 446bj2's env == 496a (e21ff752), as 527's KA-496
( clean; export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F0" VCC_DE_MDE_SIGN=clean4:0.05
  KA_Q=1.75 bash experiments/527-arc-hook-arm/run_ka_527.sh KA612b "$REF4" "$NEWB" )
# KA612c: clean4:0.05 + the key gate on the same env; run_ka_527 prints "KA FAIL" (the artifact must differ); ka612.py decides
( clean; export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F0" VCC_DE_MDE_SIGN=clean4:0.05 VCC_DE_MDE_SIGN_KEYGATE=0.5
  KA_Q=1.75 bash experiments/527-arc-hook-arm/run_ka_527.sh KA612c "$REF4" "$NEWB" )
uv run python experiments/612-key-gated-hook/ka612.py experiments/527-arc-hook-arm/work/build_KA612c/ka527-KA612c.h5ad.meta.json
```

**Pass conditions.**
- KA612a and KA612b: "KA PASS (identical)" in `experiments/527-arc-hook-arm/work/build_KA612{a,b}/runKA612{a,b}.log`,
  with driver rc 0.
- KA612c:
  - driver rc 0;
  - `ka612.py` prints KA612c PASS: config diff vs 496a exactly `de_mde_sign_keygate` (+ `de_sign_emitter_sha256` if
    the emitter moved), `de_mde_sign` unchanged ({clean4, 0.05}), q 0.5, and the output differs from e21ff752.
- If KA612b fails because 496a no longer rebuilds on the post-598 builder for a reason outside 612, rebuild KA612b's env
  on the post-598 builder BEFORE applying 612 and use that output as KA612b's reference.

### 612dd: the def_day.sh pass-through (UNAPPLIED; order 598 → 612 → 612dd; ownership agreed with [4d0207])
`patch612dd.py` adds `DEF_DAY_MDE_SIGN_KEYGATE=<q>` to the post-598 `scripts/def_day.sh` as six anchored insertions.
Nothing existing changes; the unset knob changes nothing:
- DOC: the knob is documented next to `DEF_DAY_MDE_SIGN`;
- V: validation needs `DEF_DAY_MDE_SIGN`, 0 < q < 1, and a builder that carries `VCC_DE_MDE_SIGN_KEYGATE` (612);
- E: exports `VCC_DE_MDE_SIGN_KEYGATE` for the level builds;
- C: the final meta must record `config.de_mde_sign_keygate.q` == the requested q. This is an inline check, so
  `riders_match`'s interface is untouched;
- F: `def_day.final.env` records the knob;
- A: postbuild's allowed config diff gains `de_mde_sign_keygate` when the knob is set.

It refuses a def_day.sh without 598's mark, a builder without 612's mark, and a second apply. The identity hook itself
needs no def_day edit: [4d0207]'s patch604 serves it as `DEF_DAY_MDE_SIGN=hookbank:identity:0`.

Tests T8–T11 in `test612.py`:
- the refusals, the apply, the refused second apply, and `bash -n`;
- the insert-only diff (19 lines);
- the validation block run under bash with stubs: unset → 0; 0.5 with a sign → 0; no sign, q 1, q x, or a pre-612
  builder → die (3);
- the final-meta check: match → 0; wrong q or missing key → 1.

The KA plan above does not change: 612dd is plumbing, and KA612c exercises the builder knob it exports.

## The gate quantile q on the folds (qscan612.py; work/qscan612.out)
Hybrid(q) = CLEAN4 sign where the within-target |key| rank < q, key sign elsewhere. nmae per surface; regret = nmae(q) −
the best of {key, CLEAN4, every q} on that surface.

| q | 0 (key) | 0.3 | **0.5** | **0.65** | **0.8** | 0.9 | 0.95 | 1 (CLEAN4 wherever it has a sign) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| max regret over 4 surfaces | 0.0212 | 0.0098 | **0.0033** | **0.0021** | **0.0026** | 0.0042 | 0.0058 | 0.0155 |
| mean regret | 0.0109 | 0.0046 | 0.0018 | 0.0015 | 0.0017 | 0.0025 | 0.0031 | 0.0068 |

- The regret floor is flat between q 0.5 and q 0.8. The block's own q80 protection sits on the same flat floor.
- Inside it, the differences (~0.001) are below the folds' bootstrap widths (0.004–0.008). The default stays the
  pre-specified q 0.5; q 0.65 is the measured minimax, for the record only.
- Picking q from this scan would be another forking path, and it is not worth a board read.

## The identity-augmented low half (idkg), real emitter
On 622's folds ([7772dc], 622/RESULTS.md), replacing CLEAN4 with S = z(P_C4 − 0.5) + 0.5 z(I) in the low half beats this
hybrid on 3 of 4 folds. I uses Flex-only references and fills CLEAN4's gaps with identity alone.
- h1: +0.0048 [+0.0016, +0.0086]
- hct116: +0.0031 [+0.0011, +0.0050]
- hek293t: +0.0040 [+0.0019, +0.0063]
- k562: −0.0013 (ns)

run612's `idkg` mode is that exact delivery on the H1 real emitter (idsign612.py builds S for h1 through hyb622's own
profile and identity functions: 10,777 gated genes signed, 89.5 % agreement with CLEAN4 where both sign). Prices to
follow in work/score612.json.

## Side test: the identity sign on the FILLER block (fillerid612.py; work/fillerid612.out)
Does 622's S (Flex refs, λ 0.5) predict the sign of NON-DE genes on the CPM block, which is fid's filler, better than
CLEAN4? On h1-test (genes both sources sign):

| block K | precision, all-down | CLEAN4 | identity622 | fid, CLEAN4 → identity |
|---|---:|---:|---:|---|
| 1,000 | 0.5183 | 0.5208 | **0.5234** | 0.4412 → 0.4435 |
| 2,200 | 0.5131 | 0.5186 | **0.5201** | 0.4914 → 0.4930 |
| 4,000 | 0.5088 | 0.5162 | **0.5172** | 0.5131 → 0.5143 |

Where the two disagree, identity is right 51–52 % and CLEAN4 48–49 %. The gain is real but small: +0.001 … +0.003
precision, ≈ +0.001 overall at d4cd5f's exchange rate. **It does not explain the pack's +0.014 fid edge**, which stays
unexplained (610 §7). The identity axis earns its value on GATE genes (622: +0.009 … +0.042 sign accuracy), not on the
filler.

## Delivery-faithful comparison of the four hook forms (kgh612c.py; work/kgh612c.out) — the decision table
Shape as kgh612b (480's `base + s·add`, base = BETA·consensus, MDE-law add), on the same folds. Forms:
- none;
- plain CLEAN4 (τ 0.05);
- key-gated CLEAN4 (kgC4);
- key-gated IDENTITY (kgID): 622's S = z(P_C4 − 0.5) + 0.5 z(I), Flex references, via hyb622's own functions; CLEAN4 at
  τ 0.05 where it exists, identity alone where CLEAN4 is undefined.

Δnmae vs no hook, positive = better (target bootstrap 95 %):

| surface | plain | kgC4 | **kgID** | kgID − kgC4 |
|---|---|---|---|---|
| h1 | −0.0006 (ns) | −0.0001 (ns) | **+0.0057** [+0.0015, +0.0099] | +0.0058 [+0.0016, +0.0108] |
| hct116 | +0.0149 | +0.0142 | **+0.0210** [+0.0169, +0.0251] | +0.0068 [+0.0041, +0.0096] |
| hek293t | +0.0148 | +0.0128 | **+0.0203** [+0.0160, +0.0247] | +0.0075 [+0.0046, +0.0104] |
| k562 (CLEAN4 without its own banks) | **−0.0102** | +0.0013 (ns) | +0.0003 (ns) | −0.0010 (ns) |

EV per context (raw nmae gain vs no hook), with friendly = the Orion mean and hostile = the K562 fold:

| P(hostile line) | 0 | 0.10 | 0.20 | 0.33 | 0.50 |
|---|---:|---:|---:|---:|---:|
| plain | +0.0148 | +0.0123 | +0.0098 | +0.0065 | +0.0023 |
| kgC4 | +0.0135 | +0.0123 | +0.0111 | +0.0094 | +0.0074 |
| **kgID** | **+0.0206** | **+0.0186** | **+0.0166** | **+0.0138** | **+0.0105** |

- **kgID dominates the plain hook on every surface**, including the no-tail case: it wins even at P(hostile) = 0.
- kgC4 overtakes the plain hook from P(hostile) ≈ 0.10.
- Sizes:
  - kgID vs no hook: +0.014 … +0.021 raw nmae on the cancer folds ≈ +0.004 … +0.006 overall per context; × 1.2–1.6 on
    the board (611) ≈ +0.005 … +0.009;
  - kgID vs the plain hook: +0.006 … +0.011 raw ≈ +0.0017 … +0.003 overall.
- Caveats:
  - The Orion truths are relative-only (technical blob).
  - On the 3′ surfaces the identity's references are Flex (cross-chemistry); 622 found that matches or beats mixed
    references there. D/E/F are Flex, which is cleaner.
  - The real-emitter H1 check of kgID (3 seeds) is queued, and it decides.
- Recommended ONE form: **key-gated identity hook**, i.e. 604's per-context identity artifact as the hook's source
  (`file:` for hook-only, or `programme:` with VCC_SAMELINE) + `VCC_DE_MDE_SIGN_KEYGATE=0.5`. This holds if the
  real-emitter seeds agree in sign.

## Coverage vs source (kgh612d; work/kgh612d.out): CLEAN4 forms at τ 0.05 and at τ 0 (full coverage, like kgID)

| surface | plain τ0.05 / τ0 | kgC4 τ0.05 / τ0 | kgID | kgID − kgC4 at matched τ 0 |
|---|---|---|---|---|
| h1 | −0.0006 / −0.0032 | −0.0001 / +0.0010 | +0.0057 | **+0.0048** [+0.0016, +0.0086] |
| hct116 | +0.0149 / +0.0180 | +0.0142 / +0.0179 | +0.0210 | **+0.0031** [+0.0011, +0.0050] |
| hek293t | +0.0148 / +0.0173 | +0.0128 / +0.0163 | +0.0203 | **+0.0040** [+0.0019, +0.0063] |
| k562 | −0.0102 / **−0.0139** | +0.0013 / +0.0016 | +0.0003 | −0.0013 (ns) |

- **Coverage:** with the gate, full coverage beats τ 0.05 (kgC4 +0.001 … +0.004). Without the gate, full coverage deepens
  the plain hook's K562 tail. Gate and coverage go together.
- **Source at matched coverage:** identity adds +0.003 … +0.005 on 3 of 4 folds (CIs exclude 0), ns on K562. That is the
  identity's own contribution (622 / 639).
- **The arm therefore uses τ 0**, i.e. `file:<604 artifact>:<sha>:0`. 604's quantile-mapped P passes τ 0.05 on only
  ~3,400–3,500 of 10,962 finite genes per context, which is NOT the form these numbers price.

## H1 real emitter, idkg seed 0: where it gains and loses against the plain hook (diag612b.py; work/diag612b.out)
Per-target nmae contribution of (idkg − c4h), split by uncalled-gene class. Negative = idkg better.

| class | f 1 | f 0.35 | sign accuracy base / c4h / idkg (f 0.35) |
|---|---:|---:|---|
| called | −0.00004 | −0.00009 | |
| uncalled, weak key, CLEAN4-covered | −0.00004 | +0.00000 | 0.475 / 0.566 / 0.566 |
| **uncalled, weak key, NOT CLEAN4-covered** (the identity coverage extension) | **−0.00318** | −0.00011 | 0.507 / 0.507 / **0.541** |
| **uncalled, strong key, CLEAN4-covered** | **+0.00031** | **+0.00108** | 0.657 / 0.646 / 0.657 |
| uncalled, strong key, not covered | +0.00001 | −0.00002 | |
| **total** | **−0.00294** | **+0.00087** | |

- **The identity source is real.** On weak-key uncalled genes outside CLEAN4's τ 0.05 coverage, sign accuracy rises
  0.507 → 0.541 (f 1: 0.520 → 0.541). That is worth −0.003 nmae at f 1. At f 0.35 those genes are a small share of the
  pool.
- **The key gate's cost is MAGNITUDE, not sign.** On strong-key CLEAN4-covered genes the key's sign is the better one
  (0.657 vs c4h's 0.646 at f 0.35). But the plain hook's opposite-direction add SHRINKS |lfc| exactly where two sources
  disagree, which is where the sign is least certain, and per target that helps nmae. The gate gives that shrinkage up.
- **Follow-up idea (not built):** a "disagreement shrink" mode. On strong-key genes where the hook source disagrees with
  the key, shrink |lfc| without ever flipping. It keeps the plain hook's benefit and avoids its K562-style flip tail.
  It needs an emitter change (the emitter is pinned), so it is a post-10-06 candidate only.

## ds612 / ds612b / ds612c: disagreement shrink, 521's regimes on the folds, and the call-set-aware fold model
- **ds612** (`work/ds612.out`). Shape: kgh612c's. Surfaces: 4 folds + the H1 cell-thinning ladder (c100 / c150 / c250 /
  thin).
  - KA: kgh612c's and kgh612d's 12 reference cells reproduce exactly.
  - Never-flip clamps lose to kgID on every surface.
  - Magnitude shrinks ρ·base (ρ 0.5 / 0) where CLEAN4 and/or identity disagree with the key lose −0.004 … −0.037.
  - The exact flip / shrink / grow decomposition puts kgID's whole gain in flips.
- **ds612b** (`work/ds612b.out`): the same forms at 521's pools (the truth's top round(f·n_real) genes by p).
  - The strong-half CLEAN4 shrink's value (IDw_C4s − kgID) at f 1 / 0.5 / 0.35:
    - h1 −0.0006 / −0.0009 / −0.0010;
    - HCT116 +0.0007 / −0.0014 / −0.0017;
    - HEK293T +0.0020 / +0.0013 / +0.0016;
    - K562 −0.0115 / −0.0109 / −0.0088.
  - These were all measured in the call-set-blind model, which CORRECTION 2 voids.
- **ds612c** (`work/ds612c.out`, `ds612c.json`): specified in the prereg's amendment A1 before it existed.
  - The call set is production's construction: head = top 400 by |consensus key|, block = top 2,200 of the rest by
    control expression.
  - Its h1 call set overlaps lib521's production call set by a median 0.876. The heads overlap by only 0.152, because the
    fold's key is the 3-donor consensus, not the design key.
  - KA: an empty call set equals ds612b exactly.
  - The verdict and its numbers are in CORRECTION 2.
- Try count for this family: disagreement shrink, 13 variants × 8 surfaces (ds612); 7 forms × 4 folds × 3 regimes
  (ds612b); 4 forms (ds612c). None is used to select a shipped setting.

## idkg3 and idall3 on H1's real emitter (3 seeds each): two out-of-sample tests of 613d, both PASS
Raw nmae vs c4h (the plain CLEAN4 hook), f 1 / 0.5 / 0.35; overall-eq beside:
| form | Δnmae vs c4h | 613d's registered prediction | overall vs c4h | overall vs base |
|---|---|---|---|---|
| idkg3 (KGID: v3 identity on the weak-key half) | −0.00173 / −0.00021 / +0.00005 | −0.00202 / −0.00034 / +0.00006 | +0.0009 / +0.0006 / −0.0001 | +0.0011 / +0.0009 / +0.0006 |
| **idall3 (v3 identity, ungated)** | **−0.00186 / −0.00218 / −0.00223** | −0.00160 / −0.00232 / −0.00225 | **+0.0014 / +0.0011 / +0.0006** | +0.0015 / +0.0015 / +0.0012 |
- Both predictions land within 0.0003 at every regime (`work/a3verdict612.out`, `work/idall3verdict612.out`).
- The v3 artifact (no identity-only extension) recovers v1's f 1 reach loss: idkg3 − idkg reach +0.0138.
- The key gate is not needed with v3: idall3 − idkg3 nmae −0.0001 / −0.0020 / −0.0023. 613d4 showed why: the gate costs
  on H1 and HCT116 and pays only on the CLEAN4-weakened K562 fold.
- **Outcome:**
  - idall3 replaces KGID (coordinator's rule). Its board prereg is FINAL (`PREREGISTRATION-draft-IDALL3.md`,
    cf3f8004), committed before 10-05.
  - The KEYGATE knob (patch612 / 612dd) still applies 10-06 as unset-by-default engineering. Nothing ships with it set.
- Dumps for [d4cd5f]'s 630 stacks: c4kgd, idkg3d, idall3d (seed 0). Each KA: stats == the scored seed, byte-identical.
