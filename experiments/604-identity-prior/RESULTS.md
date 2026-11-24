# 604 — delivery plumbing for the identity-augmented CLEAN4 block prior (council 600; [4d0207])

Science: 639 (d4cd5f) and 622 (7772dc). Price + prereg: d4cd5f. This experiment builds the files and the day-of patch only.
No builder change, no build, no slot. Every heavy step ran through `602/gate602.sh`.

## v3 = CURRENT (2026-10-02 13:35Z): two corrections from d4cd5f's 639d (commit 1f794019) and the 13:30Z ruling
1. **No identity-only extension (v2's fix).** v1 signed genes with no CLEAN4 prior from z(I) alone. 639d measured those
   signs right only 0.27–0.34 on HEK293T's gate genes and 0.49 on HCT116's uncalled ones. The combined score now applies
   ONLY on CLEAN4-covered genes. Uncovered genes keep CLEAN4's pinned psign with no P, exactly as without this module
   (639 block prereg (ii)).
2. **The side is sign(score), not a rank (v3's fix).** v1 and v2 quantile-mapped the SIGNED score onto CLEAN4's P. That
   keeps CLEAN4's down-count among covered genes (53.5 %), so the map undid part of the identity down-shift, and most of
   it in A (Jurkat), where the shift is largest. 639's validated rule (LOLO; 639d routing) is side = sign(zc + λ·zi),
   with zc = (P − ½)/sd centred at P = ½. v3 takes the side from that sign. The confidence |P − ½| is CLEAN4's own
   |P − ½| distribution on covered genes, assigned by the rank of |score|. The rank is stable and breaks ties by CLEAN4's
   own |P − ½|, so λ 0 returns CLEAN4 bit for bit. This is 639 prereg (iii) as reworded by d4cd5f.

| λ | v3 artifact (`work/`) | sha256 | flips vs CLEAN4 A / B / C (block) | down-fraction, covered genes (CLEAN4 0.535) | block down-fraction (CLEAN4 0.677 / 0.678 / 0.659) |
|---|---|---|---|---|---|
| 0 | `ident604v3_flex_lam0_ABC.npz` | `f65d12fc3a7d8a7557e071ea6f86dfc1e89847cd0eea859b6ac14b26ef6fc6e6` (= v2's λ 0: both are CLEAN4) | 0 / 0 / 0 | 0.535 / 0.535 / 0.535 | identical |
| **0.5** | `ident604v3_flex_lam0.5_ABC.npz` | `6cd84dc258cb609d08d93ec4c6c6ec7c7134ff8e3093eccde5ac15bbd2ac618f` | 597 / 601 / 648 (164 / 182 / 187) | **0.582** / 0.539 / 0.542 | **0.739 / 0.757 / 0.727** |
| 1 | `ident604v3_flex_lam1_ABC.npz` | `3b1e22f037aff26b7a13e06328d011757423a51a8cf4acaef7d26cd206fd6a07` | 1,102 / 1,125 / 1,199 (278 / 313 / 330) | 0.617 / 0.553 / 0.540 | 0.777 / 0.812 / 0.775 |

**KAs: `work/ka604v3.log`, 39/39 PASS.**
- λ 0: P == CLEAN4's P bit for bit on covered genes, and psign == CLEAN4 everywhere (also asserted inside
  `identity_prior`).
- For every λ and context, v3's sign == sign(zc + λ·zi) on every covered gene.
- v3 == v2 wherever the two side rules agree. They disagree on 299 / 27 / 44 genes at λ 0.5 and 522 / 111 / 31 at λ 1.
- Uncovered genes == CLEAN4's psign with P NaN.
- The confidence distribution == CLEAN4's sorted |P − ½| (max diff 5.6e-17, from mirroring).

**Superseded, still on disk and unchanged:**
- v1 (`ident604_flex_lam*_ABC.npz`, 71a63aa2 / 2a8a66eb / a9d2c26e) has the extension and the rank side. 612's KGID
  prereg and 622k's review pin it by sha.
- v2 (`ident604v2_flex_lam*_ABC.npz`, 756aec96 / 2fcb840a / f65d12fc) has no extension but keeps the rank side.

Never use v1 or v2 for a board arm. The day-of path (patch604) installs the v3 module: installed def_identity_prior.py
`a2d48da00f01…` (was fe3fbeed), and 606 ALL PASS with it (`606/work/compose606_604v3.log`).
patch604's inserted def_sameline comment still says "quantile-mapped onto CLEAN4's P". v3's exact map is defined in def_identity_prior.py's docstring. The comment was left as is on purpose, so the post-apply def_sameline sha (`d964aa03…` after 644) that the coordinator and [537757] already pinned does not move.

## Flip = BASE + octile sides only (2026-10-02 14:56Z, commit 85d4be68; 7772dc's option b; module `428f4a32…`)
- **Gap 7772dc found after the scale fix.** At λ 0.5 the arms' sides were clean, but v3's rank-mapped confidence
  reshuffled |P − ½| for ~3,800–5,500 NON-octile genes per context (max |Δ| 0.046). A τ > 0 hook could then switch on or
  off outside the octile.
- **Now.** An arm IS the base (no-flip) artifact. Only the masked octiles' SIDES are re-decided, from
  sign(zc' + λ·zi) on the flipped prior with the unflipped scale, and those genes keep the base's |P − ½|.
- **F8.** Every non-octile gene's pdown and psign is byte-identical to the base, at λ 0.5 for octiles 0 / 3 / 7, and the
  octile genes keep their base magnitudes. test604flip **29/29**.
- **The rebuilt 16** (`work/flip625_artifacts.json`): λ 0 is byte-identical to before (each arm = its octile exactly).
  At λ 0.5 the sides changed vs the base are unchanged from the scale-fixed build (206–793; 0 outside), and
  **pdown is identical outside the octile on all 16**. New λ 0.5 shas are in the manifest.
- **Full suites on the pinned module 428f4a32 (15:09Z / 15:18Z, through 602's gate):** test604 **35/35**
  (`work/test604_bb4a.log`; it installs 428f4a32) and test604f **24/24** (`work/test604f_bb4a.log`).
- **Independent check by [7772dc] (~14:57Z, on 85d4be68): clean.**
  - All 16 files match the manifest shas.
  - Outside the octile, pdown and psign are byte-identical to the base in every context.
  - Inside the octile, the pdown side agrees with psign.
  - The kept |P − ½| matches the base to ≤ 5.6e-17 (1 ULP, from rounding 1 − P) on octiles 0–3, and exactly on 4–7.

## Flip SCALE fix (2026-10-02 14:53Z, commit dcef3ca7) and the 16 readout artifacts for 625 (Amendment 3; superseded by the section above for λ 0.5)
- **Fix.** Under a flip, zc's scale sd(P − ½) is now taken from the UNFLIPPED prior. Recomputing it from the flipped
  prior rescaled zc for every covered gene. At λ 0.5 that moved NON-octile genes' sides: octile 0's arm changed 862
  genes against its 794. That would contaminate a module readout.
- **What does not move.** No mask: the same number as before, so the committed v1 / v2 / v3 artifacts are byte-identical
  (F1). λ 0: the side and the rank-mapped confidence are scale-free, so the 8 λ 0 flip artifacts are byte-identical to
  their pre-fix builds.
- **Tests.** F7 asserts that at λ 0.5 a single-octile flip changes sides ONLY inside that octile. test604flip 26/26.
- **Artifacts.** `work/flip625_artifacts.json` lists all 16 (8 single-octile masks × λ 0 / 0.5): paths, shas, the base
  (λ 0: `f65d12fc` = CLEAN4; λ 0.5: ID05 `6cd84dc2`), sides changed vs the base, and changed-outside-the-octile
  (0 everywhere).
  - λ 0 arms change exactly their octile (793–794 genes).
  - λ 0.5 arms change 206–793 genes. The weakest-prior octile 3 changes only ~210, because identity already decides
    most of its genes.

## 604f — the DAY-OF flip pass-through (UNAPPLIED; for a 10-08+ apply AFTER the frozen 10-06 chain)
`patch604f.py` makes two edits.
- **def_sameline:** `--flip-octiles <8 x 0/1>`, passed to `identity_prior(..., flip_octiles=)` and refused without an
  identity programme source.
- **def_day:** `DEF_DAY_DE_PRIOR_FLIP=<mask>`, validated (8 × 0/1) and passed to whichever def_sameline call carries an
  identity prior: the main call when `DEF_DAY_PROGRAMME_SOURCE=identity…`, the hookbank call when
  `DEF_DAY_MDE_SIGN=hookbank:identity…`. A `clean4:τ` hook is refused with a flip, since it would not see it; a flip
  with nothing to carry it is refused; final.env records it.
- It anchors on 644 v3's `${SLSKIP:+--skip-contexts "$SLSKIP"}` call and on 604's identity branch, and requires 630dd and
  644 to be applied.

**test604f 24/24 PASS** (`work/test604f.log`, 14:49Z, peak 3.8 GB, through 602's gate). It covers:
- check / apply / refusals;
- `bash -n`;
- unset and 00000000 leave the identity:0.5 artifact array-identical;
- day-of identity:0 + 10000000 and identity:0.5 + 10010001 on the relabelled bundle equal `build604 --flip`'s A/B/C board
  artifacts;
- def_day's validation in every carrier combination.

606 composes the chain + 604f ALL PASS.

## 607 review fixes (2026-10-02 ~13:55Z; module `7aeb66b4…`, post-apply def_sameline `f147f8bb…`)
- **Name collision.** A pinned reference whose name equals a bundle context's (a rehearsal reusing ids A/B/C with other
  calls) used to overwrite that context in the sister-guard pool. Its own deviation became the pin's, the pin read
  similarity 1.0 and was dropped. Pinned references are now renamed `pin:<name>` on collision. test604flip F6 shows
  the bug on the pre-fix module and its absence now.
- **Score exactly 0.** It now gets P = ½, so psign keeps the gene's own CLEAN4 sign.
- **patch604 dose parsing.** An empty or non-numeric dose is a clean SystemExit; a negative dose is refused. test604 T9.
- The committed v1 / v2 / v3 artifacts are byte-identical: test604flip F1 rebuilds them; the edge cases never occur on
  A/B/C. test604flip 23/23.
- **test604 35/35 PASS** on module 7aeb66b4 (`work/test604flip_full.log`, 14:27Z, peak 6.6 GB through 602's gate):
  the 27 earlier checks plus T9 (4 dose refusals × 2 variants: `identity:`, `identity@`, `identity:x` → "must be a
  number"; `identity@-1` → "must be >= 0"; no traceback).

## 625 octile FLIP option (2026-10-02 ~13:40Z, for [7772dc]'s module readout and the D/E/F decode)
`identity_prior(..., flip_octiles=None)`, CLI `--flip-octiles <8 x 0/1>`, and `build604.py --flip <mask> [--out-dir]`.
- **What it does.** On the covered genes of the masked octiles it applies P_C4 → 1 − P_C4 and negates the pinned psign
  where P ≠ ½. This happens BEFORE zc, so the score, the v3 map, psign and pdown all follow from the flipped prior.
  The λ 0 KA asserts the flipped prior back exactly.
- **The octiles.** They are the `octile` array of `experiments/625-module-readout/work/modules625_v8.npz` (sha256 pinned
  `c24e3304…`; genes matched BY NAME; octile 0 = strongest UP prior … 7 = strongest DOWN). They exactly partition the
  6,346 CLEAN4-covered genes (794 / 793 × 6 / 794), and no covered gene has P = ½.
- **Provenance.** It gains `flip_octiles`, `flipped_genes` and `modules625_sha256` only when the mask has a 1. None and
  "00000000" are bit-identical.
- **Artifacts.** `build604 --flip <mask>` writes `ident604v3_flex_lam<λ>_flip<mask>_ABC.npz`.

`test604flip.py` **21/21 PASS** (`work/test604flip.log`):
- F1: build604 with no mask and with 00000000 rebuilds the COMMITTED v3 artifacts byte for byte (f65d12fc / 6cd84dc2 /
  3b1e22f0).
- F2: λ 0 + 11111111 gives pdown = 1 − P_C4 and psign negated on all covered genes; uncovered genes are untouched.
- F3: λ 0 + a single octile k changes exactly that octile's genes, for all 8.
- F4: λ 0.5 + 10010001 obeys the v3 side rule on the flipped prior.
- F5: invalid masks are refused.

The installed day-of module moves to **def_identity_prior.py `60e984887bf9…`** (was a2d48da0); 606 ALL PASS
(`606/work/compose606_flip.log`). The day-of pass-through for D/E/F is a separate later patch, for a 10-08+ apply:
- `def_sameline.py --flip-octiles <mask>`, with an identity source;
- `DEF_DAY_DE_PRIOR_FLIP` in def_day, passed to the main and the hookbank def_sameline calls and recorded in final.env.

With a flip, the recipe must read its hook from the artifact (`programme:τ` / `file:`), not `clean4:τ`. That needs one
KA on a real build: λ 0 + no flip + `programme:τ` == `clean4:τ` byte for byte, with 612's key gate acting identically.

## The v1 definition (agreed 2026-10-02 12:00Z with d4cd5f and 7772dc; SUPERSEDED by v3 above for the extension and the side rule)
- I_c = log1p(ctrl CPM_c) − mean over references of log1p(ctrl CPM), on the genes c and every reference measure.
  Use only Flex references (622: 3′ references hurt on hek293t): the bundle's other contexts plus the pinned profiles
  A/B/C (challenge controls) and H1 (Arc Flex 2025). Keep one reference per identity call, never one with the context's
  own call (622e: a same-line reference kills the signal; this also de-duplicates a relabelled rehearsal bundle).
- score = z(P_C4 − 0.5) + λ·z(I), where a missing part counts as 0. Each z is taken over the genes where it is defined.
  Positive means down. The ship point is λ 0.5 (622g, Broad-16); λ 1 is the dose partner.
- P is the score QUANTILE-MAPPED onto CLEAN4's own P distribution, so the builder's |P − 0.5| logic (hook τ) behaves as
  it does today:
  - CLEAN4-covered genes get a stable rank-match.
  - Identity-only genes get the same monotone map by interpolation (λ > 0 only).
  - A gene with no score keeps CLEAN4's pinned psign.
- KA (asserted inside `identity_prior`): λ 0 reproduces CLEAN4's psign on every gene and its P on covered genes.

## Files
- `def_identity_prior.py`: the library + CLI. `patch604.py` installs it as `scripts/def_identity_prior.py`, with the refs
  sha pinned.
- `make_refs604.py` → `work/flex_refs604.npz` (+ .json): the pinned Flex reference profiles and their calls.
- `build604.py` → `work/ident604_flex_lam{0,0.5,1}_ABC.npz` (+ .json): the A/B/C board artifacts. They are built by
  the SAME function the day-of path runs.
- `patch604.py` (UNAPPLIED): edits the POST-BUNDLE `scripts/def_sameline.py` (22c6f0f0) so that
  `--programme-source identity:<lam>` gives every context its own identity-augmented prior. Unmatched contexts get w 0;
  matched contexts keep their same-line rows. `def_day.sh` needs no change: `DEF_DAY_PROGRAMME_SOURCE=identity:0.5`
  passes through, and the flags must not contain noprog (def_day.sh already enforces that).
- `test604.py`: T1–T6 on the canonical post-bundle tree (scratch copy; live repo untouched).

## How the board arm uses it (no builder change; the 458d path)
- Identity BLOCK prior: `VCC_SAMELINE=experiments/604-identity-prior/work/ident604_flex_lam0.5_ABC.npz:<sha>:0`
  (programme on, w 0, no rows).
- Identity HOOK as well (21dd4c, source-agnostic 612): `VCC_DE_MDE_SIGN=programme:0.05` (517's form reads pdown_<c>)
  + `VCC_DE_MDE_SIGN_KEYGATE=0.5`.
- Identity HOOK only, block stays CLEAN4: `VCC_DE_MDE_SIGN=file:<same npz>:<sha>:0.05` + `VCC_DE_MDE_SIGN_KEYGATE=0.5`.
- Hard gate (coordinator): per-context median prior λ from the build's audit dump ≥ 0.9 (554's floor), so 458d's 0.43
  does not repeat.

## v1 results (2026-10-02 12:17Z build; `work/build604.log`; SUPERSEDED)
Pinned references `work/flex_refs604.npz` sha `07762ab1…` (A/B/C per-cell CPM float64 + H1 Flex `ctrl_cpm` on its
10,777 genes; calls A JURKAT…, B HELA_CERVIX, C CAL33…, h1 H1_ESC). Peak RSS 4.6 GB.

| λ | artifact sha | context | flips vs CLEAN4 | identity-only genes scored | CPM-block down-fraction (CLEAN4) | hook coverage at τ 0.05 (CLEAN4 2,617) |
|---|---|---|---:|---:|---|---:|
| 0 | `a9d2c26e…` | A / B / C | 0 / 0 / 0 (**KA PASS**) | 0 | 0.677 / 0.678 / 0.659 (identical) | 2,617 |
| **0.5** | `71a63aa2…` | A / B / C | 3,127 / 2,637 / 2,678 | 4,616 | **0.702 / 0.754 / 0.719** (0.677 / 0.678 / 0.659) | 3,414 / 3,432 / 3,529 |
| 1 | `ident604_flex_lam1_ABC.npz` (.json) | A / B / C | 3,833 / 3,409 / 3,419 | 4,616 | 0.716 / 0.795 / 0.767 | see .json |

- Identity genes: 10,777 per context (H1's measured set bounds the common set). CLEAN4 covers 6,346, so the identity
  extends a scored prior to 10,962 genes per context.
- The CPM block's DOWN-fraction RISES at λ 0.5. The emitter delivers a down-skewed prior at full dose, so the 554 gate
  (median prior λ ≥ 0.9) should pass. 458d's 0.43 came from a balanced Arc programme, the opposite direction.
- **Independent verification (7772dc, 622k):** 622's own code on these pinned profiles reproduces the λ 0.5 artifact
  EXACTLY (max |ΔP| 0 on CLEAN4-covered genes; sign agreement 1.00000 and Spearman 1.000000 on all 10,962 genes; flips
  8.8 / 9.4 / 10.2 %). 622's float32 CPM cache differs from these float64 values by 1e-7 … 1e-5 relative; that is
  accumulation error on 622's side, and the KA is P-level equality.
- **Hook coverage trap (21dd4c, kgh612d):** at the builder's τ 0.05, only ~3,4xx of 10,962 genes carry a hook sign.
  Everything priced used full coverage, and with the KEY GATE full coverage beats τ 0.05 (+0.001 … +0.004), so the hook
  must run at **τ 0**: `VCC_DE_MDE_SIGN=file:<ident604_flex_lam0.5_ABC.npz>:<sha>:0` + `VCC_DE_MDE_SIGN_KEYGATE=0.5`.
  At matched full coverage the identity source still adds +0.0048 / +0.0031 / +0.0040 on h1 / HCT116 / HEK293T (CIs
  exclude 0; K562 ns).

## October wiring notes (for the coordinator; nothing here is applied)
1. **Block prior.** Set `DEF_DAY_PROGRAMME_SOURCE=identity:0.5` (after patch604). def_day.sh passes it to
   def_sameline.py unchanged, and def_sameline then writes every context's identity-augmented programme into the same
   artifact as any same-line rows.
2. **Flags.** The resolver passes `with_programme_source = ""` today (docs/october-decisions.json `sameline_flags`,
   487 §1). For the identity source it should pass **`pt`**. A matched D/E/F context then keeps its bank's per-target
   signs (450: gate signs 0.764 vs 0.551), while the identity programme replaces CLEAN4 only as the DEFAULT. `noprog`
   must stay off (def_day.sh refuses it with a programme source).
3. **Hook, at τ 0** (see the coverage trap above). Two October forms, neither needing a def_day.sh edit:
   - **(A) identity HOOK only; the block stays CLEAN4.** This is 21dd4c's arm and is run as a day-of pre-step:
     1. run def_identity.py on the D/E/F bundle;
     2. run `scripts/def_identity_prior.py --bundle <B> --identity <identity.json> --lam 0.5 --out <W>/identity_prior.npz`
        (installed by patch604);
     3. call def_day.sh with `DEF_DAY_MDE_SIGN=file:<W>/identity_prior.npz:<sha>:0`.
   - **(B) identity BLOCK + HOOK:** `DEF_DAY_PROGRAMME_SOURCE=identity:0.5` (patch604's def_sameline path) +
     `DEF_DAY_MDE_SIGN=programme:0` (517's form reads `pdown_<c>` from the VCC_SAMELINE artifact, so no sha is needed).
   - Both need a def_day pass-through for `VCC_DE_MDE_SIGN_KEYGATE` (612), which the coordinator's recipe adds.
4. **Guards already in the function.** One reference per identity call; never the context's own call (622e); the λ 0 KA
   is asserted on every run. The coordinator's hard gate is per-context median prior λ ≥ 0.9 from the build's audit dump
   (554). A day-of call that is not unanimous still yields a call string, so the guard still drops an exact same-call
   reference.

## Exact env strings for the A/B/C board arms (each ONE switch vs its base; prereg by d4cd5f / 21dd4c)
**v3 (13:35Z; CURRENT. The v1 strings below are SUPERSEDED, kept only because 612's KGID prereg pins them):**
- Identity BLOCK prior, λ 0.5 (prior-only, so it shares emission draws), 639's ID05:
  `VCC_SAMELINE=experiments/604-identity-prior/work/ident604v3_flex_lam0.5_ABC.npz:6cd84dc258cb609d08d93ec4c6c6ec7c7134ff8e3093eccde5ac15bbd2ac618f:0`
  ID10 (the dose partner) is the λ 1 file, `ident604v3_flex_lam1_ABC.npz` (`3b1e22f037aff26b7a13e06328d011757423a51a8cf4acaef7d26cd206fd6a07`).
- Identity HOOK only, key-gated, τ 0 (the block stays CLEAN4; needs 612):
  `VCC_DE_MDE_SIGN=file:experiments/604-identity-prior/work/ident604v3_flex_lam0.5_ABC.npz:6cd84dc258cb609d08d93ec4c6c6ec7c7134ff8e3093eccde5ac15bbd2ac618f:0`
  `VCC_DE_MDE_SIGN_KEYGATE=0.5`

v1 (SUPERSEDED: it signs identity-only genes, which 639d found wrong-leaning):
- `VCC_SAMELINE=experiments/604-identity-prior/work/ident604_flex_lam0.5_ABC.npz:71a63aa2525e1ab72258e5864cbfde822c6c3b549cebbad8a077d28c3ee91a11:0`
  (λ 1: `ident604_flex_lam1_ABC.npz`, `2a8a66eb3233618a47b4ce27320c1c7613dbb4f895aa4635e949a9e7f468438c`)
- Gate before any slot: per-context median prior λ ≥ 0.9 in the build's audit dump (554). Rule 10: the config diff vs
  the base is exactly the intended key(s).

## Tests: on v3 (13:15Z) `test604.py` **27/27 PASS** (`work/test604v3.log`; T5 day-of == the v3 board artifact, flips D/E/F 597 / 601 / 648 on the relabelled bundle; installed module a2d48da0). On v1 (12:38Z): 27/27, and `test604_guard.py` 5/5 (logs: `work/test604.log`)
Run on two def_sameline variants: BUNDLE (post-bundle 22c6f0f0) and FULL (+ 641 + 642, the 10-06 state).
- T1: --check, apply, and refusal of a second application.
- T2: the unset path is unchanged (all arrays identical).
- T3: the bank programme-source path is unchanged.
- T4: identity:0 == CLEAN4 exactly.
- T5: the day-of path on the relabelled bundle == the board artifact.
- T6b: a forced HCT116 same-line match keeps its 290 rows / signs / targets / w.
- T6c: the 626 sister guard drops D's own data twin (pinned A, other call), so D's programme == board A.
- T7 / T8: `identity` (bare, λ 0.5) and `identity@1` work, which is what def_day's hookbank form needs.
- The live def_sameline.py and def_day.sh were untouched.

Composition with every other pending patch: 606 (ALL PASS in the coordinator's order; 604 and 642 commute).
