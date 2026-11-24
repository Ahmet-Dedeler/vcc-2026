# 497 — October pass-throughs for the two pending bulk switches (plan and prep, 2026-09-27 06:20Z)

**Why.** 546 prices the October bulk candidate (LODO key with a percentile cap + 479's S_P2) at +14 … +20 % direction on
all 9 October-like folds; both switches are gated on their A/B/C board reads (495 on 10-03, 479's pair on 10-02). The day
script must be able to run them before those reads, so the rehearsals are done when the reads land.

**Done (committed):**
- `scripts/def_lodo_key.py --hi-pct P` (b072b2cc): hi = the P-th percentile of key / median over carried targets.
  KA-pct PASS: data/challenge, rules abc, lo 0.5, **P 96.64697522769474** reproduces 495's cap-6 artifact `dc374281` to
  3.2e-14. The rounded 96.647 (543's text) gives cap 6.000055 and moves 10 weights by 5.5e-5 (543 corrected). The old
  (0.5, 3) KA still passes, and the default-path arrays equal 460's pre-patch KA output exactly.

**Prepared (applied by `chain497.sh` after 476 lands):** `patch_defday497.py` on 476's `def_day.sh` and helpers.
- `DEF_DAY_LODO_HI_PCT=<P>` (needs `DEF_DAY_LODO=1`): `--hi-pct P`, pin id `lodo_lo0.5_pct_def`, the cap logged from
  the key's sidecar.
- `DEF_DAY_SHARPEN_PT=<tercile|stercile>:q,q,q`: the level builds export the per-target form instead of the scalar q;
  `meta_matches` checks `bulk_sharpen_pt` (mode and qs) and requires the scalar q absent; 5b's allowed diff gains
  `bulk_sharpen_q`, `bulk_sharpen_pt`.
- Both unset: same exports, argv and checks as 476's script.
- Checks at prep time: `test_helpers497.py` 9/9 on the patched helpers (476's helpers fail the two per-target cases, as
  they should); 485's `repro_fixes.py` still passes on the 497 copies (`repro/`); invalid knobs stop before setup
  (percentile without LODO, out of range, malformed per-target form).
  - Caution recorded: a first validation test passed two variables through an unquoted zsh variable, which zsh does not
    split. `env` got one malformed assignment, so no validation fired and `def_pipeline`'s setup started on
    data/challenge in a scratch dir. It stopped at once, because another builder (479's KA-T) held the gate. The builder
    stayed in KA-T's pinned state (674fafed = KA-T's recorded `builder_sha256_after`); nothing else was written.

**Rehearsals (`chain497.sh`):**
| check | known answer |
|---|---|
| (a) efp relaxed with the 497 script, neither knob | final output == `ac49d319` (474/476 (c)) |
| (b) 497dr: 460dr's recipe + `DEF_DAY_LODO_HI_PCT=96.64697522769474` + `DEF_DAY_SHARPEN_PT=stercile:2.5,2.0,1.75` | DONE, postbuild PASS, ρ ≤ 1; `gate497.py`: config diff vs 460dr only the dump, the factors and the sharpening form; package + novelty |

497dr is the October-path twin of "446bj2 + cap 6 + S_P2", scorable against 460dr: the October bulk candidate's own
rehearsal, for a slot after the 10-02/10-03 reads.

## 497e (2026-09-28, chain497e 15:48–17:39Z): every candidate October switch at once on efr — PASS
`DEF_DAY_LODO=1 DEF_DAY_LODO_HI_PCT=96.64697522769474 DEF_DAY_SHARPEN_PT=stercile:2.5,2.0,1.75
DEF_DAY_PROGRAMME_SOURCE=h1_2025_full DEF_DAY_PRIOR_GATE=le1 def_day.sh efr … 497efr 1.75 auto pt 4` (def_day e29e793a,
builder 24b95dec = 496p's, so the uncalled-gene hook's emitter pin is in; the hook itself is off).
- rc 0, DONE; postbuild PASS vs 446b with the allowed diffs; zero-gene audit 0 in 0 blocks (lfc < −5: 7, > +3: 10).
- Converged after 3 levelled builds at D 1.3017 / E 1.3685 / F 1.1549; final ρ D 0.9772, E 0.9758, F 0.9800 (all ≤ 1,
  within 0.4 % of ρ* 0.9793). Builder restored.
- LODO key under the percentile cap: hi 7.476 on efr (cap 6 on A/B/C), 3.7 % of carried targets at the cap.
- Per-target sharpening: stercile counts 100/100/100 (edges 3.45 / 4.16), no floor or tie cases.
- Prior gate le1, D: 136 design-strong targets take the programme, 163 weak ones keep CLEAN4, 1 kept on the same-line pt
  prior (threshold 2.0).
- Same-line mode switched on for all three efr contexts (auto w, pt): the mechanics compose with every other switch.
  October's decided same-line setting (w 0.5, noprog / "") is rehearsed by 510, next in the queue.
