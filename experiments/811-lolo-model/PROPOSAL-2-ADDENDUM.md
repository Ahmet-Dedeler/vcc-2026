# 811 PROPOSAL-2 ADDENDUM — VERDICT-2's binding conditions G1–G5 (frozen before the first M811 launch)

FROZEN. Written 2026-10-16 08:33–08:38Z (clock) by the 811 worker. Amends `PROPOSAL-2.md` (sha256
`8336f260d0cb7fc78e3eaf94f4dca2477b856c700e34c0405ac0a8029e3b50c7`) as `VERDICT-2.md` requires. Where they disagree, this
file wins. Nobody edits it after it is written. Its sha256 is recorded in `work/SHA256SUMS.txt`. No M811 number exists yet.

## G1. ABCp path in `m811.py`: new sha256 `e00192e644828a0d00b18ca5f927f93893ef7012661aac84a4736ca592a98394`
(The previous sha was `4818048e…`; it was never executed.) Changes, and nothing else:
- ABCp is detected as `F.fold_id.startswith("ABCp/") or F.variant == "ABCp"`, never by `F.panel`. engine801's `AbcFold`
  keeps `panel == "ABC"` there (`experiments/801-lab-folds/engine801.py:177-179`).
- `abc` (062-operator A/B/C) is `F.panel == "ABC" and not ABCp`.
- **On ABCp, `inc` = `pct+P2g`** (`E.named(F, "pct+P2g")`: the percentile cap and `p2g_oct`). It is buildable on any fold,
  because it reads only `F.kv` and `F.carried`. No 707 arm is built there, and no 804-reference KA is run; only the engine's
  own fold KAs stand (`engine801.known_answers` ABCp branch).
- "A missing reference = a failed KA" is unchanged on `ABC/*` and `L*/*` folds.
- The ABCp note, `inc_arm` and `variant` are written into `rec` before the meta file is saved.
- The KA `c0_l10` / `c0_l25` == base811 `all_l10` / `all_l25` (every array, |d| ≤ 1e-12) still applies on every L fold. It
  runs on 804's instance type and thread count (r7i.8xlarge, OPENBLAS/OMP/MKL 2 threads). Its outcome is reported in
  RESULTS.md, and a failure voids that fold.

## G2. F1 restated on the many-line increment
- The base811 arrays `src_l10.npz` / `src_l25.npz` for the 9 W-tune folds are byte-copied (`cp -c`) into the stage-1 run
  dir. Each is sha-checked against `work/run-base811/SHA256SUMS.txt` before use.
- The pairing is valid because `inc` is identical by KA: base811's and m811's `inc` both equal 804's arrays.
- **F1 fires** if, at either λ, c1's r gain over src_lX is < 30 % of c0's r gain over src_lX. Here r gain over src_lX =
  r(arm vs src_lX), the relative aggregate cosine from `agg811.py --pairs cX_lY:src_lY` on W-tune.
- If F1 fires, the class is scorer mechanics (a shared component) and there is no transfer dossier. Otherwise the class
  is transfer, provided the held-out gates pass.
- **P1's inc-based ratio** (c1's r gain over inc ÷ c0's r gain over inc) is reported beside it, labelled as P1's read. For
  reference, full removal of the many-line increment would give 3.89 / 5.86 = 66 % (λ 0.10) and 8.67 / 14.02 = 62 %
  (λ 0.25). Those are below P1's "≥ 70 %", so a P1 miss in that range is what the rival reading predicts.
- Caveat: centring also touches the source lines inside U. 706 found that centring the sources-only mean changes nothing,
  so a c1 − c0 difference is read as acting mainly on the non-source lines. That attribution is an assumption, not a
  measurement.

## G3. Selection eligibility (applied before selection)
- An arm whose maximum over the 9 W-tune folds of |dK_rel| (engine `vs_base.dK_rel` vs inc) exceeds 0.05 is ineligible as
  the winner. Rule 6 bars it from a dossier, because 803 has no external-rows hook. It is still reported on W-tune.
- The winner is the eligible arm with the highest 2.5 % bound after e(10), with PROPOSAL-2 §3.3's tie rule.
- If no arm is eligible, I report "no DE-eligible arm". Stage 2 then prices only the λ 0.10 C-pair (`c0_l10`, `c1_l10`) on
  W-conf and BULK-ABC as a mechanism read, within the same N. There is no dossier in that branch.
- Any closure written in RESULTS.md names the λ and C actually read on the decision set and closes no other configuration.

## G4. Counts and a corrected fact
- **N on BULK-ABC:** 9 (e(9) 1.494) when ABCp is the first A/B/C read of these arms. **N = 10 (e(10) 1.547) in the branch
  where ABCp is read after a 062 read of the same arms**, which is the branch this run takes if ABCp is not registered
  when stage 2 is due.
- N on W-tune is 10 (e 1.547) and on W-conf 7 (e 1.364), as in PROPOSAL-2. If the winner is a c0 arm, its W-conf read is
  the baseline's already-priced all_lX (bit-identical on L), adds no new information, and is reported that way.
- **Correction to PROPOSAL-2 §2.1 (":45-46").** The shared share of U without KOLF2 on A/B/C is NOT unmeasured.
  `work/align/align811.json` (sha256 `28b88bd2…`) has `ABC/held/{K562,HCT116,HEK293T}/all-KOLF2/shared_U` = 0.0203 / 0.0219 /
  0.0207, which is 1.8× / 4.0× / 2.0× the truth's. So C's expected effect on A/B/C is the same small size as on L, and that
  supports the revised P1.

## G5. Wear
This run's would-be wear rows are appended to `work/wear811.jsonl` (family `811-consensus`), each with its surface set's N:
W-tune 10, W-conf 7, BULK-ABC 9 (or 10 in the G4 branch). The lab lead adds them to 800's ledger. `check800.py aggregate`
stays unused.

## Stage plan (as PROPOSAL-2 §3, with the above)
1. Stage 1: `m811.py <fold> --arms c0_l10,c0_l25,c1_l10,c1_l25` on L1–L3 × 3 held lines (EC2, r7i.8xlarge, 2 threads).
   Aggregate vs inc (`--tries 10`) and vs src_lX (G2). Then apply eligibility (G3) and select.
2. Stage 2: the winner on L4–L6 × 3. The winner and its C-partner at the same λ on `ABC/*` (062). On `ABCp/*` if 801's ABCp
   dumps are registered and staged, otherwise later, under the N = 10 branch.
3. The decision rule (PROPOSAL-2 §5) is applied on ABCp only. Without ABCp, W-conf and 062-ABC are reported and the
   decision waits.
