# 811 PROPOSAL-2 — M811: the many-line consensus on the day's bank list, with and without per-line centring (C × λ = 4)

FROZEN. Written 2026-10-16 08:18–08:25Z (clock; frozen 08:25Z) by the 811 worker, in answer to `VERDICT-1.md` (REVISE). Nobody edits this
file after it is written; a change is a new numbered file. Nothing below has been priced, and nothing will be before a GO.

## 0. Provenance (sha256)
| item | path | sha256 |
|---|---|---|
| the first proposal (judged) | `PROPOSAL.md` | `cc5b8a8fdc2ee07bf790d5e25f4f6b8c3886e0f33155efe1bb303e8e70e2b7f9` |
| verdict | `VERDICT-1.md` | read 08:18Z |
| baseline record | `RESULTS.md` | `f15f0298204a5200f6869f7ad4fd754af1f5c7acbc0cdb49027edf60f50725e9` |
| baseline preregistration | `PREREGISTRATION.md` | `884405075a57ad1bf7d8b57a9dfc0072ca980ddea5f0ef7495265d6bdce0f567` |
| baseline price table | `work/table-base811.md` | `ffeee455bccbf01fe1b6de752f5715198a9169108b7b507bc2e588d8fcd0be43` |
| baseline arrays (list) | `work/run-base811/SHA256SUMS.txt` | `05643d07c5566b3aff9b3e258e32233985525b36b65d7f978362b1e5052ce02b` |
| raw-delta diagnostic | `work/align/align811.json` | `28b88bd22b7c81428be0acf18dd95191beb8f43ad6d04d0ccf167b7db1fc0786` |
| baseline driver (imported by M811) | `base811.py` | `c29f9f976d39f76ab9acc2637b9cb9a66a8db79a786c37c4159d57a2b7a32652` |
| **M811 arm module + fold driver** | **`m811.py`** | **`4818048ee37b3f180b0e7070c95980f2c8ba0481fd2a64e50e5d6b9c8f21de37`** |
| aggregator (stats800, no side effects) | `agg811.py` | `d34d1e2645c063d527019173428c8fcddeaea224334aa29248622f3e5327c9f0` |
| compact banks | `work/banks/manifest.json` | `37e343351000b9b65ceb4adc6a18004563219908cf03b805ed8b27215e57d83b` |

`m811.py` has not been executed. A bug found after GO would mean a new sha, reported before any number is read.

## 1. What changed against PROPOSAL.md (VERDICT-1 A–H)
- **A.** One bank list everywhere: `m811.DAY_GROUPS` = base811's groups minus KOLF2 (the only bank the D/E/F day would not
  have: its full rows are not on disk). The same list on W-tune, W-conf, BULK-ABC and the day path. "With KOLF2" (and with
  TeloHAEC genome-wide) is a later, separately counted try once 813 wires the full rows. Identity rule added (§2.4).
- **B.** Mechanism paragraph rewritten (§2.1): "more lines carry more target-specific direction" is the HYPOTHESIS; the
  per-target raw cosine contradicts it for HCT116/HEK293T; shared-energy ratios corrected; 706's null cited; P1 revised.
- **C.** Factor W dropped (agreement weighting: 016, 142, 271 closed; my own P2 priced it below the readability floor).
- **D–F.** Every prediction is on a surface where the arms it names are priced; per-surface N with e(N); a decision rule.
- **G.** DE-H1: the whole H1 group out; U's banks there listed with their coverage (§3.4).
- **H.** Provenance by sha256 (§0); arm module named and hashed; the A/B/C fold set named (§3.2).

## 2. Mechanism, construction, day path
### 2.1 What the baseline does and does not show
- Engine read (energy-weighted, rendered; RESULTS.md §2.2): at the same λ, the all-lines consensus beats sources-only on
  BULK-W +0.0018 / +0.0050 / +0.0105 (λ 0.10 / 0.25 / 0.50), 17/18 folds; on A/B/C positive with pdsl, literal-pds CIs span 0.
- **Per-target raw cosine (align811; unweighted mean over targets) does NOT rise with more lines for held HCT116 / HEK293T:**
  on A/B/C it falls (0.0381 → 0.0331, 0.0376 → 0.0238), on L1–L6 it is flat (HCT116 0.024–0.032 → 0.025–0.032; HEK293T
  0.024–0.029 → 0.023–0.031). Only held K562 rises (L 0.016–0.018 → 0.023–0.026; A/B/C flat 0.0278 → 0.0278).
- So "more lines carry more target-specific transferable direction" is the hypothesis this proposal tests, not a finding.
  The rival reading: the engine gain comes from strong (energy-heavy) targets and/or from a shared component that the
  aggregate cosine rewards and that pds pays for.
- Across-target shared energy, U_all relative to the truth (align811, unrounded): **3.8–7.9× on A/B/C (KOLF2 included),
  1.3–4.7× on L1–L6** (sources-only: 0.8–2.2× and 0.6–2.9×). On A/B/C without KOLF2 it is unmeasured; the run records it
  (`rec["U"]["c0"]["shared_share"]`).
- **706 found "centring U changes nothing"** on the sources-only U, whose shared energy was 0.2–1.8 %
  (`experiments/706-strong-tier-budget/RESULTS.md` §2 controls). On L folds U_all's shared share is 2.1–2.5 %, only ~1.4×
  706's upper value, so I now expect centring to move little there. **P1 is revised accordingly (§4):** the original "C cuts
  the A/B/C pds cost by ≥ 60 %" cannot be justified: a component holding ≤ 2.5 % of U's energy turns U by ≤ asin(√0.025) ≈ 9°,
  i.e. moves a λ-0.25 blend by ≤ ~0.04 of its norm, while the whole blend moves rows by λ.
- What C still buys is a **mechanism read**: if the many-line gain IS the shared component, per-line centring removes it
  (falsifier F1); if it survives centring, the gain is per-target direction (transfer).

### 2.2 Construction (`m811.py`)
- Groups (`DAY_GROUPS`): K562 {k562_gwps, viperturb_k562, jiang_K562}; HCT116; HEK293T; CD4 {rest, stim8hr, stim48hr};
  H1 {h1_pan2025, h1_vcc2025}; TeloHAEC {pilot}; DLD1; A549 {altair, jiang}; RPE1; Jurkat; BXPC3; HAP1; HT29; MCF7. The held
  line's whole group out on every fold.
- **C off:** U = base811.consensus over DAY_GROUPS (equal weight per line). On L folds this is bit-identical to the
  baseline's all_lX (KOLF2 carries 0 L targets) — a known-answer row (§3.1).
- **C on:** each line's unit rows minus that line's mean over the panel's targets it carries (lines with < 10 carried panel
  targets kept uncentred), panel columns re-zeroed, re-normalised, then the equal-weight mean.
- Blend: base811.blend (707's norm-preserving formula, non-panel genes) at λ; key, sharpening, render, ballast unchanged.
- Arms: `c0_l10`, `c0_l25`, `c1_l10`, `c1_l25`. Class: **transfer** if F1 does not fire; scorer mechanics if it does.

### 2.3 Day path
D/E/F's 300 targets → rows from exactly `DAY_GROUPS`'s banks; centring over D/E/F's own panel; blend before key/sharpen.
Nothing is fitted to A/B/C (no learned weight, λ chosen on L1–L3 only).

### 2.4 Identity rule (D/E/F)
`m811.identity_hold`: if a D/E/F context's identity call (the CCLE match procedure production runs day-of) names a bank
line (K562, HCT116, HEK293T, Jurkat, A549, HT29, MCF7, BXPC3, HAP1, DLD-1, RPE1), **that line's whole group is held out of
U for that context.** Reasons: (1) price fidelity — every priced fold held its own line's group out, so a same-line bank in U
is an unpriced construction; (2) the only same-line probe we have, VIPerturb's own K562 rows, loses −34 … −64 % r against
the composition (RESULTS.md §2.2), so a shallow same-line bank inside an equal-weight mean is a risk, not a gift; (3)
same-line data already reach production through its own same-line route. A context matching no bank line uses every group.

## 3. Protocol
### 3.1 Known answers (any failure voids the fold; the driver exits non-zero)
- `inc`, `707` == 804's arrays (`work/ref804`, SHA256SUMS `31a4ad5b…`) on every L and 062-A/B/C fold, |d| ≤ 1e-12, on 804's
  instance type and thread count (r7i.8xlarge, OPENBLAS 2): the Mac reproduces inc to 8.5e-14 but flips 707's K (RESULTS §2.1).
- `c0_l10`, `c0_l25` == baseline `all_l10`, `all_l25` on every L fold, |d| ≤ 1e-12.
- Engine fold KAs (`E.load(check=True)`) on every fold, ABCp included. A missing reference counts as a failed KA.

### 3.2 Surfaces
- **W-tune** = L1–L3 × {k562_gwps, orion_hct116, orion_hek293t} (9 folds): all 4 arms.
- **W-conf** = L4–L6 × the same (9 folds): the winner only.
- **BULK-ABC (decision set):** 801's **ABCp** folds (A/B/C held-out operators refitted with the day recipe) if they are in
  801's registry at GO time (registry sha256 and each dump's sha256 recorded before the run; inc there = the October
  incumbent pct+P2g, as ABCp is a day-recipe surface; no 804 references exist there, so only the engine KAs). **062's
  A/B/C folds** (fold522.Abc, inc = the board's ALLSW bulk cap6+P2g(abc)) are priced and reported beside it. No dossier on
  062's folds alone (810's G1); if ABCp is not ready at GO, 062's read is reported and the decision waits for ABCp (+1 read
  on the same truths, counted). On BULK-ABC both the winner and its C-partner at the winner's λ are priced (mechanism read,
  counted). ABC/kolf2 reported apart (with KOLF2 never in U it is a clean extra held line).

### 3.3 Selection (W-tune only)
The winner = the arm with the highest 2.5 % bound after the e(10) shift (literal-pds overall vs inc, `agg811.py --tries 10`
on a run dir holding only the 9 W-tune folds). Ties within 0.0002: the smaller λ, then C off.

### 3.4 DE-H1
- U on DE-H1 holds the **whole H1 group out** (h1_pan2025 and h1_vcc2025 = `data/donors/h1_2025_full`). Banks forming U
  there, with coverage of h1-test (85) / h1-trainval (126) targets: k562_gwps 82 / 114, viperturb_k562 80 / 121,
  orion_hct116 84 / 119, orion_hek293t 84 / 120, GWCD4i 80 / 111, dld1_lowmoi 28 / 36, telohaec_pilot 15 / 26, rpe1 21 / 36,
  nadig_jurkat 21 / 35, altair_a549 6 / 9, jiang lines 6 / 14.
- Trigger: the winner's |dK_rel| > 5 % on any fold (the baseline: λ 0.10 ≤ +2.3 %, λ 0.25 up to +5.9 %).
- **Limitation, stated now:** 803's engine builds its own bulk (`--bulk ps`) and has no hook for external composition
  rows. If DE-H1 is triggered, the candidate is "DE unmeasured, K moved > 5 %" and gets no dossier until such a hook exists
  (803's / the lab lead's call, not this lane's code).

### 3.5 Tries (DESIGN §8; e(N) = E[max of N normals])
| surface set | configurations of family `811-consensus` priced there | N | e(N) |
|---|---|---:|---:|
| W-tune (L1–L3) | 6 baseline + 4 M811 | 10 | 1.547 |
| W-conf (L4–L6) | 6 baseline + 1 winner | 7 | 1.364 |
| BULK-ABC (062 and ABCp: the same A/B/C truths) | 6 baseline + winner + C-partner + 1 for align811's read of the truths | 9 | 1.494 |
Every 2.5 % bound named below is after its surface's e(N) shift.

## 4. Predictions (before any number; literal pds | pdsl side by side; vs inc unless named)
**W-tune (all four arms priced):**
- c0_l10, c0_l25: known (KA = baseline on L1–L3); not predictions.
- c1_l10: overall +0.003 … +0.005 | +0.0035 … +0.0055; r +4 … +6 %.
- c1_l25: overall +0.006 … +0.009 | +0.007 … +0.010; r +11 … +15 %.
- **P1 (revised, C − C at the same λ):** |Δoverall| ≤ 0.0010 at λ 0.10 and ≤ 0.0020 at λ 0.25 (both pds variants);
  Δpds(literal) in [−0.001, +0.003]; C on keeps ≥ 70 % of C off's r gain over inc.
- Expected winner: c1_l25 or c0_l25 (W rewards dose).

**W-conf (winner):** if λ 0.25: overall +0.006 … +0.009 | +0.007 … +0.010, r +12 … +16 %, Δpds −0.008 … −0.003, Δmse raw
−0.05 … −0.065. If λ 0.10: +0.003 … +0.005 | +0.0035 … +0.0055.

**BULK-ABC (winner and C-partner):** if λ 0.25: overall −0.002 … +0.001 | −0.001 … +0.002 (baseline all_l25 −0.0009 |
+0.0006; dropping KOLF2 may add up to ~+0.001). If λ 0.10: +0.0003 … +0.0015 | +0.0010 … +0.0025. **Expected bin: NULL
(λ 0.25) or WEAK (λ 0.10); I expect no dossier from this run.**

**Mechanism read R_tgt (reported, not a decision input):** the unweighted mean over targets of the per-target rendered
cosine x/√(a·d) (frame-averaged, from the saved arrays), winner minus the baseline's src_lX at the same λ, per held line:
> 0 on K562 folds, |R_tgt| < 0.002 on HCT116 / HEK293T folds (align811's pattern).

**Falsifiers.**
- **F1 (shared component):** on W-tune, at either λ, c1's r gain over inc < 30 % of c0's r gain over inc (centring removes
  ≥ 70 %). Then the many-line gain is a shared component: class = scorer mechanics on mse, not transfer; no transfer dossier.
- **F2 (hypothesis "more lines → more direction" for October):** the decision set's winner reads NULL or HARM.

## 5. Decision rule (written before the run)
**DOSSIER only if ALL of:**
1. BULK-ABC decision set (ABCp; 062 alone never suffices): winner vs inc, 2.5 % bound > 0 after e(9), literal pds; pdsl
   reported beside it.
2. Paired: the winner ≥ its comparator on BULK-ABC by point (C-on winner: c1 vs c0 at the same λ; C-off winner: c0 vs
   the baseline's all_lX on 062's folds, i.e. dropping KOLF2 does not hurt).
3. W-conf: winner vs inc, 2.5 % bound > 0 after e(7).
4. F1 not fired.
5. Off-axis change share > 25 % on the decision set's mean → the dossier carries the PARTIAL label (702 rule 3).
6. |dK_rel| > 5 % on any fold → DE-H1 clean (each DE member within 2 × its re-roll sd); with no 803 bulk hook, no dossier.
- **Board projection:** none claimed beyond a range that includes 0 (bulk-external has one board point, 702, negative;
  W's ρ 0.19 makes W projections unusable).
**Outcomes:**
- Decision set point ≤ 0 (NULL/HARM): **family closed for October** as a W-panel effect (706's pattern); reopened only by
  the 813 lines (KOLF2 full, TeloHAEC genome-wide) as one counted try.
- Point > 0 but bound ≤ 0 (WEAK): no dossier; same reopening condition.
- F1 fires: reported as scorer mechanics; closed as transfer.
- W-conf point outside [0.7, 1.3] × W-tune point: the selection over-read; W-conf's number is the W price; the rule above
  is unchanged.

## 6. Cost
Stage 1: 9 W-tune folds × (inc, 707, 4 arms), r7i.8xlarge, ~2 min per fold, ~$0.2. Stage 2: winner on 9 W-conf folds plus
winner + partner on BULK-ABC (4 062 folds, up to 4 ABCp folds), ~$0.2. Aggregation ~$0.05. Total ≲ $1 (lane so far $0.45).
