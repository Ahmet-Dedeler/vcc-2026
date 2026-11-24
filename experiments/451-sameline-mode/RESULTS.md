# 451 — the October same-line mode (`VCC_SAMELINE`)

**Status:** code and tests written 2026-09-25 03:00–04:30Z (coordinator). The builder patch is applied by `chain451.sh`
only after two known-answer builds. The results below fill in as the chain runs.

## What it is
When the day's identity call (`def_identity.py`) names a line we hold a public CRISPRi bank for, that context reads the
bank. Evidence is in 450 (H1: a same-line bank from another lab and chemistry carries programme signs, per-target signs
and per-target direction; emitter-level H1 price of the rows at 35 % coverage: +0.08 overall in board units).

- `scripts/def_sameline.py`: identity JSON + registry → `<work>/sameline.npz`.
  - Per matched context: the programme (396's P_down rule, from an unstimulated broad bank only: `PROGRAMME_OK`),
    `programme_ok`, and for the covered panel targets the raw bank rows (unit, non-panel genes) and their signs.
  - Only a unanimous identity call switches a context on (`--force` for known-answer rows).
  - Registry:
    - HCT116 (Orion);
    - K562 (GWPS, VIPerturb, K562-essential);
    - DLD-1;
    - Jurkat and HepG2 (Nadig, essential);
    - A549 (Jiang, Myllia), HT29, MCF7, BXPC3, HAP1 (Jiang);
    - Calu-3 (Sunshine) and THP-1 (Wessels).
- Builder knob `VCC_SAMELINE="<npz>:<sha256>:<w>[:pt][:noprog]"` (`patch_builder451.py`):
  - in a matched context, the block-prior default becomes the programme, unless `noprog` or the context's
    `programme_ok` is false;
  - with w > 0, each covered target's design row and bulk row become unit(ours) + w·unit(bank row), at our row's own
    non-panel norm; panel columns and uncovered targets are untouched;
  - with `pt`, a covered target's block prior is the bank row's own sign where the bank measures the gene;
  - unset, or an unmatched context: the aliases ARE the shared matrices, and the config gains no key.
- October default (450 §9): `:<w>:pt:noprog`, with w from 450d (provisionally 2). The programme is added only if 446bja
  reads positive on the board.

## Known answers
- `def_sameline.py --force A=JURKAT…` on the challenge bundle: programme_A == 396's `psign_A` bit for bit (0 mismatches,
  P_down identical). Artifact `work/sameline_ka_A.npz` (25f91b77…).
- Unit tests (`test_sameline451.py`) on a patched copy, all PASS:
  - unset → None;
  - parse (path forms, w, pt, noprog);
  - five refusals: ja variant, w > 10, short sha, bad flag, repeated flag;
  - loader: sha, gene axis, panel order, programme_ok;
  - blend: w 0 is the identity object; uncovered rows and panel columns are untouched; norms kept; direction is the
    unit sum.
- Builder suite on the patched builder (`1dcb29a3…`): 56 passed, 0 failed, 1 skipped, through the suite's runner
  `experiments/410-pertscale-revisit/run_builder_tests.py`. The first chain attempt called the raw test file, which
  crashes on the one test stale since 481cd092; it restored b6c3f131 as designed. 417's unit tests pass.
- **KA-1 PASS (05:44Z):** unset, the patched builder rebuilds 446b byte for byte (output sha256 `9559923a…`), and the
  config carries no `sameline` key.
- **KA-2 PASS (06:17Z):** `VCC_SAMELINE=<sameline_ka_A.npz 25f91b77>:0` (A = Nadig Jurkat programme, w 0) rebuilds
  446bja byte for byte (output sha256 `10eda1d1…`). The meta records `config.sameline.contexts.A` = {covered 0,
  programme_default true, down_frac 0.7476, w 0}.
- **Committed as builder `1dcb29a3…`** (commit 6bee023a), with both known-answer hashes in the message.

## Mechanics rehearsal (`run451r.sh`)
Forced identity on `ef-bundle-r`:
- D = K562: 300/300 covered;
- E = HCT116: 300/300;
- F = A549: 6/300, rows only, because its banks are stimulated.
It builds at 448's lvl2 factors with `:1:pt` and measures how far the blend moves each context's spread.

**Result (06:44Z):**
- The knob engaged as designed:
  - D and E: programme default, 300 covered targets blended, per-target prior on;
  - F: CLEAN4 default (programme_ok false), 6 blended.
- Zero-gene audit 0 (7 cells below −5). Builder restored to `1dcb29a3` (identical).
- ntclevel ρ at 448's lvl2 factors, without → with the blend:

  | context | covered | ρ without | ρ with the blend (w 1) |
  |---|---:|---:|---:|
  | D (K562) | 300 | 0.9798 | **1.1973** |
  | E (HCT116) | 300 | 0.9798 | **1.0824** |
  | F (A549) | 6 | 0.9463 | 0.9490 |

- **The blend raises the across-perturbation spread**, because bank rows are more target-specific than the shared
  composition. At full coverage it puts a context 8–20 % past its cliff.
- So the day-of re-level must run with `VCC_SAMELINE` already exported, from lvl0 on. `scripts/def_day.sh` does this.
  452b's forced rehearsal checks that the secant loop still converges there.

## Bank quality audit (2026-09-25 20:10Z): how often each registry bank shows its own target knocked down
| bank | rows | own target measured | own < 0 | median own delta | non-zero entries |
|---|---:|---:|---:|---:|---:|
| k562_gwps | 9,787 | 7,617 | **0.990** | −0.657 | 1.00 |
| k562_essential | 2,021 | 1,835 | 0.997 | −0.909 | 1.00 |
| nadig_jurkat | 2,270 | 2,045 | 0.997 | −0.935 | 1.00 |
| nadig_hepg2 | 2,160 | 1,958 | 0.994 | −0.900 | 1.00 |
| sunshine_calu3 | 197 | 174 | 0.989 | −0.641 | 0.60 |
| orion_hct116 | 17,767 | 17,767 | 0.911 | −0.155 | 0.91 |
| wessels_thp1 | 185 | 26 | 0.885 | −0.280 | 0.69 |
| jiang_HT29 | 218 | 194 | 0.825 | −0.743 | 0.62 |
| myllia_a549 | 175 | 175 | 0.800 | −0.164 | 0.81 |
| jiang_BXPC3 / jiang_MCF7 / viperturb_k562 / jiang_A549 | 218 / 218 / 17,489 / 218 | | 0.75 / 0.73 / 0.75 / 0.70 | | |
| jiang_HAP1 | 218 | 194 | 0.619 | −0.405 | 0.61 |
| dld1_lowmoi | 4,715 | 4,715 | 0.144 | 0.000 | 0.14 (a sparse DE-lfc table: most own entries are 0, not positive) |

- The w 2 default was priced on a Pan-quality bank (gate signs 0.76).
- For the 60–80 % banks (the Jiang lines, VIPerturb), the rows are noisier and a lower w may be better. Nothing on disk
  can price that without a truth.
- The registry order already puts the strongest K562 bank first (GWPS 0.99 before VIPerturb 0.75).
- On the day, if the matched bank is one of the weaker ones, prefer w 1 (still +0.084 on H1 with a good bank; the dose
  curve is flat from 1 to 4 at Pan quality).

## Coverage on an October-like panel (`ef-bundle-r`, essential-free, realistic public coverage), 20:40Z
| forced line | bank(s) | panel targets covered | programme |
|---|---|---:|---|
| K562 | k562_gwps | **300 / 300** | yes |
| HCT116 | orion_hct116 | **300 / 300** | yes |
| DLD-1 | dld1_lowmoi | 83 / 300 | yes (sparse table) |
| HepG2 | nadig_hepg2 | 9 / 300 | yes (the main lever) |
| A549 | jiang_A549 + myllia_a549 | 6 / 300 | no (stimulated) |
| MCF7 | jiang_MCF7 | 2 / 300 | no (stimulated) |

- The mode is decisive if a D/E/F line is K562, HCT116 or DLD-1.
- For HepG2 and Jurkat-like lines it is the programme (if 446bja confirms it) plus a handful of rows.
- For the targeted-screen lines it is a handful of rows.
