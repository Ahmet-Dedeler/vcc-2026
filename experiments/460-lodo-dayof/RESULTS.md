# 460 — the LODO key's October day-of fit (`DEF_DAY_LODO=1`)

Written 2026-09-26 10:30Z, before 446bj2 (446b + the LODO key on A/B/C) is scored on 09-27. It prepares the October path
for a positive read, so no code needs to be written on the night of the read. The flag stays off unless 446bj2 transfers
(433's PREREGISTRATION-0926 bins: overall ≥ +0.001 with pds ≥ 0).

## What was added
- **`scripts/def_lodo_key.py`.** Fits 415's LODO-V1 key on any panel with lib438's arithmetic (full438's lines, h1c1
  composition), from the day's operator and CD4 banks under the October composition rules. It writes the (0.5, 3) weights
  in 410's artifact format (targets, w, key, carried, n_donors).
- **`scripts/def_day.sh` stage 2b, behind `DEF_DAY_LODO=1` (default off).**
  - Fits the key and pins it as `lodo_lo0.5_hi3_def`: a second `def_pin_panel_artifacts.py` call, additive next to the
    410 key.
  - Exports it as `VCC_PERT_SCALE` for the stack-dump build only, as 446bj2 did on A/B/C. The config diff vs the default
    path is `bulk_load` and `bulk_load_sha256`.
  - The restore trap writes the pre-pin builder back, which undoes both pins.

## Known answers
| check | result |
|---|---|
| `--ka`: data/challenge, production operator `e2f65ecf`, rules abc vs 420's `8c82bef3` | **PASS**: targets, carried, n_donors identical; max \|Δkey\| 1.0e-15, max \|Δw\| 2.0e-14 (12 s, 8.0 GB peak RSS) |
| flag off: the dump build's `VCC_*` environment (bash test of the edited lines) | `VCC_BULK_DUMP`, `VCC_H1_UNION_G` only, as before |
| flag on | adds `VCC_PERT_SCALE=lodo_lo0.5_hi3_def` |

**October rules barely move the key.** On the A/B/C panel under rules `def` (CD4 carrying sizing g 0.3443655, pinned
V5_abc deflation), w vs `8c82bef3`: Pearson 0.9999988, max |Δw| 0.0067, Spearman(key) 1.0.
- CD4 fraction 0.3432 vs 0.3461; everything else identical.
- So 446bj2's board read is a read on the October object, not on an A/B/C-only variant.

## Rehearsal (460dr)
`chain460.sh`: `DEF_DAY_LODO=1 def_day.sh` on 180's relabelled bundle (A/B/C controls and panel as D/E/F), q 1.75,
`auto`, `pt:noprog`. That is 453 with the LODO key. It is then relabelled back to A/B/C, prepped, and novelty-checked. It
never submits. Results are below once it ends.

## The 458-positive path, prepared the same way (`def_sameline.py --programme-source`, default off)
If 458 (Arc's H1 programme as every context's block prior, 09-30 slot) reads fid ≥ +0.0015, October switches every
context to Arc's programme. `scripts/def_sameline.py --programme-source h1_2025` writes that artifact on the day:
- every context gets Arc H1's P_down programme (396's rule, the day's panel excluded) with programme_ok true;
- matched contexts keep their same-line rows, signs and w; the rest get 0 rows and w 0;
- the day then runs sameline_flags `pt` (no `noprog`).

| known answer | result |
|---|---|
| KA-1, flag off: def-bundle + 460's identity vs 460's `sameline.npz` (written by the pre-edit code) | **PASS**, byte-identical (`95084748…`) |
| KA-2, flag on: data/challenge, A/B/C forced, `h1_2025` vs 458's `arch1_programme_abc.npz` | **PASS**: programme, pdown, programme_ok and zero-row arrays identical for A/B/C (235 rows, coverage 0.964). w_A is 2 vs 0 (A matches Nadig Jurkat with 0 covered rows, so w is inert) |

The `def_day.sh` pass-through (`DEF_DAY_PROGRAMME_SOURCE`) is added after 460's rehearsal ends: bash reads a running
script incrementally, and that rehearsal is executing `def_day.sh`.
- 10:30Z: `h1_2025_full` (458d's denoised Arc bank) is registered for `--programme-source`. KA-1b: with the flag off the
  output is still byte-identical (`95084748…`). KA-3: `--programme-source h1_2025_full` on A/B/C reproduces 458d's
  artifact's programme, pdown and programme_ok arrays.
- 10:43Z: `461/patch_defday461.py` (the `DEF_DAY_PROGRAMME_SOURCE` pass-through) was tested on the patched stage-2
  fragment with `def_sameline.py` replaced by an echo:
  - unset: the argument list is unchanged;
  - set with `pt`: adds `--programme-source <bank>` and logs it;
  - set with `pt:noprog`: dies, because the builder would ignore the programme;
  - `DEF_DAY_FORCE`: unchanged.
  It is applied after 460's `def_day.sh` exits; the script refuses while one runs.

## Rehearsal result (460dr, 2026-09-26 10:04Z → 11:24Z): PASS
`DEF_DAY_LODO=1 def_day.sh` on 180's relabelled A/B/C bundle, q 1.75, `auto`, `pt:noprog`, max 4 levels:
- Setup and the key: the day's fit (`d9a6b4c4`) is byte-identical to the `data/challenge` rules-`def` fit (a second
  known answer). The second pin added the `lodo_lo0.5_hi3_def` entry (5 lines), and the dump was scaled by it.
- Re-level: lvl0 ρ D 0.8489 / E 0.7793 / F 0.7417 (453: 0.8407 / 0.7703 / 0.7346; the LODO dump carries ~1 % more
  energy). **CONVERGED after 2 levelled builds** at D 1.0804 / E 1.1313 / F 1.1607 (453: 1.0859 / 1.1383 / 1.1662),
  final ρ 0.9795 / 0.9794 / 0.9795.
- Audits: zero-gene 0 in 0 blocks; `vcc prep` dry run and real prep verified 360,000 cells. The builder was restored to
  `4feae27d` (both pins undone by the trap).
- Package `submissions/460dr-october-path-lodo-abc-q1.75.vcc`, sha256 `0a911b57…`, novel.
- `gate460.py`: **ALL GATES PASS**. The config diff vs 453dr is exactly `bulk_ctx_scale` and `bulk_load_sha256`; the
  same-line arrays are identical. `work/gates460.log` is written.
- It is a board candidate only if 446bj2 reads key_ok (PREREGISTRATION.md).
- The `def_day.sh` `DEF_DAY_PROGRAMME_SOURCE` pass-through was applied at 11:26Z, after this run exited. On the live
  file: unset leaves the argument list unchanged, `pt` adds the source, and `noprog` dies.
