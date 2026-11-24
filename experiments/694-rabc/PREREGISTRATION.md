# 694 — R-ABC: the October recipe, exactly as resolved, on the real A/B/C bundle (preregistered 2026-10-15 ~05:00Z, session [17adab])

Asked by the coordinator [400b6e] for the 10-16 s1 slot (s2 as fallback). Written before any R-ABC build exists.

## What is built
`october_recipe.py command` (plain: no `--allow-pending`; docs/october-decisions.json is fully decided since 99df0055) run by
`scripts/def_day.sh` on 180's def-bundle (A/B/C relabelled D/E/F), same-panel branch: 453's re-saved operator 66610be7
pre-placed in the run dir, exactly as 524dr / 524dr2 did. The resolved env is expected to be 524dr2's env plus two DE riders:

| export | 524dr2 | R-ABC |
|---|---|---|
| DEF_DAY_LODO / LODO_HI_PCT / SHARPEN_PT / DE_MAG / MDE_SIGN | as resolved 10-06 | same values |
| DEF_DAY_DE_UPRN | off | b_k100_q20r3p10 (630, rank_up 7e6b3e52) |
| DEF_DAY_DE_DNMAG | off | d_m200 (631, dn_mag2 b9d6cfcf = v1 byte-identical, 693) |
| args | `1.75 0.5 noprog 4` | same |

The driver (`rabc694.sh`) stops NOT READY if the resolved env differs from 524dr2's recorded env by anything other than
adding exactly these two exports with these values, or if the args differ.

## Class (AGENTS 11.1c)
Scorer mechanics. Both riders are emission passes that move counts between a perturbation's carrier and ballast rows with
per-gene totals and per-cell libraries exact and 0 call flips; nothing in them is fitted to A/B/C (UPRN and DN were each
adopted from one board arm on 446bj2z's path: 630arm, DN). This read adds no new switch; it checks that the two compose
with the def_day path (tercilec + mde_sign + re-level) as they did on 446bj2z's path. Tries: UPRN 1 board arm, DN 1 board
arm; this is the first read of both on the def_day path. Held-out reading: H1 (630 / 631 RESULTS), which over-read both.

## Prediction
Reference: 524dr2 = 0.2730889 (entry dVHLws6KdMxudq0Yn5ap; pds 0.81077, mse 0.38291, nmae 0.17665, fid 0.00546,
reach 0.23432, jac 0.02844, scaled).
Additive transfer of the two board reads: UPRN +0.0009297 (630arm vs 446bj2z) + DN +0.0003736 (DN vs 630arm)
= **+0.0013033 → 0.2743922**, all of it in nmae (scaled nmae 0.17665 → ~0.18447; raw nmae −0.0047).

Both riders keep pass 1 and every per-(context, target, gene) total, so the re-level reads the same artifact spread and must
land on 524dr2's factors (D=1.0600, E=1.1113, F=1.1260). Rule 13 class: prior-only DE arm on a fixed pass 1, near-exact
(~1e-5 overall); the uncertainty is transfer across bases (UPRN's uncalled-UP carriers and DN's called-DOWN carriers on a
build whose uncalled genes carry mde_sign's CLEAN4 hook), not draw noise.

## Bins on Δ = R-ABC − 524dr2 (overall, board)
| bin | Δ | reading | consequence |
|---|---|---|---|
| SUPERADDITIVE | > +0.0016 | the riders gain more on the def_day path | ships as decided; record the interaction |
| **COMPOSES** | [+0.0010, +0.0016] | additive within the transfer error | ships as decided |
| SUBADDITIVE | [+0.0004, +0.0010) | positive, part of the gain overlaps mde_sign / tercilec | ships as decided; member split into RESULTS |
| NULL / HARM | < +0.0004 | the riders do not carry to the def_day path | flag to the coordinator: UPRN / DNMAG re-opened in october-decisions.json with this read and a member split before Oct 22 |

## Falsifiers (any one voids the composition reading; the cause is diagnosed first)
- F1 bulk: Δpds raw ≠ 0 or |Δmse raw| > 1e-9. Checked BEFORE the slot: `bulkid658.py` on the R-ABC h5ad vs 524dr2's
  package restored must report identical per-(context, target, gene) sums (rc 0), and the meta's `bulk_ctx_scale` must
  equal 524dr2's. Either failing is NOT READY (no slot).
- F2 call flips: |Δfid raw| > 1e-6 or |Δjac raw| > 1e-6 (both riders are 0-flip by construction).
- F3 reach: |Δreach raw| > 1e-5 is reported (UPRN read bit-identical, DN +9.6e-8 raw on 446bj2z's path).

## Pre-slot gates (the run log is arm1015's evidence; ready log = run log = `work/rabc694.log`)
- `[rabc] PASS: config diff vs 453dr …` (gate524's rule with UPRN/DNMAG rows) and `[rabc] PASS: config diff vs 524dr2
  == ['de_dnmag', 'de_uprn']` (exactly two keys differ from the board-scored 524dr2 meta; emitter unchanged, 0f56209d);
- `[rabc] CLIFF OK …`; postbuild all PASS; zero-gene audit 0; meta de_uprn / de_dnmag variants and module shas as above;
- `novelty rc=0`; package sha logged and matching the file; bulk identity vs 524dr2 IDENTICAL;
- terminal line `READY rabc <package sha256>` or `NOT READY <reason>`.
