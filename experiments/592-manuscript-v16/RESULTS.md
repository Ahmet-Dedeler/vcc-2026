# 592 — Manuscript Draft v16 on a copy (for the coordinator's review and republish)

2026-10-01, session "Explore contribution areas [57577c]", at the coordinator's request (scope in its message of
~14:00Z). Method: 533's (edit only the copy in `paper/`; trace every number; revision note on top; pending markers).
Nothing outside this directory was edited. The paper is private until after Nov 5.

## Files and hashes
| file | sha256 (first 16) | note |
|---|---|---|
| source: `experiments/294-writeup/PAPER.md` / `CLAIMS.md` (v15.1, artifact v33) | `9aeefb62…` / `83d3c7e7…` | copied with `cp -c` before any edit |
| `paper/PAPER.md` (v16) | `b0a6a916f82c9f90` | 10 anchored edits + §5.9 |
| `paper/CLAIMS.md` (v16) | `a36dbb460a1a1dec` | 6 anchored edits + C-94 … C-98 |
| `paper/REFERENCES.md` (v16) | `851799e4d8274641` | + CELLOSAURUS, CCLE-2025 |
| `paper/figs/MANIFEST.md` (v16) | `8d63e3dfca9976ea` | F-12's caption: the two 0.0045 noise statements corrected |
| `paper/build_html.py` | `b48336163d9d91d9` | `W` repointed to this copy; the gallery loop now `range(1, 14)` (F-12, F-13 were missing) |
| `paper.html` (built from the copy) | `4d2cb95783c8f769` | 2,850 KB, 31 TOC entries, gallery F-1 … F-13 |
| `apply592.py` | — | the edits, re-runnable on a fresh v15.1 copy (each anchor exactly once) |
| `verify592.py` → `verify592.out` | — | recomputes every new number from primary files and checks the paper states it: **ALL CHECKS PASS** |

To adopt: copy `paper/{PAPER,CLAIMS,REFERENCES}.md` and `paper/figs/MANIFEST.md` over `experiments/294-writeup/`, change
294's `build_html.py` line 27 to `range(1, 14)` (keep its own `W`), build, republish to the same artifact URL.

## 1. What entered (scope as given)
| item | where | claim | status |
|---|---|---|---|
| 460dr (10-01 s1): the day-of LODO key transfers on the October path | §3(h), T12 | C-94 | board-validated (preregistered bin) |
| 446bj2zcp (10-01 s2): the call-propensity block alone NULL | §3(l), T12 | C-95 | board (preregistered NULL bin); H1 positive |
| DE-channel read noise: 0.0045 is a LAW spread, the re-roll is ≈ 0.00025 overall | §3(k) | C-96 | board pairs + H1 + mechanism (537 §2, 547, 548) |
| delivered-dose caveat on C-89 | §3(j) amendment, abstract, plain-language summary, §4.1, §4.7 | C-89 amended, C-97 | build records (pooled) |
| 590: day-of line identification validated | §5.9 (new) | C-98 | preregistered census |
| C-82 corrected | §3(h) | C-82 amended | board (recomputed) |
| pending markers | §3(l) (446bj2zcph, 10-07); revision note (496a 10-05, 523 K pair 10-06; C-88 already carries its marker) | — | — |

## 2. Verification (`verify592.py`, primary files only)
- **Board deltas** from `data/board-entries/<id>.json` (`entry.score_*` and raw members; overall = mean of six; bulk =
  (Δpds + Δmse)/6 scaled, DE = the other four):
  - 460dr − 453dr +0.0028694; mse −0.013060 raw; pds +0.001044 raw. Against the A/B/C key's −0.014473 (446bj2 − 446b) that
    is 0.902×.
  - 460dr − 446bj2: DE −0.0000139, bulk −0.0030321.
  - DE shares against 446b: 446bj2 −0.0003882, 446bd2 −0.0004810, 453dr −0.0005950, 457dr −0.0005473 (457dr − 446b overall
    −0.0005332, bulk +0.0000141).
  - 446bj2zcp − 446bj2z: overall −0.0007942; Δ(fid + jac) −0.0001765 overall; raw jac +0.003986, fid −0.003580, reach
    −0.003353; pds identical; **mse −2.32e-08 raw** (0.6264081539 → 0.6264081307).
- **Doses** from the build metas' `output.rc_callset.prior_lambda` (min / median / max over 900 blocks):
  - 446b 0 / 1.0 / 1; 458 0 / 0.2427 / 1; 458d 0 / 0.4331 / 1; 458g 0 / 0.9526 / 1.
  - Gate sizes from `E468/RESULTS.md`: strong 214 + 217 + 217 = 648, weak 252.
  - 458d's median means ≥ 450 blocks at λ ≤ 0.433. At most 252 of them are weak-target blocks, so **≥ 198 of the 648
    strong-target blocks** are at λ ≤ 0.433.
  - 458g's strong blocks are byte-identical to 458d's (468 compare458g), so the same ≥ 198 holds for 458g.
- **Noise** (C-96) quoted from 537 §2.2–§2.6 (its `noise537b.json` and `redraw537.json` reproduce 230's 0.004501827446 to
  1e-12, per 537), 547 (13 pairs; its correction), 548 (the modern-design rms 0.000256 overall).
- **590** from its own RESULTS / `summary590.json` (this session; known answer exact).

## 3. Discrepancies in the record (for the coordinator; not edited here)
1. **"The direction stands: 458g at λ 0.95"** overstates what 458g shows. It appears at
   `experiments/458-arc-lab-programme/BOARD-RESULT-458d-458g.md:39` (addendum) and `docs/FINDINGS.md:21` (§0, 10-01 bullet).
   - 458g's pooled median (0.953) is lifted by its 252 consensus blocks (λ 1.0 in 446b's build).
   - Its 648 programme blocks are byte-identical to 458d's, and ≥ 198 of them sat at λ ≤ 0.433.
   - Proposed wording: "the harm is established for the programme as this emitter delivers it on these lines; whether a
     fully delivered programme would hurt is open (a dose-only audit rebuild with `VCC_RC_AUDIT_DUMP` would settle it)".
   - October is unaffected either way.
2. **"pds/mse identical"** for 446bj2zcp vs 446bj2z is not exact. It appears at
   `experiments/535-callprop-arm/BOARD-RESULT-446bj2zcp.md:20` ("identical to 16 significant figures") and in its table's
   mse row, and at `docs/FINDINGS.md:27`.
   - pds is identical. mse differs by −2.32e-08 raw (0.6264081539 → 0.6264081307), i.e. a few whole counts moved with the
     block's emission.
   - It makes no difference to any reading, but "identical" is not the word.
3. **The published build omits two figures.** `experiments/294-writeup/build_html.py:27` loops `range(1, 12)`, so the
   gallery stops at F-11. 533's copy had `range(1, 14)`, and v15 added F-12 (October path, single-channel entries) and F-13
   (the hook's H1 price). Artifact versions 32–33 probably lack both; the copy here builds F-1 … F-13.
4. **Attribution.** 554 (the dose guard and its calibration) is the 550–559 peer's work (session "Explore contribution
   areas [4d87df]"), not this session's (57577c, numbers 590–599). The paper cites it by experiment number only.

## 4. What I did not change
- §3(h)'s 453dr bullet keeps "the DE share −0.000595 was reach, the ordering channel" (true); the next bullet now says
  where it comes from.
- T9's table (09-27 … 09-29) is unchanged; the two 10-01 reads are a separate T12.
- C-93's "DE members inside the ordering noise" stands (−0.00039 overall, ~1.5σ of C-96).
- No figure was added. F-12's caption (`figs/MANIFEST.md:73`) did quote 0.0045 twice; both statements are corrected in the
  copy (`apply592.py` `MANIFEST_EDITS`).

## 5. Draft v16.1 (after the 593 audit, 2026-10-01 ~15:00Z)
The cloud audit (`experiments/593-cloud-v16-audit/RESULTS.md`, Grok 4.7 xhigh) confirmed every board number and the
198-block bound. It found 2 WRONG and 6 OVERSTATED statements. Each was re-checked here before writing:
- the H1 diag jsons give 83/85 targets at λ 1.0 (MED23 0.849–0.881, PHF12 0.970–0.987);
- `match590.csv`'s wrong in-CCLE calls are 6, the sixth being schraivogel8 K562 → CORL23, non-unanimous, 0.334;
- correct-unanimous q01 0.681 with 98.9 % ≥ 0.68; all correct q01 0.527 with 97.1 %;
- `noise537b.json`: board_diff_sd_raw 0.001230 = √2 × per_target_sd 0.026098 / √900 (82 targets);
- 392's meta median λ 0.9014.

`apply592b.py` makes the corrections (11 PAPER + 7 CLAIMS anchored edits, a v16.1 note on top). `verify592.py` still
passes. New shas: PAPER.md `ff0b2cd09fce109e`, CLAIMS.md `0a402441d8dd6bfa`, paper.html `e185ec55939f3828`.

**The same "85/85 on H1" statement appears outside the paper** (for the coordinator / the 550–559 peer):
`experiments/458-arc-lab-programme/BOARD-RESULT-458d-458g.md:37` and `experiments/554-dose-guard/RESULTS.md:48`.
Correct: 83/85 (MED23, PHF12 below 1.0).

