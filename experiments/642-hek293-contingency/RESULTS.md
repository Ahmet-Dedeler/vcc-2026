# 642 — the HEK293 contingency: a pinned non-CCLE reference row for day-of identity (UNAPPLIED, default OFF)

Date 2026-10-02, session [537757], at the coordinator's request (approved after the measurement; apply ~10-06 or later,
after 538 + 554 + 590 + 641). Class (AGENTS §11.1c): **day-of instrument**. Nothing is fitted to A/B/C, and it changes no
prediction unless the coordinator switches it on.

## Why

- Arc's 2026-10-01 post says the six lines span "tissue of origin, cell lineage, and disease status", and shows
  immortalized lines (e.g. RPE1, hTERT) as a type (Fig. 3). A D/E/F line may therefore be absent from CCLE.
- Today such a line can never switch on same-line rows:
  - CCLE (1,600 samples) has no HEK293, RPE1, MCF10A, BJ, HaCaT or hTERT lines;
  - 590's floor makes an absent line abstain (real-Flex HEK293: SNU-398, 5/8, 0.647).
- HEK293(T) is the one non-CCLE line with a genome-wide CRISPRi bank on disk: Orion HEK293T, 294/300 of the A/B/C panel.
- RPE1 is not worth it. Its bank is Replogle's essential screen: 0/300 on A/B/C and 8–9/300 on the essential-free
  October-like draws.

## Measurement (`match642.py`; KA PASS)

The matcher is 590's, unchanged: CCLE cache and axis, 4 feature sets × {pearson, spearman}, anchor-weighted Borda.

- **Known answer.** CCLE only reproduces 590's `match590.json` on all **423** common queries (call, unanimity, min_top1;
  0 mismatches).
- **Queries.** 590's profiles plus three real-Flex profiles: Satija VIPerturb pilot K562 / HAP1 / HEK293, Zenodo
  18460279, built by `pilot642.py`. Cells and totals equal 597's.
- **Candidate set.** One extra row, `HEK293_REF`. Variants: xu (XuCao 3′), orion (Orion HEK293T FiCS), flex (the Satija
  Flex HEK293) and nonflex (xu + orion). Positives are cross-source only.

| reference | positives (held out) | false HEK293 calls (of 423) | switch-eligible false calls |
|---|---|---:|---:|
| xu | Orion 8/8 0.889; Flex 8/8 0.783 | 1 (3′ Jurkat set, 7/8 at 0.470: below the floor) | 0 |
| orion | Xu 8/8 0.889; Flex 8/8 0.812 | 0 | 0 |
| **flex** | Orion 8/8 0.809; Xu 8/8 0.785 | **1: Satija's own Flex HAP1, 8/8 at 0.818** | **1** |
| **nonflex** (shipped) | **Flex HEK293 8/8, 0.820** (today: SNU-398 5/8, 0.647) | **0** | **0** |

**Two design rules came out of it.**

1. **Never build a non-CCLE reference from Flex data.** A same-platform reference attracts same-platform queries: the
   Flex-built row took the same study's Flex HAP1 at 8/8 and 0.818. The contexts are Flex, so the reference must not be.
   The production loader refuses a reference whose provenance lists a Flex source, or none.
2. **Freeze the method weights and make it a second stage.**
   - Re-deriving the weights over the augmented set moved 590's weights.
   - Adding a candidate also shifts Borda near-ties. 2–7 unrelated calls moved one-stage, e.g. A549 → A549H3,
     CAKI1 → SLR21.
   - So stage 1 (CCLE only) is computed exactly as today. Stage 2 can only add a HEK293 call, under its own gate:
     augmented consensus = the row, unanimous over 8 methods, min_top1 ≥ **0.75**. That is the stage-2 floor, set by the
     coordinator on 10-02 before any use. Stage 1 keeps 590's 0.65.

**Near misses: unanimity is doing real work** (`test642.json` T3; methods out of 8 with the row as top-1, non-HEK293
profiles):

| profile | methods with the row top-1 | min_top1 | outcome |
|---|---:|---:|---|
| sci-Plex4 A549 | 6 | 0.691 | blocked by unanimity |
| sci-Plex4 MCF7 | 6 | 0.665 | blocked by unanimity |
| Datlinger 3′ Jurkat | 6 | 0.461 | blocked by unanimity and the floor |
| sci-Plex2 A549 | 5 | 0.671 | blocked by unanimity |
| H1 (Arc Flex) | 4 | 0.665 | blocked by unanimity |
| Flex HAP1 (Satija) | 3 | 0.753 | blocked by unanimity |
| HCT116 anchor | 2 | 0.787 | blocked by unanimity |
| **context B (HeLa, Arc Flex, 18,400 cells)** | 2 | 0.719 | blocked by unanimity |

The true HEK293 (Flex, 486 cells) is 8/8 at 0.820. HEK293's embryonic/neural character pulls ESC and some epithelial
profiles part of the way, most of all shallow or odd chemistries. On Arc-Flex profiles the most is 4/8. Read
`row_top1_methods` and the top 5 before deciding.

## The patch (`patch_hek293_642.py`; reference `data/identity_refs/hek293_nonflex_ref642.npz`, sha256 `1b829de0…`)

The reference is built by `build_ref642.py` from XuCao 3′ (2,758 cells) + Orion HEK293T b_ctrl, as a log2-CPM mean on
the official axis. It equals the harness row to 4.8e-7.

| file | change | default |
|---|---|---|
| `def_identity.py` (590) | `NONCCLE_REFS` (pinned path + sha), `NONCCLE_FLOOR` 0.75, `load_nonccle_ref` (sha pin, Flex guard), `nonccle_stage` (stage 1's weights frozen; a gene the row lacks takes CCLE's median); in `match()`, `nonccle_ref` is added only when `DEF_DAY_ALLOW_NONCCLE_REF` names a reference | off: the output is today's |
| `def_sameline.py` (530 + 590) | `REGISTRY_NONCCLE = {HEK293: (orion_hek293t,)}`; `orion_hek293t` in BANK_DIRS, PROGRAMME_OK and STRONG (own-target knockdown visible in 0.911 of rows, equal to orion_hct116's 0.911). A context the second stage called HEK293 never uses its CCLE call, and reads orion_hek293t only with `DEF_DAY_ALLOW_NONCCLE_REF=hek293` AND a non-empty `DEF_DAY_NONCCLE_DECISION` (recorded in the provenance) | no stage-2 block: byte-identical output |
| `def_day.sh` (538 + 554 + 641) | the knob and the decision are logged at start and written to `def_day.final.env` (`%q`, so a decision with spaces round-trips) | logs only when set |

**Day-of rule.** Switch on only when all four hold:

1. `def_identity` with `DEF_DAY_ALLOW_NONCCLE_REF=hek293` calls the context HEK293: 8/8, min_top1 ≥ 0.75.
2. Read `row_top1_methods` and the top 5. A borderline case (stage-1 call unanimous and registered, or the row winning
   only on the consensus) is a stop.
3. The coordinator decides explicitly and writes the reason into `DEF_DAY_NONCCLE_DECISION`.
4. The decision is recorded in the round's RESULTS.

It is never auto-enabled.

**Open question, stated plainly.** 642 shows the line can be RECOGNISED. It does not show that Orion HEK293T's rows
transfer to an Arc-Flex HEK293 context: the 450 same-line mechanism, measured only on H1 (+0.053 … +0.079
overall-equivalent at 35 % coverage). HEK293 (Arc) vs HEK293T (Orion) is also a sister-line pairing.

## Tests (`test_hek293_642.py -v`): 6/6 OK

- First run: 180 s, peak RSS 4.86 GB. That was above the 4 GB announce level; no coordinator window was open, and it is
  reported here.
- Re-run with the 0.75 stage-2 floor: 146 s, peak RSS 3.6 GB.

- **T1, mechanics.** `--check` passes and writes nothing; the patch applies. A second apply, a pre-590/641 tree, a
  partial tree, an anchor miss and a modified reference are all refused, with nothing written.
- **T2, `def_identity --verify` on A/B/C.**
  - Knob OFF: PASS, and the JSON equals the unpatched post-bundle tree's (timing masked).
  - Knob ON: every stage-1 key is equal and `nonccle_ref` is added. A 8/8 JURKAT, B 6/8 HELA (the row is top-1 in 2/8),
    C 8/8 CAL33. None is called HEK293.
- **T3, the coordinator's pins** (through the patched module's `nonccle_stage`):
  - (1) stage 1 reproduces 590 on 423/423;
  - (2) 0 of the 423 non-HEK293 profiles are called HEK293;
  - (3) the held-out real-Flex HEK293 is called HEK293_REF, 8/8, min_top1 0.820, above the 0.75 stage-2 floor.
    In-sample sources: 0.970 / 0.972.
  - Every non-HEK293 profile with the row top-1 in ≥ 5/8 methods sits below 0.75 (0.46–0.69), so it is blocked twice.
- **T4, the Flex guard.** A reference listing a Flex source is refused, as are one with no sources and one with a
  non-pinned sha. The real reference loads.
- **T5, def_sameline.**
  - No stage-2 block: identical arrays and provenance to the unpatched tree.
  - HEK293 on B with the knob unset, no decision, or the wrong knob: B unchanged; A identical.
  - Knob + decision: B reads orion_hek293t, 294/300 covered, w 2, decision recorded; A identical.
  - HEK293 on A: JURKAT's registered call is vetoed.
  - `--force`: identical.
- **T6, def_day.sh lines.** Logged and recorded only when set; the decision round-trips through `source`.

## Reproduce

```bash
uv run python experiments/642-hek293-contingency/pilot642.py <597's rds> <spill dir>   # work/flexpilot642.npz
uv run python experiments/642-hek293-contingency/match642.py                           # the measurement table
uv run python experiments/642-hek293-contingency/build_ref642.py                       # the pinned reference
uv run python experiments/642-hek293-contingency/test_hek293_642.py -v
uv run python experiments/642-hek293-contingency/patch_hek293_642.py --root . --check  # after 538+554+590+641
```

## One sentence

A pinned non-Flex HEK293 row turns today's abstention on real-Flex HEK293 (split, 0.647) into a unanimous 0.820 call,
with no false call on 423 profiles. It stays off unless the coordinator records a decision, and whether Orion's rows then
help is the open part.
