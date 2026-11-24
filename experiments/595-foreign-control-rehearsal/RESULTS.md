# 595 — The day-of control stages on FOREIGN control cells: every registry line switches same-line ON correctly, the two non-CCLE lines stay off (HEK293T only because of 590's floor), and the only stop is the designed one (554's spread guard on a bimodal Calu-3)

2026-10-01, session "Explore contribution areas [57577c]", at the coordinator's request. Class: day-of procedure.
- **Scratch.** A worktree of HEAD `99e32a61` with the 10-05 patches applied (538 → 554 → 590; the 11 shas assert-equal
  594's final table). The stage commands are exactly the driver's / `def_day.sh`'s. No builds; `def_day.sh` never
  executed.
- **Live files untouched** (hashed before and after).

## The bundles (`build595.py`; seeded; written to this session's scratchpad, 3.8 GB, not committed)
- Every rehearsal so far used A/B/C's own controls. These two bundles hold real control cells of six other lines in the
  exact bundle format: the official 18,533-gene axis (a gene the source lacks is a zero column), 18,400 cells per
  context in 46 `ntc_id`s × 400, the CURRENT pert panel and a D/E/F manifest.
- Orion counts are integer-valued, checked.
- Bundle Y's sources hold 2,386–4,156 control cells and are resampled to 18,400 with replacement. That preserves the
  pseudobulk and the library-size quantiles, which are all these stages read.

| bundle / ctx | line (source, chemistry) | cells in pool | genes with counts | median library | **q75/q25** (554 floor 3.0) |
|---|---|---:|---:|---:|---:|
| X / D | HCT116 (X-Atlas/Orion NTC, FiCS) — registry, strong bank | 18,681 | 16,771 | 16,803 | 1.771 |
| X / E | H1 (Arc VCC 2025 controls, **10x Flex**) — not in CCLE | 18,400 | 17,117 | 20,186 | 1.575 |
| X / F | HEK293T (Orion NTC) — not in CCLE (590's false-HAP1 case) | 18,630 | 17,181 | 16,113 | 1.793 |
| Y / D | A549 (Myllia NTC + gene desert, 10x 5′ HT) — registry, weak banks | 4,156 | 15,749 | 12,220 | 1.572 |
| Y / E | THP-1 (Papalexi ECCITE controls, 10x 5′) — registry, weak bank | 2,386 | 13,099 | 8,856 | 2.100 |
| Y / F | Calu-3 (Sunshine NTC, 10x 3′, SARS-CoV-2 screen) — registry, strong bank | 3,536 | 14,661 | 11,230 | **7.532** |

## The questions
**(i) Does an unforced identity call switch same-line ON for a real registry line?** Yes, for all four, each unanimous and
above 590's floor:

| ctx | call | unanimous | **min_top1** | consensus top 5 | same-line |
|---|---|---|---:|---|---|
| X/D | HCT116_LARGE_INTESTINE | yes | **0.794** | HCT116, HCT15, SNU638, SNU-C5, SW48 | **ON**: `orion_hct116`, w 2.0, **290/300** targets |
| Y/D | A549_LUNG | yes | **0.773** | A549, A549H3, NCI-H2023, NCI-H1355, NCI-H2030 | **ON**: `jiang_A549` + `myllia_a549`, w 0.5, 22/300 |
| Y/F | CALU3_LUNG | yes | **0.765** | Calu-3, CFPAC-1, OE33, PANC03.27, SNU-213 | **ON**: `sunshine_calu3`, w 2.0, 6/300 |
| Y/E | THP1_HAEMATOPOIETIC… | yes | **0.665** | THP-1, MONO-MAC-1, MONO-MAC-6, OCI-AML3, OCI-AML2 | **ON**: `wessels_thp1`, w 0.5, 1/300 |

**(ii) Do the non-CCLE lines stay off?** Yes, by two different routes:

| ctx | call | unanimous | **min_top1** | consensus top 5 | same-line |
|---|---|---|---:|---|---|
| X/E (H1) | TERA2_TESTIS | no (7/8-type split) | 0.653 | TERA-2, NCCIT, NTERA-2 cl.D1, TERA-1, GCT27 | off: no registered bank |
| X/F (HEK293T) | **HAP1_ENGINEERED** | **yes** | **0.616** | HAP1, LS, SNU-398, COLO-320, DMS-273 | **off only because of 590's floor**: "registered and unanimous, but min top-1 correlation 0.616 < 0.65; context unchanged (590)" |

Without 590's patch the HEK293T context would have switched HAP1's rows on. This is the in-situ confirmation, in a
full-axis 18,400-cell context file, of what 590 found on pseudobulks.

**(iii) Library spreads against 554's floor (3.0).**
- Five foreign lines in four chemistries read 1.57–2.10. H1 (Arc's Flex) reads 1.575, inside A/B/C's 1.57–1.95.
- Calu-3 reads **7.53** (bimodal libraries, the infected-state screen R7b also hit). `def_preflight.py` FAILs it: "[F]
  control library q75/q25 <= 3.0 … 7.532 (q25 2704, q75 20367)", rc 1.
- With `DEF_DAY_ALLOW_WIDE_LIBS=1` the same run passes (rc 0) with the ALLOWED row. The override path behaves.

**(iv) Census.** rc 0 on both bundles with the panel's own CD4 extraction (199's A/B/C KA). Targets no production donor
carries: **0** (floor 5). The census reads only the panel, so on the current panel it cannot differ from A/B/C's. Its
override (`--max-zero-carrier`) was not exercised, because a count of 0 cannot be made to fail.

**(v) What Oct 22's log will look like.**
- A correct registry call: unanimous, min_top1 0.67–0.79 here (0.72–0.77 for A/B/C in Flex), and a
  "590 gate passed" line.
- A line absent from CCLE: either a split call naming a biological neighbour (H1 → the testicular embryonal-carcinoma
  lines, 0.653), or a unanimous call at a low score that the floor blocks (HEK293T → HAP1, 0.616).

## Anything that would stop or mislead the day
1. **A Calu-3-like library spread stops the day, by design.** Overriding it only lets the day continue: at such spreads
   the K-class emitter's steering collapses for that context (553), and 554's postbuild dose guard will then flag the
   context too. Arc's own Flex lines (A/B/C, H1) sit at 1.57–1.95, so a D/E/F trip would be a real anomaly, not a
   chemistry artefact.
2. **THP-1's correct call sits at 0.665, 0.015 above the floor.** It is a 10x 5′ profile with ~5,400 official genes
   absent (zero columns). A degraded registry-line profile can come close to the floor. The coordinator's day-of rule
   (log min_top1; read the top 5 if it is 0.60–0.70) covers it: THP-1's top 5 are all monocytic AML lines. The banks
   such a call would switch on are small (1/300).
3. **Nothing misleading found.** Every call names the true line or a biological neighbour, and same-line ON/OFF is right
   in all six contexts.
- Census coverage cannot fail on the current panel. On a new panel it is the panel, not the controls, that decides.

## Deviation from the agreed memory budget (stated, as AGENTS §9 requires)
The rehearsal itself peaked at 3.1 GB (whole process tree, sampled). **The bundle BUILD did not stay light:**
- attempt 1 reached 15.9 GB max RSS: pyarrow's default pool kept the Orion parquet list columns resident;
- attempt 2 reached 4.1 GB, then crashed on H1's nullable-string gene index;
- attempt 3 completed at **11.2 GB**: the whole-matrix CSR remap of H1's 118M nonzeros. My watchdog sampled the
  `time` wrapper, not the Python child, so it did not stop it.

The machine stayed at 79–80 % memory free throughout, with swap unchanged (2,458 MB used) and the 553 worker's ~1.9 GB
job unaffected. `build595.py` now streams Orion in 1,024-row batches with the system allocator. A rebuild should also
stream H1's rows in chunks before running under a 4 GB promise.

## Files
`build595.py` (bundles), `rehearse595.py` (stages on the patched worktree), `work/build595.json` (sources, pools,
spreads, shas of the six context files), `work/rehearse595.json` / `.log` (every stage's rc, identity top 5 and
min_top1, census, same-line provenance), `work/{preflight,identity,census}_{X,Y}.json`, `work/sameline_{X,Y}.npz.json`,
and the logs of the three build attempts.
