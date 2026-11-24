# 601: filler-call signs on H1 (council 600, [4d0207], 2026-10-02 10:55–11:05Z)

Analysis only. Reads h1-test's truth (`323/work/reference_de.parquet`) and h1-test's control cells. No emitter, no build,
no `cell-eval2`. Peak RSS was 9.6 GB (`/usr/bin/time -l`) with 26+ GiB free and no build running.

**One-sentence result:** the reference sign of non-DE high-CPM genes ("filler", which is most of our ~2,600 calls) does
follow (a) the DE set's net mass shift and (b) a per-gene generic direction shared across targets. Neither is a new
lever: (a) is worth +0.006 filler accuracy even at oracle, and (b) is the context-programme pot that 441/491 already
priced, whose best public route reaches only a third of it.

| filler sign rule (top-2,200-CPM genes outside each target's reference DE set; mean over 85 targets) | accuracy |
|---|---:|
| all-down | 0.5115 |
| compositional: −sign(net mass change of the target's reference DE set) (ORACLE: uses the truth's DE set) | 0.5176 |
| per-gene majority sign over the OTHER 84 targets of the same line (LOTO; the ceiling of any target-agnostic rule) | 0.5413 |
| the same, on the 610 genes whose majority is ≥ 60/40 | 0.6316 |
| K562 GWPS mean response sign, raw (all-down on the same genes 0.5123) | 0.5170 |
| K562 mean response projected onto H1's own control-PC span, k 20 / 50 | 0.5234 / 0.5239 |
| Jurkat / RPE1 / HCT116 / HEK293T mean response, projected k 50 | 0.5151 / 0.5185 / 0.5152 / 0.5110 |
| H1 control PC1, oriented by the donor | 0.5106 |

- Spearman(filler up-fraction, −net DE mass change) = **0.633** over targets. The compositional coupling is real, but its
  per-gene effect is tiny except on a few large-mass targets (HMGCS1, SAFB2, MED19, CYCS).
- The control-basis projection lifts a donor's generic sign by about +0.007 (K562). That is the "from controls" route of
  441/491 (+0.0055 fid-like), re-derived independently. It reaches ~1/3 of the gap from all-down to the same-line ceiling.
- So [d4cd5f]'s §3 stands: the block's best sign is the context's own programme. It is worth about +0.017 overall at
  oracle, and no zero-shot public route has been found that reaches more than about a third of it.

Files: `filler601.py` (composition), `generic601.py` (LOTO ceiling), `ctrlprog601.py` + `ctrlprog601.json` (control-basis
projection), `filler601.csv`.

## Addendum (11:10Z): the generic filler sign by target strength (`strata601.py`)
Per-gene majority sign built from one stratum's targets (leave-one-out within it), evaluated on each stratum's filler:

| majority from \ evaluated on | null (n_real ≤ 5, 16 tg) | weak (6–100, 28) | strong (> 100, 41) |
|---|---:|---:|---:|
| null | 0.5190 | 0.5182 | 0.5161 |
| weak | 0.5186 | 0.5286 | 0.5256 |
| strong | 0.5138 | 0.5263 | **0.5577** |
| all-down (same filler) | **0.4958** | 0.5109 | 0.5180 |

- Strong targets share a sizeable programme (+0.040 over all-down within the stratum). Null targets' filler leans slightly
  UP, so an all-down block sign is below chance on them (0.4958). The 84 %-down CLEAN4 block therefore loses a little
  fid on the negligible bin, where the information is weakest.
- Size: flipping null-like targets' block to neutral or up is worth about +0.008 filler accuracy on ~19 % of targets,
  so ≈ +0.0015 fid raw. That is a rider at best, and it needs a strength key (521: r 0.19–0.43 on cancer lines). Not
  proposed.

## Addendum 2 (11:12Z): panel-restricted donor programme (`panelprog601.py`): NULL
Restricting a donor's mean response to its rows for THIS panel's targets (known on the day) does not beat its
genome-wide mean on H1's filler: K562 GWPS 0.5155 vs 0.5170 (82/85 covered), HCT116 0.5082 vs 0.5129, HEK293T 0.5102 vs
0.5108; all-down 0.511–0.512; same-line ceiling 0.5413. Donor programmes carry almost none of the recipient's programme,
whichever way they are averaged. This matches 450h/462/466.
