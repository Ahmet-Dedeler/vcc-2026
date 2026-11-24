# 510 — the October recipe (as decided 2026-09-27) rehearsed end to end on efr with a forced in-mix same-line call — PASS

2026-09-28 17:39–19:26Z (`chain510.sh`, after 497e). def_day e29e793a, builder 24b95dec (496p's emitter pin in; the hook off).

## Configuration
`DEF_DAY_FORCE=D=K562…,E=HCT116…,F=A549_LUNG DEF_DAY_LODO=1 DEF_DAY_LODO_HI_PCT=96.64697522769474
DEF_DAY_SHARPEN_PT=stercile:2.5,2.0,1.75 DEF_DAY_PROGRAMME_SOURCE=h1_2025_full DEF_DAY_PRIOR_GATE=le1
def_day.sh ef-bundle-r … 510efr 1.75 0.5 "" 4` — 487's same-line setting (w 0.5, no pt, the programme on) plus every
candidate switch. The forced call puts two in-mix banks (K562, HCT116) and one sparse bank (A549) on the three contexts.

## Result
- rc 0, **DONE**; postbuild PASS vs 446b with the allowed diffs; zero-gene audit 0 in 0 blocks (lfc < −5: 6, > +3: 10).
- Converged after 3 levelled builds at D 1.2353 / E 1.3276 / F 1.1546; final ρ D 0.9730, E 0.9721, F 0.9797 (all ≤ 1;
  D and E sit 0.6–0.7 % under ρ* 0.9793, inside the loop's tolerance). Builder restored.
- Same-line: w 0.5 applied on every context (the artifact's auto w of 2.0 for the in-mix banks is overridden by the
  explicit argument, as intended); pt off; Arc's H1 programme is each context's block-prior default (down fraction 0.560);
  D and E change all 300 rows.
- Prior gate le1: D 141 strong targets on the programme / 159 weak on CLEAN4; E 133 / 167.
- LODO key under the percentile cap: hi 7.476 (the same efr panel as 497e). Per-target sharpening: stercile 100/100/100.

## Reading
The October recipe as decided on 09-27, with every switch still PENDING turned on, composes and completes on an
October-like panel with a same-line match. The two remaining engineering items for the day are 502's DE-rider
pass-through (DEF_DAY_DE_MAG, chain502 tonight) and the recipe's final switch values after the 09-30 … 10-05 reads.
