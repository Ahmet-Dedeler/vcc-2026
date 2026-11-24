# 510 — the 09-27 October recipe rehearsed on efr with a forced in-mix same-line call

Why: 487 changed the same-line setting to w 0.5 with no per-target signs. Every earlier rehearsal (452, 452b, 465, 469)
ran with `pt`. This run puts every October switch as decided on 09-27 through `def_day.sh` once, on an October-like panel,
with the same-line call forced the 452b way (D=K562 and E=HCT116, whose banks are already in our mix; F=A549, sparse).
Switches: LODO + percentile cap (P 96.64697522769474), S_P2, programme h1_2025_full, gate le1, same-line args `0.5` `""`.
Checks: DONE, postbuild PASS, every final ρ ≤ 1, the builder restored, and the audits recorded. `chain510.sh`, serial after
chain497e. It is a machinery check, not a scored arm.
