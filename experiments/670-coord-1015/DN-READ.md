# Coordinator checklist: reading DN (10-15 s1) and the 704 contingency

Prereg: experiments/631-nmae-headroom/PREREGISTRATION-DNMAG.md (FROZEN 10-02). Base: 630arm's entry QJRJK6OC6Zw1kjwnZqid.
DN's entry id: the "entry id:" line of experiments/670-coord-1015/work/submit-1015s1.log; capture in data/board-entries/<id>.json.

1. Falsifiers first (any one: no adoption, diagnose): |Δfid|, |Δjac| > 1e-6 raw; |Δreach| > 0.0014 raw; |Δpds| > 1e-6 raw
   (one pds rank flip is 3.716e-6 = 1/(3·300·299): if Δpds is exactly that, diagnose as one flip, as ALLSW's); |Δmse| > 1e-6 raw.
2. Primary: Δnmae raw (DN − 630arm). Overall ≈ −Δnmae_raw × 1.64787 / 6.
   - ≤ −0.0025: TRANSFERS → decide `de_dnmag d_m200` (October, via 690a's DEF_DAY_DE_DNMAG pass-through), written reason in the evidence.
   - −0.0025 … −0.0008: PARTIAL → A/B/C adopt; October adopts only if the build's capacity floor holds on D/E/F's own blocks (flagged).
   - −0.0008 … +0.0008: NULL → `de_dnmag off`; launch 704.
   - > +0.0008: HARM → `de_dnmag off`; close the route; launch 704.
3. 704 contingency (only on NULL/HARM): `nohup bash experiments/704-alpha125-contingency/build704.sh >
   experiments/704-alpha125-contingency/work/build704.nohup 2>&1 &` in the next free builder window (after 702 r2 or between rehearsals).
   TAG a704; run log experiments/428-zero-gene-audit/work/build_a704/runa704.log; READY line in experiments/701-kup-callset/work/READY-science.log.
   Prereg experiments/704-alpha125-contingency/ (be4e67eb). Arm it with arm1015.sh for the next free slot.
4. Record: BOARD-LEDGER entry, experiments/671-dnmag-arm/BOARD-RESULT-DN.md, FINDINGS §0 line, NEXT-SLOT-DECISION line, october-decisions.
