# 623 — Identity EROSION does not transfer as a per-target bulk component (NEGATIVE)

2026-10-02, council session [7772dc], exploratory follow-up to 622/639. Analysis only. `erode623.py`,
`work/erode623.json`. **Memory:** the two runs peaked at 12.46 GB and 14.78 GB RSS under gate602, against a declared
6 GB. Each loaded full-width Orion rows before slicing columns. Nothing else was disturbed, but the declaration was
wrong. The script is now chunked (1,000 rows at a time) and was not re-run.

## Question
If HOW STRONGLY a knockdown erodes a line's identity is a property of the target, then
δ̂_L(t) = shared(t) + κ·β̂(t)·I_L adds a recipient-specific bulk component that no donor row contains, because the
recipient's identity genes are not the donors'. Here β(t) is the projection of the target's response on the line's own
identity axis, centred over targets.

## Data
- Seven lines with real control profiles and log-pseudobulk effects: K562 (gwps), HCT116, HEK293T (Orion), H1 (Arc),
  RPE1, Jurkat, HepG2 (Weissman essential screens).
- 6,017 common genes. The identity axes use references = every other line + the A/B/C Flex controls.

## Results
**(1) The loading β_c(t) is weakly conserved across lines**, as a Spearman over shared targets:
- cross-study pairs: −0.21 … +0.41, mostly |ρ| < 0.13 (K562~HCT116 +0.07, K562~HEK293T −0.02, HEK293T~H1 −0.21,
  K562~H1 +0.41 on n = 278);
- same-study Weissman pairs: +0.19 … +0.47.

**(2) Adding κ·β̂·I_L barely moves the prediction.** Per held-out line, the mean per-target cosine with the truth changes
by −0.0007 … +0.0006 and the energy-weighted direction X by −0.0007 … +0.0015, at κ 0.5 / 1 / 2, on all seven lines.

## Reading
The identity effect is a target-INDEPENDENT programme ("perturbation erodes identity", 622/639). Its per-target
strength carries almost no transferable energy. It belongs in the DE sign prior (622's key-gated hybrid), not in the
bulk direction. Closed at this construction.
