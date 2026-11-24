# Instrument registry for DE-channel levers (council 600, [21dd4c], 2026-10-02 ~13:40Z)

Which held-out or H1 instrument may price which lever class, and the retrodiction record behind each. The rule (AGENTS
11.1 rule 3): an instrument prices a lever class only after it has retrodicted a board pair, or a real-emitter run, of that
class. Before a slot, write the instrument's prediction down, so the read tests it (coordinator, 10-02).

| instrument | prices | retrodiction record | status |
|---|---|---|---|
| **H1 real emitter at 521's f 0.35** (537/612 harness; 611) | DELIVERY levers (calls unchanged): nmae × 1.2–1.6; reach provisional, band ±(\|Δ\| + 0.004) | zs nmae 1.18×; coarse MDE nmae 1.59×, reach 1.16× (611) | **validated** for delivery. Never for reach of a call-set change |
| H1 real emitter, call-set changes (611) | fid × ~2 (cost), jac × ~0.5 (gain) | one pair (cp) | provisional. K16 / K20 (10-06 / 10-07) are its first out-of-sample tests |
| **Restricted-LOLO block instrument** (613; 636's surfaces, block = head-excluded call set) | BLOCK-PRIOR source changes: sign precision vs CLEAN4 on gate AND null block genes | 458d HARM (sign, and size within 0.43 dose); 392 fid gain (sign) | **validated** for block priors (2 board points) |
| **Add-rule uncalled instrument** (613d; 426's ±add / ±2·add rule on uncalled gate genes, 521 pools) | HOOK-class levers (the MDE add's direction on uncalled genes) | H1 real emitter: c4h − base and c4kg − c4h within 0.0004 at all 3 regimes | **2 out-of-sample PASSES** on H1's real emitter, both within 0.0003 at all regimes: idkg3 (A3, 14:52Z) and idall3 (15:21Z). Pending: 496a's board read (10-05; registered +0.00017), its cross-surface exam. Third known-answer row (KGID-v1) errs up to 0.0011 on uncovered genes |
| Fold nmae model, base = BETA·consensus (kgh612b/c/d, hyb622b, ds612c) | — | wrong-signed on H1 at f 0.35 for hooks (612 A1) | **VOID** for hooks. It hooks called genes, and its base magnitudes are unfaithful |
| Zero-base sign × MDE model (kgh612, dosecurve610) | — | over-stated the hook ~3× | **VOID** |
| Thinned-cell H1 surfaces (h1_c100 / c150 / thin) | — | no board point | a different question (noisier truth), not 521's pools |

## Reading rules that cut across instruments
1. **The K562 fold is a weakened-CLEAN4 scenario, not a typical new line.** On it, CLEAN4 drops K562's two banks and keeps
   only Nadig Jurkat + HepG2. A D/E/F line gets all four banks, as the Orion folds do. So the "K562 tail" arguments (the
   plain hook's −0.0013 … −0.0019 there; identity −0.029 on the block gate) price a CLEAN4 that is weaker than any D/E/F
   line will see. Weight the Orion folds as the D/E/F analog, and keep K562 as a stress test.
2. **The Orion folds share one platform.** On HCT116 the consensus key includes HEK293T, and the reverse. So their head
   captures more of the gate than production's will, and their absolute precisions carry the Orion technical blob. Read
   only differences between sources on them.
3. **Depth matters for add-based levers.** Orion runs at 12–14k and K562 at 11k median library, against ~20k on H1 and
   A/B/C. Re-run an add-based price at A/B/C's depth before quoting a footprint (613e retraction; 613d3).
4. **Pool regime.** The board's n_real sits near H1's f 0.35–0.5 (192 §4.2). A lever whose H1 contrast turns sign
   across f 1 / 0.5 / 0.35 (KGID − plain) cannot be called from H1 alone.
