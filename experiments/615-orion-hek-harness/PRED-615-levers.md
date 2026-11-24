# 615: lever predictions on the Orion raw HEK293T harness, registered BEFORE any lever run (and before the gate verdict)

Written 2026-10-02 ~18:45Z by [21dd4c]. At this point the reference was still being built: no 615 emitter run existed,
and the power diagnostic had not been read. Raw nmae, lever − its base, f 1 unless stated, mean of seeds 0–2.
- **Gate dependence:** these predictions are READ only if the 615 gate PASSES (TRUSTGATE.md, 6e35ff09). Otherwise the
  runs are descriptive.
- **The key, beside every row (coordinator):** 615's design key is k562_gwps alone (269 / 300 panel coverage), a
  single-line proxy. Production's key is multi-donor (3 lines) and whitened, so the call set (head) and the uncalled
  set differ from production's.
- Every lever row carries `--zero-repair -5`.

## Hook forms (base = G1's zs construction: prodem, coarse MDE, zero repair −5, `--zs dn_k250_zs`)
- **Source:** 613d (`613/work/block613d.json`, `block613d5.out`). Its HEK293T fold is 636's instrument over 480's
  gates (the add-rule model, not this emitter), with two differences from 615:
  - its truth is the THINNED Orion HEK harness, i.e. the greedy-selection blob 615 just removed (RESULTS stage 0);
  - its key is K562 + HCT116, which includes a same-platform sibling.
- 615 has a clean truth and a K562-only key, so this is a genuine cross-surface test of 613d's HEK row. It is not a
  re-reading of the same data.
- **CLEAN4 is complete here** (all four banks; none is Orion-derived), as it will be on D/E/F. So the HEK row is the D/E/F
  analog, unlike 614's weakened-CLEAN4 K562 row.

| contrast | 613d HEK293T fold, f 1 / 0.5 / 0.35 |
|---|---|
| c4h − no hook | −0.00215 / −0.00069 / −0.00076 (the plain hook helps) |
| idall3 − c4h | −0.00386 / −0.00262 / −0.00111 (identity helps more) |

**Registered (f 1, the regime 615 reads natively):**
- **R1:** c4h − no hook < 0, and idall3 − c4h < 0.
- **R2, size:** each contrast within a factor of 3 of 613d's HEK value (ratio in [1/3, 3]). That is c4h − no hook in
  [−0.0065, −0.0007] and idall3 − c4h in [−0.0116, −0.0013].
- **Reading:**
  - R1 + R2 → 613d's hook rows hold on a clean truth with a different key, for the D/E/F-like (CLEAN4-complete) case. It
    is the instrument's first test off its own gate tables.
  - R1 only → the sign transfers but 613d's size does not; the hook's price stays the H1 real-emitter number.
  - R1 fails → 613d's Orion rows were the blob's (or the sibling key's), and the hook's held-out support rests on H1
    alone.
- **How it would be built:** `run615.py hekc4h / hekidall3` (`idsign615.py`: c4h = the four banks' mean P_down at
  τ 0.05; idall3 = 604 v3's rule with HEK's own controls), 3 seeds, against G1's zs-arm stats (`g1/…_dnk250zs`).

## Not registered yet (needs the coordinator's lever list with driver commands)
- CPH, DN, K16 / K20-class call sets, 630 v3 `b_k100_q20r3p10`, 631 / v7.
- Each one gets a prediction row here, committed before its first HEK run, and a d20 note:
  - 633's d20 is a fold-model rescale of the add's counts;
  - this real-emitter harness cannot emit at 20k from a 15k truth;
  - so the depth reading for an add-based lever comes from H1 run at 15k vs 20k (h1-test thinned, as H1@11k in 614),
    not from rescaling.
