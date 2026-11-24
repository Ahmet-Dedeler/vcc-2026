# 611 — Calibrating the H1 DE instrument against the board (desk check; council 600 gate for DE dose arms)

> **THE RULE (load-bearing for every council DE dose arm):**
> - **H1 prices fid + jac and nmae.**
> - **It never prices reach for any arm that changes the call set** (composition or size). cp's reach is mis-signed in
>   all four regimes.
> - For delivery-only arms (calls unchanged), one pair (coarse MDE) matched reach at f 0.35 (1.16×). Treat that as
>   provisional and READ reach on the board inside a preregistered band of ±(H1 |Δreach| + 0.004).
> - Price at f 0.35:
>   - nmae × 1.2–1.6 (delivery / magnitude levers);
>   - fid × ~2 (cost) and jac × ~0.5 (gain) for call-set changes.
>
> **AMENDED 2026-10-04 (K20's board read; 662):**
> - **Call-set changes (composition or size) are priced at f ≈ 0.61 (0.60–0.62), NOT f 0.35.**
>   - K20's fid moved −0.0017 on the board.
>   - H1 at f ≤ 0.5 has the wrong sign; f 1 over-reads 6× and mis-signs jac.
> - **At f 0.61:**
>   - fid × 2.7 on h1-test (one value for cp and K20), × 0.8 on 618;
>   - jac × 0.43–0.49;
>   - delivery nmae × 0.92–1.37 on both panels.
> - The f 0.35 call-set lines above are void: the K16 / K20 card is graded FAIL and K16 is re-registered in
>   `662/PRED-662-K16.md`.
> - See `experiments/662-callset-regime/RESULTS.md` and 619's amendment.

2026-10-02 ~12:00Z, council session [21dd4c]. Analysis only: no build, no `vcc`, no `cell-eval2`, every process < 2 GB.
Scripts: `desk611.py` (cp from 537's `price537.json`), `desk611b.py` (zs and coarse MDE re-scored here). Outputs are
`work/desk611.out`, `work/desk611b.out` and `work/desk611b.json`.

## Question
The council's pivot rule (SYNTHESIS v4) sends DE dose arms to the board only after pricing them on an H1 instrument that
retrodicts the board pairs we already hold. Which pairs can be priced, which regime of 521's four (f 1 / 0.5 / 0.35 =
PowerSurface pools, cmp) matches the board, and on which members?

## Which board DE pairs can be priced on H1 without leakage
| board pair | what changed | H1 equivalent | priceable? |
|---|---|---|---|
| 446bj2z − 446bj2 | zs exact-U post-pass (nmae only) | 435 coarse vs coarse+zs, seeds 0–2 | **yes** |
| 433a − 417e | coarse MDE (delivery on uncalled genes) | 426 production-emitter no-MDE vs 435 coarse, seed 0 (bulk sums bit-equal) | **yes** |
| 446bj2zcp − 446bj2z | 528's call-propensity block (call-set composition) | 537 `lcap-kept` vs base, seeds 0–2 | **yes** |
| 446bja − 446b | Jurkat programme as context A's block prior | one-context prior; H1 is one context | no |
| 458d − 446b, 458g − 458d | Arc's H1 programme as block prior (± gate) | Arc's H1 programme on H1 is the line's own (leakage) | no |
| 446br2 − 446b | zero-repair floor | H1 has no zero-total gate genes | no |
| 433h − 433a | zero repair + re-level | same | no |
| K16 / K20 (10-06 / 10-07) | call-set size | 537 `cpm-n1600` / `cpm-n2000` | yes, once read |

## Board vs H1, member by member (raw Δ; board / H1 ratio)

| pair | member | board Δ | H1 f 1 | f 0.5 | **f 0.35** | cmp |
|---|---|---:|---:|---:|---:|---:|
| zs | nmae | −0.00218 | −0.00183 (1.19×) | −0.00182 (1.20×) | **−0.00185 (1.18×)** | −0.00318 (0.69×) |
| coarse MDE | nmae | −0.01679 | −0.00966 (1.74×) | −0.01101 (1.52×) | **−0.01059 (1.59×)** | −0.01263 (1.33×) |
| coarse MDE | reach | +0.00240 | −0.00549 (sign ✗) | −0.00127 (sign ✗) | **+0.00207 (1.16×) ✓** | −0.00055 (sign ✗) |
| coarse MDE | fid, jac | ≈ 0 | ≈ 0 | ≈ 0 | ≈ 0 | ≈ 0 |
| cp | fid | −0.00358 | −0.00126 (2.83×) | −0.00176 (2.03×) | **−0.00183 (1.95×)** | −0.00183 (1.95×) |
| cp | jac | +0.00399 | +0.00910 (0.44×) | +0.00887 (0.45×) | **+0.00829 (0.48×)** | +0.00918 (0.43×) |
| cp | reach | −0.00335 | +0.00956 (sign ✗) | +0.00851 (sign ✗) | **+0.00888 (sign ✗)** | +0.00860 (sign ✗) |
| cp | nmae | +0.00010 | +0.00043 | −0.00044 | −0.00013 | −0.00190 | (noise)

## Reading
1. **Use f 0.35.** It is the only regime that gets every sign right on the delivery pairs (zs, MDE), including the MDE's
   reach (+0.0021 vs +0.0024 on the board). Its ratios are the most stable: nmae 1.2–1.6×, MDE reach 1.16×. cmp
   over-reads zs (0.69×). f 1 and f 0.5 mis-sign the MDE's reach.
2. **Member multipliers at f 0.35**, for converting an H1 price into a board expectation:
   - nmae from delivery and magnitude levers: × 1.2–1.6, the same direction as FINDINGS' "magnitude 1.6–2.3×".
   - fid cost of a call-set change: × ~2 (one pair).
   - jac gain of a call-set change: × ~0.5 (one pair).
   - reach of a DELIVERY change (no call-set change): trusted at ~1× (one pair).
   - **reach of a CALL-SET change (composition or size): NOT trusted.** cp's reach is mis-signed in all four regimes.
3. **What this licenses for the pivot sweep:**
   - A rank-neutral mean-delivery arm (nmae) changes delivery, not calls. It can be priced at f 0.35 with nmae
     × 1.2–1.6 and reach × ~1. **Licensed.**
   - A reach-order arm that leaves calls unchanged (it moves only p-values of uncalled genes, i.e. the order) is
     delivery-like. Price it at f 0.35 but read it with one caveat: no board pair of that exact kind exists yet.
   - A dispersion-law arm that changes the call set (William-type) is call-set-like. fid × ~2 and jac × ~0.5 are
     usable, but **reach must be read on the board**. Preregister a reach band of ±(H1 |Δreach| + 0.004).
   - Call-set size (K16/K20) will add the first size points on 10-06/07. Re-run this table then (`desk611.py` takes
     new pairs in one line).
4. **The evidence is thin:** three pairs, one of them multi-member. Every multiplier above rests on one or two pairs. They
   are planning numbers, not calibrated constants. The sweep should prefer arms whose expected effect is ≥ 3× the
   member's band.

## One line for the coordinator
H1 at f 0.35 retrodicts the board's two delivery pairs (zs nmae 1.18×; coarse MDE nmae 1.59× and reach 1.16×, signs
right). It mis-signs reach for the one call-set pair (cp) in every regime. Price council dose arms at f 0.35 with nmae
× 1.2–1.6, fid × 2 and jac × 0.5 for call-set changes. Trust reach only for delivery changes.
