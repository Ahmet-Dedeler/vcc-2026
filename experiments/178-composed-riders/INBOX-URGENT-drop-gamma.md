# URGENT INBOX for 178 — drop gamma from the composite; read `experiments/170-gamma-mse/RESULTS.md` first

Coordinator, 2026-09-02 20:50. Experiment 170 (which the coordinator believed killed this morning)
completed a full production-staged gamma sweep with bit-identical known-answer gates
(`baseline_pre`/`candidate_pre` sha256 = a15's audit; a15's and a19's post-whitening audits
reproduced to 1e-15). Its verdict supersedes 097 §9h and your brief's first rider:

## Gamma is a TERMINAL NO-GO at every setting in [0.35, 0.80]

| gamma | mse | pds | nmae | reach | fid | jac | **TOTAL** |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.35 | −0.00918 | −0.00150 | +0.00107 | +0.00110 | +0.00119 | +0.00019 | **−0.00713** |
| **0.45** | **−0.00629** | +0.00005 | +0.00073 | +0.00075 | +0.00008 | +0.00001 | **−0.00466** |
| 0.55 | −0.00202 | +0.00046 | +0.00041 | −0.00090 | +0.00035 | +0.00006 | −0.00163 |
| 0.65 | 0 | 0 | 0 | 0 | 0 | 0 | ships |
| 0.80 | −0.00334 | −0.00083 | −0.00114 | +0.00041 | −0.00085 | −0.00013 | −0.00588 |

**Raw `mse` on the challenge panel is V-shaped with its vertex at the shipped 0.65** (0.926 / 0.881 /
**0.874** / 0.880 / 0.893 at 0.35 / 0.60 / 0.65 / 0.70 / 0.80). **Mechanism: the ρ cliff.** Lowering
gamma raises the effect's mean-square `A` by 14% (0.35 → 0.65) and with it the across-perturbation
budget `B`; a15 sits essentially on the cliff (ρ 0.99 / 0.94 / 0.88 per context), and raising `B`
past `claim` saturates `min(claim, B)` so the extra amplitude is not refunded. 097's h1-test
`mse` column assumed the arm sits ON the cliff; it does not, and h1-test cannot price `mse` at all
(its raw `mse` is 0.96–1.02, above the clamp). **Remove gamma 0.45 and gamma 0.30 from every arm.**
a19 re-attributed: gamma is 21% of its −0.0959, `lam` 79%.

## What this means for the rest of your composite — a sharper criterion

170's closing paragraph is the design rule: **"a transform that raised `X` without raising `A` would
convert the full gain … norm-preserving, not norm-adding."** Lowering gamma buys +16% alignment `X`
over [0.65, 0.35] but pays for it with amplitude the metric will not refund above the cliff.

Check each remaining rider against that criterion, and **measure `A`, `B`, `claim` and ρ per context
on the challenge panel for every arm** (170's `chal.py` does exactly this — reuse it; do NOT
linearise `min(claim, B)` around the shipped point, which is the error 170's own first pass made and
corrected):

- **Sign repair** (173): flips signs, keeps `|H|` per entry → `A` unchanged by construction,
  `X` +7.8% → **norm-preserving.** This is the rider 170's criterion endorses. Expect `mse`
  ≈ +0.0027 overall from the `X` gain alone (170: +16% X ≈ +0.0055), against 173's `pds` −0.0042
  raw. Use the **deflated** consensus as the sign source (INBOX-fourth-rider addendum 3).
- **K562-heavy weights (2,1,1)**: changes the consensus → changes `A` (Orion's effects are 1.7×
  larger, so down-weighting them LOWERS `A`, which may move ρ the *favourable* way — measure it).
- **Per-donor deflation**: removes a direction → lowers `A` slightly; on the whitened profile it is
  ~nil directly (0.11% of energy), acts through the whitener refit and the sign source.

The composite is now **sign repair (deflated source) + weights (2,1,1) + deflation**, with gamma
fixed at 0.65. Report `A`/`B`/ρ per context alongside the six members for each arm — whether an
arm stays below the cliff is now the first thing to know about it.

## Also: 170 found a small seed-free `pds` interior maximum at gamma 0.50–0.55 on the DELIVERED object

+0.00046 overall, seed-free (`cell_eval2.discrimination_score` exact). It is swamped by the `mse`
loss at those gammas, so it is not actionable — but it means 097's "pds monotone in gamma" was true
of the *intended* effect and not of the *delivered* (mass-fixed) one. Do not cite 097's pds curve.

## Operational corollary: re-tune `VCC_ALPHA` onto the knee for EVERY arm that moves `A` or `B`

From 096/171/170 the scored numerator is `s²A − 2sX + Kc − min(s²B, claim)`. With `A − B` tiny
(0.0013) the unconstrained optimum `s = X/(A−B)` is far above the knee, so **the `mse` optimum is
pinned AT the knee `s* = √(claim / B)`**, where

```
raw mse* = 1 + ( claim·A/B  −  2·X·√(claim/B)  +  Kc  −  claim ) / sumD
```

Two consequences:
1. **At the knee, `mse` improves with `X/√B`, not with `X` alone.** A rider that raises `X` and `B`
   together (rescaling) buys nothing; one that raises `X` at fixed `B` (sign repair) buys the full
   amount. That is 170's "norm-preserving" rule in one formula.
2. **Every arm must be scored at ITS OWN knee.** a15's `VCC_ALPHA = 0.982` puts a15 on the knee; a
   rider that changes `B` (weights: Orion effects are 1.7× larger, so (2,1,1) lowers `B`; deflation
   lowers `B` slightly) is off the knee at 0.982 and will read a spurious `mse` loss. For each arm,
   measure `B_new` per context on the challenge panel and set `alpha_new = 0.982 · √(B_a15 / B_new)`
   (per context if the builder allows, else the mean), THEN price. Report both the un-retuned and
   the knee-retuned `mse` so the size of this correction is on the record. `pds` is scale-invariant
   so the retune costs it nothing; `nmae`/`fid`/`reach` move with amplitude and must be read after
   the retune too.

170 §6 calls this "a joint (alpha, gamma) screen" and says it is the only open version of the gamma
question; with gamma fixed at 0.65 it is the joint (alpha, rider) screen, and it is the difference
between a fair and an unfair pricing of your composite.
