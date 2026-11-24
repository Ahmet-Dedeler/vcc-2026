# a13 — operating point on the CORRECTED shrinkage estimator

Written 2026-08-28, before any a13 build exists. Supersedes the version solved against the
pre-fix estimator, which is retracted in full: it ranked the ~16.5k zero-agreement genes by array
position and was a different estimator. Everything here imports the shipped
`agreement_shrink_vector`, `donor_cosine_key`, `pert_shrink_vector` and `MultiDonorModel` from
`scripts/make_submission_pm_cd4.py` at HEAD rather than reimplementing them.

Baseline: live board entry **a11** — overall `0.16611503`, rank 7/391, raw `mse` `0.9099893`,
`score_mse` `+0.082041`. Host: **a12** (`--keep-top 8000 --keep-spread 0.25 --tail-spread 1.0
--overdispersion 4.0`, pruned whitener), predicted raw `mse` `0.8881`.

Scripts: `mech_a13_resolve.py`, `mech_a13_bound.py`, `mech_a13_pds.py`,
`mech_a13_donorcount.py`, `mech_a13_grid.py`, `mech_a13_claim.py`, `mech_a13_floored.py`.

**Pinned to builder commit `2f7020c` (2026-08-28 20:42:36 +0300)**, with
`PERT_SHRINK_GAMMA=2.0`, `PERT_SHRINK_FLOOR=0.3`, `PERT_SHRINK_MISSING='damp'`,
`AGREEMENT_SHRINK_RANK_SCALE=2000.0`. ⚠️ An earlier pass of this solve imported the same
functions from an **earlier revision** of the file, when `PERT_SHRINK_FLOOR` was `0.0` and
`missing` did not exist. It was an import, not a reimplementation — but of a version that no
longer exists, which is the same defect one level up. Every number below is re-run against the
commit named here; re-verify the hash before building.

---

## 0. Applicability — orion cannot run this estimator, and that is now controlled

You asked whether the challenge differs from the harnesses in a way that makes the harness answer
*inapplicable* rather than merely noisy. **It does, and I nearly shipped the wrong surface.**

Production keys on **three** donors (`k562_gwps`, `orion_hct116`, `orion_hek293t`). `orion_hct116`
is the same cell line as `orion-hct116-thinned`, so that harness can only ever run a **two**-donor
key. Unanimity across three donors is a stricter gate than across two, so the key is far sparser:

| surface | donors | median non-zero genes / target |
|---|---:|---:|
| **challenge** | 3 | **1,979** |
| **h1-test** | 3 (leak-free — no donor is H1) | **1,964** |
| orion-thinned | 2 (forced) | 3,820 |

**The controlled test.** Running h1-test with orion's forced two-donor key reproduces orion's
behaviour and reverses the sign of the `pds` result:

| γ | h1-test, **3 donors** Δ`score_pds` | h1-test, **2 donors** Δ`score_pds` | orion, 2 donors |
|---:|---:|---:|---:|
| 0.25 | **+0.0266** | −0.0084 | −0.0522 |
| 0.50 | +0.0074 | −0.0223 | −0.0768 |
| 1.00 | −0.0105 | −0.0338 | −0.0806 |

So the disagreement is **donor count, not cell line**, and orion's answer is not a second opinion
— it is a measurement of a different estimator. **h1-test at three donors is the only applicable
surface, and rule 1's "take the smaller of two" does not apply here.** That is a fourth instance
of a surface being structurally incapable of exhibiting the thing it is being asked about.

---

## 1. The (γ, effect_scale) pair

Exact identity throughout, `VCC_ALPHA = 1.0796`, `A′`/`B′` computed exactly on the challenge
panel, `X′/X` transferred from h1-test-3-donor only. `Δ` columns are overall points against a12.

| γ_gene | γ_pert | effect_scale | raw `mse` | Δ`mse` | Δ`pds` | **sum** |
|---:|---:|---:|---:|---:|---:|---:|
| 0.000 | 0.0 | 0.9375 | 0.8881 | +0.00002 | +0.00000 | +0.00002 |
| 0.000 | 1.0 | 1.6318 | 0.8753 | +0.00226 | +0.00000 | +0.00226 |
| 0.125 | 1.0 | 2.1574 | 0.8407 | +0.00833 | +0.00470 | +0.01303 |
| **0.250** | **1.0** | **2.4004** | **0.8212** | **+0.01175** | **+0.00444** | **+0.01619** |
| 0.250 | 0.0 | 1.4068 | 0.8409 | +0.00829 | +0.00444 | +0.01273 |
| 0.500 | 1.0 | 2.8006 | 0.8077 | +0.01412 | +0.00124 | +0.01535 |
| 1.000 | 1.0 | 3.1766 | 0.7889 | +0.01741 | −0.00176 | +0.01566 |

**The floor changes the answer by less than the grid spacing.** Re-solving with
`floor = 0.3` applied INSIDE the exponent, as the code does, against `floor = 0.0`:

| | `mean z` | `mean z²` | `effect_scale` | raw `mse` | sum |
|---|---:|---:|---:|---:|---:|
| γ_pert=1, floor 0.0 | 0.500 | 0.333 | 2.4814 | 0.8225 | +0.01596 |
| **γ_pert=1, floor 0.3** | **0.545** | **0.351** | **2.4004** | **0.8212** | **+0.01619** |
| γ_pert=2, floor 0.3 *(shipped default γ)* | 0.351 | 0.202 | 3.1188 | 0.8281 | +0.01498 |

`effect_scale` moves **2.4814 → 2.4004**, −3.3 %, and the floor is mildly *helpful*. My `mean z²`
of 0.351 against the implementing session's 0.347 is an independent cross-check on the same
quantity. **γ_pert = 2.0, the shipped default, is worse than 1.0 at every γ_gene** — pass it
explicitly.

**Recommended: γ_gene 0.25, γ_pert 1.0, floor 0.3, `--effect-scale 2.4004`.**

**Why not the shipped default γ=1.** Raw `mse` is **monotone** in γ — it keeps falling to γ=3 and
beyond. An unbounded gain means the model is being extrapolated, so I looked for what bounds it.
Three candidates do **not**: profile clipping is 0.08–0.25 %, `n_pred` is *flat and slightly up*
(4724 → 4876 → 4983 → 4987 across γ 0…2), and peak per-gene amplitude stays ~1.0–1.3× production
at every γ because `effect_scale` compensates. **`pds` is the bound**, and it is the one member
the `mse` solve structurally cannot see. On the applicable surface `pds` *improves* at γ ≤ 0.5 and
turns negative past γ ≈ 0.7. γ=0.25 is both the joint optimum and the only setting where `pds`
gains rather than pays.

**`claim` does not move — the omission all three routes shared is now closed.** Re-emitting each
arm and re-measuring: `Ĉ_p/C_p` = 1.14 (a12) / 1.18 (recommended) / 1.29 (aggressive), so the #247
cap binds in every case and `claim` stays at `0.246529` exactly. Holding it fixed was correct, and
not by luck — a12's tail policy over-disperses past the cap and shrinkage does not undo that.

---

## 2. The per-perturbation multiplier: **same arm, and it is free of the attribution problem**

Put it in. Two reasons, one structural and one measured.

**`pds` is exactly neutral to it.** A per-perturbation scalar `z_p` rescales row `p`, and `pds`
normalises every row before taking cosines, so `z_p` cannot move it — measured identical to four
decimals at γ_pert 0 and 1, at every γ_gene. **So any `pds` movement in the arm attributes
entirely to γ_gene**, which is what would otherwise make a two-knob arm unreadable.

**The `mse` contributions are near-additive.** γ_pert 0→1 adds `+0.00226` alone and `+0.00323` on
top of γ_gene=0.25 — mild super-additivity, not redundancy, so the two knobs are carrying
different information (per-gene magnitude-gated-on-sign vs whole-vector direction agreement).

γ_pert = 1.0 is the argmax at every γ_gene; **2.0 — which is `PERT_SHRINK_GAMMA` at HEAD — is
worse everywhere**, so it must be passed explicitly rather than defaulted. `pds` remains exactly
neutral to γ_pert *and* to the floor: identical to four decimals across all five (γ_pert, floor)
combinations at every γ_gene, re-confirmed on the floored estimator.

**The residual risk is magnitude, not attribution:** `effect_scale` 2.48 is a large compensation.
It is safe on the two things I could check — clipping 0.09 %, p99 per-gene effect 1.04× production
— because α and `z` shrink first and the scale only restores what they removed. If you would
rather not spend a slot on a 2.5× scale, `γ_pert = 0` costs `−0.0032` and needs only 1.41.

---

## 3. Predicted six-member outcome vs a11, with falsifiers

| member | predicted Δ scaled | Δ overall | basis |
|---|---:|---:|---|
| **mse** | **+0.0934** | **+0.01556** | raw 0.9100 → 0.8212, exact identity |
| **pds** | **+0.0266** | **+0.00444** | h1-test-3-donor, γ_gene=0.25. **Single surface** |
| reach | +0.010 ± 0.012 | +0.00167 | net amplitude ≈ +14 % vs a11 (a11's −22 % cost −0.0114); agreement reordering should help the purity bar |
| nmae | −0.003 ± 0.008 | −0.00050 | amplitude partly restored; a11 *gained* +0.0026 on the cut |
| fid | −0.002 ± 0.006 | −0.00033 | `n_pred` ≈ 4,900, coverage slack confirmed by a11 |
| jac | −0.002 ± 0.004 | −0.00033 | union grows with `n_pred` |
| | | **+0.0205** | **overall ≈ 0.1866** |

**Falsifiers.**
1. raw `mse` outside **[0.79, 0.86]** → the exact identity or the single-surface `X′/X` transfer is
   wrong. The last two predictions landed at 0.3 % and 0.9 %; this one carries more extrapolation.
2. **`score_pds` falls at all** → the `+0.0266` is a single-surface artefact and the sign is
   wrong. This is the prediction I am least sure of and the most consequential: `pds` is our best
   member. If it fires, re-run at γ_gene = 0.125 rather than abandoning — the `mse` gain survives.
3. `n_pred` < 4,000 → shrinkage is cutting the call set, against direct measurement (4,876–4,987).
4. `score_reach` falls by more than 0.02 → the agreement reordering hurts, which would also
   undercut the key's rationale.

## 4. Is it worth shipping at all?

**Yes** — but the case now rests on `pds` as much as `mse`, and `pds` rests on one surface. The
honest bound: if the `pds` gain is real, `+0.0205`; if `pds` is flat, `+0.0162`; if `pds` moves by
the 2-donor amount (`−0.0084`), `+0.0148`. **All three are positive**, because the `mse` term
dominates and is measured on the panel itself. The downside case needs `pds` to fall by more than
`0.12` scaled, which is 4.5× the worst number any surface produced.

Smallest build that would pin what is still open: the recommended arm itself. There is no cheaper
probe — `pds` needs the reference, `reach`/`fid`/`jac` need the DE path against the real panel, and
both are exactly what a submission buys.


---

## 5. The exact flag string

Every flag name, choice and default below was read from the parser at commit `2f7020c`.
`VCC_ALPHA` is an env var (`make_submission_pm_cd4.py:348`) whose default is exactly `1.0796`,
so **do not set it** — leaving it unset is what I solved against.

```
# a13 = a12's emission + the corrected shrinkage, at the solved operating point.
# Do NOT export VCC_ALPHA; its default 1.0796 is the value A/B were measured at.
  --keep-top 8000 \
  --keep-spread 0.25 \
  --tail-spread 1.0 \
  --overdispersion 4.0 \
  --basis-prune essential \
  --shrink-key agreement \
  --shrink-toward zero \
  --shrink-gamma 0.25 \
  --shrink-pert-key donor-cos \
  --shrink-pert-gamma 1.0 \
  --shrink-pert-floor 0.3 \
  --shrink-pert-missing damp \
  --effect-scale 2.4004
```

Three notes on that string, all of them places where the default is not what I solved for:

* `--shrink-gamma 0.25` against `AGREEMENT_SHRINK_GAMMA = 1.0`.
* `--shrink-pert-gamma 1.0` against `PERT_SHRINK_GAMMA = 2.0`.
* `--shrink-pert-floor 0.3` and `--shrink-pert-missing damp` **are** the HEAD defaults, and are
  written out anyway so the arm is reproducible from the string alone if a constant moves again.

⚠️ `--shrink-key` also accepts **`agreement-pairwise`**, which appeared in the same commit. I have
**not** measured it; everything here is `agreement`. Do not substitute it.

⚠️ `--basis-prune essential` is carried from a12 on the assumption a12 used the pruned whitener.
If a12's actual invocation differs, take a12's emission flags verbatim and change only the six
shrinkage/scale flags — those six are what this document solves.
