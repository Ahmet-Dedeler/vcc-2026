# 454 — does a same-line bank know each target's effect SIZE? No: keep our magnitudes (2026-09-26 00:18Z, analytic, H1)

**Question.** The same-line mode (451) blends only direction. A covered target keeps our row norm: the cap15 composition
norm times the 410 key weight. Would the bank's own per-target effect size be a better magnitude for matched contexts
in October?

**Setup** (`mag454.py`, 2 s, 12 MB):
- h1-test truth: 323's `reference.npz`, 85 targets, bulk effect = log1p(5e4 × pseudobulk CPM) minus the control's, on
  17,992 shared non-panel genes.
- Pan 2025 hESC: 30 of the 85 covered; one pseudobulk row per target, 16,021 measured genes each.
- Our bulk rows: the harness's own (435: `native_cap15`, then `ps_bulk`, i.e. 410b's key with clip(0.3, 3)).

| predictor of the truth's per-target energy (30 covered targets) | Spearman | Pearson (log) |
|---|---:|---:|
| **ours (cap15 × 410 key), norm²** | **+0.552** | +0.548 |
| native cap15 without the key | +0.354 | +0.406 |
| **Pan's row norm²** | **+0.056** | +0.143 |
| geometric mean of ours and Pan | +0.522 | +0.577 |

Jackknife-debiased truth energy gives the same ordering (+0.552 / +0.354 / +0.056 / +0.522).
Over all 85 targets, ours reads +0.49.

**Direction vs size.**
- Pan's rows point far better than ours: median cosine to the truth 0.351 vs 0.048.
- Their norms are noise-dominated. TOR1A is null in the truth (energy 16.5 against a jackknife floor of 13.7), yet its
  Pan norm² is 22.6. MED19 has the largest true effect (324), yet its Pan norm² is only 34.9.
- Pan's direction quality rises with the true effect size: Spearman(true energy, cos(Pan, truth)) = +0.68.

**Verdict: NO magnitude lever from a pseudobulk bank.** The 410 key (cross-donor shared energy) is the better
per-target size estimate, even against a same-line bank with 0.35 cosine. 451's design (direction blend at our own norm)
is right. A same-line size estimate would need per-guide or split-half rows (a shared-energy key inside the line). Our
registry banks store one mean row per target, so that is closed until a bank ships per-guide pseudobulks.

## 454b — a strength-dependent blend weight? No (`wstrength454.py`, 00:20Z)
Setup: cos(unit(ours) + w·unit(Pan), truth) per covered target (bulk rows, non-panel genes), and the rule
w_t = w0·(s_t / median s)^g, with s = our own key-scaled norm (non-oracle).

| global w | 0 | 0.5 | 1 | 2 | 4 | 8 | 32 |
|---|---:|---:|---:|---:|---:|---:|---:|
| mean cos to truth | 0.068 | 0.197 | 0.264 | 0.307 | 0.321 | 0.323 | 0.322 |

- Both halves by our strength prefer w 4–8 (weaker half 0.286 at w 2 → 0.306 at w 8; stronger half 0.328 → 0.341).
- The per-target best w is weakly anti-correlated with our strength (Spearman −0.25).
- Every g > 0 lowers the mean cosine (w0 2: g 0 0.307, g 0.5 0.297, g 1 0.283).
- The flat w stands. The emitter-level dose (450d: w 2–4 peak, flat from 1 to 4) is the number to use. Direction
  cosine alone would favour a slightly higher w (4–8) at Pan quality. That is consistent with 450's "near-replacement
  is only 0.003 lower", and not a reason to move the October default of 2.
