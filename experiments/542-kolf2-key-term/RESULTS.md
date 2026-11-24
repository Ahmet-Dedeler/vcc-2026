# 542 — kolf2 (iPSC) as an extra truth proxy in the LODO key HURTS: NEGATIVE

2026-09-27 04:55Z, third session. Exploratory, `arms542.py` (PLOLO, lib420, 541's path, the three cancer folds only;
kolf2 cannot be held out here because it is the proxy). Peak RSS 5.5 GiB per fold.

## Question
The oracle per-target amplitude has +22 … +26 % direction headroom over the LODO key (540's side check). The key
averages only 2 leave-one-out terms per target on the cancer folds (3 in production). The essential-gene screens
(RPE1, Nadig Jurkat, K562 essential) cover 0 of the 300 panel targets and h1_2025 covers 20, so kolf2 (273 targets)
is the only broad external proxy.
- Term: t_p = ⟨F_p, T_kolf2,p⟩ / |F_p|² (F = the PP composition row with panel columns zeroed, the LODO key's
  denominator; T = kolf2's cached fold truth).
- It is rescaled to the donor terms' scale (the ratio of medians: 2.12–2.38) and averaged with them:
  (Σ_j v_j + t′) / (n_j + 1).

## Result (relative cos, cancer folds k562 / hct116 / hek293t; KA: the mean key reproduces comp420 and 540b rows, |d| 0)
| key | vs the LODO key at the same clip | vs 446bj2's key (0.5, 3) |
|---|---|---|
| LODO + kolf2, clip (0.5, 3) | −1.09 % (1/3): +0.30 / −2.09 / −1.48 | −1.09 % |
| LODO + kolf2, clip (0.5, 6) | −2.39 % (0/3): −0.74 / −4.36 / −2.08 | +0.15 % (cap 6 alone: +2.62 %) |
| kolf2 term alone, (0.5, 3) | −6.85 % (0/3) | — |
| kolf2 term alone, (0.5, 6) | −8.82 % (0/3) | — |

- The kolf2 term agrees only loosely with the LODO key (Spearman +0.44 … +0.51 over 240 targets).
- As a key on its own it loses 7–9 %. Mixed in, it dilutes the cancer donors' terms and gives back most of cap 6's
  gain.
- **Mechanism:** the key estimates how a held-out cancer line's response projects on the prediction. The other cancer
  lines are good proxies for that; an iPSC line's per-target amplitudes are not.
- More proxies help only if they respond like the target contexts. For A/B/C (Jurkat, HeLa, CAL33) and D/E/F, that
  means cancer lines. The panel-covering cancer datasets (K562, HCT116, HEK293T) are already the donors.

**Decision:** no external-proxy term. The key stays the cancer donors' LODO mean at cap 6 (540 / 495).
