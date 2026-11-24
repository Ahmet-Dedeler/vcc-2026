# 419: carriers at negative whitening gamma. Superseded before any arm ran

Status: **CLOSED, superseded** (2026-09-22 23:25Z).

- 416 (worker) ran the γ −0.3 carriers (anti-whitening) on H1 inside its own menu. The DE sum at board slopes was
  **−0.0073** (reach −0.047), against **+0.0030** for γ 0. The γ dose response peaks at 0 and turns sharply negative
  below it, so γ −0.65 is not worth an H1 slot.
- 419's own driver aborted at its known-answer row: base_s0 vs 416's base_s0 differed by 1.1e-16 in nmae (float
  noise), and the check required exact equality. The check was too strict. 416 uses |Δ| ≤ 1.1e-16. No arm ran.
- The designs (`work/carriers419.npz`, `carriers419.json`: g065 == cap15 KA, g000 == 416's g000 sha) stay for reference.

The next carrier question (416's suggestion) is the CD4 dose on the carriers. 416's 2×2 found CD4 worth about +0.006 of
DE on its own. That question is experiment 421.
