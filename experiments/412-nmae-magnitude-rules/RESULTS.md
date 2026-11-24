# 412 — oracle-free fold-change magnitude rules barely move nmae; inflating our most confident genes HURTS. NEGATIVE

Status: NEGATIVE (2026-09-22 ~19:10Z, H1, the shipped 395 construction's delivered statistics; 58 gated targets;
base nmae 0.95210 reproduces 410's H1 base395 row). Question: with the DE channel free under ballast, is there an
oracle-free magnitude design worth building for nmae (388: oracle per-target scalar −0.023; global best −0.003)?

| rule (reads only our lfc / q / p) | Δ nmae |
|---|---:|
| global ×1.25 / ×1.5 / ×2 / ×3 | **−0.0041** / −0.0035 / +0.0086 / +0.0711 |
| top-50 genes by our p ×2 / ×3 / ×5 | +0.0033 / +0.0135 / +0.0496 |
| top-100 / 200 / 400 by p ×2 | +0.0045 / +0.0061 / +0.0085 |
| floor |lfc| ≥ 0.1 / 0.2 / 0.3 on our calls | −0.0030 / +0.0436 / +0.1411 |

Our most confident genes do not carry signs reliable enough to amplify; only a mild global ×1.25 helps (≈ +0.001
overall). nmae is a sign-information problem (343/187), not a magnitude one. Reproduce: `rules412.py` (seconds).
