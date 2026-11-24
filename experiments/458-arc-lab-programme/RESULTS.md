# 458 — Arc's own H1 programme as every context's block prior (build log; board read on 09-30)

**Built 07:55Z 09-26** by `chain458.sh` through `run_arm_generic.sh`: 446b + `VCC_SAMELINE=<arch1_programme_abc.npz 32b82850…>:0`.
- **[458] PASS:** the config diff vs 446b is exactly `['sameline']`.
- CLIFF OK with render energies 0.9603 / 1.0282 / 1.0785, 446b's own: the bulk is untouched.
- Zero-gene audit 0 (17 cells below −5). Novelty rc 0.
- Package `submissions/458-cap15-partial50-arch1prog-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.vcc`,
  sha256 `e3fabbfa…`.
- Meta `config.sameline.contexts`: A, B and C each have `programme_default` true, down fraction 0.578 (Arc-H1), w 0,
  0 rows changed.
- On each block the prior's sign differs from CLEAN4 on ~40 % of genes (PREREGISTRATION addendum).

Armed for 09-30 slot 1 by `459-h1-axis-def/chain0930.sh`. It is read against 446b with `PREREGISTRATION.md`'s bins.

## 458d build (2026-09-26 11:25Z → 11:58Z, `chain458d.sh`)
`run_arm_generic.sh 458d`: 446b + `VCC_SAMELINE=arch1d_programme_abc.npz:6951dbc7…:0`.
- **[458d] PASS: config diff exactly ['sameline'].**
- Zero-gene audit 0 in 0 blocks (lfc < −5: 17; > +3: 3, max +4.11). **CLIFF OK** A 0.9603 / B 1.0282 / C 1.0785,
  identical to 458's (the bulk is untouched). Novelty rc 0.
- Package `submissions/458d-cap15-partial50-arch1dprog-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.vcc`,
  sha256 `a4201ffe…`.
- `459/chain0930b.sh` arms it as 09-30 s1 ("the pipette learned to meditate"). Reader: `read458.py prog <entry>
  UDgjOKL9gI1udAexq6ue`.

## 458e build (2026-09-26 13:21Z → 13:54Z, `chain458e.sh`)
446b + `VCC_SAMELINE=arch1e_blend_abc.npz:b67af1f5…:0` (the 0.5 CLEAN4 + Arc-full blend).
- **[458e] PASS: config diff exactly ['sameline'].** Zero-gene audit 0. CLIFF OK A 0.9603 / B 1.0282 / C 1.0785
  (bulk identical). Novelty rc 0.
- Package `submissions/458e-cap15-partial50-arch1eblend-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.vcc`,
  sha256 `8a755fd0…`.
- It is a 10-01 candidate if 458d reads positive or partial (`464/chain1001.sh`).
