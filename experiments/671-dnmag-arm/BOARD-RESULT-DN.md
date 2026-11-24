# DN board result (10-15 s1; read by the coordinator [400b6e] at 01:50Z 10-15)

Entry **wzaVLFTmlAMW3Mfw0dD8** (package 3c6a51b8…, submitted 01:12Z, upload 20 min, published ~01:49Z): **overall 0.2741292, rank 62**,
vs its base 630arm (QJRJK6OC6Zw1kjwnZqid, 0.2737556): **+0.0003736**. Prereg: experiments/631-nmae-headroom/PREREGISTRATION-DNMAG.md (frozen 10-02).

| member | raw Δ vs 630arm | scaled Δ |
|---|---:|---:|
| pds | 0 | 0 |
| mse | 0 | 0 |
| **nmae** | **−0.001646** | +0.002241 |
| fid | 0 | 0 |
| reach | +9.6e-8 | 0 |
| jac | 0 | 0 |

- **Falsifiers: all clean** (|Δfid|, |Δjac|, |Δpds|, |Δmse| = 0 ≤ 1e-6; |Δreach| 9.6e-8 ≤ 0.0014). Build fail-safes: dn_invariant_failed and
  dn_guard_failed max 0 in every block (also ur_*); config diff vs 630arm exactly ['de_dnmag']; zero-gene audit 0 in 900 blocks.
- **Bin: PARTIAL** (Δnmae raw −0.0025 … −0.0008). The prereg's own expected bin. Prediction −0.0018 … −0.0048 (centre −0.003): the read
  sits just below the range, at 0.55× the centre, i.e. A/B/C delivered about half of H1's f 0.35 price.
- **Decision (prereg §a PARTIAL):** adopt for A/B/C arms. October adopts **d_m200, flagged**, conditional on the build's capacity floor
  holding on D/E/F's own blocks (mean dn_frac ≥ 0.10 and dn_moved ≥ 200 genes/target per context); the day-of check is a postbuild
  report ([17adab]), and a context that fails it is recorded, not rebuilt (DN's harm is bounded by its own fail-safes).
- **Cost:** DN's pass doubled the build (≥ 139 min in the build stage vs 630arm's 74). [17adab]'s byte-identical v2 of dn_mag
  (experiments/693-dnmag-speed) goes in before the freeze, with its KA.
- 704 (α ×1.25 contingency) is NOT launched: it was for NULL/HARM only.
- Try count: DOWN-magnitude forms on the board 1 (this); H1 m = 2 only.
