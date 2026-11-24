# 422: sign-only carriers (shipped magnitudes, gamma-0 signs). Negative on H1; closed, patch not applied

Date: 2026-09-23 00:57–01:17Z. H1 six-member harness (416's run416.py copy), 410b's shipped bulk, three seeds, vs 416's
same-seed base. Known answer: base_s0 == 416's base_s0 (max|d| ≤ 1e-12, PASS).

**Why:** 417c on the board (γ-0 carriers) gained fid +0.0012 and reach +0.0024 raw but lost nmae +0.0089, which we read
as a magnitude-profile effect. The question: can the signs alone keep the fid/reach gain without the nmae cost?

| carriers | Δnmae (H1, not trusted) | Δfid | Δreach | Δjac | fid+reach+jac at board slopes |
|---|---:|---:|---:|---:|---:|
| g000, full γ 0 (= 417c; 416, 5 seeds) | −0.00683 | **+0.00253** | −0.0021 | +0.00025 | **+0.00113** |
| **s000** (\|cap15\| × sign(g000); 17.6 % of entries, 6.4 % of mass flip) | −0.00107 | +0.00021 ± 0.00011 | **−0.0162 ± 0.0092** | +0.00001 | **−0.00293** |
| s030 (\|cap15\| × sign(g030); 6.6 % / 1.1 % flip) | −0.00014 | +0.00004 | +0.0006 | +0.00001 | +0.00014 |

**Reading:**
- The sign swap alone removes almost all of γ 0's fid gain and costs reach.
- g000's fid/reach gain is not a sign effect. It comes from the unwhitened MAGNITUDE structure: which genes carry large
  |LFC| sets the call set and the order the Wilcoxon reads.
- Our 417c story ("better signs, worse magnitudes") is therefore wrong in its first half. The unwhitened magnitudes buy
  fid/reach and cost nmae on A/B/C.
- **Sign-only carriers are closed.** `patch_carriersign422.py` and its test stay in the repo, unapplied.

**The DE channel after 416/417c/421/422:**
- Carrier whitening: the board prefers the shipped γ 0.65 on net.
- CD4 dose: shipped is optimal.
- Sign-only: negative.
- The one untested direction is a γ between 0 and 0.65 weighed on the board (416's γ 0.3: DE +0.0015 on H1). But H1
  mis-prices nmae for carrier changes, and the board lost 0.0029 overall to nmae at γ 0, so any intermediate γ must be
  priced as nmae ≈ −(0.0029 × dose share).
