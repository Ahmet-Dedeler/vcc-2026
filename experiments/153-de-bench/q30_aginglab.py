"""Q30 — is Aginglab's reach jump headroom WE still have, or headroom we already spent?

We withdrew the emitter as the explanation for their v9->v11 (mse 1.82->0.96, reach
0.167->0.267) because our own residual artefact is only ~19% -- too small to yield +0.100.

But that compares their starting point to OURS.  Their v9 raw mse of 1.8219 is the signature of
a badly under-dispersed emitter forfeiting the #247 credit; ours at a15 is 0.8828 because a12/a15
already spent that path.  If the emitter explains them, their available artefact reduction was
the full ~74% -> ~19%, not our remaining ~19% -> 0%.

The decisive question is therefore not "what did Aginglab do" but: **has any team gained large
reach while ALREADY at good mse?**  If none has, the emitter route is spent for us and their
jump is headroom we have taken.  If one has, something else is available.
"""
import json, glob
from collections import defaultdict

hist = defaultdict(dict)
for f in sorted(glob.glob("../../data/leaderboard/*.json")):
    d = json.load(open(f))
    for r in (d if isinstance(d, list) else d.get("entries") or []):
        if isinstance(r, dict) and r.get("score_avg") is not None:
            hist[r["team_name"]][r["id"]] = r

pairs = []
for t, es in hist.items():
    seq = sorted(es.values(), key=lambda r: r.get("submission_date") or "")
    for a, b in zip(seq, seq[1:]):
        if a["id"] != b["id"]:
            pairs.append((t, a, b))

M = "expr_mse_unbiased_capped_norm"; R = "de_wilcoxon_direction_reach_raw"
big = [(t, a, b) for t, a, b in pairs if b[R] - a[R] > 0.04]
print(f"{len(pairs)} within-team pairs; {len(big)} with reach gain > +0.04\n")
print(f"{'team':22} {'reach before->after':>22} {'mse before->after':>22} {'pds Δ':>8} {'avg Δ':>8}")
for t, a, b in sorted(big, key=lambda x: -(x[2][R] - x[1][R]))[:14]:
    print(f"{t[:22]:22} {a[R]:9.4f} ->{b[R]:9.4f} {a[M]:9.3f} ->{b[M]:9.3f} "
          f"{b['pds_cosine']-a['pds_cosine']:+8.4f} {b['score_avg']-a['score_avg']:+8.4f}")

print("\nSPLIT BY STARTING mse — was the team already dispersion-tuned when it gained reach?")
already = [(t, a, b) for t, a, b in big if a[M] < 1.20]
under = [(t, a, b) for t, a, b in big if a[M] >= 1.20]
print(f"  started UNDER-dispersed (raw mse >= 1.20): {len(under)} pairs")
print(f"  started already tuned    (raw mse <  1.20): {len(already)} pairs")
for t, a, b in already:
    print(f"      {t[:24]:24} mse {a[M]:.3f}->{b[M]:.3f}  reach {a[R]:.4f}->{b[R]:.4f}  "
          f"avg {a['score_avg']:.4f}->{b['score_avg']:.4f}")
print("\nOUR OWN position for reference: a15 raw mse 0.8828, raw reach 0.1806.")
