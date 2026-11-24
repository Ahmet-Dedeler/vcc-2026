"""490: per-block dose audit over every recorded rc audit dump. How often does the feasibility fallback deliver
less than designed (prior lambda < 1, lambda = 0 = unsteered, MDE lambda < 1)?"""
import json, sys
from collections import defaultdict
import numpy as np

FILES = {
    "429a (A/B/C, plain MDE)": "experiments/429-zs-production/work/rc_audit_429a.jsonl",
    "433a (A/B/C, coarse MDE)": "experiments/433-coarse-mde-arm/work/rc_audit_433a.jsonl",
    "432 b0 (D/E/F synth)": "experiments/432-def-rehearsal-new-knobs/work/b0/rc_audit_b0.jsonl",
    "432 b1 (D/E/F synth)": "experiments/432-def-rehearsal-new-knobs/work/b1/rc_audit_b1.jsonl",
    "432 b2 (D/E/F synth)": "experiments/432-def-rehearsal-new-knobs/work/b2/rc_audit_b2.jsonl",
    "432 ef_b1 (D/E/F synth)": "experiments/432-def-rehearsal-new-knobs/work/ef_b1/rc_audit_ef_b1.jsonl",
    "448 efr lvl1 (D/E/F synth, q1.75)": "experiments/448-def-rehearsal-q175/work/efr/rc_audit_lvl1.jsonl",
}
out = {}
for name, path in FILES.items():
    rows = [json.loads(l) for l in open(path)]
    by = defaultdict(list)
    for r in rows:
        by[r["context"]].append(r)
    rec = {}
    for ctx, rr in sorted(by.items()):
        pl = np.array([r.get("prior_lambda", np.nan) for r in rr], float)
        ml = np.array([r.get("mde_lambda", np.nan) for r in rr], float)
        uns = np.array([r.get("steer_unsteered", 0.0) for r in rr], float)
        kf = np.array([r.get("prior_shrink_k", np.nan) for r in rr], float)
        K = np.array([r.get("ballast_K", np.nan) for r in rr], float)
        rec[ctx] = dict(blocks=len(rr),
                        prior_lambda_eq0=int((pl == 0).sum()), prior_lambda_lt1=int((pl < 1).sum()),
                        prior_lambda_med=float(np.nanmedian(pl)),
                        mde_lambda_eq0=int((ml == 0).sum()) if np.isfinite(ml).any() else None,
                        mde_lambda_lt1=int((ml < 1).sum()) if np.isfinite(ml).any() else None,
                        mde_missing=int((~np.isfinite(ml)).sum()),
                        unsteered=int((uns == 1).sum()),
                        kf_med=float(np.nanmedian(kf)), ballast_K_med=float(np.nanmedian(K)) if np.isfinite(K).any() else None)
    out[name] = rec
    print(name)
    for ctx, v in rec.items():
        print("  ", ctx, v)
json.dump(out, open("experiments/490-de-emitter-review/work/audit_lambda.json", "w"), indent=1)
