"""194 — assemble the screen into the deliverable tables."""
import json, sys
import numpy as np

TAG = sys.argv[1] if len(sys.argv) > 1 else "all3"
TRUTH = sys.argv[2] if len(sys.argv) > 2 else "incumbent"
SRC = sys.argv[3] if len(sys.argv) > 3 else "work/scores.json"
rows = [r for r in json.load(open(SRC))
        if r.get("swap_tag") == TAG and r["truth"] == TRUTH]
folds = ["k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"]
folds = [f for f in folds if any(r["fold"] == f for r in rows)]
order, seen = [], set()
for r in rows:
    if r["variant"] not in seen:
        seen.add(r["variant"]); order.append(r["variant"])
D = {(r["fold"], r["variant"]): r for r in rows}
base = "v0"

def col(v, f, k):
    r = D.get((f, v))
    return r[k] if r else float("nan")

print(f"=== swap={TAG} truth={TRUTH} — delta vs v0, mean over {len(folds)} PLOLO folds ===")
print(f"{'variant':<17}" + "".join(f"{f[:9]:>11}" for f in folds)
      + f"{'MEAN dpds':>11}{'sd':>9}{'d cos':>9}{'d q':>9}{'d q_w':>9}{'d sse':>10}")
out = {}
for v in order:
    dp = [col(v, f, "pds") - col(base, f, "pds") for f in folds]
    dc = np.nanmean([col(v, f, "cos") - col(base, f, "cos") for f in folds])
    dq = np.nanmean([col(v, f, "a15_sgn") - col(base, f, "a15_sgn") for f in folds])
    dw = np.nanmean([col(v, f, "a15_sgn_w") - col(base, f, "a15_sgn_w") for f in folds])
    ds = np.nanmean([col(v, f, "sse") - col(base, f, "sse") for f in folds])
    m = float(np.nanmean(dp))
    out[v] = dict(dpds=m, sd=float(np.nanstd(dp, ddof=1)), per_fold=dp,
                  dcos=float(dc), dq=float(dq), dq_w=float(dw), dsse=float(ds))
    print(f"{v:<17}" + "".join(f"{x:+11.5f}" for x in dp)
          + f"{m:+11.5f}{np.nanstd(dp, ddof=1):9.5f}{dc:+9.5f}{dq:+9.5f}{dw:+9.5f}{ds:+10.6f}")
print("\nabsolute baseline v0:")
for f in folds:
    r = D.get((f, base))
    print(f"  {f:<15} pds {r['pds']:.5f}  cos {r['cos']:.5f}  q {r['a15_sgn']:.5f}"
          f"  q_w {r['a15_sgn_w']:.5f}  sse {r['sse']:.6f}  n_targets {r['n_targets']}")
json.dump(out, open(f"work/summary-{TAG}-{TRUTH}.json", "w"), indent=2)
