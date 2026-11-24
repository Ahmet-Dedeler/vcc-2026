"""414 pre-upload check: where does a built artifact's BULK sit against each context's mse cliff?

Streams the artifact (h5py, ~10 s, < 3 GB RSS), reads A_c / B_c exactly as the scorer forms them (decomp414.decompose:
cell_eval2's `_row_weights` / `_across_pert_budget` on log1p(5e4 tot/sum) - b_ctrl), and prints rho_c = B_c / claim_c
(187's retrodicted claims) and the AMPLITUDE-ONLY raw-mse change vs 410b under 187's model at 410b's board X
(i.e. before any direction gain; the direction gain must exceed this to pay).

    uv run python experiments/414-bulk-reprice/cliffcheck414.py <artifact.h5ad> [more.h5ad ...]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from decomp414 import decompose  # noqa: E402
from analyse414 import CTX, joint187, mse_c  # noqa: E402

WORK = HERE / "work"


def main():
    dec = json.loads((WORK / "decomp414.json").read_text())
    an = json.loads((WORK / "analyse414.json").read_text())
    m = joint187.Model(joint187.F_FIT)
    claim = {c: m.t[c]["claim"] for c in CTX}
    Kc = {c: m.t[c]["Kc"] for c in CTX}
    sumD, Xb = joint187.sumD0, an["Xb"]["410b"]
    ref = {c: mse_c(dec["410b"][c]["A"], dec["410b"][c]["B"], Kc[c], claim[c], Xb, sumD) for c in CTX}
    cache, res = {}, {}
    for p in sys.argv[1:]:
        path = Path(p).resolve()
        print(f"[{path.name}]", flush=True)
        d = decompose(path, cache)
        per = {c: mse_c(d[c]["A"], d[c]["B"], Kc[c], claim[c], Xb, sumD) - ref[c] for c in CTX}
        res[path.name] = dict(per_context=d, rho={c: d[c]["B"] / claim[c] for c in CTX},
                              A_vs_410b={c: d[c]["A"] / dec["410b"][c]["A"] for c in CTX},
                              amp_only_d_raw_vs_410b=per, amp_only_d_raw_mean=sum(per.values()) / 3)
        print(f"  A_c/A_c(410b) {[round(d[c]['A'] / dec['410b'][c]['A'], 4) for c in CTX]}  rho_c "
              f"{[round(d[c]['B'] / claim[c], 4) for c in CTX]} (410b {[round(dec['410b'][c]['B'] / claim[c], 4) for c in CTX]})"
              f"  amplitude-only d raw mse vs 410b {sum(per.values()) / 3:+.5f} (per ctx "
              f"{[round(per[c], 5) for c in CTX]})", flush=True)
    out = WORK / "cliffcheck414.json"
    old = json.loads(out.read_text()) if out.exists() else {}
    old.update(res)
    out.write_text(json.dumps(old, indent=1, default=float))


if __name__ == "__main__":
    main()
