"""620d — renormalised over shared genes + the Satija Flex probes (597's pseudobulk597.npz: K562, HAP1, HEK293).

Changes vs 620b: (1) every profile is converted back to proportions (expm1(b)/TS), restricted to the common gene axis,
renormalised to sum 1 and re-logged (log1p CP50k over the SHARED genes), so feature-space totals cannot leak in
(caution from [36a81a]); (2) three new Flex probes from another lab (Satija VIPerturb pilot): K562 (known answer;
the K562 pair is dropped from its bias fit), HEK293 (near-same-line to orion_hek293t: HEK293 vs 293T), HAP1 (no real
bank profile: jiang_HAP1's b_ctrl is a placeholder, so HAP1 is a NEGATIVE here). Same bias / scoring as 620b, plus the
truth-removed negatives of 620c, all in one table with the candidate rule (top >= 0.33 AND margin >= 0.15).
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
import importlib.util

here = Path(__file__).parent
spec = importlib.util.spec_from_file_location("m620b", here / "match620b.py")
m620b = importlib.util.module_from_spec(spec); spec.loader.exec_module(m620b)  # type: ignore[union-attr]
m620 = m620b.m620
ROOT = m620.ROOT
TS = 50000.0


def main() -> None:
    prof = {n: m620.load_bank(d) for n, (d, _, _) in m620.BANKS.items() if n not in m620b.DROP}
    genes = [l.strip() for l in (ROOT / "data/challenge/gene_names.csv").read_text().splitlines()[1:] if l.strip()]
    for c, p in m620.CHALLENGE.items():
        prof[f"challenge_{c}"] = dict(zip(genes, map(float, np.load(p))))
    z = np.load(ROOT / "experiments/597-flex-hek293-floor/work/pseudobulk597.npz", allow_pickle=True)
    lines = [str(x) for x in z["lines"]]
    g597 = [str(x) for x in z["genes"]]
    for i, ln in enumerate(lines):
        cnt = z["counts"][i].astype(float)
        prof[f"satija_{ln}"] = dict(zip(g597, np.log1p(TS * cnt / cnt.sum())))
    group = {n: m620.BANKS[n][1] for n in prof if n in m620.BANKS}
    group.update({f"challenge_{c}": m620.CHALLENGE_TRUTH[c] for c in "ABC"})
    group.update({"satija_K562": "K562", "satija_HEK293": "HEK293", "satija_HAP1": "HAP1(no real bank)"})
    flex = set(m620b.FLEX) | {"satija_K562", "satija_HEK293", "satija_HAP1"}
    names = list(prof)
    cnt = Counter(g for p in prof.values() for g in p)
    axis = sorted(g for g, k in cnt.items() if k >= 0.85 * len(names))
    raw = np.array([[prof[n].get(g, np.nan) for g in axis] for n in names])
    # renormalise over the shared axis (NaN genes excluded from each profile's own total)
    P = np.expm1(raw) / TS
    P = P / np.nansum(P, axis=1, keepdims=True)
    M = np.log1p(TS * P)
    idx = {n: i for i, n in enumerate(names)}

    def rc(a, b):
        ok = np.isfinite(a) & np.isfinite(b)
        return float(np.corrcoef(np.argsort(np.argsort(a[ok])), np.argsort(np.argsort(b[ok])))[0, 1])

    def score(probe: str, drop_group: str | None) -> list[tuple[float, str]]:
        g0 = group[probe]
        use = [k for k, (f, o) in m620b.PAIRS.items() if f != probe and group.get(o) != g0
               and not (probe.startswith("satija_") and f == "viperturb_k562" and g0 == "K562")]
        bias = np.nanmean([M[idx[m620b.PAIRS[k][0]]] - M[idx[m620b.PAIRS[k][1]]] for k in use], axis=0)
        refs = [n for n in names if n not in flex and (drop_group is None or group.get(n) != drop_group)]
        R = M[[idx[r] for r in refs]]
        mu = np.nanmean(R, axis=0); sd = np.nanstd(R, axis=0)
        keep = (mu > 1.0) & (sd > 0.1) & np.isfinite(bias) & np.isfinite(M[idx[probe]])
        zx = ((M[idx[probe]] - bias)[keep] - mu[keep]) / sd[keep]
        zr = (R[:, keep] - mu[keep]) / sd[keep]
        return sorted(((rc(zx, zr[i]), r) for i, r in enumerate(refs)), reverse=True)

    rows = []
    probes = ["challenge_A", "h1_2025", "viperturb_k562", "satija_K562", "satija_HEK293", "satija_HAP1",
              "challenge_B", "challenge_C", "cd4_rest"]
    for probe in probes:
        for mode in ("full", "truth_removed"):
            g0 = group[probe]
            has_bank = any(group.get(n) == g0 for n in names if n not in flex)
            if mode == "truth_removed" and not has_bank:
                continue
            sc = score(probe, g0 if mode == "truth_removed" else None)
            top, r1 = sc[0]
            nxt = next(s for s, r in sc[1:] if group[r] != group[r1])
            margin = top - nxt
            fires = top >= 0.33 and margin >= 0.15
            correct = group[r1] == g0 or (g0 == "HEK293" and r1 == "orion_hek293t")
            kind = "positive" if (mode == "full" and has_bank) else "negative"
            rows.append({"probe": probe, "mode": mode, "kind": kind, "top": r1, "top_group": group[r1],
                         "score": round(top, 4), "margin": round(margin, 4), "rule_fires": fires,
                         "fires_correctly": bool(fires and correct)})
            print(f"{probe:15s} {mode:13s} {kind:8s} top {r1:16s}({group[r1]:10s}) {top:+.3f} margin {margin:+.3f} "
                  f"rule {'FIRES' if fires else 'quiet'}{' (correct)' if fires and correct else (' (WRONG)' if fires else '')}")
    (here / "work/match620d.json").write_text(json.dumps({"axis": len(axis), "rows": rows}, indent=1))


if __name__ == "__main__":
    main()
