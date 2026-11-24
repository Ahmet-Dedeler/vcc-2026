#!/usr/bin/env python3
"""616 S3: Δ(clean − thinned) for every instrument output, both from sandbox runs of the SAME code.
  thinned = the identity-mode run (tag `ka`; equals the committed output wherever ka616 says byte-identical);
  clean   = the clean-mode run (tag `clean`).
JSON outputs: every numeric leaf, flattened to a key path, with thinned / clean / Δ and a sign-flip flag.
Text (.out) outputs: every line whose numbers changed, thinned line above clean line.
    uv run python compare616.py  ->  work/compare616.json, work/compare616.md"""
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
SB = HERE / "work/sandbox"
FILES = ["experiments/639-identity-axis-lolo/work/decomp639.json",
         "experiments/613-block-prior-validation/work/block613.json",
         "experiments/613-block-prior-validation/work/block613d.json",
         "experiments/612-key-gated-hook/work/ds612b.json", "experiments/612-key-gated-hook/work/ds612c.json",
         "experiments/633-updose-folds/work/dose633.json",
         "experiments/613-block-prior-validation/work/block613.out",
         "experiments/613-block-prior-validation/work/block613b.out",
         "experiments/613-block-prior-validation/work/block613c.out",
         "experiments/613-block-prior-validation/work/block613d.out",
         "experiments/613-block-prior-validation/work/block613d3.out",
         "experiments/613-block-prior-validation/work/block613d4.out",
         "experiments/613-block-prior-validation/work/block613d5.out",
         "experiments/613-block-prior-validation/work/block613d6.out",
         "experiments/613-block-prior-validation/work/block613e.out",
         "experiments/613-block-prior-validation/work/block613e2.out",
         "experiments/633-updose-folds/work/dose633.out"]
NUM = re.compile(r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+(?:[eE][-+]?\d+)?")


def leaves(x, path=""):
    if isinstance(x, dict):
        for k, v in x.items():
            yield from leaves(v, f"{path}/{k}")
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from leaves(v, f"{path}[{i}]")
    elif isinstance(x, (int, float)) and not isinstance(x, bool):
        yield path, float(x)
    elif isinstance(x, bool):
        yield path, float(x)


def main():
    out, md = {}, ["# 616 — Δ(clean − thinned), every instrument output (thinned = identity sandbox run, same code)\n"]
    for rel in FILES:
        a, b = SB / "ka" / rel, SB / "clean" / rel
        if not (a.exists() and b.exists()):
            out[rel] = "missing"
            md.append(f"## {rel}\nmissing (thinned {a.exists()}, clean {b.exists()})\n")
            continue
        md.append(f"## {rel}")
        if rel.endswith(".json"):
            la, lb = dict(leaves(json.loads(a.read_text()))), dict(leaves(json.loads(b.read_text())))
            rows = []
            for k in la:
                if k in lb and la[k] != lb[k]:
                    flip = (la[k] > 0) != (lb[k] > 0) and la[k] != 0 and lb[k] != 0
                    rows.append({"key": k, "thinned": la[k], "clean": lb[k], "delta": lb[k] - la[k], "sign_flip": flip})
            out[rel] = rows
            md.append(f"{len(rows)} numeric leaves changed; sign flips: {sum(r['sign_flip'] for r in rows)}\n")
            md += [f"- `{r['key']}`: {r['thinned']:+.5g} → {r['clean']:+.5g} (Δ {r['delta']:+.5g}){'  **SIGN FLIP**' if r['sign_flip'] else ''}"
                   for r in rows[:400]]
        else:
            ta, tb = a.read_text().splitlines(), b.read_text().splitlines()
            rows = []
            for i in range(max(len(ta), len(tb))):
                x = ta[i] if i < len(ta) else ""
                y = tb[i] if i < len(tb) else ""
                if NUM.findall(x) != NUM.findall(y) or x != y:
                    rows.append({"line": i + 1, "thinned": x, "clean": y})
            out[rel] = rows
            md.append(f"{len(rows)} lines changed\n")
            for r in rows[:300]:
                md.append(f"- L{r['line']} thinned: `{r['thinned'][:220]}`")
                md.append(f"  L{r['line']} clean:   `{r['clean'][:220]}`")
        md.append("")
    # 609's statistic, directly on the truths: the DOWN share of DE calls (q < 0.05) per Orion gate table
    import numpy as np
    G480 = HERE.parents[0] / "480-nmae-instrument/work"
    md.append("## DE-call DOWN share of the gate tables (q < 0.05; 609 saw 72–74 % DOWN on Orion block gate genes)")
    out["down_share"] = {}
    for name in ("orion_hek293t", "orion_hct116"):
        r = {}
        for lab, f in (("thinned", G480 / f"gate-{name}.npz"), ("clean", HERE / f"work/gate-{name}_clean.npz")):
            if f.exists():
                z = np.load(f, allow_pickle=True)
                q, lfc = np.asarray(z["q"], float), np.asarray(z["lfc"], float)
                sig = np.isfinite(q) & (q < 0.05) & np.isfinite(lfc) & (lfc != 0)
                r[lab] = {"calls": int(sig.sum()), "down_share": float((lfc[sig] < 0).mean()),
                          "median_calls_per_target": float(np.median(sig.sum(1)))}
        out["down_share"][name] = r
        md.append(f"- {name}: " + "; ".join(f"{k} {v['down_share']:.3f} of {v['calls']:,} calls "
                                           f"(median {v['median_calls_per_target']:.0f}/target)" for k, v in r.items()))
    md.append("")
    (HERE / "work/compare616.json").write_text(json.dumps(out, indent=1))
    (HERE / "work/compare616.md").write_text("\n".join(md) + "\n")
    print("\n".join(md[:80]))


if __name__ == "__main__":
    main()
