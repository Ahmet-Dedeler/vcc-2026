"""468: which (context, target) blocks differ between 458g and 458d, and between 458g and 446b? (the pair's noise model)

458g = 458d + the prior gate. If the emitter's random draws do not depend on the block prior, every design-strong
target's 400 cells are byte-identical in the two packages, and exactly the weak third differs. The board difference
458g − 458d is then the sign change on weak targets, not an emission re-roll. If more blocks differ than the gate's
weak counts, the draws shifted, and 468's bins widen to 2 sd of the ordering noise (±0.009), as preregistered.
458g's weak targets carry 446b's prior (CLEAN4), so with prior-independent draws they should equal 446b's blocks. Then
the three board entries split the programme's effect exactly: 458g − 446b is the strong two thirds, 458d − 458g the weak
third.

Restores both packages (`scripts/restore_prediction_h5ad.py --raw`, ~12 GB each, into work/cmp/), then streams the CSR
block by block (memory-light) and compares data and indices. Deletes the restored files at the end.
usage: uv run python experiments/468-prior-gate/compare458g.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
S = ROOT / "submissions"
STEMS = {"446b": "446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5",
         "458d": "458d-cap15-partial50-arch1dprog-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5",
         "458g": "458g-cap15-partial50-arch1dprog-gate33-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5"}
CMP = HERE / "work/cmp"


def labels(f, name):
    g = f["obs"][name]
    if isinstance(g, h5py.Group):
        cats = g["categories"][:].astype(str)
        return cats[g["codes"][:]]
    return g[:].astype(str)


def main():
    CMP.mkdir(parents=True, exist_ok=True)
    paths = {}
    for k, stem in STEMS.items():
        out = CMP / f"{k}.h5ad"
        if not out.exists():
            r = subprocess.run([sys.executable, str(ROOT / "scripts/restore_prediction_h5ad.py"), str(S / f"{stem}.vcc"), str(out), "--raw"],
                               capture_output=True, text=True)
            if r.returncode != 0:
                raise SystemExit(f"restore {k} failed: {r.stderr[-400:]}")
        paths[k] = out
    fg = h5py.File(paths["458g"], "r")
    cg, tg = labels(fg, "context"), labels(fg, "target_gene")
    ipg = fg["X/indptr"][:]
    keys = np.char.add(np.char.add(cg, "|"), tg)
    starts = np.flatnonzero(np.r_[True, keys[1:] != keys[:-1]])
    ends = np.r_[starts[1:], len(keys)]
    same_as: dict[str, dict[str, set]] = {}
    for ref in ("458d", "446b"):
        fr = h5py.File(paths[ref], "r")
        assert np.array_equal(labels(fr, "context"), cg) and np.array_equal(labels(fr, "target_gene"), tg), f"row order differs vs {ref}"
        ipr = fr["X/indptr"][:]
        hit = {}
        for s_, e_ in zip(starts, ends):
            a0, a1, b0, b1 = int(ipg[s_]), int(ipg[e_]), int(ipr[s_]), int(ipr[e_])
            same = (a1 - a0 == b1 - b0 and np.array_equal(ipg[s_:e_ + 1] - a0, ipr[s_:e_ + 1] - b0)
                    and np.array_equal(fg["X/indices"][a0:a1], fr["X/indices"][b0:b1])
                    and np.array_equal(fg["X/data"][a0:a1], fr["X/data"][b0:b1]))
            hit.setdefault(str(cg[s_]), set())
            if same:
                hit[str(cg[s_])].add(str(tg[s_]))
        same_as[ref] = hit
        fr.close()
    meta = json.loads((S / f"{STEMS['458g']}.h5ad.meta.json").read_text())
    audit = meta["config"].get("de_prior_gate", {}).get("contexts", {})
    out = {}
    for c in sorted({str(x) for x in cg}):
        all_t = {str(t) for cc, t in zip(cg[starts], tg[starts]) if str(cc) == c}
        d, b = same_as["458d"][c], same_as["446b"][c]
        weak = audit.get(c, {}).get("weak_clean4"); strong = audit.get(c, {}).get("strong_programme")
        out[c] = {"blocks": len(all_t), "equal_458d": len(d), "equal_446b": len(b), "equal_both": len(d & b), "equal_neither": len(all_t - d - b),
                  "gate_weak_clean4": weak, "gate_strong_programme": strong,
                  # exact: every block equals one reference; the weak-only matches cannot exceed the gate's weak count
                  "exact_split": len(all_t - d - b) == 0 and weak is not None and len(b - d) <= weak <= len(b)}
        print(f"[compare458g] {c}: {len(all_t)} blocks; = 458d {len(d)}, = 446b {len(b)}, both {len(d & b)}, neither {len(all_t - d - b)}; "
              f"gate weak {weak} strong {strong}", flush=True)
    ok = all(v["equal_neither"] == 0 for v in out.values())
    print("[compare458g]", "EXACT SPLIT: every 458g block equals 458d's (strong) or 446b's (weak); the three entries decompose the programme"
          if ok else "some blocks equal neither reference: the draws shifted; 468's bins widen to +-0.009", flush=True)
    (HERE / "work/compare458g.json").write_text(json.dumps({"per_context": out, "exact_split": ok,
        "weak_equal_446b_only": {c: sorted(same_as["446b"][c] - same_as["458d"][c]) for c in out}}, indent=1))
    fg.close()
    for p in paths.values():
        p.unlink()


if __name__ == "__main__":
    main()
