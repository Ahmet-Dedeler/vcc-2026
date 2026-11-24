"""414: emulate 410d and 410d2 on the REAL 410c bulk before a 40-minute build.

Stand-in for the 410c effect matrix: its rendered bulk d (as in sharp410c.py; on gated genes the render is ~identity).
  410d  emulation: d * m * 0.532508 (the constant the real 410d build logged), re-rendered per context.
        KNOWN ANSWER: must reproduce the measured 410d artifact energy ratios vs 410c (x1.197 / x1.182 / x2.127) and
        SYVN1-C's blow-up (HSP90B1 = 50 % of the block's counts).
  410d2 emulation: d * m * c_rg, c_rg = the constant that preserves the energy on each context's RENDERED non-panel genes
        summed over A/B/C (patch_gain410d2.py's rule), re-rendered.
  capped variants: m clipped to [1/cap, cap] before the same render-gated normalisation (a suggested fix).
Reports the per-context render-space energy ratio vs the 410c re-render (410d2's own gate: 0.97...1.03), the max per-row
ratio and the largest single-gene mass share.

    uv run python emul410d2.py
"""
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
os.environ.setdefault("NUMBA_CACHE_DIR", str(HERE / "work" / "numba-cache"))
import numpy as np  # noqa: E402

sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "experiments/162-retrodicting-validator"))
from shared414 import rows_of  # noqa: E402
from plolo import TS, MIN_CPM, profile_from_effect  # noqa: E402

ROOT = HERE.parents[1]
ART = ROOT / "experiments/410-pertscale-revisit/work/build_c/the-kettle-applied-for-tenure.h5ad"
M = np.asarray(np.load(ROOT / "experiments/410-pertscale-revisit/work/bulk_gene_gain_g228reg_t025.npz")["m"], float)
MEASURED_410D = {"A": 1.1965, "B": 1.1818, "C": 2.1274}      # cliffcheck414 on the built artifact


def render(F, b, gate):
    out, share = np.empty_like(F), np.empty(F.shape[0])
    for i in range(F.shape[0]):
        p = profile_from_effect(F[i], b, gate, "minnorm")
        out[i] = np.log1p(TS * p) - b
        share[i] = p.max()
    return out, share


def main():
    cache = {}
    genes, R = rows_of(ART, cache)
    targets = sorted(R["A"])
    gpos = {g: i for i, g in enumerate(genes)}
    panel = np.zeros(genes.size, bool)
    panel[[gpos[t] for t in targets]] = True
    D, G = {}, {}
    for ctx in "ABC":
        b = cache["b"][ctx]
        cpm = 1e6 * np.expm1(b) / np.expm1(b).sum()
        G[ctx] = (cpm > MIN_CPM) & (b > 0)
        T = np.vstack([R[ctx][t] for t in targets])
        d = np.log1p(TS * T / T.sum(1, keepdims=True)) - b[None, :]
        d[:, panel] = 0.0
        D[ctx] = np.where(G[ctx][None, :], d, 0.0)
    base = {ctx: render(D[ctx], cache["b"][ctx], G[ctx]) for ctx in "ABC"}
    e0 = {ctx: np.einsum("ij,ij->i", base[ctx][0][:, ~panel], base[ctx][0][:, ~panel]) for ctx in "ABC"}

    def rg_const(mm):
        before = sum((D[c][:, G[c] & ~panel] ** 2).sum() for c in "ABC")
        after = sum(((D[c] * mm[None, :])[:, G[c] & ~panel] ** 2).sum() for c in "ABC")
        return float(np.sqrt(before / after))

    variants = {"410d_emul (m x 0.532508)": (M, 0.532508), "410d2_emul (render-gated)": (M, rg_const(M))}
    for cap in (3.0, 2.0, 1.5):
        mc = np.clip(M, 1 / cap, cap)
        variants[f"410d2 + m clipped to [1/{cap:g},{cap:g}]"] = (mc, rg_const(mc))
    out = {}
    syv = targets.index("SYVN1")
    for name, (mm, c) in variants.items():
        res = dict(constant=c)
        for ctx in "ABC":
            dl, share = render(D[ctx] * mm[None, :] * c, cache["b"][ctx], G[ctx])
            e = np.einsum("ij,ij->i", dl[:, ~panel], dl[:, ~panel])
            res[ctx] = dict(energy_ratio=float(e.sum() / e0[ctx].sum()), max_row_ratio=float((e / e0[ctx]).max()),
                            argmax_target=targets[int((e / e0[ctx]).argmax())], max_mass_share=float(share.max()),
                            syvn1_ratio=float(e[syv] / e0[ctx][syv]))
        out[name] = res
        print(f"{name:34s} c {c:.4f}  " + "  ".join(
            f"{ctx}: E x{res[ctx]['energy_ratio']:.3f} (row max x{res[ctx]['max_row_ratio']:.1f} {res[ctx]['argmax_target']}, "
            f"SYVN1 x{res[ctx]['syvn1_ratio']:.1f}, mass share max {res[ctx]['max_mass_share']:.3f})" for ctx in "ABC"), flush=True)
    ka = {c: out["410d_emul (m x 0.532508)"][c]["energy_ratio"] / MEASURED_410D[c] for c in "ABC"}
    print(f"KA: 410d emulation / measured 410d artifact energy ratio: {({c: round(v, 3) for c, v in ka.items()})}")
    out["ka_emul_over_measured"] = ka
    (HERE / "work" / "emul410d2.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
