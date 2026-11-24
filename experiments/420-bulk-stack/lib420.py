"""420 — shared machinery for pricing bulk STACKS vs 417a's bulk (PP + 410 key clip(0.3,3) + sharpening q 1.25).

Everything is 418's (lib418: PP composition with the composition constants as knobs, the a23 read-back with the
board-emulated pds at lam 0.5, render-space A_c/B_c) plus:
  * `realise_k`   414's ballast emulation (`ballast414.realise`) with the builder's KMAX exposed. KMAX 16 is 414's
                  function bit-for-bit (asserted by callers); for KMAX 24 / 32 the carriers of a capped block are 376 / 368
                  cells, so their total is S_C * (400 - KMAX) / 384 (S_C = 410b's median carrier total at 384 cells).
                  Blocks below the cap are realised exactly at any KMAX (K = ceil(...) makes S * p_bulk >= S_C * p_de).
  * `realised_rows` the delivered (log-pseudobulk) rows of a bulk design after ballast realisation, per frame.
  * key artifacts on the real A/B/C panel (410 at clip (0.3,3)/(0.5,3)/(0.3,inf), 415's LODO-V1 PP at (0.3,3), and
    LODO-V1 at (0.5,3) re-clipped from 415's own key vector -- written to ./work as a candidate artifact).
Nothing here writes outside experiments/420-bulk-stack/work (bytecode writing is disabled before the foreign imports).
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True          # importing 414/418/415 modules must not drop __pycache__ into their dirs
import hashlib  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
from pathlib import Path  # noqa: E402

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
W418 = ROOT / "experiments/418-bulk-composition/work"
sys.path.insert(0, str(ROOT / "experiments/418-bulk-composition"))

import lib418  # noqa: E402  (imports 162/218/228/410/169/415/414 machinery)
from lib418 import (ALPHA, E062, FOLDS, GAMMA, LAM, PROD_OPS, SHARP_Q, UNION, UNION_G, W414, W415, Donors,  # noqa: E402,F401
                    Operator, agg, agg_parts, betas_for, build_fold, cd4_for, cd4_frozen_and_partial, compose_arm,
                    compose_pp, deflate, deliver, donor_rows, endpoints, frames, key_scale, lodo_bases, lodo_key,
                    pre_normrel, score, shared_energy_key, sharpen, solve_fraction, source_residual, union_rows, weights)
import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"             # APFS clones of 418's fold caches; nothing else is read or written through plolo
from ballast414 import LB, S_C, realise as realise414  # noqa: E402
from plolo import TS, profile_from_effect  # noqa: E402

MAIN = ("k562_gwps", "orion_hct116", "orion_hek293t")
CANCER = MAIN
QS = (1.25, 1.35, 1.5, 1.75, 2.0)
KMAXES = (16, 24, 32)

# ---------------------------------------------------------------------------------------------- compositions
# name -> (m_cd4, g_h1)  [normrel cap -1.5, deflation k 5, partial g 0.5 fixed: 418 found them at their bulk optimum]
COMPS = {
    "PP": (1.0, UNION_G),
    "h1x2": (1.0, 2 * UNION_G),
    "cd4x1.5": (1.5, UNION_G),
    "stack": (1.5, 2 * UNION_G),              # the mandate's (a): H1 x2 + CD4 x1.5
    "stack_h1x1.5": (1.5, 1.5 * UNION_G),     # stack neighbours: redundancy may move the optimum of either knob
    "stack_cd4x1.25": (1.25, 2 * UNION_G),
}
COMP_ENV = {
    "PP": "",
    "h1x2": "VCC_H1_UNION_G=0.56",
    "cd4x1.5": "VCC_CD4_TARGET_F=0.51",
    "stack": "VCC_H1_UNION_G=0.56 VCC_CD4_TARGET_F=0.51",
    "stack_h1x1.5": "VCC_H1_UNION_G=0.42 VCC_CD4_TARGET_F=0.51",
    "stack_cd4x1.25": "VCC_H1_UNION_G=0.56 VCC_CD4_TARGET_F=0.425",
}
# 418's names for the single-knob compositions (known-answer rows)
COMP418 = {"PP": "PP", "h1x2": "b:h1x2", "cd4x1.5": "a:cd4x1.5"}

# ---------------------------------------------------------------------------------------------- keys
# name -> (kind, lo, hi); kind "K410" = the fold's 410-formula key (composition-free), "V1" = 415's LODO-V1 key
# computed on the PP composition and then HELD FIXED across compositions (that is what the builder's pinned artifact
# does), "V1own" = the LODO-V1 key recomputed on the arm's own composition (would need a new artifact).
KEYS = {
    "K3": ("K410", 0.3, 3.0),
    "K5_3": ("K410", 0.5, 3.0),
    "V1": ("V1", 0.3, 3.0),
    "V1_5_3": ("V1", 0.5, 3.0),
}
KEY_ENV = {
    "K3": "VCC_PERT_SCALE=wshare_lo0.3_hi3",
    "K5_3": "VCC_PERT_SCALE=wshare_lo0.5_hi3",
    "V1": "VCC_PERT_SCALE=lodo415_lo0.3_hi3",
    "V1_5_3": "VCC_PERT_SCALE=<NEW: lodo415_lo0.5_hi3 = work/pert_scale_415pp_lo0.5_hi3.npz>",
    "Kinf": "VCC_PERT_SCALE=wshare_lo0.3_inf",
}
ABC_KEYS = {  # A/B/C panel artifacts: path, pinned sha256
    "K3": (ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi3.npz",
           "d7f76a22f62ad6991d7220ce8a65218db38b2d71dcfe3c3781ec6a8f70409e0e"),
    "K5_3": (ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.5_hi3.npz",
             "387d6f1a76b866a314c633e6f879329036bee5a1f7d5bcadcad1e86cd04f1574"),
    "Kinf": (ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi1e+09.npz",
             "94c048e99aecd40004e70f559fe26bdd3b8b9aef9cb8c8cd0efbf8a72d2fe69e"),
    "V1": (ROOT / "experiments/415-lodo-key/work/pert_scale_415pp_lo0.3_hi3.npz",
           "5f1ae21f530734735b23bc0a114bb7b84988f0b62221c6a11c7a64d8176547d9"),
}
V1_53_PATH = WORK / "pert_scale_415pp_lo0.5_hi3.npz"


def sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def reclip(key, carried, lo, hi):
    """410/415's weights(): w = clip(key / median(key over carried), lo, hi) on carried targets, 1 elsewhere."""
    key = np.asarray(key, float)
    carried = np.asarray(carried, bool)
    med = float(np.median(key[carried]))
    return np.where(carried, np.clip(key / med, lo, hi), 1.0)


def make_v1_53() -> tuple[Path, str]:
    """LODO-V1 (415 PP) at clip (0.5, 3): the SAME key vector as the shipped lodo415_lo0.3_hi3 artifact, re-clipped.
    Known answer: re-clipping that key at (0.3, 3) reproduces its shipped w bit-for-bit."""
    src, want = ABC_KEYS["V1"]
    assert sha256(src) == want, "415 PP artifact sha"
    z = np.load(src)
    w03 = reclip(z["key"], z["carried"], 0.3, 3.0)
    assert np.array_equal(w03, z["w"]), "re-clipping 415's key at (0.3,3) does not reproduce its w"
    w = reclip(z["key"], z["carried"], 0.5, 3.0)
    if not V1_53_PATH.exists():
        np.savez(V1_53_PATH, targets=z["targets"], w=w, key=z["key"], carried=z["carried"], n_donors=z["n_donors"])
        V1_53_PATH.with_suffix(".json").write_text(json.dumps(dict(
            source=str(src.relative_to(ROOT)), source_sha256=want, clip=[0.5, 3.0],
            note="415's LODO-V1 PP key vector re-clipped at (0.5, 3); w = max(w_lo0.3, 0.5). Candidate artifact only: "
                 "not wired into the builder (PERT_SCALE_ARTIFACTS needs a reviewed entry with this sha).",
            w_quantiles=np.quantile(w, [0, .1, .25, .5, .75, .9, 1]).tolist(), frac_floor=float(np.mean(w[z["carried"]] <= 0.5)),
        ), indent=1))
    z2 = np.load(V1_53_PATH)
    assert np.array_equal(z2["w"], w) and np.array_equal(z2["w"], np.maximum(z["w"], 0.5))
    return V1_53_PATH, sha256(V1_53_PATH)


def abc_weights(name: str, targets) -> np.ndarray:
    if name == "V1_5_3":
        path, _ = make_v1_53()
    else:
        path, want = ABC_KEYS[name]
        got = sha256(path)
        assert got == want, f"{name}: sha {got} != pinned {want}"
    z = np.load(path)
    kidx = {str(t): i for i, t in enumerate(z["targets"])}
    return np.asarray([float(z["w"][kidx[t]]) for t in targets])


# ---------------------------------------------------------------------------------------------- ballast
def realise_k(p_de, p_bulk, kmax=16):
    """414's `realise` with KMAX exposed (bit-identical to it at kmax 16, asserted by callers)."""
    live = p_de > 1e-9
    r = float(np.min(p_bulk[live] / p_de[live]))
    sc = S_C if kmax == 16 else S_C * (400 - kmax) / 384.0
    K = int(min(kmax, max(1, int(np.ceil(sc * (1.0 / max(r, 1e-6) - 1.0) / LB)))))
    S = sc + K * LB
    want = S * p_bulk
    floor = sc * p_de
    clipped = want < floor
    tot = np.where(clipped, floor, want)
    return tot / tot.sum(), r, K, float((floor - want)[clipped].sum() / S), int(clipped.sum())


def de_profiles(carrier_design):
    """The carriers' design profile per frame (a bulk-only arm keeps the shipped carriers: PP unkeyed, unsharpened)."""
    return {ctx: [profile_from_effect(carrier_design[i], b, gate, "minnorm") for i in range(carrier_design.shape[0])]
            for ctx, (b, gate) in frames().items()}


def design_and_realised(bulk, P_DE, q=None, scale=None, kmaxes=KMAXES, zero_cols=None):
    """Per frame: the designed delivered rows and the ballast-realised delivered rows at each KMAX, plus block stats.
    `scale` = {ctx: s} multiplies the (sharpened) bulk effect per context (VCC_BULK_CTX_SCALE; sharpening restores
    each row's own norm, so s commutes with it). `zero_cols` are zeroed in every returned row (the panel)."""
    out = {}
    for ctx, (b, gate) in frames().items():
        F = bulk if q is None else sharpen(bulk, q, gate, mode="row")
        if scale is not None:
            F = F * float(scale[ctx])
        n = F.shape[0]
        des = np.empty_like(F)
        rea = {k: np.empty_like(F) for k in kmaxes}
        st = {k: dict(K=np.empty(n, int), r=np.empty(n), clip_mass=np.empty(n), clip_genes=np.empty(n, int)) for k in kmaxes}
        for i in range(n):
            p_bk = profile_from_effect(F[i], b, gate, "minnorm")
            des[i] = np.log1p(TS * p_bk) - b
            for k in kmaxes:
                p_re, r, K, cm, cg = realise_k(P_DE[ctx][i], p_bk, k)
                rea[k][i] = np.log1p(TS * p_re) - b
                st[k]["K"][i], st[k]["r"][i], st[k]["clip_mass"][i], st[k]["clip_genes"][i] = K, r, cm, cg
        if zero_cols is not None:
            des[:, zero_cols] = 0.0
            for k in kmaxes:
                rea[k][:, zero_cols] = 0.0
        out[ctx] = dict(design=des, realised=rea, stats=st)
    return out


def energy(H, nonpanel):
    """Render-space A (mean row energy) and B (mean centred row energy) over non-panel genes (abc418's render_AB)."""
    H = H[:, nonpanel]
    Hc = H - H.mean(0, keepdims=True)
    return dict(A=float((H * H).sum(1).mean()), B=float((Hc * Hc).sum(1).mean()))


def rowcos(X, Y):
    return np.einsum("ij,ij->i", X, Y) / np.maximum(np.linalg.norm(X, axis=1) * np.linalg.norm(Y, axis=1), 1e-300)
