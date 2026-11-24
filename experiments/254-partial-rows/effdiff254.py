"""254 falsifiers 2 and 3 at the EFFECT-MATRIX level (249 section 6a's exact form), split by PARTIAL row.

`effects-composed.npy` vs `effects-k5.npy` (the two staged 300 x 18,533 matrices, the same builder,
one flag apart: the arm vs cd4svd5).  The 15 partial rows must move in DIRECTION; the other 285 must
equal cd4svd5's times ONE shared scalar (the merged preserve ratio) to machine precision.  Also
vs `effects-off.npy` (h1union): 291 rows move (276 CD4-carried by the deflation, 15 partial), the 9
rows CD4 never touches are one scalar.  Non-vacuity: `||dE||/||E||`, zero-displacement rows.
KNOWN ANSWERS: `effects-off` hashes to abeb0749..., `effects-k5` to a28d5be2....

    uv run python effdiff254.py
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
S = WORK / "stage"
WANT_OFF = "abeb0749f8ecc9410fcec0d9d0cebd2569c97ca456d73f9e6b3bfd359164c315"
WANT_K5 = "a28d5be2e9cb3ed17d11a4c5f40f470893f0c968a1ae87b0b3b358321ead5d11"


def sha(a):
    return hashlib.sha256(np.ascontiguousarray(a, dtype="<f8").tobytes()).hexdigest()


def split(on, off, sel, s, tag, lines):
    n_off = np.linalg.norm(off, axis=1)
    n_on = np.linalg.norm(on, axis=1)
    d = np.linalg.norm(on - off, axis=1)
    rel = d / n_off
    cos = (on * off).sum(1) / (n_on * n_off)
    corrected = np.linalg.norm(on - s * off, axis=1) / n_off
    ratio = n_on / n_off
    lines.append(f"\n{tag}")
    lines.append(f"{'row class':18s}{'n':>5s}{'rel displ min':>15s}{'median':>10s}{'max':>10s}{'cos min':>14s}{'cos max':>14s}"
                 f"{'||dE - s E||/||E|| max':>24s}{'median':>12s}")
    res = {}
    for name, m in (("moving", sel), ("NOT moving", ~sel)):
        lines.append(f"{name:18s}{int(m.sum()):5d}{rel[m].min():15.6f}{np.median(rel[m]):10.6f}{rel[m].max():10.6f}"
                     f"{cos[m].min():14.10f}{cos[m].max():14.10f}{corrected[m].max():24.3e}{np.median(corrected[m]):12.3e}")
        res[name] = dict(n=int(m.sum()), rel_min=float(rel[m].min()), rel_median=float(np.median(rel[m])),
                         rel_max=float(rel[m].max()), cos_min=float(cos[m].min()), cos_max=float(cos[m].max()),
                         corrected_max=float(corrected[m].max()), corrected_median=float(np.median(corrected[m])),
                         ratio_min=float(ratio[m].min()), ratio_max=float(ratio[m].max()))
    sep = corrected[sel].min() / max(corrected[~sel].max(), 1e-300)
    lines.append(f"NOT-moving norm ratio: min {ratio[~sel].min()!r} max {ratio[~sel].max()!r} (s = {s!r})")
    lines.append(f"separation on the scale-corrected displacement: {sep:.1e}x (min moving {corrected[sel].min():.6f} vs "
                 f"max not-moving {corrected[~sel].max():.2e})")
    f2 = bool(corrected[~sel].max() < 1e-12 and cos[sel].max() < 1 - 1e-9 and cos[~sel].min() > 1 - 1e-12)
    rel_global = float(np.linalg.norm(on - off) / np.linalg.norm(off))
    lines.append(f"FALSIFIER 2 {'PASS' if f2 else 'FAIL'}: moving rows move in direction (cos <= {cos[sel].max():.6f}), "
                 f"not-moving rows == incumbent x one scalar to {corrected[~sel].max():.1e}; ||dE||/||E|| {rel_global:.6f}; "
                 f"moving rows carry {100*float((d[sel]**2).sum()/(d**2).sum()):.3f}% of ||dE||^2")
    res.update(separation=float(sep), rel_global=rel_global, f2=f2, s=s,
               moving_share=float((d[sel]**2).sum() / (d**2).sum()))
    return res, d, rel


def main() -> None:
    off = np.load(S / "effects-off.npy")
    k5 = np.load(S / "effects-k5.npy")
    on = np.load(S / "effects-composed.npy")
    part = np.load(S / "partial-composed.npz")
    cd4 = np.load(S / "cd4-composed.npz")
    j = {t: json.loads((S / f"stage-{t}.json").read_text()) for t in ("off", "k5", "composed", "partial", "pre")}
    targets = np.asarray(part["targets"]).astype(str)
    partial_rows = np.asarray(part["gene_votes"]).sum(1) > 0
    cd4_rows = np.asarray(cd4["gene_votes"]).sum(1) > 0
    lines = []
    ka = {"effects-off == abeb0749": j["off"]["effects_array_sha256"] == WANT_OFF == sha(off),
          "effects-pre == abeb0749": j["pre"]["effects_array_sha256"] == WANT_OFF,
          "effects-k5 == a28d5be2 (249's arm)": j["k5"]["effects_array_sha256"] == WANT_K5 == sha(k5),
          "partial rows: 15; CD4-carried rows: 276; disjoint": int(partial_rows.sum()) == 15 and int(cd4_rows.sum()) == 276
                                                                 and not np.any(partial_rows & cd4_rows)}
    for k, v in ka.items():
        lines.append(f"[{'PASS' if v else 'FAIL'}] KA {k}")
    lines.append(f"effects: off {j['off']['effects_array_sha256'][:16]}  k5 {j['k5']['effects_array_sha256'][:16]}  "
                 f"composed {j['composed']['effects_array_sha256'][:16]}  partial {j['partial']['effects_array_sha256'][:16]}  shape {on.shape}")
    pres = {t: j[t]["audit"]["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"] for t in ("off", "k5", "composed", "partial")}
    lines.append(f"merged preserve scalars: {pres}")
    s_ck = pres["composed"] / pres["k5"]
    s_po = pres["partial"] / pres["off"]
    res_ck, d_ck, rel_ck = split(on, k5, partial_rows, s_ck, "COMPOSED vs cd4svd5 (k5): the partial rows vs the rest", lines)
    res_po, _, _ = split(np.load(S / "effects-partial.npy"), off, partial_rows, s_po, "PARTIAL vs h1union (off): the partial rows vs the rest", lines)
    s_co = pres["composed"] / pres["off"]
    res_co, _, _ = split(on, off, partial_rows | cd4_rows, s_co, "COMPOSED vs h1union (off): the 291 CD4 rows (276 deflated + 15 partial) vs the 9 rows CD4 never touches", lines)
    zero_rows = int(np.sum(d_ck == 0.0))
    f3 = zero_rows == 0 and res_ck["rel_global"] > 0.001
    lines.append(f"\nNON-VACUITY (composed vs cd4svd5): ||dE||/||E|| = {res_ck['rel_global']:.6f}; rows with displacement exactly 0: "
                 f"{zero_rows}/300; per-row min {rel_ck.min():.6f}")
    lines.append(f"FALSIFIER 3 (effect level) {'PASS' if f3 else 'FAIL'}")
    order = np.argsort(-rel_ck)[:15]
    lines.append("the 15 most displaced rows (composed vs cd4svd5): " + ", ".join(f"{targets[i]} {rel_ck[i]:.3f}" for i in order))
    lines.append("partial targets: " + ", ".join(targets[partial_rows]))
    # the deflation audit and CD4 fraction must be identical between k5 and composed
    ck = j["composed"]["audit"]["cd4"]; kk = j["k5"]["audit"]["cd4"]
    same_cd4 = (ck["fraction"] == kk["fraction"] and json.dumps(ck["deflation"], sort_keys=True) == json.dumps(kk["deflation"], sort_keys=True)
                and json.dumps(ck["residual_frame"], sort_keys=True) == json.dumps(kk["residual_frame"], sort_keys=True))
    lines.append(f"[{'PASS' if same_cd4 else 'FAIL'}] CD4's fraction, frame and deflation audit identical between k5 and composed "
                 f"(fraction {ck['fraction']!r}, energy removed {ck['deflation']['energy_removed_frac']!r})")
    same_union = (j["composed"]["audit"]["h1_union"]["residual_frame"]["solved_fraction"]
                  == j["k5"]["audit"]["h1_union"]["residual_frame"]["solved_fraction"]
                  == j["off"]["audit"]["h1_union"]["residual_frame"]["solved_fraction"])
    lines.append(f"[{'PASS' if same_union else 'FAIL'}] the union's solved fraction identical across off / k5 / composed "
                 f"({j['off']['audit']['h1_union']['residual_frame']['solved_fraction']!r})")
    pb = j["composed"]["audit"]["cd4_partial"]
    same_pfrac = pb["residual_frame"]["solved_fraction"] == j["partial"]["audit"]["cd4_partial"]["residual_frame"]["solved_fraction"]
    lines.append(f"[{'PASS' if same_pfrac else 'FAIL'}] the partial's solved fraction identical with and without the deflation "
                 f"({pb['residual_frame']['solved_fraction']!r}; shipped all-rows rule would solve {pb['residual_frame']['shipped_median_rule_would_solve']!r})")
    ok = all(ka.values()) and res_ck["f2"] and res_po["f2"] and res_co["f2"] and f3 and same_cd4 and same_union and same_pfrac
    out = dict(ka=ka, composed_vs_k5=res_ck, partial_vs_off=res_po, composed_vs_off=res_co, zero_rows=zero_rows, f3=f3,
               preserve=pres, same_cd4=same_cd4, same_union=same_union, same_partial_fraction=same_pfrac,
               partial_targets=targets[partial_rows].tolist(), ok=ok)
    (WORK / "effdiff254.json").write_text(json.dumps(out, indent=1))
    (WORK / "effdiff254.log").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nEFFDIFF {'PASS' if ok else 'FAIL'}")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
