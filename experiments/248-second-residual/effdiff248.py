"""248 -- the rider's geometry at the EFFECT-MATRIX level (222 section 6a's exact form), on the
staged 300 x 18,533 matrices: `effects-cum1.npy` (and cum4) vs `effects-off.npy` (the shipped
h1union path through the same builder).  Rows union 2 carries must move in DIRECTION; rows it does
not must equal the incumbent's times ONE shared scalar -- the merged preserve ratio -- to machine
precision.  Non-vacuity: ||dE||/||E|| and the number of rows with exactly zero displacement.
Also the two unions' overlap and the frames the builder wrote.  No package is built, so this is
the level at which the arm's profile split is stated.

    uv run python effdiff248.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
S = HERE / "work" / "stage"
lines = []


def P(t=""):
    print(t)
    lines.append(t)


def main() -> None:
    off = np.load(S / "effects-off.npy")
    j_off = json.loads((S / "stage-off.json").read_text())
    s_off = j_off["audit"]["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"]
    out = {}
    for tag in ("cum1", "cum4"):
        on = np.load(S / f"effects-{tag}.npy")
        j_on = json.loads((S / f"stage-{tag}.json").read_text())
        u2 = np.load(S / f"union2-{tag}.npz")
        u1 = np.load(S / f"union-{tag}.npz")
        carry2 = np.asarray(u2["target_condition_count"]) > 0
        carry1 = np.asarray(u1["target_condition_count"]) > 0
        targets = np.asarray(u2["targets"]).astype(str)
        s_on = j_on["audit"]["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"]
        s = s_on / s_off
        fr2 = j_on["audit"]["union2"]["residual_frame"]
        fr1 = j_on["audit"]["h1_union"]["residual_frame"]
        P(f"=== {tag}: {j_on['effects_array_sha256'][:16]} vs off {j_off['effects_array_sha256'][:16]}  "
          f"(union 2 = {j_on['union2_sources']} @ {j_on['union2_g']}) ===")
        P(f"union 2 carries {int(carry2.sum())}/300 (coverage {fr2['coverage_frac']:.3f}); union 1 carries {int(carry1.sum())}/300; "
          f"both {int((carry1 & carry2).sum())}; union 2 only {int((carry2 & ~carry1).sum())}")
        P(f"frames: cd4 fraction {j_on['audit']['cd4']['fraction']!r} (off {j_off['audit']['cd4']['fraction']!r}); "
          f"union 1 solved {fr1['solved_fraction']!r} (off {j_off['audit']['h1_union']['residual_frame']['solved_fraction']!r}); "
          f"union 2 solved {fr2['solved_fraction']!r}, g_per_unit {fr2['g_per_unit_fraction']:.6f}, "
          f"shipped-median rule would solve {fr2['shipped_median_rule_would_solve']!r}, median carrying {fr2['residual_median_row_norm_over_carrying_rows']:.4f}, "
          f"median all rows {fr2['residual_median_row_norm_over_all_rows']:.4f}")
        P(f"merged preserve scalar: off {s_off!r} -> {s_on!r}  ratio s = {s!r}  (written identically under cd4 and h1_union and union2: "
          f"{j_on['audit']['h1_union']['merged_preserve']['post_add_global_preserve_scale'] == s_on == j_on['audit']['union2']['merged_preserve']['post_add_global_preserve_scale']})")
        d = on - off
        rel = np.linalg.norm(d, axis=1) / np.linalg.norm(off, axis=1)
        cos = np.einsum("ij,ij->i", on, off) / (np.linalg.norm(on, axis=1) * np.linalg.norm(off, axis=1))
        corr = np.linalg.norm(on - s * off, axis=1) / np.linalg.norm(off, axis=1)
        for label, sel in (("union-2-carried", carry2), ("NOT carried by union 2", ~carry2)):
            P(f"  {label:24s} n {int(sel.sum()):3d}  rel displacement min/median/max {rel[sel].min():.6f}/{np.median(rel[sel]):.6f}/{rel[sel].max():.6f}  "
              f"cosine(on, off) min/max {cos[sel].min():.10f}/{cos[sel].max():.10f}  ||on - s off||/||off|| max {corr[sel].max():.2e}")
        sep = corr[carry2].min() / max(corr[~carry2].max(), 1e-300)
        share = float((d[carry2] ** 2).sum() / (d ** 2).sum())
        P(f"  scale-corrected separation carried-min / non-carried-max = {sep:.1e}x; carried rows hold {100*share:.3f}% of ||dE||^2; "
          f"||dE||/||E|| {np.linalg.norm(d)/np.linalg.norm(off):.6f}; rows with exactly zero displacement {int((rel == 0).sum())}/300")
        top = np.argsort(-rel)[:5]
        P(f"  most displaced targets: " + ", ".join(f"{targets[i]} {rel[i]:.3f}" for i in top))
        out[tag] = dict(sha=j_on["effects_array_sha256"], off_sha=j_off["effects_array_sha256"], carry2=int(carry2.sum()),
                        carry1=int(carry1.sum()), both=int((carry1 & carry2).sum()), s=s, separation=float(sep),
                        share_carried=share, rel_total=float(np.linalg.norm(d) / np.linalg.norm(off)),
                        zero_rows=int((rel == 0).sum()), frames=dict(cd4=j_on["audit"]["cd4"]["fraction"], u1=fr1, u2=fr2),
                        carried=dict(rel_min=float(rel[carry2].min()), rel_median=float(np.median(rel[carry2])), rel_max=float(rel[carry2].max()),
                                     cos_min=float(cos[carry2].min()), cos_max=float(cos[carry2].max())),
                        not_carried=dict(rel_min=float(rel[~carry2].min()), rel_max=float(rel[~carry2].max()),
                                         cos_min=float(cos[~carry2].min()), corr_max=float(corr[~carry2].max())))
    (HERE / "work/effdiff248.json").write_text(json.dumps(out, indent=1))
    (HERE / "work/effdiff248.log").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
