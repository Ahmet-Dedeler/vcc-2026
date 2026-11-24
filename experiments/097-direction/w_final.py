"""The headline table: every candidate arm on both harnesses, with the amplitude control.

`w_scale.py` showed `W` is not exactly amplitude-invariant through the floor+gate nonlinearity,
so each arm's raw `sqrt(W/W_base)` is reported alongside a SCALE-CORRECTED version: interpolate
the pure-rescale curve `W(s)` at the arm's own effective scale `s_eff = med||dh||_arm /
med||dh||_base`, and divide it out.  What survives is the part of the arm an amplitude move
could not have produced.
"""
from __future__ import annotations

import json

import numpy as np
from w_arms import key_normrel, soft
from wkit import ROOT, Ctx, fixed, uniform

SCALES = (0.25, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.25, 1.5, 2.0)


def main():
    out = {}
    for name in ("h1-test", "orion"):
        c = Ctx(name)
        b = c.run(uniform)
        Wb, nb = b["W"], np.linalg.norm(b["E"], axis=1)
        dhb = float(np.median(np.linalg.norm(c.dh_exact(b["E"]), axis=1)))

        # pure-rescale reference curve, for the amplitude control
        ss, ww, dd = [], [], []
        for s in SCALES:
            dh = c.dh_exact(b["E"] * s)
            ss.append(s)
            ww.append(c.XAW(dh)[2] / Wb)
            dd.append(float(np.median(np.linalg.norm(dh, axis=1))) / dhb)
        dd, ww = np.asarray(dd), np.asarray(ww)

        def corr(seff):
            return float(np.interp(seff, dd, ww))

        keys = key_normrel(c)
        arms = {
            "a1 uniform (base)": ("full", uniform),
            "a3 k562 share 0.35 (orion argmax)":
                ("full", fixed({"k562_gwps": 0.35, "orion_hek293t": 0.65})),
            "a3 k562 share 0.60 (h1-test argmax)":
                ("full", fixed({"k562_gwps": 0.60, "orion_hek293t": 0.40})),
            "a4 basal exp(10*sim)":
                ("full", fixed({n: float(np.exp(10.0 * v))
                                for n, v in c.model.basal_weights(c.b_ctrl).items()})),
            "a6 normrel^-1.25": ("full", soft(c, keys, -1.25)),
            "a6 normrel^-1.50": ("full", soft(c, keys, -1.50)),
            "a6 normrel^-1.75": ("full", soft(c, keys, -1.75)),
            "a6 normrel^-1.25 dir-only": ("dir", soft(c, keys, -1.25)),
            "a6 normrel^-1.50 dir-only": ("dir", soft(c, keys, -1.50)),
            "a6 normrel^-1.75 dir-only": ("dir", soft(c, keys, -1.75)),
        }
        rows = []
        print(f"\n########## {name}   W_base={Wb:.8g}  med||dh||_base={dhb:.5f}")
        print(f"  {'arm':<38}{'X':>12}{'A':>12}{'W':>12}{'mult':>8}"
              f"{'s_eff':>8}{'scale-only':>12}{'corrected':>11}")
        for lbl, (mode, wf) in arms.items():
            E = c.effects(wf)
            if mode == "dir":
                E = E * (nb / np.maximum(np.linalg.norm(E, axis=1), 1e-30))[:, None]
            dh = c.dh_exact(E)
            X, A, W = c.XAW(dh)
            seff = float(np.median(np.linalg.norm(dh, axis=1))) / dhb
            m, ms = np.sqrt(W / Wb), np.sqrt(corr(seff))
            rows.append(dict(arm=lbl, X=X, A=A, W=W, mult=float(m), s_eff=seff,
                             scale_only=float(ms), corrected=float(m / ms)))
            print(f"  {lbl:<38}{X:>12.6g}{A:>12.6g}{W:>12.6g}{m:>8.4f}"
                  f"{seff:>8.4f}{ms:>12.4f}{m/ms:>11.4f}")
        out[name] = dict(W_base=Wb, dh_base=dhb, arms=rows,
                         scale_curve={"dh_ratio": dd.tolist(), "W_ratio": ww.tolist()})
        del c

    # the three-donor arms, each measured against its own harness's two-donor base
    for tag, harness, extra in (("h1-test+orion_hct116", "h1-test", "orion_hct116"),
                                ("orion+h1_2025", "orion", "h1_2025")):
        c = Ctx(harness, extra_donors=(extra,))
        b = c.run(uniform)
        nb = np.linalg.norm(b["E"], axis=1)
        keys = key_normrel(c)
        Wb, dhb = out[harness]["W_base"], out[harness]["dh_base"]
        dd = np.asarray(out[harness]["scale_curve"]["dh_ratio"])
        ww = np.asarray(out[harness]["scale_curve"]["W_ratio"])
        print(f"\n########## {tag}   (vs that harness's 2-donor base W={Wb:.8g})")
        rows = []
        for lbl, (mode, wf) in {"a7 3-donor uniform": ("full", uniform),
                                "a7 3-donor + normrel^-1.50": ("full", soft(c, keys, -1.5)),
                                "a7 3-donor + normrel^-1.50 dir-only":
                                    ("dir", soft(c, keys, -1.5))}.items():
            E = c.effects(wf)
            if mode == "dir":
                E = E * (nb / np.maximum(np.linalg.norm(E, axis=1), 1e-30))[:, None]
            dh = c.dh_exact(E)
            X, A, W = c.XAW(dh)
            seff = float(np.median(np.linalg.norm(dh, axis=1))) / dhb
            m, ms = np.sqrt(W / Wb), np.sqrt(float(np.interp(seff, dd, ww)))
            rows.append(dict(arm=lbl, X=X, A=A, W=W, mult=float(m), s_eff=seff,
                             scale_only=float(ms), corrected=float(m / ms)))
            print(f"  {lbl:<38}{X:>12.6g}{A:>12.6g}{W:>12.6g}{m:>8.4f}"
                  f"{seff:>8.4f}{ms:>12.4f}{m/ms:>11.4f}")
        out[tag] = dict(arms=rows)
        del c

    (ROOT / "experiments/097-direction/work/final.json").write_text(json.dumps(out, indent=2))
    print("\nwrote work/final.json")


if __name__ == "__main__":
    main()
