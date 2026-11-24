"""374: Jurkat blend knee on 341's saved source-frame predictions.

Reads 341's saved predictions.npz per line (no refit), blends
crossline_diag -> same_line_full under 346's fixed-norm rule, scores each
beta with 341/334's own endpoint functions, intervals via core341.interval.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/341-known-line-calibration"))
sys.path.insert(0, str(ROOT / "experiments/334-text-conditioned-transfer"))
import core341 as C341
import core334 as C334

LINES = ["k562_gwps", "orion_hct116", "orion_hek293t", "rpe1",
         "nadig_jurkat", "nadig_hepg2"]
BETAS = [0, 0.1, 0.25, 0.35, 0.5, 0.75, 1.0]
WORK341 = ROOT / "experiments/341-known-line-calibration/work"


def blend(d, p, beta):
    nd = np.linalg.norm(d, axis=1, keepdims=True)
    np_ = np.linalg.norm(p, axis=1, keepdims=True)
    assert (nd > 0).all() and (np_ > 0).all()
    m = (1 - beta) * d / nd + beta * p / np_
    nm = np.linalg.norm(m, axis=1, keepdims=True)
    return m * (nd / nm)


def endpoints(pred, truth):
    e = C334.endpoints(pred, truth)
    s = C341.sign_endpoints(pred, truth)
    return {"weighted_sign": s["weighted_sign"],
            "top64_weighted_sign": s["top64_weighted_sign"],
            "retrieval": e["retrieval"], "cosine": e["cosine"]}


def main():
    out = HERE / "work"
    out.mkdir(parents=True, exist_ok=True)
    # ---- load ----
    data = {}
    for line in LINES:
        z = np.load(WORK341 / line / "predictions.npz")
        data[line] = {k: np.asarray(z[k], dtype=np.float64)
                      for k in ("truth", "crossline_diag", "same_line_full",
                                "crossline_full", "labels_shuffled_full")}
        assert list(z["targets"]) == json.loads(
            (ROOT / "experiments/341-known-line-calibration/selection.json")
            .read_text())["targets"]
    results = json.loads((WORK341 / "results.json").read_text())

    # ---- known-answer gate ----
    print("== known-answer gate (recompute from saved predictions) ==")
    gate = {}
    for arm, key in (("same_line_full", "same_line_full"),
                     ("crossline_diag", "crossline_diag")):
        m = {ln: endpoints(data[ln][key], data[ln]["truth"]) for ln in LINES}
        for metric in ("weighted_sign", "retrieval", "cosine",
                       "top64_weighted_sign"):
            macro = float(np.mean([v[metric].mean() for v in m.values()]))
            rec = results["means"][arm][metric]
            diff = macro - rec
            gate[f"{arm}/{metric}"] = {"mine": macro, "record": rec,
                                       "diff": diff}
            print(f"{arm}/{metric}: mine={macro:.15f} rec={rec:.15f} "
                  f"diff={diff:.3e}")
    # Jurkat contrasts
    j = "nadig_jurkat"
    ej_s = endpoints(data[j]["same_line_full"], data[j]["truth"])
    ej_d = endpoints(data[j]["crossline_diag"], data[j]["truth"])
    ej_f = endpoints(data[j]["crossline_full"], data[j]["truth"])
    for name, delta in (
            ("jurkat wsign same-crossdiag",
             float(ej_s["weighted_sign"].mean() - ej_d["weighted_sign"].mean())),
            ("jurkat wsign same-crossfull",
             float(ej_s["weighted_sign"].mean() - ej_f["weighted_sign"].mean())),
            ("pooled retr same-crossdiag",
             float(np.mean([endpoints(data[ln]["same_line_full"],
                                      data[ln]["truth"])["retrieval"].mean()
                             - endpoints(data[ln]["crossline_diag"],
                                         data[ln]["truth"])["retrieval"].mean()
                             for ln in LINES])))):
        print(f"{name}: {delta:.8f}")
    print("record jurkat vs diag : 0.03380218 ; vs full: 0.02800431 ; "
          "pooled retr vs diag: -0.04829784")
    gate_ok = all(abs(v["diff"]) < 1e-9 for v in gate.values())
    print("GATE:", "PASS" if gate_ok else "FAIL")
    if not gate_ok:
        (out / "curve.json").write_text(json.dumps(
            {"gate": gate, "gate_ok": False}, indent=2))
        print("stopping: known-answer gate failed")
        sys.exit(1)

    # beta=0 identity + alternate-form equivalence
    d0 = data[j]["crossline_diag"]
    p0 = data[j]["same_line_full"]
    b0 = blend(d0, p0, 0.0)
    print("beta0 max abs diff:", float(np.abs(b0 - d0).max()))
    assert float(np.abs(b0 - d0).max()) < 1e-12
    # alternate form: scale ridge first, blend, rescale
    beta = 0.35
    nd = np.linalg.norm(d0, axis=1, keepdims=True)
    np_ = np.linalg.norm(p0, axis=1, keepdims=True)
    alt = (1 - beta) * d0 + beta * p0 * (nd / np_)
    alt = alt * (nd / np.linalg.norm(alt, axis=1, keepdims=True))
    print("alt-form max abs diff:",
          float(np.abs(alt - blend(d0, p0, beta)).max()))
    assert float(np.abs(alt - blend(d0, p0, beta)).max()) < 1e-14
    # beta=1 norm check
    b1 = blend(d0, p0, 1.0)
    print("beta1 row-norm max err:",
          float(np.abs(np.linalg.norm(b1, axis=1)
                       - np.linalg.norm(d0, axis=1)).max()))

    # ---- curves ----
    curve = {"betas": BETAS, "jurkat": {}, "pooled5": {}, "per_line": {}}
    jur_base = endpoints(d0, data[j]["truth"])
    jur_metrics = {}
    for beta in BETAS:
        pred = blend(d0, p0, beta)
        jur_metrics[beta] = endpoints(pred, data[j]["truth"])
    # pooled-five baseline + metrics (other five lines)
    others = [ln for ln in LINES if ln != j]
    base5 = {m: np.stack([endpoints(data[ln]["crossline_diag"],
                                   data[ln]["truth"])[m] for ln in others])
             for m in ("weighted_sign", "top64_weighted_sign",
                       "retrieval", "cosine")}
    met5 = {}
    for beta in BETAS:
        met5[beta] = {m: np.stack(
            [endpoints(blend(data[ln]["crossline_diag"],
                             data[ln]["same_line_full"], beta),
                       data[ln]["truth"])[m] for ln in others])
            for m in base5}
    for beta in BETAS:
        row = {}
        for m in ("weighted_sign", "top64_weighted_sign", "retrieval",
                  "cosine"):
            v = jur_metrics[beta][m]
            row[m] = {"mean": float(v.mean())}
            if beta != 0:
                delta = (v - jur_base[m])[None, :]
                iv = C341.interval(delta)
                row[m].update({"delta": float(iv["mean"]),
                               "target_ci95": iv["target_ci95"]})
        curve["jurkat"][str(beta)] = row
    for beta in BETAS:
        row = {}
        for m in base5:
            v = met5[beta][m]
            row[m] = {"mean": float(v.mean())}
            if beta != 0:
                iv = C341.interval(v - base5[m])
                row[m].update({"delta": float(iv["mean"]),
                               "target_ci95": iv["target_ci95"],
                               "crossed_ci95": iv["crossed_ci95"]})
        curve["pooled5"][str(beta)] = row
    # per-line means for appendix
    for ln in LINES:
        dd = data[ln]["crossline_diag"]
        pp = data[ln]["same_line_full"]
        tt = data[ln]["truth"]
        base = endpoints(dd, tt)
        lr = {}
        for beta in BETAS:
            e = endpoints(blend(dd, pp, beta), tt)
            lr[str(beta)] = {m: {"mean": float(e[m].mean()),
                                 "delta": float(e[m].mean()
                                                - base[m].mean())}
                             for m in base}
        curve["per_line"][ln] = lr
    # shuffled floor on Jurkat
    sh = endpoints(data[j]["labels_shuffled_full"], data[j]["truth"])
    curve["jurkat_shuffled_floor"] = {m: float(sh[m].mean())
                                      for m in sh}
    curve["gate"] = gate
    curve["gate_ok"] = True
    (out / "curve.json").write_text(json.dumps(curve, indent=2))
    # CSV for the report
    import csv
    with (out / "curve.csv").open("w", newline="") as h:
        w = csv.writer(h)
        w.writerow(["scope", "beta", "wsign_mean", "wsign_delta",
                    "wsign_lo", "wsign_hi", "top64_mean", "top64_delta",
                    "retr_mean", "retr_delta", "retr_lo", "retr_hi",
                    "cos_mean", "cos_delta", "pds_raw_price"])
        for scope, blk, xkey in (("jurkat", curve["jurkat"], None),
                                 ("pooled5", curve["pooled5"], None)):
            for beta in BETAS:
                r = blk[str(beta)]
                if beta == 0:
                    w.writerow([scope, beta, f"{r['weighted_sign']['mean']:.6f}",
                                "0", "", "", f"{r['top64_weighted_sign']['mean']:.6f}",
                                "0", f"{r['retrieval']['mean']:.6f}", "0", "",
                                "", f"{r['cosine']['mean']:.6f}", "0", "0"])
                else:
                    w.writerow([scope, beta,
                                f"{r['weighted_sign']['mean']:.6f}",
                                f"{r['weighted_sign']['delta']:+.6f}",
                                f"{r['weighted_sign']['target_ci95'][0]:+.6f}",
                                f"{r['weighted_sign']['target_ci95'][1]:+.6f}",
                                f"{r['top64_weighted_sign']['mean']:.6f}",
                                f"{r['top64_weighted_sign']['delta']:+.6f}",
                                f"{r['retrieval']['mean']:.6f}",
                                f"{r['retrieval']['delta']:+.6f}",
                                f"{r['retrieval']['target_ci95'][0]:+.6f}",
                                f"{r['retrieval']['target_ci95'][1]:+.6f}",
                                f"{r['cosine']['mean']:.6f}",
                                f"{r['cosine']['delta']:+.6f}",
                                f"{r['retrieval']['delta'] * 1.009:+.6f}"])
    print("wrote work/curve.json + work/curve.csv")
    print("jurkat shuffled floor:", curve["jurkat_shuffled_floor"])


if __name__ == "__main__":
    main()
