#!/usr/bin/env python3
"""553 — EMISSION-side diagnostics of a surface's arms (item 5): what was emitted, never how it scores. Reads only the
run553 outputs (run553.json, diag / stats / metrics files of the predicted blocks) — the predicted blocks' own statistics
against the CONTROL cells — so it is allowed on a sealed surface (RESERVE.md section 6).

Per arm x seed (zs form unless --pre): the ballast K histogram; the prediction's own call set per target (q < 0.05, the
target's own gene excluded) — size, up/down split, and its overlap with the incumbent's (`cpm`, same seed) call set; the
emitter's call-set record (prescribed / delivered / missed / unsuppressible); the MDE hook's coverage (signed complement
genes) and the prior's flips; library per cell; and two identity checks against `cpm` at the same seed that a DE-only
knob must pass: the per-gene block sums (pds / mse inputs) and pass 1 (ballast K per target) are equal.
usage: uv run python experiments/553-r7b-harness/emitdiag553.py <surface> <mode,mode,...> [--seeds 0,1,2]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import lib553 as S  # noqa: E402

TAG = "base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse"


def files(surface, mode, seed, zs, tag=""):
    d = HERE / "work" / f"{surface}_{mode}{tag}"
    t = TAG.format(seed=seed) + ("_dnk250zs" if zs else "")
    return d / f"stats-{t}.npz", d / f"diag-{t}.json", d / f"diag-{TAG.format(seed=seed)}.json", d / "run553.json"


def _q(dg, key, q):
    v = np.array([np.nan if d.get(key) is None else float(d[key]) for d in dg], dtype=float)
    return float(np.nanquantile(v, q)) if np.isfinite(v).any() else None


def calls(st, genes_keep, targets):
    own = np.equal.outer(np.asarray(targets).astype(str), genes_keep)
    sig = (np.asarray(st["q"]) < 0.05) & ~own
    return sig, np.asarray(st["lfc"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("surface")
    ap.add_argument("modes")
    ap.add_argument("--seeds", default="0,1,2")
    ap.add_argument("--pre", action="store_true")
    ap.add_argument("--tag", default="", help="run directory suffix, e.g. -flex")
    a = ap.parse_args()
    spec = S.surface_spec(a.surface)
    cr = S.load_controlref(spec, need_cs=False)
    genes_keep = cr["genes"].astype(str)[np.asarray(cr["keep"], bool)]
    targets = cr["targets"].astype(str)
    seeds = [int(s) for s in a.seeds.split(",")]
    zs = not a.pre
    modes = ["cpm"] + [m for m in a.modes.split(",") if m and m != "cpm"]
    base = {}
    out = {}
    for mode in modes:
        rows = []
        for s in seeds:
            fs, fd, fd1, fr = files(a.surface, mode, s, zs, a.tag)
            if not fs.exists():
                continue
            st = np.load(fs, allow_pickle=True)
            dg1 = json.loads(fd1.read_text())               # pass-1 (pre-zs) block records: K, call set, MDE, prior
            # the run's own target rows (a run with --exclude-targets has fewer rows than the surface)
            sig, lfc = calls(st, genes_keep, st["targets"].astype(str) if "targets" in st.files else targets)
            n = sig.sum(1)
            up = (sig & (lfc > 0)).sum(1)
            r = dict(seed=s, n_called_median=float(np.median(n)), n_called_q10_q90=[float(np.quantile(n, .1)), float(np.quantile(n, .9))],
                     n_called_total=int(n.sum()), up_frac=float(up.sum() / max(n.sum(), 1)),
                     K_hist={int(k): int(v) for k, v in zip(*np.unique([d["ballast_K"] for d in dg1], return_counts=True))},
                     prescribed_median=float(np.median([d["call_set_prescribed"] for d in dg1])),
                     delivered_median=float(np.median([d["call_set_delivered"] for d in dg1])),
                     missed_mean=float(np.mean([d["call_set_missed"] for d in dg1])),
                     unsuppressible_mean=float(np.mean([d["call_set_unsuppressible"] for d in dg1])),
                     mde_sign_cover_median=float(np.median([d.get("mde_n_sign_cover", 0) for d in dg1])),
                     mde_sign_flip_median=float(np.median([d.get("mde_n_sign_flip", 0) for d in dg1])),
                     prior_flips_median=float(np.median([d["prior_n_flip_target"] for d in dg1])),
                     per_cell_total_median=float(np.median([d["per_cell_total"] for d in dg1])),
                     n_de_universe=int(dg1[0]["n_genes_de_universe"]),
                     # the steering doses (the CLEAN4 block prior's lambda; the MDE add's, which carries the hook)
                     prior_lambda_min_median=[_q(dg1, "prior_lambda", 0), _q(dg1, "prior_lambda", .5)],
                     mde_lambda_min_median=[_q(dg1, "mde_lambda", 0), _q(dg1, "mde_lambda", .5)],
                     n_mde_lambda_below_half=int(sum(1 for d in dg1 if d.get("mde_lambda") is not None and d["mde_lambda"] < .5)),
                     sub_gate_counts_per_cell_median=_q(dg1, "sub_gate_counts_per_cell", .5),
                     unresolved_total=int(sum(int(d.get("unresolved", 0) or 0) for d in dg1)))
            if mode == "cpm":
                base[s] = dict(sig=sig, sums=np.asarray(st["sums"]), K=[d["ballast_K"] for d in dg1])
            elif s in base:
                b = base[s]
                inter = (sig & b["sig"]).sum(1)
                r.update(calls_new_vs_cpm_median=float(np.median(n - inter)),
                         calls_lost_vs_cpm_median=float(np.median(b["sig"].sum(1) - inter)),
                         sums_equal_cpm=bool(np.array_equal(np.asarray(st["sums"]), b["sums"])),
                         pass1_K_equal_cpm=bool([d["ballast_K"] for d in dg1] == b["K"]))
            rows.append(r)
        fr0 = files(a.surface, mode, seeds[0], zs, a.tag)[3]
        rec = json.loads(fr0.read_text()) if fr0.exists() else {}
        out[mode] = dict(seeds=rows, hook_signs=rec.get("hook_signs"), hook_signs_down=rec.get("hook_signs_down"),
                         key_sha256=rec.get("key_sha256"), classes=rec.get("classes"), bulk_key=rec.get("bulk_key"))
    S.write_json(HERE / "work" / f"emitdiag553_{a.surface}{a.tag}.json", out)
    print(f"| arm | seeds | called / target (median) | up frac | K hist (s{seeds[0]}) | prescribed / delivered | missed | "
          f"unsuppressible | hook signs | calls new vs cpm | sums == cpm | pass-1 K == cpm | prior λ min / med | MDE λ min / med |")
    print("|---|---|---:|---:|---|---|---:|---:|---:|---:|---|---|---|---|")
    def lam(rs, k):
        v = [r[k] for r in rs if r[k][0] is not None]
        return "n/a" if not v else f"{min(x[0] for x in v):.3f} / {np.mean([x[1] for x in v]):.3f}"

    for mode, o in out.items():
        rs = o["seeds"]
        if not rs:
            continue
        f = lambda k: np.mean([r[k] for r in rs])  # noqa: E731
        is_base = mode == "cpm"
        new_calls = "" if is_base else f"{f('calls_new_vs_cpm_median'):.0f}"
        sums_eq = "" if is_base else str(all(r.get("sums_equal_cpm", False) for r in rs))
        k_eq = "" if is_base else str(all(r.get("pass1_K_equal_cpm", False) for r in rs))
        hook = (o["hook_signs"] or {}).get("weak", "")
        print(f"| {mode} | {len(rs)} | {f('n_called_median'):.0f} | {f('up_frac'):.3f} | {rs[0]['K_hist']} | "
              f"{f('prescribed_median'):.0f} / {f('delivered_median'):.0f} | {f('missed_mean'):.1f} | "
              f"{f('unsuppressible_mean'):.1f} | {hook} | {new_calls} | {sums_eq} | {k_eq} | {lam(rs, 'prior_lambda_min_median')} | "
              f"{lam(rs, 'mde_lambda_min_median')} |")

if __name__ == "__main__":
    main()
