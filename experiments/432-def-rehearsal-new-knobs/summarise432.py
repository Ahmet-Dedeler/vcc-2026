"""432: collect every stage's wall clock, peak RSS and outcome from the rehearsal's logs into one table.

Sources: the driver's `timings.json` per work dir (its own stages), `/usr/bin/time -l` tails of the wrapper logs
(`logs/repin-panel-support.log`, `logs/build.log`), the ntclevel/audit outputs and the rc audit dumps.
usage: uv run python summarise432.py            (writes work/summary432.json and prints markdown)
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

W = Path(__file__).resolve().parent / "work"


def time_tail(log: Path) -> dict:
    if not log.exists():
        return {}
    t = log.read_text(errors="replace")
    real = re.findall(r"^\s*([\d.]+) real", t, re.M)
    rss = re.findall(r"(\d+)\s+maximum resident set size", t)
    return dict(seconds=float(real[-1]) if real else None, peak_rss_gib=round(int(rss[-1]) / 2**30, 2) if rss else None)


def rc_dump(path: Path) -> dict:
    if not path.exists():
        return {}
    rows = [json.loads(l) for l in path.read_text().splitlines() if l.strip()]
    out = {"blocks": len(rows)}
    for ctx in sorted({r["context"] for r in rows}):
        R = [r for r in rows if r["context"] == ctx]
        d = {}
        for k in ("mde_lambda", "prior_lambda", "zero_repaired", "zs_n_moved", "ballast_K", "clip_frac", "mde_overfill_frac"):
            v = np.array([r[k] for r in R if k in r], dtype=float)
            if v.size:
                d[k] = dict(n=int(v.size), median=float(np.median(v)), q10=float(np.quantile(v, .1)),
                            q90=float(np.quantile(v, .9)), min=float(v.min()), max=float(v.max()),
                            frac_zero=float(np.mean(v == 0)), frac_lt1=float(np.mean(v < 1)), total=float(v.sum()))
        out[ctx] = d
    return out


def main() -> None:
    S = {}
    for wd in sorted(p for p in W.iterdir() if p.is_dir()):
        row = {}
        tj = wd / "timings.json"
        if tj.exists():
            row["driver"] = {k: {kk: v.get(kk) for kk in ("seconds", "peak_rss_gib", "returncode")}
                             for k, v in json.loads(tj.read_text()).items() if not k.startswith("_")}
        for st in ("repin-panel-support", "build", "prep"):
            tt = time_tail(wd / "logs" / f"{st}.log")
            if tt:
                row[st] = tt
        nt = wd / f"{wd.name}.ntclevel.json"
        if nt.exists():
            row["ntclevel"] = {c: {k: v[k] for k in ("Cp46", "A", "B", "claim", "rho", "cliff", "scale")}
                               for c, v in json.loads(nt.read_text())["per_context"].items()}
            row["ntclevel_time"] = time_tail(wd / f"{wd.name}.ntclevel.out")
        rd = wd / f"rc_audit_{wd.name}.jsonl"
        if rd.exists():
            row["rc_dump"] = rc_dump(rd)
        au = wd / "audit428.out"
        if au.exists():
            row["audit428_head"] = au.read_text().splitlines()[:4]
        meta = wd / f"{wd.name}.h5ad.meta.json"
        if meta.exists():
            m = json.loads(meta.read_text())
            row["meta"] = dict(sha256=m["output"].get("sha256"), peak_rss_gib=m.get("peak_rss_gib"),
                               bytes=m["output"].get("bytes"), n_rows=m["output"].get("n_rows"))
        if row:
            S[wd.name] = row
    # (c) did the re-level move the energies as intended? B_c(level) / B_c(unlevelled) against s_c^2, and the levelled
    # build's own cliff (should read exp(-sd * Phi^-1(q)) = 1.0105 at sd 0.01, i.e. EV-scale 1.0000)
    for pre in ("", "ef_"):
        n1, n2 = W / f"{pre}b1/{pre}b1.ntclevel.json", W / f"{pre}b2/{pre}b2.ntclevel.json"
        if n1.exists() and n2.exists():
            p1, p2 = json.loads(n1.read_text())["per_context"], json.loads(n2.read_text())["per_context"]
            S[f"{pre}relevel_check"] = {c: dict(scale_asked=p1[c]["scale"], B_ratio=p2[c]["B"] / p1[c]["B"],
                                                 scale_sq=p1[c]["scale"] ** 2,
                                                 B_ratio_over_scale_sq=p2[c]["B"] / p1[c]["B"] / p1[c]["scale"] ** 2,
                                                 A_ratio=p2[c]["A"] / p1[c]["A"], cliff_after=p2[c]["cliff"],
                                                 ev_scale_after=p2[c]["scale"], rho_before=p1[c]["rho"], rho_after=p2[c]["rho"])
                                        for c in p1}
            print(f"\n{pre}re-level check (b1 -> b2):")
            for c, v in S[f"{pre}relevel_check"].items():
                print(f"  {c}: scale asked {v['scale_asked']:.4f}  B ratio {v['B_ratio']:.5f} vs s^2 {v['scale_sq']:.5f} "
                      f"(ratio {v['B_ratio_over_scale_sq']:.5f})  rho {v['rho_before']:.4f} -> {v['rho_after']:.4f}  "
                      f"cliff after {v['cliff_after']:.4f}  EV-scale after {v['ev_scale_after']:.4f}")
    (W / "summary432.json").write_text(json.dumps(S, indent=1, default=float))
    print("| work dir | stage | seconds | peak RSS GiB | rc |\n|---|---|---:|---:|---:|")
    for name, row in S.items():
        for st, v in row.get("driver", {}).items():
            print(f"| {name} | driver:{st} | {v['seconds']} | {v['peak_rss_gib']} | {v['returncode']} |")
        for st in ("repin-panel-support", "build", "prep"):
            if st in row:
                print(f"| {name} | {st} | {row[st]['seconds']} | {row[st]['peak_rss_gib']} | |")
    for name, row in S.items():
        if "rc_dump" in row:
            print(f"\n{name} rc dump ({row['rc_dump']['blocks']} blocks):")
            for ctx, d in row["rc_dump"].items():
                if ctx == "blocks":
                    continue
                print(f"  {ctx}: " + "; ".join(f"{k} med {v['median']:.4f} q10 {v['q10']:.4f} min {v['min']:.4f} "
                                             f"frac0 {v['frac_zero']:.3f}" for k, v in d.items()
                                             if k in ("mde_lambda", "zero_repaired", "zs_n_moved")))
        if "ntclevel" in row:
            print(f"\n{name} ntclevel: " + "; ".join(f"{c}: Cp46 {v['Cp46']:.4f} rho {v['rho']:.4f} cliff {v['cliff']:.4f} "
                                                  f"scale {v['scale']:.4f}" for c, v in row["ntclevel"].items()))


if __name__ == "__main__":
    main()
