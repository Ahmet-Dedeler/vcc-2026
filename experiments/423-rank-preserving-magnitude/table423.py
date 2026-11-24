"""RESULTS table rows from work/analysis423[_bps].json (after analyse423.py)."""
import json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sfx = sys.argv[1] if len(sys.argv) > 1 else ""
rows = json.loads((HERE / f"work/analysis423{sfx}.json").read_text())
order = ["all_k080", "all_k120", "all_k140", "all_k160", "all_k200", "q50_k140", "q80_k140", "q80_k200", "up_k140", "dn_k140",
         "dn_k160", "dn_k180", "dn_k200", "dn_k250", "dn_k300", "dn_k400", "q80dn_k200", "q80dn_k250",
         "dn_k180_zspc", "dn_k200_zspc", "dn_k250_zspc",
         "all_k140_zs", "dn_k200_zs", "dn_k250_zs", "dn_k300_zs", "dn_k500_zs",
         "all_k140_ex", "all_k200u120_ex", "dn_k200_ex", "dn_k250_ex"]
by = {r["variant"]: r for r in rows}
mode = lambda v: v.split("_")[-1] if v.split("_")[-1] in ("zs", "ex", "zspc") else "pc"
print("| variant | mode | seeds | Δnmae (± sd over seeds) | Δpds | Δmse | Δfid | Δreach | Δjac | total → overall | feasible (\\|dA\\|-w) | genes moved / target | genes with ΔU / target |")
print("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
for v in order + [v for v in by if v not in order]:
    if v not in by:
        continue
    r = by[v]; n = len(r["seeds"])
    sd = f" ± {r['sd_nmae']:.5f}" if r["sd_nmae"] is not None else ""
    b = "**" if v in ("dn_k250_zs",) else ""
    f = lambda x, d=5: "0" if abs(x) < 10 ** -(d + 1) else f"{x:+.{d}f}"
    print(f"| {b}{v}{b} | {mode(v)} | {n} | {b}{r['d_nmae']:+.5f}{sd}{b} | {f(r['d_pds'])} | {f(r['d_mse'], 6)} | {f(r['d_fid'], 6)} | {f(r['d_reach'])} | "
          f"{f(r['d_jac'], 6)} | {r['ovr']:+.5f} | {r['frac_massw_mean']:.3f} | {r['n_moved_mean']:.0f} | {r['n_genes_U_changed_mean']:.0f} |")
