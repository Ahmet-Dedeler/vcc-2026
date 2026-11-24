"""643: print the per-construct table and the preregistered verdict from work/k562_f16_643.json."""
import json
from pathlib import Path
d = json.loads((Path(__file__).parent / "work/k562_f16_643.json").read_text())
r = d["result"]
print(json.dumps(r, indent=1))
print(f"{'construct':14s} {'partner':14s} {'n':>2s} {'cells':>5s} {'p.cells':>7s} {'own':>8s} {'partner':>8s} {'D':>8s} {'nullD':>8s} {'K562kd':>7s} {'H1kd':>6s} {'H1 drep':>8s}")
for x in sorted(d["rows"], key=lambda x: x.get("D", 9)):
    if not x.get("n_cand"):
        print(f"{x['construct']:14s} {x['partner']:14s} no expressed candidates")
        continue
    print(f"{x['construct']:14s} {x['partner']:14s} {x['n_cand']:2d} {x['k562_cells']:5d} {x['partner_cells']:7d} {x['own']:+8.4f} {x['partner_val']:+8.4f} "
          f"{x['D']:+8.4f} {x['null_D_mean']:+8.4f} {x['k562_kd_own'] if x['k562_kd_own'] is None else round(x['k562_kd_own'],3):>7} {x['h1_kd']:+6.2f} {x['h1_delta_repaired']:+8.4f}")
