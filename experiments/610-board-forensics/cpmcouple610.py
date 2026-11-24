"""610 cpmcouple: is the null-gene sign channel of fid (p_null) driven by per-cell CPM coupling to each perturbation's net
library-mass change, and can we reach it? H1 truth (323 reference) and 435's production-emitter stats (lib521.load_state).
Three parts, run in order (copied from the 21dd4c scratch scripts cpmcouple.py / cpmcouple2.py / cpmcouple3.py).
Analysis only, < 2 GB, seconds. Output: work/cpmcouple610.out."""
import runpy, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
for k in (1, 2, 3):
    print(f"\n### part {k}")
    sys.argv = [str(HERE / f"work/cpmcouple610_part{k}.py")]
    runpy.run_path(sys.argv[0], run_name="__main__")
