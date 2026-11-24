"""656 reader (PREREGISTRATION 28d99181): 450's board units (read450 SLOPE on harness raw members at f 1, vs 435's base, paired
by seed) for w 0.5 / 1 / 2 x p 0.2 / 0.3 with zero repair; Delta(w) = overall(w) - overall(0.5) per seed, mean and sd; the
decision rule; per-arm zero-total counts. usage: uv run python experiments/656-sameline-w-zr/read656.py"""
import glob, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L  # noqa: E402
SLOPE = {"pds": 2.2138, "mse": -1.051, "nmae": -1.6532, "fid": 3.3645, "reach": 1.1292, "jac": 2.7263}   # read450's
BASE = ROOT / "experiments/435-zs-on-coarse/work/h1"
TAG = "base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse"
S = L.load_state(0); OWN = S["own"]; P35 = L.pools(S, 0.35) & ~OWN; P1 = L.pools(S, 1.0) & ~OWN
out, res = [], {}
say = lambda s="": (print(s), out.append(s))  # noqa: E731


def overall(metrics, base):
    return sum(SLOPE[m] * (metrics[m] - base[m]) for m in SLOPE) / 6.0


for p in ("0.2", "0.3"):
    say(f"\n== bank flip p {p} (gate signs ~{0.66 if p == '0.2' else 0.61}); board units vs 435's base")
    ov = {}
    for w in ("0.5", "1", "2"):
        d = HERE / f"work/h1c_w{w}_clean4_f{p}s0"
        vals, zt, raws = [], [], {}
        for s in (0, 1, 2):
            m = glob.glob(str(d / f"metrics-base_s{s}_bps*zr-5*.json"))
            st = glob.glob(str(d / f"stats-base_s{s}_bps*zr-5*.npz"))
            b = BASE / f"metrics-{TAG.format(s=s)}.json"
            if len(m) == 1:
                raws[s] = json.loads(Path(m[0]).read_text())["raw"]
                if b.exists():
                    vals.append(overall(raws[s], json.loads(b.read_text())["raw"]))
            if len(st) == 1:
                l = np.load(st[0])["lfc"]; z = (l < -30) & ~OWN
                zt.append((int(z.sum()), int((z & P1).sum()), int((z & P35).sum())))
        ov[w] = raws
        say(f"  w {w:>3s}: overall {np.mean(vals) if vals else float('nan'):+.5f} seeds {[round(v, 5) for v in vals]}   "
            f"non-own zero totals (all, in f1, in f0.35) per seed {zt}")
        res[f"{p}|{w}"] = dict(seeds=vals, zero_totals=zt)
    for w in ("1", "2"):
        common = sorted(set(ov[w]) & set(ov["0.5"]))                 # base-free: base cancels in Delta(w)
        dd = [overall(ov[w][s], ov["0.5"][s]) for s in common]
        n = len(dd)
        if dd:
            mu, sd = float(np.mean(dd)), (float(np.std(dd, ddof=1)) if n > 1 else float("nan"))
            res[f"{p}|delta{w}"] = dict(mean=mu, sd=sd, seeds=dd)
            say(f"  Delta(w {w} - w 0.5): mean {mu:+.5f}, per-seed sd {sd:.5f}, seeds {[round(x, 5) for x in dd]}  "
                f"-> {'BEATS 0.5 by > 2 sd' if (n > 1 and mu > 0 and mu > 2 * sd) else 'does not beat 0.5'}")
q = lambda p, w: res.get(f"{p}|delta{w}", {})  # noqa: E731
beats = {w: all(q(p, w).get("mean", -1) > 0 and q(p, w).get("mean", -1) > 2 * q(p, w).get("sd", np.inf) for p in ("0.2", "0.3")) for w in ("1", "2")}
p2only = {w: (q("0.2", w).get("mean", -1) > 2 * q("0.2", w).get("sd", np.inf)) and not beats[w] for w in ("1", "2")}
if any(len(q(p, w).get("seeds", [])) < 3 for p in ("0.2", "0.3") for w in ("1", "2")):
    verdict = "INCOMPLETE (fewer than 3 paired seeds for some arm): no verdict"
elif any(beats.values()):
    wstar = max([w for w in beats if beats[w]], key=lambda w: q("0.3", w)["mean"])
    verdict = f"OCTOBER VALUE MOVES: sameline_w 0.5 -> {wstar} (beats 0.5 at both p by > 2 sd)"
elif any(p2only.values()):
    verdict = "RECORDED: a higher w beats 0.5 at p 0.2 only (a per-context w needs builder support); October stays 0.5"
else:
    verdict = "0.5 STANDS (no w beats it at both p by > 2 per-seed sd), now with zero repair and 3 seeds"
say(f"\nVERDICT: {verdict}")
res["verdict"] = verdict
(HERE / "work/read656.json").write_text(json.dumps(res, indent=1, default=float))
(HERE / "work/read656.out").write_text("\n".join(out) + "\n")
