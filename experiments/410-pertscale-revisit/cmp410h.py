import json, sys
from pathlib import Path
W = Path(__file__).resolve().parent / "work" / "h1"
base = json.loads((W / "metrics-base395.json").read_text())["raw"]
SL = {"pds": 2.2138, "mse": -1.051, "nmae": -1.6532, "fid": 3.3645, "reach": 1.1292, "jac": 2.7263}  # board raw->scaled slopes
print(f"{'arm':16s} " + " ".join(f"{m:>9s}" for m in SL) + "   board-slope overall")
print(f"{'base395':16s} " + " ".join(f"{base[m]:9.5f}" for m in SL))
for f in sorted(W.glob("metrics-*.json")):
    arm = f.stem[len("metrics-"):]
    if arm == "base395": continue
    raw = json.loads(f.read_text())["raw"]
    d = {m: raw[m] - base[m] for m in SL}
    ov = sum(SL[m] * d[m] for m in SL) / 6
    print(f"{arm:16s} " + " ".join(f"{d[m]:+9.5f}" for m in SL) + f"   {ov:+.5f}")
