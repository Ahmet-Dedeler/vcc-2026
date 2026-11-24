"""Small helpers for scripts/def_day.sh, kept out of inline `python3 -c` strings (bash brace-expands `{a, b}` in them:
452's first rehearsal logged an empty rho line three times because of it).

usage:
  def_day_helpers.py rho <build.ntclevel.json>            -> "D 0.6861 E 0.6206 F 0.5911"
  def_day_helpers.py past_cliff <build.ntclevel.json>     -> exit 0 if any context's rho > 1, else 1
  def_day_helpers.py same <D=..,E=..> <D=..,E=..>         -> exit 0 if the factor strings agree within 5e-5
  def_day_helpers.py inbound <D=..,E=..,F=..>             -> exit 0 if every factor is inside the builder bound [0.7, 1.4]
  def_day_helpers.py nctx <sameline.npz>                  -> the number of contexts the same-line artifact switches on
"""
import json
import sys


def factors(s: str) -> dict[str, float]:
    return {k: float(v) for k, v in (p.split("=") for p in s.split(","))}


def main() -> int:
    cmd, *args = sys.argv[1:]
    if cmd == "rho":
        per = json.load(open(args[0]))["per_context"]
        print(" ".join(f"{c} {per[c]['rho']:.4f}" for c in sorted(per)))
        return 0
    if cmd == "past_cliff":
        per = json.load(open(args[0]))["per_context"]
        return 0 if max(v["rho"] for v in per.values()) > 1.0 else 1
    if cmd == "same":
        a, b = factors(args[0]), factors(args[1])
        return 0 if a.keys() == b.keys() and all(abs(a[k] - b[k]) < 5e-5 for k in a) else 1
    if cmd == "inbound":
        return 0 if all(0.7 <= v <= 1.4 for v in factors(args[0]).values()) else 1
    if cmd == "nctx":
        import numpy as np  # noqa: PLC0415
        print(len(np.load(args[0])["contexts"]))
        return 0
    raise SystemExit(f"unknown command {cmd!r}")


if __name__ == "__main__":
    sys.exit(main())
