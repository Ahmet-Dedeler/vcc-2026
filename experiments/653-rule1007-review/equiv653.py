"""653: v3 vs v4 rule() over every bin combination (1,920 states). v4 runs with --arc-hook so only R1 can differ; reports
every differing state's class. Usage: uv run python experiments/653-rule1007-review/equiv653.py"""
import importlib.util, itertools, collections
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
def mod(p, n):
    s = importlib.util.spec_from_file_location(n, ROOT / p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
v3 = mod("experiments/537-arm-pipeline/rule1007_draft537.py", "v3")
v4 = mod("experiments/653-rule1007-review/rule1007_v4_653.py", "v4")
diff = collections.Counter(); n = 0; ex = {}
for cp, res, h, k, c, cph, k20s in itertools.product(["POS", "PARTIAL", "NULL", "NEG"], [False, True],
        ["TRANSFERS", "PARTIAL", "FAILS", "UNKNOWN"], ["POS", "MID", "NEG", "UNKNOWN"], ["OK", "FAILS"],
        ["UNKNOWN", "SHIP", "NO"], [False, True]):
    if res and cp in ("POS", "PARTIAL"):
        continue
    B = {"CP": cp, "CP_rescue": res, "H": h, "K": k, "C": c, "CPH": cph}
    a, b = v3.rule(dict(B), k20s), v4.rule(dict(B), k20s, arc_hook=True)
    n += 1
    if a[0] != b[0] or a[1] != b[1]:
        diff["10-07 lists"] += 1; ex.setdefault("10-07", (B, k20s, a, b))
    if a[2] != b[2]:
        why = "k20 scored on 10-08 (R1)" if (k20s and cph != "UNKNOWN") else ("CPH read, k20 none (R1 safe mode)" if cph != "UNKNOWN" else "other")
        diff[why] += 1; ex.setdefault(why, (B, k20s, a[2], b[2]))
print(f"{n} states; differences: {dict(diff)}")
for k, v in ex.items():
    print(" ", k, "e.g.", v)

# With --scored = the ACTUAL 10-07 picks (the heads of the 10-07 lists computed with 10-07's own knowledge, k20 not yet
# scored), v4's 10-08 list must equal v3's wherever v3 is called with 10-07's knowledge, and must differ only where v3 is
# called with K20 scored on 10-07 (R1).
agree = bad = r1 = 0; r1ex = None
for cp, res, h, k, c, cph in itertools.product(["POS", "PARTIAL", "NULL", "NEG"], [False, True],
        ["TRANSFERS", "PARTIAL", "FAILS", "UNKNOWN"], ["POS", "MID", "NEG", "UNKNOWN"], ["OK", "FAILS"], ["SHIP", "NO"]):
    if res and cp in ("POS", "PARTIAL"):
        continue
    B = {"CP": cp, "CP_rescue": res, "H": h, "K": k, "C": c, "CPH": cph}
    s1, s2, _ = v3.rule(dict(B, CPH="UNKNOWN"), False)          # what 10-07 actually took
    picks = [s1[0], s2[0]]
    v3_ok = v3.rule(dict(B), False)[2]                          # v3 called with 10-07's knowledge (k20 none)
    v4_d8 = v4.rule(dict(B), "K20" in picks, scored=picks, arc_hook=True)[2]
    agree += v3_ok == v4_d8
    bad += v3_ok != v4_d8
    if "K20" in picks:                                          # the operator passes K20's id on 10-08 (it WAS scored)
        v3_r1 = v3.rule(dict(B), True)[2]
        if v3_r1 != v4_d8:
            r1 += 1; r1ex = r1ex or (B, picks, v3_r1, v4_d8)
print(f"--scored = actual 10-07 picks: v4 == v3(10-07 knowledge) in {agree} states, differs in {bad}")
print(f"R1: v3 called with K20's id on 10-08 drops or mis-lists an arm in {r1} states, e.g. {r1ex}")
