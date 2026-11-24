p = "experiments/253-read-first-v2/READ-FIRST.md"
s = open(p, encoding="utf-8").read()
for old, new in [("0.21287742 (LB:66);\n", "0.21287742 (LB:66).\n"),
                 ("from a15 0.1740865 (BL:28) (ids: SL):", "from a15 0.1740865 (BL:28; ids in SL):"),
                 ("(BR222:49). tested by the 237 slot", "(BR222:49). Tested by the 237 slot")]:
    assert s.count(old) == 1, old
    s = s.replace(old, new)
open(p, "w", encoding="utf-8").write(s)
print(len(s.encode("utf-8")))
