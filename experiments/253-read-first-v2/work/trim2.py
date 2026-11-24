p = "experiments/253-read-first-v2/READ-FIRST.md"
s = open(p, encoding="utf-8").read()
edits = [
 (" 4th 0.21092928 at 15:33Z (LB:86-87).", ""),
 ("; none beside a build, nor `vcc prep` (RB:475-476)", ""),
 ("On trial in the 237 slot (§6); 243 does not price `nmae` (R243:292)", "tested by the 237 slot (§6); 243 does not price it (R243:292)"),
 ("254 CD4 partial rows and the composed `cd4svd5+partial` build, no submit", "254 partial rows + composed `cd4svd5+partial` build, no submit"),
 ("3/4 closed by 088/213/164 and 232/234, 5 worthless, 2 open-low", "3/4 closed, 5 worthless, 2 open-low"),
 (" (ids: SL:1-77):", " (ids: SL):"),
 ("235 `axis_sha256`, `REFUSING TO BUILD`, A/B/C axis", "235 `axis_sha256` guard, A/B/C axis"),
 ("Prereg before the build: ", "Prereg: "),
 ("`--emitter rc_callset` (RB:147,208)", "`--emitter rc_callset` (RB:208)"),
 (" Trap −0.0429 on `pds` alone (RB:255).", ""),
 ("9/10 statistics bit-exact, 10th 3.3e-16 (R221:8)", "9/10 bit-exact, 10th 3.3e-16 (R221:8)"),
]
for old, new in edits:
    assert s.count(old) == 1, old
    s = s.replace(old, new)
open(p, "w", encoding="utf-8").write(s)
print(len(s.encode("utf-8")), "bytes")
