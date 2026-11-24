import sys
p = "experiments/253-read-first-v2/READ-FIRST.md"
s = open(p, encoding="utf-8").read()
edits = [
 ("(board fetched 16:53Z)", "(16:53Z)"),
 ("; entry ids and package shas in SL:1-77:", " (ids: SL:1-77):"),
 ("; donors 1:1:1 (RB:216)", ""),
 (" Knee: `mse` −0.00226 (BR180:32).", ""),
 (" a23 mass fix −0.00073 (RB:213); a21 resample `reach` −0.0185 (BL:30).", ""),
 ("never the floor, which mis-reads 1% amplitude as 0.7–1.7% cosine (R243:289-290)", "never the floor (R243:289-290)"),
 ("a21 2/4 DE signs wrong, +0.050 vs −0.0185 (A:556-559)", "a21 2/4 DE signs wrong (A:556-559)"),
 ("unpriced/priced 7.4× (R162:17). PLOLO→board", "PLOLO→board"),
 (", on `cos_agg`, a23 tag, 3 cancer folds, per-fold mean, measured cosine (R243:44-47,286-290)", ", on `cos_agg`, 3 cancer folds, measured cosine (R243:44-47,286-290)"),
 ("hct116 depth-confounded, 093 +0.0202 → −0.00291 (A:592-598)", "hct116 depth-confounded (A:592-598)"),
 (" as `Reviewer 2 Was Right About Everything` (E237", " (E237"),
 ("243 predicts +0.0022…+0.0038; the −0.004783 `nmae` leg is the whole gap to 229 (AD:9-16). Board raw `mse` bins: ≤ 0.80768 above every model; –0.81002 dose; –0.81135 coverage; –0.81392 new-source (point 0.81244); ≥ 0.81392 below every model (AD:23-29).",
  "243: +0.0022…+0.0038, differing from 229 by the −0.004783 `nmae` leg (AD:9-16). Raw `mse` bins: ≤0.80768 above all models; –0.81002 dose; –0.81135 coverage; –0.81392 new-source, point 0.81244; ≥0.81392 below all (AD:23-29)."),
 (" as `The Control Group Has Opinions` after", " after"),
 (", raw `mse` ≈ 0.8183 (P249:27-40)", " (P249:27-40)"),
 ("244's rehearsal on an essential-free pessimistic bundle: 1,692.7 s of stages, 29.7 min to a prepped artifact (RB:35-36)", "244 rehearsal: 29.7 min to a prepped artifact (RB:35-36)"),
 ("arms asserted by value; before it the driver dropped both arms at 14/14 PASS (RB:129-135; R240:23-24)", "arms asserted by value (RB:129-135)"),
 (" Still open from 242: harness rows (252), monitoring rows (R247:385-388).", ""),
 ("; `n_carried < 250` ⇒ re-fit forbidden; captured fraction outside [0.05, 0.19] unmeasured (R251:394-395)", " and its census trigger (R251:394-395)"),
 ("; CD4 union 174 (RB:432)", ""),
 ("; fitter `scripts/build_union_whitener.py` (RB:160)", ""),
 ("; basis registry MB:793-806", ""),
 ("; `202-rc_callset` is not a valid `--incumbent` (RB:96)", ""),
]
for old, new in edits:
    n = s.count(old)
    assert n == 1, (n, old)
    s = s.replace(old, new)
open(p, "w", encoding="utf-8").write(s)
print(len(s.encode("utf-8")), "bytes")
