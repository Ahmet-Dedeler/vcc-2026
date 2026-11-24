"""181 — audit of the 096/170/180 `mse` identity against the board's own alpha ladder.

The per-context knee rule (180) and 170's rho-cliff both rest on
    raw_mse(s) = 1 + mean_c[ s^2 A_c + Kc_c - 2 s X - min(s^2 B_c, claim_c) ] / sumD
with claim_c capped by the NTC-proxied C_p.  We have THREE board points that differ
mainly in VCC_ALPHA.  If the identity is right, it must reproduce their raw mse.

Known-answer rows: (i) a15/a23's own raw mse must come out at 0.8723-0.8738;
(ii) the per-context knees must reproduce 180's 0.9823 / 1.0129 / 1.0479.
"""
import numpy as np

# --- constants, all from primary sources
SUMD = 0.229141            # 096 board-solved
CLAMP_B, CLAMP_S = 0.987862879, 0.950162566   # 096 (d)
S_AT_982 = 0.90960         # 164 S5 / 171 S5: VCC_ALPHA 0.982 -> s
K_ALPHA = S_AT_982/0.982

# 180 S5 table, measured on a23's effect matrix at VCC_ALPHA 0.982 (s=0.9096)
CTX = {
 'A': dict(A=0.234830, B=0.233593, claim=0.233754, Kc=0.005084),
 'B': dict(A=0.224778, B=0.223624, claim=0.237923, Kc=0.005228),
 'C': dict(A=0.234941, B=0.233728, claim=0.266129, Kc=0.005919),
}
X_AT_982 = 0.017932        # 180 S6 / 170 gate 8

def raw_mse(alpha, Xraw=None, sumd=SUMD, claim_scale=1.0):
    s = K_ALPHA*alpha
    Xr = X_AT_982/S_AT_982 if Xraw is None else Xraw
    nums=[]
    for c,v in CTX.items():
        Ar = v['A']/S_AT_982**2; Br = v['B']/S_AT_982**2; cl = v['claim']*claim_scale
        nums.append(s*s*Ar - 2*s*Xr + v['Kc'] - min(s*s*Br, cl))
    return 1 + np.mean(nums)/sumd

def knees(claim_scale=1.0):
    out={}
    for c,v in CTX.items():
        Br=v['B']/S_AT_982**2
        out[c] = np.sqrt(v['claim']*claim_scale/Br)/K_ALPHA
    return out

print('KA1  per-context knee alphas :', {k:round(v,5) for k,v in knees().items()}, ' (180: 0.98234 / 1.01291 / 1.04786)')
print('KA2  raw mse at alpha 0.982  :', round(raw_mse(0.982),6), ' (a23 board raw 0.872343, a15 0.873794)')
print()

# --- the board's own alpha ladder
BOARD = [   # (name, alpha, score_mse, note)
 ('a11', 0.838, 0.082041, 'a10 geometry + alpha 1.0796->0.838'),
 ('a14', 0.982, 0.110615, 'a12 alpha 1.0108->0.982, pure alpha pair'),
 ('a12', 1.0108, 0.091163, 'a11 alpha 0.838->1.0108 + tier-8000 spread'),
 ('a15', 0.982, 0.120057, 'a14 + CD4 f=0.10'),
 ('a23', 0.982, 0.121581, 'a15 + minnorm mass fix'),
]
print(f"{'arm':5} {'alpha':>7} {'board raw':>10} {'model raw':>10} {'resid':>9}")
for n,a,sm,note in BOARD:
    raw = CLAMP_B - sm*CLAMP_S
    m = raw_mse(a)
    print(f'{n:5} {a:7.4f} {raw:10.6f} {m:10.6f} {m-raw:+9.6f}   {note}')

print()
print('--- the shape of the modelled curve vs the three board alphas ---')
for a in [0.80,0.838,0.86,0.90,0.94,0.96,0.982,1.00,1.0108,1.03,1.05,1.08]:
    print(f'  alpha {a:6.4f}  raw {raw_mse(a):.6f}   score {max(0,(CLAMP_B-raw_mse(a))/CLAMP_S):.6f}')
amin = min(np.arange(0.80,1.10,0.0002), key=raw_mse)
print(f'  modelled argmin alpha = {amin:.4f}  raw {raw_mse(amin):.6f}')

# --- sensitivity: what if the 400-cell NTC proxy over-states C_p because the
#     reference groups hold ~461 cells (manifest: 138,400 gt cells / 300 perts)?
print()
print('--- sensitivity of the knee to a biased C_p proxy (claim scales with C_p where the cap binds) ---')
for f in [0.80, 0.868, 0.90, 1.00, 1.10, 1.20]:
    k = knees(f); am = min(np.arange(0.70,1.30,0.0005), key=lambda a: raw_mse(a, claim_scale=f))
    print(f'  C_p x {f:4.3f} -> knees {k["A"]:.4f}/{k["B"]:.4f}/{k["C"]:.4f}   global argmin alpha {am:.4f}')

# --- sensitivity of the KNEE GAIN to the unmeasurable per-context X_c and sumD_c
print()
print('--- knee gain (180: +0.000838 overall) under per-context X / sumD uncertainty ---')
def knee_gain(Xmult=1.0, sumd=SUMD):
    tot=0.0
    Xr = X_AT_982/S_AT_982*Xmult
    for c,v in CTX.items():
        Ar=v['A']/S_AT_982**2; Br=v['B']/S_AT_982**2
        s0=S_AT_982; s1=np.sqrt(v['claim']/Br)
        n0=s0*s0*Ar-2*s0*Xr+v['Kc']-min(s0*s0*Br,v['claim'])
        n1=s1*s1*Ar-2*s1*Xr+v['Kc']-min(s1*s1*Br,v['claim'])
        tot += ((n0-n1)/sumd)/CLAMP_S
    return tot/3/6
print(f'  X x1.00, sumD 0.229141 : {knee_gain():+.6f}   (180 reports +0.000838)')
for xm in [0.5,0.75,1.0,1.5,2.0]:
    print(f'  X x{xm:4.2f}                : {knee_gain(xm):+.6f}')
for sd in [0.19,0.229141,0.27]:
    print(f'  sumD {sd:.5f}            : {knee_gain(sumd=sd):+.6f}')
