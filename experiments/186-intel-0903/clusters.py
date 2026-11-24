import json, numpy as np
new = json.load(open('board/leaderboard_2026-09-03.json'))['entries']
E = sorted(new, key=lambda e: -(e['score_avg'] or -9))[:150]
keys = ['pds_cosine','expr_mse_unbiased_capped_norm','de_wilcoxon_lfc_nmae','de_wilcoxon_direction_fidelity_yield_raw','de_wilcoxon_direction_reach_raw','de_wilcoxon_sig_jaccard']
X = np.array([[e[k] for k in keys] for e in E])
# relative distance: normalize each col by its std within top150
sd = X.std(0); Z = X/sd
print("Near-identical pairs among top 150 (all six members within 0.05 sd each):")
seen=set()
for i in range(len(E)):
    for j in range(i+1,len(E)):
        if np.all(np.abs(Z[i]-Z[j])<0.05):
            a,b=E[i],E[j]
            print(f"  {a['team_name'][:18]:<18} `{str(a['model_name'])[:22]:<22}` <-> {b['team_name'][:18]:<18} `{str(b['model_name'])[:22]:<22}` pds {a['pds_cosine']:.4f}/{b['pds_cosine']:.4f} reach {a['de_wilcoxon_direction_reach_raw']:.4f}/{b['de_wilcoxon_direction_reach_raw']:.4f} mse {a['expr_mse_unbiased_capped_norm']:.3f}/{b['expr_mse_unbiased_capped_norm']:.3f} dates {a['submission_date'][:13]} {b['submission_date'][:13]}")
print("\nLoose clusters (within 0.15 sd on all six):")
for i in range(len(E)):
    for j in range(i+1,len(E)):
        if np.all(np.abs(Z[i]-Z[j])<0.15) and not np.all(np.abs(Z[i]-Z[j])<0.05):
            a,b=E[i],E[j]
            print(f"  {a['team_name'][:18]:<18} `{str(a['model_name'])[:22]:<22}` <-> {b['team_name'][:18]:<18} `{str(b['model_name'])[:22]:<22}` pds {a['pds_cosine']:.4f}/{b['pds_cosine']:.4f} reach {a['de_wilcoxon_direction_reach_raw']:.4f}/{b['de_wilcoxon_direction_reach_raw']:.4f} mse {a['expr_mse_unbiased_capped_norm']:.3f}/{b['expr_mse_unbiased_capped_norm']:.3f}")
# reach vs pds regression on pds>0.70
allE=[e for e in new if e['pds_cosine']>0.70]
p=np.array([e['pds_cosine'] for e in allE]); r=np.array([e['de_wilcoxon_direction_reach_raw'] for e in allE])
A=np.vstack([p,np.ones_like(p)]).T; coef,res,_,_=np.linalg.lstsq(A,r,rcond=None)
pred=A@coef; resid=r-pred; s=resid.std()
print(f"\nreach ~ pds on n={len(allE)} (pds>0.70): slope {coef[0]:.3f} intercept {coef[1]:.3f} resid sd {s:.4f}")
rows=[(e['team_name'],e['model_name'],e['pds_cosine'],e['de_wilcoxon_direction_reach_raw'],(e['de_wilcoxon_direction_reach_raw']-(coef[0]*e['pds_cosine']+coef[1]))/s,e['expr_mse_unbiased_capped_norm'],e['score_avg']) for e in allE]
rows.sort(key=lambda t:-t[4])
print("Top 25 positive reach residuals (z):")
for t in rows[:25]: print(f"  {t[0][:20]:<20} `{str(t[1])[:24]:<24}` pds {t[2]:.4f} reach {t[3]:.4f} z {t[4]:+.2f} mse {t[5]:.2f} score {t[6]:.4f}")
print("Ours:", [t for t in rows if t[0]=='Ahmet Dedeler'])
