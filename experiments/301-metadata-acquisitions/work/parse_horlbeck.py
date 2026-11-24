import openpyxl, pickle, sys
wb=openpyxl.load_workbook('horlbeck2016_supp3_hCRISPRi-v2.xlsx', read_only=True)
print('sheets', wb.sheetnames, flush=True)
out={}
for sn in wb.sheetnames:
    ws=wb[sn]
    rows=ws.iter_rows(values_only=True)
    hdr=None; n=0; recs=[]
    for r in rows:
        if hdr is None:
            hdr=r; print(sn,'header',hdr, flush=True); continue
        recs.append(r); n+=1
    print(sn, n, 'rows', flush=True)
    out[sn]=(hdr,recs)
pickle.dump(out, open('horlbeck_supp3.pkl','wb'))
print('done', flush=True)
