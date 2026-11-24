"""181 — fetch Ensembl genomic spans for the h1-test DE-universe gene axis, so the
'long-gene / nuclear-fraction technical axis' identification (175 addendum, AGENTS.md
harness warning) can be tested as a correlation instead of a top-30 eyeball."""
import json, time, urllib.request
import numpy as np, pandas as pd

S = np.load('experiments/177-resample-emitter/de177.KA_a15_raw.npz', allow_pickle=True)
genes = [str(g) for g in S['genes']]
meta = pd.read_parquet('data/external/orion/metadata/gene_metadata.parquet')
name2ens = dict(zip(meta.gene_name, meta.ensembl_id))
ids = [name2ens[g] for g in genes if g in name2ens]
print(f'{len(ids)} / {len(genes)} genes map to an ensembl_id in the Orion metadata', flush=True)

out = {}
URL = 'https://rest.ensembl.org/lookup/id'
for i in range(0, len(ids), 900):
    chunk = ids[i:i+900]
    req = urllib.request.Request(URL, data=json.dumps({'ids': chunk}).encode(),
                                 headers={'Content-Type': 'application/json'})
    r = {}
    for attempt in range(4):
        try:
            r = json.load(urllib.request.urlopen(req, timeout=180)); break
        except Exception as e:
            print('retry', i, repr(e)[:80], flush=True); time.sleep(5)
    for k, v in r.items():
        if v and 'start' in v:
            out[k] = dict(span=int(v['end'])-int(v['start'])+1, biotype=v.get('biotype'),
                          chrom=str(v.get('seq_region_name')))
    print(f'  {i+len(chunk)} / {len(ids)}  resolved {len(out)}', flush=True)

rows = [dict(gene=g, ensembl=name2ens[g], **out[name2ens[g]]) for g in genes
        if g in name2ens and name2ens[g] in out]
pd.DataFrame(rows).to_csv('experiments/181-audit-0902/work/gene_spans.csv', index=False)
print('wrote', len(rows), 'rows')
