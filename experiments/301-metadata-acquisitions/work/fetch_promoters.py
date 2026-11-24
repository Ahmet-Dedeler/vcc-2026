"""301 item 2: TSS +/- 500 bp promoter windows, GRCh38, for the 18,080 genes measured in the
VCC 2025 h5ads. Coordinates+strand from Ensembl BioMart (biomart_genes.tsv); sequence from the
Ensembl REST POST /sequence/region/human endpoint (GRCh38, current release)."""
import csv, gzip, json, sys, time, urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor

FLANK = 500
GOOD = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}

symbols = [l.strip() for l in open("~/Code/virtual-cell-challenge-2026/data/external/vcc2025/gene_names.csv") if l.strip()]
want = set(symbols)

rows = []
with open("biomart_genes.tsv") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for r in rd:
        if r["Gene name"] in want and r["Chromosome/scaffold name"] in GOOD:
            rows.append(r)
# one row per symbol: prefer protein_coding, then the longest gene
best = {}
for r in rows:
    s = r["Gene name"]
    key = (r["Gene type"] == "protein_coding", int(r["Gene end (bp)"]) - int(r["Gene start (bp)"]))
    if s not in best or key > best[s][0]:
        best[s] = (key, r)
print(f"symbols in var: {len(want)}; with a GRCh38 primary-assembly BioMart row: {len(best)}", flush=True)

regions = []
for s, (_, r) in sorted(best.items()):
    strand = int(r["Strand"])
    tss = int(r["Gene start (bp)"]) if strand == 1 else int(r["Gene end (bp)"])
    lo = max(1, tss - FLANK)
    hi = tss + FLANK
    regions.append((s, r["Gene stable ID"], r["Chromosome/scaffold name"], strand, tss, lo, hi))
print(f"regions: {len(regions)}", flush=True)

URL = "https://rest.ensembl.org/sequence/region/human?coord_system_version=GRCh38"
def fetch(batch):
    body = json.dumps({"regions": [f"{c}:{lo}..{hi}:1" for _, _, c, _, _, lo, hi in batch]}).encode()
    for attempt in range(6):
        try:
            req = urllib.request.Request(URL, data=body, headers={
                "Content-Type": "application/json", "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as resp:
                return json.loads(resp.read())
        except Exception as e:
            if attempt == 5:
                raise
            time.sleep(2 ** attempt)

B = 50
batches = [regions[i:i + B] for i in range(0, len(regions), B)]
out = gzip.open("promoters_tss1000_grch38.fa.gz", "wt")
done = 0
with ThreadPoolExecutor(max_workers=4) as ex:
    for batch, seqs in zip(batches, ex.map(fetch, batches)):
        assert len(seqs) == len(batch), (len(seqs), len(batch))
        for (sym, eid, chrom, strand, tss, lo, hi), s in zip(batch, seqs):
            seq = s["seq"].upper()
            out.write(f">{sym}|{eid}|{chrom}|{strand}|{tss}|{lo}|{hi}|{len(seq)}\n{seq}\n")
        done += len(batch)
        if done % 2000 < B:
            print(f"  {done}/{len(regions)}", flush=True)
out.close()
print("WROTE promoters_tss1000_grch38.fa.gz", flush=True)
