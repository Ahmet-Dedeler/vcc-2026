"""425: gene TSS table from GENCODE v44 basic (gene records), for CRISPRi cis-neighbour analysis.

Writes work/tss_v44.tsv: gene_name, gene_id, chrom, strand, tss (1-based gene start on +, gene end on -), gene_type.
Duplicate names (PAR copies, readthrough clones) keep the first protein-coding record, else the first record.
"""
import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GTF = ROOT / "experiments/384-promoterai-gene-priors/data/gencode.v44.basic.annotation.gtf.gz"
OUT = Path(__file__).resolve().parent / "work/tss_v44.tsv"

rows = {}
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] != "gene":
            continue
        attrs = dict((kv.strip().split(" ", 1)[0], kv.strip().split(" ", 1)[1].strip('"'))
                     for kv in f[8].split(";") if kv.strip())
        name, gtype = attrs.get("gene_name"), attrs.get("gene_type")
        tss = int(f[3]) if f[6] == "+" else int(f[4])
        rec = (name, attrs.get("gene_id"), f[0], f[6], tss, gtype)
        if name not in rows or (rows[name][5] != "protein_coding" and gtype == "protein_coding"):
            rows[name] = rec
OUT.write_text("gene_name\tgene_id\tchrom\tstrand\ttss\tgene_type\n" +
               "".join("\t".join(map(str, r)) + "\n" for r in rows.values()))
print(len(rows), "genes ->", OUT)
