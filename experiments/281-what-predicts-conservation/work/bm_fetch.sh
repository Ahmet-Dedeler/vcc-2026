#!/bin/zsh
# Ensembl BioMart (useast mirror, GRCh38 current release) — homologs page only.
# Terms: Ensembl data are released under the Apache 2.0 licence / no restriction on use (https://www.ensembl.org/info/about/legal/disclaimer.html)
H=https://useast.ensembl.org/biomart/martservice
q() {  # $1 = filter block, $2 = attribute block, $3 = out
  Q="<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query virtualSchemaName=\"default\" formatter=\"TSV\" header=\"1\" uniqueRows=\"1\" count=\"\" datasetConfigVersion=\"0.6\"><Dataset name=\"hsapiens_gene_ensembl\" interface=\"default\">$1$2</Dataset></Query>"
  curl -s --max-time 900 -G --data-urlencode "query=$Q" "$H" -o "$3"
  echo "$3 $(wc -l < $3) lines"
}
for sp in scerevisiae celegans dmelanogaster drerio mmusculus; do
  q "<Filter name=\"with_${sp}_homolog\" excluded=\"0\"/>" "<Attribute name=\"ensembl_gene_id\"/>" "bm_orth_${sp}.tsv"
done
q "<Filter name=\"with_hsapiens_paralog\" excluded=\"0\"/>" "<Attribute name=\"ensembl_gene_id\"/><Attribute name=\"hsapiens_paralog_associated_gene_name\"/>" "bm_paralog.tsv"
echo DONE
