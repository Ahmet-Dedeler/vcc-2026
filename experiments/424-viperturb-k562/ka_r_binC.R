# 424 known answer for the streaming parser: R's own readRDS of binC, group sums for NO-TARGET + 12 labels.
suppressMessages(library(Matrix))
obj <- readRDS("work/dl/genome_wide_binC.RDS")
rna <- attr(obj, "assays")[["RNA"]]
cnt <- attr(rna, "layers")[["counts"]]
meta <- attr(obj, "meta.data")
feats <- attr(attr(rna, "features"), "dimnames")[[1]]
cells <- attr(attr(rna, "cells"), "dimnames")[[1]]
stopifnot(identical(rownames(meta), cells), nrow(cnt) == length(feats))
rm(obj); invisible(gc())
labs <- c("NO-TARGET", "PSKH1", "PTMS", "RBM18", "SIN3B", "SMARCA5", "SND1", "STAT6", "TAF15", "TEAD4", "TMED3", "ZNF706", "SPATA20")
labs <- labs[labs %in% meta$gene]
res <- sapply(labs, function(l) Matrix::rowSums(cnt[, which(meta$gene == l), drop = FALSE]))
rownames(res) <- feats
write.csv(res, "work/conv/ka_r_binC.csv")
cat("labels:", paste(labs, collapse = ","), "\n")
cat("nnz", length(cnt@x), "sum", sum(cnt@x), "\n")
