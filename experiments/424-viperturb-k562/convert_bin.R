# 424: one VIPerturb-seq genome-wide bin (Seurat v5 RDS, Zenodo 18460279) -> per-target count sums.
#
# Reads the object WITHOUT SeuratObject (S4 slots via attr()), takes the RNA assay's raw `counts` layer(s) and the
# per-cell `gene` label from meta.data, and writes, for every label (targets and NO-TARGET):
#   sums.f64       (n_groups x n_features) float64, row-major: summed raw UMI counts over the label's cells
#   groups.txt     label per row of sums;  ncells.txt  cells per label
#   features.txt   feature (gene symbol) per column
#   cells.tsv      per cell: barcode, gene, guide, sample, nCount_RNA (meta), libsize (colSums of counts),
#                  target_count (counts of the cell's own target gene; NA when the target is not a feature)
#   info.json      shapes, layer names, meta columns, timings
# Group sums are all the donor table needs (b = log1p(TS * P / sum P) of summed counts, scripts/build_donor_table.py).
#   Rscript convert_bin.R <in.RDS> <outdir>
suppressMessages(library(Matrix))
args <- commandArgs(TRUE)
inp <- args[1]; out <- args[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
el <- function() round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
obj <- readRDS(inp)
cat(sprintf("read %s in %ds; class %s\n", inp, el(), paste(class(obj), collapse = ",")))
cat("object attributes:", paste(names(attributes(obj)), collapse = ", "), "\n")
assays <- attr(obj, "assays")
cat("assays:", paste(names(assays), collapse = ", "), "\n")
rna <- assays[["RNA"]]
cat("RNA assay class:", class(rna), " attributes:", paste(names(attributes(rna)), collapse = ", "), "\n")
layers <- attr(rna, "layers")
cat("layers:", paste(names(layers), collapse = ", "), "\n")
meta <- attr(obj, "meta.data")
cat(sprintf("meta.data: %d x %d; columns %s\n", nrow(meta), ncol(meta), paste(colnames(meta), collapse = ", ")))
feat_lm <- attr(rna, "features"); cell_lm <- attr(rna, "cells")
features <- rownames(feat_lm); cells_all <- rownames(cell_lm)
cat(sprintf("features %d, cells %d; LogMap feature cols %s; cell cols %s\n", length(features), length(cells_all),
            paste(colnames(feat_lm), collapse = ","), paste(colnames(cell_lm), collapse = ",")))
cnames <- grep("^counts", names(layers), value = TRUE)
stopifnot(length(cnames) >= 1)
cat("counts layers:", paste(cnames, collapse = ", "), "\n")
rm(obj); invisible(gc())

groups <- NULL; sums <- NULL; ncell <- NULL
cells_tab <- list()
for (ln in cnames) {
  cnt <- layers[[ln]]
  # the layer's own features / cells from the LogMaps (Seurat v5 keeps layer dimnames there)
  fsel <- if (ln %in% colnames(feat_lm)) features[as.logical(feat_lm[, ln])] else features
  csel <- if (ln %in% colnames(cell_lm)) cells_all[as.logical(cell_lm[, ln])] else cells_all
  if (!is.null(rownames(cnt))) fsel <- rownames(cnt)
  if (!is.null(colnames(cnt))) csel <- colnames(cnt)
  stopifnot(nrow(cnt) == length(fsel), ncol(cnt) == length(csel), identical(fsel, features))
  cat(sprintf("layer %s: %d x %d, class %s, nnz %.0f, x storage %s, all-integer %s (%ds)\n", ln, nrow(cnt), ncol(cnt),
              class(cnt), as.numeric(length(cnt@x)), typeof(cnt@x), all(cnt@x == round(cnt@x)), el()))
  m <- meta[csel, , drop = FALSE]
  stopifnot(identical(rownames(m), csel))
  g <- as.character(m$gene)
  lev <- sort(unique(g)); code <- match(g, lev)
  ind <- sparseMatrix(i = seq_along(code), j = code, x = 1, dims = c(length(code), length(lev)))
  S <- as.matrix(cnt %*% ind)                   # features x groups
  cat(sprintf("group sums done: %d groups (%ds)\n", length(lev), el()))
  lib <- colSums(cnt)
  ti <- match(g, features)
  tc <- rep(NA_real_, length(g))
  ok <- which(!is.na(ti))
  if (length(ok)) tc[ok] <- cnt[cbind(ti[ok], ok)]
  cells_tab[[ln]] <- data.frame(barcode = csel, gene = g,
                                guide = if ("guide" %in% colnames(m)) as.character(m$guide) else NA,
                                sample = if ("sample" %in% colnames(m)) as.character(m$sample) else NA,
                                nCount_RNA = if ("nCount_RNA" %in% colnames(m)) m$nCount_RNA else NA,
                                libsize = lib, target_count = tc, stringsAsFactors = FALSE)
  nc <- as.numeric(table(factor(g, levels = lev)))
  if (is.null(groups)) {
    groups <- lev; sums <- S; ncell <- nc
  } else {                                        # several counts layers: merge by label
    all <- sort(union(groups, lev))
    S2 <- matrix(0, nrow(S), length(all)); S2[, match(groups, all)] <- sums; S2[, match(lev, all)] <- S2[, match(lev, all)] + S
    n2 <- numeric(length(all)); n2[match(groups, all)] <- ncell; n2[match(lev, all)] <- n2[match(lev, all)] + nc
    groups <- all; sums <- S2; ncell <- n2
  }
  rm(cnt, S, ind); invisible(gc())
}
con <- file(file.path(out, "sums.f64"), "wb"); writeBin(as.vector(sums), con, size = 8); close(con)   # column-major f x g
writeLines(groups, file.path(out, "groups.txt"))
writeLines(as.character(ncell), file.path(out, "ncells.txt"))
writeLines(features, file.path(out, "features.txt"))
ct <- do.call(rbind, cells_tab)
write.table(ct, file.path(out, "cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
info <- sprintf('{"input": "%s", "n_features": %d, "n_groups": %d, "n_cells": %d, "layers": "%s", "meta_columns": "%s", "seconds": %d, "layout": "sums.f64 = float64, numpy shape (n_groups, n_features) C-order"}',
                basename(inp), length(features), length(groups), nrow(ct), paste(cnames, collapse = ","),
                paste(colnames(meta), collapse = ","), el())
writeLines(info, file.path(out, "info.json"))
cat(sprintf("wrote %s (%ds)\n", out, el()))
