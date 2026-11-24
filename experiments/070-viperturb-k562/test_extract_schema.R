#!/usr/bin/env Rscript

# Base-R regression for the exact manifest and cell-metadata schema observed in
# genome_wide_filtered.rds. It deliberately evaluates only pure helper
# definitions from the extractor, so Seurat and the 3.6 GB RDS are not needed.

file_arguments <- commandArgs(trailingOnly = FALSE)
file_arguments <- file_arguments[grepl("^--file=", file_arguments)]
if (length(file_arguments) != 1L) {
  stop("cannot identify test script path", call. = FALSE)
}
test_path <- normalizePath(sub("^--file=", "", file_arguments[[1L]]))
source_path <- file.path(dirname(test_path), "extract_viperturb.R")
expressions <- parse(source_path)
environment <- new.env(parent = baseenv())
wanted <- c(
  "EXPECTED_SOURCE_PERTURBATIONS",
  "is_control",
  "source_catalogue_from_manifest",
  "map_target_values",
  "choose_target_column"
)
for (expression in expressions) {
  if (
    is.call(expression) &&
      identical(expression[[1L]], as.name("<-")) &&
      is.symbol(expression[[2L]]) &&
      as.character(expression[[2L]]) %in% wanted
  ) {
    eval(expression, envir = environment)
  }
}
stopifnot(all(wanted %in% ls(environment)))

manifest <- data.frame(
  gene = c("NO-TARGET", "NO-TARGET", "NO-TARGET", "A", "B", "NOT_FILTERED"),
  cell_count = c(7949L, 7949L, 7949L, 11L, 17L, 20L),
  file = c("partA", "partB", "partC", "partA", "partB", "partC"),
  is_in_filtered_obj = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
  stringsAsFactors = FALSE
)
catalogue <- environment$source_catalogue_from_manifest(
  manifest,
  expected_perturbations = 2L
)
stopifnot(
  identical(catalogue$targets, c("A", "B")),
  identical(catalogue$control_labels, "NO-TARGET"),
  identical(catalogue$control_rows, 3L)
)

# Match the observed Seurat v5 object[[]] schema: RNA count summaries, GDO
# summaries, guide/gene/sample character fields, and one NO-TARGET label.
metadata <- data.frame(
  orig.ident = rep("SeuratProject", 5L),
  nCount_RNA = c(1000, 1100, 1200, 1300, 1400),
  nFeature_RNA = c(500L, 510L, 520L, 530L, 540L),
  nCount_GDO = c(20, 21, 22, 23, 24),
  guide = c("NO-TARGETg0", "NO-TARGETg1", "Ag0", "Ag1", "Bg0"),
  gene = c("NO-TARGET", "NO-TARGET", "A", "A", "B"),
  sample = c("FXB1A01_L1", "FXB1B01_L1", "FXB1A01_L1", "FXB1B01_L1", "FXB1C01_L1"),
  nFeature_GDO = c(2L, 2L, 3L, 3L, 4L),
  row.names = paste0("cell", seq_len(5L)),
  stringsAsFactors = FALSE
)
assignment <- environment$choose_target_column(
  metadata,
  catalogue$targets,
  minimum_target_coverage = 2L
)
stopifnot(
  identical(assignment$column, "gene"),
  identical(sum(assignment$control), 2L),
  identical(unique(assignment$mapped[!is.na(assignment$mapped)]), c("A", "B"))
)

duplicated_perturbation <- rbind(
  manifest,
  data.frame(
    gene = "A",
    cell_count = 9L,
    file = "partC",
    is_in_filtered_obj = TRUE
  )
)
duplicate_error <- tryCatch(
  {
    environment$source_catalogue_from_manifest(
      duplicated_perturbation,
      expected_perturbations = 3L
    )
    ""
  },
  error = function(error) conditionMessage(error)
)
stopifnot(grepl("duplicated perturbations: A", duplicate_error, fixed = TRUE))

cat("3 R schema regression checks passed\n")
