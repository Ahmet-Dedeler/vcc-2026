#!/usr/bin/env Rscript

# Synthetic-only oracle for the frozen Python EmptyDrops primitives.
#
# This script intentionally calls the pinned R package implementations. It is
# never used on GSE249595 itself; the production pass remains the streaming
# Python implementation after this oracle agrees on controlled matrices.

suppressPackageStartupMessages({
    library(BiocParallel)
    library(DropletUtils)
    library(jsonlite)
    library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
    stop("usage: r_oracle.R INPUT.mtx OUTPUT.json")
}

expected <- c(DropletUtils = "1.31.2", edgeR = "4.99.0")
observed <- c(
    DropletUtils = as.character(packageVersion("DropletUtils")),
    edgeR = as.character(packageVersion("edgeR"))
)
if (!identical(observed, expected)) {
    stop(
        sprintf(
            "pinned package mismatch: observed %s; expected %s",
            paste(names(observed), observed, collapse = ", "),
            paste(names(expected), expected, collapse = ", ")
        )
    )
}

matrix <- as(readMM(args[[1L]]), "dgCMatrix")
lower <- 100L
barcode <- barcodeRanks(
    matrix,
    lower = lower,
    exclude.from = 50L,
    window = 1,
    gradient.threshold = -1,
    BPPARAM = SerialParam()
)
ambient <- ambientProfileEmpty(
    matrix,
    lower = lower,
    round = FALSE,
    good.turing = TRUE,
    BPPARAM = SerialParam()
)
totals <- unname(Matrix::colSums(matrix))
tested <- totals > lower
observed.data <- DropletUtils:::.compute_multinom_prob_data(
    matrix[, tested, drop = FALSE],
    prop = ambient,
    alpha = Inf,
    BPPARAM = SerialParam()
)

bh.input <- c(0.9, 0.01, 0.01, 0.2, 1, 0.049, 0.5, 0)
payload <- list(
    versions = as.list(observed),
    totals = totals,
    knee = S4Vectors::metadata(barcode)$knee,
    inflection = S4Vectors::metadata(barcode)$inflection,
    ambient = unname(ambient),
    tested = which(tested) - 1L,
    observed_data_log_probability = unname(observed.data),
    bh_input = bh.input,
    bh_adjusted = p.adjust(bh.input, method = "BH")
)
write_json(payload, args[[2L]], auto_unbox = TRUE, digits = NA, pretty = TRUE)
