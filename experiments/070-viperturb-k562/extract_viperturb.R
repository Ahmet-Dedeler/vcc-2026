#!/usr/bin/env Rscript

# Extract a compact VCC-panel donor table from the official VIPerturb Seurat
# RDS, while estimating normalization statistics over every perturbation in the
# filtered source object. No whole-cell matrix export is produced.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
  library(jsonlite)
})

EXPECTED_SOURCE_MD5 <- "ac52c41233fc7a5a5c3810e0a075007f"
EXPECTED_MANIFEST_MD5 <- "c095c8f2d16136753c22ef5d9b30a22b"
EXPECTED_FILTERED_COVERAGE <- 124L
EXPECTED_SOURCE_PERTURBATIONS <- 6724L
TARGET_SUM <- 50000
LOW_CELL_THRESHOLD <- 30L
CONTROL_LABEL <- "__CONTROL__"

parse_args <- function(values) {
  result <- list()
  index <- 1L
  while (index <= length(values)) {
    key <- values[[index]]
    if (!startsWith(key, "--")) {
      stop(sprintf("unexpected positional argument %s", key), call. = FALSE)
    }
    if (index == length(values)) {
      stop(sprintf("missing value for %s", key), call. = FALSE)
    }
    result[[substring(key, 3L)]] <- values[[index + 1L]]
    index <- index + 2L
  }
  result
}

required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(value)) {
    stop(sprintf("--%s is required", name), call. = FALSE)
  }
  value
}

read_unique_column <- function(path, column) {
  table <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!column %in% colnames(table)) {
    stop(sprintf("%s lacks column %s", path, column), call. = FALSE)
  }
  values <- trimws(as.character(table[[column]]))
  values <- values[nzchar(values)]
  if (length(values) == 0L || anyDuplicated(values)) {
    stop(sprintf("%s column %s is empty or duplicated", path, column), call. = FALSE)
  }
  values
}

is_control <- function(values) {
  text <- tolower(trimws(as.character(values)))
  grepl(
    "(^|[^a-z0-9])(non?.?target|non?.?targeting|ntc|negative|control|safe.?harbou?r|nt)([^a-z0-9]|$)",
    text,
    perl = TRUE
  )
}

source_catalogue_from_manifest <- function(
  manifest,
  expected_perturbations = EXPECTED_SOURCE_PERTURBATIONS
) {
  filtered <- manifest[manifest$is_in_filtered_obj, , drop = FALSE]
  labels <- trimws(as.character(filtered$gene))
  if (length(labels) == 0L || anyNA(labels) || any(!nzchar(labels))) {
    stop("filtered source catalogue has empty perturbation labels", call. = FALSE)
  }
  control_rows <- is_control(labels)
  if (!any(control_rows)) {
    stop("filtered source catalogue has no control rows", call. = FALSE)
  }
  control_labels <- unique(labels[control_rows])
  if (length(control_labels) != 1L) {
    stop(
      sprintf(
        "filtered source catalogue has ambiguous control labels: %s",
        paste(control_labels, collapse = ",")
      ),
      call. = FALSE
    )
  }
  targets <- labels[!control_rows]
  duplicate_targets <- unique(targets[duplicated(targets)])
  if (length(targets) == 0L || length(duplicate_targets) > 0L) {
    stop(
      sprintf(
        "filtered source catalogue has empty or duplicated perturbations: %s",
        paste(utils::head(duplicate_targets, 10L), collapse = ",")
      ),
      call. = FALSE
    )
  }
  if (length(targets) != expected_perturbations) {
    stop(
      sprintf(
        "filtered source perturbation count %d != %d",
        length(targets),
        expected_perturbations
      ),
      call. = FALSE
    )
  }
  list(
    targets = targets,
    control_labels = control_labels,
    control_rows = as.integer(sum(control_rows))
  )
}

map_target_values <- function(values, target_set) {
  raw <- trimws(as.character(values))
  raw[is.na(raw)] <- ""
  unique_values <- unique(raw)
  mapping <- rep(NA_character_, length(unique_values))
  names(mapping) <- unique_values
  exact <- unique_values %in% target_set
  mapping[exact] <- unique_values[exact]
  unresolved <- unique_values[
    !is.na(unique_values) &
      nzchar(unique_values) &
      !exact &
      !is_control(unique_values)
  ]
  for (value in unresolved) {
    tokens <- unique(unlist(strsplit(value, "[,;|:_[:space:]]+", perl = TRUE)))
    hits <- intersect(tokens, target_set)
    if (length(hits) == 1L) {
      mapping[[value]] <- hits[[1L]]
    }
  }
  unname(mapping[raw])
}

choose_target_column <- function(
  metadata,
  targets,
  explicit = NULL,
  minimum_target_coverage = 70L
) {
  if (!is.null(explicit)) {
    if (!explicit %in% colnames(metadata)) {
      stop(sprintf("requested target column %s is absent", explicit), call. = FALSE)
    }
    candidates <- explicit
  } else {
    candidates <- colnames(metadata)[vapply(
      metadata,
      function(column) is.character(column) || is.factor(column),
      logical(1)
    )]
  }
  preferred <- c(
    "Gene_assignment", "gene_target", "guide_target", "target_gene",
    "perturbation", "gene", "Guide_assignment", "guide_ID"
  )
  audits <- vector("list", length(candidates))
  mapped_values <- vector("list", length(candidates))
  for (index in seq_along(candidates)) {
    column <- candidates[[index]]
    raw <- metadata[[column]]
    mapped <- map_target_values(raw, targets)
    mapped_values[[index]] <- mapped
    audits[[index]] <- data.frame(
      column = column,
      target_coverage = length(unique(mapped[!is.na(mapped)])),
      target_cells = sum(!is.na(mapped)),
      control_cells = sum(is_control(raw)),
      preference = match(column, preferred, nomatch = length(preferred) + 1L),
      stringsAsFactors = FALSE
    )
  }
  audit <- do.call(rbind, audits)
  if (is.null(audit) || nrow(audit) == 0L) {
    stop("no character/factor metadata columns can encode perturbations", call. = FALSE)
  }
  order_index <- order(
    -audit$target_coverage,
    -audit$control_cells,
    audit$preference,
    audit$column
  )
  winner <- order_index[[1L]]
  if (
    audit$target_coverage[[winner]] < minimum_target_coverage ||
      audit$control_cells[[winner]] < 1L
  ) {
    stop(
      sprintf(
        paste0(
          "no metadata column identifies >=%d filtered source targets and ",
          "controls; best=%s coverage=%d controls=%d"
        ),
        minimum_target_coverage,
        audit$column[[winner]],
        audit$target_coverage[[winner]],
        audit$control_cells[[winner]]
      ),
      call. = FALSE
    )
  }
  list(
    column = audit$column[[winner]],
    mapped = mapped_values[[winner]],
    control = is_control(metadata[[audit$column[[winner]]]]),
    audit = audit[order_index, , drop = FALSE]
  )
}

get_expression_plan <- function(object, assay_name) {
  assay <- object[[assay_name]]
  layers <- tryCatch(SeuratObject::Layers(assay), error = function(error) character())
  count_layers <- layers[grepl("^counts($|\\.)", layers)]
  data_layers <- layers[grepl("^data($|\\.)", layers)]
  if (length(count_layers) > 0L) {
    return(list(
      mode = "layers",
      layers = count_layers,
      source = paste0("layer:", count_layers),
      kind = "counts"
    ))
  }
  if (length(data_layers) > 0L) {
    return(list(
      mode = "layers",
      layers = data_layers,
      source = paste0("layer:", data_layers),
      kind = "data"
    ))
  }
  matrix <- tryCatch(
    SeuratObject::GetAssayData(object, assay = assay_name, layer = "counts"),
    error = function(error) NULL
  )
  if (!is.null(matrix) && nrow(matrix) > 0L && ncol(matrix) > 0L) {
    return(list(
      mode = "matrix",
      matrix = matrix,
      source = "GetAssayData:counts",
      kind = "counts"
    ))
  }
  matrix <- tryCatch(
    SeuratObject::GetAssayData(object, assay = assay_name, layer = "data"),
    error = function(error) NULL
  )
  if (!is.null(matrix) && nrow(matrix) > 0L && ncol(matrix) > 0L) {
    return(list(
      mode = "matrix",
      matrix = matrix,
      source = "GetAssayData:data",
      kind = "data"
    ))
  }
  for (slot_name in c("counts", "data")) {
    matrix <- tryCatch(
      SeuratObject::GetAssayData(object, assay = assay_name, slot = slot_name),
      error = function(error) NULL
    )
    if (!is.null(matrix) && nrow(matrix) > 0L && ncol(matrix) > 0L) {
      return(list(
        mode = "matrix",
        matrix = matrix,
        source = paste0("slot:", slot_name),
        kind = slot_name
      ))
    }
  }
  stop("could not obtain a non-empty RNA counts or data matrix", call. = FALSE)
}

expression_piece <- function(object, assay_name, plan, index) {
  if (plan$mode == "matrix") {
    if (index != 1L) {
      stop("single-matrix expression plan requested more than once", call. = FALSE)
    }
    return(plan$matrix)
  }
  SeuratObject::LayerData(
    object,
    assay = assay_name,
    layer = plan$layers[[index]]
  )
}

expression_piece_count <- function(plan) {
  if (plan$mode == "matrix") {
    1L
  } else {
    length(plan$layers)
  }
}

expression_piece_name <- function(plan, index) {
  if (plan$mode == "matrix") {
    plan$source
  } else {
    plan$source[[index]]
  }
}

resolve_features <- function(assay, matrix, official) {
  candidates <- list(rownames = rownames(matrix))
  feature_metadata <- tryCatch(assay[[]], error = function(error) NULL)
  if (!is.null(feature_metadata)) {
    for (column in intersect(
      c("gene_name", "gene_symbol", "symbol", "feature_name"),
      colnames(feature_metadata)
    )) {
      values <- as.character(feature_metadata[[column]])
      if (length(values) == nrow(matrix)) {
        candidates[[paste0("feature_metadata:", column)]] <- values
      }
    }
  }
  scores <- vapply(
    candidates,
    function(values) length(intersect(official, values)),
    integer(1)
  )
  source <- names(scores)[which.max(scores)]
  labels <- candidates[[source]]
  official_rows <- match(official, labels)
  keep <- !is.na(official_rows)
  duplicate_labels <- unique(labels[duplicated(labels)])
  list(
    source = source,
    genes = official[keep],
    rows = as.integer(official_rows[keep]),
    positions = which(keep),
    overlap = sum(keep),
    duplicate_labels = duplicate_labels
  )
}

accumulate_expression <- function(
  object,
  assay_name,
  plan,
  selected_names,
  selected_labels,
  official,
  group_levels,
  min_genes,
  group_chunk_size
) {
  group_sums <- matrix(0, nrow = length(official), ncol = length(group_levels))
  group_counts <- integer(length(group_levels))
  feature_seen <- logical(length(official))
  cell_seen <- logical(length(selected_names))
  feature_sources <- character()
  duplicate_features <- character()
  for (piece_index in seq_len(expression_piece_count(plan))) {
    matrix_piece <- expression_piece(object, assay_name, plan, piece_index)
    if (is.null(colnames(matrix_piece)) || is.null(rownames(matrix_piece))) {
      stop("expression matrix piece lacks row or column names", call. = FALSE)
    }
    selected_positions <- which(selected_names %in% colnames(matrix_piece))
    if (length(selected_positions) == 0L) {
      rm(matrix_piece)
      invisible(gc())
      next
    }
    if (any(cell_seen[selected_positions])) {
      duplicated <- selected_names[selected_positions[cell_seen[selected_positions]]]
      stop(
        sprintf(
          "selected cells occur in multiple expression layers: %s",
          paste(utils::head(duplicated, 10L), collapse = ",")
        ),
        call. = FALSE
      )
    }
    feature_map <- resolve_features(object[[assay_name]], matrix_piece, official)
    if (feature_map$overlap < min_genes) {
      stop(
        sprintf(
          "%s maps only %d official VCC genes",
          expression_piece_name(plan, piece_index),
          feature_map$overlap
        ),
        call. = FALSE
      )
    }
    group_codes <- match(selected_labels[selected_positions], group_levels)
    if (anyNA(group_codes)) {
      stop("selected cells contain an unmapped group", call. = FALSE)
    }
    piece_names <- selected_names[selected_positions]
    cell_columns <- match(piece_names, colnames(matrix_piece))
    for (
      group_lo in seq.int(1L, length(group_levels), by = group_chunk_size)
    ) {
      group_hi <- min(group_lo + group_chunk_size - 1L, length(group_levels))
      in_chunk <- group_codes >= group_lo & group_codes <= group_hi
      if (!any(in_chunk)) {
        next
      }
      chunk_codes <- group_codes[in_chunk] - group_lo + 1L
      indicator <- Matrix::sparseMatrix(
        i = seq_along(chunk_codes),
        j = chunk_codes,
        x = 1,
        dims = c(length(chunk_codes), group_hi - group_lo + 1L)
      )
      selected_matrix <- matrix_piece[
        feature_map$rows,
        cell_columns[in_chunk],
        drop = FALSE
      ]
      if (!inherits(selected_matrix, "Matrix")) {
        selected_matrix <- methods::as(selected_matrix, "dgCMatrix")
      }
      piece_sums <- as.matrix(selected_matrix %*% indicator)
      group_sums[feature_map$positions, group_lo:group_hi] <-
        group_sums[
          feature_map$positions,
          group_lo:group_hi,
          drop = FALSE
        ] + piece_sums
      group_counts[group_lo:group_hi] <-
        group_counts[group_lo:group_hi] +
        tabulate(chunk_codes, nbins = group_hi - group_lo + 1L)
      rm(selected_matrix, piece_sums, indicator)
    }
    feature_seen[feature_map$positions] <- TRUE
    cell_seen[selected_positions] <- TRUE
    feature_sources <- c(
      feature_sources,
      sprintf(
        "%s=%s",
        expression_piece_name(plan, piece_index),
        feature_map$source
      )
    )
    duplicate_features <- union(
      duplicate_features,
      feature_map$duplicate_labels
    )
    rm(matrix_piece)
    invisible(gc())
  }
  if (!all(cell_seen)) {
    missing <- selected_names[!cell_seen]
    stop(
      sprintf(
        "%d selected cells are absent from every expression layer: %s",
        length(missing),
        paste(utils::head(missing, 10L), collapse = ",")
      ),
      call. = FALSE
    )
  }
  if (!any(feature_seen)) {
    stop("no official VCC features were accumulated", call. = FALSE)
  }
  list(
    genes = official[feature_seen],
    group_sums = group_sums[feature_seen, , drop = FALSE],
    group_counts = group_counts,
    feature_sources = feature_sources,
    duplicate_features = duplicate_features
  )
}

write_effect_table <- function(path, targets, genes, effects) {
  connection <- gzfile(path, open = "wt", compression = 6)
  on.exit(close(connection), add = TRUE)
  writeLines(paste(c("target_gene", genes), collapse = "\t"), connection)
  for (index in seq_along(targets)) {
    values <- formatC(
      effects[index, ],
      digits = 17L,
      format = "g",
      flag = "#"
    )
    writeLines(paste(c(targets[[index]], values), collapse = "\t"), connection)
  }
}

write_control_table <- function(path, genes, control) {
  connection <- gzfile(path, open = "wt", compression = 6)
  on.exit(close(connection), add = TRUE)
  writeLines("gene\tcontrol_profile", connection)
  for (index in seq_along(genes)) {
    value <- formatC(control[[index]], digits = 9L, format = "g", flag = "#")
    writeLines(paste(genes[[index]], value, sep = "\t"), connection)
  }
}

write_gene_value_table <- function(path, genes, column, values) {
  connection <- gzfile(path, open = "wt", compression = 6)
  on.exit(close(connection), add = TRUE)
  writeLines(paste("gene", column, sep = "\t"), connection)
  for (index in seq_along(genes)) {
    value <- formatC(values[[index]], digits = 17L, format = "g", flag = "#")
    writeLines(paste(genes[[index]], value, sep = "\t"), connection)
  }
}

profile_effect <- function(group_sums, column, control_profile, kind, group_count) {
  values <- group_sums[, column]
  if (kind == "counts") {
    library <- sum(values)
    if (!is.finite(library) || library <= 0) {
      stop(sprintf("source group %d has an empty count library", column), call. = FALSE)
    }
    return(log1p(TARGET_SUM * values / library) - control_profile)
  }
  values / group_count - control_profile
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  input <- normalizePath(required_arg(args, "input"), mustWork = TRUE)
  manifest_path <- normalizePath(required_arg(args, "manifest"), mustWork = TRUE)
  genes_path <- normalizePath(required_arg(args, "genes-csv"), mustWork = TRUE)
  panel_path <- normalizePath(required_arg(args, "panel-csv"), mustWork = TRUE)
  output <- normalizePath(required_arg(args, "output"), mustWork = FALSE)
  target_column <- args[["target-column"]]
  min_genes <- if (is.null(args[["min-genes"]])) {
    7000L
  } else {
    as.integer(args[["min-genes"]])
  }
  low_threshold <- if (is.null(args[["low-cell-threshold"]])) {
    LOW_CELL_THRESHOLD
  } else {
    as.integer(args[["low-cell-threshold"]])
  }
  group_chunk_size <- if (is.null(args[["group-chunk-size"]])) {
    256L
  } else {
    as.integer(args[["group-chunk-size"]])
  }
  if (
    is.na(min_genes) || min_genes < 1L ||
      is.na(low_threshold) || low_threshold < 1L ||
      is.na(group_chunk_size) || group_chunk_size < 1L
  ) {
    stop("min-genes, low-cell-threshold, and group-chunk-size must be positive", call. = FALSE)
  }
  verified_source_md5 <- if (is.null(args[["verified-source-md5"]])) {
    EXPECTED_SOURCE_MD5
  } else {
    args[["verified-source-md5"]]
  }
  verified_manifest_md5 <- if (is.null(args[["verified-manifest-md5"]])) {
    EXPECTED_MANIFEST_MD5
  } else {
    args[["verified-manifest-md5"]]
  }
  if (verified_source_md5 != EXPECTED_SOURCE_MD5) {
    stop("source MD5 is not the frozen official checksum", call. = FALSE)
  }
  if (verified_manifest_md5 != EXPECTED_MANIFEST_MD5) {
    stop("manifest MD5 is not the frozen official checksum", call. = FALSE)
  }
  if (dir.exists(output) || file.exists(output)) {
    stop(sprintf("refusing to overwrite %s", output), call. = FALSE)
  }

  official <- read_unique_column(genes_path, "gene_name")
  panel <- read_unique_column(panel_path, "target_gene")
  manifest <- utils::read.table(
    manifest_path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  expected_manifest_columns <- c("gene", "cell_count", "file", "is_in_filtered_obj")
  if (!identical(colnames(manifest), expected_manifest_columns)) {
    stop("unexpected genome_wide_manifest.txt columns", call. = FALSE)
  }
  manifest$is_in_filtered_obj <- toupper(as.character(manifest$is_in_filtered_obj)) == "TRUE"
  source_catalogue <- source_catalogue_from_manifest(manifest)
  source_targets <- source_catalogue$targets
  filtered_targets <- panel[
    panel %in% manifest$gene[
      manifest$is_in_filtered_obj
    ]
  ]
  filtered_manifest_rows <- manifest[
    manifest$gene %in% filtered_targets,
    ,
    drop = FALSE
  ]
  if (anyDuplicated(filtered_manifest_rows$gene)) {
    stop("manifest contains duplicate filtered VCC target rows", call. = FALSE)
  }
  if (length(filtered_targets) != EXPECTED_FILTERED_COVERAGE) {
    stop(
      sprintf(
        "manifest filtered VCC coverage %d != %d",
        length(filtered_targets),
        EXPECTED_FILTERED_COVERAGE
      ),
      call. = FALSE
    )
  }
  filtered_manifest <- manifest[
    match(filtered_targets, manifest$gene),
    ,
    drop = FALSE
  ]

  message("READ_RDS ", input)
  object <- readRDS(input)
  if (!inherits(object, "Seurat")) {
    stop("genome_wide_filtered.rds is not a Seurat object", call. = FALSE)
  }
  assays <- Seurat::Assays(object)
  assay_name <- if ("RNA" %in% assays) {
    "RNA"
  } else if (length(assays) == 1L) {
    assays[[1L]]
  } else {
    stop(sprintf("RNA assay absent; available assays: %s", paste(assays, collapse = ",")), call. = FALSE)
  }
  metadata <- object[[]]
  if (nrow(metadata) == 0L || is.null(rownames(metadata))) {
    stop("Seurat metadata is empty or lacks cell names", call. = FALSE)
  }
  assignment <- choose_target_column(
    metadata,
    source_targets,
    target_column,
    minimum_target_coverage = length(source_targets)
  )
  mapped <- assignment$mapped
  controls <- assignment$control
  raw_assignment <- trimws(as.character(metadata[[assignment$column]]))
  unmapped_non_controls <- unique(
    raw_assignment[
      !controls &
        !is.na(raw_assignment) &
        nzchar(raw_assignment) &
        is.na(mapped)
    ]
  )
  if (length(unmapped_non_controls) > 0L) {
    stop(
      sprintf(
        "selected metadata column %s has %d unexpected non-control labels: %s",
        assignment$column,
        length(unmapped_non_controls),
        paste(utils::head(unmapped_non_controls, 10L), collapse = ",")
      ),
      call. = FALSE
    )
  }
  selected <- controls | !is.na(mapped)
  selected_names <- rownames(metadata)[selected]
  selected_labels <- ifelse(controls[selected], CONTROL_LABEL, mapped[selected])
  observed_targets <- unique(selected_labels[selected_labels != CONTROL_LABEL])
  missing_targets <- setdiff(source_targets, observed_targets)
  if (length(missing_targets) > 0L) {
    stop(
      sprintf(
        "selected metadata column %s misses %d filtered source targets: %s",
        assignment$column,
        length(missing_targets),
        paste(missing_targets, collapse = ",")
      ),
      call. = FALSE
    )
  }

  group_levels <- c(CONTROL_LABEL, source_targets)
  expression <- get_expression_plan(object, assay_name)
  accumulated <- accumulate_expression(
    object,
    assay_name,
    expression,
    selected_names,
    selected_labels,
    official,
    group_levels,
    min_genes,
    group_chunk_size
  )
  group_sums <- accumulated$group_sums
  group_counts <- accumulated$group_counts
  output_genes <- accumulated$genes
  if (group_counts[[1L]] < 1L || any(group_counts[-1L] < 1L)) {
    stop("control or filtered source target has zero selected cells", call. = FALSE)
  }
  if (any(!is.finite(group_sums))) {
    stop("non-finite pseudobulk sums", call. = FALSE)
  }

  if (expression$kind == "counts") {
    control_library <- sum(group_sums[, 1L])
    if (!is.finite(control_library) || control_library <= 0) {
      stop("control pseudobulk contains an empty library", call. = FALSE)
    }
    control_profile <- log1p(TARGET_SUM * group_sums[, 1L] / control_library)
    effect_definition <- "log1p(50000 * group_count_sum / library) minus control"
  } else {
    control_profile <- group_sums[, 1L] / group_counts[[1L]]
    effect_definition <- "mean Seurat data-layer expression minus control mean"
  }

  source_sum <- numeric(length(output_genes))
  source_sum_squares <- numeric(length(output_genes))
  panel_effects <- matrix(
    0,
    nrow = length(filtered_targets),
    ncol = length(output_genes)
  )
  panel_index <- match(source_targets, filtered_targets)
  for (source_index in seq_along(source_targets)) {
    effect <- profile_effect(
      group_sums,
      source_index + 1L,
      control_profile,
      expression$kind,
      group_counts[[source_index + 1L]]
    )
    if (any(!is.finite(effect))) {
      stop(sprintf("non-finite effect for %s", source_targets[[source_index]]), call. = FALSE)
    }
    source_sum <- source_sum + effect
    source_sum_squares <- source_sum_squares + effect * effect
    panel_row <- panel_index[[source_index]]
    if (!is.na(panel_row)) {
      panel_effects[panel_row, ] <- effect
    }
  }
  source_center <- source_sum / length(source_targets)
  source_variance <- pmax(
    source_sum_squares / length(source_targets) - source_center * source_center,
    0
  )
  source_norms <- numeric(length(source_targets))
  for (source_index in seq_along(source_targets)) {
    effect <- profile_effect(
      group_sums,
      source_index + 1L,
      control_profile,
      expression$kind,
      group_counts[[source_index + 1L]]
    )
    source_norms[[source_index]] <- sqrt(sum((effect - source_center)^2))
  }
  source_median_norm <- median(source_norms)
  if (!is.finite(source_median_norm) || source_median_norm <= 0) {
    stop("source-wide median centered norm is not positive", call. = FALSE)
  }
  panel_center <- colMeans(panel_effects)
  panel_centered <- sweep(
    panel_effects,
    2L,
    panel_center,
    FUN = "-",
    check.margin = TRUE
  )
  panel_norms <- sqrt(rowSums(panel_centered * panel_centered))
  rm(panel_centered)
  panel_median_norm <- median(panel_norms)
  if (!is.finite(panel_median_norm) || panel_median_norm <= 0) {
    stop("panel-centered median norm is not positive", call. = FALSE)
  }

  source_actual_counts <- group_counts[-1L]
  actual_counts <- source_actual_counts[match(filtered_targets, source_targets)]
  target_audit <- data.frame(
    target_gene = filtered_targets,
    manifest_cell_count = as.integer(filtered_manifest$cell_count),
    observed_cell_count = as.integer(actual_counts),
    count_difference = as.integer(actual_counts - filtered_manifest$cell_count),
    low_cell_manifest = filtered_manifest$cell_count < low_threshold,
    low_cell_observed = actual_counts < low_threshold,
    stringsAsFactors = FALSE
  )

  parent <- dirname(output)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".exp070-compact-", tmpdir = parent)
  dir.create(temporary)
  completed <- FALSE
  on.exit({
    if (!completed && dir.exists(temporary)) {
      unlink(temporary, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)

  writeLines(output_genes, file.path(temporary, "genes.txt"))
  writeLines(filtered_targets, file.path(temporary, "perts.txt"))
  writeLines(source_targets, file.path(temporary, "source_perts.txt"))
  write_effect_table(
    file.path(temporary, "effects.tsv.gz"),
    filtered_targets,
    output_genes,
    panel_effects
  )
  write_control_table(
    file.path(temporary, "control.tsv.gz"),
    output_genes,
    control_profile
  )
  write_gene_value_table(
    file.path(temporary, "source_center.tsv.gz"),
    output_genes,
    "source_center",
    source_center
  )
  write_gene_value_table(
    file.path(temporary, "source_variance.tsv.gz"),
    output_genes,
    "source_variance",
    source_variance
  )
  write_gene_value_table(
    file.path(temporary, "panel_center.tsv.gz"),
    output_genes,
    "panel_center",
    panel_center
  )
  utils::write.csv(
    target_audit,
    file.path(temporary, "target_audit.csv"),
    row.names = FALSE,
    quote = TRUE
  )
  utils::write.csv(
    assignment$audit,
    file.path(temporary, "metadata_column_audit.csv"),
    row.names = FALSE,
    quote = TRUE
  )
  structure_lines <- c(
    sprintf("class=%s", paste(class(object), collapse = ",")),
    sprintf("assays=%s", paste(assays, collapse = ",")),
    sprintf("selected_assay=%s", assay_name),
    sprintf("layers=%s", paste(tryCatch(SeuratObject::Layers(object[[assay_name]]), error = function(error) character()), collapse = ",")),
    sprintf("expression_source=%s", paste(expression$source, collapse = ",")),
    sprintf("expression_kind=%s", expression$kind),
    sprintf("target_column=%s", assignment$column),
    sprintf("manifest_control_labels=%s", paste(source_catalogue$control_labels, collapse = ",")),
    sprintf("manifest_control_rows=%d", source_catalogue$control_rows),
    sprintf("metadata_columns=%s", paste(colnames(metadata), collapse = ","))
  )
  writeLines(structure_lines, file.path(temporary, "object_audit.txt"))
  writeLines(capture.output(sessionInfo()), file.path(temporary, "r_session_info.txt"))

  artifact_names <- c(
    "genes.txt", "perts.txt", "source_perts.txt", "effects.tsv.gz",
    "control.tsv.gz", "source_center.tsv.gz", "source_variance.tsv.gz",
    "panel_center.tsv.gz", "target_audit.csv", "metadata_column_audit.csv",
    "object_audit.txt", "r_session_info.txt"
  )
  artifact_md5 <- as.list(unname(tools::md5sum(file.path(temporary, artifact_names))))
  names(artifact_md5) <- artifact_names
  metadata_output <- list(
    experiment = "070-viperturb-k562",
    dataset = "Genome-wide single-cell perturbation screens with VIPerturb-seq",
    zenodo_record = "18460279",
    source_file = "genome_wide_filtered.rds",
    source_bytes = 3611371127,
    source_md5 = verified_source_md5,
    manifest_file = "genome_wide_manifest.txt",
    manifest_md5 = verified_manifest_md5,
    manifest_control_labels = as.list(source_catalogue$control_labels),
    manifest_control_rows = unname(source_catalogue$control_rows),
    manifest_filtered_vcc_coverage = length(filtered_targets),
    filtered_targets = filtered_targets,
    normalization_scope = "all_filtered_source_perturbations",
    n_source_perturbations = length(source_targets),
    source_perturbation_order_md5 = unname(
      tools::md5sum(file.path(temporary, "source_perts.txt"))
    ),
    source_median_centered_norm = unname(source_median_norm),
    panel_median_centered_norm = unname(panel_median_norm),
    selected_assay = assay_name,
    expression_source = paste(expression$source, collapse = ","),
    expression_kind = expression$kind,
    effect_definition = effect_definition,
    target_sum = if (expression$kind == "counts") TARGET_SUM else NULL,
    target_column = assignment$column,
    control_label = CONTROL_LABEL,
    n_control_cells = unname(group_counts[[1L]]),
    n_perturbations = nrow(panel_effects),
    n_genes = ncol(panel_effects),
    feature_label_source = as.list(accumulated$feature_sources),
    duplicate_feature_labels = as.list(accumulated$duplicate_features),
    low_cell_threshold = low_threshold,
    manifest_low_cell_targets = as.list(
      filtered_targets[filtered_manifest$cell_count < low_threshold]
    ),
    observed_low_cell_targets = as.list(
      filtered_targets[actual_counts < low_threshold]
    ),
    manifest_observed_count_mismatches = sum(target_audit$count_difference != 0L),
    matrix_layout = list(
      effects = "gzip TSV; raw control-relative panel rows by official VCC gene",
      control = "gzip TSV; one control profile value per official VCC gene",
      source_center = "gzip TSV; source-wide mean raw effect by gene",
      source_variance = "gzip TSV; source-wide population variance by gene",
      panel_center = "gzip TSV; panel-only mean raw effect by gene"
    ),
    artifacts_md5 = artifact_md5
  )
  jsonlite::write_json(
    metadata_output,
    file.path(temporary, "meta.json"),
    auto_unbox = TRUE,
    pretty = TRUE,
    digits = NA,
    null = "null"
  )
  if (!file.rename(temporary, output)) {
    stop(sprintf("failed to atomically publish compact output %s", output), call. = FALSE)
  }
  completed <- TRUE
  message(
    sprintf(
      paste0(
        "WROTE_COMPACT panel_targets=%d source_targets=%d genes=%d ",
        "controls=%d low_manifest=%d low_observed=%d"
      ),
      nrow(panel_effects),
      length(source_targets),
      ncol(panel_effects),
      group_counts[[1L]],
      sum(filtered_manifest$cell_count < low_threshold),
      sum(actual_counts < low_threshold)
    )
  )
}

main()
