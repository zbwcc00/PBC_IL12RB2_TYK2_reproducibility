suppressPackageStartupMessages({
  library(BiocParallel)
  library(Matrix)
  library(scuttle)
  library(SingleCellExperiment)
  library(SingleR)
})

arguments <- commandArgs(trailingOnly = TRUE)
set.seed(38003)
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_argument), winslash = "/")
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/")
work_root <- file.path(repo_root, "work", "single_cell", "hra008003")
requested_samples <- unique(arguments)
matrix_root <- file.path(work_root, "reconstruction")
adaptive_root <- file.path(work_root, "qc")
soupx_root <- file.path(work_root, "soupx")
output_root <- file.path(work_root, "annotation")
reference_root <- file.path(repo_root, "data", "raw", "single_cell", "references", "reference_pseudobulk")
dir.create(file.path(output_root, "sample_metadata"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "sample_predictions"), recursive = TRUE, showWarnings = FALSE)

sample_name <- function(donor_id) {
  matches <- regexec("^(control|pbc)([0-9]+)$", donor_id, ignore.case = TRUE)
  pieces <- regmatches(donor_id, matches)[[1]]
  if (length(pieces) != 3L) {
    stop("Unexpected donor_id: ", donor_id)
  }
  paste0(tolower(pieces[[2]]), "_liver", pieces[[3]])
}

order_table <- read.delim(
  file.path(matrix_root, "01_liver_BAM_processing_order.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
order_table <- order_table[order(order_table$processing_order), ]
order_table$sample_id <- vapply(order_table$donor_id, sample_name, character(1))
if (length(requested_samples) > 0L) {
  unknown_samples <- setdiff(requested_samples, order_table$sample_id)
  if (length(unknown_samples) > 0L) {
    stop("Unknown requested samples: ", paste(unknown_samples, collapse = ", "))
  }
  order_table <- order_table[order_table$sample_id %in% requested_samples, ]
}

reference_counts <- as(
  readMM(gzfile(file.path(reference_root, "CELLxGENE_reference_pseudobulk_counts.mtx.gz"))),
  "dgCMatrix"
)
reference_features <- read.delim(
  gzfile(file.path(reference_root, "CELLxGENE_reference_features.tsv.gz")),
  stringsAsFactors = FALSE
)
reference_profiles <- read.delim(
  file.path(reference_root, "CELLxGENE_reference_profiles.tsv"),
  stringsAsFactors = FALSE
)
rownames(reference_counts) <- reference_features$feature_id
colnames(reference_counts) <- reference_profiles$index
reference <- SingleCellExperiment(
  assays = list(counts = reference_counts),
  colData = reference_profiles
)
reference <- logNormCounts(reference)

sample_feature_tables <- list()
for (sample_id in order_table$sample_id) {
  qc_path <- file.path(
    adaptive_root,
    "objects",
    paste0(sample_id, "_adaptive_called_cells.rds")
  )
  if (!file.exists(qc_path)) {
    stop("Adaptive QC object is missing for requested sample: ", qc_path)
  }
  qc_object <- readRDS(qc_path)
  sample_feature_tables[[sample_id]] <- data.frame(
    feature_id = rownames(qc_object$counts),
    gene_symbol = qc_object$feature_data$Symbol,
    stringsAsFactors = FALSE
  )
  rm(qc_object)
  gc()
}
common_feature_ids <- Reduce(
  intersect,
  lapply(sample_feature_tables, function(table) table$feature_id)
)
query_features <- sample_feature_tables[[1]][
  match(common_feature_ids, sample_feature_tables[[1]]$feature_id),
  ,
  drop = FALSE
]

excluded_targets <- c(
  "IL12RB2", "IL12RB1", "JAK2", "STAT4", "IFNG",
  "FCRL3", "FCRL3-AS1", "AL356276.1", "CD226", "DOK6"
)
nuisance_rows <- grepl(
  "^(MT-|RPS|RPL|HB[ABDEGQZ]|MALAT1$)",
  toupper(query_features$gene_symbol)
)
shared_rows <- query_features$feature_id %in% rownames(reference)
eligible_rows <- shared_rows & !nuisance_rows & !query_features$gene_symbol %in% excluded_targets
classification_ids <- query_features$feature_id[eligible_rows]
if (length(classification_ids) < 5000L) {
  stop("Fewer than 5000 target-independent shared genes are available for annotation")
}

trained <- trainSingleR(
  ref = reference[classification_ids, ],
  labels = reference$author_cell_type,
  test.genes = classification_ids,
  genes = "de",
  de.method = "classic",
  de.n = 50,
  aggr.ref = TRUE,
  BPPARAM = SerialParam()
)

pure_tnk_reference_labels <- c("lrNK", "cNK", "CD4T", "CD8T", "Tcell")
pure_tnk_reference <- reference[, reference$author_cell_type %in% pure_tnk_reference_labels]
trained_pure_tnk <- trainSingleR(
  ref = pure_tnk_reference[classification_ids, ],
  labels = pure_tnk_reference$author_cell_type,
  test.genes = classification_ids,
  genes = "de",
  de.method = "classic",
  de.n = 50,
  aggr.ref = TRUE,
  BPPARAM = SerialParam()
)

author_to_coarse <- vapply(
  split(reference_profiles$coarse_label, reference_profiles$author_cell_type),
  function(values) {
    unique_values <- unique(values)
    if (length(unique_values) != 1L) {
      stop("Reference author label maps to multiple coarse labels")
    }
    unique_values[[1]]
  },
  character(1)
)
fine_label_map <- c(
  lrNK = "lrNK",
  cNK = "cNK",
  CD4T = "CD4_T",
  CD8T = "CD8_T",
  Tcell = "T_other",
  `CD3T-lrNK` = "TNK_mixed",
  `CD8T-cNK` = "TNK_mixed",
  NKT = "TNK_mixed"
)

metadata_rows <- list()
summary_rows <- list()
for (row_index in seq_len(nrow(order_table))) {
  record <- order_table[row_index, ]
  sample_id <- record$sample_id
  cat("Annotating", sample_id, "\n")
  corrected_path <- file.path(soupx_root, sample_id, "SoupX_corrected_counts.rds")
  adaptive_path <- file.path(
    adaptive_root,
    "objects",
    paste0(sample_id, "_adaptive_called_cells.rds")
  )
  if (!file.exists(corrected_path) || !file.exists(adaptive_path)) {
    stop("Required post-QC input is missing for ", sample_id)
  }

  corrected_counts <- readRDS(corrected_path)
  adaptive <- readRDS(adaptive_path)
  sample_features <- sample_feature_tables[[sample_id]]
  if (!identical(rownames(corrected_counts), sample_features$feature_id)) {
    stop("Corrected-count and adaptive-QC feature rows differ within ", sample_id)
  }
  if (!all(classification_ids %in% rownames(corrected_counts))) {
    stop("Target-independent common annotation genes are incomplete for ", sample_id)
  }
  final_metadata <- adaptive$cell_data[adaptive$cell_data$passes_final_qc, , drop = FALSE]
  final_metadata <- final_metadata[match(colnames(corrected_counts), final_metadata$barcode), , drop = FALSE]
  if (!identical(final_metadata$barcode, colnames(corrected_counts))) {
    stop("Final singlet metadata do not align to corrected counts for ", sample_id)
  }

  query <- SingleCellExperiment(
    assays = list(counts = corrected_counts[classification_ids, , drop = FALSE])
  )
  query <- logNormCounts(query)
  predicted <- classifySingleR(
    test = query,
    trained = trained,
    prune = FALSE,
    BPPARAM = SerialParam()
  )
  pruned <- pruneScores(predicted)
  reference_label <- as.character(predicted$labels)
  reference_coarse_label <- unname(author_to_coarse[reference_label])
  conservative_label <- unname(fine_label_map[reference_label])
  conservative_label[is.na(conservative_label)] <- reference_coarse_label[is.na(conservative_label)]
  conservative_label[pruned] <- "uncertain"

  tnk_candidate_labels <- c(
    pure_tnk_reference_labels,
    "CD3T-lrNK",
    "CD8T-cNK",
    "NKT"
  )
  tnk_candidate <- reference_label %in% tnk_candidate_labels
  tnk_refined_label <- rep(NA_character_, ncol(query))
  tnk_refined_delta <- rep(NA_real_, ncol(query))
  tnk_refined_pruned <- rep(NA, ncol(query))
  tnk_predicted <- NULL
  if (sum(tnk_candidate) > 0L) {
    tnk_predicted <- classifySingleR(
      test = query[, tnk_candidate, drop = FALSE],
      trained = trained_pure_tnk,
      prune = FALSE,
      BPPARAM = SerialParam()
    )
    tnk_pruned <- pruneScores(tnk_predicted)
    tnk_indices <- which(tnk_candidate)
    tnk_labels <- unname(fine_label_map[as.character(tnk_predicted$labels)])
    tnk_labels[tnk_pruned] <- "TNK_unresolved"
    tnk_refined_label[tnk_indices] <- tnk_labels
    tnk_refined_delta[tnk_indices] <- as.numeric(tnk_predicted$delta.next)
    tnk_refined_pruned[tnk_indices] <- as.logical(tnk_pruned)
  }

  raw_il12rb2_rows <- which(sample_features$gene_symbol == "IL12RB2")
  raw_final_counts <- adaptive$counts[, adaptive$cell_data$passes_final_qc, drop = FALSE]
  raw_final_counts <- raw_final_counts[, match(colnames(corrected_counts), colnames(raw_final_counts)), drop = FALSE]
  raw_il12rb2 <- Matrix::colSums(raw_final_counts[raw_il12rb2_rows, , drop = FALSE])
  corrected_il12rb2 <- Matrix::colSums(corrected_counts[raw_il12rb2_rows, , drop = FALSE])

  sample_metadata <- data.frame(
    cell_id = paste(sample_id, colnames(corrected_counts), sep = "::"),
    barcode = colnames(corrected_counts),
    sample_id = sample_id,
    run = record$runAcc,
    donor_id = record$donor_id,
    condition = record$condition,
    total_UMI_raw = final_metadata$total_UMI,
    detected_genes = final_metadata$detected_genes,
    percent_mt = final_metadata$percent_mt,
    scDblFinder_score = final_metadata$scDblFinder_score,
    total_UMI_SoupX = Matrix::colSums(corrected_counts),
    raw_IL12RB2_UMI = raw_il12rb2,
    SoupX_IL12RB2_UMI = corrected_il12rb2,
    reference_label = reference_label,
    reference_coarse_label = reference_coarse_label,
    annotation_delta = as.numeric(predicted$delta.next),
    annotation_pruned = as.logical(pruned),
    conservative_label = conservative_label,
    tnk_refined_label = tnk_refined_label,
    tnk_refined_delta = tnk_refined_delta,
    tnk_refined_pruned = tnk_refined_pruned,
    stringsAsFactors = FALSE
  )

  metadata_connection <- gzfile(
    file.path(output_root, "sample_metadata", paste0(sample_id, "_cell_metadata.tsv.gz")),
    "wt"
  )
  write.table(sample_metadata, metadata_connection, sep = "\t", quote = FALSE, row.names = FALSE)
  close(metadata_connection)
  saveRDS(
    list(primary = predicted, pure_tnk_sensitivity = tnk_predicted),
    file.path(output_root, "sample_predictions", paste0(sample_id, "_SingleR_predictions.rds")),
    compress = "gzip"
  )

  metadata_rows[[sample_id]] <- sample_metadata
  sample_summary <- as.data.frame(table(sample_metadata$conservative_label), stringsAsFactors = FALSE)
  names(sample_summary) <- c("conservative_label", "n_cells")
  sample_summary$sample_id <- sample_id
  sample_summary$donor_id <- record$donor_id
  sample_summary$condition <- record$condition
  summary_rows[[sample_id]] <- sample_summary[, c(
    "sample_id", "donor_id", "condition", "conservative_label", "n_cells"
  )]

  rm(corrected_counts, adaptive, raw_final_counts, query, predicted, tnk_predicted)
  gc()
}

all_metadata <- do.call(rbind, metadata_rows)
rownames(all_metadata) <- NULL
all_summary <- do.call(rbind, summary_rows)
rownames(all_summary) <- NULL

metadata_connection <- gzfile(file.path(output_root, "all_processed_cell_metadata.tsv.gz"), "wt")
write.table(all_metadata, metadata_connection, sep = "\t", quote = FALSE, row.names = FALSE)
close(metadata_connection)
write.table(
  all_summary,
  file.path(output_root, "annotation_cell_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  data.frame(
    feature_id = classification_ids,
    gene_symbol = query_features$gene_symbol[match(classification_ids, query_features$feature_id)]
  ),
  file.path(output_root, "classification_gene_universe.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  data.frame(
    parameter = c(
      "reference",
      "reference_profiles",
      "classification_genes",
      "excluded_targets",
      "nuisance_exclusion",
      "annotation_method",
      "confidence_rule",
      "mixed_TNK_rule",
      "pure_TNK_sensitivity",
      "random_seed"
    ),
    value = c(
      "CELLxGENE liver normal/PBC/PSC donor-celltype pseudobulks reused from matched project reference",
      ncol(reference),
      length(classification_ids),
      paste(excluded_targets, collapse = ","),
      "MT-, RPS, RPL, HB family, MALAT1",
      "SingleR trained once on target-independent shared genes",
      "SingleR pruneScores; pruned cells labeled uncertain",
      "CD3T-lrNK, CD8T-cNK, and NKT retained as TNK_mixed in the primary label",
      "TNK candidates reclassified against pure lrNK/cNK/CD4T/CD8T/Tcell profiles; pruned calls labeled TNK_unresolved and used only for sensitivity analyses",
      38003
    )
  ),
  file.path(output_root, "annotation_parameter_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
capture.output(sessionInfo(), file = file.path(output_root, "R_sessionInfo.txt"))
cat("Annotation complete for", nrow(order_table), "samples and", nrow(all_metadata), "cells\n")
print(with(all_metadata, table(sample_id, conservative_label)))
