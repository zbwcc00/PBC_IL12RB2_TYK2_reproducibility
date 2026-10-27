suppressPackageStartupMessages({
  library(BiocParallel)
  library(DropletUtils)
  library(Matrix)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(scuttle)
})

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 6L) {
  stop("Usage: Rscript script.R <matrix_dir> <output_root> <sample> <donor> <condition> <seed>")
}

matrix_dir <- arguments[[1]]
output_root <- arguments[[2]]
sample_id <- arguments[[3]]
donor_id <- arguments[[4]]
condition <- arguments[[5]]
seed <- as.integer(arguments[[6]])

dir.create(file.path(output_root, "objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "cell_metrics"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "sample_summary"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "plots"), recursive = TRUE, showWarnings = FALSE)

object_path <- file.path(output_root, "objects", paste0(sample_id, "_adaptive_called_cells.rds"))
metrics_path <- file.path(output_root, "cell_metrics", paste0(sample_id, "_called_cell_metrics.tsv.gz"))
summary_path <- file.path(output_root, "sample_summary", paste0(sample_id, "_summary.tsv"))
threshold_path <- file.path(output_root, "sample_summary", paste0(sample_id, "_thresholds.tsv"))
plot_path <- file.path(output_root, "plots", paste0(sample_id, "_QC_thresholds.png"))

set.seed(seed)
sce <- read10xCounts(
  samples = matrix_dir,
  col.names = TRUE,
  type = "sparse",
  BPPARAM = SerialParam()
)
count_matrix <- counts(sce)
gene_symbols <- rowData(sce)$Symbol
if (is.null(gene_symbols)) {
  gene_symbols <- rownames(sce)
}

total_counts <- Matrix::colSums(count_matrix)
detected_genes <- Matrix::colSums(count_matrix > 0)
mitochondrial_rows <- grepl("^MT-", toupper(gene_symbols))
mitochondrial_counts <- Matrix::colSums(count_matrix[mitochondrial_rows, , drop = FALSE])
mitochondrial_percent <- ifelse(total_counts > 0, 100 * mitochondrial_counts / total_counts, 0)

barcode_ranks <- barcodeRanks(count_matrix, lower = 100)
rank_metadata <- metadata(barcode_ranks)
empty_drops <- emptyDrops(
  count_matrix,
  lower = 100,
  niters = 10000,
  BPPARAM = SerialParam()
)

emptydrops_called <- !is.na(empty_drops$FDR) & empty_drops$FDR <= 0.01
knee_value <- as.numeric(rank_metadata$knee)
inflection_value <- as.numeric(rank_metadata$inflection)
knee_called <- rep(FALSE, length(total_counts))
if (length(knee_value) == 1L && is.finite(knee_value)) {
  knee_called <- total_counts >= knee_value
}
cell_called <- emptydrops_called | knee_called
if (sum(cell_called) < 100L) {
  stop("Fewer than 100 cell-called barcodes for ", sample_id)
}

called_total <- total_counts[cell_called]
called_genes <- detected_genes[cell_called]
called_mt <- mitochondrial_percent[cell_called]

low_total <- isOutlier(called_total, nmads = 3, type = "lower", log = TRUE)
low_genes <- isOutlier(called_genes, nmads = 3, type = "lower", log = TRUE)
high_total <- isOutlier(called_total, nmads = 4, type = "higher", log = TRUE)
high_genes <- isOutlier(called_genes, nmads = 4, type = "higher", log = TRUE)
high_mt <- isOutlier(called_mt, nmads = 3, type = "higher", log = FALSE)
passes_core_qc <- !(low_total | low_genes | high_mt)
if (sum(passes_core_qc) < 100L) {
  stop("Fewer than 100 barcodes pass adaptive core QC for ", sample_id)
}

called_indices <- which(cell_called)
core_indices <- called_indices[passes_core_qc]
doublet_input <- sce[, core_indices]
row_keep <- Matrix::rowSums(counts(doublet_input) > 0) > 0
doublet_work <- doublet_input[row_keep, ]
doublet_work$sample_id <- sample_id
doublet_work <- scDblFinder(
  doublet_work,
  clusters = FALSE,
  BPPARAM = SerialParam(),
  verbose = FALSE
)

doublet_score <- rep(NA_real_, sum(cell_called))
doublet_class <- rep(NA_character_, sum(cell_called))
doublet_score[passes_core_qc] <- doublet_work$scDblFinder.score
doublet_class[passes_core_qc] <- as.character(doublet_work$scDblFinder.class)
passes_final_qc <- passes_core_qc & !is.na(doublet_class) & doublet_class == "singlet"

call_source <- ifelse(
  emptydrops_called[cell_called] & knee_called[cell_called],
  "emptyDrops_FDR01_and_knee",
  ifelse(emptydrops_called[cell_called], "emptyDrops_FDR01", "knee_rescue")
)
cell_metrics <- data.frame(
  barcode = colnames(sce)[cell_called],
  sample_id = sample_id,
  donor_id = donor_id,
  condition = condition,
  total_UMI = called_total,
  detected_genes = called_genes,
  mitochondrial_UMI = mitochondrial_counts[cell_called],
  percent_mt = called_mt,
  emptyDrops_FDR = empty_drops$FDR[cell_called],
  emptyDrops_Limited = empty_drops$Limited[cell_called],
  cell_call_source = call_source,
  low_total_flag = as.logical(low_total),
  low_genes_flag = as.logical(low_genes),
  high_total_review_flag = as.logical(high_total),
  high_genes_review_flag = as.logical(high_genes),
  high_mt_flag = as.logical(high_mt),
  passes_core_qc = passes_core_qc,
  scDblFinder_score = doublet_score,
  scDblFinder_class = doublet_class,
  passes_final_qc = passes_final_qc,
  stringsAsFactors = FALSE
)

threshold_value <- function(flag, direction) {
  thresholds <- attr(flag, "thresholds")
  if (is.null(thresholds) || !direction %in% names(thresholds)) {
    return(NA_real_)
  }
  as.numeric(thresholds[[direction]])
}

threshold_table <- data.frame(
  metric = c(
    "barcodeRanks_knee_total_UMI",
    "barcodeRanks_inflection_total_UMI",
    "adaptive_total_UMI_lower",
    "adaptive_total_UMI_upper_review",
    "adaptive_detected_genes_lower",
    "adaptive_detected_genes_upper_review",
    "adaptive_percent_mt_upper"
  ),
  value = c(
    knee_value,
    inflection_value,
    threshold_value(low_total, "lower"),
    threshold_value(high_total, "higher"),
    threshold_value(low_genes, "lower"),
    threshold_value(high_genes, "higher"),
    threshold_value(high_mt, "higher")
  ),
  use = c(
    "cell_call_rescue",
    "review_only",
    "hard_core_QC",
    "review_only_doublet_context",
    "hard_core_QC",
    "review_only_doublet_context",
    "hard_core_QC"
  ),
  stringsAsFactors = FALSE
)

summary_table <- data.frame(
  sample_id = sample_id,
  donor_id = donor_id,
  condition = condition,
  raw_corrected_barcodes = ncol(sce),
  raw_total_UMI = sum(total_counts),
  emptyDrops_tested = sum(!is.na(empty_drops$FDR)),
  emptyDrops_FDR01 = sum(emptydrops_called),
  knee_rescue_total = sum(knee_called),
  union_cell_called = sum(cell_called),
  low_total_flagged = sum(low_total),
  low_genes_flagged = sum(low_genes),
  high_mt_flagged = sum(high_mt),
  passes_core_qc = sum(passes_core_qc),
  scDblFinder_doublets = sum(doublet_class == "doublet", na.rm = TRUE),
  final_singlets = sum(passes_final_qc),
  final_singlet_UMI = sum(called_total[passes_final_qc]),
  final_singlet_median_UMI = median(called_total[passes_final_qc]),
  final_singlet_median_genes = median(called_genes[passes_final_qc]),
  final_singlet_median_percent_mt = median(called_mt[passes_final_qc]),
  stringsAsFactors = FALSE
)

connection <- gzfile(metrics_path, open = "wt")
write.table(cell_metrics, connection, sep = "\t", quote = FALSE, row.names = FALSE)
close(connection)
write.table(summary_table, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(threshold_table, threshold_path, sep = "\t", quote = FALSE, row.names = FALSE)

raw_called_counts <- as(count_matrix[, called_indices, drop = FALSE], "dgCMatrix")
result <- list(
  counts = raw_called_counts,
  feature_data = as.data.frame(rowData(sce)),
  cell_data = cell_metrics,
  thresholds = threshold_table,
  summary = summary_table,
  provenance = list(
    emptyDrops_lower = 100,
    emptyDrops_niters = 10000,
    emptyDrops_FDR = 0.01,
    cell_call_rule = "emptyDrops FDR <= 0.01 union barcodeRanks knee",
    core_QC_rule = "sample-adaptive 3 MAD lower total UMI and genes plus 3 MAD upper percent.mt",
    upper_UMI_gene_rule = "4 MAD upper flags retained for review; scDblFinder makes doublet call",
    doublet_method = "scDblFinder per sample",
    seed = seed
  )
)
saveRDS(result, object_path, compress = "gzip")

png(plot_path, width = 2200, height = 1800, res = 220)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
hist(log10(called_total), breaks = 60, main = paste(sample_id, "called-barcode UMI"), xlab = "log10(total UMI)")
abline(v = log10(threshold_value(low_total, "lower")), col = "#D7301F", lwd = 2, lty = 2)
hist(log10(called_genes), breaks = 60, main = "Detected genes", xlab = "log10(detected genes)")
abline(v = log10(threshold_value(low_genes, "lower")), col = "#D7301F", lwd = 2, lty = 2)
hist(called_mt, breaks = 60, main = "Mitochondrial fraction", xlab = "percent.mt")
abline(v = threshold_value(high_mt, "higher"), col = "#D7301F", lwd = 2, lty = 2)
plot(
  log10(called_total),
  log10(called_genes),
  pch = 16,
  cex = 0.35,
  col = ifelse(passes_final_qc, "#238B45", "#BDBDBD"),
  xlab = "log10(total UMI)",
  ylab = "log10(detected genes)",
  main = "Final singlets (green)"
)
dev.off()

sink(file.path(output_root, "R_sessionInfo.txt"))
sessionInfo()
sink()

cat(jsonlite::toJSON(summary_table, dataframe = "rows", auto_unbox = TRUE, pretty = TRUE), "\n")
