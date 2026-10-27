suppressPackageStartupMessages({
  library(BiocParallel)
  library(DropletUtils)
  library(igraph)
  library(Matrix)
  library(scater)
  library(scran)
  library(SingleCellExperiment)
  library(SoupX)
})

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 5L) {
  stop("Usage: Rscript script.R <matrix_dir> <adaptive_qc_rds> <output_root> <sample> <seed>")
}

matrix_dir <- arguments[[1]]
adaptive_qc_path <- arguments[[2]]
output_root <- arguments[[3]]
sample_id <- arguments[[4]]
seed <- as.integer(arguments[[5]])

sample_output <- file.path(output_root, sample_id)
dir.create(sample_output, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(sample_output, "SoupX_status.json")

write_status <- function(status, message, extra = list()) {
  payload <- c(
    list(
      sample_id = sample_id,
      status = status,
      message = message,
      timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    ),
    extra
  )
  writeLines(
    jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = TRUE, null = "null"),
    status_path
  )
}

write_status("started", "SoupX workflow started")
set.seed(seed)

tryCatch({
  raw_sce <- read10xCounts(
    samples = matrix_dir,
    col.names = TRUE,
    type = "sparse",
    BPPARAM = SerialParam()
  )
  raw_counts <- as(counts(raw_sce), "dgCMatrix")
  qc_result <- readRDS(adaptive_qc_path)
  final_barcodes <- qc_result$cell_data$barcode[qc_result$cell_data$passes_final_qc]

  missing_barcodes <- setdiff(final_barcodes, colnames(raw_counts))
  if (length(missing_barcodes) > 0L) {
    stop(length(missing_barcodes), " final singlet barcodes were absent from the raw matrix")
  }
  if (length(final_barcodes) < 100L) {
    stop("Fewer than 100 final singlets available for SoupX")
  }

  cell_counts <- raw_counts[, final_barcodes, drop = FALSE]
  cluster_sce <- SingleCellExperiment(list(counts = cell_counts))
  cluster_sce <- logNormCounts(cluster_sce)
  variance_fit <- modelGeneVar(cluster_sce, BPPARAM = SerialParam())

  gene_symbols <- rowData(raw_sce)$Symbol
  if (is.null(gene_symbols)) {
    gene_symbols <- rownames(raw_counts)
  }
  excluded_targets <- c("IL12RB2", "FCRL3", "FCRL3-AS1", "AL356276.1", "CD226", "DOK6")
  nuisance_rows <- grepl(
    "^(MT-|RPS|RPL|HB[ABDEGQZ]|MALAT1$)",
    toupper(gene_symbols)
  )
  eligible_rows <- !nuisance_rows & !gene_symbols %in% excluded_targets
  eligible_variance <- variance_fit[eligible_rows, , drop = FALSE]
  hvg_names <- getTopHVGs(eligible_variance, n = min(3000L, nrow(eligible_variance)))
  if (length(hvg_names) < 500L) {
    stop("Fewer than 500 target-independent HVGs available for preliminary clustering")
  }

  cluster_sce <- runPCA(
    cluster_sce,
    subset_row = hvg_names,
    ncomponents = 30,
    BSPARAM = BiocSingular::IrlbaParam(),
    BPPARAM = SerialParam()
  )
  graph <- buildSNNGraph(cluster_sce, use.dimred = "PCA", k = 20)
  clusters <- factor(membership(cluster_louvain(graph)))
  names(clusters) <- colnames(cluster_sce)
  if (nlevels(clusters) < 2L) {
    stop("Preliminary target-independent clustering produced fewer than two clusters")
  }

  soup_channel <- SoupChannel(tod = raw_counts, toc = cell_counts, calcSoupProfile = TRUE)
  soup_channel <- setClusters(soup_channel, clusters)
  soup_channel <- autoEstCont(soup_channel, doPlot = FALSE)
  corrected_counts <- adjustCounts(soup_channel, roundToInt = TRUE)

  contamination_table <- data.frame(
    barcode = rownames(soup_channel$metaData),
    sample_id = sample_id,
    preliminary_cluster = as.character(clusters[rownames(soup_channel$metaData)]),
    contamination_fraction = soup_channel$metaData$rho,
    stringsAsFactors = FALSE
  )
  write.table(
    contamination_table,
    file.path(sample_output, "barcode_contamination.tsv.gz"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  saveRDS(corrected_counts, file.path(sample_output, "SoupX_corrected_counts.rds"), compress = "gzip")
  saveRDS(
    list(
      clusters = clusters,
      hvg_genes = hvg_names,
      excluded_target_genes = excluded_targets,
      soup_profile = soup_channel$soupProfile,
      contamination = contamination_table,
      seed = seed
    ),
    file.path(sample_output, "SoupX_provenance.rds"),
    compress = "gzip"
  )

  summary_table <- data.frame(
    sample_id = sample_id,
    raw_barcodes = ncol(raw_counts),
    final_singlets = ncol(cell_counts),
    preliminary_clusters = nlevels(clusters),
    median_contamination_fraction = median(contamination_table$contamination_fraction),
    min_contamination_fraction = min(contamination_table$contamination_fraction),
    max_contamination_fraction = max(contamination_table$contamination_fraction),
    raw_final_singlet_UMI = sum(cell_counts),
    corrected_final_singlet_UMI = sum(corrected_counts),
    stringsAsFactors = FALSE
  )
  write.table(
    summary_table,
    file.path(sample_output, "SoupX_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  sink(file.path(sample_output, "R_sessionInfo.txt"))
  sessionInfo()
  sink()
  write_status(
    "complete",
    "SoupX autoEstCont and count correction completed",
    as.list(summary_table[1, ])
  )
  cat(jsonlite::toJSON(summary_table, dataframe = "rows", auto_unbox = TRUE, pretty = TRUE), "\n")
}, error = function(error) {
  write_status(
    "failed",
    conditionMessage(error),
    list(fallback_contamination_fraction_applied = FALSE)
  )
  stop(error)
})
