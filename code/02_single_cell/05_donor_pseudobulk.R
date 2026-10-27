suppressPackageStartupMessages({
  library(edgeR)
  library(Matrix)
})

arguments <- commandArgs(trailingOnly = TRUE)
set.seed(39003)
if (length(arguments) != 0L) stop("Usage: Rscript script.R")
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_argument), winslash = "/")
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/")
work_root <- file.path(repo_root, "work", "single_cell", "hra008003")
matrix_root <- file.path(work_root, "reconstruction")
adaptive_root <- file.path(work_root, "qc")
soupx_root <- file.path(work_root, "soupx")
annotation_root <- file.path(work_root, "annotation")
output_root <- file.path(work_root, "pseudobulk")
dir.create(file.path(output_root, "all_gene_results"), recursive = TRUE, showWarnings = FALSE)

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

target_genes <- c(
  "IL12RB2", "FCRL3", "AL356276.1", "CD226", "DOK6",
  "IL12RB1", "JAK2", "STAT4", "IFNG"
)
candidate_chain <- c("IL12RB2", "FCRL3", "AL356276.1", "CD226", "DOK6")
display_name <- c(
  IL12RB2 = "IL12RB2",
  FCRL3 = "FCRL3",
  AL356276.1 = "FCRL3-AS1 (AL356276.1)",
  CD226 = "CD226",
  DOK6 = "DOK6",
  IL12RB1 = "IL12RB1",
  JAK2 = "JAK2",
  STAT4 = "STAT4",
  IFNG = "IFNG"
)

sample_feature_tables <- list()
for (sample_id in order_table$sample_id) {
  qc_object <- readRDS(file.path(
    adaptive_root,
    "objects",
    paste0(sample_id, "_adaptive_called_cells.rds")
  ))
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
feature_table <- sample_feature_tables[[1]][
  match(common_feature_ids, sample_feature_tables[[1]]$feature_id),
  ,
  drop = FALSE
]
target_feature_ids <- lapply(target_genes, function(gene) {
  feature_table$feature_id[feature_table$gene_symbol == gene]
})
names(target_feature_ids) <- target_genes
missing_targets <- names(target_feature_ids)[lengths(target_feature_ids) == 0L]
if (length(missing_targets) > 0L) {
  stop("Target genes absent from feature universe: ", paste(missing_targets, collapse = ", "))
}
all_target_ids <- unique(unlist(target_feature_ids, use.names = FALSE))

schemes <- list(
  primary_conservative = "conservative_label",
  pure_TNK_sensitivity = "tnk_refined_label"
)
primary_labels <- c("lrNK", "cNK", "CD4_T", "CD8_T", "T_other", "TNK_mixed")
sensitivity_labels <- c("lrNK", "cNK", "CD4_T", "CD8_T", "T_other", "TNK_unresolved")

pseudobulk_columns <- list()
pseudobulk_metadata <- list()
capture_rows <- list()

for (row_index in seq_len(nrow(order_table))) {
  record <- order_table[row_index, ]
  sample_id <- record$sample_id
  cat("Aggregating", sample_id, "\n")
  corrected_counts <- readRDS(file.path(soupx_root, sample_id, "SoupX_corrected_counts.rds"))
  sample_features <- sample_feature_tables[[sample_id]]
  if (!identical(rownames(corrected_counts), sample_features$feature_id)) {
    stop("Corrected-count and adaptive-QC feature rows differ within ", sample_id)
  }
  if (!all(common_feature_ids %in% rownames(corrected_counts))) {
    stop("Common feature universe is incomplete for ", sample_id)
  }
  corrected_counts <- corrected_counts[common_feature_ids, , drop = FALSE]
  metadata <- read.delim(
    gzfile(file.path(annotation_root, "sample_metadata", paste0(sample_id, "_cell_metadata.tsv.gz"))),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  metadata <- metadata[match(colnames(corrected_counts), metadata$barcode), , drop = FALSE]
  if (!identical(metadata$barcode, colnames(corrected_counts))) {
    stop("Annotation metadata and corrected counts are misaligned for ", sample_id)
  }
  stopifnot(identical(rownames(corrected_counts), feature_table$feature_id))
  total_cell_umi <- Matrix::colSums(corrected_counts)
  depth_filters <- list(
    all_final_singlets = rep(TRUE, ncol(corrected_counts)),
    at_least_500_UMI = total_cell_umi >= 500
  )
  target_cell_counts <- lapply(target_feature_ids, function(ids) {
    Matrix::colSums(corrected_counts[ids, , drop = FALSE])
  })

  for (scheme in names(schemes)) {
    label_column <- schemes[[scheme]]
    allowed_labels <- if (scheme == "primary_conservative") primary_labels else sensitivity_labels
    base_labels <- as.character(metadata[[label_column]])
    base_labels[base_labels == "NA" | base_labels == ""] <- NA_character_
    base_labels[!base_labels %in% allowed_labels] <- NA_character_

    for (cell_filter in names(depth_filters)) {
      labels <- base_labels
      labels[!depth_filters[[cell_filter]]] <- NA_character_
      observed_labels <- intersect(allowed_labels, unique(labels[!is.na(labels)]))
      if (length(observed_labels) == 0L) {
        next
      }

      group_matrix <- sparseMatrix(
        i = which(!is.na(labels)),
        j = match(labels[!is.na(labels)], observed_labels),
        x = 1,
        dims = c(ncol(corrected_counts), length(observed_labels)),
        dimnames = list(colnames(corrected_counts), observed_labels)
      )
      aggregated <- corrected_counts %*% group_matrix

      for (label_index in seq_along(observed_labels)) {
        current_label <- observed_labels[[label_index]]
        key <- paste(scheme, cell_filter, current_label, sep = "__")
        cell_indices <- which(labels == current_label)
        pseudobulk_columns[[key]][[sample_id]] <- aggregated[, label_index, drop = FALSE]
        pseudobulk_metadata[[key]][[sample_id]] <- data.frame(
          sample_id = sample_id,
          donor_id = record$donor_id,
          condition = record$condition,
          n_cells = length(cell_indices),
          total_UMI = sum(total_cell_umi[cell_indices]),
          stringsAsFactors = FALSE
        )

        for (target_gene in target_genes) {
          target_values <- target_cell_counts[[target_gene]][cell_indices]
          capture_rows[[length(capture_rows) + 1L]] <- data.frame(
            scheme = scheme,
            cell_filter = cell_filter,
            cell_type = current_label,
            sample_id = sample_id,
            donor_id = record$donor_id,
            condition = record$condition,
            gene_symbol = target_gene,
            display_gene = unname(display_name[target_gene]),
            n_cells = length(cell_indices),
            positive_cells = sum(target_values > 0),
            positive_fraction = mean(target_values > 0),
            target_UMI = sum(target_values),
            total_UMI = sum(total_cell_umi[cell_indices]),
            target_CPM = if (sum(total_cell_umi[cell_indices]) > 0) {
              1e6 * sum(target_values) / sum(total_cell_umi[cell_indices])
            } else {
              NA_real_
            },
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  rm(corrected_counts, metadata, group_matrix, aggregated)
  gc()
}

capture_table <- do.call(rbind, capture_rows)
rownames(capture_table) <- NULL
write.table(
  capture_table,
  file.path(output_root, "01_donor_celltype_target_capture.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

exact_permutation <- function(values, conditions) {
  valid <- is.finite(values) & !is.na(conditions)
  values <- values[valid]
  conditions <- conditions[valid]
  n_pbc <- sum(conditions == "PBC")
  if (n_pbc == 0L || n_pbc == length(values)) {
    return(c(effect = NA_real_, p_value = NA_real_, assignments = 0))
  }
  assignments <- combn(seq_along(values), n_pbc)
  observed <- mean(values[conditions == "PBC"]) - mean(values[conditions == "control"])
  null_effects <- apply(assignments, 2, function(indices) {
    mean(values[indices]) - mean(values[-indices])
  })
  c(
    effect = observed,
    p_value = mean(abs(null_effects) >= abs(observed) - 1e-12),
    assignments = ncol(assignments)
  )
}

effect_ci <- function(values, conditions) {
  control_values <- values[conditions == "control" & is.finite(values)]
  pbc_values <- values[conditions == "PBC" & is.finite(values)]
  if (length(control_values) < 2L || length(pbc_values) < 2L) {
    return(c(effect = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }
  effect <- mean(pbc_values) - mean(control_values)
  test <- tryCatch(
    t.test(pbc_values, control_values, var.equal = FALSE),
    error = function(error) NULL
  )
  if (is.null(test)) {
    return(c(effect = effect, ci_low = NA_real_, ci_high = NA_real_))
  }
  c(
    effect = effect,
    ci_low = unname(test$conf.int[[1]]),
    ci_high = unname(test$conf.int[[2]])
  )
}

availability_rows <- list()
target_result_rows <- list()
donor_logcpm_rows <- list()
permutation_rows <- list()
ci_rows <- list()

for (key in names(pseudobulk_columns)) {
  key_parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  scheme <- key_parts[[1]]
  cell_filter <- key_parts[[2]]
  cell_type <- key_parts[[3]]
  columns <- pseudobulk_columns[[key]]
  metadata_rows <- pseudobulk_metadata[[key]]
  donor_metadata <- do.call(rbind, metadata_rows)
  rownames(donor_metadata) <- NULL
  counts_matrix <- do.call(cbind, columns)
  colnames(counts_matrix) <- names(columns)
  donor_metadata <- donor_metadata[match(colnames(counts_matrix), donor_metadata$sample_id), , drop = FALSE]

  for (minimum_cells in c(20L, 10L)) {
    eligible <- donor_metadata$n_cells >= minimum_cells
    eligible_counts <- counts_matrix[, eligible, drop = FALSE]
    eligible_metadata <- donor_metadata[eligible, , drop = FALSE]
    n_control <- sum(eligible_metadata$condition == "control")
    n_pbc <- sum(eligible_metadata$condition == "PBC")
    exploratory_estimable <- n_control >= 2L && n_pbc >= 2L
    primary_donor_support <- n_control >= 3L && n_pbc >= 3L
    availability_rows[[length(availability_rows) + 1L]] <- data.frame(
      scheme = scheme,
      cell_filter = cell_filter,
      cell_type = cell_type,
      minimum_cells = minimum_cells,
      eligible_control_donors = n_control,
      eligible_PBC_donors = n_pbc,
      eligible_cells = sum(eligible_metadata$n_cells),
      exploratory_estimable = exploratory_estimable,
      primary_donor_support = primary_donor_support,
      stringsAsFactors = FALSE
    )
    if (!exploratory_estimable) {
      next
    }

    condition <- factor(eligible_metadata$condition, levels = c("control", "PBC"))
    design <- model.matrix(~condition)
    initial_dge <- DGEList(counts = eligible_counts)
    passes_filter <- filterByExpr(initial_dge, design = design, min.count = 5)
    target_present <- rownames(eligible_counts) %in% all_target_ids & Matrix::rowSums(eligible_counts) > 0
    keep <- passes_filter | target_present
    dge <- initial_dge[keep, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge)

    fit_result <- tryCatch({
      dge <- estimateDisp(dge, design, robust = TRUE)
      fit <- glmQLFit(dge, design, robust = TRUE)
      test <- glmQLFTest(fit, coef = "conditionPBC")
      result <- topTags(test, n = Inf, sort.by = "none")$table
      result$feature_id <- rownames(result)
      result$gene_symbol <- feature_table$gene_symbol[match(result$feature_id, feature_table$feature_id)]
      result$passes_filterByExpr <- passes_filter[match(result$feature_id, rownames(eligible_counts))]
      result$forced_prespecified_target <- !result$passes_filterByExpr & result$feature_id %in% all_target_ids
      result
    }, error = function(error) {
      warning("edgeR failed for ", key, " minimum_cells=", minimum_cells, ": ", conditionMessage(error))
      NULL
    })
    if (is.null(fit_result)) {
      next
    }

    fit_result$scheme <- scheme
    fit_result$cell_filter <- cell_filter
    fit_result$cell_type <- cell_type
    fit_result$minimum_cells <- minimum_cells
    fit_result$n_control_donors <- n_control
    fit_result$n_PBC_donors <- n_pbc
    fit_result$primary_donor_support <- primary_donor_support
    fit_result$contrast <- "PBC_minus_control"
    safe_key <- gsub("[^A-Za-z0-9_]+", "_", key)
    write.table(
      fit_result,
      file.path(output_root, "all_gene_results", paste0(safe_key, "_min", minimum_cells, "_edgeR.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    target_result <- fit_result[fit_result$gene_symbol %in% target_genes, , drop = FALSE]
    target_result$display_gene <- unname(display_name[target_result$gene_symbol])
    missing_fit_targets <- setdiff(target_genes, unique(target_result$gene_symbol))
    if (length(missing_fit_targets) > 0L) {
      placeholder <- fit_result[rep(NA_integer_, length(missing_fit_targets)), , drop = FALSE]
      placeholder$feature_id <- NA_character_
      placeholder$gene_symbol <- missing_fit_targets
      placeholder$passes_filterByExpr <- FALSE
      placeholder$forced_prespecified_target <- FALSE
      placeholder$scheme <- scheme
      placeholder$cell_filter <- cell_filter
      placeholder$cell_type <- cell_type
      placeholder$minimum_cells <- minimum_cells
      placeholder$n_control_donors <- n_control
      placeholder$n_PBC_donors <- n_pbc
      placeholder$primary_donor_support <- primary_donor_support
      placeholder$contrast <- "PBC_minus_control"
      placeholder$display_gene <- unname(display_name[missing_fit_targets])
      target_result <- rbind(target_result, placeholder)
    }
    target_result_rows[[length(target_result_rows) + 1L]] <- target_result

    log_cpm <- cpm(dge, log = TRUE, prior.count = 0.5)
    capture_subset <- capture_table[
      capture_table$scheme == scheme &
        capture_table$cell_filter == cell_filter &
        capture_table$cell_type == cell_type &
        capture_table$sample_id %in% eligible_metadata$sample_id,
      ,
      drop = FALSE
    ]
    for (target_gene in target_genes) {
      ids <- intersect(target_feature_ids[[target_gene]], rownames(log_cpm))
      gene_logcpm <- if (length(ids) == 0L) {
        rep(NA_real_, ncol(log_cpm))
      } else if (length(ids) == 1L) {
        as.numeric(log_cpm[ids, ])
      } else {
        log2(Matrix::colSums(2^log_cpm[ids, , drop = FALSE]) + 0.5)
      }
      donor_rows <- data.frame(
        scheme = scheme,
        cell_filter = cell_filter,
        cell_type = cell_type,
        minimum_cells = minimum_cells,
        sample_id = eligible_metadata$sample_id,
        donor_id = eligible_metadata$donor_id,
        condition = eligible_metadata$condition,
        gene_symbol = target_gene,
        display_gene = unname(display_name[target_gene]),
        logCPM = gene_logcpm,
        stringsAsFactors = FALSE
      )
      donor_logcpm_rows[[length(donor_logcpm_rows) + 1L]] <- donor_rows

      target_capture <- capture_subset[capture_subset$gene_symbol == target_gene, , drop = FALSE]
      target_capture <- target_capture[match(eligible_metadata$sample_id, target_capture$sample_id), , drop = FALSE]
      metrics <- list(
        logCPM = gene_logcpm,
        log2_CPM_plus_0.5 = log2(target_capture$target_CPM + 0.5),
        positive_fraction = target_capture$positive_fraction
      )
      for (metric_name in names(metrics)) {
        values <- metrics[[metric_name]]
        permutation <- exact_permutation(values, eligible_metadata$condition)
        confidence <- effect_ci(values, eligible_metadata$condition)
        permutation_rows[[length(permutation_rows) + 1L]] <- data.frame(
          scheme = scheme,
          cell_filter = cell_filter,
          cell_type = cell_type,
          minimum_cells = minimum_cells,
          gene_symbol = target_gene,
          display_gene = unname(display_name[target_gene]),
          metric = metric_name,
          effect_PBC_minus_control = permutation[["effect"]],
          exact_PValue = permutation[["p_value"]],
          exact_assignments = permutation[["assignments"]],
          n_control_donors = n_control,
          n_PBC_donors = n_pbc,
          primary_donor_support = primary_donor_support,
          stringsAsFactors = FALSE
        )
        ci_rows[[length(ci_rows) + 1L]] <- data.frame(
          scheme = scheme,
          cell_filter = cell_filter,
          cell_type = cell_type,
          minimum_cells = minimum_cells,
          gene_symbol = target_gene,
          display_gene = unname(display_name[target_gene]),
          metric = metric_name,
          effect_PBC_minus_control = confidence[["effect"]],
          CI_low = confidence[["ci_low"]],
          CI_high = confidence[["ci_high"]],
          n_control_donors = n_control,
          n_PBC_donors = n_pbc,
          primary_donor_support = primary_donor_support,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

availability <- do.call(rbind, availability_rows)
rownames(availability) <- NULL
write.table(
  availability,
  file.path(output_root, "02_pseudobulk_contrast_availability.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

target_results <- do.call(rbind, target_result_rows)
rownames(target_results) <- NULL
target_results$BH_IL12RB2_across_celltype_tests <- NA_real_
target_results$BH_candidate_chain_across_tests <- NA_real_
target_results$BH_pathway_family_across_tests <- NA_real_
for (current_scheme in unique(target_results$scheme)) {
  for (current_filter in unique(target_results$cell_filter)) {
    for (minimum_cells in unique(target_results$minimum_cells)) {
      group <- target_results$scheme == current_scheme &
        target_results$cell_filter == current_filter &
        target_results$minimum_cells == minimum_cells
      il12_group <- group & target_results$gene_symbol == "IL12RB2"
      chain_group <- group & target_results$gene_symbol %in% candidate_chain
      pathway_group <- group & target_results$gene_symbol %in% target_genes
      target_results$BH_IL12RB2_across_celltype_tests[il12_group] <- p.adjust(
        target_results$PValue[il12_group],
        method = "BH"
      )
      target_results$BH_candidate_chain_across_tests[chain_group] <- p.adjust(
        target_results$PValue[chain_group],
        method = "BH"
      )
      target_results$BH_pathway_family_across_tests[pathway_group] <- p.adjust(
        target_results$PValue[pathway_group],
        method = "BH"
      )
    }
  }
}
write.table(
  target_results,
  file.path(output_root, "03_prespecified_target_edgeR_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

donor_logcpm <- do.call(rbind, donor_logcpm_rows)
rownames(donor_logcpm) <- NULL
write.table(
  donor_logcpm,
  file.path(output_root, "04_donor_target_logCPM.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

permutation_results <- do.call(rbind, permutation_rows)
rownames(permutation_results) <- NULL
permutation_results$BH_IL12RB2_across_celltype_tests <- NA_real_
for (current_scheme in unique(permutation_results$scheme)) {
  for (current_filter in unique(permutation_results$cell_filter)) {
    for (minimum_cells in unique(permutation_results$minimum_cells)) {
      for (metric_name in unique(permutation_results$metric)) {
        group <- permutation_results$scheme == current_scheme &
          permutation_results$cell_filter == current_filter &
          permutation_results$minimum_cells == minimum_cells &
          permutation_results$metric == metric_name &
          permutation_results$gene_symbol == "IL12RB2"
        permutation_results$BH_IL12RB2_across_celltype_tests[group] <- p.adjust(
          permutation_results$exact_PValue[group],
          method = "BH"
        )
      }
    }
  }
}
write.table(
  permutation_results,
  file.path(output_root, "05_exact_donor_permutation_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

confidence_results <- do.call(rbind, ci_rows)
rownames(confidence_results) <- NULL
write.table(
  confidence_results,
  file.path(output_root, "06_donor_effect_95CI.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    parameter = c(
      "independent_unit",
      "primary_label_scheme",
      "sensitivity_label_scheme",
      "minimum_cells_primary",
      "minimum_cells_sensitivity",
      "minimum_donors_primary_support",
      "minimum_donors_exploratory_fit",
      "cell_depth_sensitivity",
      "all_gene_filter",
      "prespecified_low_expression_rule",
      "differential_model",
      "exact_test",
      "confidence_interval",
      "multiple_testing",
      "FCRL3_AS1_alias"
    ),
    value = c(
      "donor; one liver sample per donor",
      "conservative SingleR labels; mixed TNK retained and not forced into pure subsets",
      "pure-TNK-only SingleR reclassification; target genes excluded; used only as sensitivity",
      20,
      10,
      "at least 3 control and 3 PBC donors",
      "at least 2 control and 2 PBC donors; explicitly flagged exploratory",
      "all final singlets are primary; cells with at least 500 SoupX-corrected UMI are a fixed sensitivity analysis",
      "edgeR filterByExpr with min.count=5",
      "predeclared targets with any pseudobulk count are retained even if filterByExpr fails; flag preserved",
      "edgeR TMM quasi-likelihood negative binomial GLM; PBC minus control",
      "all donor-label assignments enumerated for logCPM, log2(CPM+0.5), and positive fraction",
      "Welch donor-level mean difference CI; descriptive complement to edgeR",
      "all-gene FDR within each fit plus BH across IL12RB2 cell-type tests, candidate chain, and pathway family",
      "ENSG00000227217 is labeled AL356276.1 in this Cell Ranger reference"
    )
  ),
  file.path(output_root, "parameter_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
capture.output(sessionInfo(), file = file.path(output_root, "R_sessionInfo.txt"))
cat("Donor-level pseudobulk and low-expression robust statistics completed\n")
print(target_results[target_results$gene_symbol == "IL12RB2", c(
  "scheme", "cell_filter", "cell_type", "minimum_cells", "logFC", "PValue", "FDR",
  "BH_IL12RB2_across_celltype_tests", "n_control_donors", "n_PBC_donors",
  "primary_donor_support", "passes_filterByExpr", "forced_prespecified_target"
)])
