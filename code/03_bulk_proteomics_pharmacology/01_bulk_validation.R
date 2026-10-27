options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(limma)
  library(lmtest)
  library(patchwork)
  library(sandwich)
})

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_argument), winslash = "/")
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/")
input_dir <- file.path(repo_root, "data", "raw", "bulk", "input")
result_dir <- file.path(repo_root, "work", "bulk", "validation", "results")
figure_dir <- file.path(repo_root, "work", "bulk", "validation", "figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20261014)

read_expression <- function(path) {
  data <- read.delim(gzfile(path), check.names = FALSE, quote = "", comment.char = "")
  genes <- data[[1]]
  matrix <- as.matrix(data[, -1, drop = FALSE])
  storage.mode(matrix) <- "double"
  rownames(matrix) <- genes
  matrix
}

limma_contrast <- function(expression, group, contrast_text) {
  group <- droplevels(factor(group))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  contrast <- makeContrasts(contrasts = contrast_text, levels = design)
  fit <- lmFit(expression, design)
  fit <- contrasts.fit(fit, contrast)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  table <- topTable(fit, number = Inf, sort.by = "P")
  table$gene_symbol <- rownames(table)
  table$rank_by_p <- seq_len(nrow(table))
  table[, c("gene_symbol", setdiff(names(table), "gene_symbol"))]
}

hedges_g <- function(case, control) {
  case <- case[is.finite(case)]
  control <- control[is.finite(control)]
  n_case <- length(case)
  n_control <- length(control)
  pooled_sd <- sqrt(((n_case - 1) * var(case) + (n_control - 1) * var(control)) /
    (n_case + n_control - 2))
  if (!is.finite(pooled_sd) || pooled_sd == 0) {
    return(c(g = NA_real_, se = NA_real_))
  }
  correction <- 1 - 3 / (4 * (n_case + n_control) - 9)
  g <- correction * (mean(case) - mean(control)) / pooled_sd
  se <- sqrt((n_case + n_control) / (n_case * n_control) + g^2 / (2 * (n_case + n_control - 2)))
  c(g = g, se = se)
}

permutation_mean_test <- function(values, labels, case_label, control_label, max_exact = 100000, b = 100000) {
  keep <- labels %in% c(case_label, control_label) & is.finite(values)
  values <- values[keep]
  labels <- droplevels(factor(labels[keep], levels = c(control_label, case_label)))
  n_case <- sum(labels == case_label)
  n_total <- length(values)
  observed <- mean(values[labels == case_label]) - mean(values[labels == control_label])
  combinations <- choose(n_total, n_case)
  if (is.finite(combinations) && combinations <= max_exact) {
    indices <- combn(n_total, n_case)
    case_sums <- colSums(matrix(values[indices], nrow = n_case))
    permuted <- case_sums / n_case - (sum(values) - case_sums) / (n_total - n_case)
    method <- "exact_enumeration"
  } else {
    permuted <- replicate(b, {
      shuffled <- sample(labels)
      mean(values[shuffled == case_label]) - mean(values[shuffled == control_label])
    })
    method <- paste0("monte_carlo_", b)
  }
  c(
    observed_difference = observed,
    p_two_sided = (sum(abs(permuted) >= abs(observed)) + 1) / (length(permuted) + 1),
    p_direction_case_greater = (sum(permuted >= observed) + 1) / (length(permuted) + 1),
    permutations = length(permuted),
    method = method
  )
}

compare_vector <- function(values, labels, case_label, control_label, dataset, contrast) {
  keep <- labels %in% c(case_label, control_label) & is.finite(values)
  values <- values[keep]
  labels <- labels[keep]
  case <- values[labels == case_label]
  control <- values[labels == control_label]
  effect <- hedges_g(case, control)
  permutation <- permutation_mean_test(values, labels, case_label, control_label)
  data.frame(
    dataset = dataset,
    contrast = contrast,
    case_label = case_label,
    control_label = control_label,
    n_case = length(case),
    n_control = length(control),
    mean_case = mean(case),
    mean_control = mean(control),
    median_case = median(case),
    median_control = median(control),
    mean_difference = mean(case) - mean(control),
    hedges_g = effect[["g"]],
    hedges_g_se = effect[["se"]],
    wilcoxon_p = suppressWarnings(wilcox.test(case, control, exact = FALSE)$p.value),
    permutation_p_two_sided = as.numeric(permutation[["p_two_sided"]]),
    permutation_p_direction_case_greater = as.numeric(permutation[["p_direction_case_greater"]]),
    permutation_method = as.character(permutation[["method"]]),
    stringsAsFactors = FALSE
  )
}

signature_score <- function(expression, genes) {
  covered <- intersect(genes, rownames(expression))
  if (length(covered) < 2) {
    return(list(score = rep(NA_real_, ncol(expression)), genes = covered))
  }
  z <- t(scale(t(expression[covered, , drop = FALSE])))
  list(score = colMeans(z, na.rm = TRUE), genes = covered)
}

write_table <- function(data, name) {
  write.table(data, file.path(result_dir, name), sep = "\t", quote = FALSE, row.names = FALSE)
}

save_png <- function(plot, path, width, height, dpi = 300) {
  grDevices::png(path, width = width * dpi, height = height * dpi, res = dpi, type = "cairo-png")
  print(plot)
  grDevices::dev.off()
}

bind_rows_fill <- function(...) {
  frames <- list(...)
  columns <- unique(unlist(lapply(frames, names)))
  frames <- lapply(frames, function(frame) {
    missing <- setdiff(columns, names(frame))
    for (column in missing) {
      frame[[column]] <- NA
    }
    frame[, columns, drop = FALSE]
  })
  do.call(rbind, frames)
}

message("GSE119600: limma and signatures")
expr_119600 <- read_expression(file.path(input_dir, "GSE119600_gene_log2expr.tsv.gz"))
meta_119600 <- read.delim(file.path(input_dir, "GSE119600_metadata.tsv"), check.names = FALSE)
stopifnot(identical(colnames(expr_119600), meta_119600$sample_id))
group_119600 <- ifelse(meta_119600$source_name_ch1 == "Primary biliary cholangitis, adult", "PBC", "Control")
de_119600 <- limma_contrast(expr_119600, factor(group_119600, levels = c("Control", "PBC")), "PBC-Control")
write_table(de_119600, "01_GSE119600_limma_PBC_vs_control.tsv")

target_comparisons <- compare_vector(
  expr_119600["IL12RB2", ], group_119600, "PBC", "Control", "GSE119600_whole_blood", "PBC_vs_control"
)
target_limma <- de_119600[de_119600$gene_symbol == "IL12RB2", ]
target_comparisons$limma_log2FC <- target_limma$logFC
target_comparisons$limma_t <- target_limma$t
target_comparisons$limma_p <- target_limma$P.Value
target_comparisons$limma_fdr <- target_limma$adj.P.Val
target_comparisons$limma_rank <- target_limma$rank_by_p

signatures <- list(
  lrNK = c("CXCR6", "CD69", "ITGAE", "XCL1", "XCL2", "KLRD1", "NKG7"),
  cNK = c("FCGR3A", "S1PR5", "FGFBP2", "KLF2"),
  IL12_IFNG = c("IL12RB1", "IL12RB2", "STAT4", "IFNG", "TBX21", "STAT1", "IRF1", "CXCL9", "CXCL10"),
  cytotoxic_NK = c("KLRD1", "NKG7", "GNLY", "PRF1", "GZMB", "FCGR3A", "TYROBP")
)
signature_matrix <- matrix(nrow = 0, ncol = ncol(expr_119600), dimnames = list(NULL, colnames(expr_119600)))
signature_coverage <- data.frame()
for (signature_name in names(signatures)) {
  scored <- signature_score(expr_119600, signatures[[signature_name]])
  signature_matrix <- rbind(signature_matrix, scored$score)
  rownames(signature_matrix)[nrow(signature_matrix)] <- signature_name
  signature_coverage <- rbind(
    signature_coverage,
    data.frame(
      signature = signature_name,
      n_prespecified = length(signatures[[signature_name]]),
      n_covered = length(scored$genes),
      covered_genes = paste(scored$genes, collapse = ";")
    )
  )
}
signature_tests <- limma_contrast(signature_matrix, factor(group_119600, levels = c("Control", "PBC")), "PBC-Control")
signature_tests$BH_FDR_within_signatures <- p.adjust(signature_tests$P.Value, method = "BH")
write_table(signature_coverage, "02_signature_coverage.tsv")
write_table(signature_tests, "03_GSE119600_signature_tests.tsv")
write.table(
  data.frame(sample_id = colnames(signature_matrix), group = group_119600, t(signature_matrix), check.names = FALSE),
  file.path(result_dir, "04_GSE119600_signature_scores.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

signature_correlations <- do.call(rbind, lapply(rownames(signature_matrix), function(signature_name) {
  test <- cor.test(expr_119600["IL12RB2", ], signature_matrix[signature_name, ], method = "spearman", exact = FALSE)
  data.frame(signature = signature_name, rho = unname(test$estimate), p = test$p.value)
}))
signature_correlations$BH_FDR <- p.adjust(signature_correlations$p, method = "BH")
write_table(signature_correlations, "05_IL12RB2_signature_correlations.tsv")

message("GSE119600: xCell immune enrichment")
xcell_status <- "not_run"
xcell_scores <- tryCatch({
  xcell_environment <- new.env(parent = emptyenv())
  utils::data("xCell.data", package = "xCell", envir = xcell_environment)
  xcell_data <- xcell_environment$xCell.data
  scores <- xCell::xCellAnalysis(
    expr_119600,
    signatures = xcell_data$signatures,
    genes = xcell_data$genes,
    spill = xcell_data$spill.array,
    rnaseq = FALSE,
    scale = TRUE,
    save.raw = FALSE,
    parallel.sz = 1
  )
  xcell_status <<- "success"
  scores
}, error = function(error) {
  xcell_status <<- paste("failed", conditionMessage(error))
  NULL
})
writeLines(xcell_status, file.path(result_dir, "06_xcell_status.txt"))

if (!is.null(xcell_scores)) {
  write.table(
    data.frame(cell_type = rownames(xcell_scores), xcell_scores, check.names = FALSE),
    file.path(result_dir, "07_GSE119600_xcell_scores.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  xcell_tests <- limma_contrast(xcell_scores, factor(group_119600, levels = c("Control", "PBC")), "PBC-Control")
  write_table(xcell_tests, "08_GSE119600_xcell_group_tests.tsv")
  xcell_correlations <- do.call(rbind, lapply(rownames(xcell_scores), function(cell_type) {
    test <- suppressWarnings(cor.test(expr_119600["IL12RB2", ], xcell_scores[cell_type, ], method = "spearman", exact = FALSE))
    data.frame(cell_type = cell_type, rho = unname(test$estimate), p = test$p.value)
  }))
  xcell_correlations$BH_FDR <- p.adjust(xcell_correlations$p, method = "BH")
  write_table(xcell_correlations, "09_IL12RB2_xcell_correlations.tsv")

  desired <- c("NK cells", "NKT", "CD4+ T-cells", "CD8+ T-cells", "B-cells", "Monocytes", "ImmuneScore")
  available <- intersect(desired, rownames(xcell_scores))
  adjusted_data <- data.frame(
    IL12RB2 = as.numeric(expr_119600["IL12RB2", ]),
    group = factor(group_119600, levels = c("Control", "PBC"))
  )
  for (cell_type in available) {
    adjusted_data[[make.names(cell_type)]] <- as.numeric(scale(xcell_scores[cell_type, ]))
  }
  formula_text <- paste("IL12RB2 ~ group", paste(make.names(available), collapse = " + "), sep = if (length(available)) " + " else "")
  adjusted_model <- lm(as.formula(formula_text), data = adjusted_data)
  adjusted_test <- coeftest(adjusted_model, vcov. = vcovHC(adjusted_model, type = "HC3"))
  adjusted_table <- data.frame(
    term = rownames(adjusted_test),
    estimate = adjusted_test[, 1],
    robust_se = adjusted_test[, 2],
    t_value = adjusted_test[, 3],
    p_value = adjusted_test[, 4],
    row.names = NULL,
    check.names = FALSE
  )
  write_table(adjusted_table, "10_IL12RB2_composition_adjusted_HC3.tsv")
}

message("GSE159676: independent liver replication")
expr_159676 <- read_expression(file.path(input_dir, "GSE159676_gene_log2expr.tsv.gz"))
meta_159676 <- read.delim(file.path(input_dir, "GSE159676_metadata.tsv"), check.names = FALSE)
stopifnot(identical(colnames(expr_159676), meta_159676$sample_id))
group_159676 <- ifelse(grepl("Primary_biliary", meta_159676$title), "PBC", "Control")
de_159676 <- limma_contrast(expr_159676, factor(group_159676, levels = c("Control", "PBC")), "PBC-Control")
write_table(de_159676, "11_GSE159676_limma_PBC_vs_control.tsv")
comparison_159676 <- compare_vector(
  expr_159676["IL12RB2", ], group_159676, "PBC", "Control", "GSE159676_liver", "PBC_vs_control"
)
target_159676 <- de_159676[de_159676$gene_symbol == "IL12RB2", ]
comparison_159676$limma_log2FC <- target_159676$logFC
comparison_159676$limma_t <- target_159676$t
comparison_159676$limma_p <- target_159676$P.Value
comparison_159676$limma_fdr <- target_159676$adj.P.Val
comparison_159676$limma_rank <- target_159676$rank_by_p
target_comparisons <- bind_rows_fill(target_comparisons, comparison_159676)

message("GSE79850: liver outcome association")
expr_79850 <- read_expression(file.path(input_dir, "GSE79850_gene_log2expr.tsv.gz"))
meta_79850 <- read.delim(file.path(input_dir, "GSE79850_metadata.tsv"), check.names = FALSE)
stopifnot(identical(colnames(expr_79850), meta_79850$sample_id))
group_79850 <- meta_79850$sample_type
case_control_79850 <- ifelse(group_79850 == "Control", "Control", "PBC")
de_79850_case <- limma_contrast(expr_79850, factor(case_control_79850, levels = c("Control", "PBC")), "PBC-Control")
risk_group_79850 <- gsub(" ", "_", group_79850[group_79850 != "Control"])
de_79850_risk <- limma_contrast(
  expr_79850[, group_79850 != "Control", drop = FALSE],
  factor(risk_group_79850, levels = c("Low_Risk", "High_Risk")),
  "High_Risk-Low_Risk"
)
write_table(de_79850_case, "12_GSE79850_limma_PBC_vs_control.tsv")
write_table(de_79850_risk, "13_GSE79850_limma_high_vs_low_risk.tsv")
comparison_79850_case <- compare_vector(
  expr_79850["IL12RB2", ], case_control_79850, "PBC", "Control", "GSE79850_liver_NanoString", "PBC_vs_control"
)
comparison_79850_risk <- compare_vector(
  expr_79850["IL12RB2", ], group_79850, "High Risk", "Low Risk", "GSE79850_liver_NanoString", "high_vs_low_risk"
)
for (index in seq_len(2)) {
  source <- if (index == 1) de_79850_case else de_79850_risk
  target <- source[source$gene_symbol == "IL12RB2", ]
  object_name <- if (index == 1) "comparison_79850_case" else "comparison_79850_risk"
  object <- get(object_name)
  object$limma_log2FC <- target$logFC
  object$limma_t <- target$t
  object$limma_p <- target$P.Value
  object$limma_fdr <- target$adj.P.Val
  object$limma_rank <- target$rank_by_p
  assign(object_name, object)
}
target_comparisons <- bind_rows_fill(target_comparisons, comparison_79850_case, comparison_79850_risk)
raw_79850 <- 2^expr_79850["IL12RB2", ] - 1
floor_table <- aggregate(
  list(at_floor = raw_79850 <= 1.000001),
  list(group = group_79850),
  function(value) c(n = length(value), at_floor = sum(value), fraction = mean(value))
)
floor_table <- do.call(data.frame, floor_table)
write_table(floor_table, "14_GSE79850_IL12RB2_floor_audit.tsv")

message("GSE138888: DNT and UDCA-response sensitivity")
expr_138888 <- read_expression(file.path(input_dir, "GSE138888_gene_log2expr.tsv.gz"))
meta_138888 <- read.delim(file.path(input_dir, "GSE138888_analysis_metadata.tsv"), check.names = FALSE)
stopifnot(identical(colnames(expr_138888), meta_138888$sample_name))
group_138888 <- meta_138888$group
comparisons_138888 <- rbind(
  compare_vector(expr_138888["IL12RB2", ], group_138888, "Treatment_naive", "Control", "GSE138888_DNT", "treatment_naive_vs_control"),
  compare_vector(expr_138888["IL12RB2", ], group_138888, "UDCA_nonresponder", "UDCA_responder", "GSE138888_DNT", "UDCA_nonresponder_vs_responder")
)
target_comparisons <- bind_rows_fill(target_comparisons, comparisons_138888)

message("GSE304352: liver disease-control sensitivity")
counts_304352 <- read_expression(file.path(input_dir, "GSE304352_gene_counts.tsv.gz"))
meta_304352 <- read.delim(file.path(input_dir, "GSE304352_analysis_metadata.tsv"), check.names = FALSE)
stopifnot(identical(colnames(counts_304352), meta_304352$sample_name))
group_304352 <- factor(meta_304352$group, levels = c("CHB", "AIH", "PBC", "PBC_AIH"))
dge <- DGEList(counts = round(counts_304352), group = group_304352)
keep <- filterByExpr(dge, group = group_304352)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
design_304352 <- model.matrix(~ 0 + group_304352)
colnames(design_304352) <- levels(group_304352)
dge <- estimateDisp(dge, design_304352, robust = TRUE)
fit_304352 <- glmQLFit(dge, design_304352, robust = TRUE)
for (contrast_name in c("PBC-CHB", "PBC-AIH")) {
  contrast <- makeContrasts(contrasts = contrast_name, levels = design_304352)
  result <- topTags(glmQLFTest(fit_304352, contrast = contrast), n = Inf, sort.by = "PValue")$table
  result$gene_symbol <- rownames(result)
  result$rank_by_p <- seq_len(nrow(result))
  write_table(result[, c("gene_symbol", setdiff(names(result), "gene_symbol"))], paste0("15_GSE304352_edgeR_", gsub("-", "_vs_", contrast_name), ".tsv"))
}
log_cpm_304352 <- cpm(calcNormFactors(DGEList(counts = round(counts_304352))), log = TRUE, prior.count = 2)
comparison_304352_chb <- compare_vector(
  log_cpm_304352["IL12RB2", ], as.character(group_304352), "PBC", "CHB", "GSE304352_liver_RNAseq", "PBC_vs_CHB"
)
comparison_304352_aih <- compare_vector(
  log_cpm_304352["IL12RB2", ], as.character(group_304352), "PBC", "AIH", "GSE304352_liver_RNAseq", "PBC_vs_AIH"
)
target_comparisons <- bind_rows_fill(target_comparisons, comparison_304352_chb, comparison_304352_aih)

write_table(target_comparisons, "16_IL12RB2_all_cohort_comparisons.tsv")

message("Figures")
plot_data <- rbind(
  data.frame(dataset = "GSE119600 whole blood", group = group_119600, value = expr_119600["IL12RB2", ]),
  data.frame(dataset = "GSE159676 liver", group = group_159676, value = expr_159676["IL12RB2", ]),
  data.frame(dataset = "GSE79850 liver NanoString", group = group_79850, value = expr_79850["IL12RB2", ]),
  data.frame(dataset = "GSE138888 DNT", group = group_138888, value = expr_138888["IL12RB2", ])
)
plot_data$dataset <- factor(plot_data$dataset, levels = unique(plot_data$dataset))
p_box <- ggplot(plot_data, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.12, size = 1.5, alpha = 0.7) +
  facet_wrap(~dataset, scales = "free", ncol = 2) +
  labs(x = NULL, y = "IL12RB2 expression (log2 scale)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(figure_dir, "01_IL12RB2_cross_cohort_boxplots.pdf"), p_box, width = 9, height = 7)
save_png(p_box, file.path(figure_dir, "01_IL12RB2_cross_cohort_boxplots.png"), width = 9, height = 7)

forest_data <- target_comparisons[is.finite(target_comparisons$hedges_g), ]
forest_data$label <- paste(forest_data$dataset, forest_data$contrast, sep = ": ")
forest_data$lower <- forest_data$hedges_g - 1.96 * forest_data$hedges_g_se
forest_data$upper <- forest_data$hedges_g + 1.96 * forest_data$hedges_g_se
p_forest <- ggplot(forest_data, aes(x = hedges_g, y = reorder(label, hedges_g))) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(size = 2.5, color = "#B2182B") +
  labs(x = "Hedges g (positive = higher in case group)", y = NULL) +
  theme_bw(base_size = 10)
ggsave(file.path(figure_dir, "02_IL12RB2_effect_forest.pdf"), p_forest, width = 9, height = 5.5)
save_png(p_forest, file.path(figure_dir, "02_IL12RB2_effect_forest.png"), width = 9, height = 5.5)

if (!is.null(xcell_scores)) {
  xcell_tests <- read.delim(file.path(result_dir, "08_GSE119600_xcell_group_tests.tsv"), check.names = FALSE)
  relevant <- grepl("NK|NKT|T-cell|T cell|Th1|Treg|B-cell|B cell|Mono|Macro|DC|Dendritic|Neutro", xcell_tests$gene_symbol, ignore.case = TRUE)
  xcell_plot <- xcell_tests[relevant, ]
  xcell_plot <- xcell_plot[order(xcell_plot$P.Value), ]
  xcell_plot <- head(xcell_plot, 25)
  xcell_plot$direction <- ifelse(xcell_plot$logFC > 0, "Higher in PBC", "Lower in PBC")
  p_xcell <- ggplot(xcell_plot, aes(x = logFC, y = reorder(gene_symbol, logFC), color = direction, size = -log10(P.Value))) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    geom_point() +
    scale_color_manual(values = c("Higher in PBC" = "#B2182B", "Lower in PBC" = "#2166AC")) +
    labs(x = "xCell score difference (PBC - control)", y = NULL, color = NULL, size = "-log10(P)") +
    theme_bw(base_size = 10)
  ggsave(file.path(figure_dir, "03_GSE119600_xcell_differences.pdf"), p_xcell, width = 8, height = 7)
  save_png(p_xcell, file.path(figure_dir, "03_GSE119600_xcell_differences.png"), width = 8, height = 7)
}

session <- capture.output(sessionInfo())
writeLines(session, file.path(result_dir, "99_sessionInfo.txt"))
message("Bulk validation complete")
