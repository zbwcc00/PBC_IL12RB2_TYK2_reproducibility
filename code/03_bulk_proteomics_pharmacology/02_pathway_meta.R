options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(lmtest)
  library(sandwich)
})

set.seed(20261014)

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_argument), winslash = "/")
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/")
input_dir <- file.path(repo_root, "data", "raw", "bulk", "input")
prior_result_dir <- file.path(repo_root, "work", "bulk", "validation", "results")
output_dir <- file.path(repo_root, "work", "bulk", "pathway_meta", "results")
figure_dir <- file.path(repo_root, "work", "bulk", "pathway_meta", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(data, name) {
  write.table(data, file.path(output_dir, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

read_expression <- function(path) {
  data <- read.delim(gzfile(path), check.names = FALSE, quote = "", comment.char = "")
  genes <- make.unique(as.character(data[[1]]))
  expression <- as.matrix(data[, -1, drop = FALSE])
  storage.mode(expression) <- "double"
  rownames(expression) <- genes
  expression
}

signatures <- list(
  IL12_STAT4 = c("IL12RB1", "TYK2", "JAK2", "STAT4", "TBX21", "IFNG"),
  IFNG_RESPONSE = c("IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1", "IRF1", "GBP1", "GBP2", "TAP1", "CXCL9", "CXCL10")
)

trimmed_signatures <- list(
  IL12_STAT4_no_IFNG = setdiff(signatures$IL12_STAT4, "IFNG"),
  IFNG_RESPONSE_no_CXCL9_CXCL10 = setdiff(signatures$IFNG_RESPONSE, c("CXCL9", "CXCL10"))
)

signature_manifest <- do.call(rbind, lapply(names(signatures), function(signature_name) {
  data.frame(
    signature = signature_name,
    gene_symbol = signatures[[signature_name]],
    role = "primary_prespecified_non_circular",
    IL12RB2_included = FALSE,
    stringsAsFactors = FALSE
  )
}))
write_tsv(signature_manifest, "00_signature_manifest.tsv")

robust_standardize <- function(values) {
  center <- median(values, na.rm = TRUE)
  spread <- mad(values, center = center, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(spread) || spread <= 0) {
    spread <- sd(values, na.rm = TRUE)
  }
  if (!is.finite(spread) || spread <= 0) {
    return(rep(0, length(values)))
  }
  standardized <- (values - center) / spread
  pmax(-3, pmin(3, standardized))
}

rank_normalize <- function(values) {
  n <- sum(is.finite(values))
  output <- rep(NA_real_, length(values))
  output[is.finite(values)] <- qnorm((rank(values[is.finite(values)], ties.method = "average") - 0.5) / n)
  output
}

score_signature <- function(expression, genes, minimum_genes = 4L) {
  covered <- intersect(genes, rownames(expression))
  if (length(covered) < minimum_genes) {
    return(list(
      robust_z = rep(NA_real_, ncol(expression)),
      rank_z = rep(NA_real_, ncol(expression)),
      covered = covered,
      eligible = FALSE
    ))
  }
  robust_matrix <- t(apply(expression[covered, , drop = FALSE], 1, robust_standardize))
  rank_matrix <- t(apply(expression[covered, , drop = FALSE], 1, rank_normalize))
  list(
    robust_z = colMeans(robust_matrix, na.rm = TRUE),
    rank_z = colMeans(rank_matrix, na.rm = TRUE),
    covered = covered,
    eligible = TRUE
  )
}

hedges_g <- function(case, control) {
  case <- case[is.finite(case)]
  control <- control[is.finite(control)]
  n_case <- length(case)
  n_control <- length(control)
  pooled_sd <- sqrt(((n_case - 1) * var(case) + (n_control - 1) * var(control)) /
    (n_case + n_control - 2))
  if (!is.finite(pooled_sd) || pooled_sd <= 0) {
    return(c(g = NA_real_, variance = NA_real_, se = NA_real_))
  }
  correction <- 1 - 3 / (4 * (n_case + n_control) - 9)
  effect <- correction * (mean(case) - mean(control)) / pooled_sd
  variance <- (n_case + n_control) / (n_case * n_control) + effect^2 / (2 * (n_case + n_control - 2))
  c(g = effect, variance = variance, se = sqrt(variance))
}

permutation_test <- function(values, labels, case_label, control_label, max_exact = 100000L, monte_carlo = 100000L) {
  keep <- labels %in% c(case_label, control_label) & is.finite(values)
  values <- values[keep]
  labels <- labels[keep]
  n_case <- sum(labels == case_label)
  n_total <- length(values)
  observed <- mean(values[labels == case_label]) - mean(values[labels == control_label])
  combinations <- choose(n_total, n_case)
  if (is.finite(combinations) && combinations <= max_exact) {
    indices <- combn(n_total, n_case)
    case_sums <- colSums(matrix(values[indices], nrow = n_case))
    null <- case_sums / n_case - (sum(values) - case_sums) / (n_total - n_case)
    method <- "exact_enumeration"
  } else {
    null <- replicate(monte_carlo, {
      indices <- sample.int(n_total, n_case, replace = FALSE)
      mean(values[indices]) - mean(values[-indices])
    })
    method <- paste0("monte_carlo_", monte_carlo)
  }
  c(
    mean_difference = observed,
    p_two_sided = (sum(abs(null) >= abs(observed)) + 1) / (length(null) + 1),
    p_case_greater = (sum(null >= observed) + 1) / (length(null) + 1),
    n_assignments = length(null),
    method = method
  )
}

compare_score <- function(values, labels, case_label, control_label, dataset, tissue, contrast, signature, score_method, run_permutation = TRUE) {
  keep <- labels %in% c(case_label, control_label) & is.finite(values)
  values <- values[keep]
  labels <- labels[keep]
  case <- values[labels == case_label]
  control <- values[labels == control_label]
  effect <- hedges_g(case, control)
  permutation <- if (run_permutation) {
    permutation_test(values, labels, case_label, control_label)
  } else {
    c(mean_difference = mean(case) - mean(control), p_two_sided = NA, p_case_greater = NA, n_assignments = NA, method = "not_run")
  }
  data.frame(
    dataset = dataset,
    tissue = tissue,
    contrast = contrast,
    signature = signature,
    score_method = score_method,
    case_label = case_label,
    control_label = control_label,
    n_case = length(case),
    n_control = length(control),
    mean_case = mean(case),
    mean_control = mean(control),
    mean_difference = as.numeric(permutation[["mean_difference"]]),
    hedges_g = effect[["g"]],
    variance = effect[["variance"]],
    se = effect[["se"]],
    wilcoxon_p = suppressWarnings(wilcox.test(case, control, exact = FALSE)$p.value),
    permutation_p_two_sided = as.numeric(permutation[["p_two_sided"]]),
    permutation_p_case_greater = as.numeric(permutation[["p_case_greater"]]),
    permutation_method = as.character(permutation[["method"]]),
    stringsAsFactors = FALSE
  )
}

random_effect_mkh <- function(effect, variance, model_name) {
  keep <- is.finite(effect) & is.finite(variance) & variance > 0
  effect <- effect[keep]
  variance <- variance[keep]
  k <- length(effect)
  if (k < 2) {
    return(data.frame(model = model_name, k = k, estimate = NA, se_mKH = NA, ci_low_mKH = NA, ci_high_mKH = NA, p_mKH = NA, tau2_DL = NA, Q = NA, Q_p = NA, I2_percent = NA))
  }
  fixed_weights <- 1 / variance
  fixed_effect <- sum(fixed_weights * effect) / sum(fixed_weights)
  q <- sum(fixed_weights * (effect - fixed_effect)^2)
  df <- k - 1
  c_term <- sum(fixed_weights) - sum(fixed_weights^2) / sum(fixed_weights)
  tau2 <- max(0, (q - df) / c_term)
  random_weights <- 1 / (variance + tau2)
  pooled <- sum(random_weights * effect) / sum(random_weights)
  hk_q <- sum(random_weights * (effect - pooled)^2) / df
  modified_hk_q <- max(1, hk_q)
  se <- sqrt(modified_hk_q / sum(random_weights))
  critical <- qt(0.975, df = df)
  statistic <- pooled / se
  data.frame(
    model = model_name,
    k = k,
    estimate = pooled,
    se_mKH = se,
    ci_low_mKH = pooled - critical * se,
    ci_high_mKH = pooled + critical * se,
    p_mKH = 2 * pt(-abs(statistic), df = df),
    tau2_DL = tau2,
    Q = q,
    Q_p = pchisq(q, df = df, lower.tail = FALSE),
    I2_percent = if (q > 0) max(0, (q - df) / q) * 100 else 0,
    positive_effects = sum(effect > 0),
    negative_effects = sum(effect < 0),
    stringsAsFactors = FALSE
  )
}

expr_119600 <- read_expression(file.path(input_dir, "GSE119600_gene_log2expr.tsv.gz"))
meta_119600 <- read.delim(file.path(input_dir, "GSE119600_metadata.tsv"), check.names = FALSE)
group_119600 <- ifelse(meta_119600$source_name_ch1 == "Primary biliary cholangitis, adult", "PBC", "Control")

expr_159676 <- read_expression(file.path(input_dir, "GSE159676_gene_log2expr.tsv.gz"))
meta_159676 <- read.delim(file.path(input_dir, "GSE159676_metadata.tsv"), check.names = FALSE)
group_159676 <- ifelse(grepl("Primary_biliary", meta_159676$title), "PBC", "Control")

expr_79850 <- read_expression(file.path(input_dir, "GSE79850_gene_log2expr.tsv.gz"))
meta_79850 <- read.delim(file.path(input_dir, "GSE79850_metadata.tsv"), check.names = FALSE)
group_79850 <- meta_79850$sample_type
case_control_79850 <- ifelse(group_79850 == "Control", "Control", "PBC")

expr_138888 <- read_expression(file.path(input_dir, "GSE138888_gene_log2expr.tsv.gz"))
meta_138888 <- read.delim(file.path(input_dir, "GSE138888_analysis_metadata.tsv"), check.names = FALSE)
group_138888 <- meta_138888$group

counts_304352 <- read_expression(file.path(input_dir, "GSE304352_gene_counts.tsv.gz"))
meta_304352 <- read.delim(file.path(input_dir, "GSE304352_analysis_metadata.tsv"), check.names = FALSE)
group_304352 <- meta_304352$group
expr_304352 <- cpm(calcNormFactors(DGEList(counts = round(counts_304352))), log = TRUE, prior.count = 2)

stopifnot(
  identical(colnames(expr_119600), meta_119600$sample_id),
  identical(colnames(expr_159676), meta_159676$sample_id),
  identical(colnames(expr_79850), meta_79850$sample_id),
  identical(colnames(expr_138888), meta_138888$sample_name),
  identical(colnames(expr_304352), meta_304352$sample_name)
)

datasets <- list(
  GSE119600_whole_blood = list(expression = expr_119600, labels = group_119600, tissue = "whole_blood", case = "PBC", control = "Control", contrast = "PBC_vs_healthy_control", primary = TRUE),
  GSE159676_liver = list(expression = expr_159676, labels = group_159676, tissue = "liver", case = "PBC", control = "Control", contrast = "PBC_vs_healthy_control", primary = TRUE),
  GSE79850_liver_NanoString = list(expression = expr_79850, labels = case_control_79850, tissue = "liver", case = "PBC", control = "Control", contrast = "PBC_vs_healthy_control", primary = TRUE),
  GSE138888_DNT = list(expression = expr_138888, labels = group_138888, tissue = "double_negative_T_cells", case = "Treatment_naive", control = "Control", contrast = "treatment_naive_PBC_vs_healthy_control", primary = TRUE),
  GSE304352_liver_PBC_vs_CHB = list(expression = expr_304352, labels = group_304352, tissue = "liver", case = "PBC", control = "CHB", contrast = "PBC_vs_CHB", primary = FALSE),
  GSE304352_liver_PBC_vs_AIH = list(expression = expr_304352, labels = group_304352, tissue = "liver", case = "PBC", control = "AIH", contrast = "PBC_vs_AIH", primary = FALSE)
)

coverage_rows <- list()
score_rows <- list()
effect_rows <- list()
score_cache <- list()

for (dataset_name in names(datasets)) {
  dataset <- datasets[[dataset_name]]
  score_cache[[dataset_name]] <- list()
  for (signature_name in names(signatures)) {
    scored <- score_signature(dataset$expression, signatures[[signature_name]])
    coverage_rows[[length(coverage_rows) + 1L]] <- data.frame(
      dataset = dataset_name,
      tissue = dataset$tissue,
      signature = signature_name,
      n_prespecified = length(signatures[[signature_name]]),
      n_covered = length(scored$covered),
      coverage_fraction = length(scored$covered) / length(signatures[[signature_name]]),
      covered_genes = paste(scored$covered, collapse = ";"),
      missing_genes = paste(setdiff(signatures[[signature_name]], scored$covered), collapse = ";"),
      eligible_minimum_four_genes = scored$eligible,
      stringsAsFactors = FALSE
    )
    score_cache[[dataset_name]][[signature_name]] <- scored
    if (!scored$eligible) {
      next
    }
    score_rows[[length(score_rows) + 1L]] <- data.frame(
      dataset = dataset_name,
      sample_id = colnames(dataset$expression),
      group = dataset$labels,
      signature = signature_name,
      robust_z = scored$robust_z,
      rank_z = scored$rank_z,
      stringsAsFactors = FALSE
    )
    for (score_method in c("robust_z", "rank_z")) {
      effect_rows[[length(effect_rows) + 1L]] <- compare_score(
        scored[[score_method]], dataset$labels, dataset$case, dataset$control,
        dataset_name, dataset$tissue, dataset$contrast, signature_name, score_method
      )
    }
  }
}

coverage_table <- do.call(rbind, coverage_rows)
score_table <- do.call(rbind, score_rows)
effect_table <- do.call(rbind, effect_rows)
write_tsv(coverage_table, "01_signature_coverage_by_dataset.tsv")
write_tsv(score_table, "02_signature_scores_by_sample.tsv")
write_tsv(effect_table, "03_all_cohort_effects.tsv")

primary_names <- names(datasets)[vapply(datasets, function(dataset) isTRUE(dataset$primary), logical(1))]
primary_effects <- effect_table[effect_table$dataset %in% primary_names, , drop = FALSE]
meta_rows <- list()
for (signature_name in names(signatures)) {
  for (score_method in c("robust_z", "rank_z")) {
    subset <- primary_effects[primary_effects$signature == signature_name & primary_effects$score_method == score_method, ]
    meta <- random_effect_mkh(subset$hedges_g, subset$variance, "random_effect_modified_Hartung_Knapp")
    meta$signature <- signature_name
    meta$score_method <- score_method
    meta$datasets <- paste(subset$dataset, collapse = ";")
    meta_rows[[length(meta_rows) + 1L]] <- meta
  }
}
meta_table <- do.call(rbind, meta_rows)
meta_table$BH_FDR_four_primary_tests <- p.adjust(meta_table$p_mKH, method = "BH")
write_tsv(meta_table, "04_primary_random_effect_meta.tsv")

lodo_rows <- list()
for (signature_name in names(signatures)) {
  subset <- primary_effects[primary_effects$signature == signature_name & primary_effects$score_method == "robust_z", ]
  for (left_out in subset$dataset) {
    retained <- subset[subset$dataset != left_out, ]
    meta <- random_effect_mkh(retained$hedges_g, retained$variance, "leave_one_dataset_out_mKH")
    meta$signature <- signature_name
    meta$left_out_dataset <- left_out
    lodo_rows[[length(lodo_rows) + 1L]] <- meta
  }
}
lodo_table <- do.call(rbind, lodo_rows)
write_tsv(lodo_table, "05_leave_one_dataset_out.tsv")

loo_gene_rows <- list()
for (signature_name in names(signatures)) {
  for (left_out_gene in signatures[[signature_name]]) {
    reduced_genes <- setdiff(signatures[[signature_name]], left_out_gene)
    effects <- list()
    for (dataset_name in primary_names) {
      dataset <- datasets[[dataset_name]]
      scored <- score_signature(dataset$expression, reduced_genes, minimum_genes = 4L)
      if (!scored$eligible) {
        next
      }
      effects[[length(effects) + 1L]] <- compare_score(
        scored$robust_z, dataset$labels, dataset$case, dataset$control,
        dataset_name, dataset$tissue, dataset$contrast, signature_name, "robust_z", run_permutation = FALSE
      )
    }
    effects <- do.call(rbind, effects)
    meta <- random_effect_mkh(effects$hedges_g, effects$variance, "leave_one_gene_out_mKH")
    meta$signature <- signature_name
    meta$left_out_gene <- left_out_gene
    meta$covered_in_k_datasets <- nrow(effects)
    loo_gene_rows[[length(loo_gene_rows) + 1L]] <- meta
  }
}
loo_gene_table <- do.call(rbind, loo_gene_rows)
write_tsv(loo_gene_table, "06_leave_one_gene_out.tsv")

trimmed_rows <- list()
for (signature_name in names(trimmed_signatures)) {
  effects <- list()
  for (dataset_name in primary_names) {
    dataset <- datasets[[dataset_name]]
    scored <- score_signature(dataset$expression, trimmed_signatures[[signature_name]], minimum_genes = 4L)
    if (!scored$eligible) {
      next
    }
    effects[[length(effects) + 1L]] <- compare_score(
      scored$robust_z, dataset$labels, dataset$case, dataset$control,
      dataset_name, dataset$tissue, dataset$contrast, signature_name, "robust_z", run_permutation = FALSE
    )
  }
  effects <- do.call(rbind, effects)
  meta <- random_effect_mkh(effects$hedges_g, effects$variance, "driver_trimmed_mKH")
  meta$signature <- signature_name
  meta$genes <- paste(trimmed_signatures[[signature_name]], collapse = ";")
  trimmed_rows[[length(trimmed_rows) + 1L]] <- meta
}
trimmed_table <- do.call(rbind, trimmed_rows)
write_tsv(trimmed_table, "07_inflammatory_driver_trimmed_meta.tsv")

clinical_rows <- list()
for (signature_name in names(signatures)) {
  score_79850 <- score_cache$GSE79850_liver_NanoString[[signature_name]]
  score_138888 <- score_cache$GSE138888_DNT[[signature_name]]
  if (score_79850$eligible) {
    clinical_rows[[length(clinical_rows) + 1L]] <- compare_score(
      score_79850$robust_z, group_79850, "High Risk", "Low Risk",
      "GSE79850_liver_NanoString", "liver", "high_vs_low_risk", signature_name, "robust_z"
    )
  }
  if (score_138888$eligible) {
    clinical_rows[[length(clinical_rows) + 1L]] <- compare_score(
      score_138888$robust_z, group_138888, "UDCA_nonresponder", "UDCA_responder",
      "GSE138888_DNT", "double_negative_T_cells", "UDCA_nonresponder_vs_responder", signature_name, "robust_z"
    )
  }
}
clinical_table <- do.call(rbind, clinical_rows)
clinical_table$BH_FDR_four_exploratory_tests <- p.adjust(clinical_table$permutation_p_two_sided, method = "BH")
write_tsv(clinical_table, "08_clinical_outcome_exploratory_effects.tsv")

composition_rows <- list()
xcell_path <- file.path(prior_result_dir, "07_GSE119600_xcell_scores.tsv")
if (file.exists(xcell_path)) {
  xcell_data <- read.delim(xcell_path, check.names = FALSE)
  rownames(xcell_data) <- xcell_data$cell_type
  xcell_matrix <- as.matrix(xcell_data[, -1, drop = FALSE])
  storage.mode(xcell_matrix) <- "double"
  desired <- c("NK cells", "NKT", "CD4+ T-cells", "CD8+ T-cells", "B-cells", "Monocytes", "ImmuneScore")
  available <- intersect(desired, rownames(xcell_matrix))
  for (signature_name in names(signatures)) {
    scored <- score_cache$GSE119600_whole_blood[[signature_name]]
    model_data <- data.frame(
      score = scored$robust_z,
      group = factor(group_119600, levels = c("Control", "PBC"))
    )
    for (cell_type in available) {
      model_data[[make.names(cell_type)]] <- as.numeric(scale(xcell_matrix[cell_type, colnames(expr_119600)]))
    }
    formula_text <- paste("score ~ group", paste(make.names(available), collapse = " + "), sep = if (length(available)) " + " else "")
    model <- lm(as.formula(formula_text), data = model_data)
    test <- coeftest(model, vcov. = vcovHC(model, type = "HC3"))
    term <- "groupPBC"
    composition_rows[[length(composition_rows) + 1L]] <- data.frame(
      signature = signature_name,
      term = term,
      estimate = test[term, 1],
      robust_se_HC3 = test[term, 2],
      t_value = test[term, 3],
      p_value = test[term, 4],
      adjusted_cell_scores = paste(available, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
}
composition_table <- do.call(rbind, composition_rows)
if (!is.null(composition_table)) {
  composition_table$BH_FDR_two_signatures <- p.adjust(composition_table$p_value, method = "BH")
  write_tsv(composition_table, "09_GSE119600_composition_adjusted_pathways.tsv")
}

primary_plot <- primary_effects[primary_effects$score_method == "robust_z", ]
primary_plot$label <- paste(primary_plot$dataset, primary_plot$tissue, sep = " | ")
primary_plot$ci_low <- primary_plot$hedges_g - qnorm(0.975) * primary_plot$se
primary_plot$ci_high <- primary_plot$hedges_g + qnorm(0.975) * primary_plot$se
meta_plot <- meta_table[meta_table$score_method == "robust_z", ]
meta_plot$dataset <- "Random-effects mHK"
meta_plot$tissue <- "pooled"
meta_plot$label <- paste(meta_plot$dataset, meta_plot$tissue, sep = " | ")
meta_plot$hedges_g <- meta_plot$estimate
meta_plot$ci_low <- meta_plot$ci_low_mKH
meta_plot$ci_high <- meta_plot$ci_high_mKH
plot_data <- rbind(
  primary_plot[, c("signature", "label", "hedges_g", "ci_low", "ci_high")],
  meta_plot[, c("signature", "label", "hedges_g", "ci_low", "ci_high")]
)
plot_data$is_pooled <- plot_data$label == "Random-effects mHK | pooled"
p_forest <- ggplot(plot_data, aes(x = hedges_g, y = reorder(label, hedges_g), color = is_pooled)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18) +
  geom_point(aes(size = is_pooled)) +
  facet_wrap(~signature, scales = "free_y") +
  scale_color_manual(values = c(`FALSE` = "#2166AC", `TRUE` = "#B2182B"), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 2.2, `TRUE` = 3.5), guide = "none") +
  labs(x = "Hedges g (positive = higher pathway score in PBC)", y = NULL) +
  theme_bw(base_size = 10)
ggsave(file.path(figure_dir, "01_pathway_random_effect_forest.pdf"), p_forest, width = 10, height = 6.5)
png(file.path(figure_dir, "01_pathway_random_effect_forest.png"), width = 3000, height = 1950, res = 300, type = "cairo-png")
print(p_forest)
dev.off()

decision_rows <- do.call(rbind, lapply(names(signatures), function(signature_name) {
  robust <- meta_table[meta_table$signature == signature_name & meta_table$score_method == "robust_z", ]
  rank <- meta_table[meta_table$signature == signature_name & meta_table$score_method == "rank_z", ]
  lodo <- lodo_table[lodo_table$signature == signature_name, ]
  loo_gene <- loo_gene_table[loo_gene_table$signature == signature_name, ]
  data.frame(
    signature = signature_name,
    robust_meta_effect = robust$estimate,
    robust_meta_p_mKH = robust$p_mKH,
    robust_meta_FDR = robust$BH_FDR_four_primary_tests,
    rank_meta_effect = rank$estimate,
    rank_meta_p_mKH = rank$p_mKH,
    same_direction_two_scoring_methods = sign(robust$estimate) == sign(rank$estimate),
    all_four_primary_effects_positive = robust$positive_effects == 4,
    all_leave_one_dataset_pooled_estimates_positive = all(lodo$estimate > 0, na.rm = TRUE),
    all_leave_one_gene_pooled_estimates_positive = all(loo_gene$estimate > 0, na.rm = TRUE),
    decision = ifelse(
      robust$estimate > 0 && rank$estimate > 0 && robust$positive_effects >= 3,
      "retain_as_cross_cohort_directional_pathway_support",
      "do_not_promote_stop_expansion"
    ),
    stringsAsFactors = FALSE
  )
}))
write_tsv(decision_rows, "10_prespecified_decision_summary.tsv")

capture.output(sessionInfo(), file = file.path(output_dir, "99_sessionInfo.txt"))
print(meta_table, row.names = FALSE)
print(decision_rows, row.names = FALSE)
