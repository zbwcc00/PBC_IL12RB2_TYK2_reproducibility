suppressPackageStartupMessages({
  library(data.table)
  library(susieR)
  library(coloc)
})

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_argument), winslash = "/")
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/")
work_root <- file.path(repo_root, "work", "genetics")
result_dir <- file.path(work_root, "discovery", "results")
pbc_input_dir <- file.path(work_root, "fine_mapping", "input")
pbc_ld_dir <- file.path(work_root, "fine_mapping", "ld")
pair_input_dir <- file.path(work_root, "signal_level", "input")
pair_ld_dir <- file.path(work_root, "signal_level", "ld")
fit_dir <- file.path(work_root, "signal_level", "fits")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

pbc_n_case <- 8021
pbc_n_control <- 16489
pbc_n <- pbc_n_case + pbc_n_control
pbc_case_fraction <- pbc_n_case / pbc_n
analysis_filter <- Sys.getenv("ANALYSIS_FILTER", unset = "")
run_pbc_full <- Sys.getenv("RUN_PBC_FULL", unset = "1") != "0"
run_pairs <- Sys.getenv("RUN_PAIRS", unset = "1") != "0"
pair_L <- as.integer(Sys.getenv("PAIR_L", unset = "10"))
output_suffix <- Sys.getenv("OUTPUT_SUFFIX", unset = "")

pair_result_path <- function(stem) {
  file.path(result_dir, paste0(stem, output_suffix, ".tsv"))
}

read_ld <- function(prefix) {
  vars_path <- paste0(prefix, ".unphased.vcor1.bin.vars")
  matrix_path <- paste0(prefix, ".unphased.vcor1.bin")
  variants <- scan(vars_path, what = character(), quiet = TRUE)
  variant_n <- length(variants)
  connection <- file(matrix_path, "rb")
  on.exit(close(connection))
  values <- readBin(connection, what = numeric(), n = variant_n * variant_n,
                    size = 4, endian = "little")
  if (length(values) != variant_n * variant_n) {
    stop("Incomplete LD matrix: ", prefix)
  }
  matrix <- matrix(values, nrow = variant_n, ncol = variant_n, byrow = TRUE)
  matrix <- (matrix + t(matrix)) / 2
  diag(matrix) <- 1
  dimnames(matrix) <- list(variants, variants)
  list(variants = variants, matrix = matrix)
}

read_statistics <- function(path) fread(path)

align_statistics_to_ld <- function(statistics, variants) {
  index <- match(variants, statistics$rsid)
  if (anyNA(index)) {
    stop("LD variants missing from statistics: ", sum(is.na(index)))
  }
  aligned <- statistics[index]
  if (!identical(as.character(aligned$rsid), variants)) {
    stop("Failed to align statistics to LD order")
  }
  aligned
}

run_diagnostic <- function(z, R, n) {
  eigen_decomposition <- attr(R, "eigen")
  if (is.null(eigen_decomposition)) {
    eigen_decomposition <- eigen(R, symmetric = TRUE)
  }
  minimum_eigenvalue <- min(eigen_decomposition$values)
  negative_eigenvalue_n <- sum(eigen_decomposition$values < -1e-6)
  attr(R, "eigen") <- eigen_decomposition
  s_hat <- susieR::estimate_s_rss(z, R, n = n, method = "null-mle")
  kriging <- susieR::kriging_rss(z, R, n = n, s = s_hat)
  conditional <- kriging$conditional_dist
  list(
    R = R,
    s_hat = s_hat,
    minimum_eigenvalue = minimum_eigenvalue,
    negative_eigenvalue_n = negative_eigenvalue_n,
    kriging_outlier_n = sum(conditional$logLR > 2 & abs(conditional$z) > 2),
    kriging_max_abs_z_difference = max(abs(conditional$z_std_diff)),
    kriging_max_logLR = max(conditional$logLR)
  )
}

run_susie <- function(dataset, L) {
  fit <- susieR::susie_rss(
    z = dataset$beta / sqrt(dataset$varbeta),
    R = dataset$LD,
    n = dataset$N,
    L = L,
    max_iter = 500,
    refine = FALSE,
    estimate_residual_variance = FALSE
  )
  getFromNamespace("annotate_susie", "coloc")(fit, dataset$snp, dataset$LD)
}

credible_set_rows <- function(fit, metadata) {
  if (is.null(fit$sets$cs) || length(fit$sets$cs) == 0) {
    return(data.table())
  }
  rbindlist(lapply(seq_along(fit$sets$cs), function(index) {
    members <- fit$sets$cs[[index]]
    component <- fit$sets$cs_index[index]
    member_pip <- fit$pip[members]
    purity <- fit$sets$purity[index, , drop = FALSE]
    cbind(
      as.data.table(metadata),
      data.table(
      susie_component = component,
      credible_set_label = names(fit$sets$cs)[index],
      credible_set_size = length(members),
      credible_set_coverage = fit$sets$coverage[index],
      minimum_abs_correlation = purity$min.abs.corr,
      mean_abs_correlation = purity$mean.abs.corr,
      median_abs_correlation = purity$median.abs.corr,
      top_snp = names(member_pip)[which.max(member_pip)],
      top_snp_pip = max(member_pip),
      credible_set_snps = paste(names(member_pip), collapse = ";")
      )
    )
  }), fill = TRUE)
}

fit_summary <- function(fit) {
  credible_set_n <- if (is.null(fit$sets$cs)) 0L else length(fit$sets$cs)
  top_index <- which.max(fit$pip)
  list(
    converged = isTRUE(fit$converged),
    credible_set_n = credible_set_n,
    top_snp = names(fit$pip)[top_index],
    top_pip = fit$pip[top_index]
  )
}

get_component_cs <- function(fit, component) {
  index <- match(component, fit$sets$cs_index)
  if (is.na(index)) integer() else fit$sets$cs[[index]]
}

classification <- function(pp_h3, pp_h4) {
  shared_fraction <- pp_h4 / (pp_h3 + pp_h4)
  if (is.na(pp_h4)) return("not_estimable")
  if (pp_h4 >= 0.8 && shared_fraction >= 0.9) return("strong")
  if (pp_h4 >= 0.5) return("candidate")
  if (pp_h4 >= 0.2) return("exploratory")
  "no_support"
}

pbc_summary_rows <- list()
pbc_cs_rows <- list()
if (run_pbc_full) {
  pbc_audit <- fread(file.path(pbc_input_dir, "pbc_locus_allele_alignment_audit.tsv"))
  for (row_index in seq_len(nrow(pbc_audit))) {
    record <- pbc_audit[row_index]
    locus_id <- record$locus_id
    message("PBC_LOCUS_START\t", locus_id, "\t", record$pbc_locus)
    statistics <- read_statistics(record$data_path)
    ld <- read_ld(file.path(pbc_ld_dir, locus_id))
    statistics <- align_statistics_to_ld(statistics, ld$variants)
    diagnostic <- run_diagnostic(statistics$pbc_z_alt, ld$matrix, pbc_n)
    dataset <- list(
      beta = statistics$pbc_beta_alt,
      varbeta = statistics$pbc_se^2,
      snp = statistics$rsid,
      position = statistics$position_grch37_1based,
      MAF = pmin(statistics$reference_alt_freq, 1 - statistics$reference_alt_freq),
      LD = diagnostic$R,
      N = pbc_n,
      type = "cc",
      s = pbc_case_fraction
    )
    for (L in c(3L, 10L)) {
      start_time <- proc.time()[["elapsed"]]
      fit <- run_susie(dataset, L = L)
      elapsed_seconds <- proc.time()[["elapsed"]] - start_time
      summary <- fit_summary(fit)
      pbc_summary_rows[[length(pbc_summary_rows) + 1L]] <- data.table(
        locus_id = locus_id,
        pbc_locus = record$pbc_locus,
        lead_rsid = record$lead_rsid,
        variant_n = nrow(statistics),
        model_L = L,
        converged = summary$converged,
        credible_set_n = summary$credible_set_n,
        top_snp = summary$top_snp,
        top_pip = summary$top_pip,
        rss_s_hat = diagnostic$s_hat,
        ld_minimum_eigenvalue = diagnostic$minimum_eigenvalue,
        ld_negative_eigenvalue_n = diagnostic$negative_eigenvalue_n,
        kriging_outlier_n = diagnostic$kriging_outlier_n,
        kriging_max_abs_z_difference = diagnostic$kriging_max_abs_z_difference,
        kriging_max_logLR = diagnostic$kriging_max_logLR,
        elapsed_seconds = elapsed_seconds
      )
      cs <- credible_set_rows(
        fit,
        list(locus_id = locus_id, pbc_locus = record$pbc_locus, model_L = L)
      )
      if (nrow(cs)) pbc_cs_rows[[length(pbc_cs_rows) + 1L]] <- cs
      saveRDS(fit, file.path(fit_dir, paste0(locus_id, "__PBC__L", L, ".rds")), compress = "xz")
      message("PBC_LOCUS_DONE\t", locus_id, "\tL", L, "\tCS=", summary$credible_set_n)
    }
    fwrite(rbindlist(pbc_summary_rows, fill = TRUE),
           file.path(result_dir, "19_pbc_susie_multisignal_finemapping.tsv"), sep = "\t")
    fwrite(rbindlist(pbc_cs_rows, fill = TRUE),
           file.path(result_dir, "19b_pbc_susie_credible_set_details.tsv"), sep = "\t")
    rm(statistics, ld, diagnostic, dataset)
    invisible(gc())
  }
  fwrite(rbindlist(pbc_summary_rows, fill = TRUE),
         file.path(result_dir, "19_pbc_susie_multisignal_finemapping.tsv"), sep = "\t")
  fwrite(rbindlist(pbc_cs_rows, fill = TRUE),
         file.path(result_dir, "19b_pbc_susie_credible_set_details.tsv"), sep = "\t")
}

if (!run_pairs) {
  message("ALL_COMPLETE\tPBC_only")
  quit(save = "no", status = 0)
}

pair_audit <- fread(file.path(pair_input_dir, "signal_aware_allele_alignment_audit.tsv"))
if (nzchar(analysis_filter)) {
  pair_audit <- pair_audit[analysis_id %in% strsplit(analysis_filter, ",", fixed = TRUE)[[1]]]
}
if (nrow(pair_audit) == 0) stop("No signal-aware analyses selected")

pair_summary_rows <- list()
primary_coloc_rows <- list()
sensitivity_coloc_rows <- list()
pair_cs_rows <- list()
priors <- c(1e-5, 5e-6, 1e-6)

for (row_index in seq_len(nrow(pair_audit))) {
  record <- pair_audit[row_index]
  analysis_id <- record$analysis_id
  message("PAIR_START\t", analysis_id, "\t", record$gene_symbol, "\t", record$sample_group)
  statistics <- read_statistics(record$data_path)
  ld <- read_ld(file.path(pair_ld_dir, analysis_id))
  statistics <- align_statistics_to_ld(statistics, ld$variants)

  pbc_diagnostic <- run_diagnostic(statistics$pbc_z_alt, ld$matrix, pbc_n)
  eqtl_diagnostic <- run_diagnostic(statistics$eqtl_z_alt, pbc_diagnostic$R, record$sample_size)
  R <- eqtl_diagnostic$R

  pbc_dataset <- list(
    beta = statistics$pbc_beta_alt,
    varbeta = statistics$pbc_se^2,
    snp = statistics$rsid,
    position = statistics$position_grch37_1based,
    MAF = pmin(statistics$reference_alt_freq, 1 - statistics$reference_alt_freq),
    LD = R,
    N = pbc_n,
    type = "cc",
    s = pbc_case_fraction
  )
  eqtl_dataset <- list(
    beta = statistics$eqtl_beta_alt,
    varbeta = statistics$eqtl_se^2,
    snp = statistics$rsid,
    position = statistics$eqtl_position_grch38_1based,
    MAF = statistics$eqtl_maf,
    LD = R,
    N = record$sample_size,
    type = "quant"
  )

  pbc_fit <- run_susie(pbc_dataset, L = pair_L)
  eqtl_fit <- run_susie(eqtl_dataset, L = pair_L)
  pbc_fit_summary <- fit_summary(pbc_fit)
  eqtl_fit_summary <- fit_summary(eqtl_fit)
  pair_summary_rows[[length(pair_summary_rows) + 1L]] <- data.table(
    analysis_id = analysis_id,
    pbc_locus = record$pbc_locus,
    gene_symbol = record$gene_symbol,
    gene_id = record$gene_id,
    sample_group = record$sample_group,
    sample_size = record$sample_size,
    variant_n = nrow(statistics),
    pbc_converged = pbc_fit_summary$converged,
    pbc_credible_set_n = pbc_fit_summary$credible_set_n,
    pbc_top_snp = pbc_fit_summary$top_snp,
    pbc_top_pip = pbc_fit_summary$top_pip,
    pbc_rss_s_hat = pbc_diagnostic$s_hat,
    pbc_kriging_outlier_n = pbc_diagnostic$kriging_outlier_n,
    eqtl_converged = eqtl_fit_summary$converged,
    eqtl_credible_set_n = eqtl_fit_summary$credible_set_n,
    eqtl_top_snp = eqtl_fit_summary$top_snp,
    eqtl_top_pip = eqtl_fit_summary$top_pip,
    eqtl_rss_s_hat = eqtl_diagnostic$s_hat,
    eqtl_kriging_outlier_n = eqtl_diagnostic$kriging_outlier_n,
    ld_minimum_eigenvalue = min(pbc_diagnostic$minimum_eigenvalue, eqtl_diagnostic$minimum_eigenvalue),
    ld_negative_eigenvalue_n = max(pbc_diagnostic$negative_eigenvalue_n, eqtl_diagnostic$negative_eigenvalue_n)
  )

  pbc_cs <- credible_set_rows(
    pbc_fit,
    list(analysis_id = analysis_id, trait = "PBC", gene_symbol = record$gene_symbol,
         sample_group = record$sample_group)
  )
  eqtl_cs <- credible_set_rows(
    eqtl_fit,
    list(analysis_id = analysis_id, trait = "eQTL", gene_symbol = record$gene_symbol,
         sample_group = record$sample_group)
  )
  if (nrow(pbc_cs)) pair_cs_rows[[length(pair_cs_rows) + 1L]] <- pbc_cs
  if (nrow(eqtl_cs)) pair_cs_rows[[length(pair_cs_rows) + 1L]] <- eqtl_cs

  saveRDS(list(pbc = pbc_fit, eqtl = eqtl_fit),
          file.path(fit_dir, paste0(analysis_id, "__", record$gene_symbol, "__", record$sample_group,
                                    "__L", pair_L, output_suffix, ".rds")),
          compress = "xz")

  for (p12 in priors) {
    coloc_result <- suppressWarnings(coloc::coloc.susie(
      pbc_fit,
      eqtl_fit,
      p1 = 1e-4,
      p2 = 1e-4,
      p12 = p12
    ))
    if (!is.list(coloc_result) || is.null(coloc_result$summary) || nrow(coloc_result$summary) == 0) {
      sensitivity_coloc_rows[[length(sensitivity_coloc_rows) + 1L]] <- data.table(
        analysis_id = analysis_id,
        pbc_locus = record$pbc_locus,
        gene_symbol = record$gene_symbol,
        sample_group = record$sample_group,
        p12 = p12,
        status = "no_credible_set_pair"
      )
      next
    }
    summary_table <- as.data.table(copy(coloc_result$summary))
    posterior_table <- as.data.table(coloc_result$results)
    for (signal_row in seq_len(nrow(summary_table))) {
      result <- summary_table[signal_row]
      posterior_column <- names(posterior_table)[signal_row + 1L]
      top_posterior_index <- which.max(posterior_table[[posterior_column]])
      top_shared_snp <- posterior_table$snp[top_posterior_index]
      top_shared_snp_posterior <- posterior_table[[posterior_column]][top_posterior_index]
      pbc_members <- get_component_cs(pbc_fit, result$idx1)
      eqtl_members <- get_component_cs(eqtl_fit, result$idx2)
      overlap_n <- length(intersect(pbc_members, eqtl_members))
      union_n <- length(union(pbc_members, eqtl_members))
      output_row <- data.table(
        analysis_id = analysis_id,
        pbc_locus = record$pbc_locus,
        gene_symbol = record$gene_symbol,
        gene_id = record$gene_id,
        sample_group = record$sample_group,
        sample_size = record$sample_size,
        variant_n = nrow(statistics),
        pbc_signal = result$idx1,
        eqtl_signal = result$idx2,
        p12 = p12,
        PP.H0 = result$PP.H0.abf,
        PP.H1 = result$PP.H1.abf,
        PP.H2 = result$PP.H2.abf,
        PP.H3 = result$PP.H3.abf,
        PP.H4 = result$PP.H4.abf,
        shared_fraction_H4_over_H3_H4 = result$PP.H4.abf / (result$PP.H3.abf + result$PP.H4.abf),
        top_pbc_snp = result$hit1,
        top_eqtl_snp = result$hit2,
        top_shared_snp = top_shared_snp,
        top_shared_snp_posterior = top_shared_snp_posterior,
        pbc_credible_set_size = length(pbc_members),
        eqtl_credible_set_size = length(eqtl_members),
        credible_set_overlap_n = overlap_n,
        credible_set_jaccard = if (union_n == 0) NA_real_ else overlap_n / union_n,
        evidence_class = classification(result$PP.H3.abf, result$PP.H4.abf),
        status = "estimated"
      )
      sensitivity_coloc_rows[[length(sensitivity_coloc_rows) + 1L]] <- output_row
      if (identical(p12, 1e-6)) {
        primary_coloc_rows[[length(primary_coloc_rows) + 1L]] <- output_row
      }
    }
  }
  message("PAIR_DONE\t", analysis_id, "\tPBC_CS=", pbc_fit_summary$credible_set_n,
          "\teQTL_CS=", eqtl_fit_summary$credible_set_n)
  fwrite(rbindlist(pair_summary_rows, fill = TRUE),
         pair_result_path("20_onek1k_susie_concordance"), sep = "\t")
  fwrite(rbindlist(pair_cs_rows, fill = TRUE),
         pair_result_path("20b_signal_aware_susie_credible_set_details"), sep = "\t")
  fwrite(rbindlist(primary_coloc_rows, fill = TRUE),
         pair_result_path("21_coloc_susie_signal_pairs"), sep = "\t")
  fwrite(rbindlist(sensitivity_coloc_rows, fill = TRUE),
         pair_result_path("22_coloc_susie_prior_sensitivity"), sep = "\t")
  rm(statistics, ld, pbc_diagnostic, eqtl_diagnostic, R, pbc_dataset, eqtl_dataset,
     pbc_fit, eqtl_fit)
  invisible(gc())
}

fwrite(rbindlist(pair_summary_rows, fill = TRUE),
       pair_result_path("20_onek1k_susie_concordance"), sep = "\t")
fwrite(rbindlist(pair_cs_rows, fill = TRUE),
       pair_result_path("20b_signal_aware_susie_credible_set_details"), sep = "\t")
fwrite(rbindlist(primary_coloc_rows, fill = TRUE),
       pair_result_path("21_coloc_susie_signal_pairs"), sep = "\t")
fwrite(rbindlist(sensitivity_coloc_rows, fill = TRUE),
       pair_result_path("22_coloc_susie_prior_sensitivity"), sep = "\t")

message("ALL_COMPLETE\tpairs=", nrow(pair_audit))
