suppressPackageStartupMessages(library(data.table))

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_argument), winslash = "/")
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/")
work_root <- file.path(repo_root, "work", "genetics")
result_dir <- file.path(work_root, "discovery", "results")
fit_dir <- file.path(work_root, "signal_level", "fits")
pair_input_dir <- file.path(work_root, "signal_level", "input")

credible_set_rows <- function(fit, metadata) {
  if (is.null(fit$sets$cs) || length(fit$sets$cs) == 0) return(data.table())
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

pbc_summary <- fread(file.path(result_dir, "19_pbc_susie_multisignal_finemapping.tsv"))
pbc_rows <- list()
for (index in seq_len(nrow(pbc_summary))) {
  record <- pbc_summary[index]
  fit <- readRDS(file.path(fit_dir, paste0(record$locus_id, "__PBC__L", record$model_L, ".rds")))
  rows <- credible_set_rows(
    fit,
    list(locus_id = record$locus_id, pbc_locus = record$pbc_locus, model_L = record$model_L)
  )
  if (nrow(rows)) pbc_rows[[length(pbc_rows) + 1L]] <- rows
}
fwrite(rbindlist(pbc_rows, fill = TRUE),
       file.path(result_dir, "19b_pbc_susie_credible_set_details.tsv"), sep = "\t")

pair_audit <- fread(file.path(pair_input_dir, "signal_aware_allele_alignment_audit.tsv"))
pair_rows <- list()
for (index in seq_len(nrow(pair_audit))) {
  record <- pair_audit[index]
  fit_path <- file.path(
    fit_dir,
    paste0(record$analysis_id, "__", record$gene_symbol, "__", record$sample_group, ".rds")
  )
  fits <- readRDS(fit_path)
  for (trait in c("PBC", "eQTL")) {
    fit <- if (trait == "PBC") fits$pbc else fits$eqtl
    rows <- credible_set_rows(
      fit,
      list(
        analysis_id = record$analysis_id,
        pbc_locus = record$pbc_locus,
        trait = trait,
        gene_symbol = record$gene_symbol,
        sample_group = record$sample_group,
        model_L = 10L
      )
    )
    if (nrow(rows)) pair_rows[[length(pair_rows) + 1L]] <- rows
  }
}
fwrite(rbindlist(pair_rows, fill = TRUE),
       file.path(result_dir, "20b_signal_aware_susie_credible_set_details.tsv"), sep = "\t")

review_summary <- fread(file.path(result_dir, "20_onek1k_susie_concordance_l3_review.tsv"))
review_rows <- list()
for (index in seq_len(nrow(review_summary))) {
  record <- review_summary[index]
  fit_path <- file.path(
    fit_dir,
    paste0(record$analysis_id, "__", record$gene_symbol, "__", record$sample_group,
           "__L3_L3_nonconverged_review.rds")
  )
  fits <- readRDS(fit_path)
  for (trait in c("PBC", "eQTL")) {
    fit <- if (trait == "PBC") fits$pbc else fits$eqtl
    rows <- credible_set_rows(
      fit,
      list(
        analysis_id = record$analysis_id,
        pbc_locus = record$pbc_locus,
        trait = trait,
        gene_symbol = record$gene_symbol,
        sample_group = record$sample_group,
        model_L = 3L
      )
    )
    if (nrow(rows)) review_rows[[length(review_rows) + 1L]] <- rows
  }
}
fwrite(rbindlist(review_rows, fill = TRUE),
       file.path(result_dir, "20b_signal_aware_susie_credible_set_details_l3_review.tsv"), sep = "\t")

message("CREDIBLE_SET_DETAILS_REBUILT")
