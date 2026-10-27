suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_argument), winslash = "/")
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/")
work_root <- file.path(repo_root, "work", "genetics")
project <- file.path(work_root, "external_signal_coloc")
manifest_path <- file.path(work_root, "external_eqtl", "manifests", "external_signal_level_manifest.tsv")
cs_detail_path <- file.path(work_root, "external_eqtl", "results", "external_eqtl_susie_credible_set_details.tsv")
result_dir <- file.path(project, "results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

resolve_repository_path <- function(path) {
  if (grepl("^([A-Za-z]:[\\\\/]|/)", path)) return(path)
  file.path(repo_root, path)
}

read_gzip_tsv <- function(path) {
  fread(resolve_repository_path(path))
}

classification <- function(pp_h3, pp_h4) {
  if (is.na(pp_h4)) return("not_estimable")
  shared_fraction <- pp_h4 / (pp_h3 + pp_h4)
  if (pp_h4 >= 0.8 && shared_fraction >= 0.9) return("strong")
  if (pp_h4 >= 0.5) return("candidate")
  if (pp_h4 >= 0.2) return("exploratory")
  "no_support"
}

posterior_mass <- function(log_bf, prior) {
  variant_n <- length(log_bf) - 1L
  prior_vector <- c(rep(prior, variant_n), 1 - variant_n * prior)
  prior_vector[prior_vector <= 0] <- 1e-16
  log_weight <- log_bf + log(prior_vector)
  log_weight <- log_weight - max(log_weight)
  weight <- exp(log_weight)
  weight / sum(weight)
}

locus_config <- list(
  "1q23.1" = list(
    fit = file.path(work_root, "signal_level", "fits", "PBC03__PBC__L10.rds"),
    map = file.path(work_root, "signal_level", "input", "SA05__FCRL3__B_naive.tsv.gz")
  ),
  "1p31.3" = list(
    fit = file.path(work_root, "signal_level", "fits", "PBC04__PBC__L10.rds"),
    map = file.path(work_root, "signal_level", "input", "SA09__IL12RB2__NK.tsv.gz")
  ),
  "18q22.2" = list(
    fit = file.path(work_root, "signal_level", "fits", "PBC06__PBC__L10.rds"),
    map = file.path(work_root, "signal_level", "input", "SA13__CD226__CD4_TCM.tsv.gz")
  )
)

manifest <- fread(manifest_path)
manifest <- manifest[extraction_status == "success"]
cs_detail <- fread(cs_detail_path)
priors <- c(1e-6, 5e-6, 1e-5)

fit_cache <- list()
map_cache <- list()
for (locus in names(locus_config)) {
  fit <- readRDS(locus_config[[locus]]$fit)
  if (is.null(fit$sets$cs) || length(fit$sets$cs) == 0L) {
    stop("PBC fit has no credible sets: ", locus)
  }
  fit_cache[[locus]] <- fit

  map <- read_gzip_tsv(locus_config[[locus]]$map)
  map[, variant_grch38 := paste0(
    "chr", chromosome_grch37, "_", eqtl_position_grch38_1based, "_",
    eqtl_ref_original, "_", eqtl_alt_original
  )]
  if (map[, anyDuplicated(rsid)] > 0L) stop("Duplicate rsID in build map: ", locus)
  if (map[, anyDuplicated(variant_grch38)] > 0L) stop("Duplicate GRCh38 variant in build map: ", locus)
  map_cache[[locus]] <- map
}

qc_rows <- list()
result_rows <- list()

for (row_index in seq_len(nrow(manifest))) {
  record <- manifest[row_index]
  locus <- record$locus
  fit <- fit_cache[[locus]]
  pbc_components <- fit$sets$cs_index
  pbc_bf <- fit$lbf_variable[pbc_components, , drop = FALSE]
  rownames(pbc_bf) <- paste0("PBC_L", pbc_components)

  external <- fread(resolve_repository_path(record$lbf_selected_path))
  external <- external[molecular_trait_id == record$ensembl_gene_id]
  component_ids <- as.integer(strsplit(as.character(record$external_cs_components), ";", fixed = TRUE)[[1]])
  component_columns <- paste0("lbf_variable", component_ids)
  if (!all(component_columns %in% names(external))) {
    stop("Missing external log-BF component: ", record$dataset_id, " ", record$gene_symbol)
  }

  pair_cs <- cs_detail[
    dataset_id == record$dataset_id & gene_symbol == record$gene_symbol
  ]
  base_map <- map_cache[[locus]][, .(
    variant = variant_grch38,
    rsid,
    pbc_beta_alt,
    pbc_se,
    pbc_pvalue,
    reference_allele,
    alternate_allele,
    position_grch37_1based,
    position_grch38_1based = eqtl_position_grch38_1based
  )]
  variant_map <- copy(base_map)
  conflicts <- variant_map[, .(rsid_n = uniqueN(rsid)), by = variant][rsid_n > 1L, variant]
  if (length(conflicts)) variant_map <- variant_map[!variant %in% conflicts]
  setkey(variant_map, variant)

  external[, rsid := variant_map[external$variant, rsid]]
  mapped <- external[!is.na(rsid) & nzchar(rsid)]
  duplicate_rsids <- mapped[, .N, by = rsid][N > 1L, rsid]
  if (length(duplicate_rsids)) mapped <- mapped[!rsid %in% duplicate_rsids]
  mapped <- mapped[rsid %in% colnames(pbc_bf)]

  external_bf <- t(as.matrix(mapped[, ..component_columns]))
  rownames(external_bf) <- paste0("external_L", component_ids)
  colnames(external_bf) <- mapped$rsid
  pbc_bf_with_null <- cbind(pbc_bf, null = 0)
  external_bf_with_null <- cbind(external_bf, null = 0)
  common_rsids <- setdiff(intersect(colnames(pbc_bf_with_null), colnames(external_bf_with_null)), "null")

  pbc_coverages <- sapply(seq_len(nrow(pbc_bf_with_null)), function(index) {
    posterior <- posterior_mass(pbc_bf_with_null[index, ], 1e-4)
    sum(posterior[match(common_rsids, colnames(pbc_bf_with_null))], na.rm = TRUE)
  })
  external_coverages <- sapply(seq_len(nrow(external_bf_with_null)), function(index) {
    posterior <- posterior_mass(external_bf_with_null[index, ], 1e-4)
    sum(posterior[match(common_rsids, colnames(external_bf_with_null))], na.rm = TRUE)
  })

  qc_rows[[length(qc_rows) + 1L]] <- data.table(
    analysis_order = record$analysis_order,
    locus = locus,
    dataset_id = record$dataset_id,
    study_id = record$study_id,
    study_label = record$study_label,
    sample_group = record$sample_group,
    sample_size = record$sample_size,
    ancestry = record$ancestry,
    gene_symbol = record$gene_symbol,
    gene_id = record$ensembl_gene_id,
    external_variant_n = nrow(external),
    exact_build_allele_mapped_n = nrow(mapped),
    ambiguous_rsid_dropped_n = length(duplicate_rsids),
    common_variant_n = length(common_rsids),
    pbc_component_n = nrow(pbc_bf),
    external_component_n = nrow(external_bf),
    minimum_pbc_signal_overlap_posterior = min(pbc_coverages),
    minimum_external_signal_overlap_posterior = min(external_coverages),
    overlap_qc = ifelse(
      length(common_rsids) >= 100 & min(pbc_coverages) >= 0.5 & min(external_coverages) >= 0.5,
      "pass", "review"
    )
  )

  if (length(common_rsids) == 0L) next
  for (p12 in priors) {
    coloc_result <- suppressWarnings(coloc::coloc.bf_bf(
      pbc_bf_with_null,
      external_bf_with_null,
      p1 = 1e-4,
      p2 = 1e-4,
      p12 = p12,
      overlap.min = 0.5,
      trim_by_posterior = TRUE
    ))
    if (!is.list(coloc_result) || is.null(coloc_result$summary)) next
    summary_table <- as.data.table(coloc_result$summary)
    posterior_table <- as.data.table(coloc_result$results)
    for (signal_row in seq_len(nrow(summary_table))) {
      summary_record <- summary_table[signal_row]
      posterior_column <- names(posterior_table)[signal_row + 1L]
      top_index <- which.max(posterior_table[[posterior_column]])
      top_shared_snp <- posterior_table$snp[top_index]
      top_shared_posterior <- posterior_table[[posterior_column]][top_index]
      pbc_component <- pbc_components[summary_record$idx1]
      external_component <- component_ids[summary_record$idx2]
      pbc_variant <- base_map[rsid == top_shared_snp][1L]
      external_variant <- pair_cs[component == external_component & rsid == top_shared_snp][1L]
      result_rows[[length(result_rows) + 1L]] <- data.table(
        analysis_order = record$analysis_order,
        locus = locus,
        dataset_id = record$dataset_id,
        study_id = record$study_id,
        study_label = record$study_label,
        sample_group = record$sample_group,
        sample_size = record$sample_size,
        ancestry = record$ancestry,
        gene_symbol = record$gene_symbol,
        gene_id = record$ensembl_gene_id,
        pbc_signal = pbc_component,
        external_signal = external_component,
        p12 = p12,
        common_variant_n = length(common_rsids),
        pbc_overlap_posterior = pbc_coverages[summary_record$idx1],
        external_overlap_posterior = external_coverages[summary_record$idx2],
        PP.H0 = summary_record$PP.H0.abf,
        PP.H1 = summary_record$PP.H1.abf,
        PP.H2 = summary_record$PP.H2.abf,
        PP.H3 = summary_record$PP.H3.abf,
        PP.H4 = summary_record$PP.H4.abf,
        shared_fraction_H4_over_H3_H4 = summary_record$PP.H4.abf /
          (summary_record$PP.H3.abf + summary_record$PP.H4.abf),
        top_pbc_snp = summary_record$hit1,
        top_external_snp = summary_record$hit2,
        top_shared_snp = top_shared_snp,
        top_shared_snp_posterior = top_shared_posterior,
        pbc_beta_alt = ifelse(nrow(pbc_variant), pbc_variant$pbc_beta_alt, NA_real_),
        pbc_alt = ifelse(nrow(pbc_variant), pbc_variant$alternate_allele, NA_character_),
        external_variant_grch38 = ifelse(nrow(external_variant), external_variant$variant, NA_character_),
        external_beta_alt = ifelse(nrow(external_variant), external_variant$beta, NA_real_),
        evidence_class = classification(summary_record$PP.H3.abf, summary_record$PP.H4.abf),
        status = ifelse(is.na(summary_record$PP.H4.abf), "insufficient_overlap", "estimated")
      )
    }
  }
  message("DONE\t", record$dataset_id, "\t", record$gene_symbol, "\tcommon=", length(common_rsids))
}

qc <- rbindlist(qc_rows, fill = TRUE)
results <- rbindlist(result_rows, fill = TRUE)
setorder(qc, analysis_order)
setorder(results, analysis_order, p12, -PP.H4)
fwrite(qc, file.path(result_dir, "01_external_signal_level_dataset_qc.tsv"), sep = "\t")
fwrite(results, file.path(result_dir, "02_external_signal_pair_prior_sensitivity.tsv"), sep = "\t")
fwrite(results[p12 == 1e-6], file.path(result_dir, "03_external_signal_pair_primary.tsv"), sep = "\t")
