suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_argument), winslash = "/")
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/")
work_root <- file.path(repo_root, "work", "genetics")
region_path <- file.path(repo_root, "data", "raw", "genetics", "finngen", "r12", "chirbil_prim", "regions", "02_1p31.3_IL12RB2.tsv.gz")
manifest_path <- file.path(work_root, "external_eqtl", "manifests", "external_signal_level_manifest.tsv")
result_dir <- file.path(work_root, "external_signal_coloc", "results")
ld_path <- file.path(result_dir, "12_IL12RB2_SiSu42_lead_LD.tsv")

resolve_repository_path <- function(path) {
  if (grepl("^([A-Za-z]:[\\\\/]|/)", path)) return(path)
  file.path(repo_root, path)
}

classification <- function(pp_h3, pp_h4) {
  if (is.na(pp_h4)) return("not_estimable")
  shared_fraction <- pp_h4 / (pp_h3 + pp_h4)
  if (pp_h4 >= 0.8 && shared_fraction >= 0.9) return("strong")
  if (pp_h4 >= 0.5) return("candidate")
  if (pp_h4 >= 0.2) return("exploratory")
  "no_support"
}

softmax <- function(log_weight) {
  shifted <- log_weight - max(log_weight)
  weight <- exp(shifted)
  weight / sum(weight)
}

credible_set <- function(variant, posterior, coverage = 0.95) {
  order_index <- order(posterior, decreasing = TRUE)
  ordered_posterior <- posterior[order_index]
  crossing <- which(cumsum(ordered_posterior) >= coverage)[1L]
  selected <- order_index[seq_len(crossing)]
  data.table(
    variant = variant[selected],
    posterior = posterior[selected],
    rank = seq_along(selected),
    cumulative_posterior = cumsum(ordered_posterior[seq_len(crossing)])
  )
}

finngen <- fread(cmd = paste("gzip -dc", shQuote(region_path)))
setnames(finngen, "#chrom", "chromosome")
finngen[, variant := paste0("chr", chromosome, "_", pos, "_", ref, "_", alt)]
finngen[, variant_colon := paste(chromosome, pos, ref, alt, sep = ":")]
finngen <- finngen[
  is.finite(beta) & is.finite(sebeta) & sebeta > 0 & is.finite(pval) &
    !is.na(ref) & !is.na(alt)
]
if (finngen[, anyDuplicated(variant)] > 0L) stop("Duplicate exact GRCh38 REF/ALT variants in FinnGen input")

prior_sd_values <- c(0.15, 0.20, 0.30)
credible_rows <- list()
for (prior_sd in prior_sd_values) {
  variance <- finngen$sebeta^2
  prior_variance <- prior_sd^2
  log_abf <- 0.5 * (
    log(variance / (variance + prior_variance)) +
      (finngen$beta / finngen$sebeta)^2 * prior_variance / (variance + prior_variance)
  )
  posterior <- softmax(log_abf)
  cs <- credible_set(finngen$variant, posterior, coverage = 0.95)
  cs <- merge(
    cs,
    finngen[, .(
      variant,
      variant_colon,
      chromosome,
      position_grch38_1based = pos,
      ref,
      alt,
      rsids,
      nearest_genes,
      beta_alt = beta,
      se = sebeta,
      pvalue = pval,
      af_alt
    )],
    by = "variant",
    all.x = TRUE,
    sort = FALSE
  )
  setorder(cs, rank)
  cs[, `:=`(
    prior_sd_log_odds = prior_sd,
    prior_variance = prior_variance,
    finngen_model = "single-causal-variant Wakefield ABF; regional posterior conditional on one association",
    genome_build = "GRCh38",
    coordinate_convention = "1-based VCF-style REF/ALT"
  )]
  credible_rows[[as.character(prior_sd)]] <- cs
  if (abs(prior_sd - 0.20) < 1e-12) {
    finngen[, `:=`(log_abf = log_abf, posterior_conditional = posterior)]
  }
}

credible_output <- rbindlist(credible_rows, use.names = TRUE, fill = TRUE, idcol = "prior_sd_key")
setcolorder(credible_output, c(
  "prior_sd_key", "prior_sd_log_odds", "prior_variance", "rank", "variant",
  "variant_colon", "chromosome", "position_grch38_1based", "ref", "alt", "rsids",
  "nearest_genes", "beta_alt", "se", "pvalue", "af_alt", "posterior",
  "cumulative_posterior", "finngen_model", "genome_build", "coordinate_convention"
))
fwrite(
  credible_output,
  file.path(result_dir, "16_il12rb2_finngen_single_causal_abf_95cs.tsv"),
  sep = "\t"
)

manifest <- fread(manifest_path)
targets <- manifest[locus == "1p31.3" & gene_symbol == "IL12RB2" & extraction_status == "success"]
if (nrow(targets) != 3L) stop("Expected exactly three external IL12RB2 datasets")

priors_p12 <- c(1e-6, 5e-6, 1e-5)
result_rows <- list()
qc_rows <- list()
overlap_rows <- list()

finngen_cs <- credible_rows[["0.2"]]
finngen_cs_variants <- finngen_cs$variant
lead_ld <- fread(ld_path)
lead_ld[, variant := paste0(
  "chr", chromosome, "_", position_grch38_1based, "_", ref, "_", alt
)]
if (lead_ld[, anyDuplicated(variant)] > 0L) stop("Duplicate exact variants in lead LD audit")
setkey(lead_ld, variant)

for (row_index in seq_len(nrow(targets))) {
  record <- targets[row_index]
  external_path <- resolve_repository_path(record$lbf_selected_path)
  external_path <- gsub("\\\\", "/", external_path)
  external <- fread(external_path)
  external <- external[molecular_trait_id == record$ensembl_gene_id]
  component_column <- paste0("lbf_variable", record$external_cs_components)
  if (length(component_column) != 1L || !component_column %in% names(external)) {
    stop("Unexpected external component configuration for ", record$dataset_id)
  }
  if (external[, anyDuplicated(variant)] > 0L) stop("Duplicate external variant: ", record$dataset_id)
  external[, external_log_bf := get(component_column)]
  external <- external[is.finite(external_log_bf)]
  external[, posterior_conditional := softmax(external_log_bf)]

  common_variants <- intersect(finngen$variant, external$variant)
  if (length(common_variants) < 100L) stop("Insufficient exact variant overlap: ", record$dataset_id)
  finngen_overlap <- finngen[match(common_variants, variant)]
  external_overlap <- external[match(common_variants, variant)]
  if (!identical(finngen_overlap$variant, external_overlap$variant)) stop("Variant order mismatch")

  finngen_overlap_coverage <- sum(finngen_overlap$posterior_conditional)
  external_overlap_coverage <- sum(external_overlap$posterior_conditional)
  finngen_bf <- matrix(
    finngen_overlap$log_abf,
    nrow = 1L,
    dimnames = list("FinnGen_single_causal_ABF", common_variants)
  )
  external_bf <- matrix(
    external_overlap$external_log_bf,
    nrow = 1L,
    dimnames = list(paste0("external_L", record$external_cs_components), common_variants)
  )
  finngen_bf <- cbind(finngen_bf, null = 0)
  external_bf <- cbind(external_bf, null = 0)

  external_cs <- credible_set(
    external$variant,
    external$posterior_conditional,
    coverage = 0.95
  )
  shared_cs <- intersect(finngen_cs_variants, external_cs$variant)
  overlap_rows[[length(overlap_rows) + 1L]] <- data.table(
    dataset_id = record$dataset_id,
    study_id = record$study_id,
    study_label = record$study_label,
    sample_group = record$sample_group,
    ancestry = record$ancestry,
    gene_symbol = record$gene_symbol,
    finngen_95cs_variant_n = nrow(finngen_cs),
    external_95cs_variant_n = nrow(external_cs),
    shared_95cs_variant_n = length(shared_cs),
    finngen_posterior_mass_in_external_95cs = sum(
      finngen$posterior_conditional[match(external_cs$variant, finngen$variant)],
      na.rm = TRUE
    ),
    external_posterior_mass_in_finngen_95cs = sum(
      external$posterior_conditional[match(finngen_cs_variants, external$variant)],
      na.rm = TRUE
    ),
    shared_95cs_variants = paste(shared_cs, collapse = ";"),
    finngen_model = "single-causal-variant Wakefield ABF, prior SD=0.20 log-odds",
    external_model = paste0("SuSiE component ", record$external_cs_components)
  )

  qc_rows[[length(qc_rows) + 1L]] <- data.table(
    dataset_id = record$dataset_id,
    study_id = record$study_id,
    study_label = record$study_label,
    sample_group = record$sample_group,
    sample_size = record$sample_size,
    ancestry = record$ancestry,
    gene_symbol = record$gene_symbol,
    finngen_region_variant_n = nrow(finngen),
    external_variant_n = nrow(external),
    exact_grch38_ref_alt_common_variant_n = length(common_variants),
    finngen_conditional_posterior_covered_by_overlap = finngen_overlap_coverage,
    external_conditional_posterior_covered_by_overlap = external_overlap_coverage,
    coordinate_match = "exact GRCh38 position+REF+ALT",
    finngen_prior_sd_log_odds = 0.20,
    finngen_assumption = "one causal variant in region; no signed Finnish LD matrix available",
    external_component = record$external_cs_components,
    overlap_qc = ifelse(
      finngen_overlap_coverage >= 0.5 & external_overlap_coverage >= 0.5,
      "pass", "review"
    ),
    analysis_label = "BF-BF colocalization; not coloc.susie"
  )

  for (p12 in priors_p12) {
    coloc_result <- suppressWarnings(coloc::coloc.bf_bf(
      finngen_bf,
      external_bf,
      p1 = 1e-4,
      p2 = 1e-4,
      p12 = p12,
      overlap.min = 0.5,
      trim_by_posterior = TRUE
    ))
    summary_row <- as.data.table(coloc_result$summary)[1L]
    posterior_table <- as.data.table(coloc_result$results)
    posterior_column <- names(posterior_table)[2L]
    top_index <- which.max(posterior_table[[posterior_column]])
    top_shared_variant <- posterior_table$snp[top_index]
    top_finngen <- finngen[variant == top_shared_variant][1L]
    top_external <- external[variant == top_shared_variant][1L]
    external_top_ld <- lead_ld[summary_row$hit2]
    shared_top_ld <- lead_ld[top_shared_variant]
    result_rows[[length(result_rows) + 1L]] <- data.table(
      dataset_id = record$dataset_id,
      study_id = record$study_id,
      study_label = record$study_label,
      sample_group = record$sample_group,
      sample_size = record$sample_size,
      ancestry = record$ancestry,
      gene_symbol = record$gene_symbol,
      p1 = 1e-4,
      p2 = 1e-4,
      p12 = p12,
      common_variant_n = length(common_variants),
      finngen_overlap_posterior = finngen_overlap_coverage,
      external_overlap_posterior = external_overlap_coverage,
      PP.H0 = summary_row$PP.H0.abf,
      PP.H1 = summary_row$PP.H1.abf,
      PP.H2 = summary_row$PP.H2.abf,
      PP.H3 = summary_row$PP.H3.abf,
      PP.H4 = summary_row$PP.H4.abf,
      shared_fraction_H4_over_H3_H4 = summary_row$PP.H4.abf /
        (summary_row$PP.H3.abf + summary_row$PP.H4.abf),
      finngen_top_variant = summary_row$hit1,
      external_top_variant = summary_row$hit2,
      external_top_r2_to_finngen_region_top = ifelse(
        nrow(external_top_ld), external_top_ld$r2_to_finngen_region_top, NA_real_
      ),
      external_top_r2_to_discovery_prespecified_lead = ifelse(
        nrow(external_top_ld), external_top_ld$r2_to_discovery_prespecified_lead, NA_real_
      ),
      top_shared_variant = top_shared_variant,
      top_shared_variant_posterior = posterior_table[[posterior_column]][top_index],
      top_shared_rsid = ifelse(nrow(top_finngen), top_finngen$rsids, NA_character_),
      top_shared_finngen_beta_alt = ifelse(nrow(top_finngen), top_finngen$beta, NA_real_),
      top_shared_finngen_alt = ifelse(nrow(top_finngen), top_finngen$alt, NA_character_),
      top_shared_external_log_bf = ifelse(nrow(top_external), top_external$external_log_bf, NA_real_),
      top_shared_r2_to_finngen_region_top = ifelse(
        nrow(shared_top_ld), shared_top_ld$r2_to_finngen_region_top, NA_real_
      ),
      top_shared_r2_to_discovery_prespecified_lead = ifelse(
        nrow(shared_top_ld), shared_top_ld$r2_to_discovery_prespecified_lead, NA_real_
      ),
      evidence_class = classification(summary_row$PP.H3.abf, summary_row$PP.H4.abf),
      finngen_model = "single-causal-variant Wakefield ABF, prior SD=0.20 log-odds",
      external_model = paste0("SuSiE component ", record$external_cs_components),
      method = "coloc.bf_bf",
      interpretation_guardrail = "BF-BF colocalization; not coloc.susie because a complete signed Finnish LD matrix was unavailable"
    )
  }
}

results <- rbindlist(result_rows, use.names = TRUE, fill = TRUE)
qc <- rbindlist(qc_rows, use.names = TRUE, fill = TRUE)
overlap <- rbindlist(overlap_rows, use.names = TRUE, fill = TRUE)
setorder(results, dataset_id, p12)
setorder(qc, dataset_id)
setorder(overlap, dataset_id)

fwrite(results, file.path(result_dir, "17_il12rb2_finngen_external_eqtl_bf_bf_sensitivity.tsv"), sep = "\t")
fwrite(qc, file.path(result_dir, "18_il12rb2_finngen_external_eqtl_bf_bf_qc.tsv"), sep = "\t")
fwrite(overlap, file.path(result_dir, "19_il12rb2_finngen_external_eqtl_95cs_overlap.tsv"), sep = "\t")

print(results[p12 == 1e-6, .(
  dataset_id,
  study_label,
  sample_group,
  common_variant_n,
  PP.H3,
  PP.H4,
  shared_fraction_H4_over_H3_H4,
  top_shared_variant,
  top_shared_rsid,
  evidence_class
)])
print(overlap[, .(
  dataset_id,
  finngen_95cs_variant_n,
  external_95cs_variant_n,
  shared_95cs_variant_n,
  finngen_posterior_mass_in_external_95cs,
  external_posterior_mass_in_finngen_95cs
)])
