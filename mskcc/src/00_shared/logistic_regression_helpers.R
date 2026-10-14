# ============================================================================
# Shared helper: sourced by the analysis scripts, never run on its own.
# ============================================================================
###################
# Logistic regression of mutation status against a variable of interest
# Assumes mutation_mtx is pre-joined with clinical covariates (see calling script)
###################

impact_logistic_regression <- function(
    variable_of_interest = "BMI_NEAREST_RESULT",
    extra_covariates = c(),
    mutation_mtx,
    seq_bmi_primary_mut,
    single_sex_cancers,
    gene_in_panel,
    min_mutations = 10) {
  timing <- Sys.time()

  # ── validate required columns exist in mutation_mtx ───────────────────────
  required_cols <- c(
    variable_of_interest, "age", "sex", "ancestry",
    "stage", "smoking", "prior_treatment"
  )
  missing_cols <- setdiff(required_cols, colnames(mutation_mtx))
  if (length(missing_cols) > 0) {
    stop(
      "mutation_mtx is missing required columns: ", paste(missing_cols, collapse = ", "),
      "\nAvailable columns (first 20): ", paste(head(colnames(mutation_mtx), 20), collapse = ", ")
    )
  }

  # Sanitise extra_covariates: unique, present in matrix, not duplicating variable_of_interest
  extra_covariates <- intersect(
    unique(setdiff(as.character(extra_covariates), variable_of_interest)),
    colnames(mutation_mtx)
  )
  if (length(setdiff(as.character(extra_covariates), colnames(mutation_mtx))) > 0) {
    message(
      "Dropping extra_covariates not found in mutation_mtx: ",
      paste(setdiff(extra_covariates, colnames(mutation_mtx)), collapse = ", ")
    )
  }

  # Fixed clinical covariates (always present after the pre-join)
  clinical_cols <- c(
    "age", "sex", "ancestry", "cancer", "stage",
    "smoking", "prior_treatment", "sample_type"
  )

  # Pre-split by cancer type once
  cancer_df_list <- split(seq_bmi_primary_mut, seq_bmi_primary_mut$CANCER_TYPE_DETAILED)

  results <- list()
  i <- 1L

  for (cancer_df in cancer_df_list) {
    cancer <- cancer_df$CANCER_TYPE_DETAILED[[1]]
    message("Cancer: ", cancer)

    is_single_sex <- cancer %in% single_sex_cancers$CANCER_TYPE_DETAILED

    # Genes observed in this cancer (drop NA = unsequenced samples)
    genes <- unique(na.omit(cancer_df$Hugo_Symbol))

    for (gene in genes) {
      # ── panel filtering ────────────────────────────────────────────────────
      if (!gene %in% gene_in_panel$gene_name) {
        message("Skipping ", gene, ": not found in gene_in_panel for ", cancer)
        next
      }
      panels <- toupper(gene_in_panel$panels[gene_in_panel$gene_name == gene][[1]][[1]])

      tsbs <- unique(cancer_df$Tumor_Sample_Barcode[
        cancer_df$GENE_PANEL %in% panels &
          cancer_df$CANCER_TYPE_DETAILED == cancer
      ])

      df <- mutation_mtx[tsbs, , drop = FALSE]

      # ── minimum mutation guard ─────────────────────────────────────────────
      if (!gene %in% colnames(df) || sum(df[, gene], na.rm = TRUE) < min_mutations) next

      # ── build analysis data frame ──────────────────────────────────────────
      adf <- data.frame(
        mut = df[, gene],
        test_variable = df[, variable_of_interest],
        age = df[, "age"],
        sex = df[, "sex"],
        ancestry = df[, "ancestry"],
        stage = df[, "stage"],
        smoking = df[, "smoking"],
        prior_treatment = df[, "prior_treatment"],
        stringsAsFactors = FALSE
      )

      # append extra covariates
      for (cov in extra_covariates) adf[[cov]] <- df[, cov]

      # ── drop single-level columns (can happen during downsampling) ─────────
      single_level <- names(adf)[sapply(adf, function(x) length(unique(x)) == 1)]
      single_level <- setdiff(single_level, c("mut", "sex")) # sex handled separately
      if (length(single_level) > 0) {
        message(
          "Dropping single-level columns for ", gene, " / ", cancer, ": ",
          paste(single_level, collapse = ", ")
        )
        adf <- adf[, !names(adf) %in% single_level, drop = FALSE]
      }

      # ── decide whether sex enters the model ───────────────────────────────
      use_sex <- !is_single_sex && length(unique(adf$sex)) > 1

      partial_vars <- c("test_variable", "age", if (use_sex) "sex", extra_covariates)
      # Self-reported race is deliberately NOT a covariate. Genetic ancestry
      # (`ancestry`) is already in the model, and race was the only covariate in
      # this study fit here and nowhere else -- every other adjusted model in the
      # manuscript adjusts for age, sex, stage, ancestry and smoking without it.
      # Keeping it made the gene-level Methods formula unstateable alongside the
      # rest of the paper's.
      full_vars <- c(
        "test_variable", "age", if (use_sex) "sex",
        "ancestry", "stage", "smoking", "prior_treatment",
        extra_covariates
      )
      # remove any that were dropped above
      partial_vars <- intersect(partial_vars, names(adf))
      full_vars <- intersect(full_vars, names(adf))

      # ── fit models ────────────────────────────────────────────────────────
      fit_base <- glm(mut ~ test_variable, data = adf, family = "binomial")
      fit_partial <- glm(reformulate(partial_vars, "mut"), data = adf, family = "binomial")
      fit_full <- glm(reformulate(full_vars, "mut"),
        data = adf, family = "binomial",
        na.action = "na.omit"
      )

      # ── extract row 2 coefficients (test_variable) ────────────────────────
      extract_coef <- function(fit, n_mut, n_total) {
        s <- summary(fit)$coefficients
        list(
          cancer       = cancer,
          gene         = gene,
          panels       = paste(panels, collapse = ", "),
          sample_size  = n_total,
          num_mutated  = n_mut,
          estimate     = s[2, 1],
          se           = s[2, 2],
          pvalue       = s[2, 4]
        )
      }

      n_total <- nrow(adf)
      n_mut <- sum(adf$mut, na.rm = TRUE)
      results[["base"]][[i]] <- extract_coef(fit_base, n_mut, n_total)
      results[["partial"]][[i]] <- extract_coef(fit_partial, n_mut, n_total)
      results[["full"]][[i]] <- extract_coef(fit_full, n_mut, n_total)
      i <- i + 1L
    }
  }

  message("Total time: ", round(difftime(Sys.time(), timing, units = "mins"), 2), " mins")

  # ── collate results ────────────────────────────────────────────────────────
  collate <- function(lst, model_label) {
    bind_rows(lst) %>%
      as.data.frame() %>%
      mutate(
        model         = model_label,
        test_variable = variable_of_interest,
        padj          = p.adjust(pvalue, method = "BH"),
        odds_ratio    = exp(estimate),
        ci_lower      = exp(estimate - 1.96 * se),
        ci_upper      = exp(estimate + 1.96 * se)
      ) %>%
      arrange(padj)
  }

  list(
    base_results    = collate(results[["base"]], "base"),
    mv_results      = collate(results[["partial"]], "partially_adjusted"),
    mv_model_stage  = collate(results[["full"]], "fully_adjusted")
  )
}


save_impact_logistic_regression <- function(
    impact_logistic_regression_results_obj,
    results_dir = "./code_outputs/results/tables/bmi/impact/gene_cancer_bmi_associations",
    analysis_name) {
  out_dir <- file.path(results_dir, analysis_name)
  create_dir_if_needed(out_dir)

  # Only the fully adjusted model is written. impact_logistic_regression() fits
  # all three nested models the Methods describe (unadjusted, + age/sex, and the
  # full covariate set) and returns them in `res`, but the manuscript reports the
  # fully adjusted fit alone -- it is what Figure 4a's volcano and Supplementary
  # Table 18 are built from -- so the other two are not saved.
  res <- impact_logistic_regression_results_obj
  write_tsv(
    res$mv_model_stage,
    file.path(out_dir, paste0(analysis_name, "_mv_model_stage_results.tsv"))
  )
}
