# =============================================================================
# RPPAnalyzeR :: 06_pipeline.R
# One-function pipeline: import → QC → DE → visualize → export
# =============================================================================

#' Run the complete RPPA analysis pipeline
#'
#' A single entry-point that chains import → QC → differential expression →
#' visualizations → Excel export. Results and plots are written to `output_dir`.
#'
#' @param xlsx_path Character. Path to the RPPA Core Excel file.
#' @param output_dir Character. Directory for all outputs (created if absent).
#' @param min_qc_score Numeric. Minimum antibody QC score (default 0.8).
#' @param remove_caution Logical. Exclude "C" (caution) antibodies (default FALSE).
#' @param remove_low_protein Logical. Exclude low-protein samples (default TRUE).
#' @param ctrl_timepoint Character. Baseline timepoint label (default "0h").
#' @param p_threshold Numeric. DE significance threshold (default 0.05).
#' @param fc_threshold Numeric. |Log2FC| threshold (default 0.5).
#' @param key_proteins Character vector. Proteins for time-course plot.
#'   Defaults to top 6 significant proteins at 24h.
#' @param verbose Logical. Progress messages (default TRUE).
#'
#' @return Invisibly returns a list with `rppa`, `de_results`, and plot paths.
#' @export
#'
#' @examples
#' \dontrun{
#' results <- run_pipeline(
#'   xlsx_path  = "01_Jonathan_Coloff__Vipin_Rawat.xlsx",
#'   output_dir = "rppa_analysis_output"
#' )
#' }
run_pipeline <- function(xlsx_path,
                          output_dir         = "rppa_output",
                          min_qc_score       = 0.8,
                          remove_caution     = FALSE,
                          remove_low_protein = TRUE,
                          ctrl_timepoint     = "0h",
                          p_threshold        = 0.05,
                          fc_threshold       = 0.5,
                          key_proteins       = NULL,
                          verbose            = TRUE) {

  msg <- function(...) if (verbose) message("\n[Pipeline] ", ...)

  # ── Setup output directory ───────────────────────────────────────────────
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir)

  msg("════════════════════════════════════")
  msg("  RPPAnalyzeR — Full Pipeline")
  msg("════════════════════════════════════")
  msg("Input : ", xlsx_path)
  msg("Output: ", output_dir)

  # ── Step 1: Import ───────────────────────────────────────────────────────
  msg("STEP 1 / 5 — Importing data ...")
  rppa <- import_rppa(xlsx_path, verbose = verbose)

  # ── Step 2: QC ──────────────────────────────────────────────────────────
  msg("STEP 2 / 5 — Quality control ...")
  rppa <- run_qc(rppa,
                 min_qc_score       = min_qc_score,
                 remove_caution     = remove_caution,
                 remove_low_protein = remove_low_protein,
                 verbose            = verbose)
  print_qc_summary(rppa)

  # ── Step 3: Differential Expression ─────────────────────────────────────
  msg("STEP 3 / 5 — Differential expression ...")
  de_all <- diff_expression_all(rppa,
                                 ctrl_timepoint = ctrl_timepoint,
                                 p_threshold    = p_threshold,
                                 fc_threshold   = fc_threshold)

  n_sig <- sum(de_all$Significant, na.rm = TRUE)
  msg("  Total significant hits across all timepoints: ", n_sig)

  # ── Step 4: Visualizations ───────────────────────────────────────────────
  msg("STEP 4 / 5 — Generating plots ...")
  plot_paths <- list()

  # 4a. Sample QC bar chart
  p_sampleqc <- plot_sample_qc(rppa)
  path_sqc <- file.path(plots_dir, "01_sample_QC.png")
  save_plot(p_sampleqc, path_sqc, width = 10, height = 5)
  plot_paths$sample_qc <- path_sqc

  # 4b. Antibody QC histogram
  p_abqc <- plot_antibody_qc(rppa)
  path_abqc <- file.path(plots_dir, "02_antibody_QC.png")
  save_plot(p_abqc, path_abqc, width = 8, height = 5)
  plot_paths$antibody_qc <- path_abqc

  # 4c. Volcano plots for each time point
  timepoints <- setdiff(unique(rppa$sample_info$Timepoint), ctrl_timepoint)
  for (tp in timepoints) {
    de_tp <- dplyr::filter(de_all, Timepoint == tp)
    if (nrow(de_tp) == 0) next
    p_vol <- plot_volcano(de_tp, timepoint = tp,
                          p_threshold = p_threshold, fc_threshold = fc_threshold)
    path_vol <- file.path(plots_dir, paste0("03_volcano_0h_vs_", tp, ".png"))
    save_plot(p_vol, path_vol, width = 8, height = 7)
    plot_paths[[paste0("volcano_", tp)]] <- path_vol
  }

  # 4d. Heatmap (top 50 most variable proteins)
  expr  <- get_expression_matrix(rppa, "l4_chm")
  mat   <- expr$matrix
  vars  <- apply(mat, 2, var, na.rm = TRUE)
  top50 <- names(sort(vars, decreasing = TRUE))[1:min(50, length(vars))]
  top50_clean <- stringr::str_remove(top50, "-[RMG]-[VCE]$")
  path_hm <- file.path(plots_dir, "04_heatmap_top50.png")
  plot_heatmap(rppa, proteins = top50_clean,
               title    = "Top 50 Most Variable Proteins (L4 CHM)",
               filename = path_hm, width = 10, height = 14)
  plot_paths$heatmap <- path_hm
  msg("  Heatmap saved.")

  # 4e. Time-course for key proteins
  if (is.null(key_proteins)) {
    de_24h <- dplyr::filter(de_all, Timepoint == "24h")
    key_proteins <- top_proteins(de_24h, n = 6) %>% dplyr::pull(Protein)
    key_proteins <- stringr::str_remove(key_proteins, "-[RMG]-[VCE]$")
  }
  if (length(key_proteins) > 0) {
    p_tc <- plot_timecourse(rppa, proteins = key_proteins)
    path_tc <- file.path(plots_dir, "05_timecourse_key_proteins.png")
    save_plot(p_tc, path_tc, width = 9, height = 6)
    plot_paths$timecourse <- path_tc
  }

  # ── Step 5: Export ────────────────────────────────────────────────────────
  msg("STEP 5 / 5 — Exporting results ...")
  excel_path <- file.path(output_dir, "RPPAnalyzeR_results.xlsx")
  export_rppa_excel(rppa, de_results = de_all, output_path = excel_path)

  csv_path <- file.path(output_dir, "DE_significant_hits.csv")
  export_de_csv(de_all, output_path = csv_path, significant_only = TRUE)

  # ── Summary ───────────────────────────────────────────────────────────────
  msg("════════════════════════════════════")
  msg("  Pipeline complete!")
  msg("  Outputs in: ", normalizePath(output_dir))
  msg("  Excel     : ", basename(excel_path))
  msg("  Sig. hits : ", basename(csv_path), " (", n_sig, " rows)")
  msg("  Plots     : ", length(plot_paths), " files in /plots/")
  msg("════════════════════════════════════")

  invisible(list(
    rppa        = rppa,
    de_results  = de_all,
    plot_paths  = plot_paths,
    excel_path  = excel_path,
    csv_path    = csv_path
  ))
}
