# =============================================================================
# RPPAnalyzeR :: 02_qc.R
# Quality control filtering and summary reporting
#
# KEY DESIGN NOTE:
# The QC sheet uses human-readable names ("Cyclin B1") while the expression
# sheets use heatmap names ("Cyclin-B1-R-V"). These NEVER match directly.
# The shared key is Antigen_ID (e.g. "AGID00100") which is identical in both
# sheets and in the same order. All QC filtering uses Antigen_ID.
# =============================================================================

#' Run QC checks and filter the RPPA dataset
#'
#' @param rppa An `rppa_data` object from [import_rppa()].
#' @param min_qc_score Numeric. Minimum antibody QC score (default 0.8).
#' @param remove_caution Logical. Remove "C" validation antibodies (default FALSE).
#' @param remove_low_protein Logical. Remove samples with Total Protein < -3 (default TRUE).
#' @param verbose Logical. Print QC summary (default TRUE).
#' @return A filtered `rppa_data` list with a `qc_report` element added.
#' @export
run_qc <- function(rppa,
                   min_qc_score       = 0.8,
                   remove_caution     = FALSE,
                   remove_low_protein = TRUE,
                   verbose            = TRUE) {

  msg <- function(...) if (verbose) message("[RPPAnalyzeR::QC] ", ...)

  # ── 1. Antibody QC score filtering ───────────────────────────────────────
  # Match QC sheet to expression sheet via Antigen_ID (the only shared key).
  # QC sheet:        Antigen_ID column (e.g. "AGID00100")
  # Expression sheet: Antigen_ID in rppa$metadata (same values, same order)
  # Heatmap names:   rppa$metadata$Heatmap_Label (e.g. "Cyclin-B1-R-V")
  #                  — these ARE the column names in l4_log2 etc.

  ab_qc    <- rppa$antibody_qc    # has Antigen_ID and QC_Score
  meta     <- rppa$metadata       # has Antigen_ID and Heatmap_Label

  # Join QC scores onto metadata via Antigen_ID
  meta_qc <- merge(
    meta[, c("Antigen_ID", "Heatmap_Label", "Validation_Status")],
    ab_qc[, c("Antigen_ID", "QC_Score")],
    by   = "Antigen_ID",
    all.x = TRUE
  )

  pass_ab  <- meta_qc[!is.na(meta_qc$QC_Score) & meta_qc$QC_Score >= min_qc_score, ]
  fail_ab  <- meta_qc[!is.na(meta_qc$QC_Score) & meta_qc$QC_Score <  min_qc_score, ]
  no_score <- meta_qc[is.na(meta_qc$QC_Score), ]

  msg(nrow(pass_ab), " / ", nrow(meta_qc),
      " antibodies pass QC score >= ", min_qc_score)

  if (nrow(fail_ab) > 0)
    msg("  Removed (low QC): ",
        paste(fail_ab$Heatmap_Label, collapse = ", "))

  if (nrow(no_score) > 0)
    msg("  Note: ", nrow(no_score),
        " antibodies had no QC score — kept by default")

  # ── 2. Validation status filtering ───────────────────────────────────────
  caution_labels <- meta_qc$Heatmap_Label[
    !is.na(meta_qc$Validation_Status) & meta_qc$Validation_Status == "C"
  ]

  if (remove_caution && length(caution_labels) > 0) {
    msg("Removing ", length(caution_labels), " 'Caution' (C) antibodies")
    pass_ab <- pass_ab[pass_ab$Validation_Status != "C" |
                         is.na(pass_ab$Validation_Status), ]
  } else if (length(caution_labels) > 0) {
    msg(length(caution_labels),
        " antibodies have Validation_Status='C'. ",
        "Set remove_caution=TRUE to exclude them.")
  }

  # Heatmap_Label = the column names in the expression sheets
  keep_proteins <- c(pass_ab$Heatmap_Label, no_score$Heatmap_Label)
  keep_proteins <- keep_proteins[!is.na(keep_proteins)]

  # ── 3. Sample QC filtering ───────────────────────────────────────────────
  samp_qc   <- rppa$sample_qc
  fail_samp <- samp_qc[!is.na(samp_qc$Low_Protein) & samp_qc$Low_Protein, ]

  if (nrow(fail_samp) > 0) {
    msg("WARNING: ", nrow(fail_samp),
        " sample(s) have Total Protein < -3:")
    msg("  ", paste(fail_samp$Sample_Description, collapse = ", "))
    if (remove_low_protein) msg("  -> Removing these samples.")
  } else {
    msg("All ", nrow(samp_qc), " samples pass Total Protein QC \u2713")
  }

  keep_samples <- samp_qc$Sample_Name
  if (remove_low_protein && nrow(fail_samp) > 0)
    keep_samples <- setdiff(keep_samples, fail_samp$Sample_Name)

  # ── 4. Filter each expression sheet ─────────────────────────────────────
  filter_sheet <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(df)

    # sample metadata columns (non-protein columns)
    meta_cols <- c("Order", "Sample_Source", "Category_1", "Category_2",
                   "Category_3", "Sample", "Sample_Name", "Sample_Description",
                   "Sample_Type", "Timepoint", "Replicate")
    meta_cols <- intersect(meta_cols, colnames(df))
    prot_cols <- setdiff(colnames(df), meta_cols)

    # keep only proteins that passed QC (matched by Heatmap_Label = column name)
    keep_pc <- intersect(prot_cols, keep_proteins)

    # keep only samples that passed QC
    df_filt <- df[df$Sample_Name %in% keep_samples, c(meta_cols, keep_pc),
                  drop = FALSE]
    df_filt
  }

  rppa_filtered <- rppa
  rppa_filtered$l4_log2   <- filter_sheet(rppa$l4_log2)
  rppa_filtered$l4_linear <- filter_sheet(rppa$l4_linear)
  rppa_filtered$l4_chm    <- filter_sheet(rppa$l4_chm)
  rppa_filtered$l3_log2   <- filter_sheet(rppa$l3_log2)
  rppa_filtered$l2_log2   <- filter_sheet(rppa$l2_log2)

  # ── 5. QC report ─────────────────────────────────────────────────────────
  qc_report <- list(
    n_antibodies_input   = nrow(meta_qc),
    n_antibodies_kept    = length(keep_proteins),
    n_antibodies_removed = nrow(fail_ab),
    antibodies_removed   = fail_ab$Heatmap_Label,
    caution_antibodies   = caution_labels,
    n_samples_input      = nrow(samp_qc),
    n_samples_kept       = length(keep_samples),
    low_protein_samples  = fail_samp$Sample_Description,
    settings = list(
      min_qc_score       = min_qc_score,
      remove_caution     = remove_caution,
      remove_low_protein = remove_low_protein
    )
  )

  rppa_filtered$qc_report <- qc_report

  msg("QC complete. ",
      length(keep_proteins), " antibodies x ",
      length(keep_samples), " samples retained.")

  rppa_filtered
}


#' Print a formatted QC summary
#'
#' @param rppa An `rppa_data` object that has been through [run_qc()].
#' @export
print_qc_summary <- function(rppa) {
  if (is.null(rppa$qc_report)) {
    message("No QC report found. Run run_qc() first.")
    return(invisible(NULL))
  }
  r <- rppa$qc_report
  cat("==========================================\n")
  cat("  RPPAnalyzeR - QC Summary\n")
  cat("==========================================\n")
  cat(sprintf("  Antibodies  : %d input -> %d kept (%d removed)\n",
              r$n_antibodies_input, r$n_antibodies_kept, r$n_antibodies_removed))
  cat(sprintf("  Samples     : %d input -> %d kept\n",
              r$n_samples_input, r$n_samples_kept))
  if (length(r$low_protein_samples) > 0) {
    cat("  Low-protein samples:\n")
    cat("   ", paste(r$low_protein_samples, collapse = "\n    "), "\n")
  } else {
    cat("  All samples passed protein QC \u2713\n")
  }
  cat(sprintf("  Caution (C) antibodies: %d\n", length(r$caution_antibodies)))
  cat("==========================================\n")
  invisible(rppa)
}
