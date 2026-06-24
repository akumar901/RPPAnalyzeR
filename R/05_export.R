# =============================================================================
# RPPAnalyzeR :: 05_export.R
# Export important sheets and analysis results to Excel and CSV
# =============================================================================

#' Export all important RPPA sheets to a tidy Excel workbook
#'
#' Writes the most information-rich sheets to a single Excel file:
#'   - L4_log2 (primary data)
#'   - L4_CHM (heatmap data)
#'   - L4_linear (for bar graphs)
#'   - All time-point CHM comparisons (with Log2FC and P values)
#'   - Antibody QC scores
#'   - Sample QC metrics
#'   - Antibody metadata
#'   - DE summary (if provided)
#'
#' @param rppa An `rppa_data` object (optionally after [run_qc()]).
#' @param de_results Output from [diff_expression_all()] (optional).
#' @param output_path Character. Output .xlsx file path.
#' @param overwrite Logical. Overwrite if file exists (default TRUE).
#'
#' @return Invisibly returns the output path.
#' @export
export_rppa_excel <- function(rppa,
                               de_results  = NULL,
                               output_path = "RPPAnalyzeR_output.xlsx",
                               overwrite   = TRUE) {

  if (!requireNamespace("openxlsx", quietly = TRUE))
    stop("Package 'openxlsx' is required.")

  if (file.exists(output_path) && !overwrite)
    stop("File exists. Set overwrite=TRUE to replace.")

  wb <- openxlsx::createWorkbook()

  # ── Styles ────────────────────────────────────────────────────────────────
  header_style <- openxlsx::createStyle(
    fontColour = "#FFFFFF", bgFill = "#2171B5",
    textDecoration = "bold", halign = "CENTER", border = "Bottom"
  )
  sig_up_style <- openxlsx::createStyle(bgFill = "#FFCCCC")  # up = light red
  sig_dn_style <- openxlsx::createStyle(bgFill = "#CCFFCC")  # down = light green
  flag_style   <- openxlsx::createStyle(bgFill = "#FFD700", fontColour = "#7B3F00")

  .write_sheet <- function(wb, sheet_name, data, header_row = TRUE) {
    openxlsx::addWorksheet(wb, sheetName = sheet_name)
    openxlsx::writeData(wb, sheet = sheet_name, x = data,
                        headerStyle = header_style, borders = "rows",
                        borderColour = "#D0D0D0")
    openxlsx::setColWidths(wb, sheet = sheet_name,
                           cols = seq_len(ncol(data)),
                           widths = "auto")
    openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  }

  # ── 1. L4 log2 ────────────────────────────────────────────────────────────
  message("[RPPAnalyzeR::Export] Writing L4 (log2) ...")
  .write_sheet(wb, "L4_log2", rppa$l4_log2)

  # ── 2. L4 CHM ─────────────────────────────────────────────────────────────
  message("[RPPAnalyzeR::Export] Writing L4 CHM ...")
  .write_sheet(wb, "L4_CHM", rppa$l4_chm)

  # ── 3. L4 linear ──────────────────────────────────────────────────────────
  message("[RPPAnalyzeR::Export] Writing L4 linear ...")
  .write_sheet(wb, "L4_linear", rppa$l4_linear)

  # ── 4. Time-point CHM comparisons ─────────────────────────────────────────
  for (tp_name in names(rppa$chm_timepoints)) {
    sheet_label <- paste0("CHM_", tp_name)
    if (nchar(sheet_label) > 31) sheet_label <- substr(sheet_label, 1, 31)
    message("[RPPAnalyzeR::Export] Writing ", sheet_label, " ...")
    chm_data <- rppa$chm_timepoints[[tp_name]]
    .write_sheet(wb, sheet_label, chm_data)

    # Highlight significant proteins
    if ("Significant" %in% colnames(chm_data) &&
        "Direction" %in% colnames(chm_data)) {
      sig_up_rows <- which(chm_data$Direction == "Up") + 1
      sig_dn_rows <- which(chm_data$Direction == "Down") + 1
      if (length(sig_up_rows) > 0)
        openxlsx::addStyle(wb, sheet_label, sig_up_style,
                           rows = sig_up_rows, cols = 1:ncol(chm_data),
                           gridExpand = TRUE)
      if (length(sig_dn_rows) > 0)
        openxlsx::addStyle(wb, sheet_label, sig_dn_style,
                           rows = sig_dn_rows, cols = 1:ncol(chm_data),
                           gridExpand = TRUE)
    }
  }

  # ── 5. Antibody QC ────────────────────────────────────────────────────────
  message("[RPPAnalyzeR::Export] Writing Antibody QC ...")
  .write_sheet(wb, "Antibody_QC", rppa$antibody_qc)

  # ── 6. Sample QC ─────────────────────────────────────────────────────────
  message("[RPPAnalyzeR::Export] Writing Sample QC ...")
  .write_sheet(wb, "Sample_QC", rppa$sample_qc)

  # Highlight low-protein samples
  low_rows <- which(rppa$sample_qc$Low_Protein) + 1
  if (length(low_rows) > 0)
    openxlsx::addStyle(wb, "Sample_QC", flag_style,
                       rows = low_rows, cols = 1:ncol(rppa$sample_qc),
                       gridExpand = TRUE)

  # ── 7. Antibody metadata ─────────────────────────────────────────────────
  message("[RPPAnalyzeR::Export] Writing Antibody Metadata ...")
  .write_sheet(wb, "Antibody_Metadata", rppa$metadata)

  # ── 8. DE results (if provided) ───────────────────────────────────────────
  if (!is.null(de_results)) {
    message("[RPPAnalyzeR::Export] Writing DE results ...")
    .write_sheet(wb, "DE_All_Timepoints", de_results)
    # Colour significant rows
    sig_up_rows <- which(de_results$Direction == "Up") + 1
    sig_dn_rows <- which(de_results$Direction == "Down") + 1
    if (length(sig_up_rows) > 0)
      openxlsx::addStyle(wb, "DE_All_Timepoints", sig_up_style,
                         rows = sig_up_rows, cols = 1:ncol(de_results),
                         gridExpand = TRUE)
    if (length(sig_dn_rows) > 0)
      openxlsx::addStyle(wb, "DE_All_Timepoints", sig_dn_style,
                         rows = sig_dn_rows, cols = 1:ncol(de_results),
                         gridExpand = TRUE)
  }

  # ── Save ─────────────────────────────────────────────────────────────────
  openxlsx::saveWorkbook(wb, output_path, overwrite = overwrite)
  message("[RPPAnalyzeR::Export] Saved to: ", output_path)
  invisible(output_path)
}


#' Export DE results to a CSV
#'
#' @param de_results Output from [diff_expression()] or [diff_expression_all()].
#' @param output_path Character. Path for CSV output.
#' @param significant_only Logical. Write only significant hits. Default FALSE.
#' @export
export_de_csv <- function(de_results,
                           output_path      = "DE_results.csv",
                           significant_only = FALSE) {
  dat <- de_results
  if (significant_only) dat <- dplyr::filter(dat, Significant)
  readr_check <- requireNamespace("readr", quietly = TRUE)
  if (readr_check) {
    readr::write_csv(dat, output_path)
  } else {
    write.csv(dat, output_path, row.names = FALSE)
  }
  message("[RPPAnalyzeR::Export] DE CSV saved to: ", output_path,
          " (", nrow(dat), " rows)")
  invisible(output_path)
}


#' Save a ggplot to file
#'
#' @param plot A ggplot object.
#' @param filename Character. Output file path (.png, .pdf, .svg).
#' @param width Numeric. Width in inches (default 10).
#' @param height Numeric. Height in inches (default 7).
#' @param dpi Integer. Resolution for raster outputs (default 300).
#' @export
save_plot <- function(plot, filename, width = 10, height = 7, dpi = 300) {
  ggplot2::ggsave(filename, plot = plot, width = width,
                  height = height, dpi = dpi)
  message("[RPPAnalyzeR::Export] Plot saved: ", filename)
  invisible(filename)
}
