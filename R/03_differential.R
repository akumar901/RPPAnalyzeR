# =============================================================================
# RPPAnalyzeR :: 03_differential.R
# Differential expression analysis across time points
# =============================================================================

#' Extract the expression matrix from an RPPA data sheet
#'
#' Returns a numeric matrix (samples × proteins) from the L4_log2 sheet,
#' with timepoint and replicate annotations.
#'
#' @param rppa An `rppa_data` object.
#' @param sheet Which sheet to use: "l4_log2" (default), "l3_log2", "l4_chm", "l4_linear".
#' @return A list with `matrix` (numeric, samples × proteins) and `sample_info` data frame.
#' @export
get_expression_matrix <- function(rppa, sheet = "l4_log2") {
  dat <- rppa[[sheet]]
  meta_cols <- c("Order", "Sample Source", "Category_1", "Category_2",
                 "Category_3", "Sample", "Sample_Name", "Sample_Description",
                 "Sample_Type", "Timepoint", "Replicate")
  meta_cols <- intersect(meta_cols, colnames(dat))
  prot_cols <- setdiff(colnames(dat), meta_cols)

  mat <- as.matrix(dat[, prot_cols])
  mode(mat) <- "numeric"
  rownames(mat) <- dat$Sample_Name

  list(
    matrix      = mat,
    sample_info = dat[, meta_cols]
  )
}


#' Differential expression analysis: one time point vs. 0h baseline
#'
#' For each protein, runs a two-sample t-test (or uses the pre-computed
#' p-values from the CHM sheets) comparing a treatment time point to
#' the 0h control.
#'
#' @param rppa An `rppa_data` object.
#' @param timepoint Character. The treatment time point label to test, e.g. "2h", "4h".
#'   Must match values in `rppa$sample_info$Timepoint`.
#' @param ctrl_timepoint Character. Control time point (default "0h").
#' @param p_threshold Numeric. Significance threshold (default 0.05).
#' @param fc_threshold Numeric. |Log2FC| threshold (default 0.5).
#' @param use_chm_sheet Logical. Use pre-computed p-values from CHM sheets
#'   (more accurate for small n). Default TRUE.
#'
#' @return A tibble with columns: Protein, Mean_Ctrl, Mean_Treat, Log2FC,
#'   Pvalue, Padj (BH), Significant, Direction.
#' @export
diff_expression <- function(rppa,
                            timepoint      = "2h",
                            ctrl_timepoint = "0h",
                            p_threshold    = 0.05,
                            fc_threshold   = 0.5,
                            use_chm_sheet  = TRUE) {

  # ── Option A: use pre-computed CHM timepoint sheet ───────────────────────
  if (use_chm_sheet && !is.null(rppa$chm_timepoints)) {
    # Find matching CHM sheet
    tp_key <- names(rppa$chm_timepoints)
    match_idx <- grep(gsub("h", "", timepoint), tp_key, ignore.case = TRUE)[1]

    if (!is.na(match_idx)) {
      chm <- rppa$chm_timepoints[[match_idx]]
      result <- chm %>%
        dplyr::select(Protein, Mean_Ctrl, Mean_Treat, Log2FC, Pvalue) %>%
        dplyr::filter(!is.na(Pvalue)) %>%
        dplyr::mutate(
          Padj        = p.adjust(Pvalue, method = "BH"),
          Significant = Pvalue < p_threshold & abs(Log2FC) > fc_threshold,
          Direction   = dplyr::case_when(
            Significant & Log2FC > 0 ~ "Up",
            Significant & Log2FC < 0 ~ "Down",
            TRUE                     ~ "NS"
          ),
          Timepoint   = timepoint
        ) %>%
        dplyr::arrange(Pvalue)
      return(result)
    }
  }

  # ── Option B: compute from expression matrix ─────────────────────────────
  expr <- get_expression_matrix(rppa, "l4_log2")
  mat  <- expr$matrix
  sinfo <- expr$sample_info

  ctrl_idx  <- which(sinfo$Timepoint == ctrl_timepoint)
  treat_idx <- which(sinfo$Timepoint == timepoint)

  if (length(ctrl_idx) == 0)  stop("Control timepoint '", ctrl_timepoint, "' not found.")
  if (length(treat_idx) == 0) stop("Treatment timepoint '", timepoint, "' not found.")

  ctrl_mat  <- mat[ctrl_idx,  , drop = FALSE]
  treat_mat <- mat[treat_idx, , drop = FALSE]

  result <- purrr::map_dfr(colnames(mat), function(prot) {
    ctrl_vals  <- as.numeric(ctrl_mat[, prot])
    treat_vals <- as.numeric(treat_mat[, prot])
    if (all(is.na(ctrl_vals)) || all(is.na(treat_vals))) {
      return(tibble::tibble(Protein = prot, Mean_Ctrl = NA_real_,
                            Mean_Treat = NA_real_, Log2FC = NA_real_,
                            Pvalue = NA_real_))
    }
    ttest <- tryCatch(
      t.test(treat_vals, ctrl_vals, var.equal = FALSE),
      error = function(e) list(p.value = NA_real_)
    )
    tibble::tibble(
      Protein    = prot,
      Mean_Ctrl  = mean(ctrl_vals,  na.rm = TRUE),
      Mean_Treat = mean(treat_vals, na.rm = TRUE),
      Log2FC     = mean(treat_vals, na.rm = TRUE) - mean(ctrl_vals, na.rm = TRUE),
      Pvalue     = ttest$p.value
    )
  })

  result %>%
    dplyr::filter(!is.na(Pvalue)) %>%
    dplyr::mutate(
      Padj        = p.adjust(Pvalue, method = "BH"),
      Significant = Pvalue < p_threshold & abs(Log2FC) > fc_threshold,
      Direction   = dplyr::case_when(
        Significant & Log2FC > 0 ~ "Up",
        Significant & Log2FC < 0 ~ "Down",
        TRUE                     ~ "NS"
      ),
      Timepoint   = timepoint
    ) %>%
    dplyr::arrange(Pvalue)
}


#' Run differential expression for all time points vs. baseline
#'
#' Convenience wrapper that loops over all non-baseline time points and
#' returns a combined results table.
#'
#' @param rppa An `rppa_data` object.
#' @param ctrl_timepoint Character. Baseline time point label (default "0h").
#' @param p_threshold Numeric. Default 0.05.
#' @param fc_threshold Numeric. Default 0.5.
#'
#' @return A tibble combining results for all time points, with a Timepoint column.
#' @export
diff_expression_all <- function(rppa,
                                ctrl_timepoint = "0h",
                                p_threshold    = 0.05,
                                fc_threshold   = 0.5) {
  timepoints <- setdiff(unique(rppa$sample_info$Timepoint), ctrl_timepoint)
  message("[RPPAnalyzeR] Running DE for time points: ",
          paste(timepoints, collapse = ", "))

  purrr::map_dfr(timepoints, function(tp) {
    res <- tryCatch(
      diff_expression(rppa, timepoint = tp,
                      ctrl_timepoint = ctrl_timepoint,
                      p_threshold    = p_threshold,
                      fc_threshold   = fc_threshold),
      error = function(e) {
        message("  Warning: could not compute DE for ", tp, ": ", e$message)
        NULL
      }
    )
    res
  })
}


#' Get top N significant proteins at a given time point
#'
#' @param de_result Output from [diff_expression()] or [diff_expression_all()].
#' @param n Integer. Number of top proteins to return (default 20).
#' @param timepoint Character. Filter by timepoint (optional).
#' @param direction Character. "Up", "Down", or "both" (default).
#'
#' @return A filtered and sorted tibble.
#' @export
top_proteins <- function(de_result, n = 20, timepoint = NULL, direction = "both") {
  res <- de_result %>% dplyr::filter(Significant)

  if (!is.null(timepoint)) res <- dplyr::filter(res, Timepoint == timepoint)

  if (direction == "Up")   res <- dplyr::filter(res, Direction == "Up")
  if (direction == "Down") res <- dplyr::filter(res, Direction == "Down")

  res %>%
    dplyr::arrange(Pvalue) %>%
    dplyr::slice_head(n = n)
}
