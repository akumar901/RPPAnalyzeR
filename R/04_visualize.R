# =============================================================================
# RPPAnalyzeR :: 04_visualize.R
# Publication-ready plots: heatmap, volcano, time-course, QC bar chart
# =============================================================================

#' Plot a clustered heatmap from RPPA L4 CHM data
#'
#' @param rppa An `rppa_data` object.
#' @param proteins Character vector of protein names to include.
#'   If NULL, uses all proteins (may be slow for 497).
#' @param scale_rows Logical. Scale rows (proteins) to z-score. Default TRUE.
#' @param annotation_colors Named list for annotation colors.
#' @param title Character. Plot title.
#' @param filename Character. If provided, saves PNG to this path.
#' @param ... Additional arguments passed to pheatmap::pheatmap().
#'
#' @return Invisibly returns the pheatmap object.
#' @export
plot_heatmap <- function(rppa,
                         proteins         = NULL,
                         scale_rows       = TRUE,
                         annotation_colors = NULL,
                         title            = "RPPA L4 CHM Heatmap",
                         filename         = NULL,
                         ...) {

  if (!requireNamespace("pheatmap", quietly = TRUE))
    stop("Package 'pheatmap' is required. Install with: install.packages('pheatmap')")

  dat   <- rppa$l4_chm
  sinfo <- rppa$sample_info

  # Build expression matrix
  meta_cols <- c("Order", "Sample Source", "Category_1", "Category_2",
                 "Category_3", "Sample", "Sample_Name", "Sample_Description",
                 "Sample_Type", "Timepoint", "Replicate")
  meta_cols  <- intersect(meta_cols, colnames(dat))
  prot_cols  <- setdiff(colnames(dat), meta_cols)

  if (!is.null(proteins)) {
    prot_cols <- prot_cols[stringr::str_remove(prot_cols, "-[RMG]-[VCE]$") %in% proteins |
                           prot_cols %in% proteins]
  }

  mat <- t(as.matrix(dat[, prot_cols]))
  mode(mat) <- "numeric"
  colnames(mat) <- dat$Sample_Description

  # Sample annotation
  ann_col <- data.frame(
    Timepoint = dat$Timepoint,
    row.names  = dat$Sample_Description
  )

  # Timepoint color scale
  tps <- unique(dat$Timepoint)
  tp_colors <- setNames(
    RColorBrewer::brewer.pal(max(3, length(tps)), "Set2")[seq_along(tps)],
    tps
  )
  if (is.null(annotation_colors)) {
    annotation_colors <- list(Timepoint = tp_colors)
  }

  # Color scale: green–black–red (classic RPPA)
  heat_colors <- colorRampPalette(c("#1B7837", "black", "#C2182B"))(100)

  args <- list(
    mat                   = mat,
    scale                 = if (scale_rows) "row" else "none",
    cluster_rows          = TRUE,
    cluster_cols          = FALSE,   # keep time order
    annotation_col        = ann_col,
    annotation_colors     = annotation_colors,
    color                 = heat_colors,
    fontsize_row          = 7,
    fontsize_col          = 9,
    main                  = title,
    border_color          = NA,
    show_rownames         = nrow(mat) <= 60,
    filename              = if (!is.null(filename)) filename else NA
  )
  extra <- list(...)
  args[names(extra)] <- extra

  ph <- do.call(pheatmap::pheatmap, args)
  invisible(ph)
}


#' Volcano plot for a single time point
#'
#' @param de_result Output from [diff_expression()].
#' @param timepoint Character. Timepoint label (for title).
#' @param p_threshold Numeric. Significance line (default 0.05).
#' @param fc_threshold Numeric. FC line (default 0.5).
#' @param label_top Integer. Number of top proteins to label (default 15).
#' @param title Character. Plot title. Auto-generated if NULL.
#'
#' @return A ggplot object.
#' @export
plot_volcano <- function(de_result,
                         timepoint    = NULL,
                         p_threshold  = 0.05,
                         fc_threshold = 0.5,
                         label_top    = 15,
                         title        = NULL) {

  if (!requireNamespace("ggrepel", quietly = TRUE))
    stop("Package 'ggrepel' is required.")

  if (!is.null(timepoint)) de_result <- dplyr::filter(de_result, Timepoint == timepoint)

  de_result <- de_result %>%
    dplyr::mutate(
      neg_log10_p = -log10(Pvalue + 1e-10),
      Label       = dplyr::if_else(
        dplyr::row_number() <= label_top & Significant,
        Protein, ""
      ),
      Color = dplyr::case_when(
        Significant & Direction == "Up"   ~ "Upregulated",
        Significant & Direction == "Down" ~ "Downregulated",
        TRUE                              ~ "NS"
      )
    )

  tp_label <- if (!is.null(timepoint)) timepoint else unique(de_result$Timepoint)[1]
  plot_title <- title %||% paste0("Volcano Plot: 0h vs ", tp_label)

  ggplot2::ggplot(de_result,
                  ggplot2::aes(x = Log2FC, y = neg_log10_p,
                               color = Color, label = Label)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.8) +
    ggrepel::geom_text_repel(
      size = 3, max.overlaps = 20, segment.size = 0.3,
      box.padding = 0.4, show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = -log10(p_threshold),
                        linetype = "dashed", color = "grey40", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = c(-fc_threshold, fc_threshold),
                        linetype = "dashed", color = "grey40", linewidth = 0.5) +
    ggplot2::scale_color_manual(
      values = c("Upregulated" = "#C2182B", "Downregulated" = "#1B7837", "NS" = "grey60")
    ) +
    ggplot2::labs(
      title    = plot_title,
      subtitle = paste0("n = ", nrow(de_result), " proteins | ",
                        "Sig: |Log2FC| > ", fc_threshold,
                        " & P < ", p_threshold),
      x        = "Log2 Fold Change (treat - ctrl)",
      y        = "-log10(P value)",
      color    = NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      legend.position = "top",
      panel.grid.minor = ggplot2::element_blank()
    )
}


#' Time-course line plot for selected proteins
#'
#' Plots mean ± SE expression across all time points for a set of proteins.
#'
#' @param rppa An `rppa_data` object.
#' @param proteins Character vector. Protein names to plot.
#' @param sheet Character. Which expression sheet to use (default "l4_log2").
#' @param timepoint_order Character vector. Order of timepoints for x-axis.
#'
#' @return A ggplot object.
#' @export
plot_timecourse <- function(rppa,
                            proteins,
                            sheet           = "l4_log2",
                            timepoint_order = c("0h", "2h", "4h", "8h", "24h")) {

  expr  <- get_expression_matrix(rppa, sheet)
  mat   <- expr$matrix
  sinfo <- expr$sample_info

  # Match protein names — accepts EITHER full name ("Akt_pS473-R-V")
  # OR short name without suffix ("Akt_pS473"). Tries exact match first,
  # then falls back to suffix-stripped match.
  avail_prots  <- colnames(mat)
  stripped_map <- setNames(avail_prots,
                           stringr::str_remove(avail_prots, "-[RMG]-[VCE]$"))

  resolved <- sapply(proteins, function(p) {
    if (p %in% avail_prots) return(p)          # exact full name match
    hit <- stripped_map[p]
    if (!is.na(hit)) return(unname(hit))        # stripped name match
    return(NA_character_)
  })
  names(resolved) <- proteins
  resolved <- resolved[!is.na(resolved)]

  if (length(resolved) == 0)
    stop("None of the requested proteins found in expression matrix.\n",
         "Available (first 10): ", paste(head(avail_prots, 10), collapse = ", "))

  # Build long data frame
  long <- purrr::map_dfr(seq_along(resolved), function(i) {
    prot_col <- resolved[i]
    prot_nm  <- names(resolved)[i]
    vals     <- as.numeric(mat[, prot_col])
    tibble::tibble(
      Sample      = rownames(mat),
      Protein     = prot_nm,
      Expression  = vals,
      Timepoint   = sinfo$Timepoint
    )
  }) %>%
    dplyr::mutate(
      Timepoint = factor(Timepoint, levels = timepoint_order)
    ) %>%
    dplyr::filter(!is.na(Timepoint))

  # Summary stats
  summ <- long %>%
    dplyr::group_by(Protein, Timepoint) %>%
    dplyr::summarise(
      Mean = mean(Expression, na.rm = TRUE),
      SE   = sd(Expression,   na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    )

  ggplot2::ggplot(summ, ggplot2::aes(x = Timepoint, y = Mean,
                                     color = Protein, group = Protein)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = Mean - SE, ymax = Mean + SE),
      width = 0.2, linewidth = 0.6
    ) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(
      title    = "Protein Time-Course",
      subtitle = "Mean ± SE across replicates (L4 log2)",
      x        = "Time point",
      y        = "L4 log2 expression",
      color    = "Protein"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title     = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}


#' Bar chart of sample total protein content (QC)
#'
#' @param rppa An `rppa_data` object.
#' @return A ggplot object.
#' @export
plot_sample_qc <- function(rppa) {
  sq <- rppa$sample_qc %>%
    dplyr::mutate(
      Sample_Description = factor(Sample_Description,
                                  levels = Sample_Description),
      Flag = dplyr::if_else(Low_Protein, "Low (<-3)", "Pass")
    )

  ggplot2::ggplot(sq, ggplot2::aes(x = Sample_Description,
                                    y = Total_Protein,
                                    fill = Flag)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = -3, color = "red",
                        linetype = "dashed", linewidth = 0.8) +
    ggplot2::scale_fill_manual(values = c("Pass" = "#3182BD", "Low (<-3)" = "#E41A1C")) +
    ggplot2::labs(
      title    = "Sample QC: Total Protein Content",
      subtitle = "Red dashed line = threshold (log2 < -3); samples below are unreliable",
      x        = "Sample",
      y        = "Total Protein (log2)",
      fill     = "QC Status"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      plot.title    = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    )
}


#' Plot antibody QC score distribution
#'
#' @param rppa An `rppa_data` object.
#' @param threshold Numeric. Minimum QC threshold line (default 0.8).
#' @return A ggplot object.
#' @export
plot_antibody_qc <- function(rppa, threshold = 0.8) {
  aq <- rppa$antibody_qc %>% dplyr::filter(!is.na(QC_Score))

  ggplot2::ggplot(aq, ggplot2::aes(x = QC_Score)) +
    ggplot2::geom_histogram(binwidth = 0.01, fill = "#2171B5",
                            color = "white", alpha = 0.85) +
    ggplot2::geom_vline(xintercept = threshold, color = "red",
                        linetype = "dashed", linewidth = 0.9) +
    ggplot2::annotate("text", x = threshold - 0.005, y = Inf,
                      label = paste0("Threshold = ", threshold),
                      hjust = 1, vjust = 1.5, color = "red", size = 3.5) +
    ggplot2::labs(
      title    = "Antibody QC Score Distribution",
      subtitle = paste0("n = ", nrow(aq), " antibodies | ",
                        sum(aq$QC_Score >= threshold), " pass threshold"),
      x        = "QC Score (0 = worst, 1 = best)",
      y        = "Count"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}


# Null coalescing operator (if not available in older R)
`%||%` <- function(a, b) if (!is.null(a)) a else b
