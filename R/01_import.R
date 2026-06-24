# =============================================================================
# RPPAnalyzeR :: 01_import.R
# Tested against: 01_Jonathan_Coloff__Vipin_Rawat.xlsx
# Sheet structure verified row-by-row before writing this code.
# =============================================================================

#' Import the complete RPPA Excel workbook
#'
#' @param path Character. Full path to the .xlsx file from the RPPA Core.
#' @param verbose Logical. Print progress messages (default TRUE).
#' @return A named list: l4_log2, l4_linear, l4_chm, l3_log2, l2_log2,
#'   chm_timepoints, antibody_qc, sample_qc, metadata, sample_info.
#' @export
#' @examples
#' \dontrun{
#' rppa <- import_rppa("path/to/your_rppa_file.xlsx")
#' names(rppa)
#' }
import_rppa <- function(path, verbose = TRUE) {
  if (!file.exists(path)) stop("File not found: ", path)

  msg <- function(...) if (verbose) message("[RPPAnalyzeR] ", ...)

  msg("Reading workbook: ", basename(path))
  sheets <- readxl::excel_sheets(path)
  msg("Found ", length(sheets), " sheets: ", paste(sheets, collapse = ", "))

  msg("Importing L4 (log2) — primary analysis sheet ...")
  l4_log2   <- .read_expression_sheet(path, "L4 (log_2)")

  msg("Importing L4 (linear) ...")
  l4_linear <- .read_expression_sheet(path, "L4 (linear)")

  msg("Importing L4 (CHM) ...")
  l4_chm    <- .read_expression_sheet(path, "L4 (CHM)")

  msg("Importing L3 (log2) ...")
  l3_log2   <- .read_expression_sheet(path, "L3 (log_2)")

  msg("Importing L2 (log2) ...")
  l2_log2   <- .read_expression_sheet(path, "L2 (log_2)")

  tp_sheets <- grep("CHM.*heat map|CHM.*vs", sheets, value = TRUE)
  msg("Importing ", length(tp_sheets), " time-point CHM sheets ...")
  chm_timepoints <- lapply(tp_sheets, function(s) .read_chm_timepoint(path, s))
  names(chm_timepoints) <- .clean_tp_names(tp_sheets)

  msg("Importing Antibody QC Scores ...")
  antibody_qc <- .read_antibody_qc(path)

  msg("Importing Sample QC metrics ...")
  sample_qc <- .read_sample_qc(path)

  msg("Extracting antibody metadata ...")
  metadata <- .extract_metadata(path)

  sample_info <- l4_log2[, c("Sample_Name", "Sample_Description",
                              "Sample_Type", "Timepoint", "Replicate"),
                         drop = FALSE]

  msg("Done. All sheets imported successfully.")

  structure(
    list(
      l4_log2        = l4_log2,
      l4_linear      = l4_linear,
      l4_chm         = l4_chm,
      l3_log2        = l3_log2,
      l2_log2        = l2_log2,
      chm_timepoints = chm_timepoints,
      antibody_qc    = antibody_qc,
      sample_qc      = sample_qc,
      metadata       = metadata,
      sample_info    = sample_info
    ),
    class = "rppa_data",
    file  = path
  )
}

# =============================================================================
# Internal helpers — verified against actual file structure
# =============================================================================

# Expression sheets (L2/L3/L4):
#   Row 1        : empty
#   Rows 2-9     : metadata (cols 1-8 empty; col 9 = label; cols 10-506 = values)
#   Row 10       : column headers (Order, Sample Source, ..., 497 protein names)
#   Rows 11-25   : 15 sample data rows
#   Row 26+      : may contain stray values (e.g. 'RI' in L4 linear) — ignored

#' @keywords internal
.read_expression_sheet <- function(path, sheet_name) {

  # --- Step 1: get protein names from row 10, cols 10-506 ---
  hdr <- readxl::read_excel(
    path, sheet = sheet_name,
    skip = 9, n_max = 1,          # skip rows 1-9, read only row 10
    col_names = FALSE,
    .name_repair = "minimal"
  )
  # cols 1-9 are the 9 fixed metadata columns; proteins start at col 10
  protein_names <- as.character(unlist(hdr[1, 10:ncol(hdr)]))
  # fix any NA names
  bad <- is.na(protein_names) | trimws(protein_names) == "" | protein_names == "NA"
  protein_names[bad] <- paste0("Protein_", which(bad))
  # ensure unique
  protein_names <- make.unique(protein_names, sep = "_dup")

  fixed_names <- c("Order", "Sample_Source", "Category_1", "Category_2",
                   "Category_3", "Sample", "Sample_Name",
                   "Sample_Description", "Sample_Type")

  all_colnames <- c(fixed_names, protein_names)

  # --- Step 2: read exactly 15 data rows (rows 11-25) ---
  dat <- readxl::read_excel(
    path, sheet = sheet_name,
    skip = 10, n_max = 15,        # skip rows 1-10, read rows 11-25 only
    col_names = FALSE,
    .name_repair = "minimal"
  )

  # trim columns to match our names (handles any trailing extra cols)
  n_use <- min(ncol(dat), length(all_colnames))
  dat   <- dat[, 1:n_use, drop = FALSE]
  colnames(dat) <- all_colnames[1:n_use]

  # remove any fully empty rows
  dat <- dat[rowSums(!is.na(dat)) > 1, , drop = FALSE]

  # --- Step 3: type conversions and derived columns ---
  dat$Order             <- suppressWarnings(as.numeric(dat$Order))
  dat$Sample_Name       <- as.character(dat$Sample_Name)
  dat$Sample_Description <- as.character(dat$Sample_Description)
  dat$Sample_Type       <- as.character(dat$Sample_Type)

  # extract timepoint and replicate from Sample_Description
  # format is e.g. "0h-1- wer" or "24h-3-wo-ser"
  dat$Timepoint <- trimws(gsub("-.*", "", dat$Sample_Description))
  dat$Replicate <- gsub("^[^-]+-([0-9]+).*", "\\1", dat$Sample_Description)

  # convert protein columns to numeric individually
  prot_present <- intersect(protein_names, colnames(dat))
  for (pc in prot_present) {
    dat[[pc]] <- suppressWarnings(as.numeric(dat[[pc]]))
  }

  tibble::as_tibble(dat)
}


# CHM timepoint sheets:
#   Row 1  : headers — col 1=None, cols 2-4=ctrl sample names,
#             cols 5-7=treat sample names, cols 8-10=None (FC formula headers),
#             col 11='Log 2 FC', col 12='P value'
#   Rows 2+ : protein data — cols 8-10 are Excel formulas (come through as NA
#             in readxl since they reference other cells); cols 11-12 are
#             pre-computed values

#' @keywords internal
.read_chm_timepoint <- function(path, sheet_name) {

  dat <- readxl::read_excel(
    path, sheet = sheet_name,
    skip = 1,               # skip the header row (row 1)
    col_names = FALSE,
    .name_repair = "minimal"
  )

  # Confirmed column structure (12 cols):
  # 1=Protein, 2-4=Ctrl values, 5-7=Treat values,
  # 8-10=FC formulas (NA in readxl), 11=Log2FC, 12=Pvalue
  colnames(dat) <- c("Protein",
                     "Ctrl_1", "Ctrl_2", "Ctrl_3",
                     "Treat_1", "Treat_2", "Treat_3",
                     "FC_1", "FC_2", "FC_3",
                     "Log2FC", "Pvalue")

  dat <- dat %>%
    dplyr::filter(!is.na(Protein), as.character(Protein) != "NA") %>%
    dplyr::mutate(
      dplyr::across(c(Ctrl_1, Ctrl_2, Ctrl_3,
                      Treat_1, Treat_2, Treat_3,
                      Log2FC, Pvalue),
                    ~ suppressWarnings(as.numeric(.x))),
      Mean_Ctrl   = rowMeans(cbind(Ctrl_1, Ctrl_2, Ctrl_3),  na.rm = TRUE),
      Mean_Treat  = rowMeans(cbind(Treat_1, Treat_2, Treat_3), na.rm = TRUE),
      # Use pre-computed Log2FC if available; recalculate if NA
      Log2FC      = dplyr::if_else(is.na(Log2FC),
                                   Mean_Treat - Mean_Ctrl, Log2FC),
      Significant = !is.na(Pvalue) & Pvalue < 0.05 & abs(Log2FC) > 0.5,
      Direction   = dplyr::case_when(
        Significant & Log2FC > 0 ~ "Up",
        Significant & Log2FC < 0 ~ "Down",
        TRUE                     ~ "NS"
      )
    ) %>%
    dplyr::select(Protein, Ctrl_1, Ctrl_2, Ctrl_3,
                  Treat_1, Treat_2, Treat_3,
                  Mean_Ctrl, Mean_Treat, Log2FC, Pvalue,
                  Significant, Direction) %>%
    dplyr::arrange(dplyr::desc(abs(Log2FC)))

  dat
}


# Antibody QC sheet:
#   Rows 1-2 : empty
#   Row 3    : header (Row | Antibody Lab ID | Antibody Name | Antigen ID | QC Score)
#   Rows 4+  : 497 antibody rows

#' @keywords internal
.read_antibody_qc <- function(path) {

  dat <- readxl::read_excel(
    path, sheet = "Antibody QC Scores",
    skip = 3, n_max = 497,         # skip 3 header rows, read 497 antibodies
    col_names = FALSE,
    .name_repair = "minimal"
  )
  colnames(dat) <- c("Row", "Lab_ID", "Antibody_Name", "Antigen_ID", "QC_Score_Raw")

  dat <- dat %>%
    dplyr::filter(!is.na(Row)) %>%
    dplyr::mutate(
      Row      = suppressWarnings(as.numeric(Row)),
      Lab_ID   = suppressWarnings(as.numeric(Lab_ID)),
      # QC score is formatted as "0.886 (0.8)" — extract the first number
      QC_Score = suppressWarnings(as.numeric(
        gsub("^([0-9.]+).*", "\\1", trimws(as.character(QC_Score_Raw)))
      ))
    ) %>%
    dplyr::filter(!is.na(Row))

  tibble::as_tibble(dat)
}


# Sample QC sheet:
#   Rows 1-2 : empty
#   Row 3    : header
#   Rows 4+  : 15 sample rows

#' @keywords internal
.read_sample_qc <- function(path) {

  dat <- readxl::read_excel(
    path, sheet = "Sample QC metrics",
    skip = 3, n_max = 15,          # skip 3 header rows, read 15 samples
    col_names = FALSE,
    .name_repair = "minimal"
  )
  colnames(dat) <- c("Order", "Sample_Source", "Category_1", "Category_2",
                     "Category_3", "Sample", "Sample_Name",
                     "Sample_Description", "Sample_Type", "Total_Protein")

  dat <- dat %>%
    dplyr::filter(!is.na(Order)) %>%
    dplyr::mutate(
      Order         = suppressWarnings(as.numeric(Order)),
      Total_Protein = suppressWarnings(as.numeric(Total_Protein)),
      Low_Protein   = !is.na(Total_Protein) & Total_Protein < -3,
      Timepoint     = trimws(gsub("-.*", "", as.character(Sample_Description)))
    )

  tibble::as_tibble(dat)
}


# Metadata: rows 2-9 of L4 (log_2)
#   Each row: cols 1-8 empty | col 9 = label | cols 10-506 = 497 protein values
#   We read each row individually using skip + n_max=1

#' @keywords internal
.extract_metadata <- function(path) {

  # fixed clean names for the 8 metadata rows (in order)
  known_labels <- c("Antibody_Name", "Antibody_Origin", "Gene_Name",
                    "Validation_Status", "Slide_ID", "Heatmap_Label",
                    "Antigen_ID", "Slide_Number")

  # get protein count from header row
  hdr <- readxl::read_excel(
    path, sheet = "L4 (log_2)",
    skip = 9, n_max = 1,
    col_names = FALSE, .name_repair = "minimal"
  )
  n_proteins <- ncol(hdr) - 9   # first 9 cols are sample metadata

  # read each of the 8 metadata rows one at a time
  # skip=1 -> reads file row 2; skip=2 -> reads file row 3; etc.
  meta_list <- lapply(seq_len(8), function(i) {
    row_raw <- readxl::read_excel(
      path, sheet = "L4 (log_2)",
      skip = i, n_max = 1,
      col_names = FALSE, .name_repair = "minimal"
    )
    # find the label column (first non-NA) then take everything after it
    row_vec  <- as.character(unlist(row_raw[1, ]))
    non_na   <- which(!is.na(row_vec) & row_vec != "NA")
    if (length(non_na) == 0) return(rep(NA_character_, n_proteins))
    label_col <- non_na[1]
    vals <- row_vec[(label_col + 1):length(row_vec)]
    # pad or trim to exactly n_proteins
    length(vals) <- n_proteins
    vals
  })

  # build data frame — one col per metadata field, one row per protein
  meta <- as.data.frame(
    stats::setNames(meta_list, known_labels),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- NULL

  tibble::as_tibble(meta)
}


#' @keywords internal
.clean_tp_names <- function(sheet_names) {
  stringr::str_to_lower(
    make.unique(
      gsub("_+", "_",
           gsub("[^a-zA-Z0-9]", "_",
                gsub("L4 CHM for heat map\\s*|L4 CHM\\s*", "",
                     sheet_names))),
      sep = "_"
    )
  )
}
