# RPPAnalyzeR <img src="man/figures/logo.png" align="right" height="140" alt="RPPAnalyzeR logo"/>

> A complete R analysis pipeline for MD Anderson RPPA Core Excel output

**Author:** Amar Kumar  
**Version:** 0.1.0  
**License:** MIT  

<!-- badges: start -->
[![R-CMD-check](https://github.com/akumar901/RPPAnalyzeR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/akumar901/RPPAnalyzeR/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

> **Note:** The CI badge may show a build failure as the GitHub Actions environment does not have all RPPA dependencies pre-installed. The package installs and runs correctly locally following the instructions below.

---

## What is RPPAnalyzeR?

**RPPAnalyzeR** is an R package that takes the Excel workbook produced by the
[MD Anderson Functional Proteomics RPPA Core](https://www.mdanderson.org/research/research-resources/core-facilities/functional-proteomics-rppa-core.html)
and turns it into a complete, reproducible analysis — with one function call or
step by step.

It was built and fully tested on a **MCF7 human breast cancer cell serum-starvation
time-course experiment** (15 samples × 497 antibodies, Set208) but works
with any MD Anderson RPPA Core output file.

---

## Background — What is RPPA?

**Reverse Phase Protein Array (RPPA)** is a high-throughput antibody-based
proteomics platform. Lysates from your samples are printed onto nitrocellulose
slides and probed with hundreds of validated antibodies simultaneously.

The MD Anderson RPPA Core returns results as a multi-sheet Excel workbook
containing raw signal, normalised values at multiple levels, quality control
metrics, and pairwise comparison sheets.

### About this dataset

| Parameter | Value |
|---|---|
| Cell line | MCF7 human breast cancer cells |
| Experiment type | Serum starvation time-course |
| Condition | Cells switched from serum-containing → serum-free media |
| Time points | 0h (baseline), 2h, 4h, 8h, 24h |
| Replicates | 3 per time point |
| Total samples | 15 |
| Antibodies | 497 |
| RPPA Set | Set208 |
| Analyst | Amar Kumar |

**Why serum starvation?** Removing serum (growth factors) from cell culture
media activates major signalling pathways as cells respond to nutrient
deprivation. This is a classic model for studying PI3K/Akt/mTOR, MAPK/ERK,
and cell cycle regulation dynamics in real time.

---

## Understanding the Input File

The MD Anderson RPPA Core delivers a single `.xlsx` file containing 12 sheets.
Here is what each sheet contains and whether you need it:

| Sheet | Contents | Use |
|---|---|---|
| **L4 (log_2)** | Fully normalised log2 values — loading + batch corrected | ✅ PRIMARY — use for all analysis |
| **L4 (CHM)** | L4 + antibody median centering | ✅ Use for heatmaps |
| **L4 (linear)** | L4 converted to linear scale | ✅ Use for bar graphs only |
| **L4 CHM 0h vs 2h/4h/8h/24h** | Pairwise comparisons with Log2FC and pre-computed P values | ✅ Use for quick DE results |
| **Antibody QC Scores** | QC score per antibody (0–1 scale, min = 0.8) | ✅ Filter low-quality antibodies |
| **Sample QC metrics** | Total protein content per sample | ✅ Flag failed samples |
| **L3 (log_2)** | Loading-normalised only (intermediate step) | ⚠️ Intermediate only |
| **L2 (log_2)** | Raw RPPASPACE output, no normalisation | ⚠️ Raw only |
| **Report Software Versions** | Pipeline version info | ℹ️ Reference only |

### Why is the L4 sheet so disorganised?

The MD Anderson Excel file has a non-standard structure that confuses
standard Excel readers:

```
Row 1        : Completely empty
Rows 2-9     : Antibody metadata — OFFSET to columns 9-506 only
               (columns 1-8 are empty in these rows)
               Col 9  = metadata label (e.g. "Antibody Name")
               Cols 10-506 = values for each of the 497 antibodies
Row 10       : TRUE column headers (Order, Sample Name, ..., 497 protein names)
Rows 11-25   : 15 sample data rows (your actual expression values)
Row 26+      : May contain stray values (e.g. "RI" annotation in L4 linear)
               — these must be ignored
```

**RPPAnalyzeR handles this automatically.** It uses `skip` and `n_max`
arguments in `readxl` to read exactly the rows it needs and ignores everything
else.

### The antibody naming problem

The RPPA Core uses THREE different naming systems across sheets:

| Sheet | Name format | Example |
|---|---|---|
| Antibody QC Scores | Human-readable | `Cyclin B1` |
| L4 expression sheets | Heatmap label with suffix | `Cyclin-B1-R-V` |
| Metadata rows | Short antibody name | `Cyclin-B1` |

These **never match directly** — `"Cyclin B1" != "Cyclin-B1-R-V"`.
The only shared key across all sheets is the **Antigen ID**
(e.g. `AGID00024`), which is consistent everywhere.
RPPAnalyzeR joins all sheets via Antigen ID to ensure correct QC filtering.

### The antibody suffix key

Protein column names in L4 sheets follow the pattern:
`AntibodyName-Species-ValidationStatus`

| Code | Meaning |
|---|---|
| **R** | Rabbit antibody |
| **M** | Mouse antibody |
| **G** | Goat antibody |
| **V** | Validated — performs well in all RPPA assays |
| **C** | Use with Caution — mostly reliable, some edge cases |
| **Q** | Tissue-reactive — detects non-specific components |

Example: `Akt_pS473-R-V` = phospho-Akt Ser473, rabbit antibody, validated.

### The normalisation levels explained

```
L2 (log2)   Raw RPPASPACE curve-fitting output. No normalisation.
    |
    v
L3 (log2)   + Loading normalisation
              (bidirectional median centering: by antibody, then by sample)
    |
    v  
L4 (log2)   + Set-to-Set batch correction         <-- USE THIS for analysis
    |
    v
L4 (CHM)    + Antibody median centering            <-- USE THIS for heatmaps
```

Always use L4 for analysis. L2 and L3 are included for transparency only.

---

## Installation

### Prerequisites

- R version ≥ 4.0.0
- RStudio (recommended)

**Find your R version:** In RStudio Console type `R.version$version.string`

### Install dependencies

```r
install.packages(c(
  "readxl", "dplyr", "tidyr", "ggplot2", "pheatmap",
  "stringr", "tibble", "purrr", "openxlsx",
  "ggrepel", "RColorBrewer", "scales", "writexl",
  "remotes", "rmarkdown"
))
```

### Install RPPAnalyzeR from GitHub

```r
remotes::install_github("akumar901/RPPAnalyzeR")
```

### Install from local folder (without GitHub)

```r
install.packages(
  "/path/to/RPPAnalyzeR",
  repos = NULL,
  type  = "source"
)
```

---

## Usage

### Option 1 — Full pipeline in one line

```r
library(RPPAnalyzeR)

results <- run_pipeline(
  xlsx_path  = "your_rppa_file.xlsx",
  output_dir = "rppa_output"
)
```

### Option 2 — Step by step

```r
library(RPPAnalyzeR)

# 1. Import all 12 sheets
rppa <- import_rppa("your_rppa_file.xlsx")

# 2. Quality control
rppa <- run_qc(rppa, min_qc_score = 0.8)
print_qc_summary(rppa)

# 3. Differential expression — all time points vs 0h baseline
de_all <- diff_expression_all(rppa)

# Top 10 hits at 24h
top_proteins(de_all, n = 10, timepoint = "24h")

# 4. Plots
plot_sample_qc(rppa)                    # sample total protein bar chart
plot_antibody_qc(rppa)                  # QC score histogram
plot_volcano(de_all, timepoint = "24h") # volcano plot
plot_heatmap(rppa)                      # clustered heatmap top 50 proteins
plot_timecourse(rppa, proteins = c(     # time-course line plot
  "Akt_pS473-R-V",
  "p70-S6K_pT389-R-V",
  "MAPK_pT202_Y204-R-C"
))

# 5. Export
export_rppa_excel(rppa, de_results = de_all,
                  output_path = "results.xlsx")
export_de_csv(de_all, significant_only = TRUE)
```

### Option 3 — RMarkdown report

Open `RPPAnalyzeR_Analysis.Rmd`, update the two path lines at the top,
and click **Knit**. This generates a single self-contained HTML report
with all plots, tables, and results embedded.

---

## Known Issues and Solutions

These are real issues encountered during development — documented here
so other users don't have to debug them.

### Issue 1 — `Can't rename columns that don't exist`

**Cause:** `readxl` with `col_names = FALSE` assigns internal names
(`...1`, `...2`) to columns. Renaming by column name string then fails
because the original names no longer exist.

**Fix:** Assign column names by position using a hardcoded vector of
9 fixed metadata names, not by matching strings.

### Issue 2 — `Can't transform a data frame with duplicate names`

**Cause:** A stray value `"RI"` (Reference Interval annotation) sits
in row 30, column 316 of the `L4 (linear)` sheet — 5 rows below the
actual data. When all rows were read, this triggered column count
mismatches leading to recycled (duplicate) column names.

**Fix:** Use `n_max = 15` when reading expression sheets to read
exactly the 15 sample rows and nothing beyond.

### Issue 3 — `Column 1 must be named / Empty name found at location 1`

**Cause:** The metadata rows (rows 2–9) have columns 1–8 completely
empty. `readxl` detects the used range as starting at column 9,
returning a tibble where the label column becomes column 1 instead
of column 9. Transposing this matrix then produced empty column names.

**Fix:** Read each of the 8 metadata rows individually using
`skip = i, n_max = 1` and build the data frame directly with
hardcoded clean column names — no transpose needed.

### Issue 4 — `run_qc()` silently removes almost all proteins

**Cause:** The Antibody QC sheet uses human-readable names
(`"Cyclin B1"`) while the expression sheets use heatmap-label names
(`"Cyclin-B1-R-V"`). These share zero exact matches, so name-based
joining drops every protein.

**Fix:** Join QC scores to expression sheet proteins via **Antigen ID**
(e.g. `AGID00024`) — the only field that is identical across all sheets.
All 497 Antigen IDs match perfectly.

### Issue 5 — `None of the requested proteins found in expression matrix`

**Cause:** `plot_timecourse()` was building a lookup map using
suffix-stripped names as keys (e.g. `"Akt_pS473"`) but then looking
up full names (e.g. `"Akt_pS473-R-V"`) — a key mismatch.

**Fix:** `plot_timecourse()` now tries exact match first, then
falls back to suffix-stripped match, accepting protein names in
either format.

### Issue 6 — Old package version persists after reinstall

**Cause:** R keeps loaded packages in memory for the session even
after reinstalling. Knitting an RMD uses the in-memory version.

**Fix:** After reinstalling, always go to
**Session → Restart R** in RStudio before knitting.

### Issue 7 — Low-QC antibodies may appear in volcano plots and DE results

**Cause:** The MD Anderson CHM pairwise comparison sheets contain
pre-computed Log2FC and P values for **all 497 antibodies**, including
any that failed QC filtering (e.g. B7-H4, QC score < 0.8). The
`diff_expression()` function reads these pre-computed values directly
from the CHM sheets, meaning a QC-failed antibody could theoretically
appear as a significant hit in a volcano plot even though it was removed
from the expression matrix by `run_qc()`.

**Fix:** `diff_expression()` now cross-references the list of
QC-passing proteins from `rppa$l4_log2` and filters DE results to
exclude any antibody that did not pass QC. This ensures complete
consistency between the QC-filtered expression data and all downstream
DE analysis and visualisations. Download the fixed `03_differential.R`
from the repository and reinstall the package.

**Note:** In practice this only affects antibodies with a QC score
below 0.8. In the tested dataset (Set208), only B7-H4 was removed,
and its low QC score means it is very unlikely to produce a strong
DE signal. However the fix ensures rigorous consistency for all
future datasets.

### Issue 8 — Antibodies with "Q" validation code incorrectly excluded from DE results

**Cause:** The protein name matching regex used in `diff_expression()`
to cross-reference QC-passing proteins was `-[RMG]-[VCE]$`, which
handles Validated (V) and Caution (C) antibody suffixes but missed
the **Q (tissue-reactive)** suffix. This caused antibodies like
`Snail-M-Q` to be incorrectly excluded from DE results even though
they passed the QC score threshold, because `"Snail"` (CHM sheet name)
did not match `"Snail-M-Q"` (expression sheet name).

**Impact:** In the tested dataset (Set208), `Snail` was incorrectly
absent from the 2h and 8h DE results when the Issue 7 fix was first
applied.

**Fix:** The regex was updated to `-[RMG]-[VCEQq]$` to correctly
handle all three validation codes:
- **V** — Validated
- **C** — Use with Caution
- **Q** — Tissue-reactive

**Important note for users:** If you modify the protein name matching
regex in `diff_expression()`, ensure it covers all three validation
codes (V, C, Q). Missing any one code will silently exclude the
affected antibodies from DE results without any error message,
potentially affecting your significant hits list.

---

## Troubleshooting Tips

**Find your Mac username:**
```r
Sys.getenv("USER")
```

**Confirm the fix is in your installed version:**
```r
grep("stripped_map", deparse(body(plot_timecourse)))
# Should return line numbers — if character(0), old version is still loaded
```

**Find where R installed the package:**
```r
find.package("RPPAnalyzeR")
```

**Check what proteins are available after QC:**
```r
expr <- get_expression_matrix(rppa, "l4_log2")
head(colnames(expr$matrix), 20)
```

**Check protein names before plotting:**
```r
# Search for a protein
colnames(expr$matrix)[grep("Cyclin", colnames(expr$matrix))]
```

---

## Output Files

Running the full pipeline produces:

```
output_dir/
├── RPPAnalyzeR_results.xlsx        All key sheets + DE results (colour-coded)
├── DE_significant_hits.csv         Significant proteins across all time points
└── plots/
    ├── 01_sample_QC.png            Bar chart of total protein per sample
    ├── 02_antibody_QC.png          Antibody QC score histogram
    ├── 03_volcano_0h_vs_2h.png     Volcano plot: 0h vs 2h
    ├── 03_volcano_0h_vs_4h.png     Volcano plot: 0h vs 4h
    ├── 03_volcano_0h_vs_8h.png     Volcano plot: 0h vs 8h
    ├── 03_volcano_0h_vs_24h.png    Volcano plot: 0h vs 24h
    ├── 04_heatmap_top50.png        Clustered heatmap top 50 proteins
    ├── 05_timecourse_mTOR_PI3K.png mTOR and PI3K pathway time-course
    ├── 05_timecourse_MAPK_ERK.png  MAPK and ERK pathway time-course
    └── 05_timecourse_CellCycle.png Cell cycle and apoptosis time-course
```

---

## License

MIT © Amar Kumar — see [LICENSE](LICENSE) for details.

You are free to use, modify, share, and build on this package.
If you use it in a publication, a citation is appreciated.

---

## Citation

```
Kumar A (2025). RPPAnalyzeR: A complete R analysis pipeline for
MD Anderson RPPA Core output. R package version 0.1.0.
https://github.com/akumar901/RPPAnalyzeR
```

Please also acknowledge the MD Anderson RPPA Core:

> *The Functional Proteomics RPPA Core is supported by MD Anderson
> Cancer Center Support Grant # 5 P30 CA016672-40.*

---

## Contact

**Amar Kumar**  
Email: amarcompbio@gmail.com  
Issues and pull requests welcome:
[https://github.com/akumar901/RPPAnalyzeR/issues](https://github.com/akumar901/RPPAnalyzeR/issues)
