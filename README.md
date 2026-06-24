# RPPAnalyzeR <img src="man/figures/logo.png" align="right" height="100" alt="" />

> A complete R analysis pipeline for MD Anderson RPPA Core Excel output

**Author:** Amar Kumar  
**Version:** 0.1.0  
**License:** MIT

<!-- badges: start -->
[![R-CMD-check](https://github.com/AmarKumar/RPPAnalyzeR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AmarKumar/RPPAnalyzeR/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

---

## What does this package do?

**RPPAnalyzeR** turns the messy, multi-sheet Excel workbook from the [MD Anderson Functional Proteomics RPPA Core](https://www.mdanderson.org/research/research-resources/core-facilities/functional-proteomics-rppa-core.html) into a complete, reproducible analysis — in one function call or step by step.

It was built and tested on a **MCF7 serum-starvation time-course experiment** (15 samples × 497 antibodies, Set208) but works with any MD Anderson RPPA Core output file.

### What it handles

- Imports **all 12 sheets** from the RPPA Core Excel file (L2/L3/L4 log2, linear, CHM, pairwise time-point comparison sheets, both QC sheets)
- Parses the **disorganized 9-row metadata header** that MD Anderson uses — automatically
- **QC filtering** by antibody score threshold and sample total protein
- **Differential expression** across all time points vs. baseline (uses pre-computed p-values from CHM sheets or runs fresh t-tests)
- Publication-ready **volcano plots**, **clustered heatmaps**, **time-course plots**, and **QC charts**
- **Excel export** with coloured significant hits, all key sheets, and antibody metadata

---

## Installation

```r
# Install from GitHub (one-time setup)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("AmarKumar/RPPAnalyzeR")
```

Install dependencies first if needed:

```r
install.packages(c(
  "readxl", "dplyr", "tidyr", "ggplot2", "pheatmap",
  "stringr", "tibble", "purrr", "openxlsx",
  "ggrepel", "RColorBrewer", "scales", "writexl"
))
```

---

## Quickstart — one line

```r
library(RPPAnalyzeR)

results <- run_pipeline(
  xlsx_path  = "01_Jonathan_Coloff__Vipin_Rawat.xlsx",
  output_dir = "rppa_output"
)
```

This single call runs the entire pipeline and writes to `rppa_output/`:

```
rppa_output/
├── RPPAnalyzeR_results.xlsx        ← all key sheets + DE results, colour-coded
├── DE_significant_hits.csv         ← significant proteins across all time points
└── plots/
    ├── 01_sample_QC.png            ← bar chart of total protein per sample
    ├── 02_antibody_QC.png          ← QC score distribution histogram
    ├── 03_volcano_0h_vs_2h.png     ← volcano plot: 0h vs 2h
    ├── 03_volcano_0h_vs_4h.png     ← volcano plot: 0h vs 4h
    ├── 03_volcano_0h_vs_8h.png     ← volcano plot: 0h vs 8h
    ├── 03_volcano_0h_vs_24h.png    ← volcano plot: 0h vs 24h
    ├── 04_heatmap_top50.png        ← clustered heatmap, top 50 variable proteins
    └── 05_timecourse_key_proteins.png  ← time-course for top DE proteins
```

---

## Step-by-step usage

### Step 1 — Import all sheets

```r
library(RPPAnalyzeR)

rppa <- import_rppa("01_Jonathan_Coloff__Vipin_Rawat.xlsx")
names(rppa)
```

Returns a named list:

| Element | Contents |
|---|---|
| `rppa$l4_log2` | **Primary data** — fully normalised log2 (loading + batch corrected) |
| `rppa$l4_chm` | L4 antibody-median-centred — input for heatmaps |
| `rppa$l4_linear` | L4 converted to linear scale — for bar graphs |
| `rppa$l3_log2` | L3 loading-normalised only |
| `rppa$l2_log2` | L2 raw un-normalised |
| `rppa$chm_timepoints` | Named list of pairwise CHM sheets (Log2FC + p-values) |
| `rppa$antibody_qc` | QC score table for all 497 antibodies |
| `rppa$sample_qc` | Total protein and flags for all samples |
| `rppa$metadata` | Antibody gene names, species (R/M/G), validation (V/C/Q) |
| `rppa$sample_info` | Sample descriptions, timepoints, replicate IDs |

### Step 2 — Quality control

```r
rppa <- run_qc(
  rppa,
  min_qc_score       = 0.8,    # MD Anderson minimum — scores below this are removed
  remove_caution     = FALSE,  # keep C antibodies (flagged but not removed)
  remove_low_protein = TRUE    # remove samples with Total Protein log2 < -3
)

print_qc_summary(rppa)
```

### Step 3 — Differential expression

```r
# Single time point
de_24h <- diff_expression(rppa, timepoint = "24h")

# All time points at once
de_all <- diff_expression_all(rppa)

# Top hits at 24h
top_proteins(de_all, n = 10, timepoint = "24h")

# Only upregulated hits
top_proteins(de_all, n = 10, timepoint = "24h", direction = "Up")
```

### Step 4 — Plots

```r
# Sample QC bar chart
plot_sample_qc(rppa)

# Antibody QC score distribution
plot_antibody_qc(rppa)

# Volcano plot — 0h vs 24h, label top 15 proteins
plot_volcano(de_all, timepoint = "24h", label_top = 15)

# Clustered heatmap — top 50 most variable proteins
plot_heatmap(rppa)

# Custom protein set heatmap
plot_heatmap(rppa, proteins = c("Akt_pS473", "p70-S6K_pT389", "MAPK_pT202_Y204",
                                 "4E-BP1", "mTOR_pS2448", "S6_pS240_S244"))

# Time-course plot — mean ± SE across replicates
plot_timecourse(rppa, proteins = c(
  "p70-S6K_pT389",
  "Akt_pS473",
  "MAPK_pT202_Y204",
  "4E-BP1",
  "eEF2K"
))
```

### Step 5 — Export

```r
# Full Excel workbook — all key sheets, highlighted DE hits
export_rppa_excel(
  rppa,
  de_results  = de_all,
  output_path = "RPPAnalyzeR_results.xlsx"
)

# Significant hits only as CSV
export_de_csv(de_all, significant_only = TRUE,
              output_path = "significant_proteins.csv")

# Save any plot
p <- plot_volcano(de_all, timepoint = "24h")
save_plot(p, "volcano_24h.png", width = 8, height = 7)
```

---

## Understanding the MD Anderson normalisation levels

The package imports all four levels — always use **L4** for analysis:

```
L2 (log2)   Raw RPPASPACE output. No normalisation.
    ↓
L3 (log2)   + Loading normalisation (bidirectional median centering by antibody then by sample)
    ↓
L4 (log2)   + Set-to-Set batch correction                    ← USE THIS for analysis
    ↓
L4 (CHM)    + Antibody median centering                      ← USE THIS for heatmaps
```

L2 and L3 are included for transparency and troubleshooting only.

---

## Antibody label key

In heatmaps, antibody names follow: `AntibodyName-Species-ValidationStatus`

| Code | Meaning |
|---|---|
| R | Rabbit antibody |
| M | Mouse antibody |
| G | Goat antibody |
| V | Validated — performs well in all RPPA assays |
| C | Use with Caution — mostly reliable, some edge cases |
| Q | "Tissue reactive" — detects non-specific components in tissue samples |

Example: `Akt_pS473-R-V` = Akt phospho-Ser473, rabbit, validated.

---

## Why the L4 sheet looks disorganized

The MD Anderson Excel file has a specific structure:

- Rows 1 is blank
- Rows 2–9 contain antibody metadata (gene names, slide IDs, etc.) shifted to column 9 onwards
- Row 10 is the actual column header
- Rows 11 onwards are the 15 sample data rows

`RPPAnalyzeR` knows this format and parses it automatically. You never need to manually rearrange anything.

---

## Citing this package

If you use RPPAnalyzeR in your work, please also cite the MD Anderson RPPA Core:

> *The Functional Proteomics RPPA Core is supported by MD Anderson Cancer Center Support Grant # 5 P30 CA016672-40.*

---

## Contact

Amar Kumar — [your.email@example.com]

Issues and pull requests welcome: [https://github.com/AmarKumar/RPPAnalyzeR/issues](https://github.com/AmarKumar/RPPAnalyzeR/issues)
