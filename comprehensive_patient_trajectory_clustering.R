# =============================================================================
# COMPREHENSIVE PATIENT TRAJECTORY CLUSTERING WORKFLOW
# =============================================================================
# Author : Comprehensive Analysis Pipeline
# Created: 2024
# Purpose: Complete patient clustering workflow combining baseline clustering,
#          trajectory analysis, and data preparation utilities
# 
# This script combines functionality from:
# - Patient_Trajectory_Clustering.R (main analysis)
# - create_mrn_ngs.R (MRN/NGS extraction)
# - convert_sheet2_long.R (wide-to-long conversion)
# =============================================================================

# =============================================================================
# SECTION 1: SETUP AND CONFIGURATION
# =============================================================================

# Load Required Packages -----------------------------------------------------
required_packages <- c("data.table", "tidyverse", "pheatmap", "RColorBrewer", 
                      "cluster", "NbClust", "readxl", "lubridate", "ggplot2", 
                      "reshape2", "survival", "survminer")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# =============================================================================
# SECTION 2: UTILITY FUNCTIONS
# =============================================================================

# Function 1: Extract MRN and NGS Dates (from create_mrn_ngs.R)
extract_mrn_ngs <- function(excel_file, sheet_name = "Sheet2", output_file = "MRN_NGS_dates.csv") {
  message("Extracting MRN and NGS date columns...")
  
  df <- as.data.table(read_excel(excel_file, sheet = sheet_name, .name_repair = 'minimal'))
  
  mrn_idx <- which(tolower(names(df)) == 'mrn')
  ngs_idx <- which(tolower(names(df)) == 'ngs date')
  
  if (length(ngs_idx) == 0) {
    stop('No NGS date columns found')
  }
  
  sel_idx <- c(mrn_idx, ngs_idx)
  sub <- df[, ..sel_idx]
  
  # Rename columns
  if (length(mrn_idx) > 0) {
    setnames(sub, 1, 'MRN')
  } else {
    sub[, MRN := NA]
    setcolorder(sub, c('MRN', setdiff(names(sub), 'MRN')))
  }
  
  ngs_names <- paste0('NGS_date_', seq_len(ncol(sub) - 1))
  setnames(sub, 2:ncol(sub), ngs_names)
  
  fwrite(sub, output_file)
  message(paste('Saved', nrow(sub), 'rows and', ncol(sub), 'cols to', output_file))
  
  return(sub)
}

# Function 2: Convert Wide to Long Format (from convert_sheet2_long.R)
convert_wide_to_long <- function(excel_file, sheet_name = "Sheet2", output_file = "CMML_Serial_long.csv") {
  message(paste('Converting Sheet2 from', excel_file, 'to long format...'))
  
  dt <- as.data.table(read_excel(excel_file, sheet = sheet_name, .name_repair = 'minimal'))
  
  # Identify MRN and repeated NGS date columns
  mrn_idx <- which(tolower(names(dt)) == 'mrn')
  if (length(mrn_idx) == 0) stop('MRN column not found')
  
  ngs_idx <- which(tolower(names(dt)) == 'ngs date')
  if (length(ngs_idx) == 0) stop('No "NGS date" columns found')
  
  ngs_idx <- sort(ngs_idx)
  seg_end <- c(ngs_idx[-1] - 1, ncol(dt))
  
  segments <- vector('list', length(ngs_idx))
  
  for (i in seq_along(ngs_idx)) {
    col_range <- ngs_idx[i]:seg_end[i]
    seg <- dt[, c(mrn_idx, col_range), with = FALSE]
    setnames(seg, 1:2, c('MRN', 'NGS_date'))
    
    # Ensure unique gene column names
    gnam <- names(seg)
    if (length(gnam) > 2) {
      gnam[3:length(gnam)] <- make.unique(gnam[3:length(gnam)], sep = '_dup')
      setnames(seg, gnam)
    }
    
    # Drop rows with NA NGS_date
    seg <- seg[!is.na(NGS_date) & NGS_date != '']
    segments[[i]] <- seg
    message(paste('Segment', i, ':', nrow(seg), 'rows,', ncol(seg), 'cols'))
  }
  
  long_dt <- rbindlist(segments, use.names = TRUE, fill = TRUE)
  setcolorder(long_dt, c('MRN', 'NGS_date', setdiff(names(long_dt), c('MRN', 'NGS_date'))))
  
  fwrite(long_dt, output_file)
  message(paste('Saved long-form data to', output_file, '(', nrow(long_dt), 'rows)'))
  
  return(long_dt)
}

# Function 3: Normalize for Cosine Similarity
normalize_l2 <- function(x) {
  nrm <- sqrt(sum(x^2))
  if (nrm == 0) return(x)
  x / nrm
}

# =============================================================================
# SECTION 3: CONFIGURATION AND DATA LOADING
# =============================================================================

# User-Configurable Parameters -----------------------------------------------
clinical_excel <- "CMML Project_2.xlsx"   # workbook with MRN column
clinical_sheet <- "Sheet1"                # sheet containing baseline clinical data
serial_sheet   <- "Sheet2"                # sheet containing serial NGS data
mutation_file  <- "CMML_Sheet2.csv"       # wide mutation data (baseline genes)
output_prefix  <- "baseline"              # prefix for output files
set.seed(123)

message("=== STARTING COMPREHENSIVE PATIENT CLUSTERING ANALYSIS ===")

# -----------------------------
# Step 1 : Load and Prepare Data
# -----------------------------
message("Loading clinical data …")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
library(readxl)

if (file.exists(clinical_excel)) {
  message("Loading clinical data from Excel …")
  df_clinical <- read_excel(clinical_excel, sheet = clinical_sheet, .name_repair = "minimal")
  df_clinical <- as.data.table(df_clinical)
} else {
  stop(paste0("Clinical Excel file (", clinical_excel, ") not found."))
}

# Detect diagnosis / subtype column ------------------------------------------------
diag_candidates <- grep("diagnosis", names(df_clinical), ignore.case = TRUE, value = TRUE)
diag_col <- if (length(diag_candidates) > 0) diag_candidates[1] else NULL

# Prepare survival columns if present -----------------------------------------
if (all(c("OS months", "Survival status") %in% names(df_clinical))) {
  df_clinical$`OS months` <- suppressWarnings(as.numeric(df_clinical$`OS months`))
  # Flexible parsing of event indicator
  stat_col <- df_clinical$`Survival status`
  event_parsed <- rep(NA_integer_, length(stat_col))
  # Case 1: already numeric 0/1
  if (is.numeric(stat_col)) {
    event_parsed <- ifelse(is.na(stat_col), NA_integer_, ifelse(stat_col > 0, 1L, 0L))
  } else {
    stat_chr <- tolower(trimws(as.character(stat_col)))
    event_parsed[stat_chr %in% c("dead", "deceased", "expired", "1", "yes", "y")] <- 1L
    event_parsed[stat_chr %in% c("alive", "0", "no", "n", "living", "censored", "unknown", "na", "")] <- 0L
  }
  df_clinical$event <- event_parsed
}

# Progression-free survival columns
if ("Progression Free Survival Months" %in% names(df_clinical)) {
  df_clinical$`Progression Free Survival Months` <- suppressWarnings(as.numeric(df_clinical$`Progression Free Survival Months`))
  # Event indicator: try "Leukemia progression" or "Progression of MDS"
  if ("Leukemia progression" %in% names(df_clinical)) {
    df_clinical$pfs_event <- ifelse(tolower(trimws(df_clinical$`Leukemia progression`)) %in% c("yes", "y", "1", "true"), 1, 0)
  } else if ("Progression of MDS" %in% names(df_clinical)) {
    df_clinical$pfs_event <- ifelse(tolower(trimws(df_clinical$`Progression of MDS`)) %in% c("yes", "y", "1", "true"), 1, 0)
  }
}

# BMT dates column (Bone Marrow Transplant) ----------------------------------
if ("BMT dates" %in% names(df_clinical)) {
  message("Processing BMT dates for trajectory analysis filtering...")
  df_clinical$bmt_date_parsed <- df_clinical$`BMT dates`
  
  # Parse BMT dates (handle Excel serial numbers and various date formats)
  bmt_vec <- as.character(df_clinical$`BMT dates`)
  bmt_parsed <- rep(NA_real_, length(bmt_vec))
  
  # Check for Excel serial numbers (numeric format)
  is_serial <- grepl("^[0-9]+(\\.[0-9]+)?$", bmt_vec) & !is.na(bmt_vec) & bmt_vec != ""
  if (any(is_serial)) {
    bmt_parsed[is_serial] <- as.numeric(bmt_vec[is_serial])
  }
  
  # Convert Excel serials to dates
  bmt_dates <- as.Date(bmt_parsed, origin = "1899-12-30")
  
  # Parse remaining dates using lubridate for various formats
  if (any(!is_serial & !is.na(bmt_vec) & bmt_vec != "")) {
    non_serial_idx <- which(!is_serial & !is.na(bmt_vec) & bmt_vec != "")
    parsed_dates <- suppressWarnings(lubridate::parse_date_time(bmt_vec[non_serial_idx], 
                                                               orders = c("ymd", "mdy", "dmy", "ymd HMS", "mdy HMS", "dmy HMS")))
    bmt_dates[non_serial_idx] <- as.Date(parsed_dates)
  }
  
  # Handle placeholder Excel serial 0 (1899-12-30) as missing
  bmt_dates[bmt_dates == as.Date("1899-12-30")] <- NA
  
  df_clinical$bmt_date_parsed <- bmt_dates
  
  n_bmt_patients <- sum(!is.na(df_clinical$bmt_date_parsed))
  message(paste("Found", n_bmt_patients, "patients with BMT dates"))
} else {
  message("No 'BMT dates' column found - all NGS data will be included in trajectory analysis")
  df_clinical$bmt_date_parsed <- NA
}

message("Loading mutation data …")
df_genes_raw <- fread(mutation_file, header = FALSE, stringsAsFactors = FALSE)

# The first row contains headers → extract and clean
headers <- df_genes_raw[1, ] %>% unlist() %>% trimws()
df_genes      <- df_genes_raw[-1, ]
setnames(df_genes, headers)

# Align patient rows between clinical and mutation data -----------------------
# Assumption: rows correspond 1-to-1 across files. Warn otherwise.
if (nrow(df_genes) != nrow(df_clinical)) {
  warning("Row count mismatch between clinical and gene files – verify alignment!")
}

# -----------------------------
# Build Binary Mutation Matrix
# -----------------------------
message("Constructing binary mutation matrix …")
all_cols <- colnames(df_genes)

unique_genes <- sort(unique(sapply(all_cols, function(col) {
  gene <- strsplit(col, "\\.")[[1]][1]
  if (!startsWith(gene, "p.") && !startsWith(gene, "c.") && !(gene %in% c("NGS date"))) {
    return(trimws(gene))
  } else {
    return(NA)
  }
})))
unique_genes <- unique_genes[!is.na(unique_genes)]

mutation_matrix <- sapply(unique_genes, function(gene) {
  gene_cols <- grep(paste0("^", gene, "(\\.|$)"), colnames(df_genes), value = TRUE)
  mut_cols  <- gene_cols[!grepl("pos|VAF|NM|TYPE", gene_cols, ignore.case = TRUE)]
  if (length(mut_cols) == 0) return(rep(0, nrow(df_genes)))
  gene_dat <- as.data.frame(df_genes[, ..mut_cols])
  gene_dat[] <- lapply(gene_dat, function(x) as.numeric(gsub("[^0-9.]", "", x)))
  gene_dat[is.na(gene_dat)] <- 0
  apply(gene_dat, 1, max, na.rm = TRUE)
})
mutation_matrix <- as.data.frame(mutation_matrix)

# remove genes with no variance ------------------------------------------------
variances <- apply(mutation_matrix, 2, var, na.rm = TRUE)
keep_cols <- !(variances == 0 | is.na(variances))
mutation_matrix <- mutation_matrix[, keep_cols, drop = FALSE]

# -----------------------------
# Normalise for Cosine Similarity
# -----------------------------
normalize_l2 <- function(x) {
  nrm <- sqrt(sum(x^2))
  if (nrm == 0) return(x)
  x / nrm
}
mutation_mat_norm <- t(apply(as.matrix(mutation_matrix), 1, normalize_l2))

# Track which patients retain non-zero mutation vectors
row_norms <- rowSums(mutation_mat_norm^2)
kept_idx <- which(row_norms > 0)
mutation_mat_norm <- mutation_mat_norm[kept_idx, , drop = FALSE]

message(paste("Patients with ≥1 mutation:", nrow(mutation_mat_norm), "/", nrow(mutation_matrix)))

# -----------------------------
# Determine Optimal k via Silhouette
# -----------------------------
sil_width <- sapply(2:10, function(k) {
  fit <- kmeans(mutation_mat_norm, centers = k, nstart = 25)
  d   <- as.matrix(dist(mutation_mat_norm, method = "euclidean"))
  sil <- silhouette(fit$cluster, d)
  mean(sil[, 3])
})
opt_k <- which.max(sil_width) + 1

message(paste("Optimal k (baseline):", opt_k))

# -----------------------------
# Final K-means Clustering
# -----------------------------
final_fit <- kmeans(mutation_mat_norm, centers = opt_k, nstart = 25)

# Map cluster labels back to ALL patients (NA for all-zero rows) --------------
baseline_cluster <- rep(NA_integer_, nrow(mutation_matrix))
baseline_cluster[kept_idx] <- final_fit$cluster

# Append to clinical data -----------------------------------------------------
df_clinical$baseline_cluster <- baseline_cluster

# Save results ----------------------------------------------------------------
message("Saving baseline clustering output …")
fwrite(df_clinical, paste0("CMML_clinical_", output_prefix, "_clustered.csv"))

# Silhouette plot -------------------------------------------------------------
sil_df <- data.frame(k = 2:10, sil = sil_width)

png(paste0("silhouette_analysis_", output_prefix, ".png"), width = 1200, height = 700)
plot(sil_df$k, sil_df$sil, type = "b", pch = 19, col = "steelblue",
     main = "Silhouette Analysis – Baseline Clustering",
     xlab = "Number of Clusters (k)", ylab = "Mean Silhouette Width")
abline(v = opt_k, col = "red", lty = 2)
text(opt_k + 0.5, max(sil_width), labels = paste("Optimal k =", opt_k), col = "red")
dev.off()

# ---------------------------------------------------------------------------
# PCA plot of baseline mutation profiles -------------------------------------
# ---------------------------------------------------------------------------
if (nrow(mutation_mat_norm) >= 2) {
  pca_base <- prcomp(mutation_mat_norm, center = TRUE, scale. = FALSE)
  pca_df <- data.frame(PC1 = pca_base$x[, 1], PC2 = pca_base$x[, 2],
                       idx = kept_idx,
                       Cluster = factor(final_fit$cluster))
  # Attach diagnosis if available
  if (!is.null(diag_col)) {
    pca_df$Diagnosis <- df_clinical[[diag_col]][pca_df$idx]
  } else {
    pca_df$Diagnosis <- factor(final_fit$cluster)
  }
  # Remove outliers via Z-score (PC1 & PC2)
  z1 <- as.numeric(scale(pca_df$PC1))
  z2 <- as.numeric(scale(pca_df$PC2))
  keep_filt <- ifelse(is.na(z1) | is.na(z2), FALSE, abs(z1) <= 3 & abs(z2) <= 3)
  pca_df_filt <- pca_df[keep_filt, ]
  p_pca_base <- ggplot(pca_df_filt, aes(x = PC1, y = PC2, colour = Cluster, shape = factor(Diagnosis))) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 12) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background  = element_rect(fill = "white", colour = NA)) +
    labs(title = "PCA – Baseline Mutation Profiles", colour = "Baseline Cluster", shape = "Diagnosis")
  ggsave("pca_plot_baseline_clusters.png", p_pca_base, width = 8, height = 6, dpi = 300, bg = "white")
}

# -----------------------------
# Patient Heatmap by Cluster ---------------------------------------------------
row_annot <- data.frame(Cluster = baseline_cluster)
rownames(row_annot) <- 1:nrow(mutation_matrix)
row_annot$Cluster <- ifelse(is.na(row_annot$Cluster), "None", as.character(row_annot$Cluster))
row_annot$Cluster <- factor(row_annot$Cluster)
row_annot$Cluster <- droplevels(row_annot$Cluster)

# Load viridis for color palette if not present
if (!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis")
library(viridis)

# Get all unique actual annotation values (not just levels)
cluster_values <- sort(unique(as.character(row_annot$Cluster)))
cluster_palette <- viridis::viridis(length(cluster_values))
names(cluster_palette) <- cluster_values
annotation_colors <- list(
  Cluster = cluster_palette,
  Serial = c(Serial = "#1b9e77", `Not Serial` = "#d95f02")
)
missing_clusters <- setdiff(unique(as.character(row_annot$Cluster)), names(annotation_colors$Cluster))
if (length(missing_clusters) > 0) {
  annotation_colors$Cluster[missing_clusters] <- "#bdbdbd" # gray color for missing annotation
}

# Use viridis color scale, higher font size, and PDF output
pheatmap(
  as.matrix(mutation_matrix),
  color = viridis(100, option = "D", direction = 1, alpha = 1, begin = 0, end = 1),
  cluster_rows = FALSE,             # Don't cluster patients, only columns
  cluster_cols = TRUE,              # Cluster genes
  show_rownames = FALSE,
  show_colnames = TRUE,             # Show gene names
  fontsize = 10,                    # Larger font for better readability
  annotation_row = row_annot,       # Patient clusters
  annotation_colors = annotation_colors,
  legend = TRUE,
  legend_breaks = c(0, 1),          # Show scale from 0 (no mutation) to 1 (mutated)
  legend_labels = c("No Mutation", "Mutation"),
  main = "Baseline Mutation Matrix by Patient Cluster",
  filename = paste0("all_samples_", output_prefix, "_mutation_heatmap.pdf"), # Save as PDF
  width = max(8, min(20, ncol(mutation_matrix) * 0.3 + 2)),
  height = max(6, min(15, nrow(mutation_matrix) * 0.08 + 3))
)

message("Baseline patient clustering completed ✔️")

# =============================================================================
# Step 2 : LOAD SERIAL NGS DATA & BUILD PATIENT TRAJECTORIES
# =============================================================================

# Load serial NGS data --------------------------------------------------------
# Prefer long-form CSV produced by convert_sheet2_long.R (one row per MRN × date).
# If the CSV is absent, fall back to reading Sheet2 and converting on-the-fly.

serial_long_csv <- "CMML_Serial_long.csv"

if (file.exists(serial_long_csv)) {
  message("Loading serial NGS data from long CSV …")
  df_serial_raw <- fread(serial_long_csv, stringsAsFactors = FALSE)
} else if (file.exists(clinical_excel)) {
  message("Long CSV not found – reading Sheet2 and converting in memory …")
  wide_dt <- as.data.table(read_excel(clinical_excel, sheet = serial_sheet, .name_repair = "minimal"))
  mrn_idx  <- which(tolower(names(wide_dt)) == "mrn")
  ngs_idx  <- which(tolower(names(wide_dt)) == "ngs date")
  if (length(mrn_idx) == 0 || length(ngs_idx) == 0) stop("MRN or NGS date cols missing in Sheet2.")
  ngs_idx  <- sort(ngs_idx)
  seg_end  <- c(ngs_idx[-1] - 1, ncol(wide_dt))
  segments <- vector("list", length(ngs_idx))
  for (i in seq_along(ngs_idx)) {
    rng <- ngs_idx[i]:seg_end[i]
    seg <- wide_dt[, c(mrn_idx, rng), with = FALSE]
    setnames(seg, 1:2, c("MRN", "NGS_date"))
    if (ncol(seg) > 2) {
      gnam <- names(seg)
      gnam[3:length(gnam)] <- make.unique(gnam[3:length(gnam)], sep = "_dup")
      setnames(seg, gnam)
    }
    seg <- seg[!is.na(NGS_date) & NGS_date != ""]
    segments[[i]] <- seg
  }
  df_serial_raw <- rbindlist(segments, use.names = TRUE, fill = TRUE)
  fwrite(df_serial_raw, serial_long_csv)
  message("Converted wide Sheet2 to long CSV (", serial_long_csv, ")")
} else {
  stop("Neither long CSV nor Excel workbook found – cannot load serial NGS data.")
}

# Clean column names -----------------------------------------------------------
orig_names <- names(df_serial_raw)
clean_names <- gsub("[^A-Za-z0-9]", "", orig_names)   # keep alphanum only
clean_names <- make.unique(clean_names, sep = "_")
colnames(df_serial_raw) <- clean_names

# Identify key columns --------------------------------------------------------
# Prefer explicit MRN column for patient identifier
if ("MRN" %in% names(df_serial_raw)) {
  id_candidates <- "MRN"
} else {
  id_candidates <- names(df_serial_raw)[grepl("Patient|Subject|ID", names(df_serial_raw), ignore.case = TRUE)]
}

# detect date and VAF columns (after id detection)
date_candidates <- names(df_serial_raw)[grepl("NGS_date|NGSdate|NGSDate|Date", names(df_serial_raw), ignore.case = TRUE)]
vaf_candidates  <- names(df_serial_raw)[grepl("VAF", clean_names, ignore.case = TRUE)]

if (length(id_candidates) == 0) {
  warning("No obvious patient identifier column found – generating PatientID placeholder.")
  df_serial_raw$PatientID <- paste0("Patient_", 1:nrow(df_serial_raw))
  id_col <- "PatientID"
} else {
  id_col <- id_candidates[1]
}

if (length(date_candidates) == 0) {
  stop("No date column detected (e.g., containing 'NGSdate'). Cannot proceed with trajectories.")
} else {
  date_col <- date_candidates[1]
}

if (length(vaf_candidates) == 0) {
  stop("No VAF columns detected in serial data. Check sheet structure.")
}

message(paste("Using columns – id:", id_col, ", date:", date_col, ", VAF cols:", length(vaf_candidates)))

# Melt to long format ----------------------------------------------------------
dt_serial <- as.data.table(df_serial_raw[, c(id_col, date_col, vaf_candidates), with = FALSE])
setnames(dt_serial, c(id_col, date_col), c("Patient", "Date"))

long_serial <- melt(dt_serial,
                    id.vars   = c("Patient", "Date"),
                    variable.name = "Gene",
                    value.name   = "VAF",
                    na.rm = TRUE)

# Ensure long_serial is a data.table
long_serial <- as.data.table(long_serial)

# Cast data types --------------------------------------------------------------
if (!requireNamespace("lubridate", quietly = TRUE)) install.packages("lubridate")

# Parse Date column (supports Excel serials stored as char/numeric or common date strings)
long_serial[, Date := {
  vec <- as.character(Date)
  is_serial <- grepl("^[0-9]+(\\.[0-9]+)?$", vec)
  out <- rep(NA_real_, length(vec))
  if (any(is_serial)) {
    out[is_serial] <- as.numeric(vec[is_serial])
  }
  res_dates <- as.Date(out, origin = "1899-12-30")
  if (any(!is_serial)) {
    parsed <- suppressWarnings(lubridate::parse_date_time(vec[!is_serial], orders = c("ymd", "mdy", "dmy")))
    res_dates[!is_serial] <- as.Date(parsed)
  }
  res_dates
}]
long_serial[, VAF  := as.numeric(as.character(VAF))]

# Treat placeholder Excel serial 0 (1899-12-30) as missing
long_serial[Date == as.Date("1899-12-30"), Date := NA]

# Identify serial status BEFORE filtering -----------------------------------
serial_counts_all <- long_serial[!is.na(Date), uniqueN(Date), by = Patient]
serial_mrns <- serial_counts_all[V1 >= 2]$Patient

# Keep only rows with numeric VAF ------------------------------------------------
valid_dates_flag <- !is.na(long_serial$Date)
long_serial <- long_serial[!is.na(VAF) & !is.na(Date)]

# =============================================================================
# BMT FILTERING: Exclude post-BMT NGS data from trajectory analysis
# =============================================================================
if (any(!is.na(df_clinical$bmt_date_parsed))) {
  message("Applying BMT date filtering to exclude post-BMT NGS data...")
  
  # Create patient-BMT date mapping
  # Find the patient ID column used in serial data
  patient_id_candidates <- names(df_clinical)[grepl("MRN|Patient|Subject|ID", names(df_clinical), ignore.case = TRUE)]
  
  if (length(patient_id_candidates) > 0) {
    # Use the first matching patient ID column
    patient_id_col <- patient_id_candidates[1]
    bmt_subset <- df_clinical[!is.na(bmt_date_parsed)]
    bmt_mapping <- data.table(Patient = bmt_subset[[patient_id_col]], 
                             BMT_Date = bmt_subset$bmt_date_parsed)
  } else {
    # Fallback to row numbers if no obvious patient ID column
    bmt_mapping <- data.table(Patient = 1:nrow(df_clinical), 
                             BMT_Date = df_clinical$bmt_date_parsed)
    bmt_mapping <- bmt_mapping[!is.na(BMT_Date)]
  }
  
  if (nrow(bmt_mapping) > 0) {
    # Merge BMT dates with serial data
    long_serial_with_bmt <- merge(long_serial, bmt_mapping, by = "Patient", all.x = TRUE)
    
    # Count rows before filtering
    n_before <- nrow(long_serial)
    
    # Filter: keep only NGS dates that are before BMT date (or patients without BMT)
    long_serial <- long_serial_with_bmt[is.na(BMT_Date) | Date < BMT_Date]
    
    # Remove the BMT_Date column
    long_serial[, BMT_Date := NULL]
    
    n_after <- nrow(long_serial)
    n_filtered <- n_before - n_after
    
    message(paste("Filtered out", n_filtered, "post-BMT NGS data points"))
    message(paste("Retained", n_after, "pre-BMT NGS data points for trajectory analysis"))
    
    # Report per-patient filtering
    if (n_filtered > 0) {
      filtered_patients <- unique(long_serial_with_bmt[!is.na(BMT_Date) & Date >= BMT_Date]$Patient)
      message(paste("BMT filtering applied to", length(filtered_patients), "patients"))
    }
  } else {
    message("No valid BMT dates found for filtering")
  }
} else {
  message("No BMT dates available - including all NGS data in trajectory analysis")
}

# Identify patients with ≥2 timepoints -----------------------------------------
pat_time_counts <- long_serial[, .N, by = .(Patient, Date)]  # one row per patient-date
multi_tp_pats   <- pat_time_counts[, .N, by = Patient][N >= 2]$Patient

message(paste("Patients with ≥2 NGS timepoints (after BMT filtering):", length(multi_tp_pats)))

if (length(multi_tp_pats) == 0) {
  stop("No patients with multiple NGS timepoints after BMT filtering. Trajectory analysis skipped.")
}

long_serial <- long_serial[Patient %in% multi_tp_pats]

# Compute mean VAF per patient-date (aggregate across genes) -------------------
pat_date_vaf <- long_serial[, .(mean_vaf = mean(VAF, na.rm = TRUE)), by = .(Patient, Date)]

# ENHANCED: Identify new mutations and compute VAF slopes -------------------------
message("Computing VAF slopes for new mutations appearing over time...")

# Function to compute VAF slope for new mutations
compute_new_mutation_slopes <- function(patient_data) {
  # Sort by date
  setorder(patient_data, Date)
  
  if (nrow(patient_data) < 2) return(NULL)
  
  # Get unique dates
  dates <- unique(patient_data$Date)
  if (length(dates) < 2) return(NULL)
  
  # Identify baseline (first timepoint) mutations
  baseline_genes <- patient_data[Date == dates[1] & VAF > 0]$Gene
  
  # Track new mutations appearing at later timepoints
  new_mutation_slopes <- c()
  
  for (i in 2:length(dates)) {
    current_date <- dates[i]
    current_genes <- patient_data[Date == current_date & VAF > 0]$Gene
    
    # Find genes that weren't present at baseline
    new_genes <- setdiff(current_genes, baseline_genes)
    
    if (length(new_genes) > 0) {
      # For each new gene, compute slope from appearance to current timepoint
      for (gene in new_genes) {
        gene_data <- patient_data[Gene == gene & Date >= current_date]
        if (nrow(gene_data) >= 2) {
          # Compute slope using linear regression
          time_numeric <- as.numeric(gene_data$Date - min(gene_data$Date))
          if (length(unique(time_numeric)) > 1) {
            slope <- lm(VAF ~ time_numeric, data = gene_data)$coefficients[2]
            if (!is.na(slope)) {
              new_mutation_slopes <- c(new_mutation_slopes, slope)
            }
          }
        }
      }
    }
  }
  
  return(new_mutation_slopes)
}

# Enhanced function to get detailed new mutation information
get_new_mutation_details <- function(patient_data) {
  setorder(patient_data, Date)
  
  if (nrow(patient_data) < 2) return(NULL)
  
  dates <- unique(patient_data$Date)
  if (length(dates) < 2) return(NULL)
  
  baseline_genes <- patient_data[Date == dates[1] & VAF > 0]$Gene
  new_mutation_details <- data.table()
  
  for (i in 2:length(dates)) {
    current_date <- dates[i]
    current_genes <- patient_data[Date == current_date & VAF > 0]$Gene
    new_genes <- setdiff(current_genes, baseline_genes)
    
    if (length(new_genes) > 0) {
      for (gene in new_genes) {
        gene_data <- patient_data[Gene == gene & Date >= current_date]
        if (nrow(gene_data) >= 2) {
          time_numeric <- as.numeric(gene_data$Date - min(gene_data$Date))
          if (length(unique(time_numeric)) > 1) {
            slope <- lm(VAF ~ time_numeric, data = gene_data)$coefficients[2]
            if (!is.na(slope)) {
              new_mutation_details <- rbind(new_mutation_details,
                data.table(Gene = gene, 
                          First_Appearance_Date = current_date,
                          VAF_Slope = slope,
                          Initial_VAF = gene_data$VAF[1],
                          Latest_VAF = gene_data$VAF[nrow(gene_data)]))
            }
          }
        }
      }
    }
  }
  
  return(new_mutation_details)
}

# Get detailed new mutation information for all patients
new_mutation_detailed <- long_serial[, {
  details <- get_new_mutation_details(.SD)
  if (!is.null(details) && nrow(details) > 0) {
    details
  } else {
    data.table(Gene = character(0), First_Appearance_Date = as.Date(character(0)),
               VAF_Slope = numeric(0), Initial_VAF = numeric(0), Latest_VAF = numeric(0))
  }
}, by = Patient]

# Compute new mutation features for each patient
new_mutation_features <- long_serial[, {
  slopes <- compute_new_mutation_slopes(.SD)
  if (length(slopes) > 0) {
    list(
      n_new_mutations = as.numeric(length(slopes)),
      mean_new_mut_slope = as.numeric(mean(slopes, na.rm = TRUE)),
      max_new_mut_slope = as.numeric(max(slopes, na.rm = TRUE)),
      min_new_mut_slope = as.numeric(min(slopes, na.rm = TRUE)),
      sd_new_mut_slope = as.numeric(sd(slopes, na.rm = TRUE))
    )
  } else {
    list(
      n_new_mutations = 0.0,
      mean_new_mut_slope = 0.0,
      max_new_mut_slope = 0.0,
      min_new_mut_slope = 0.0,
      sd_new_mut_slope = 0.0
    )
  }
}, by = Patient]

# Enhanced trajectory features combining mean VAF and new mutation slopes -----------
traj_features <- pat_date_vaf[, {
  # sort dates
  setorder(.SD, Date)
  vafs <- mean_vaf
  if (length(vafs) < 2) {
    NULL  # skip patients with <2 pts (shouldn't happen)
  } else {
    time_idx <- seq_along(vafs)
    list(
      n_timepoints = length(vafs),
      mean_vaf     = mean(vafs, na.rm = TRUE),
      sd_vaf       = sd(vafs, na.rm = TRUE),
      slope_vaf    = {
        if (length(unique(time_idx)) > 1) {
          coef(lm(vafs ~ time_idx))[2]
        } else {
          0
        }
      },
      max_vaf      = max(vafs, na.rm = TRUE),
      min_vaf      = min(vafs, na.rm = TRUE),
      range_vaf    = max(vafs, na.rm = TRUE) - min(vafs, na.rm = TRUE)
    )
  }
}, by = Patient]

# Merge with new mutation features
traj_features <- merge(traj_features, new_mutation_features, by = "Patient", all.x = TRUE)

# Replace missing values with 0 for new mutation features
traj_features[is.na(n_new_mutations), n_new_mutations := 0]
traj_features[is.na(mean_new_mut_slope), mean_new_mut_slope := 0]
traj_features[is.na(max_new_mut_slope), max_new_mut_slope := 0]
traj_features[is.na(min_new_mut_slope), min_new_mut_slope := 0]
traj_features[is.na(sd_new_mut_slope), sd_new_mut_slope := 0]

message(paste("Trajectory features computed for", nrow(traj_features), "patients"))

# =============================================================================
# SECTION 4: TRAJECTORY CLUSTERING
# =============================================================================

# Prepare trajectory feature matrix for clustering
traj_matrix <- as.matrix(traj_features[, -1])  # exclude Patient column
rownames(traj_matrix) <- traj_features$Patient

# Remove features with zero variance
traj_variances <- apply(traj_matrix, 2, var, na.rm = TRUE)
keep_traj_cols <- !(traj_variances == 0 | is.na(traj_variances))
traj_matrix <- traj_matrix[, keep_traj_cols, drop = FALSE]

# Scale features for clustering
traj_matrix_scaled <- scale(traj_matrix)

# Handle any remaining NaN/Inf values
traj_matrix_scaled[is.nan(traj_matrix_scaled) | is.infinite(traj_matrix_scaled)] <- 0

message("Performing trajectory clustering...")

# Determine optimal k for trajectory clustering
if (nrow(traj_matrix_scaled) >= 3) {
  traj_sil_width <- sapply(2:min(8, nrow(traj_matrix_scaled)-1), function(k) {
    fit <- kmeans(traj_matrix_scaled, centers = k, nstart = 25)
    d <- as.matrix(dist(traj_matrix_scaled, method = "euclidean"))
    sil <- silhouette(fit$cluster, d)
    mean(sil[, 3])
  })
  
  traj_opt_k <- which.max(traj_sil_width) + 1
  message(paste("Optimal k (trajectory):", traj_opt_k))
  
  # Final trajectory clustering
  traj_final_fit <- kmeans(traj_matrix_scaled, centers = traj_opt_k, nstart = 25)
  
  # Add trajectory cluster to features
  traj_features$trajectory_cluster <- traj_final_fit$cluster
  
  # Save trajectory clustering results
  fwrite(traj_features, "trajectory_clustering_results.csv")
  
  # Plot trajectory silhouette analysis
  traj_sil_df <- data.frame(k = 2:length(traj_sil_width)+1, sil = traj_sil_width)
  
  png("trajectory_silhouette_analysis.png", width = 1200, height = 700)
  plot(traj_sil_df$k, traj_sil_df$sil, type = "b", pch = 19, col = "darkgreen",
       main = "Silhouette Analysis – Trajectory Clustering",
       xlab = "Number of Clusters (k)", ylab = "Mean Silhouette Width")
  abline(v = traj_opt_k, col = "red", lty = 2)
  text(traj_opt_k + 0.5, max(traj_sil_width), labels = paste("Optimal k =", traj_opt_k), col = "red")
  dev.off()
  
} else {
  warning("Insufficient patients for trajectory clustering")
  traj_features$trajectory_cluster <- 1
}

# =============================================================================
# SECTION 5: COMBINED ANALYSIS AND VISUALIZATION
# =============================================================================

# Add serial status to clinical data
serial_status_vec <- rep("Not Serial", nrow(df_clinical))
if ("MRN" %in% names(df_clinical)) {
  serial_indices <- which(df_clinical$MRN %in% serial_mrns)
} else {
  # Use row numbers as fallback
  patient_id_candidates <- names(df_clinical)[grepl("MRN|Patient|Subject|ID", names(df_clinical), ignore.case = TRUE)]
  if (length(patient_id_candidates) > 0) {
    serial_indices <- which(df_clinical[[patient_id_candidates[1]]] %in% serial_mrns)
  } else {
    serial_indices <- integer(0)
  }
}
serial_status_vec[serial_indices] <- "Serial"
df_clinical$serial_status <- serial_status_vec

# Update row annotation for heatmaps
row_annot$Serial <- serial_status_vec

# Generate per-cluster heatmaps with Serial annotation
unique_clusters <- sort(unique(baseline_cluster[!is.na(baseline_cluster)]))

for (cluster_id in unique_clusters) {
  cluster_patients <- which(baseline_cluster == cluster_id)
  
  if (length(cluster_patients) > 1) {
    cluster_matrix <- mutation_matrix[cluster_patients, , drop = FALSE]
    cluster_annot <- row_annot[cluster_patients, , drop = FALSE]
    
    # Remove genes with no variation in this cluster
    cluster_variances <- apply(cluster_matrix, 2, var, na.rm = TRUE)
    cluster_keep_cols <- !(cluster_variances == 0 | is.na(cluster_variances))
    cluster_matrix <- cluster_matrix[, cluster_keep_cols, drop = FALSE]
    
    if (ncol(cluster_matrix) > 0) {
      pheatmap(
        as.matrix(cluster_matrix),
        color = viridis(100, option = "D", direction = 1, alpha = 1, begin = 0, end = 1),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = FALSE,
        show_colnames = TRUE,
        fontsize = 8,
        annotation_row = cluster_annot,
        annotation_colors = annotation_colors,
        legend = TRUE,
        legend_breaks = c(0, 1),
        legend_labels = c("No Mutation", "Mutation"),
        main = paste("Cluster", cluster_id, "Mutation Matrix"),
        filename = paste0("cluster_", cluster_id, "_mutation_heatmap.pdf"),
        width = max(6, min(15, ncol(cluster_matrix) * 0.4 + 2)),
        height = max(4, min(12, nrow(cluster_matrix) * 0.15 + 2))
      )
    }
  }
}

# Save final clinical data with all annotations
fwrite(df_clinical, "CMML_clinical_final_annotated.csv")

# Save detailed new mutation analysis
if (nrow(new_mutation_detailed) > 0) {
  fwrite(new_mutation_detailed, "new_mutations_detailed_analysis.csv")
}

message("=== COMPREHENSIVE PATIENT CLUSTERING ANALYSIS COMPLETED ===")
message("Output files generated:")
message("- CMML_clinical_baseline_clustered.csv")
message("- CMML_clinical_final_annotated.csv")
message("- trajectory_clustering_results.csv")
message("- new_mutations_detailed_analysis.csv")
message("- silhouette_analysis_baseline.png")
message("- trajectory_silhouette_analysis.png")
message("- pca_plot_baseline_clusters.png")
message("- all_samples_baseline_mutation_heatmap.pdf")
message("- cluster_*_mutation_heatmap.pdf (per cluster)")