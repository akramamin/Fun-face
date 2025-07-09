# Comprehensive Patient Trajectory Clustering Workflow

This R script provides a complete analysis pipeline for patient trajectory clustering in CMML (Chronic Myelomonocytic Leukemia) research, combining baseline clustering, trajectory analysis, and data preparation utilities.

## Overview

The script combines functionality from multiple sources:
- Patient_Trajectory_Clustering.R (main analysis)
- create_mrn_ngs.R (MRN/NGS extraction)
- convert_sheet2_long.R (wide-to-long conversion)

## Features

### 1. Data Preparation
- **MRN and NGS Date Extraction**: Extract Medical Record Numbers and NGS sequencing dates
- **Wide-to-Long Format Conversion**: Convert serial NGS data from wide to long format
- **BMT Date Filtering**: Exclude post-bone marrow transplant NGS data from trajectory analysis
- **Data Alignment**: Ensure proper alignment between clinical and mutation data

### 2. Baseline Clustering
- **Binary Mutation Matrix Construction**: Build mutation presence/absence matrix
- **Cosine Similarity Normalization**: Normalize mutation vectors for clustering
- **Optimal k Selection**: Use silhouette analysis to determine optimal number of clusters
- **K-means Clustering**: Perform baseline mutation profile clustering

### 3. Trajectory Analysis
- **Serial NGS Data Processing**: Handle multiple timepoints per patient
- **VAF Slope Computation**: Calculate variant allele frequency trends over time
- **New Mutation Detection**: Identify mutations appearing after baseline
- **Enhanced Trajectory Features**: Compute comprehensive trajectory metrics including:
  - Mean VAF across timepoints
  - VAF variability (standard deviation)
  - VAF slope (trend over time)
  - VAF range (max - min)
  - New mutation slopes and counts

### 4. Advanced Analysis
- **Trajectory Clustering**: Cluster patients based on mutation trajectory patterns
- **Combined Analysis**: Integrate baseline and trajectory clustering results
- **Statistical Validation**: Silhouette analysis for both clustering approaches

### 5. Visualization and Output
- **PCA Plots**: Principal component analysis of mutation profiles
- **Heatmaps**: Comprehensive mutation matrix visualizations
- **Silhouette Plots**: Clustering validation visualizations
- **Per-cluster Analysis**: Individual cluster mutation patterns

## Input Files Required

1. **Clinical Excel File** (`CMML Project_2.xlsx`):
   - Sheet1: Clinical data with patient information
   - Sheet2: Serial NGS data (wide format)

2. **Mutation Data** (`CMML_Sheet2.csv`):
   - Wide format mutation data for baseline clustering

## Output Files Generated

### Data Files
- `CMML_clinical_baseline_clustered.csv` - Clinical data with baseline clusters
- `CMML_clinical_final_annotated.csv` - Final annotated clinical data with all clusters
- `trajectory_clustering_results.csv` - Trajectory clustering results and features
- `new_mutations_detailed_analysis.csv` - Detailed analysis of new mutations over time
- `CMML_Serial_long.csv` - Long format serial NGS data (if converted)

### Visualization Files
- `silhouette_analysis_baseline.png` - Baseline clustering validation
- `trajectory_silhouette_analysis.png` - Trajectory clustering validation
- `pca_plot_baseline_clusters.png` - PCA visualization of baseline mutation profiles
- `all_samples_baseline_mutation_heatmap.pdf` - Overall mutation heatmap with clusters
- `cluster_*_mutation_heatmap.pdf` - Individual cluster mutation heatmaps

## Key Enhancements

### Enhanced BMT Filtering
- Robust date parsing for Excel serial numbers and various date formats
- Automatic exclusion of post-BMT NGS data from trajectory analysis
- Detailed reporting of filtering effects

### New Mutation Analysis
- Detection of mutations appearing after baseline timepoint
- Computation of VAF slopes for new mutations
- Tracking of mutation evolution patterns

### Improved Clustering
- Enhanced trajectory feature computation
- Better handling of missing data and edge cases
- Comprehensive validation through silhouette analysis

### Advanced Visualizations
- Per-cluster mutation heatmaps with serial status annotation
- Improved color schemes using viridis palette
- Better annotation and legends

## Usage

```r
# Run the complete analysis
source("comprehensive_patient_trajectory_clustering.R")
```

### Configuration
Modify the following parameters at the top of the script as needed:

```r
clinical_excel <- "CMML Project_2.xlsx"   # Excel workbook path
clinical_sheet <- "Sheet1"                # Clinical data sheet
serial_sheet   <- "Sheet2"                # Serial NGS data sheet
mutation_file  <- "CMML_Sheet2.csv"       # Mutation data file
output_prefix  <- "baseline"              # Output file prefix
```

## Dependencies

The script automatically installs and loads required R packages:
- data.table
- tidyverse
- pheatmap
- RColorBrewer
- cluster
- NbClust
- readxl
- lubridate
- ggplot2
- reshape2
- survival
- survminer
- viridis

## Key Functions

### Utility Functions
- `extract_mrn_ngs()`: Extract MRN and NGS date columns
- `convert_wide_to_long()`: Convert wide format to long format
- `normalize_l2()`: L2 normalization for cosine similarity

### Analysis Functions
- `compute_new_mutation_slopes()`: Calculate VAF slopes for new mutations
- `get_new_mutation_details()`: Get detailed new mutation information

## Notes

- The script handles missing data gracefully
- BMT date filtering is optional and automatic based on data availability
- All visualizations are saved as high-resolution files
- The analysis pipeline is robust to various data format variations
- Comprehensive logging provides detailed progress information

## Error Handling

The script includes robust error handling for:
- Missing input files
- Data format inconsistencies
- Insufficient data for clustering
- Date parsing errors
- Missing annotations

For questions or issues, please refer to the inline documentation within the script.