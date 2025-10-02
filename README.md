# Anaconda

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/Anaconda)](https://cran.r-project.org/package=Anaconda) 
[![CRAN_Release_Date](https://www.r-pkg.org/badges/ago/Anaconda)](https://cran.r-project.org/package=Anaconda)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/Anaconda)](https://cran.r-project.org/package=Anaconda)

> **tArgeted differeNtial and globAl enriChment analysis of taxOnomic raNk by shareD Asvs**

<div align="center">
  <img src="/Stenger_2022b.png" width="250" alt="Anaconda package logo">
  <br><br>
  
  **DOI:** [10.13140/RG.2.2.11117.67048](https://doi.org/10.13140/RG.2.2.11117.67048)
</div>

---

## Official Publication

ðŸ“– **[Download paper from PLOS ONE](https://doi.org/10.1371/journal.pone.0311986)**

### Citation

**Primary citation:**
```
Stenger PL., LÃ©opold A., Dinh K., Mournet P., Robert N., Drouin J., Wamejonengo J., 
Russet S., Ibanez T., Maggia L., Carriconde F. (2025). Advancing biomonitoring of eDNA 
studies with the Anaconda R package: Integrating soil and One Health perspectives in the 
face of evolving traditional agriculture practices. PLOS ONE, 20(1), e0311986.
https://doi.org/10.1371/journal.pone.0311986
```

**Package citation:**
```
Stenger P-L. The R Package "Anaconda": Targeted Differential and Global Enrichment 
Analysis of Taxonomic Rank by Shared Asvs. CRAN; 2022. pp. 1â€“28. 
http://doi.org/10.13140/RG.2.2.11117.67048
```

---

## Installation

### From CRAN (Recommended)
```r
install.packages("Anaconda")
```

### From GitHub (Development version)
```r
install.packages("devtools")
library(devtools)
install_github("PLStenger/Anaconda")
library("Anaconda")

# If you encounter "Error in fetch(key) : lazy-load database" 
# Run:
.rs.restartR()
```

---

## What is Anaconda?

**Anaconda** is an R package designed for **targeted differential and global enrichment analysis** of taxonomic ranks using shared ASVs (Amplicon Sequence Variants). It supports high-throughput eDNA sequencing analysis of **fungi**, **bacteria**, and **metazoan**.

### Analysis Workflow

The package operates in **two main steps**:

#### **Step I: Targeted Differential Analysis** 
- Uses QIIME2 data with **DESeq2 algorithm**
- Estimates variance-mean dependence in count/abundance ASV data
- Tests for differential represented ASVs using **negative binomial distribution**

#### **Step II: Global Enrichment Analysis**
- **Taxon Mann-Whitney U test** analysis from targeted analysis results
- Uses continuous significance measures (fold-change or -log(p-value))
- Identifies taxonomic ranks significantly enriched with up/down-represented ASVs

---

## Key Features

- ðŸ§¬ **Multi-kingdom support**: Fungi, Bacteria, Metazoan
- ðŸ“ˆ **Statistical robustness**: Negative binomial modeling for count data
- ðŸŒ³ **Taxonomic hierarchy**: Analysis across all taxonomic ranks
- ðŸŽ¨ **Rich visualizations**: Heatmaps, PCA plots, MA plots, taxonomic trees
- ðŸ”— **Database integration**: FunGuild and Bactotraits databases
- ðŸ“ **Automated organization**: Creates structured output folders and files

---

## Required External Databases

> **Note**: These files are too large for the R package and must be downloaded separately.

### Download Links

| Kingdom | File | Download Link |
|---------|------|---------------|
| **Metazoan** | `taxonomy_all_metazoan_QIIME2_and_NCBI_format.txt` | [Download](https://www.dropbox.com/s/qpnkhmvskardlt7/taxonomy_all_metazoan_QIIME2_and_NCBI_format.txt?dl=0) |
| **Bacteria** | `taxonomy_all_bacteria_QIIME2_and_NCBI_format.txt` | [Download](https://www.dropbox.com/s/hzu31gm5ivvh3cu/taxonomy_all_bacteria_QIIME2_and_NCBI_format.txt?dl=0) |
| **Fungi** | `taxonomy_all_fungi_QIIME2_and_NCBI_format.txt` | [Download](https://www.dropbox.com/s/yulal8i56ewcv45/taxonomy_all_fungi_QIIME2_and_NCBI_format.txt?dl=0) |

### Optional Database Files

| File | Purpose | Download Link |
|------|---------|---------------|
| `ncbitaxon_ontology.obo` | Custom database creation | [Download](https://www.dropbox.com/s/m74cxd46wk31w53/ncbitaxon_ontology.obo?dl=0) |

---

## Input Files Required

Before running Anaconda, prepare these **4 essential files** from your QIIME2 pipeline:

| File | Description |
|------|-------------|
| `ASV.tsv` | ASV abundance table for all samples |
| `taxonomy.tsv` | Taxonomy-ASV mapping (rarefied dataset) |
| `taxonomy_RepSeq.tsv` | Taxonomy-ASV mapping (representative sequences) |
| `SampleSheet_comparison.txt` | **User-created** sample metadata file |

### SampleSheet_comparison.txt Format

```txt
Sample_Name    Input_File         Condition
F1             input_F1.txt       F
F2             input_F2.txt       F
F3             input_F3.txt       F
LF1            input_LF1.txt      LF
LF2            input_LF2.txt      LF
SF1            input_SF1.txt      SF
```

---

## Quick Start Guide

### Basic Setup
```r
library(Anaconda)

# Save original directory
original_dir <- getwd()

# Place your 4 input files in the working directory
```

### Choose Your Kingdom
```r
# Run analysis for specific kingdom
Fungi()          # For fungal analysis
# OR
Bacteria()       # For bacterial analysis
```

### Organize Files
```r
# Move into kingdom directory
setwd("Fungi")  # or "Bacteria"
kingdom <- getwd()

# Move input files to kingdom folder
move_files()
```

### Targeted Analysis (Step I)
```r
# Prepare taxonomy and create input files
taxo <- get_input_files()
setwd("01_Targeted_analysis")
targeted_analysis_dir <- getwd()

# Import sample information
target_file <- target_file()
samplesInfo <- samplesInfo()

# Set minimum ASV threshold
threshold <- 1

# Create differential ASV abundance object
dasva <- get_dasva(fitType="parametric")

# Quality control plots
plotDispASVs(dasva)      # Dispersion plot
plotSparsityASV(dasva)   # Sparsity plot

# PCA analysis
data <- PCA_data_dasva()
# Use 'data' object with ggplot2 for custom PCA plots

# Heatmap of top 75 most abundant ASVs
log2.norm.counts <- heatmap_data_dasva()
colnames(log2.norm.counts) <- NULL
heatmap_condition_df <- heatmap_condition()
pheatmap(log2.norm.counts, annotation_col = heatmap_condition_df)
```

### Differential Analysis
```r
# Compare conditions (example: Forest vs Long Fallow)
res_forest_vs_long_fallow <- results(
  dasva, 
  contrast = c("condition", "F", "LF"), 
  pAdjustMethod = "BH", 
  alpha = 0.01
)

# Export results with taxonomy
write.table(
  merge(data.frame(res_forest_vs_long_fallow), taxo, 
        by.x="row.names", by.y="Feature.ID"),
  "Table_res_forest_vs_long_fallow_taxonomy.txt", 
  sep="\t"
)

# Add functional annotations (Fungi only)
res_forest_vs_long_fallow_guilds <- funguild_input_targeted(res_forest_vs_long_fallow)
get_funguilds_targeted(res_forest_vs_long_fallow_guilds)

# Add Bactotraits (Bacteria only)
get_bactotraits_targeted(res_forest_vs_long_fallow)

# Create MA plot
plotMA.dasva(res_forest_vs_long_fallow, alpha=0.01)
```

### Global Analysis (Step II)
```r
# Return to kingdom directory
setwd(kingdom)

# Create annotation database
database_fungi_creation()     # For fungi
# OR
database_bacteria_creation()  # For bacteria

# Set database path
taxon_Annotations <- "database_fungi_package_all.tab"  # Adjust for bacteria

# Move to global analysis directory
setwd("02_Global_analysis")

# Prepare input for MWU test
input_res_forest_vs_long_fallow <- input_global_analysis(res_forest_vs_long_fallow)
write.table(
  input_res_forest_vs_long_fallow, 
  "input_res_forest_vs_long_fallow.txt", 
  sep=",", quote = FALSE, row.names=FALSE
)

# Set paths
input <- "input_res_forest_vs_long_fallow.txt"
taxon_Database <- file.path(original_dir, "Working_scripts/ncbitaxon_ontology.obo")

# Run Mann-Whitney U test
taxon_mwuStats_res <- taxon_mwuStats(
  input, taxon_Database, taxon_Annotations, 
  TR, perlPath="perl", 
  largest=0.1, smallest=1, clusterCutHeight=0.1
)

# Generate taxonomic tree plot
taxon_mwuPlot(
  input, taxon_Annotations, taxon_Division, 
  absValue= -log(0.05,10), 
  level1=0.1, level2=0.05, level3=0.01, 
  txtsize=1.2, treeHeight=0.5
)
```

### Enhanced Fungi Analysis (Optional)
```r
# Add FunGuild information to taxonomic analysis
taxon_list <- taxon_mwu_list(input, taxon_Annotations, taxon_Division)
taxon_list_drawer <- get_taxon_list_drawer(taxon_list)
funguilds <- get_funguilds(taxon_list_drawer)
link_guilds <- get_link_guilds(taxon_list, funguilds)

# Plot with guild information
taxon_mwuPlot_guilds(
  input, taxon_Annotations, taxon_Division, 
  absValue= -log(0.05,10), 
  level1=0.1, level2=0.05, level3=0.01, 
  txtsize=1.2, treeHeight=0.5
)
```

---

## System Requirements

- **R version**: â‰¥ 3.5.0
- **Python**: â‰¥ 2.7 (for FunGuild functionality)
- **Perl**: Required for global analysis (recommend [Strawberry Perl](https://strawberryperl.com))

---

## Getting Started

1. **Download** the main script `Anaconda.R`
2. **Follow** the step-by-step workflow in the script
3. **Download** required taxonomy databases before analysis
4. **Ensure** all system requirements are installed

---

## Scientific Background

The Anaconda package adapts differential expression analysis methods (originally from Wright et al., 2015) for taxonomic analysis of high-throughput sequencing data. Instead of analyzing gene expression differences, Anaconda focuses on:

- **ASV abundance differences** between experimental conditions
- **Taxonomic enrichment** across hierarchical levels
- **Functional annotations** through established databases

This approach enables researchers to identify not just individual species differences, but broader taxonomic patterns that may indicate ecological shifts or functional changes in microbial communities.

---

## Contributing

Found a bug or have suggestions? Please visit our [GitHub repository](https://github.com/PLStenger/Anaconda) to:
- Report issues
- Suggest improvements  
- Contribute code
- Access development protocols

---

## License

This package is distributed under standard CRAN licensing terms. See the DESCRIPTION file for details.

---

<div align="center">
  <strong>Happy analyzing! ðŸ§¬ðŸ”¬</strong>
</div>
