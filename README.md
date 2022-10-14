[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/Anaconda)](https://cran.r-project.org/package=Anaconda) [![CRAN_Release_Date](https://www.r-pkg.org/badges/ago/Anaconda)](https://cran.r-project.org/package=Anaconda)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/Anaconda)](https://cran.r-project.org/package=Anaconda)

# Anaconda
ANACONDA : tArgeted differeNtial and globAl enriChment analysis of taxOnomic raNk by shareD Asvs
	
DOI: 10.13140/RG.2.2.11117.67048


## Install from CRAN:
`install.packages("Anaconda")`

## Install from Github: (much more up to date)

    install.packages("devtools")
    library(devtools)
    install_github("PLStenger/Anaconda")
    library("Anaconda")

    # If "Error in fetch(key) : lazy-load database "
    # Run:
    .rs.restartR()
    
 
## Why this name ?

ANACONDA: tArgeted differeNtial and globAl enriChment analysis of taxOnomic raNk by shareD Asvs

## What is it ?
 
Targeted differential and global enrichment analysis of taxonomic rank by shared ASVs (Amplicon Sequence Variant), for high-throughput eDNA sequencing of fungi, bacteria, and metazoan.
    Actually works in two steps: I) Targeted differential analysis from QIIME2 data and II) Global analysis by Taxon Mann-Whitney U test analysis from targeted analysis (I)
    (I) Estimate variance-mean dependence in count/abundance ASVs data from high-throughput sequencing assays and test for differential represented ASVs based on a model using the negative binomial distribution.
    (II) NCBITaxon_MWU uses continuous measure of significance (such as fold-change or -log(p-value)) to identify NCBITaxon that are significantly enriches with either up- or down-represented ASVs. If the measure is binary (0 or 1) the script will perform a typical 'NCBITaxon enrichment' analysis based Fisher's exact test: it will show NCBITaxon over-represented among the ASVs that have 1 as their measure. On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated ASVs. No colors are shown for binary measure analysis. The tree on the plot is hierarchical clustering of NCBITaxon based on shared ASVs. Categories with no branch length between them are subsets of each other. The fraction next to the category name indicates the fraction of 'good' ASVs in it; 'good' ASVs are the ones exceeding the arbitrary absValue cutoff (option in taxon_mwuPlot()). For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only. The original idea was for genes differential expression analysis from Wright et al (2015) <doi:10.1186/s12864-015-1540-2>; adapted here for taxonomic analysis.
    The 'Anaconda' package makes it possible to carry out these analyses by automatically creating several graphs and tables and storing them in specially created subfolders. You will need your QIIME2 pipeline output for each kingdom (eg; Fungi and/or Bacteria and/or Metazoan): i) taxonomy.tsv, ii) taxonomy_RepSeq.tsv, iii) ASV.tsv and iv) SampleSheet_comparison.txt (the latter being created by you).    
