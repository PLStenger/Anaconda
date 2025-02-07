# ANACONDA
# tArgeted differeNtial and globAl enriChment analysis of taxOnomic raNk by shareD Asvs

rm(list=ls())

# Adapt this line to your own directory file
setwd("~/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/17_Meta-analysis_soil_ultramafism/14_Anaconda")

# Integrate your working directory into "original_dir" object
original_dir= getwd()

library(Anaconda)
library(ggrepel)
library(pheatmap)
library(rafalib)

# Add python
# Add perl (strawberry)

# Add taxonomy_all_bacteria_QIIME2_and_NCBI_format.txt
# Add taxonomy_all_fungi_QIIME2_and_NCBI_format.txt
# Add ncbitaxon_ontology.obo

# In Working_scripts folder

#source(file.path(original_dir, paste("Working_scripts/Anaconda.functions.R", sep="")))

#library(DESeq2)
#library(ggplot2)
#library(ggrepel)
#library(pheatmap)
#library(rafalib)
#library(RColorBrewer)
#library(gplots)
#library(VennDiagram)
#library(ape)
#library(lookup)
#library(plyr)
#library(dplyr)
#library(data.table)

###################################################################################################
#
#     I - DESeq2 analysis from QIIME2 data
#
###################################################################################################

####################################################################################################
# Set parameters for targeted analysis
####################################################################################################

# #Files needed in your root:

## Install Python + Perl (https://strawberryperl.com)
##
## In the package

## IN 'ROOTS' :
## Anaconda.R

## IN 'Working_scripts' FOLDER :
## Anaconda.functions.R
## BactoTraits_database
## Guilds_v1.1.py
## mwu_a_NCBITaxon.pl
## mwu_b_NCBITaxon.pl
## ncbitaxon_ontology.obo
## search_taxonomic_drawer.sh
## taxonomy_all_bacteria_QIIME2_and_NCBI_format_duplicated.txt
## taxonomy_all_bacteria_QIIME2_and_NCBI_format.txt
## taxonomy_all_fungi_QIIME2_and_NCBI_format.txt

## FILES you need to add in the ROOTS (not in 'Input_Fungi' or 'Input_Bacteria' folder) :
## taxonomy.tsv
## taxonomy_RepSeq.tsv
## SampleSheet_comparison.txt
## ASV.tsv

# - taxonomy.tsv is in Diversity_in_Mare_yam_crop/05_QIIME2/Paired_end/V4/export/taxonomy/taxonomy_reads-per-batch_RarRepSeq
# - taxonomy_RepSeq.tsv (named change, normally it's also named "taxonomy.tsv") is in Diversity_in_Mare_yam_crop/05_QIIME2/Paired_end/V4/export/taxonomy/taxonomy_reads-per-batch_RepSeq
# - ASV.tsv is in Diversity_in_Mare_yam_crop/05_QIIME2/Paired_end/V4/export/subtables/RarTable-all/
# - SampleSheet_comparison.txt is manually made

# Is your data are fungi or bacteria?
# run this if it's Fungi :
Fungi() # this function create a new folder named Fungi and set your working directory into this folder
setwd("Fungi")
# Or run this if it's bacteria :
Bacteria() # this function create a new folder named Bacteria and set your working directory into this folder
setwd("Bacteria")

# For remember the path
kingdom <- getwd()

# Move the file in the good folders
move_files()

# Created sub directory "Targeted_analysis" if not already exist
# Then, create one file by condition into it, and then upload the taxonomy file
taxo <- get_input_files()
setwd("01_Targeted_analysis")

# Integrate path to current directory containing these count files into "targeted_analysis_dir" object
targeted_analysis_dir= getwd()

####################################################################################################
# Imports conditions informations
####################################################################################################

target_file <- target_file()
samplesInfo <- samplesInfo()

threshold=1 # minimum of ASVs value across samples

####################################################################################################
# Creates dASVa object
####################################################################################################

# dasva = dASVa = Differential ASV abundance
# Fit a Gamma-Poisson Generalized Linear Model
# https://bioconductor.org/packages/release/bioc/html/glmGamPoi.html
# dispersion estimates for Negative Binomial distributed data.
dasva <- get_dasva(fitType="parametric")
#dasva <- get_dasva(fitType="local")
#dasva <- get_dasva(fitType="mean")

####################################################################################################
# Dispersion plot
####################################################################################################

plotDispASVs(dasva)

##################################################
# Save it into PDF file
pdf("plot_01_dispersion_plot.pdf")
plotDispASVs(dasva)
dev.off()
##################################################

####################################################################################################
# Sparsity plot
####################################################################################################

plotSparsityASV(dasva)

##################################################
# Save it into PDF file
pdf("plot_02_sparsity_plot.pdf")
plotSparsityASV(dasva)
dev.off()
##################################################

####################################################################################################
# PCA
####################################################################################################

data <- PCA_data_dasva()

ggplot(data, aes(PC1, PC2, color=condition)) +
  geom_point(size=2.5) +
  xlab(paste0("PC1: ",round(100*attr(data, "percentVar"))[1],"% variance")) +
  ylab(paste0("PC2: ",round(100*attr(data, "percentVar"))[2],"% variance")) +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0)) +
  theme_bw()+
  #scale_color_manual(values=c("dodgerblue4", "firebrick4", "forestgreen")) +
  geom_text_repel(aes(label=target_file$label))

##################################################
# Save it into PDF file
pdf("plot_03_PCA_plot.pdf")
ggplot(data, aes(PC1, PC2, color=condition)) +
  geom_point(size=2.5) +
  xlab(paste0("PC1: ",round(100*attr(data, "percentVar"))[1],"% variance")) +
  ylab(paste0("PC2: ",round(100*attr(data, "percentVar"))[2],"% variance")) +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0)) +
  theme_bw()+
  #scale_color_manual(values=c("dodgerblue4", "firebrick4", "forestgreen")) +
  geom_text_repel(aes(label=target_file$label))
dev.off()
##################################################

####################################################################################################
# Heatmap
####################################################################################################

log2.norm.counts <- heatmap_data_dasva()
colnames(log2.norm.counts) <- NULL
heatmap_condition_df <- heatmap_condition()

#Draw heatmap
pheatmap(log2.norm.counts,
         clustering_distance_cols = "euclidean",
         clustering_method = "average",
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         annotation_col=heatmap_condition_df,
         show_rownames=TRUE,
         show_colnames=TRUE,
         fontsize_row=5,
         fontsize_col=10,
         fontsize=8,
         labels_col=rownames(heatmap_condition_df),
         main=paste("Clustered heatmap of 75 most abundant ASVs\n","euclidean"," distance with ","average", " clustering method",sep=""))

##################################################
# Save it into PDF file
pheatmap(log2.norm.counts,
         filename = "plot_04_pheatmap_log2_norm_counts.pdf",
        # width = 15,
         #height = 15,
         clustering_distance_cols = "euclidean",
         clustering_method = "average",
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         annotation_col=heatmap_condition_df,
         show_rownames=TRUE,
         show_colnames=TRUE,
         fontsize_row=5,
         fontsize_col=10,
         fontsize=8,
         labels_col=rownames(heatmap_condition_df),
         main=paste("Clustered heatmap of 75 most abundant ASVs\n","euclidean"," distance with ","average", " clustering method",sep=""))
dev.off()
##################################################

log2.norm.counts_taxo <- heatmap_taxo()

pheatmap(log2.norm.counts_taxo,
         clustering_distance_cols = "euclidean",
         clustering_method = "average",
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         #annotation_col=dat3,
         show_rownames=TRUE,
         show_colnames=TRUE,
         fontsize_row=8,
         fontsize_col=10,
         fontsize=8,
        # annotation_colors = ann_colors,
         labels_col=colnames(log2.norm.counts_taxo),
         main=paste("Clustered heatmap of 75 most abundant ASVs\n","euclidean"," distance with ","average", " clustering method",sep=""))

##################################################
# Save it into PDF file
pheatmap(log2.norm.counts_taxo,
         filename = "plot_05_pheatmap_log2_norm_counts_with_taxonomy.pdf",
         width = 15,
         height = 15,
         clustering_distance_cols = "euclidean",
         clustering_method = "average",
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         #annotation_col=dat3,
         show_rownames=TRUE,
         show_colnames=TRUE,
         fontsize_row=8,
         fontsize_col=10,
         fontsize=8,
         # annotation_colors = ann_colors,
         labels_col=colnames(log2.norm.counts_taxo),
         main=paste("Clustered heatmap of 75 most abundant ASVs\n","euclidean"," distance with ","average", " clustering method",sep=""))
dev.off()
##################################################

write.table(log2.norm.counts, "Table_01_log2.norm.counts.txt")


####################################################################################################
# Clustering
####################################################################################################

hh <- hclust(dist(t(log2.norm.counts), method = "euclidean"), method = "average")

# For fungi
myplclust(hh, labels=rownames(heatmap_condition_df), lab.col=c(rep("forestgreen",4),rep("olivedrab3",3),rep("darkgoldenrod4",4)))

# For bacteria
myplclust(hh, labels=rownames(heatmap_condition_df), lab.col=c("dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "forestgreen", "forestgreen", "forestgreen", "forestgreen", "forestgreen"))

##################################################
# Save it into PDF file
pdf("plot_06_clustering.pdf")
myplclust(hh, labels=rownames(heatmap_condition_df), lab.col=c("dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "forestgreen", "forestgreen", "forestgreen", "forestgreen", "forestgreen"))
dev.off()
##################################################

####################################################################################################
# Heatmap sample to sample
####################################################################################################


heatmap_samples_matrix <- function (nothing)
{
  
  sampleFiles <- grep("input_",list.files(targeted_analysis_dir),value=TRUE)
  sampleCondition <- samplesInfo$V3
  sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles)
  sampleTable <- data.frame(sampleName = sampleFiles,
                            fileName = sampleFiles,
                            condition = sampleCondition)
  
  dasva_raw <- dasva_raw_input(sampleTable = sampleTable,
                               directory = targeted_analysis_dir,
                               design= ~ condition)
  
  #Trims too low represented ASVs
  dasva_raw <- dasva_raw[ rowSums(counts(dasva_raw)) > threshold, ] ; dasva_raw
  
  rld <- rlog(dasva_raw, blind=FALSE)
  distsRL <- dist(t(assay(rld)))
  mat <- as.matrix(distsRL)
  # rownames(mat) <- colnames(mat) <- with(colData(dasva), paste(Id, condition , sep=' : '))
  rownames(mat)  <- gsub(".txt", "", rownames(mat))
  rownames(mat)  <- gsub("input_", "", rownames(mat))
  colnames(mat)  <- gsub(".txt", "", colnames(mat))
  colnames(mat)  <- gsub("input_", "", colnames(mat))
  return(mat)
}

heatmap_samples_hclust <- function (nothing)
{
  
  sampleFiles <- grep("input_",list.files(targeted_analysis_dir),value=TRUE)
  sampleCondition <- samplesInfo$V3
  sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles)
  sampleTable <- data.frame(sampleName = sampleFiles,
                            fileName = sampleFiles,
                            condition = sampleCondition)
  
  dasva_raw <- dasva_raw_input(sampleTable = sampleTable,
                               directory = targeted_analysis_dir,
                               design= ~ condition)
  
  #Trims too low represented ASVs
  dasva_raw <- dasva_raw[ rowSums(counts(dasva_raw)) > threshold, ] ; dasva_raw
  
  rld <- rlog(dasva_raw, blind=FALSE)
  distsRL <- dist(t(assay(rld)))
  hc <- hclust(distsRL)
  
  return(hc)
}

mat <- heatmap_samples_matrix()
hc <- heatmap_samples_hclust()

#  You need to adapt the number of color like above

#install.packages("gplots")
library(gplots)
library(RColorBrewer)

heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace='none',
          col=rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)),
          margin=c(10, 10),
          colRow=c("dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "forestgreen", "forestgreen", "forestgreen", "forestgreen", "forestgreen"),
          colCol=c("dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "forestgreen", "forestgreen", "forestgreen", "forestgreen", "forestgreen"))

dev.off()


heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace='none',
          col=rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)),
          margin=c(10, 10),
          colRow=c("dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "forestgreen", "forestgreen", "forestgreen", "forestgreen", "forestgreen"),
          colCol=c("dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "forestgreen", "forestgreen", "forestgreen", "forestgreen", "forestgreen"))

dev.off()



##################################################
# Save it into PDF file
pdf("plot_07_heatmap_samplebysample.pdf")
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace='none',
          col=rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)),
          margin=c(10, 10),
          colRow=c("dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "forestgreen", "forestgreen", "forestgreen", "forestgreen", "forestgreen"),
          colCol=c("dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "firebrick4", "forestgreen", "forestgreen", "forestgreen", "forestgreen", "forestgreen"))
dev.off()
##################################################

####################################################################################################
# Differential ASV abundance (DASVA) analysis
####################################################################################################


res <- results(dasva, pAdjustMethod = "BH", alpha=0.01)
res

# Comparisons
res_noumea_vs_voh <- results(dasva,contrast=c("condition","Noumea","Voh"), pAdjustMethod = "BH", alpha=0.01)
res_poya_vs_voh <- results(dasva,contrast=c("condition","Poya","Voh"), pAdjustMethod = "BH", alpha=0.01)
res_noumea_vs_poya <- results(dasva,contrast=c("condition","Noumea","Poya"), pAdjustMethod = "BH", alpha=0.01)

# Export differential ASVs lists
write.table(res_noumea_vs_voh,"Table_02_01_res_noumea_vs_voh.txt", sep="\t")
write.table(res_poya_vs_voh,"Table_02_01_res_poya_vs_voh.txt", sep="\t")
write.table(res_noumea_vs_poya,"Table_02_01_res_noumea_vs_poya.txt", sep="\t")

# Adding the taxonomy
write.table(merge(data.frame(res_noumea_vs_voh), taxo, by.x="row.names", by.y="Feature.ID"),"Table_02_02_res_noumea_vs_voh_taxonomy.txt", sep="\t")
write.table(merge(data.frame(res_poya_vs_voh), taxo, by.x="row.names", by.y="Feature.ID"),"Table_02_02_res_poya_vs_voh_taxonomy.txt", sep="\t")
write.table(merge(data.frame(res_noumea_vs_poya), taxo, by.x="row.names", by.y="Feature.ID"),"Table_02_02_res_noumea_vs_poya_taxonomy.txt", sep="\t")

# Adding FunGuilds for FUNGI
res_forest_vs_long_fallow_guilds <- funguild_input_targeted(res_forest_vs_long_fallow)
get_funguilds_targeted(res_forest_vs_long_fallow_guilds)

res_forest_vs_short_fallow_guilds <- funguild_input_targeted(res_forest_vs_short_fallow)
get_funguilds_targeted(res_forest_vs_short_fallow)

res_short_fallow_vs_long_fallow_guilds <- funguild_input_targeted(res_short_fallow_vs_long_fallow)
get_funguilds_targeted(res_short_fallow_vs_long_fallow_guilds)


# Adding BactoTraits for BACTERIA
get_bactotraits_targeted(res_forest_vs_long_fallow)
write.table(get_bactotraits_targeted(res_forest_vs_long_fallow),"Table_02_03_res_forest_vs_long_fallow_taxonomy_BactoTraits.txt", sep="\t")

get_bactotraits_targeted(res_forest_vs_short_fallow)
write.table(get_bactotraits_targeted(res_forest_vs_short_fallow),"Table_02_03_res_forest_vs_short_fallow_fallow_taxonomy_BactoTraits.txt", sep="\t")

get_bactotraits_targeted(res_short_fallow_vs_long_fallow)
write.table(get_bactotraits_targeted(res_short_fallow_vs_long_fallow),"Table_02_03_res_short_fallow_vs_long_fallow_taxonomy_BactoTraits.txt", sep="\t")


####################################################################################################
# MA plots
####################################################################################################


plotMA.dasva(res_noumea_vs_poya, main=expression(bold(paste("Noumea ", italic("vs. "), "Poya"))), ylim=c(-20,20),alpha=0.01)

##################################################
# Save them into PDF files
pdf("plot_08_plotMA.dasva_Forest_vs_Long_fallow.pdf")
plotMA.dasva(res_forest_vs_long_fallow, main="Forest vs Long fallow", ylim=c(-20,20),alpha=0.01)
dev.off()
pdf("plot_08_plotMA.dasva_Forest_vs_Short_fallow.pdf")
plotMA.dasva(res_forest_vs_short_fallow, main="Forest vs Short fallow", ylim=c(-20,20),alpha=0.01)
dev.off()
pdf("plot_08_plotMA.dasva_Short_fallow_vs_Long_fallow.pdf")
plotMA.dasva(res_short_fallow_vs_long_fallow, main="Short fallow vs Long fallow", ylim=c(-20,20),alpha=0.01)
dev.off()
##################################################

#Selects results with p-value <= 0.05 et |log2FC| >= 2
res_forest_vs_long_fallow_subet <- res_forest_vs_long_fallow[ which(res_forest_vs_long_fallow$padj<=0.05 & abs(res_forest_vs_long_fallow$log2FoldChange) >= 2),]
res_forest_vs_short_fallow_subet <- res_forest_vs_short_fallow[ which(res_forest_vs_short_fallow$padj<=0.05 & abs(res_forest_vs_short_fallow$log2FoldChange) >= 2),]
res_short_fallow_vs_long_fallow_subet <- res_short_fallow_vs_long_fallow[ which(res_short_fallow_vs_long_fallow$padj<=0.05 & abs(res_short_fallow_vs_long_fallow$log2FoldChange) >= 2),]

write.table(res_forest_vs_long_fallow_subet,"Table_03_res_forest_vs_long_fallow_subset.txt", sep="\t")
write.table(res_forest_vs_short_fallow_subet,"Table_03_res_forest_vs_short_fallow_subset.txt", sep="\t")
write.table(res_short_fallow_vs_long_fallow_subet,"Table_03_res_short_fallow_vs_long_fallow_subset.txt", sep="\t")

####################################################################################################
# Venn diagrams FUNGI
####################################################################################################

length(rownames(res_forest_vs_long_fallow_subet))
length(rownames(res_forest_vs_short_fallow_subet))
length(rownames(res_short_fallow_vs_long_fallow_subet))

temp1 <- intersect(rownames(res_forest_vs_long_fallow_subet), rownames(res_forest_vs_short_fallow_subet)) ; length(temp1)
temp2 <- intersect(rownames(res_forest_vs_long_fallow_subet), rownames(res_short_fallow_vs_long_fallow_subet)) ; length(temp2)
temp3 <- intersect(rownames(res_forest_vs_short_fallow_subet), rownames(res_short_fallow_vs_long_fallow_subet)) ; length(temp3)
temp123 <- intersect(temp1,temp2); length(temp123)

library(VennDiagram)

dev.off()
draw.triple.venn(length(rownames(res_forest_vs_long_fallow_subet)),length(rownames(res_forest_vs_short_fallow_subet)),length(rownames(res_short_fallow_vs_long_fallow_subet)),length(temp1),length(temp3),length(temp2),length(temp123),
                 category=c(paste("Forest\n vs\n Long fallow\nn=",length(rownames(res_forest_vs_long_fallow_subet)),sep=""),paste("Forest\n vs\n Short fallow\nn=",length(rownames(res_forest_vs_short_fallow_subet)),sep=""),paste("Short fallow\n vs\n long fallow\nn=",length(rownames(res_short_fallow_vs_long_fallow_subet)),sep="")),
                 fill=c("forestgreen","chartreuse4","darkgreen"),
                 cat.col=c("forestgreen","chartreuse4","darkgreen"), cat.fontfamily=c("sans","sans","sans"),margin=0.05,cat.dist=c(0.07,0.07,0.05))


# Save it into PDF file
pdf("plot_09_Venn_diagram.pdf")
draw.triple.venn(length(rownames(res_forest_vs_long_fallow_subet)),length(rownames(res_forest_vs_short_fallow_subet)),length(rownames(res_short_fallow_vs_long_fallow_subet)),length(temp1),length(temp3),length(temp2),length(temp123),
                 category=c(paste("Forest\n vs\n Long fallow\nn=",length(rownames(res_forest_vs_long_fallow_subet)),sep=""),paste("Forest\n vs\n Short fallow\nn=",length(rownames(res_forest_vs_short_fallow_subet)),sep=""),paste("Short fallow\n vs\n long fallow\nn=",length(rownames(res_short_fallow_vs_long_fallow_subet)),sep="")),
                 fill=c("forestgreen","chartreuse4","darkgreen"),
                 cat.col=c("forestgreen","chartreuse4","darkgreen"), cat.fontfamily=c("sans","sans","sans"),margin=0.05,cat.dist=c(0.07,0.07,0.05))
dev.off()

####################################################################################################
# Example of how extract some specific ASVs
####################################################################################################

################################################################
# Extract the 13 Specific ASVs from forest vs short fallow :
list <- rownames(res_forest_vs_short_fallow_subet)
tp1 <- list[-pmatch(temp1,list)]
tp2 <- tp1[-pmatch(temp3,tp1)]
write.table(tp2,"Table_04_forest_vs_short_fallow_subet_13_ASVs_specific.txt", sep="\t")

################################################################

df <- data.frame(
  x = c(1,2),
  y = c(length(which(res_forest_vs_long_fallow_subet$log2FoldChange >=0)), length(which(res_forest_vs_long_fallow_subet$log2FoldChange <0))),
  grp = c("More quantities", "Less quantities")
)

ggplot(data = df, aes(x, y, group = grp)) +
  geom_col(aes(fill = grp)) +
  ylim(0, 20) +
  geom_text(aes(label = y), position = position_stack(vjust = 0.5), size = 10) +
  scale_fill_manual(values=c("firebrick3", "dodgerblue3")) +
  theme(legend.position = "none") +
  ggtitle("Forest vs long fallow") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))



df <- data.frame(
  x = c(1,2),
  y = c(length(which(res_forest_vs_short_fallow_subet$log2FoldChange >=0)), length(which(res_forest_vs_short_fallow_subet$log2FoldChange <0))),
  grp = c("More quantities", "Less quantities")
)

ggplot(data = df, aes(x, y, group = grp)) +
  geom_col(aes(fill = grp)) +
  ylim(0, 20) +
  geom_text(aes(label = y), position = position_stack(vjust = 0.5), size = 10) +
  ggtitle("Forest vs short fallow") +
  scale_fill_manual(values=c("firebrick3", "dodgerblue3")) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))



df <- data.frame(
  x = c(1,2),
  y = c(length(which(res_short_fallow_vs_long_fallow_subet$log2FoldChange >=0)), length(which(res_short_fallow_vs_long_fallow_subet$log2FoldChange <0))),
  grp = c("More quantities", "Less quantities")
)

ggplot(data = df, aes(x, y, group = grp)) +
  geom_col(aes(fill = grp)) +
  ylim(0, 20) +
  geom_text(aes(label = y), position = position_stack(vjust = 0.5), size = 10) +
  scale_fill_manual(values=c("firebrick3", "dodgerblue3")) +
  theme(legend.position = "none") +
  ggtitle("Short fallow vs long fallow") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))






####################################################################################################
# Venn diagrams BACTERIA
####################################################################################################

length(rownames(res_Pa_5_vs_Pa_10_subset)) # 11
length(rownames(res_Pf_5_vs_Pf_10_subset)) # 7
length(rownames(res_Pf_5_vs_Pa_5_subset)) # 8
length(rownames(res_Pf_10_vs_Pa_10_subset)) # 18

temp12 <- intersect(rownames(res_Pa_5_vs_Pa_10_subset), rownames(res_Pf_5_vs_Pf_10_subset)) ; length(temp12)
temp13 <- intersect(rownames(res_Pa_5_vs_Pa_10_subset), rownames(res_Pf_5_vs_Pa_5_subset)) ; length(temp13)
temp14 <- intersect(rownames(res_Pa_5_vs_Pa_10_subset), rownames(res_Pf_10_vs_Pa_10_subset)) ; length(temp14)
temp23 <- intersect(rownames(res_Pf_5_vs_Pf_10_subset), rownames(res_Pf_5_vs_Pa_5_subset)) ; length(temp23)
temp24 <- intersect(rownames(res_Pf_5_vs_Pf_10_subset), rownames(res_Pf_10_vs_Pa_10_subset)) ; length(temp24)
temp34 <- intersect(rownames(res_Pf_5_vs_Pa_5_subset), rownames(res_Pf_10_vs_Pa_10_subset)) ; length(temp34)

temp123 <- intersect(temp12,temp13); length(temp123)
temp124 <- intersect(temp12,temp24); length(temp124)
temp134 <- intersect(temp13,temp34); length(temp134)
temp234 <- intersect(temp23,temp34); length(temp234)

temp1234 <- intersect(temp123,temp234); length(temp1234)


dev.off()


venn.plot <- draw.quad.venn(
  area1 = length(rownames(res_Pa_5_vs_Pa_10_subset)),
  area2 = length(rownames(res_Pf_5_vs_Pf_10_subset)),
  area3 = length(rownames(res_Pf_5_vs_Pa_5_subset)),
  area4 = length(rownames(res_Pf_10_vs_Pa_10_subset)) ,
  n12 = length(temp12),
  n13 = length(temp13),
  n14 = length(temp14),
  n23 = length(temp23),
  n24 = length(temp24),
  n34 = length(temp34),
  n123 = length(temp123),
  n124 = length(temp124),
  n134 = length(temp134),
  n234 = length(temp234),
  n1234 = length(temp1234),
  category = c("Pa 5 vs Pa 10", "Pf 5 vs Pf 10", "Pf 5 vs Pa 5", "Pf 10 vs Pa 10"),
  fill = c("orange", "red", "green", "blue"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("orange", "red", "green", "blue")
)



##################################################
# Save it into PDF file
pdf("plot_09_Venn_diagram.pdf")

venn.plot <- draw.quad.venn(
  area1 = length(rownames(res_Pa_5_vs_Pa_10_subset)),
  area2 = length(rownames(res_Pf_5_vs_Pf_10_subset)),
  area3 = length(rownames(res_Pf_5_vs_Pa_5_subset)),
  area4 = length(rownames(res_Pf_10_vs_Pa_10_subset)) ,
  n12 = length(temp12),
  n13 = length(temp13),
  n14 = length(temp14),
  n23 = length(temp23),
  n24 = length(temp24),
  n34 = length(temp34),
  n123 = length(temp123),
  n124 = length(temp124),
  n134 = length(temp134),
  n234 = length(temp234),
  n1234 = length(temp1234),
  category = c("Pa 5 vs Pa 10", "Pf 5 vs Pf 10", "Pf 5 vs Pa 5", "Pf 10 vs Pa 10"),
  fill = c("orange", "red", "green", "blue"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("orange", "red", "green", "blue")
)

dev.off()
##################################################


###################################################################################################
#
#     II - Global analysis by Taxon_MWU analysis from targeted analysis (I)
#
###################################################################################################

####################################################################################################
# Set parameters for global analysis
####################################################################################################

# NCBITaxon_MWU uses continuous measure of significance (such as fold-change or -log(p-value)) to identify NCBITaxon that are significantly enriches with either up- or down-represented ASVs. The advantage - no need to impose arbitrary significance cutoff.
# If the measure is binary (0 or 1) the script will perform a typical "NCBITaxon enrichment" analysis based Fisher's exact test: it will show NCBITaxon over-represented among the ASVs that have 1 as their measure.
# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated ASVs. No colors are shown for binary measure analysis.
# The tree on the plot is hierarchical clustering of NCBITaxon based on shared ASVs. Categories with no branch length between them are subsets of each other.
# The fraction next to GO category name indicates the fraction of "good" ASVs in it; "good" ASVs being the ones exceeding the arbitrary absValue cutoff (option in taxon_mwuPlot). For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only.
# Stretch the plot manually to match tree to text

# original idea, scripts, and text for genes from Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu
# adapted for ASVs by Pierre-Louis Stenger, IAC, Noumea, New-Caledonia; august 2021; Pierrelouis.stenger@gmail.com

####################################################################################################
# Edit these to match your data file names:

# If you re-start analysis, you will be need to run this :
#res_forest_vs_long_fallow <- read.table("Table_02_res_forest_vs_long_fallow.txt", sep="\t")
#res_forest_vs_short_fallow <- read.table("Table_02_res_forest_vs_short_fallow.txt", sep="\t")
#res_short_fallow_vs_long_fallow <- read.table("Table_02_res_short_fallow_vs_long_fallow.txt", sep="\t")

setwd(kingdom)

#source(file.path(original_dir, paste("Working_scripts/Anaconda.functions.R", sep="")))

#########################################################################################################
# Database creation
#########################################################################################################

database_fungi_creation_RepSeq()
database_bacteria_creation()
setwd("02_Global_analysis")

#########################################################################################################
# Input files creation for each condition
#########################################################################################################

input_res_noumea_vs_poya <- input_global_analysis(res_noumea_vs_poya)
write.table(input_res_noumea_vs_poya, "input_res_noumea_vs_poya.txt", sep=",", quote = FALSE, row.names=FALSE)

#res_noumea_vs_voh <- read.table("Table_02_01_res_noumea_vs_voh.txt", header=T)
input_res_noumea_vs_voh <- input_global_analysis(res_noumea_vs_voh)
write.table(input_res_noumea_vs_voh, "input_res_noumea_vs_voh.txt", sep=",", quote = FALSE, row.names=FALSE)

#res_poya_vs_voh <- read.table("Table_02_01_res_poya_vs_voh.txt", header=T)
input_res_poya_vs_voh <- input_global_analysis(res_poya_vs_voh)
write.table(input_res_poya_vs_voh, "input_res_poya_vs_voh.txt", sep=",", quote = FALSE, row.names=FALSE)



####################################################################################################
# Launch analysis
####################################################################################################

input="input_res_noumea_vs_poya.txt" 
input="input_res_noumea_vs_voh.txt" 
input="input_res_poya_vs_voh.txt" 
# If FUNGI
taxon_Annotations="database_fungi_package_all.tab"
# If BACTERIA
taxon_Annotations="database_bacteria_package_all.tab"
#taxon_Database=file.path(original_dir, paste("Working_scripts/ncbitaxon_ontology.obo", sep="")) # download from http://www.obofoundry.org/ontology/ncbitaxon.html and then modified for being read by the package (follow the XXXX protocol for create your updated database).
#taxon_Database=system.file("extdata", "ncbitaxon_ontology.obo", package="Anaconda")
taxon_Database=file.path(original_dir, paste("Working_scripts/ncbitaxon_ontology.obo", sep="")) # download from http://www.obofoundry.org/ontology/ncbitaxon.html and then modified for being read by the package (follow the XXXX protocol for create your updated database).

taxon_Division="TR" # TR = taxonomic Rank, don't change this

#source(file.path(original_dir, paste("Working_scripts/Anaconda.functions.R", sep="")))

taxon_mwuStats_res <- taxon_mwuStats(input, taxon_Database, taxon_Annotations, taxon_Division,
               perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
               largest=0.1,  # a NCBITaxon terms will not be considered if it contains more than this fraction of the total number of ASVs #0.1
               smallest=1,   # a NCBITaxon terms contain at least this many ASVs to be considered #5
               clusterCutHeight=0.1 )

# Save results, if you want to redo downstream analysis without redo the previous step
#save(taxon_mwuStats_res, file = "input_res_forest_vs_short_fallowl_largest_0.1_smallest_1.rda")

# For loading results
#load("taxon_mwuStats_res_forest_vs_short_fallow_input_all_package_all_largest_0.1_smallest_1.rda")

dev.off()

taxon_mwuPlot(input,taxon_Annotations,taxon_Division,
                     absValue= -log(0.05,10),  # ASVs with the measure value exceeding this will be counted as "good ASVs". Specify absValue=0.5 if you are doing Fisher's exact test for standard NCBITaxon terms enrichment. # -log(0.05,10)
                     level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all NCBITaxon terms containing ASVs exceeding the absValue. # -log(0.05,10)
                     level2=0.05, # FDR cutoff to print in regular (not italic) font.
                     level3=0.01, # FDR cutoff to print in large bold font.
                     txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                     treeHeight=0.5 # height of the hierarchical clustering tree
                     #	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)

##################################################
# Save it to a PDF file
pdf(file = "Figure_01_taxon_mwuPlot_input_res_forest_vs_short_fallow.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 6)
taxon_mwuPlot(input,taxon_Annotations,taxon_Division,
              absValue= -log(0.05,10),  # ASVs with the measure value exceeding this will be counted as "good ASVs". Specify absValue=0.5 if you are doing Fisher's exact test for standard NCBITaxon terms enrichment. # -log(0.05,10)
              level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all NCBITaxon terms containing ASVs exceeding the absValue. # -log(0.05,10)
              level2=0.05, # FDR cutoff to print in regular (not italic) font.
              level3=0.01, # FDR cutoff to print in large bold font.
              txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
              treeHeight=0.5 # height of the hierarchical clustering tree
              #	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)
dev.off()
##################################################

####################################################################################################
# With guilds
####################################################################################################

#source(file.path(original_dir, paste("Working_scripts/Anaconda.functions.R", sep="")))

taxon_list <- taxon_mwu_list(input,taxon_Annotations,taxon_Division)#, absValue= -log(0.05,10))
taxon_list
taxon_list_drawer <- get_taxon_list_drawer(taxon_list)
taxon_list_drawer


# If you prefer make an input from your own file from your own Funguild analysis,
# Don't run funguilds <- get_funguilds(taxon_list_drawer)
# and run this
# link_guilds <- read.table("YOUR_FILE.txt", header=F, sep="\t")
# were OUR_FILE.txt is your own funguild results
# See examples files for help

funguilds <- get_funguilds(taxon_list_drawer)
funguilds

link_guilds <- get_link_guilds(taxon_list, funguilds)
link_guilds


taxon_mwuPlot_guilds(input,taxon_Annotations,taxon_Division,
                     absValue= -log(0.05,10),  # ASVs with the measure value exceeding this will be counted as "good ASVs". Specify absValue=0.5 if you are doing Fisher's exact test for standard NCBITaxon terms enrichment. # -log(0.05,10)
                     level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all NCBITaxon terms containing ASVs exceeding the absValue. # -log(0.05,10)
                     level2=0.05, # FDR cutoff to print in regular (not italic) font.
                     level3=0.01, # FDR cutoff to print in large bold font.
                     txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                     treeHeight=0.5 # height of the hierarchical clustering tree
                     #	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)

##################################################
# Save it to a PDF file
pdf(file = "Figure_02_taxon_mwuPlot_guilds_res_forest_vs_short_fallow.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 6)
taxon_mwuPlot_guilds(input,taxon_Annotations,taxon_Division,
                     absValue= -log(0.05,10),  # ASVs with the measure value exceeding this will be counted as "good ASVs". Specify absValue=0.5 if you are doing Fisher's exact test for standard NCBITaxon terms enrichment. # -log(0.05,10)
                     level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all NCBITaxon terms containing ASVs exceeding the absValue. # -log(0.05,10)
                     level2=0.05, # FDR cutoff to print in regular (not italic) font.
                     level3=0.01, # FDR cutoff to print in large bold font.
                     txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                     treeHeight=0.5 # height of the hierarchical clustering tree
                     #	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)
dev.off()
##################################################

