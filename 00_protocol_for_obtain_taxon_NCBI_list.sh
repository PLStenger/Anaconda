# ncbitaxon.owl and ncbitaxon.obo files
# from http://www.obofoundry.org/ontology/ncbitaxon.html
# Find from publication : https://www.sciencedirect.com/science/article/pii/S1532046407000809

###################################################################################################################################
# file ncbitaxon.obo 
###################################################################################################################################

mv ncbitaxon.obo > ncbitaxon.txt

sed -e ':a' -e 'N' -e '$!ba' -e 's/\nname:/XXXXname:/g' /Users/pierre-louisstenger/Desktop/ncbitaxon.txt > /Users/pierre-louisstenger/Desktop/ncbitaxon_02.txt

grep "id: NCBITaxon" /Users/pierre-louisstenger/Desktop/ncbitaxon_02.txt > /Users/pierre-louisstenger/Desktop/ncbitaxon_03.txt

# grep -c "id: " /Users/pierre-louisstenger/Desktop/ncbitaxon_03.txt 
# 2398056
# 2 398 056

sed $'s/XXXX/\t/g' /Users/pierre-louisstenger/Desktop/ncbitaxon_03.txt > /Users/pierre-louisstenger/Desktop/ncbitaxon_04.txt
sed $'s/id: //g' /Users/pierre-louisstenger/Desktop/ncbitaxon_04.txt > /Users/pierre-louisstenger/Desktop/ncbitaxon_05.txt
sed $'s/name: //g' /Users/pierre-louisstenger/Desktop/ncbitaxon_05.txt > /Users/pierre-louisstenger/Desktop/ncbitaxon_06.txt

#awk '{print $2"\t"$1}' /Users/pierre-louisstenger/Desktop/ncbitaxon_06.txt > /Users/pierre-louisstenger/Desktop/ncbitaxon_list.txt 

awk '{first = $1; $1=""; print $0, first}' /Users/pierre-louisstenger/Desktop/ncbitaxon_06.txt > /Users/pierre-louisstenger/Desktop/ncbitaxon_list.txt 

sed $'s/ NCBITaxon:/\tNCBITaxon:/g' /Users/pierre-louisstenger/Desktop/ncbitaxon_list.txt > /Users/pierre-louisstenger/Desktop/ncbitaxon_list_02.txt

sed 's/ //' /Users/pierre-louisstenger/Desktop/ncbitaxon_list_02.txt > /Users/pierre-louisstenger/Desktop/ncbitaxon_list_03.txt

mv /Users/pierre-louisstenger/Desktop/ncbitaxon_list_03.txt > /Users/pierre-louisstenger/Desktop/ncbitaxon_list.txt

###################################################################################################################################
# file taxonomy.tsv (from RefTaxo folder) 
###################################################################################################################################

mv taxonomy.tsv > taxonomy_all_fungi.txt

# taxonomy_all_fungi.txt to taxonomy_all_fungi_02.txt by deleting with search tool all "xx__unasigned".
# Then delete first line

awk '{print $NF}' FS=__ /Users/pierre-louisstenger/Desktop/taxonomy_all_fungi_02.txt > /Users/pierre-louisstenger/Desktop/taxonomy_all_fungi_03.txt
sed $'s/"//g' /Users/pierre-louisstenger/Desktop/taxonomy_all_fungi_03.txt > /Users/pierre-louisstenger/Desktop/taxonomy_all_fungi_04.txt

paste -d '\t' /Users/pierre-louisstenger/Desktop/taxonomy_all_fungi_02.txt /Users/pierre-louisstenger/Desktop/taxonomy_all_fungi_04.txt > /Users/pierre-louisstenger/Desktop/taxonomy_all_fungi_QIIME2_and_NCBI_format.txt 

# Deal with the "carriage return" problem with textMate and Excel edition
# Deleted all"_" in the NCBI corresponding column (edition with excel for seleciton of the column).

# Ok for the vlookup !
# /Users/pierre-louisstenger/Desktop/taxonomy_all_fungi_QIIME2_and_NCBI_format.txt 

###################################################################################################################################
# Vlookup script
###################################################################################################################################


#!/usr/bin/env line

awk '
    # { sub(/\r$/,"") }    # uncomment to remove Windows style line-endings.
    NR==FNR{a[$1]          # hash $1 of genes file to a
    next
}
($1 in a) {                # lookup from transcriptome
    print
}' /Users/pierre-louisstenger/Desktop/taxonomy_all_fungi_QIIME2_and_NCBI_format.txt /Users/pierre-louisstenger/Desktop/ncbitaxon_list.txt > /Users/pierre-louisstenger/Desktop/taxonomy_all_fungi_QIIME2_and_NCBI_format_final.txt


###################################################################################################################################
# Vlookup excel between taxonomy_all_fungi_QIIME2_and_NCBI_format_final.txt and /Users/pierre-louisstenger/Desktop/taxonomy_all_fungi_QIIME2_and_NCBI_format.txt
###################################################################################################################################

# --> Put the thrid column into the taxonomy_all_fungi_QIIME2_and_NCBI_format.txt file in order to didn't expect duplications

# Then change into database_fungi_RGBOA.tab


awk '{print $1"\t"$2}' /Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/05_Mare_ignames/Diversity_in_Mare_yam_crop/09_RBGOA/database_fungi_RGBOA_02.tab > /Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/05_Mare_ignames/Diversity_in_Mare_yam_crop/09_RBGOA/database_fungi_RGBOA_03.tab 

# Delete duplicates
awk '!seen[$0]++' /Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/05_Mare_ignames/Diversity_in_Mare_yam_crop/09_RBGOA/database_fungi_RGBOA_03.tab > /Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/05_Mare_ignames/Diversity_in_Mare_yam_crop/09_RBGOA/database_fungi_RGBOA_04.tab


sed $'s/NCBITaxon/go/g' /Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/05_Mare_ignames/Diversity_in_Mare_yam_crop/09_RBGOA/ncbitaxon.obo > /Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/05_Mare_ignames/Diversity_in_Mare_yam_crop/09_RBGOA/ncbitaxon_GO.obo




sed 's/ncbitaxon/go/g' /Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/05_Mare_ignames/Diversity_in_Mare_yam_crop/09_RBGOA/ncbitaxon_GO.obo > /Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/05_Mare_ignames/Diversity_in_Mare_yam_crop/09_RBGOA/ncbitaxon_GO_02.obo

sed 's/ncbi_taxonomy/biological_process/g' /Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/05_Mare_ignames/Diversity_in_Mare_yam_crop/09_RBGOA/ncbitaxon_GO_02.obo > /Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/05_Mare_ignames/Diversity_in_Mare_yam_crop/09_RBGOA/ncbitaxon_GO_03.obo

sed 's/id: ontology: go:/id: GO:/g' /Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/05_Mare_ignames/Diversity_in_Mare_yam_crop/09_RBGOA/ncbitaxon_GO_03.obo > /Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/05_Mare_ignames/Diversity_in_Mare_yam_crop/09_RBGOA/ncbitaxon_GO_04.obo


sed 's/is_a: ontology: go:/is_a: GO:/g' /Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/05_Mare_ignames/Diversity_in_Mare_yam_crop/09_RBGOA/ncbitaxon_GO_04.obo > /Users/pierre-louisstenger/Documents/PostDoc_02_MetaBarcoding_IAC/02_Data/05_Mare_ignames/Diversity_in_Mare_yam_crop/09_RBGOA/ncbitaxon_GO_05.obo







> gomwuStats(input, goDatabase, goAnnotations, goDivision,
+            perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
+            largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes #0.1
+            smallest=2,   # a GO category should contain at least this many genes to be considered #5
+            clusterCutHeight=0.1 )
ncbitaxon_GO_05.obo database_fungi_RGBOA_06.tab res_forest_vs_short_fallow_taxo_ok_04.txt BP largest=0.1 smallest=2cutHeight=0.1

Run parameters:

largest GO category as fraction of all genes (largest)  : 0.1
         smallest GO category as # of genes (smallest)  : 2
                clustering threshold (clusterCutHeight) : 0.1

-----------------
retrieving GO hierarchy, reformatting data...

-------------
go_reformat:
Genes with GO annotations, but not listed in measure table: 0

Terms without defined level (old ontology?..): 0
-------------
-------------
go_nrify:
21452 categories, 18639 genes; size range 2-1863.9
	22 too broad
	18450 too small
	2980 remaining

removing redundancy:

calculating GO term similarities based on shared genes...
2662 non-redundant GO categories of good size
-------------

Secondary clustering:
calculating similarities....
Continuous measure of interest: will perform MWU test
52  GO terms at 10% FDR
Warning message:
In for (i in (1L:cols)[do]) { :
  closing unused connection 3 (dissim0_BP_database_fungi_RGBOA_06.tab)
  
sed 's/biological_process/taxonomy_rank/g' /Users/pierre-louisstenger/Desktop/09_RBGOA/ncbitaxon_GO_05.obo > /Users/pierre-louisstenger/Desktop/09_RBGOA/ncbitaxon_GO_08.obo







 setwd("~/Desktop/09_RBGOA")
> # Edit these to match your data file names: 
> input="res_forest_vs_short_fallow_taxo_ok_04.txt" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
> goAnnotations="database_fungi_RGBOA_06.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
> goDatabase="ncbitaxon_GO_08.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
> #goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
> 
> goDivision="TR" # either MF, or BP, or CC
> source("gomwu.functions_03.R")
> 
> gomwuStats(input, goDatabase, goAnnotations, goDivision,
+            perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
+            largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes #0.1
+            smallest=2,   # a GO category should contain at least this many genes to be considered #5
+            clusterCutHeight=0.1 )
ncbitaxon_GO_08.obo database_fungi_RGBOA_06.tab res_forest_vs_short_fallow_taxo_ok_04.txt TR largest=0.1 smallest=2cutHeight=0.1

Run parameters:

largest NCBITaxon category as fraction of all ASVs (largest)  : 0.1
         smallest NCBITaxon category as # of ASVs (smallest)  : 2
                clustering threshold (clusterCutHeight) : 0.1

-----------------
retrieving NCBITaxon terms hierarchy, reformatting data...

-------------
go_reformat:
ASVs with NCBITaxon terms annotations, but not listed in measure table: 0

Terms without defined level (old ontology?..): 0
-------------
-------------
go_nrify:
21452 NCBITaxon terms, 18639 ASVs; size range 2-1863.9
	22 too broad
	18450 too small
	2980 remaining

removing redundancy:

calculating NCBITaxon term similarities based on shared ASVs...
2662 non-redundant NCBITaxon terms of good size
-------------

Secondary clustering:
calculating similarities....
Continuous measure of interest: will perform MWU test
52  NCBITaxon terms at 10% FDR
  
  
sed 's/GO:/NCBITaxon:/g' /Users/pierre-louisstenger/Desktop/09_RBGOA/ncbitaxon_GO_08.obo > /Users/pierre-louisstenger/Desktop/09_RBGOA/ncbitaxon_GO_09.obo

sed 's/GO:/NCBITaxon:/g' /Users/pierre-louisstenger/Desktop/09_RBGOA/database_fungi_RGBOA_06.tab > /Users/pierre-louisstenger/Desktop/09_RBGOA/database_fungi_RGBOA_07.tab  
  
---> OK !!!


> ######################################################################################################
> ######################################################################################################
> setwd("~/Desktop/09_RBGOA")
> # Edit these to match your data file names: 
> input="res_forest_vs_short_fallow_taxo_ok_04.txt" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
> taxon_Annotations="database_fungi_RGBOA_07.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
> taxon_Database="ncbitaxon_GO_09.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
> taxon_Division="TR" # either MF, or BP, or CC
> source("gomwu.functions_04.R")
> taxon_mwuStats(input, taxon_Database, taxon_Annotations, taxon_Division,
+            perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
+            largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes #0.1
+            smallest=2,   # a GO category should contain at least this many genes to be considered #5
+            clusterCutHeight=0.1 )
ncbitaxon_GO_09.obo database_fungi_RGBOA_07.tab res_forest_vs_short_fallow_taxo_ok_04.txt TR largest=0.1 smallest=2cutHeight=0.1

Run parameters:

largest NCBITaxon category as fraction of all ASVs (largest)  : 0.1
         smallest NCBITaxon category as number of ASVs (smallest)  : 2
                clustering threshold (clusterCutHeight) : 0.1

-----------------
retrieving NCBITaxon terms hierarchy, reformatting data...

-------------
go_reformat:
ASVs with NCBITaxon terms annotations, but not listed in measure table: 0

Terms without defined level (old ontology?..): 0
-------------
-------------
go_nrify:
21452 NCBITaxon terms, 18639 ASVs; size range 2-1863.9
	22 too broad
	18450 too small
	2980 remaining

removing redundancy:

calculating NCBITaxon term similarities based on shared ASVs...
2662 non-redundant NCBITaxon terms of good size
-------------

Secondary clustering:
calculating similarities....
Continuous measure of interest: will perform MWU test
52  NCBITaxon terms at 10% FDR
  
  
  
  
> taxon_mwuPlot(input,taxon_Annotations,taxon_Division,
+           absValue= -log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.5 if you are doing Fisher's exact test for standard GO enrichment. # -log(0.05,10)
+           level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. # -log(0.05,10)
+           level2=0.05, # FDR cutoff to print in regular (not italic) font.
+           level3=0.01, # FDR cutoff to print in large bold font.
+           txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
+           treeHeight=0.5 # height of the hierarchical clustering tree
+           #	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
+ )
NCBITaxon terms dispayed:  52 
"Good ASVs" accounted for:  9 out of 13 ( 69% )
                                      pval direction       color
Clitopilus                    4.832016e-05         1  firebrick1
Lasiodiplodia                 3.110810e-02         0 dodgerblue2
Cystolepiota                  4.832016e-05         1  firebrick1
Exophiala                     1.520382e-05         0 dodgerblue2
Metarhizium                   6.908487e-03         0 dodgerblue2
Lycoperdon                    2.247101e-02         0 dodgerblue2
Exobasidiaceae                5.507175e-02         1  lightcoral
Ceratobasidiaceae             4.270245e-03         0 dodgerblue2
Neurospora                    1.000000e-15         0 dodgerblue2
Sordariaceae                  9.505784e-05         0 dodgerblue2
Veronaeopsis                  1.000000e-15         0 dodgerblue2
Trechisporales                4.691952e-02         1  firebrick1
Botryosphaeria                7.175982e-06         0 dodgerblue2
Daldinia                      1.031289e-02         0 dodgerblue2
Minimedusa                    1.000000e-15         0 dodgerblue2
Staphylotrichum               1.001144e-15         0 dodgerblue2
Veronaea                      3.478299e-14         1  firebrick1
Amauroderma                   6.507747e-04         0 dodgerblue2
Sarocladiaceae                2.329664e-03         0 dodgerblue2
Trichoglossum                 5.376352e-08         1  firebrick1
Auriculariales incertae sedis 2.338168e-09         1  firebrick1
Ochroconis                    5.433471e-02         0    skyblue2
Hawksworthiomyces             1.000000e-15         0 dodgerblue2
Amorosiaceae                  1.000000e-15         1  firebrick1
Tetracladium                  1.001187e-15         0 dodgerblue2
Dactylonectria                1.001190e-15         0 dodgerblue2
Sebacinaceae                  1.492247e-02         0 dodgerblue2
Sebacina                      2.944700e-05         0 dodgerblue2
Verruconis                    1.000793e-15         1  firebrick1
Volutella                     3.493030e-14         1  firebrick1
Chaetomella                   1.001119e-15         0 dodgerblue2
Chaetomellales                9.505784e-05         0 dodgerblue2
Acrocalymma                   1.000000e-15         0 dodgerblue2
Morosphaeriaceae              1.044724e-15         0 dodgerblue2
Wallemiomycotina              1.308824e-03         0 dodgerblue2
Basidioascus                  1.000000e-15         0 dodgerblue2
Geminibasidiales              5.574759e-11         0 dodgerblue2
Spizellomyces                 1.068433e-07         0 dodgerblue2
Spizellomycetaceae            2.568113e-04         0 dodgerblue2
Chytridiomycota               3.180827e-03         0 dodgerblue2
Spizellomycetales             2.329664e-03         0 dodgerblue2
Rhizophlyctidales             1.068433e-07         0 dodgerblue2
Rhizophlyctis                 1.000000e-15         0 dodgerblue2
Calcarisporiella              1.000000e-15         0 dodgerblue2
Saksenaea                     5.658529e-11         0 dodgerblue2
Mortierella                   1.071656e-05         0 dodgerblue2
Modicella                     1.000000e-15         1  firebrick1
Gongronella                   1.000793e-15         1  firebrick1
Fungi incertae sedis          1.614907e-03         0 dodgerblue2
Glomerales                    1.065358e-04         0 dodgerblue2
Claroideoglomus               1.068433e-07         0 dodgerblue2
Claroideoglomeraceae          1.181661e-06         0 dodgerblue2
Warning messages:
1: In plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab = "",  :
  the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'
2: In plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab = "",  :
  the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1' 
  
  
  
  
  
  
  

sed $'s/GO:/NCBITaxon:/g' /Users/pierre-louisstenger/Desktop/09_RBGOA/database_fungi_RGBOA_06.tab > /Users/pierre-louisstenger/Desktop/09_RBGOA/database_fungi_RGBOA_07.tab

sed 's/biological_process/taxonomy_rank/g' /Users/pierre-louisstenger/Desktop/09_RBGOA/ncbitaxon_GO_06.obo > /Users/pierre-louisstenger/Desktop/09_RBGOA/ncbitaxon_GO_07.obo











# UPDATE 31/08/21

> ######################################################################################################
> ######################################################################################################
> setwd("~/Desktop/09_RBGOA")
> # Edit these to match your data file names: 
> input="res_forest_vs_short_fallow_taxo_ok_04.txt" # two columns of comma-separated values: asv id, continuous measure of significance. To perform standard NCBITaxon terms enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
> taxon_Annotations="database_fungi_RGBOA_07.tab" # two-column, tab-delimited, one line per gene, multiple NCBITaxon terms separated by semicolon.
> taxon_Database="ncbitaxon_GO_09.obo" # download from http://www.obofoundry.org/ontology/ncbitaxon.html and then modified for being read by the package (follow the XXXX protocol for create your updated database).
> taxon_Division="TR" # TR = taxonomic Rank
> 
> source("NCBITaxon_mwu.functions.R")
> 
> taxon_mwuStats(input, taxon_Database, taxon_Annotations, taxon_Division,
+            perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
+            largest=0.1,  # a NCBITaxon terms will not be considered if it contains more than this fraction of the total number of ASVs #0.1
+            smallest=2,   # a NCBITaxon terms contain at least this many ASVs to be considered #5
+            clusterCutHeight=0.1 )
ncbitaxon_GO_09.obo database_fungi_RGBOA_07.tab res_forest_vs_short_fallow_taxo_ok_04.txt TR largest=0.1 smallest=2cutHeight=0.1

Run parameters:

largest NCBITaxon category as fraction of all ASVs (largest)  : 0.1
         smallest NCBITaxon category as number of ASVs (smallest)  : 2
                clustering threshold (clusterCutHeight) : 0.1

-----------------
retrieving NCBITaxon terms hierarchy, reformatting data...

-------------
go_reformat:
ASVs with NCBITaxon terms annotations, but not listed in measure table: 0

Terms without defined level (old ontology?..): 0
-------------
-------------
go_nrify:
21452 NCBITaxon terms, 18639 ASVs; size range 2-1863.9
	22 too broad
	18450 too small
	2980 remaining

removing redundancy:

calculating NCBITaxon term similarities based on shared ASVs...
2662 non-redundant NCBITaxon terms of good size
-------------

Secondary clustering:
calculating similarities....
Continuous measure of interest: will perform MWU test
52  NCBITaxon terms at 10% FDR




> taxon_mwuPlot(input,taxon_Annotations,taxon_Division,
+           absValue= -log(0.05,10),  # ASVs with the measure value exceeding this will be counted as "good ASVs". Specify absValue=0.5 if you are doing Fisher's exact test for standard NCBITaxon terms enrichment. # -log(0.05,10)
+           level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all NCBITaxon terms containing ASVs exceeding the absValue. # -log(0.05,10)
+           level2=0.05, # FDR cutoff to print in regular (not italic) font.
+           level3=0.01, # FDR cutoff to print in large bold font.
+           txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
+           treeHeight=0.5 # height of the hierarchical clustering tree
+           #	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
+ )
NCBITaxon terms dispayed:  52 
"Good ASVs" accounted for:  9 out of 13 ( 69% )
                                      pval direction       color
Clitopilus                    4.832016e-05         1  firebrick1
Lasiodiplodia                 3.110810e-02         0 dodgerblue2
Cystolepiota                  4.832016e-05         1  firebrick1
Exophiala                     1.520382e-05         0 dodgerblue2
Metarhizium                   6.908487e-03         0 dodgerblue2
Lycoperdon                    2.247101e-02         0 dodgerblue2
Exobasidiaceae                5.507175e-02         1  lightcoral
Ceratobasidiaceae             4.270245e-03         0 dodgerblue2
Neurospora                    1.000000e-15         0 dodgerblue2
Sordariaceae                  9.505784e-05         0 dodgerblue2
Veronaeopsis                  1.000000e-15         0 dodgerblue2
Trechisporales                4.691952e-02         1  firebrick1
Botryosphaeria                7.175982e-06         0 dodgerblue2
Daldinia                      1.031289e-02         0 dodgerblue2
Minimedusa                    1.000000e-15         0 dodgerblue2
Staphylotrichum               1.001144e-15         0 dodgerblue2
Veronaea                      3.478299e-14         1  firebrick1
Amauroderma                   6.507747e-04         0 dodgerblue2
Sarocladium                   2.329664e-03         0 dodgerblue2
Trichoglossum                 5.376352e-08         1  firebrick1
Auriculariales incertae sedis 2.338168e-09         1  firebrick1
Ochroconis                    5.433471e-02         0    skyblue2
Hawksworthiomyces             1.000000e-15         0 dodgerblue2
Amorosiaceae                  1.000000e-15         1  firebrick1
Tetracladium                  1.001187e-15         0 dodgerblue2
Dactylonectria                1.001190e-15         0 dodgerblue2
Sebacinaceae                  1.492247e-02         0 dodgerblue2
Sebacina                      2.944700e-05         0 dodgerblue2
Verruconis                    1.000793e-15         1  firebrick1
Volutella                     3.493030e-14         1  firebrick1
Chaetomella                   1.001119e-15         0 dodgerblue2
Chaetomellales                9.505784e-05         0 dodgerblue2
Acrocalymma                   1.000000e-15         0 dodgerblue2
Morosphaeriaceae              1.044724e-15         0 dodgerblue2
Wallemiomycotina              1.308824e-03         0 dodgerblue2
Basidioascus                  1.000000e-15         0 dodgerblue2
Geminibasidiaceae             5.574759e-11         0 dodgerblue2
Spizellomyces                 1.068433e-07         0 dodgerblue2
Spizellomycetaceae            2.568113e-04         0 dodgerblue2
Chytridiomycota               3.180827e-03         0 dodgerblue2
Spizellomycetales             2.329664e-03         0 dodgerblue2
Rhizophlyctidales             1.068433e-07         0 dodgerblue2
Rhizophlyctis                 1.000000e-15         0 dodgerblue2
Calcarisporiellales           1.000000e-15         0 dodgerblue2
Saksenaeaceae                 5.658529e-11         0 dodgerblue2
Mortierella                   1.071656e-05         0 dodgerblue2
Modicella                     1.000000e-15         1  firebrick1
Gongronella                   1.000793e-15         1  firebrick1
Fungi incertae sedis          1.614907e-03         0 dodgerblue2
Glomerales                    1.065358e-04         0 dodgerblue2
Claroideoglomus               1.068433e-07         0 dodgerblue2
Claroideoglomeraceae          1.181661e-06         0 dodgerblue2
Warning messages:
1: In plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab = "",  :
  the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'
2: In plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab = "",  :
  the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'
  
  
  
  
  
  
######################### AVEC LES ASVS AU LIEU DE function LOC 



 # Edit these to match your data file names: 
> input="res_forest_vs_short_fallow_taxo_ok_04_02.txt" # two columns of comma-separated values: asv id, continuous measure of significance. To perform standard NCBITaxon terms enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
> taxon_Annotations="database_fungi_RGBOA_07_02.tab" # two-column, tab-delimited, one line per gene, multiple NCBITaxon terms separated by semicolon.
> taxon_Database="ncbitaxon_GO_09.obo" # download from http://www.obofoundry.org/ontology/ncbitaxon.html and then modified for being read by the package (follow the XXXX protocol for create your updated database).
> taxon_Division="TR" # TR = taxonomic Rank
> 
> source("NCBITaxon_mwu.functions.R")
> 
> taxon_mwuStats(input, taxon_Database, taxon_Annotations, taxon_Division,
+                perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
+                largest=0.1,  # a NCBITaxon terms will not be considered if it contains more than this fraction of the total number of ASVs #0.1
+                smallest=2,   # a NCBITaxon terms contain at least this many ASVs to be considered #5
+                clusterCutHeight=0.1 )
ncbitaxon_GO_09.obo database_fungi_RGBOA_07_02.tab res_forest_vs_short_fallow_taxo_ok_04_02.txt TR largest=0.1 smallest=2cutHeight=0.1

Run parameters:

largest NCBITaxon category as fraction of all ASVs (largest)  : 0.1
         smallest NCBITaxon category as number of ASVs (smallest)  : 2
                clustering threshold (clusterCutHeight) : 0.1

-----------------
retrieving NCBITaxon terms hierarchy, reformatting data...

-------------
go_reformat:
ASVs with NCBITaxon terms annotations, but not listed in measure table: 0

Terms without defined level (old ontology?..): 0
-------------
-------------
go_nrify:
21452 NCBITaxon terms, 18639 ASVs; size range 2-1863.9
	22 too broad
	18450 too small
	2980 remaining

removing redundancy:

calculating NCBITaxon term similarities based on shared ASVs...
2662 non-redundant NCBITaxon terms of good size
-------------

Secondary clustering:
calculating similarities....
Continuous measure of interest: will perform MWU test
52  NCBITaxon terms at 10% FDR
> library(ape)
> 
> # Plotting results
> #quartz()
> taxon_mwuPlot(input,taxon_Annotations,taxon_Division,
+           absValue= -log(0.05,10),  # ASVs with the measure value exceeding this will be counted as "good ASVs". Specify absValue=0.5 if you are doing Fisher's exact test for standard NCBITaxon terms enrichment. # -log(0.05,10)
+           level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all NCBITaxon terms containing ASVs exceeding the absValue. # -log(0.05,10)
+           level2=0.05, # FDR cutoff to print in regular (not italic) font.
+           level3=0.01, # FDR cutoff to print in large bold font.
+           txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
+           treeHeight=0.5 # height of the hierarchical clustering tree
+           #	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
+ )
NCBITaxon terms dispayed:  52 
"Good ASVs" accounted for:  9 out of 13 ( 69% )
                                      pval direction       color
Clitopilus                    4.832016e-05         1  firebrick1
Lasiodiplodia                 3.110810e-02         0 dodgerblue2
Cystolepiota                  4.832016e-05         1  firebrick1
Exophiala                     1.520382e-05         0 dodgerblue2
Metarhizium                   6.908487e-03         0 dodgerblue2
Lycoperdon                    2.247101e-02         0 dodgerblue2
Exobasidiaceae                5.507175e-02         1  lightcoral
Ceratobasidiaceae             4.270245e-03         0 dodgerblue2
Neurospora                    1.000000e-15         0 dodgerblue2
Sordariaceae                  9.505784e-05         0 dodgerblue2
Veronaeopsis                  1.000000e-15         0 dodgerblue2
Trechisporales                4.691952e-02         1  firebrick1
Botryosphaeria                7.175982e-06         0 dodgerblue2
Daldinia                      1.031289e-02         0 dodgerblue2
Minimedusa                    1.000000e-15         0 dodgerblue2
Staphylotrichum               1.001144e-15         0 dodgerblue2
Veronaea                      3.478299e-14         1  firebrick1
Amauroderma                   6.507747e-04         0 dodgerblue2
Sarocladium                   2.329664e-03         0 dodgerblue2
Trichoglossum                 5.376352e-08         1  firebrick1
Auriculariales incertae sedis 2.338168e-09         1  firebrick1
Ochroconis                    5.433471e-02         0    skyblue2
Hawksworthiomyces             1.000000e-15         0 dodgerblue2
Amorosiaceae                  1.000000e-15         1  firebrick1
Tetracladium                  1.001187e-15         0 dodgerblue2
Dactylonectria                1.001190e-15         0 dodgerblue2
Sebacinaceae                  1.492247e-02         0 dodgerblue2
Sebacina                      2.944700e-05         0 dodgerblue2
Verruconis                    1.000793e-15         1  firebrick1
Volutella                     3.493030e-14         1  firebrick1
Chaetomella                   1.001119e-15         0 dodgerblue2
Chaetomellaceae               9.505784e-05         0 dodgerblue2
Acrocalymma                   1.000000e-15         0 dodgerblue2
Morosphaeriaceae              1.044724e-15         0 dodgerblue2
Wallemiomycotina              1.308824e-03         0 dodgerblue2
Basidioascus                  1.000000e-15         0 dodgerblue2
Geminibasidiaceae             5.574759e-11         0 dodgerblue2
Spizellomyces                 1.068433e-07         0 dodgerblue2
Spizellomycetaceae            2.568113e-04         0 dodgerblue2
Chytridiomycota               3.180827e-03         0 dodgerblue2
Spizellomycetales             2.329664e-03         0 dodgerblue2
Rhizophlyctidales             1.068433e-07         0 dodgerblue2
Rhizophlyctis                 1.000000e-15         0 dodgerblue2
Calcarisporiella              1.000000e-15         0 dodgerblue2
Saksenaea                     5.658529e-11         0 dodgerblue2
Mortierella                   1.071656e-05         0 dodgerblue2
Modicella                     1.000000e-15         1  firebrick1
Gongronella                   1.000793e-15         1  firebrick1
Fungi incertae sedis          1.614907e-03         0 dodgerblue2
Glomerales                    1.065358e-04         0 dodgerblue2
Claroideoglomus               1.068433e-07         0 dodgerblue2
Claroideoglomeraceae          1.181661e-06         0 dodgerblue2
Warning messages:
1: In plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab = "",  :
  the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'
2: In plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab = "",  :
  the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'
> # manually rescale the plot so the tree matches the text 
> # if there are too many categories displayed, try make it more stringent with level1=0.01,level2=0.005,level3=0.001.  















> # Test 01/09/2021
> 
> # Edit these to match your data file names: 
> input="res_forest_vs_short_fallow_taxo_ok_04_02.txt" # two columns of comma-separated values: asv id, continuous measure of significance. To perform standard NCBITaxon terms enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
> taxon_Annotations="database_fungi_RGBOA_07_02.tab" # two-column, tab-delimited, one line per gene, multiple NCBITaxon terms separated by semicolon.
> taxon_Database="ncbitaxon_GO_09.obo" # download from http://www.obofoundry.org/ontology/ncbitaxon.html and then modified for being read by the package (follow the XXXX protocol for create your updated database).
> taxon_Division="TR" # TR = taxonomic Rank
> 
> source("NCBITaxon_mwu.functions.R")
> 
> taxon_mwuStats(input, taxon_Database, taxon_Annotations, taxon_Division,
+                perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
+                largest=0.1,  # a NCBITaxon terms will not be considered if it contains more than this fraction of the total number of ASVs #0.1
+                smallest=2,   # a NCBITaxon terms contain at least this many ASVs to be considered #5
+                clusterCutHeight=0.1 )
ncbitaxon_GO_09.obo database_fungi_RGBOA_07_02.tab res_forest_vs_short_fallow_taxo_ok_04_02.txt TR largest=0.1 smallest=2cutHeight=0.1

Run parameters:

largest NCBITaxon category as fraction of all ASVs (largest)  : 0.1
         smallest NCBITaxon category as number of ASVs (smallest)  : 2
                clustering threshold (clusterCutHeight) : 0.1

-----------------
retrieving NCBITaxon terms hierarchy, reformatting data...

-------------
taxon_reformat:
ASVs with NCBITaxon terms annotations, but not listed in measure table: 0

Terms without defined level (old ontology?..): 0
-------------
-------------
taxon_nrify:
21452 NCBITaxon terms, 18639 ASVs; size range 2-1863.9
	22 too broad
	18450 too small
	2980 remaining

removing redundancy:

calculating NCBITaxon term similarities based on shared ASVs...
2662 non-redundant NCBITaxon terms of good size
-------------

Secondary clustering:
calculating similarities....
Continuous measure of interest: will perform MWU test
52  NCBITaxon terms at 10% FDR
> library(ape)
> 
> # Plotting results
> #quartz()
> taxon_mwuPlot(input,taxon_Annotations,taxon_Division,
+               absValue= -log(0.05,10),  # ASVs with the measure value exceeding this will be counted as "good ASVs". Specify absValue=0.5 if you are doing Fisher's exact test for standard NCBITaxon terms enrichment. # -log(0.05,10)
+               level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all NCBITaxon terms containing ASVs exceeding the absValue. # -log(0.05,10)
+               level2=0.05, # FDR cutoff to print in regular (not italic) font.
+               level3=0.01, # FDR cutoff to print in large bold font.
+               txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
+               treeHeight=0.5 # height of the hierarchical clustering tree
+               #	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
+ )
NCBITaxon terms dispayed:  52 
"Good ASVs" accounted for:  9 out of 13 ( 69% )
                                      pval direction       color
Clitopilus                    4.832016e-05         1  firebrick1
Lasiodiplodia                 3.110810e-02         0 dodgerblue2
Cystolepiota                  4.832016e-05         1  firebrick1
Exophiala                     1.520382e-05         0 dodgerblue2
Metarhizium                   6.908487e-03         0 dodgerblue2
Lycoperdon                    2.247101e-02         0 dodgerblue2
Exobasidiaceae                5.507175e-02         1  lightcoral
Ceratobasidiaceae             4.270245e-03         0 dodgerblue2
Neurospora                    1.000000e-15         0 dodgerblue2
Sordariaceae                  9.505784e-05         0 dodgerblue2
Veronaeopsis                  1.000000e-15         0 dodgerblue2
Trechisporales                4.691952e-02         1  firebrick1
Botryosphaeria                7.175982e-06         0 dodgerblue2
Daldinia                      1.031289e-02         0 dodgerblue2
Minimedusa                    1.000000e-15         0 dodgerblue2
Staphylotrichum               1.001144e-15         0 dodgerblue2
Veronaea                      3.478299e-14         1  firebrick1
Amauroderma                   6.507747e-04         0 dodgerblue2
Sarocladium                   2.329664e-03         0 dodgerblue2
Trichoglossum                 5.376352e-08         1  firebrick1
Auriculariales incertae sedis 2.338168e-09         1  firebrick1
Ochroconis                    5.433471e-02         0    skyblue2
Hawksworthiomyces             1.000000e-15         0 dodgerblue2
Amorosiaceae                  1.000000e-15         1  firebrick1
Tetracladium                  1.001187e-15         0 dodgerblue2
Dactylonectria                1.001190e-15         0 dodgerblue2
Sebacinaceae                  1.492247e-02         0 dodgerblue2
Sebacina                      2.944700e-05         0 dodgerblue2
Verruconis                    1.000793e-15         1  firebrick1
Volutella                     3.493030e-14         1  firebrick1
Chaetomella                   1.001119e-15         0 dodgerblue2
Chaetomellales                9.505784e-05         0 dodgerblue2
Acrocalymma                   1.000000e-15         0 dodgerblue2
Morosphaeriaceae              1.044724e-15         0 dodgerblue2
Wallemiomycotina              1.308824e-03         0 dodgerblue2
Basidioascus                  1.000000e-15         0 dodgerblue2
Geminibasidiaceae             5.574759e-11         0 dodgerblue2
Spizellomyces                 1.068433e-07         0 dodgerblue2
Spizellomycetaceae            2.568113e-04         0 dodgerblue2
Chytridiomycota               3.180827e-03         0 dodgerblue2
Spizellomycetales             2.329664e-03         0 dodgerblue2
Rhizophlyctidales             1.068433e-07         0 dodgerblue2
Rhizophlyctis                 1.000000e-15         0 dodgerblue2
Calcarisporiellaceae          1.000000e-15         0 dodgerblue2
Saksenaeaceae                 5.658529e-11         0 dodgerblue2
Mortierella                   1.071656e-05         0 dodgerblue2
Modicella                     1.000000e-15         1  firebrick1
Gongronella                   1.000793e-15         1  firebrick1
Fungi incertae sedis          1.614907e-03         0 dodgerblue2
Glomerales                    1.065358e-04         0 dodgerblue2
Claroideoglomus               1.068433e-07         0 dodgerblue2
Claroideoglomeraceae          1.181661e-06         0 dodgerblue2
Warning messages:
1: In plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab = "",  :
  the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'
2: In plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab = "",  :
  the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'
  
  