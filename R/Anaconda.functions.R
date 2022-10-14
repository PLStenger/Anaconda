#' @title database_metazoan_creation
#'
#' @description Create a Database for metazoan kingdom for Global analysis by Taxon_MWU analysis from targeted analysis. Please, run setwd("02_Global_analysis") after this function.
#'
#' @param nothing It's important not to write anything between the brackets, the database will create itself.
#'
#' @return A data frame file named database_metazoan_package_all.tab created from the taxonomy_all_metazoan_QIIME2_and_NCBI_format.txt file and your own taxonomy_RepSeq.tsv file.
#' database_metazoan_creation()
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' # It is important not to write anything between the brackets, the database will create itself.
#' \dontrun{database_metazoan_creation()}
#' # Please, run setwd("02_Global_analysis") after this function.
#' @export
database_metazoan_creation <- function (nothing)
{

  oldwd <- getwd()
  on.exit(setwd(oldwd))

  # Created sub directory "Targeted_analysis" if not already exist
  if (!dir.exists("02_Global_analysis")) {dir.create("02_Global_analysis")}
  setwd("02_Global_analysis")

  file.copy(from = file.path(kingdom, paste("taxonomy_RepSeq.tsv", sep="")), to ="taxonomy_RepSeq.tsv",
            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
  file.remove(file.path(kingdom, paste("taxonomy_RepSeq.tsv", sep="")))

  all_taxon <- read.table(file.path(original_dir, paste("Working_scripts/taxonomy_all_metazoan_QIIME2_and_NCBI_format.txt", sep="")), sep="\t")
  #all_taxon <- read.table(system.file("extdata", "taxonomy_all_metazoan_QIIME2_and_NCBI_format.txt", package="Anaconda"), sep="\t")

  taxo <- read.table("taxonomy_RepSeq.tsv", sep="\t", header=T)
  taxo_database <- merge(taxo,all_taxon, by.x="Taxon", by.y="V2", all.x=T)


  # Delete the duplicated rows that didn't have assigned if the duplicated ID have assignation
  i <- taxo_database$V3!=""
  taxo_database <- taxo_database[i | !taxo_database$Feature.ID %in% taxo_database$Feature.ID[i],]


  #taxo_database_df <- taxo_database[-which(duplicated(taxo_database)), ]
  taxon_Annotations <- data.frame(taxo_database$Feature.ID, taxo_database$V3)
  write.table(taxon_Annotations, file="database_metazoan_package_all.tab", sep="\t", col.names = F, row.names = F, quote = F)
}


#' @title database_bacteria_creation
#'
#' @description Create a Database for Bacteria kingdom for Global analysis by Taxon_MWU analysis from targeted analysis. Please, run setwd("02_Global_analysis") after this function.
#'
#' @param nothing It's important not to write anything between the brackets, the database will create itself.
#'
#' @return A data frame file named database_bacteria_package_all.tab created from the taxonomy_all_bacteria_QIIME2_and_NCBI_format.txt file and your own taxonomy_RepSeq.tsv file.
#' database_bacteria_creation()
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' # It is important not to write anything between the brackets, the database will create itself.
#' \dontrun{database_bacteria_creation()}
#' # Please, run setwd("02_Global_analysis") after this function.
#' @export
database_bacteria_creation <- function (nothing)
{

  oldwd <- getwd()
  on.exit(setwd(oldwd))

  # Created sub directory "Targeted_analysis" if not already exist
  if (!dir.exists("02_Global_analysis")) {dir.create("02_Global_analysis")}
  setwd("02_Global_analysis")

  file.copy(from = file.path(kingdom, paste("taxonomy_RepSeq.tsv", sep="")), to ="taxonomy_RepSeq.tsv",
            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
  file.remove(file.path(kingdom, paste("taxonomy_RepSeq.tsv", sep="")))

  all_taxon <- read.table(file.path(original_dir, paste("Working_scripts/taxonomy_all_bacteria_QIIME2_and_NCBI_format.txt", sep="")), sep="\t")
  #all_taxon <- read.table(system.file("extdata", "taxonomy_all_bacteria_QIIME2_and_NCBI_format.txt", package="Anaconda"), sep="\t")

  taxo <- read.table("taxonomy_RepSeq.tsv", sep="\t", header=T)
  taxo_database <- merge(taxo,all_taxon, by.x="Taxon", by.y="V2", all.x=T)


  # Delete the duplicated rows that didn't have assigned if the duplicated ID have assignation
  i <- taxo_database$V3!=""
  taxo_database <- taxo_database[i | !taxo_database$Feature.ID %in% taxo_database$Feature.ID[i],]


  #taxo_database_df <- taxo_database[-which(duplicated(taxo_database)), ]
  taxon_Annotations <- data.frame(taxo_database$Feature.ID, taxo_database$V3)
  write.table(taxon_Annotations, file="database_bacteria_package_all.tab", sep="\t", col.names = F, row.names = F, quote = F)
}

#' @title get_bactotraits_targeted
#'
#' @description Obtain Bacterial Traits for Bacteria kingdom for targeted analysis
#'
#' @param x Object from the Differential ASV abundance (DASVA) analysis
#'
#' @return A data frame file with Bacterial Traits informations for Bacteria kingdom for targeted analysis from the Differential ASV abundance (DASVA) analysis
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{get_bactotraits_targeted(res_forest_vs_long_fallow)}
#' @export
get_bactotraits_targeted <- function(x) {

  #file.copy(from = file.path(original_dir, paste("Working_scripts/BactoTraits_database.txt", sep="")), to ="BactoTraits_database.txt",
  #          overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

  file.copy(from = system.file("extdata", "BactoTraits_database.txt", package="Anaconda"), to ="BactoTraits_database.txt",
            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

  BactoTraits <- read.table("BactoTraits_database.txt", sep="\t", header=T)
  res_taxo <- merge(data.frame(x), taxo, by.x="row.names", by.y="Feature.ID")
  res_taxo_DT <- res_taxo
  res_taxo_02 <- setDT(res_taxo_DT)[, paste0("Taxon", 1:7) := tstrsplit(Taxon, ";")]
  res_taxo_03 <- data.frame(res_taxo_02$Taxon6, res_taxo_02$Row.names)
  colnames(res_taxo_03) <- c("Genus", "Row.names")
  res_taxo_03$Genus <- sub(" g__", "", res_taxo_03$Genus)
  BT <- data.frame(BactoTraits$Genus, BactoTraits$Heterotroph, BactoTraits$Autotroph, BactoTraits$Organotroph, BactoTraits$Lithotroph, BactoTraits$Chemotroph, BactoTraits$Phototroph, BactoTraits$Copiotroph_Diazotroph, BactoTraits$Methylotroph, BactoTraits$Oligotroph)
  colnames(BT) <- c("Genus", "Heterotroph", "Autotroph", "Organotroph", "Lithotroph", "Chemotroph", "Phototroph", "Copiotroph_Diazotroph", "Methylotroph", "Oligotroph")
  BT_link <- merge(res_taxo_03, BT, by.x="Genus", by.y="Genus", all.x=T)
  res_taxo_final <- merge(res_taxo, BT_link, by.x="Row.names", by.y="Row.names")
  unlink("BactoTraits_database.txt")

  return(res_taxo_final)

}

#' @title get_funguilds_targeted
#'
#' @description Obtain Fungi Guilds for Fungi kingdom for targeted analysis
#'
#' @param x Object from the funguild_input_targeted() output.
#'
#' @return A data frame file with Fungi Guilds informations for Fungi kingdom for targeted analysis from the Differential ASV abundance (DASVA) analysis
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{get_funguilds_targeted(res_forest_vs_long_fallow_guilds)}
#' @export
get_funguilds_targeted <- function(x) {

  # file.copy(from = file.path(original_dir, paste("Working_scripts/Guilds_v1.1.py", sep="")), to ="Guilds_v1.1.py",
  #            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

  file.copy(from = system.file("extdata", "Guilds_v1.1.py", package="Anaconda"), to ="Guilds_v1.1.py",
            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

  name = deparse(substitute(x))
  write.table(x, file=paste0("input_", name,"_targeted.txt", sep=""), sep="\t", row.names = F)
  latin = readLines(paste0("input_", name,"_targeted.txt", sep=""),-1)
  latin[1] = "ID\ttaxonomy"
  writeLines(latin,paste0("input_", name,"_targeted.txt", sep=""))

  # to deal with https://stackoverflow.com/questions/17309288/importerror-no-module-named-requests
  # if
  #FunGuild v1.1 Beta
  #Traceback (most recent call last):
  #  File "Guilds_v1.1.py", line 119, in <module>
  #  import requests
  #ImportError: No module named requests
  system('python -m pip install --user requests')

  # https://stackoverflow.com/questions/41638558/how-to-call-python-script-from-r-with-arguments
  #system('python Guilds_v1.1.py -otu taxon_list_drawer_input.txt -m -u')
  system(paste0('python Guilds_v1.1.py -otu input_',name,'_targeted.txt -m -u', sep=""))

  ## CLEAN taxon_list_drawer_input.guilds.txt
  # dat <- read.table(paste0('python Guilds_v1.1.py -otu input_',name,'_targeted.guilds.txt -m -u', sep=""), sep="\t", header=T)

  # https://stackoverflow.com/questions/26289681/r-regex-find-last-occurrence-of-delimiter
  # funguilds <- data.frame(sub(".*[__]", "", dat$taxonomy,),  dat$Guild, dat$Trophic.Mode,stringsAsFactors = FALSE)
  # colnames(funguilds) <- c("Taxon",  "Guild", "Trophic_Mode")


  delfiles <- dir(path=targeted_analysis_dir ,pattern="*matched.txt")
  file.remove(file.path(targeted_analysis_dir, delfiles))


  Files <- list.files(pattern="*.guilds.txt")
  newName <- sub("input_", "Table_02_03_", Files)
  newName_02 <- sub("guilds.", "", newName)
  file.rename(Files, newName_02)

  unlink("Guilds_v1.1.py")

}

#' @title funguild_input_targeted
#'
#' @description Prepare Object for Fungi Guilds for Fungi kingdom for targeted analysis
#'
#' @param x Object from the Differential ASV abundance (DASVA) analysis
#'
#' @return An Object used for Fungi Guilds informations for Fungi kingdom for targeted analysis from the Differential ASV abundance (DASVA) analysis
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{get_funguilds_targeted(res_forest_vs_long_fallow_guilds)}
#' @export
funguild_input_targeted <- function (x)
{
  taxon_for_guild <- data.frame(merge(data.frame(x), taxo, by.x="row.names", by.y="Feature.ID"))$Taxon
  taxon_for_guild_02 <- data.frame(rep(1:length(taxon_for_guild), 1), taxon_for_guild )
  colnames(taxon_for_guild_02) <- c("ID", "taxonomy")

  return(taxon_for_guild_02)
}

#' @title input_global_analysis
#'
#' @description Input files creation for each condition for Global analysis by Taxon_MWU analysis from targeted analysis (I)
#'
#' @param x Object from the Differential ASV abundance (DASVA) analysis
#'
#' @return Input Object for Global analysis by Taxon_MWU analysis from targeted analysis (I)
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{input_global_analysis(res_forest_vs_long_fallow)}
#' @export
input_global_analysis <- function (x)
{
  pre_input <- format_input(x)
  taxo_repseq <- read.table("taxonomy_RepSeq.tsv", sep="\t", header=T)
  pre_input_all <- merge(pre_input, taxo_repseq, by.x="ASV_ID", by.y="Feature.ID", all.y=T)
  pre_input_all$Taxon <- NULL
  pre_input_all$Confidence <- NULL
  pre_input_all$logP[is.na(pre_input_all$logP)] <- 0
  return(pre_input_all)
}

#' @title database_fungi_creation_RepSeq
#'
#' @description Create a Database for Fungi kingdom for Global analysis by Taxon_MWU analysis from targeted analysis. Please, run setwd("02_Global_analysis") after this function.
#'
#' @param nothing It's important not to write anything between the brackets, the database will create itself.
#'
#' @return A data frame file named database_fungi_package_all.tab created from the taxonomy_all_bacteria_QIIME2_and_NCBI_format.txt file and your own taxonomy_RepSeq.tsv file.
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' # It is important not to write anything between the brackets, the database will create itself.
#' \dontrun{database_fungi_creation_RepSeq()}
#' # Please, run setwd("02_Global_analysis") after this function.
#' @export
database_fungi_creation_RepSeq <- function (nothing)
{

  oldwd <- getwd()
  on.exit(setwd(oldwd))

  # Created sub directory "Targeted_analysis" if not already exist
  if (!dir.exists("02_Global_analysis")) {dir.create("02_Global_analysis")}
  setwd("02_Global_analysis")

  file.copy(from = file.path(kingdom, paste("taxonomy_RepSeq.tsv", sep="")), to ="taxonomy_RepSeq.tsv",
            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
  file.remove(file.path(kingdom, paste("taxonomy_RepSeq.tsv", sep="")))

  all_taxon <- read.table(file.path(original_dir, paste("Working_scripts/taxonomy_all_fungi_QIIME2_and_NCBI_format.txt", sep="")), sep="\t")
  #all_taxon <-  read.table(system.file("extdata", "taxonomy_all_fungi_QIIME2_and_NCBI_format.txt", package="Anaconda"), sep="\t")

  taxo <- read.table("taxonomy_RepSeq.tsv", sep="\t", header=T)
  taxo_database <- merge(taxo,all_taxon, by.x="Taxon", by.y="V2", all.x=T, all.y=F)
  taxo_database_df <- taxo_database[-which(duplicated(taxo_database)), ]
  taxon_Annotations <- data.frame(taxo_database_df$Feature.ID, taxo_database_df$V3)
  write.table(taxon_Annotations, file="database_fungi_package_all.tab", sep="\t", col.names = F, row.names = F, quote = F)
}

#' @title database_fungi_creation
#'
#' @description Create a Database for Fungi kingdom for Global analysis by Taxon_MWU analysis from targeted analysis only from rarefied ASVs. Please, run setwd("02_Global_analysis") after this function.
#'
#' @param nothing It's important not to write anything between the brackets, the database will create itself.
#'
#' @return A data frame file named database_fungi_package_all.tab created from the taxonomy_all_bacteria_QIIME2_and_NCBI_format.txt file and your own taxonomy.tsv file.
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' # It is important not to write anything between the brackets, the database will create itself.
#' \dontrun{database_fungi_creation()}
#' # Please, run setwd("02_Global_analysis") after this function.
#' @export
database_fungi_creation <- function (nothing)
{

  oldwd <- getwd()
  on.exit(setwd(oldwd))

  # Created sub directory "Targeted_analysis" if not already exist
  if (!dir.exists("02_Global_analysis")) {dir.create("02_Global_analysis")}
  setwd("02_Global_analysis")

  file.copy(from = file.path(kingdom, paste("taxonomy.tsv", sep="")), to ="taxonomy.tsv",
            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
  file.remove(file.path(kingdom, paste("taxonomy.tsv", sep="")))

  #all_taxon <- read.table(file.path(original_dir, paste("Working_scripts/taxonomy_all_fungi_QIIME2_and_NCBI_format.txt", sep="")), sep="\t")
  all_taxon <-  read.table(system.file("extdata", "taxonomy_all_fungi_QIIME2_and_NCBI_format.txt", package="Anaconda"), sep="\t")

  taxo <- read.table("taxonomy.tsv", sep="\t", header=T)
  taxo_database <- merge(taxo,all_taxon, by.x="Taxon", by.y="V2", all.x=T, all.y=F)
  taxo_database_df <- taxo_database[-which(duplicated(taxo_database)), ]
  taxon_Annotations <- data.frame(taxo_database_df$Feature.ID, taxo_database_df$V3)
  write.table(taxon_Annotations, file="database_fungi_package_all.tab", sep="\t", col.names = F, row.names = F, quote = F)
}

#' @title target_file
#'
#' @description Imports conditions information from your SampleSheet_comparison.txt file, with focus on iput files.
#'
#' @param nothing It's important not to write anything between the brackets, comparisons will create themselves.
#'
#' @return a data.frame with conditions information from your SampleSheet_comparison.txt file
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{target_file <- target_file()}
#' @export
target_file <- function (nothing)
{

  sampleSheet="SampleSheet_comparison.txt"

  samplesInfo <- read.table(sampleSheet, header=F) ; head(samplesInfo) ; dim(samplesInfo)
  target_file <- read.table(sampleSheet, header=F) ; head(samplesInfo) ; dim(samplesInfo)
  colnames(target_file) <- c("label","files","condition") ; head(target_file) ; dim(target_file)

  #threshold=1 # minimum of ASVs value across samples
  return(target_file)
}

#' @title samplesInfo
#'
#' @description Imports conditions information from your SampleSheet_comparison.txt file, with focus on samplesInfo.
#'
#' @param nothing It's important not to write anything between the brackets, comparisons will create themselves.
#'
#' @return a data.frame with conditions information from your SampleSheet_comparison.txt file, with focus on samplesInfo.
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{samplesInfo <- samplesInfo()}
#' @export
samplesInfo <- function (nothing)
{

  sampleSheet="SampleSheet_comparison.txt"

  samplesInfo <- read.table(sampleSheet, header=F) ; head(samplesInfo) ; dim(samplesInfo)
  target_file <- read.table(sampleSheet, header=F) ; head(samplesInfo) ; dim(samplesInfo)
  colnames(target_file) <- c("label","files","condition") ; head(target_file) ; dim(target_file)

  #threshold=1 # minimum of ASVs value across samples
  return(samplesInfo)
}

#' @title move_files
#'
#' @description Move the file in the good folders. Depending on the previous Kingdom selection (e.g., Fungi 'Fungi()', Bacteria 'Bacteria()', etc.)
#'
#' @param nothing It's important not to write anything between the brackets, files will move in the good folders, depending of the selected Kingdom before.
#'
#' @return Move the file in the good folders.
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{move_files()}
#' @export
move_files <- function (nothing)
{
  file.copy(from = file.path(original_dir, paste("taxonomy.tsv", sep="")), to ="taxonomy.tsv",
            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
  file.remove(file.path(original_dir, paste("taxonomy.tsv", sep="")))

  file.copy(from = file.path(original_dir, paste("taxonomy_RepSeq.tsv", sep="")), to ="taxonomy_RepSeq.tsv",
            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
  file.remove(file.path(original_dir, paste("taxonomy_RepSeq.tsv", sep="")))

  #file.copy(from = file.path(original_dir, paste("SampleSheet_comparison.txt", sep="")), to ="SampleSheet_comparison.txt",
  #          overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
  #file.remove(file.path(original_dir, paste("SampleSheet_comparison.txt", sep="")))

  file.copy(from = file.path(original_dir, paste("ASV.tsv", sep="")), to ="ASV.tsv",
            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
  file.remove(file.path(original_dir, paste("ASV.tsv", sep="")))

}

#' @title Bacteria
#'
#' @description This function create a new folder named Bacteria and set your working directory into this folder. Please, run setwd("Bacteria") after this function.
#'
#' @param nothing It's important not to write anything between the brackets, a new folder named Bacteria will be created and your working directory will be set into this folder, depending of the selected Kingdom.
#'
#' @return A new folder named Bacteria will be created and your working directory will be set into this folder, depending of the selected Kingdom.
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{Bacteria()}
#' # Please, run setwd("Bacteria") after this function.
#' @export
Bacteria <- function (nothing)
{
  if (!dir.exists("Bacteria")) {dir.create("Bacteria")}
  oldwd <- getwd()
  on.exit(setwd(oldwd))
  setwd("Bacteria")
  Bacteria= getwd()
}

#' @title Fungi
#'
#' @description This function create a new folder named Fungi and set your working directory into this folder. Please, run setwd("Fungi") after this function.
#'
#' @param nothing It's important not to write anything between the brackets, a new folder named Fungi will be created and your working directory will be set into this folder, depending of the selected Kingdom.
#'
#' @return A new folder named Fungi will be created and your working directory will be set into this folder, depending of the selected Kingdom.
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{Fungi()}
#' # plaese, run setwd("Fungi") after this function.
#' @export
Fungi <- function (nothing)
{
  if (!dir.exists("Fungi")) {dir.create("Fungi")}
  oldwd <- getwd()
  on.exit(setwd(oldwd))
  setwd("Fungi")
  Fungi= getwd()
}

#' @title Metazoan
#'
#' @description This function create a new folder named Metazoan and set your working directory into this folder. Please, run setwd("Metazoan") after this function.
#'
#' @param nothing It's important not to write anything between the brackets, a new folder named Metazoan will be created and your working directory will be set into this folder, depending of the selected Kingdom.
#'
#' @return A new folder named Metazoan will be created and your working directory will be set into this folder, depending of the selected Kingdom.
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{Metazoan()}
#' # plaese, run setwd("Metazoan") after this function.
#' @export
Metazoan <- function (nothing)
{
  if (!dir.exists("Metazoan")) {dir.create("Metazoan")}
  oldwd <- getwd()
  on.exit(setwd(oldwd))
  setwd("Metazoan")
  Metazoan= getwd()
}


#' @title format_input
#'
#' @description Apply logP on both positive and negative ASVs FC
#'
#' @param x Object from the Differential ASV abundance (DASVA) analysis
#'
#' @return an input for the input_global_analysis() function
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{format_input(x)}
#' @export
format_input <- function(x)
{
  dat <- data.frame(x)
  dat$ASV_ID <- rownames(dat)

  dat_x <- merge(taxo,dat, by.x="Feature.ID", by.y="ASV_ID", all.x=F)

  dat_x[is.na(dat_x)] <- 0
  dat <- dat_x

  dat_negativ <- dat[ which(dat$log2FoldChange < 0), ]
  dat_negativ_02 <- data.frame(dat_negativ$Feature.ID, log10(dat_negativ$padj))
  colnames(dat_negativ_02) <- c("ASV_ID", "logP")

  dat_positiv <- dat[ which(dat$log2FoldChange >= 0), ]
  dat_positiv_02 <- data.frame(dat_positiv$Feature.ID, -log10(dat_positiv$padj))
  colnames(dat_positiv_02) <- c("ASV_ID", "logP")

  input <- rbind(dat_negativ_02, dat_positiv_02)
  is.na(input)<-sapply(input, is.infinite)
  input[is.na(input)]<-0
  return(input)
}

#' @title heatmap_samples_hclust
#'
#' @description Adapt hclust for heatmap sample to sample analysis
#'
#' @param nothing It's important not to write anything between the brackets, all inputs will be adapted automatically.
#'
#' @return hclust object for the heatmap.2() function
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{hc <- heatmap_samples_hclust()}
#' @export
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

#' @title heatmap_samples_matrix
#'
#' @description Adapt samples matrix for heatmap sample to sample analysis
#'
#' @param nothing It's important not to write anything between the brackets, all inputs will be adapted automatically.
#'
#' @return samples matrix object for the heatmap.2() function
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{mat <- heatmap_samples_matrix()}
#' @export
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

#' @title heatmap_taxo
#'
#' @description Adding taxonomy in the pheatmap plot, instead of ASVs codes
#'
#' @param nothing It's important not to write anything between the brackets, all inputs will be adapted automatically.
#'
#' @return log2.norm.counts_taxo used fro adding taxonomy in the pheatmap plot, instead of ASVs codes
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{log2.norm.counts_taxo <- heatmap_taxo()}
#' @export
heatmap_taxo <- function (nothing)
{
  colnames(log2.norm.counts) <- row.names(heatmap_condition_df)
  dat <- merge(log2.norm.counts, taxo, by.x="row.names", by.y="Feature.ID", sort = TRUE)

  dat$ASV_code <- NULL
  dat$Confidence <- NULL
  dat$NCBITaxon <- NULL
  dat$Taxon.1 <- NULL
  dat$Row.names <- NULL
  dat$length_taxon <- 1:length(dat$Taxon)
  dat$length_taxon <- as.numeric(dat$length_taxon)

  # In order to deal with duplicated taxon name
  row.names(dat) <- c(paste(dat$Taxon, dat$length_taxon))
  dat$length_taxon <- NULL
  dat$Taxon <- NULL
  return(dat)
}

#' @title heatmap_condition
#'
#' @description For Clustering step. Fill directly the annotation_col variable of the pheatmap() function
#'
#' @param nothing It's important not to write anything between the brackets, all inputs will be adapted automatically.
#'
#' @return Fill directly the annotation_col variable of the pheatmap() function
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @export
heatmap_condition <- function (nothing)
{
  df <- as.data.frame(colData(dasva)["condition"])
  row.names(df) <- gsub('input_', '', row.names(df))
  row.names(df) <- gsub('.txt', '', row.names(df))

  df <- data.matrix(df, rownames.force=NA)

  return(df)
}

#' @title heatmap_data_dasva
#'
#' @description For Clustering step. Create the log2.norm.counts object.
#'
#' @param nothing It's important not to write anything between the brackets, all inputs will be adapted automatically.
#'
#' @return Create the log2.norm.counts object.
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @export
heatmap_data_dasva <- function (nothing)
{
  #Selects only most abundant ASVs (here 75)
  nCounts <- counts(dasva, normalized=TRUE)
  select <- order(rowMeans(nCounts),decreasing=TRUE)[1:75]

  #Selects corresponding norm counts
  nt <- normTransform(dasva)
  log2.norm.counts <- assay(nt)[select,]

  #Gets the metadata
  log2.norm.counts <- data.frame(log2.norm.counts)
  names(log2.norm.counts) <- gsub("input_", "", names(log2.norm.counts))
  names(log2.norm.counts) <- gsub(".txt", "", names(log2.norm.counts))

  colnames(log2.norm.counts) <- colData(dasva)$Id

  return(log2.norm.counts)
}

#' @title PCA_data_dasva
#'
#' @description Compute the PCA (Pincipal Component Analysis) data.
#'
#' @param nothing It's important not to write anything between the brackets, all inputs will be adapted automatically.
#'
#' @return data. The PCA (Pincipal Component Analysis) data.
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @export
PCA_data_dasva <- function (nothing)
{

  sampleSheet="SampleSheet_comparison.txt"

  samplesInfo <- read.table(sampleSheet, header=F) ; head(samplesInfo) ; dim(samplesInfo)
  target_file <- read.table(sampleSheet, header=F) ; head(samplesInfo) ; dim(samplesInfo)
  colnames(target_file) <- c("label","files","condition") ; head(target_file) ; dim(target_file)

  #threshold=1 # minimum of ASVs value across samples

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
  data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)

  return(data)
}

#' @title plotSparsityASV
#'
#' @description Create a plot of Sparsity ASV
#'
#' @param x Corresponding to the DASVA (Differential ASV abundance) object
#' @param normalized normalized
#' @param ... ...
#' @return A plot of Sparsity ASV
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{plotSparsityASV(dasva)}
#' @export
plotSparsityASV <- function(x, normalized=TRUE, ...) {
  if (is(x, "DESeqDataSet")) {
    x <- counts(x, normalized=normalized)
  }
  rs <- rowSums(x)
  rmx <- apply(x, 1, max)
  plot(rs[rs > 0], (rmx/rs)[rs > 0], log="x", ylim=c(0,1), xlab="sum of counts per ASV",
       ylab="max count / sum", main="Concentration of counts over total sum of counts", ...)
}

#' @title plotDispASVs
#'
#' @description Create a plot of Dispersion ASV
#'
#' @param object Corresponding to the DASVA (Differential ASV abundance) object
#'
#' @param ymin ymin
#' @param CV CV
#' @param genecol genecol
#' @param fitcol fitcol
#' @param finalcol finalcol
#' @param legend legend
#' @param xlab xlab
#' @param ylab ylab
#' @param log log
#' @param cex cex
#' @param ... ...
#' @return A plot of Dispersion ASV
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{plotDispASVs(dasva)}
#' @export
plotDispASVs <- function( object, ymin, CV=FALSE,
                          genecol = "black", fitcol = "red", finalcol = "dodgerblue",
                          legend=TRUE, xlab, ylab, log = "xy", cex = 0.45, ... )
{
  if (missing(xlab)) xlab <- "mean of normalized counts"
  if (missing(ylab)) {
    if (CV) {
      ylab <- "coefficient of variation"
    } else {
      ylab <- "dispersion"
    }
  }

  px = mcols(object)$baseMean
  sel = (px>0)
  px = px[sel]

  # transformation of dispersion into CV or not
  f <- if (CV) sqrt else I

  py = f(mcols(object)$dispGeneEst[sel])
  if(missing(ymin))
    ymin = 10^floor(log10(min(py[py>0], na.rm=TRUE))-0.1)

  plot(px, pmax(py, ymin), xlab=xlab, ylab=ylab,
       log=log, pch=ifelse(py<ymin, 6, 20), col=genecol, cex=cex, ... )

  # use a circle over outliers
  pchOutlier <- ifelse(mcols(object)$dispOutlier[sel],1,16)
  cexOutlier <- ifelse(mcols(object)$dispOutlier[sel],2*cex,cex)
  lwdOutlier <- ifelse(mcols(object)$dispOutlier[sel],2,1)
  if (!is.null(dispersions(object))) {
    points(px, f(dispersions(object)[sel]), col=finalcol, cex=cexOutlier,
           pch=pchOutlier, lwd=lwdOutlier)
  }

  if (!is.null(mcols(object)$dispFit)) {
    points(px, f(mcols(object)$dispFit[sel]), col=fitcol, cex=cex, pch=16)
  }

  if (legend) {
    legend("bottomright",c("ASVs","fitted","final"),pch=16,
           col=c(genecol,fitcol,finalcol),bg="white")
  }
}


#' @title get_dasva
#'
#' @description Creates the DASVA object. Fit a Gamma-Poisson Generalized Linear Model, dispersion estimates for Negative Binomial distributed data, "parametric", local" or "mean"
#'
#' @param fitType Fit a Gamma-Poisson Generalized Linear Model, dispersion estimates for Negative Binomial distributed data, "parametric", local" or "mean"
#'
#' @return DASVA object
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{dasva <- get_dasva(fitType="parametric")
#' dasva <- get_dasva(fitType="local")
#' dasva <- get_dasva(fitType="mean")}
#' @export
get_dasva <- function (fitType="")
{

  sampleSheet="SampleSheet_comparison.txt"

  samplesInfo <- read.table(sampleSheet, header=F) ; head(samplesInfo) ; dim(samplesInfo)
  target_file <- read.table(sampleSheet, header=F) ; head(samplesInfo) ; dim(samplesInfo)
  colnames(target_file) <- c("label","files","condition") ; head(target_file) ; dim(target_file)

  #threshold=1 # minimum of ASVs value across samples

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

  #dasva object
  dasva <- suppressMessages(DESeq(dasva_raw, fitType=fitType))
  #dasva <- suppressMessages(DESeq(dasva_raw, fitType="parametric"))
  #options(warn=-1)
  return(dasva)
}

#' @title dasva_raw_input
#'
#' @description Used in heatmap_samples_hclust(), heatmap_samples_matrix(), PCA_data_dasva() and get_dasva() functions.
#'
#' @param sampleTable Depending of the heatmap_samples_hclust(), heatmap_samples_matrix(), PCA_data_dasva() and get_dasva() functions.
#'
#' @param directory directory
#' @param design design
#' @param ignoreRank ignoreRank
#' @param ... ...
#' @return object
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#'   \dontrun{dasva_raw <- dasva_raw_input(sampleTable = sampleTable,
#'   directory = targeted_analysis_dir,
#'   design= ~ condition)}
#' @export
dasva_raw_input <- function (sampleTable, directory = ".", design, ignoreRank = FALSE,
                             ...)
{
  if (missing(design))
    stop("design is missing")
  l <- lapply(as.character(sampleTable[, 2]), function(fn) read.table(file.path(directory,
                                                                                fn), fill = TRUE))
  if (!all(sapply(l, function(a) all(a$V1 == l[[1]]$V1))))
    stop("ASV IDs (first column) differ between files.")
  tbl <- sapply(l, function(a) a[, ncol(a)])
  colnames(tbl) <- sampleTable[, 1]
  rownames(tbl) <- l[[1]]$V1
  rownames(sampleTable) <- sampleTable[, 1]
  oldSpecialNames <- c("no_feature", "ambiguous", "too_low_aQual",
                       "not_aligned", "alignment_not_unique")
  specialRows <- (substr(rownames(tbl), 1, 1) == "_") | rownames(tbl) %in%
    oldSpecialNames
  tbl <- tbl[!specialRows, , drop = FALSE]
  object <- DESeqDataSetFromMatrix(countData = tbl, colData = sampleTable[,
                                                                          -(1:2), drop = FALSE], design = design, ignoreRank, ...)
  return(object)
}


#' @title get_input_files
#'
#' @description Created sub directory "Targeted_analysis" if not already exist. Then, create one file by condition into it, and then upload the taxonomy file. Please, run setwd("01_Targeted_analysis") after this function.
#'
#' @param nothing It's important not to write anything between the brackets, all inputs will be adapted automatically.
#'
#' @return taxo
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{taxo <- get_input_files()}
#' # please, run setwd("01_Targeted_analysis") after this function.
#' @export
get_input_files <- function (nothing)
{

  taxo <- read.table("taxonomy.tsv", sep="\t", header=T)
  #taxo$Taxon <- gsub(";", "__", taxo$Taxon)
  ASV <- read.table("ASV.tsv", sep="\t", header=T)
  ##dat_samples <- merge(ASV, taxo, by.x="ASV_ID", by.y="Feature.ID")
  ##dat_samples$ASV_ID <- NULL
  ##dat_samples$Confidence <- NULL

  # I did this extrat step in order to deal with the Error in estimateSizeFactorsForMatrix --> every ASVs can contains at least one zero and so cannot compute log geometric means
  # see https://help.galaxyproject.org/t/error-with-deseq2-every-gene-contains-at-least-one-zero/564
  ASV_02 <- ASV
  ASV_02$ASV_ID <- NULL
  ASV_02 <- ASV_02 + 1
  ASV_02$ASV_ID <- ASV$ASV_ID


  # https://stackoverflow.com/questions/5620885/how-does-one-reorder-columns-in-a-data-frame
  ##arrange df vars by position
  ## vars must be a named vector, e.g. c("var.name""=1)
  arrange.vars <- function(data, vars){
    ##stop if not a data.frame (but should work for matrices as well)
    stopifnot(is.data.frame(data))

    ##sort out inputs
    data.nms <- names(data)
    var.nr <- length(data.nms)
    var.nms <- names(vars)
    var.pos <- vars
    ##sanity checks
    stopifnot( !any(duplicated(var.nms)),
               !any(duplicated(var.pos)) )
    stopifnot( is.character(var.nms),
               is.numeric(var.pos) )
    stopifnot( all(var.nms %in% data.nms) )
    stopifnot( all(var.pos > 0),
               all(var.pos <= var.nr) )

    ##prepare output
    out.vec <- character(var.nr)
    out.vec[var.pos] <- var.nms
    out.vec[-var.pos] <- data.nms[ !(data.nms %in% var.nms) ]
    stopifnot( length(out.vec)==var.nr )

    ##re-arrange vars by position
    data <- data[ , out.vec]
    return(data)
  }

  ASV_02 <- arrange.vars(ASV_02, c("ASV_ID"=1))


  dat <- ASV_02

  oldwd <- getwd()
  on.exit(setwd(oldwd))

  # Created sub directory "Targeted_analysis" if not already exist
  if (!dir.exists("01_Targeted_analysis")) {dir.create("01_Targeted_analysis")}
  setwd("01_Targeted_analysis")

  file.copy(from = file.path(original_dir, paste("SampleSheet_comparison.txt", sep="")), to ="SampleSheet_comparison.txt",
            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
  file.remove(file.path(original_dir, paste("SampleSheet_comparison.txt", sep="")))

  targeted_analysis_dir= getwd()     #path to current directory containing sample sheet and count files from step 04

  # Create one file by condition
  # https://stackoverflow.com/questions/42532940/how-to-create-multiple-text-files-for-each-column-in-a-dataframe-and-keep-the-fi
  for (i in names(dat)) {
    raw_file <- paste("input_", i, ".txt", sep = "")
    write.table(cbind(dat[1], dat[i]), raw_file, sep = "\t", quote = F, row.names = F, col.names = F, append = F)
  }

  system("rm input_ASV_ID.txt")
  unlink("input_ASV_ID.txt")

  sampleSheet="SampleSheet_comparison.txt"

  samplesInfo <- read.table(sampleSheet, header=F) ; head(samplesInfo) ; dim(samplesInfo)
  target_file <- read.table(sampleSheet, header=F) ; head(samplesInfo) ; dim(samplesInfo)
  colnames(target_file) <- c("label","files","condition") ; head(target_file) ; dim(target_file)

  return(taxo)

}


#' @title plotPCA.san
#'
#' @description Custom plotPCA function to plot PC1 et PC3
#'
#' @param object An object use for the PCA
#' @param intgroup intgroup
#' @param ntop ntop
#' @param returnData returnData
#' @return A PCA
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{plotPCA.san(object)}
#' @export
plotPCA.san <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE)
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])

  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }

  d <- data.frame(PC1 = pca$x[, 1], PC3 = pca$x[, 3], group = group,
                  intgroup.df, name = colData(rld)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:3]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC3", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3)

}

#' @title plotMA.dasva
#'
#' @description Custom MA plots for the Differential ASV abundance (DASVA) analysis. defining a new function to plot all ASVs and not only log2FoldChange > 2
#'
#' @param object Object from the Differential ASV abundance (DASVA) analysis
#' @param ... ...
#' @param alpha alpha
#' @param main main
#' @param xlab xlab
#' @param ylim ylim
#' @param MLE MLE
#' @param verbose verbose
#' @return A MA plot
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{plotMA.dasva(rXXX, main="XXX", ylim=c(-20,20))}
#' @export
plotMA.dasva <- function(object, alpha, main="", xlab="mean of normalized counts", ylim, MLE=FALSE, verbose=TRUE, ...) {
  if (missing(alpha)) {
    alpha <- if (is.null(metadata(object)$alpha)) {
      0.1
    } else {
      metadata(object)$alpha
    }
  }
  df <- if (MLE) {
    # test if MLE is there
    if (is.null(object$lfcMLE)) {
      stop("lfcMLE column is not present: you should first run results() with addMLE=TRUE")
    }
    data.frame(mean = object$baseMean,
               lfc = object$lfcMLE,
               isDE = ifelse(is.na(object$padj), FALSE, object$padj < alpha))
  } else {
    data.frame(mean = object$baseMean,
               lfc = object$log2FoldChange,
               isDE = ifelse(is.na(object$padj), FALSE, object$padj < alpha))
    #isDE = ifelse(object$padj < alpha & abs(object$log2FoldChange) >= 2, TRUE, FALSE))

  }
  if(verbose){print(dim(df[df$isDE==TRUE,]))}
  if (missing(ylim)) {
    plotMA(df, main=main, xlab=xlab, ...)
  } else {
    plotMA(df, main=main, xlab=xlab, ylim=ylim, ...)
  }
}


#' @title clusteringGOs
#'
#' @description clusteringGOs from DESeq2 analysis pipeline
#'
#' @param gen2go from DESeq2 analysis pipeline
#'
#' @param div div
#' @param cutHeight cutHeight
#' @return a clustering GO
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{clusteringGOs()}
#' @export
clusteringGOs=function(gen2go,div,cutHeight) {
  inname=paste("dissim0_",div,"_",gen2go,sep="")
  diss=read.table(inname,sep="\t",header=T,check.names=F)
  row.names(diss)=names(diss)
  hc=hclust(as.dist(diss),method="complete")
  cc=cutree(hc,h=cutHeight)
  outname=paste("cl_",inname,sep="")
  write.csv(cc,file=outname,quote=F)
}

#' @title mwuTest
#'
#' @description Mann-Whitney U Test from RBGOA
#'
#' @param gotable from gomwuStats from RBGOA
#' @param Alternative from gomwuStats from RBGOA
#'
#' @return mwuTest
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{mwuTest()}
#' @export
mwuTest=function(gotable,Alternative) {
  gos=gotable
  terms=levels(gos$term)
  gos$seq=as.character(gos$seq)
  nrg=gos[!duplicated(gos$seq),5]
  names(nrg)=gos[!duplicated(gos$seq),4]
  #	nrg=nrg+rnorm(nrg,sd=0.01) # to be able to do exact wilcox test
  rnk=rank(nrg)
  names(rnk)=names(nrg)
  pvals=c();drs=c();nams=c();levs=c();nseqs=c()
  for (t in terms){
    got=gos[gos$term==t,]
    got=got[!duplicated(got$seq),]
    ngot=gos[gos$term!=t,]
    ngot=ngot[!duplicated(ngot$seq),]
    ngot=ngot[!(ngot$seq %in% got$seq),]
    sgo.yes=got$seq
    n1=length(sgo.yes)
    sgo.no=ngot$seq
    n2=length(sgo.no)
    wi=wilcox.test(nrg[sgo.yes],nrg[sgo.no],alternative=Alternative)	# removed correct=FALSE
    r1=sum(rnk[sgo.yes])/n1
    r0=sum(rnk[sgo.no])/n2
    dr=r1-r0
    drs=append(drs,round(dr,0))
    levs=append(levs,got$lev[1])
    nams=append(nams,as.character(got$name[1]))
    pvals=append(pvals,wi$p.value)
    nseqs=append(nseqs,n1)
  }
  res=data.frame(cbind("delta.rank"=drs,"pval"=pvals,"level"=levs,nseqs))
  res=cbind(res,"term"=as.character(terms),"name"=nams)
  res$pval=as.numeric(as.character(res$pval))
  res$delta.rank=as.numeric(as.character(res$delta.rank))
  res$level=as.numeric(as.character(res$level))
  res$nseqs=as.numeric(as.character(res$nseqs))
  return(res)
}

#' @title fisherTest
#'
#' @description Ficher Test from RBGOA
#'
#' @param gotable from gomwuStats from RBGOA
#'
#' @return fisherTest
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{fisherTest()}
#' @export
fisherTest=function(gotable) {
  gos=gotable
  terms=levels(gos$term)
  gos$seq=as.character(gos$seq)
  pft=c();nam=c();lev=c();nseqs=c()
  for (t in terms) {
    got=gos[gos$term==t,]
    got=got[!duplicated(got$seq),]
    ngot=gos[gos$term!=t,]
    ngot=ngot[!duplicated(ngot$seq),]
    ngot=ngot[!(ngot$seq %in% got$seq),]
    go.sig=sum(got$value)
    go.ns=length(got[,1])-go.sig
    ngo.sig=sum(ngot$value)
    ngo.ns=length(ngot[,1])-ngo.sig
    sig=c(go.sig,ngo.sig) # number of significant genes belonging and not belonging to the tested GO category
    ns=c(go.ns,ngo.ns) # number of not-significant genes belonging and not belonging to the tested GO category
    mm=matrix(c(sig,ns),nrow=2,dimnames=list(ns=c("go","notgo"),sig=c("go","notgo")))
    ff=fisher.test(mm,alternative="greater")
    pft=append(pft,ff$p.value)
    nam=append(nam,as.character(got$name[1]))
    lev=append(lev,got$lev[1])
    nseqs=append(nseqs,length(got[,1]))
  }
  res=data.frame(cbind("delta.rank"=rep(0),"pval"=pft,"level"=lev,nseqs,"term"=terms,"name"=nam))
  res[,1]=as.numeric(as.character(res[,1]))
  res[,2]=as.numeric(as.character(res[,2]))
  res[,3]=as.numeric(as.character(res[,3]))
  res$nseqs=as.numeric(as.character(res$nseqs))
  return(res)
}

#' @title taxon_mwuStats
#'
#' @description mwuStats from RBGOA adapted for taxonomic analysis
#'
#' @param input input
#' @param goDatabase goDatabase
#' @param goAnnotations goAnnotations
#' @param goDivision goDivision
#' @param Alternative Alternative
#' @param adjust.multcomp adjust.multcomp
#' @param clusterCutHeight clusterCutHeight
#' @param largest largest
#' @param smallest smallest
#' @param perlPath perlPath
#' @param verbose verbose
#'
#' @return Statistical analysis for taxonomic rank
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{taxon_mwuStats(input, ...)}
#' @export
taxon_mwuStats=function(input,goDatabase,goAnnotations, goDivision, Alternative="t", adjust.multcomp="BH", clusterCutHeight=0.25,largest=0.1,smallest=5,perlPath="perl", verbose=TRUE){

  # file.copy(from = file.path(original_dir, paste("Working_scripts/mwu_a_NCBITaxon.pl", sep="")), to ="mwu_a_NCBITaxon.pl",
  #            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
  file.copy(from = system.file("extdata", "mwu_a_NCBITaxon.pl", package="Anaconda"), to ="mwu_a_NCBITaxon.pl",
            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

  #file.copy(from = file.path(original_dir, paste("Working_scripts/mwu_b_NCBITaxon.pl", sep="")), to ="mwu_b_NCBITaxon.pl",
  #          overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

  file.copy(from = system.file("extdata", "mwu_b_NCBITaxon.pl", package="Anaconda"), to ="mwu_b_NCBITaxon.pl",
            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

  extraOptions=paste("largest=",largest," smallest=",smallest,"cutHeight=",clusterCutHeight,sep="")

  system(paste(perlPath,"./mwu_a_NCBITaxon.pl",goDatabase,goAnnotations,input,goDivision,extraOptions))
  clusteringGOs(goAnnotations,goDivision,clusterCutHeight)
  system(paste(perlPath,"./mwu_b_NCBITaxon.pl",goAnnotations,input,goDivision))

  inname=paste(goDivision,"_",input,sep="")
  rsq=read.table(inname,sep="\t",header=T)
  rsq$term=as.factor(rsq$term)

  mwut.t=TRUE
  if (length(levels(as.factor(rsq$value)))==2) {
    if(verbose){cat("Binary classification detected; will perform Fisher's test\n");}
    mwut.t=F
    rr=fisherTest(rsq)
  } else {
    if(verbose){cat("Continuous measure of interest: will perform MWU test\n"); }
    rr=mwuTest(rsq,Alternative)
  }

  if (adjust.multcomp=="shuffle"){
    if(verbose){cat("shuffling values to calculate FDR, 5 reps\n")}
    reps=5
    spv=c()
    for (i in 1:reps) {
      if(verbose){print(paste("replicate",i))}
      rsqq=rsq
      rsqq$value=sample(rsq$value)
      if (mwut.t==TRUE) { rs=mwuTest(rsqq,Alternative) } else { rs=fisherTest(rsqq) }
      spv=append(spv,rs$pval)
    }
    fdr=c()
    for (p in rr$pval){
      fdr=append(fdr,(sum(spv<=p)/reps)/sum(rr$pval<=p))
    }
    fdr[fdr>1]=1
  } else {
    fdr=p.adjust(rr$pval,method=adjust.multcomp)
  }

  if(verbose){cat(paste(sum(fdr<0.1)," NCBITaxon terms at 10% FDR\n"))}
  rr$p.adj=fdr
  fname=paste("MWU_",inname,sep="")
  write.table(rr,fname,row.names=F)
}



#' @title taxon_mwuPlot
#'
#' @description taxon mwuPlot for taxonomic analysis
#'
#' @param inFile inFile - results object from the DASVA analysis
#' @param goAnnotations parallel to goAnnotations from gomwuStats from RBGOA. Here, "database_bacteria_package_all.tab" if Bacteria, "database_fungi_package_all.tab" if Fungi
#' @param goDivision parallel to goAnnotations from gomwuStats from RBGOA. Here, "TR" = taxonomic Rank, don't change this
#' @param level1 level1
#' @param level2 level2
#' @param level3 level3
#' @param absValue absValue
#' @param adjusted adjusted
#' @param txtsize txtsize
#' @param font.family font.family
#' @param treeHeight treeHeight
#' @param colors colors
#' @param verbose verbose
#'
#' @return taxon mwuPlot and goods "Table_02_taxon_mwuPlot.txt"
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{taxon_mwuPlot(input,...)}
#' @export
taxon_mwuPlot=function(inFile,goAnnotations,goDivision,level1=0.1,level2=0.05,level3=0.01,absValue=-log(0.05,10),adjusted=TRUE,txtsize=1,font.family="sans",treeHeight=0.5,colors=NULL, verbose=TRUE) {

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  input=inFile
  in.mwu=paste("MWU",goDivision,input,sep="_")
  in.dissim=paste("dissim",goDivision,goAnnotations,sep="_")

  cutoff=-log(level1,10)
  pv=read.table(in.mwu,header=T)
  row.names(pv)=pv$term
  in.raw=paste(goDivision,input,sep="_")
  rsq=read.table(in.raw,sep="\t",header=T)
  rsq$term=as.factor(rsq$term)

  if (adjusted==TRUE) { pvals=pv$p.adj } else { pvals=pv$pval }
  heat=data.frame(cbind("pval"=pvals))
  row.names(heat)=pv$term
  heat$pval=-log(heat$pval+1e-15,10)
  heat$direction=0
  heat$direction[pv$delta.rank>0]=1
  if (cutoff>0) {
    goods=subset(heat,pval>=cutoff)
  } else {
    goods.names=unique(rsq$term[abs(rsq$value)>=absValue])
    goods=heat[row.names(heat) %in% goods.names,]
  }

  if (is.null(colors) | length(colors)<4 ) {
    colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
    if (sum(goods$direction)==nrow(goods) | sum(goods$direction)==0) {
      colors=c("black","black","grey50","grey50")
    }
  }
  goods.names=row.names(goods)

  # reading and subsetting dissimilarity matrix
  diss=read.table(in.dissim,sep="\t",header=T,check.names=F)
  row.names(diss)=names(diss)
  diss.goods=diss[goods.names,goods.names]

  # how many genes out of what we started with we account for with our best categories?
  good.len=c();good.genes=c()
  for (g in goods.names) {
    sel=rsq[rsq$term==g,]
    pass=abs(sel$value)>=absValue
    sel=sel[pass,]
    good.genes=append(good.genes,as.character(sel$seq))
    good.len=append(good.len,nrow(sel))
  }
  ngenes=length(unique(good.genes))

  ################### HERE TO DELETE GENES NUMBERS

  #hist(rsq$value)
  totSum=length(unique(rsq$seq[abs(rsq$value)>=absValue]))
  # row.names(goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
  row.names(goods)=paste(pv[pv$term %in% goods.names,]$name,sep="") # modifier #############################################################################################
  #row.names(heat)=paste(good.len,"/",pv$nseqs," ",pv$name,sep="")
  row.names(heat)=paste(pv$name,sep="") # modifier ####################################
  #  row.names(diss.goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
  row.names(diss.goods)=paste(pv[pv$term %in% goods.names,]$name,sep="")

  # clustering terms better than cutoff
  GO.categories=as.dist(diss.goods)
  cl.goods=hclust(GO.categories,method="average")
  labs=cl.goods$labels[cl.goods$order] # saving the labels to order the plot
  goods=goods[labs,]
  labs=sub(" activity","",labs)

  old.par <- par( no.readonly = TRUE )

  plots=layout(matrix(c(1,2,3),1,3,byrow=T),c(treeHeight,3,1),TRUE)

  par(mar = c(2,2,0.85,0))
  plot(as.phylo(cl.goods),show.tip.label=FALSE,cex=0.0000001)
  step=100
  left=1
  top=step*(2+length(labs))

  par(mar = c(0,0,0.3,0))
  plot(c(1:top)~c(1:top),type="n",axes=F,xlab="",ylab="")
  ii=1
  goods$color=1
  goods$color[goods$direction==1 & goods$pval>cutoff]=colors[4]
  goods$color[goods$direction==0 & goods$pval>cutoff]=colors[3]
  goods$color[goods$direction==1 & goods$pval>(-log(level2,10))]=colors[2]
  goods$color[goods$direction==0 & goods$pval>(-log(level2,10))]=colors[1]
  goods$color[goods$direction==1 & goods$pval>(-log(level3,10))]=colors[2]
  goods$color[goods$direction==0 & goods$pval>(-log(level3,10))]=colors[1]
  for (i in length(labs):1) {
    ypos=top-step*ii
    ii=ii+1
    if (goods$pval[i]> -log(level3,10)) {
      text(left,ypos,labs[i],font=2,cex=1*txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
    } else {
      if (goods$pval[i]>-log(level2,10)) {
        text(left,ypos,labs[i],font=1,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
      } else {
        #			if (goods$pval[i]>cutoff) {
        #				text(left,ypos,labs[i],font=3,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
        #		} else {
        text(left,ypos,labs[i],font=3,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
        #}
      }
    }
  }

  par(mar = c(3,1,1,0))

  plot(c(1:top)~c(1:top),type="n",axes=F,xlab="",ylab="")
  text(left,top-step*2,paste("p < ",level3,sep=""),font=2,cex=1* txtsize,adj=c(0,0),family=font.family)
  text(left,top-step*3,paste("p < ",level2,sep=""),font=1,cex=0.8* txtsize,adj=c(0,0),family=font.family)
  text(left,top-step*4,paste("p < ",10^(-cutoff),sep=""),font=3,col="grey50",cex=0.8* txtsize,adj=c(0,0),family=font.family)

  if(verbose){cat(paste("NCBITaxon terms dispayed: ",length(goods.names)),"\n")}
  if(verbose){cat(paste("\"Good ASVs\" accounted for:  ", ngenes," out of ",totSum, " ( ",round(100*ngenes/totSum,0), "% )","\n",sep=""))}
  par(old.par)
  goods$pval=10^(-1*goods$pval)

  write.table(goods, "Table_02_taxon_mwuPlot.txt")


  return(goods)
}

#' @title taxon_mwu_list
#'
#' @description taxon Mann-Whitney U list for taxonomic analysis
#'
#' @param inFile inFile - results object from the DASVA analysis
#' @param goAnnotations parallel to goAnnotations from gomwuStats from RBGOA. Here, "database_bacteria_package_all.tab" if Bacteria, "database_fungi_package_all.tab" if Fungi
#' @param goDivision parallel to goAnnotations from gomwuStats from RBGOA. Here, "TR" = taxonomic Rank, don't change this
#' @param level1 level1
#' @param level2 level2
#' @param level3 level3
#' @param absValue absValue
#' @param adjusted adjusted
#' @param txtsize txtsize
#' @param font.family font.family
#' @param treeHeight treeHeight
#' @param colors colors
#'
#' @return List for the statistical analysis for taxonomic rank
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{taxon_list <- taxon_mwu_list(input, ...)}
#' @export
taxon_mwu_list=function(inFile,goAnnotations,goDivision,level1=0.1,level2=0.05,level3=0.01,absValue=-log(0.05,10),adjusted=TRUE,txtsize=1,font.family="sans",treeHeight=0.5,colors=NULL) {

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  input=inFile
  in.mwu=paste("MWU",goDivision,input,sep="_")
  in.dissim=paste("dissim",goDivision,goAnnotations,sep="_")

  cutoff=-log(level1,10)
  pv=read.table(in.mwu,header=T)
  row.names(pv)=pv$term
  in.raw=paste(goDivision,input,sep="_")
  rsq=read.table(in.raw,sep="\t",header=T)
  rsq$term=as.factor(rsq$term)

  if (adjusted==TRUE) { pvals=pv$p.adj } else { pvals=pv$pval }
  heat=data.frame(cbind("pval"=pvals))
  row.names(heat)=pv$term
  heat$pval=-log(heat$pval+1e-15,10)
  heat$direction=0
  heat$direction[pv$delta.rank>0]=1
  if (cutoff>0) {
    goods=subset(heat,pval>=cutoff)
  } else {
    goods.names=unique(rsq$term[abs(rsq$value)>=absValue])
    goods=heat[row.names(heat) %in% goods.names,]
  }

  if (is.null(colors) | length(colors)<4 ) {
    colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
    if (sum(goods$direction)==nrow(goods) | sum(goods$direction)==0) {
      colors=c("black","black","grey50","grey50")
    }
  }
  goods.names=row.names(goods)

  # reading and subsetting dissimilarity matrix
  diss=read.table(in.dissim,sep="\t",header=T,check.names=F)
  row.names(diss)=names(diss)
  diss.goods=diss[goods.names,goods.names]

  # how many genes out of what we started with we account for with our best categories?
  good.len=c();good.genes=c()
  for (g in goods.names) {
    sel=rsq[rsq$term==g,]
    pass=abs(sel$value)>=absValue
    sel=sel[pass,]
    good.genes=append(good.genes,as.character(sel$seq))
    good.len=append(good.len,nrow(sel))
  }
  ngenes=length(unique(good.genes))

  ################### HERE TO DELETE GENES NUMBERS

  #hist(rsq$value)
  totSum=length(unique(rsq$seq[abs(rsq$value)>=absValue]))
  # row.names(goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
  row.names(goods)=paste(pv[pv$term %in% goods.names,]$name,sep="") # modifier #############################################################################################
  #row.names(heat)=paste(good.len,"/",pv$nseqs," ",pv$name,sep="")
  row.names(heat)=paste(pv$name,sep="") # modifier ####################################
  #  row.names(diss.goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
  row.names(diss.goods)=paste(pv[pv$term %in% goods.names,]$name,sep="")

  # clustering terms better than cutoff
  GO.categories=as.dist(diss.goods)
  cl.goods=hclust(GO.categories,method="average")
  labs=cl.goods$labels[cl.goods$order] # saving the labels to order the plot
  goods=goods[labs,]
  labs=sub(" activity","",labs)

  old.par <- par( no.readonly = TRUE )


  # cat(paste("NCBITaxon terms dispayed: ",length(goods.names)),"\n")
  # cat(paste("\"Good ASVs\" accounted for:  ", ngenes," out of ",totSum, " ( ",round(100*ngenes/totSum,0), "% )","\n",sep=""))
  # par(old.par)
  # goods$pval=10^(-1*goods$pval)
  #return(goods)

  write.table(data.frame(labs), "Table_03_taxon_mwu_list.txt")


  return(data.frame(labs))

}


#' @title taxon_mwuPlot_guilds
#'
#' @description taxon Mann-Whitney U Plot with Fungi Guilds added
#'
#' @param inFile inFile - results object from the DASVA analysis
#' @param goAnnotations parallel to goAnnotations from gomwuStats from RBGOA. Here, "database_bacteria_package_all.tab" if Bacteria, "database_fungi_package_all.tab" if Fungi
#' @param goDivision parallel to goAnnotations from gomwuStats from RBGOA. Here, "TR" = taxonomic Rank, don't change this
#' @param level1 level1
#' @param level2 level2
#' @param level3 level3
#' @param absValue absValue
#' @param adjusted adjusted
#' @param txtsize txtsize
#' @param font.family font.family
#' @param treeHeight treeHeight
#' @param colors colors
#' @param verbose verbose
#'
#' @return List for the statistical analysis for taxonomic rank
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{taxon_mwuPlot_guilds(input, ...)}
#' @export
taxon_mwuPlot_guilds=function(inFile,goAnnotations,goDivision,level1=0.1,level2=0.05,level3=0.01,absValue=-log(0.05,10),adjusted=TRUE,txtsize=1,font.family="sans",treeHeight=0.5,colors=NULL, verbose=TRUE) {

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  input=inFile
  in.mwu=paste("MWU",goDivision,input,sep="_")
  in.dissim=paste("dissim",goDivision,goAnnotations,sep="_")

  cutoff=-log(level1,10)
  pv=read.table(in.mwu,header=T)
  row.names(pv)=pv$term
  in.raw=paste(goDivision,input,sep="_")
  rsq=read.table(in.raw,sep="\t",header=T)
  rsq$term=as.factor(rsq$term)

  if (adjusted==TRUE) { pvals=pv$p.adj } else { pvals=pv$pval }
  heat=data.frame(cbind("pval"=pvals))
  row.names(heat)=pv$term
  heat$pval=-log(heat$pval+1e-15,10)
  heat$direction=0
  heat$direction[pv$delta.rank>0]=1
  if (cutoff>0) {
    goods=subset(heat,pval>=cutoff)
  } else {
    goods.names=unique(rsq$term[abs(rsq$value)>=absValue])
    goods=heat[row.names(heat) %in% goods.names,]
  }

  if (is.null(colors) | length(colors)<4 ) {
    colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
    if (sum(goods$direction)==nrow(goods) | sum(goods$direction)==0) {
      colors=c("black","black","grey50","grey50")
    }
  }
  goods.names=row.names(goods)

  # reading and subsetting dissimilarity matrix
  diss=read.table(in.dissim,sep="\t",header=T,check.names=F)
  row.names(diss)=names(diss)
  diss.goods=diss[goods.names,goods.names]

  # how many genes out of what we started with we account for with our best categories?
  good.len=c();good.genes=c()
  for (g in goods.names) {
    sel=rsq[rsq$term==g,]
    pass=abs(sel$value)>=absValue
    sel=sel[pass,]
    good.genes=append(good.genes,as.character(sel$seq))
    good.len=append(good.len,nrow(sel))
  }
  ngenes=length(unique(good.genes))

  ################### HERE TO DELETE GENES NUMBERS

  #hist(rsq$value)
  totSum=length(unique(rsq$seq[abs(rsq$value)>=absValue]))
  # row.names(goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
  row.names(goods)=paste(pv[pv$term %in% goods.names,]$name,sep="") # modifier #############################################################################################
  #row.names(heat)=paste(good.len,"/",pv$nseqs," ",pv$name,sep="")
  row.names(heat)=paste(pv$name,sep="") # modifier ####################################
  #  row.names(diss.goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
  row.names(diss.goods)=paste(pv[pv$term %in% goods.names,]$name,sep="")

  # clustering terms better than cutoff
  GO.categories=as.dist(diss.goods)
  cl.goods=hclust(GO.categories,method="average")
  labs=cl.goods$labels[cl.goods$order] # saving the labels to order the plot
  goods=goods[labs,]
  labs=sub(" activity","",labs)

  old.par <- par( no.readonly = TRUE )

  #plots=layout(matrix(c(1,1,1,2,
  #                      3,4,5,6,
  #                      3,4,5,6,
  #                      3,4,5,6),c(4,3.5,
  #                                        treeHeight,3,1,3), byrow=T))

  plots=layout(matrix(c(1,1,1,1,2,2,2,2,2,2,
                        3,3,3,3,4,4,4,4,4,4,
                        3,3,3,3,4,4,4,4,4,4,
                        3,3,3,3,4,4,4,4,4,4,
                        3,3,3,3,4,4,4,4,4,4,
                        5,6,7,7,8,8,8,8,8,8,
                        5,6,7,7,8,8,8,8,8,8,
                        5,6,7,7,8,8,8,8,8,8,
                        5,6,7,7,8,8,8,8,8,8,
                        5,6,7,7,8,8,8,8,8,8,
                        5,6,7,7,8,8,8,8,8,8,
                        5,6,7,7,8,8,8,8,8,8,
                        5,6,7,7,8,8,8,8,8,8,
                        5,6,7,7,8,8,8,8,8,8), nrow=14, byrow=T))

  #plots=layout(matrix(c(1,2,3),1,3,byrow=T),c(treeHeight,3,1),TRUE)
  #####  plots=layout(matrix(c(1, 2, 3, 3,
  #####                        4, 5, 6, 6,
  #####                        4, 5, 6, 6), nrow=3, byrow=TRUE))
  #####
  #layout.show(n=6)

  step=100
  left=1
  top=step*(2+length(labs))

  #par(mar = c(1,1,1,1))
  #plot(c(1:top)~c(1:top),type="n",axes=F,xlab="",ylab="")
  #text(left,top-step*2,paste("p < ",level3,sep=""),font=2,cex=1* txtsize,adj=c(0,0),family=font.family)
  #text(left,top-step*3,paste("p < ",level2,sep=""),font=1,cex=0.8* txtsize,adj=c(0,0),family=font.family)
  #text(left,top-step*4,paste("p < ",10^(-cutoff),sep=""),font=3,col="grey50",cex=0.8* txtsize,adj=c(0,0),family=font.family)

  par(mar = c(0,0,0.3,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.39, y = 0.5, paste("Taxa enrichment"),
       cex = 2.5, col = "black", family="serif", font=2, adj=0.5)


  par(mar = c(0,0,0.3,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.39, y = 0.5, paste("Fungal guilds"),
       cex = 2.5, col = "black", family="serif", font=2, adj=0.5)


  par(mar = c(0,0,0.3,0))
  plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="", xlab="", axes=FALSE)
  legend(0, 1, legend=c("Over-represented ASVs in the first condition",
                        "Under-represented ASVs in the first condition"),
         col=c("firebrick1",
               "dodgerblue2"), pch=19, cex=1, box.lty=0, y.intersp=0.55)

  par(mar = c(0,0,0.3,0))
  plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="", xlab="", axes=FALSE)
  legend(0,1, legend=c("Pathotroph-Saprotroph-Symbiotroph",
                       "Pathotroph-Saprotroph",
                       "Saprotroph-Symbiotroph",
                       "Pathotroph-Symbiotroph",
                       "Saprotroph",
                       "Symbiotroph",
                       "Pathotroph",
                       "Unassigned",
                       "NA"),
         col=c("black",
               "steelblue2",
               "plum",
               "hotpink",
               "forestgreen",
               "goldenrod3",
               "brown4",
               "cornsilk4",
               "azure4"), pch=19, cex=1, box.lty=0, y.intersp=0.55)

  par(mar = c(3,1,1,0))

  plot(c(1:top)~c(1:top),type="n",axes=F,xlab="",ylab="")
  text(left,top-step*2,paste("p < ",level3,sep=""),font=2,cex=1* txtsize,adj=c(0,0),family=font.family)
  text(left,top-step*3,paste("p < ",level2,sep=""),font=1,cex=0.8* txtsize,adj=c(0,0),family=font.family)
  text(left,top-step*4,paste("p < ",10^(-cutoff),sep=""),font=3,col="grey50",cex=0.8* txtsize,adj=c(0,0),family=font.family)


  par(mar = c(2,2,0.85,0))
  plot(as.phylo(cl.goods),show.tip.label=F,cex=0.0000001)
  #plot(as.dendrogram(cl.goods))
  step=100
  left=1
  top=step*(2+length(labs))

  par(mar = c(0,0,0.3,0))
  plot(c(1:top)~c(1:top),type="n",axes=F,xlab="",ylab="")
  ii=1
  goods$color=1
  goods$color[goods$direction==1 & goods$pval>cutoff]=colors[4]
  goods$color[goods$direction==0 & goods$pval>cutoff]=colors[3]
  goods$color[goods$direction==1 & goods$pval>(-log(level2,10))]=colors[2]
  goods$color[goods$direction==0 & goods$pval>(-log(level2,10))]=colors[1]
  goods$color[goods$direction==1 & goods$pval>(-log(level3,10))]=colors[2]
  goods$color[goods$direction==0 & goods$pval>(-log(level3,10))]=colors[1]
  for (i in length(labs):1) {
    ypos=top-step*ii
    ii=ii+1
    if (goods$pval[i]> -log(level3,10)) {
      text(left,ypos,labs[i],font=2,cex=1*txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
    } else {
      if (goods$pval[i]>-log(level2,10)) {
        text(left,ypos,labs[i],font=1,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
      } else {
        #			if (goods$pval[i]>cutoff) {
        #				text(left,ypos,labs[i],font=3,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
        #		} else {
        text(left,ypos,labs[i],font=3,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
        #}
      }
    }
  }

  ###labs_02 <- merge(data.frame(labs), link_guilds, by.x="labs", by.y="V1")
  # From qed answer in https://stackoverflow.com/questions/17878048/merge-two-data-frames-while-keeping-the-original-row-order
  merge_sameord = function(x, y, ...) {
    UseMethod('merge_sameord')
  }

  merge_sameord.data.frame = function(x, y, ...) {
    rstr = paste(sample(c(0:9, letters, LETTERS), 12, replace=TRUE), collapse='')
    x[, rstr] = 1:nrow(x)
    res = merge(x, y, all.x=TRUE, sort=FALSE, ...)
    res = res[order(res[, rstr]), ]
    res[, rstr] = NULL
    res
  }




  labs_02 <- merge_sameord(data.frame(labs), link_guilds, by.x="labs", by.y="V1")
  # labs_02 <- data.table::fintersect(setDT(labs$labs), setDT(link_guilds))
  labs_02$V3 <- sapply(labs_02$V3,function(x) {x <- gsub("Pathotroph-Saprotroph-Symbiotroph","black",x)})
  labs_02$V3 <- sapply(labs_02$V3,function(x) {x <- gsub("Pathotroph-Saprotroph","steelblue2",x)})
  labs_02$V3 <- sapply(labs_02$V3,function(x) {x <- gsub("Saprotroph-Symbiotroph","plum",x)})
  labs_02$V3 <- sapply(labs_02$V3,function(x) {x <- gsub("Pathotroph-Symbiotroph","hotpink",x)})
  labs_02$V3 <- sapply(labs_02$V3,function(x) {x <- gsub("Saprotroph","forestgreen",x)})
  labs_02$V3 <- sapply(labs_02$V3,function(x) {x <- gsub("Symbiotroph","goldenrod3",x)})
  labs_02$V3 <- sapply(labs_02$V3,function(x) {x <- gsub("Pathotroph","brown4",x)})
  labs_02$V3 <- sapply(labs_02$V3,function(x) {x <- gsub("Unassigned","cornsilk4",x)})
  labs_02$V3 <- sapply(labs_02$V3,function(x) {x <- gsub("NA","azure4",x)})

  ### ICI CHANGER SAPRO EN VERT, ETC...

  #par(mar = c(2,2,0.85,0))
  par(mar = c(0,0,0.3,0))
  #plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="y label", xlab="x lablel")
  plot(c(1:top)~c(1:top),type="n",axes=F,xlab="",ylab="")
  ii=1
  goods$color=1
  goods$color[goods$direction==1 & goods$pval>cutoff]=colors[4]
  goods$color[goods$direction==0 & goods$pval>cutoff]=colors[3]
  goods$color[goods$direction==1 & goods$pval>(-log(level2,10))]=colors[2]
  goods$color[goods$direction==0 & goods$pval>(-log(level2,10))]=colors[1]
  goods$color[goods$direction==1 & goods$pval>(-log(level3,10))]=colors[2]
  goods$color[goods$direction==0 & goods$pval>(-log(level3,10))]=colors[1]
  for (i in length(labs):1) {
    ypos=top-step*ii
    ii=ii+1
    if (goods$pval[i]> -log(level3,10)) {
      # text(left,ypos,paste(labs[i], labs_02$V2[i]),font=2,cex=1*txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
      #text(left,ypos,expression(labs[i]*phantom(labs_02$V2[i])),col=goods$color[i],font=2,cex=1*txtsize,adj=c(0,0),family=font.family)
      #text(left,ypos,expression(phantom(labs[i])*labs_02$V2[i]),col="black",font=2,cex=1*txtsize,adj=c(0,0),family=font.family)

      ### FOR PUTTING NCBI TAXON and GUILDS Together
      ### text(left,ypos,bquote(bold(.(labs_02$labs[i]))*'  Guild  '~phantom(.( labs_02$V2[i]))), col=goods$color[i],cex=1*txtsize,adj=c(0,0),family=font.family)
      ### text(left,ypos,bquote(phantom(.(labs_02$labs[i]))*'  Guild  '~.( labs_02$V2[i])), col="black",cex=1*txtsize,adj=c(0,0),family=font.family)

      text(left,ypos,paste(labs_02$V2[i]),font=2,cex=1*txtsize,col=labs_02$V3[i],adj=c(0,0),family=font.family)

      #  text(left,ypos,labs_02$V2[i],font=2,cex=1*txtsize,col=goods$color[i],adj=c(0,0),family=font.family) + text(left,ypos,labs[i],font=2,cex=1*txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
    } else {
      if (goods$pval[i]>-log(level2,10)) {
        # text(left,ypos,labs_02$V2[i],font=1,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
        ### FOR PUTTING NCBI TAXON and GUILDS Together
        ###   text(left,ypos,bquote(.(labs_02$labs[i]) ~ "    " phantom(.( labs_02$V2[i]))),col=goods$color[i],font=1,cex=0.8*txtsize,adj=c(0,0),family=font.family)
        ###   text(left,ypos,bquote(phantom(.(labs_02$labs[i])) ~ "    " .( labs_02$V2[i])),col="black",font=1,cex=0.8*txtsize,adj=c(0,0),family=font.family)

        text(left,ypos,labs_02$V2[i],font=2,cex=1* txtsize,col=labs_02$V3[i],adj=c(0,0),family=font.family)

      } else {
        #			if (goods$pval[i]>cutoff) {
        #				text(left,ypos,labs[i],font=3,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
        #		} else {
        #  text(left,ypos,labs_02$V2[i],font=3,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
        ### FOR PUTTING NCBI TAXON and GUILDS Together
        ###   text(left,ypos,bquote(italic(.(labs_02$labs[i]) ~ "    " ~ phantom(.( labs_02$V2[i])))),col=goods$color[i],cex=0.8*txtsize,adj=c(0,0),family=font.family)
        ###   text(left,ypos,bquote(phantom(.(labs_02$labs[i])) ~ "    " ~ .( labs_02$V2[i])),col="black",cex=0.8*txtsize,adj=c(0,0),family=font.family)

        text(left,ypos,labs_02$V2[i],font=2,cex=1* txtsize,col=labs_02$V3[i],adj=c(0,0),family=font.family)

        #}
      }
    }
  }


if(verbose){cat(paste("NCBITaxon terms dispayed: ",length(goods.names)),"\n")}
if(verbose){cat(paste("\"Good ASVs\" accounted for:  ", ngenes," out of ",totSum, " ( ",round(100*ngenes/totSum,0), "% )","\n",sep=""))}
  par(old.par)
  goods$pval=10^(-1*goods$pval)

  write.table(labs_02, "Table_04_taxon_mwuPlot_guilds.txt")


  return(labs_02)

}

#' @title get_taxon_list_drawer
#'
#' @description get taxonomic list drawer
#'
#' @param taxon_list object from taxon_mwu_list() function
#'
#' @return taxon_list_drawer Object and "taxon_list_drawer_input.txt" file
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{taxon_list_drawer <- get_taxon_list_drawer(taxon_list)}
#' @export
get_taxon_list_drawer <- function(taxon_list) {

  taxon_list$labs<-gsub(" ", "_", taxon_list$labs)
  names(taxon_list) <- NULL

  write.table(taxon_list, "taxon_list.txt", row.names=FALSE, quote = FALSE)

  # file.copy(from = file.path(original_dir, paste("Working_scripts/search_taxonomic_drawer.sh", sep="")), to ="search_taxonomic_drawer.sh",
  #            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

  file.copy(from = system.file("extdata", "search_taxonomic_drawer.sh", package="Anaconda"), to ="search_taxonomic_drawer.sh",
            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

  # file.copy(from = file.path(kingdom, paste("taxonomy.tsv", sep="")), to ="taxonomy.tsv",
  #            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

  file.copy(from = file.path(kingdom, paste("taxonomy.tsv", sep="")), to ="taxonomy.tsv",
            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

  system('bash search_taxonomic_drawer.sh')
  # If Taxon have NAs, that's corresponding to "holes" in the taxonomic.tsv file used in QIIME2 (limited to current knownledge and current worldwide database).

  taxon_list_drawer <- read.table("taxon_list_drawer.txt")
  taxon_list_drawer <- data.frame(rep(1:length(taxon_list_drawer$V1), 1), taxon_list_drawer)
  #colnames(taxon_list_drawer) <- c("ID", "taxonomy")
  names(taxon_list_drawer) <- NULL
  write.table(taxon_list_drawer, "taxon_list_drawer_input.txt", sep="\t", row.names=FALSE, quote = TRUE)
  system("sed -i '' $'1s/^/ID\ttaxonomy\\\n/' taxon_list_drawer_input.txt")

  taxon_list_drawer <- read.table("taxon_list_drawer_input.txt", sep="\t", header=T)

  return(taxon_list_drawer)
}

#' @title get_funguilds
#'
#' @description get Fungi Guilds from taxon_list_drawer Object
#'
#' @param taxon_list_drawer object from get_taxon_list_drawer() function
#'
#' @return funguilds Object
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{funguilds <- get_funguilds(taxon_list_drawer)}
#' @export
get_funguilds <- function(taxon_list_drawer) {

  #file.copy(from = file.path(original_dir, paste("Working_scripts/Guilds_v1.1.py", sep="")), to ="Guilds_v1.1.py",
  #          overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

  file.copy(from = system.file("extdata", "Guilds_v1.1.py", package="Anaconda"), to ="Guilds_v1.1.py",
            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

  taxon_list_drawer <- read.table("taxon_list_drawer_input.txt", sep="\t", header=T)

  # to deal with https://stackoverflow.com/questions/17309288/importerror-no-module-named-requests
  # if
  #FunGuild v1.1 Beta
  #Traceback (most recent call last):
  #  File "Guilds_v1.1.py", line 119, in <module>
  #  import requests
  #ImportError: No module named requests
  system('python -m pip install --user requests')

  # https://stackoverflow.com/questions/41638558/how-to-call-python-script-from-r-with-arguments
  system('python Guilds_v1.1.py -otu taxon_list_drawer_input.txt -m -u')

  ## CLEAN taxon_list_drawer_input.guilds.txt
  dat <- read.table("taxon_list_drawer_input.guilds_matched.txt", sep="\t", header=T)

  # https://stackoverflow.com/questions/26289681/r-regex-find-last-occurrence-of-delimiter
  funguilds <- data.frame(sub(".*__", "", dat$taxonomy,),  dat$Guild, dat$Trophic.Mode,stringsAsFactors = FALSE)
  # funguilds <- data.frame(sub(".*[__]", "", dat$taxonomy,),  dat$Guild, dat$Trophic.Mode,stringsAsFactors = FALSE)
  colnames(funguilds) <- c("Taxon",  "Guild", "Trophic_Mode")

  unlink("Guilds_v1.1.py")

  return(funguilds)
}


#' @title get_link_guilds
#'
#' @description get link guilds from taxon_list and funguilds Objects
#'
#' @param taxon_list object from taxon_mwu_list() function
#' @param funguilds object from get_funguilds() function
#'
#' @return link_guilds Object
#' @import ggrepel pheatmap lookup plyr ape DESeq2 ggplot2 stats utils data.table RColorBrewer rafalib
#' @importFrom graphics layout legend par points text
#' @importFrom methods is
#' @examples
#' \dontrun{link_guilds <- get_link_guilds(taxon_list, funguilds)}
#' @export
get_link_guilds <- function(taxon_list, funguilds) {
  # MATCH taxon_list with dat_02
  funguilds$Taxon_02 <- data.frame(sub("_", " ", funguilds$Taxon))
  funguilds_02 <- data.frame(funguilds, funguilds$Taxon_02)
  colnames(funguilds_02) <- c("Taxon", "Guild", "Trophic_Mode", "Taxon_02", "Taxon_03")
  funguilds_03 <- data.frame(funguilds_02$Taxon_03, funguilds_02$Guild, funguilds_02$Trophic_Mode)
  colnames(funguilds_03) <- c("Taxon", "Guild", "Trophic_Mode")
  link <- merge(taxon_list, funguilds_03, by.x="labs", by.y="Taxon", all=TRUE)
  link$Guild[is.na(link$Guild)] <- "Unassigned"
  link$Trophic_Mode[is.na(link$Trophic_Mode)] <- "Unassigned"
  colnames(link) <- c("V1", "V2", "V3")
  link_guilds <- link
  return(link_guilds)
}


