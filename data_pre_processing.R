## ---------------------------
##
## Script name: data_pre_processing.R
##
## Purpose of script: pre-process data that will be used in: https://github.com/laurent2207/TiO2-scripts.
## Gene-expression data: GEO:GSE42069
## GO-term genelists: Apoptopic process (GO:0006915), Inflammatory response (GO:0006954)
## GO-term genelists: Cellular response to DNA damage stimulus (GO:0006974), Cellular response to oxidative stress (GO:0034599)
## Pathway databases: WikiPathways-20190910-gmt-Homo_sapiens, KEGG-v7.0-entrezgene, Reactome-v7.0-entrez
##
## Author: Laurent Winckers, Martina Kutmon
##
## Date Created: 2020-03-23
##
## Session info:
## R version 3.6.1 (2019-07-05)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 18362)
## Packages: biomaRt_2.42.0
##
## ---------------------------

##### Set up environment #####

rm(list=ls())
options(stringsAsFactors = F)

### set working directroy
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

### install packages 
library(biomaRt)

### link functions
source("./functions/GO_annotation.R")

##### Gene-expression data pre-processing #####

### Load in gene-expression data, previously analysed using ArrayAnalysis (arrayanalysis.org)
# Caco2
caco2_low <- read.table("./data-input/TiO2_24hrs_10_Caco2.txt", sep = "\t", header = T)
caco2_high <- read.table("./data-input/TiO2_24hrs_100_Caco2.txt", sep = "\t", header = T)
# SAE
SAE_low <- read.table("./data-input/TiO2_24hrs_10_SAE.txt", sep = "\t", header = T)
SAE_high <- read.table("./data-input/TiO2_24hrs_100_SAE.txt", sep = "\t", header = T)
# THP1
THP1_low <- read.table("./data-input/TiO2_24hrs_10_THP1.txt", sep = "\t", header = T)
THP1_high <- read.table("./data-input/TiO2_24hrs_100_THP1.txt", sep = "\t", header = T)

### select necessary columns, remove columns that are not needed further down the process
# Caco2
caco2_low <- caco2_low[c(1,2,3,6)]
caco2_high <- caco2_high[c(1,2,3,6)]
# SAE
SAE_low <- SAE_low[c(1,2,3,6)]
SAE_high <- SAE_high[c(1,2,3,6)]
# THP1
THP1_low <- THP1_low[c(1,2,3,6)]
THP1_high <- THP1_high[c(1,2,3,6)]

### change column names adressing respective cell line and either 10 ug/ml, low (L) or 100 ug/ml, high (H) concentration
# Caco2
colnames(caco2_low)[c(2,3,4)] <- c("caco2_L_logFC", "caco2_L_FC", "caco2_L_pval")
colnames(caco2_high)[c(2,3,4)] <- c("caco2_H_logFC", "caco2_H_FC", "caco2_H_pval")
# SAE
colnames(SAE_low)[c(2,3,4)] <- c("SAE_L_logFC", "SAE_L_FC", "SAE_L_pval")
colnames(SAE_high)[c(2,3,4)] <- c("SAE_H_logFC", "SAE_H_FC", "SAE_H_pval")
# THP1
colnames(THP1_low)[c(2,3,4)] <- c("THP1_L_logFC", "THP1_L_FC", "THP1_L_pval")
colnames(THP1_high)[c(2,3,4)] <- c("THP1_H_logFC", "THP1_H_FC", "THP1_H_pval")

### select ensembl IDs from one of the datasets
ids <- as.data.frame(caco2_low$ENSG_ID)
colnames(ids) <- "ENSG_ID"

### map entrezgene IDs to hgnc symbols
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

genes <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id'), 
  filters = 'ensembl_gene_id',
  values = ids$ENSG_ID,
  mart = ensembl
)

### remove NAs from gene file
## 11324
genes2 <- genes[!(is.na(genes$entrezgene_id)),]

### remove rows without HGNC-symbol
## 11298
genes3 <- genes2[-which(genes2$hgnc_symbol == ""),]

### merge datasets with annotated gene identifiers
colnames(genes3)[1] <- "ENSG_ID"
bound <- merge(genes3, caco2_low, by = "ENSG_ID")
bound <- merge(bound, caco2_high, by = "ENSG_ID")
bound <- merge(bound, SAE_low, by = "ENSG_ID")
bound <- merge(bound, SAE_high, by = "ENSG_ID")
bound <- merge(bound, THP1_low, by = "ENSG_ID")
bound <- merge(bound, THP1_high, by = "ENSG_ID")

### create seperate files for each cell line
caco2 <- bound[,c(1:9)]
SAE <- bound[,c(1:3,10:15)]
THP1 <- bound[,c(1:3,16:21)]

### save annotated gene-expression files for each cell line and combined
write.table(bound, "./data-output/merged_TiO2.txt", sep = "\t", quote = F, row.names = F)
write.table(caco2, "./data-output/merged_caco2.txt", sep = "\t", quote = F, row.names = F)
write.table(SAE, "./data-output/merged_SAE.txt", sep = "\t", quote = F, row.names = F)
write.table(THP1, ".//data-output/merged_THP1.txt", sep = "\t", quote = F, row.names = F)


##### GO-term genelists pre-processing #####

rm(list=ls())

### load in GO-term genelists, only GO.TERM and SYMBOL column 
path = paste0(getwd(), "/data-input")

ls_GO <- list.files(path = path, pattern = "GO00")

for (i in 1:length(ls_GO)){
  assign(ls_GO[i], read.table(paste0(path,"/", ls_GO[i]), header = T, sep = "\t")[c("GO.TERM", "SYMBOL")])
}

### combine them in one list to loop over them
ls_GO <- mget(ls(pattern = "GO00"))

### remove .txt at end of dataframe variable names
ls_GO <- setNames(ls_GO, gsub(pattern = ".txt", replacement = "", x = names(ls_GO)))

### map HGNC-symbols to Entrezgene IDs using the GO_annotation function
for (i in 1:length(ls_GO)){
  GO_annotation(fileName = paste0("ann_", names(ls_GO[i])), 
                data = ls_GO[[i]], 
                values = ls_GO[[i]]$SYMBOL,
                attributes = c('hgnc_symbol', 'entrezgene_id'),
                filter = "hgnc_symbol",
                path = paste0(getwd(), "/data-output/"))
}


##### Combining pathway databases #####

rm(list=ls())

### load gmt files and transpose the obtained lists to data frames
path = paste0(getwd(), "/data-input")
gmtFile <- list.files(path = path, pattern = ".gmt")
for (i in 1:length(gmtFile)){
  assign(gmtFile[i], read.gmt(paste0(path,"/",gmtFile[i])))
}

gmtFile <- mget(ls(pattern = ".gmt"))
gmtFile <- lapply(gmtFile, function(x){plyr::ldply(x,data.frame)})

### combine pathway database files and clean it
databases <- do.call(rbind, gmtFile)
colnames(databases) <- c("pathway", "entrezgene")
databases$entrezgene <- as.numeric(databases$entrezgene)

### save file of combined pathway databases
write.table(databases, "./data-output/pw_databases.txt", quote = F, row.names = F, sep = "\t")
