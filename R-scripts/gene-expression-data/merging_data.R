### Merging gene expression datasets obtained from GEO:GSE42069 
### Pre-processed and analysed with arrayanalysis.org
### Laurent Winckers, Maastricht University - BiGCaT
### 14-11-2019

##### SET UP ENVIRONMENT #####

# set working directroy
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# install packages 
source("./functions/autoInstallPackages.R")
using("biomaRt")

# Caco2
caco2_low <- read.table("./data-input/TiO2_24hrs_10_Caco2.txt", sep = "\t", header = T)
caco2_high <- read.table("./data-input/TiO2_24hrs_100_Caco2.txt", sep = "\t", header = T)
# SAE
SAE_low <- read.table("./data-input/TiO2_24hrs_10_SAE.txt", sep = "\t", header = T)
SAE_high <- read.table("./data-input/TiO2_24hrs_100_SAE.txt", sep = "\t", header = T)
# THP1
THP1_low <- read.table("./data-input/TiO2_24hrs_10_THP1.txt", sep = "\t", header = T)
THP1_high <- read.table("./data-input/TiO2_24hrs_100_THP1.txt", sep = "\t", header = T)

##### CLEAN DATA FRAMES #####

# Caco2
caco2_low <- caco2_low[c(1,2,6)]
caco2_high <- caco2_high[c(1,2,6)]
# SAE
SAE_low <- SAE_low[c(1,2,6)]
SAE_high <- SAE_high[c(1,2,6)]
# THP1
THP1_low <- THP1_low[c(1,2,6)]
THP1_high <- THP1_high[c(1,2,6)]

# Caco2
colnames(caco2_low)[c(2,3)] <- c(paste0(substring(list.files()[1], first = 1, last = 19),"-",colnames(caco2_low)[2]), paste0(substring(list.files()[1], first = 1, last = 19),"-",colnames(caco2_low)[3]))
colnames(caco2_high)[c(2,3)] <- c(paste0(substring(list.files()[4], first = 1, last = 20),"-",colnames(caco2_high)[2]), paste0(substring(list.files()[4], first = 1, last = 20),"-",colnames(caco2_high)[3]))
# SAE
colnames(SAE_low)[c(2,3)] <- c(paste0(substring(list.files()[2], first = 1, last = 17),"-",colnames(SAE_low)[2]), paste0(substring(list.files()[2], first = 1, last = 17),"-",colnames(SAE_low)[3]))
colnames(SAE_high)[c(2,3)] <- c(paste0(substring(list.files()[5], first = 1, last = 18),"-",colnames(SAE_high)[2]), paste0(substring(list.files()[5], first = 1, last = 18),"-",colnames(SAE_high)[3]))
# THP1
colnames(THP1_low)[c(2,3)] <- c(paste0(substring(list.files()[3], first = 1, last = 18),"-",colnames(THP1_low)[2]), paste0(substring(list.files()[3], first = 1, last = 18),"-",colnames(THP1_low)[3]))
colnames(THP1_high)[c(2,3)] <- c(paste0(substring(list.files()[6], first = 1, last = 19),"-",colnames(THP1_high)[2]), paste0(substring(list.files()[6], first = 1, last = 19),"-",colnames(THP1_high)[3]))

##### GENE ID ANNOTATION #####

ids <- as.data.frame(caco2_low$ENSG_ID)
colnames(ids) <- "ENSG_ID"

# map entrezgene IDs to hgnc symbols
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

genes <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id'), 
  filters = 'ensembl_gene_id',
  values = ids$ENSG_ID,
  mart = ensembl
)

# check which ENSG ids are not mapped
mismatch <- as.data.frame(caco2_low$ENSG_ID[!(caco2_low$ENSG_ID %in% genes$ensembl_gene_id)])
colnames(mismatch) <- "mismatch"
mismatch <- as.data.frame(mismatch[grepl("ENSG.*", mismatch$mismatch),])

##### CLEAN DATA FRAMES 2 AND SAVE #####

colnames(genes)[1] <- "ENSG_ID"
bound <- merge(genes, caco2_low, by = "ENSG_ID")
bound <- merge(bound, caco2_high, by = "ENSG_ID")
bound <- merge(bound, SAE_low, by = "ENSG_ID")
bound <- merge(bound, SAE_high, by = "ENSG_ID")
bound <- merge(bound, THP1_low, by = "ENSG_ID")
bound <- merge(bound, THP1_high, by = "ENSG_ID")

caco2 <- bound[,c(1:7)]
colnames(caco2)[c(4:7)] <- c("low-logFC", "low-P.Value", "high-logFC", "high-P.Value")
SAE <- bound[,c(1:3,8:11)]
colnames(SAE)[c(4:7)] <- c("low-logFC", "low-P.Value", "high-logFC", "high-P.Value")
THP1 <- bound[,c(1:3,12:15)]
colnames(THP1)[c(4:7)] <- c("low-logFC", "low-P.Value", "high-logFC", "high-P.Value")

write.table(bound, "./data-output/merged_TiO2.txt", sep = "\t", quote = F, row.names = F)
write.table(caco2, "./data-output/merged_caco2.txt", sep = "\t", quote = F, row.names = F)
write.table(SAE, "./data-output/merged_SAE.txt", sep = "\t", quote = F, row.names = F)
write.table(THP1, ".//data-output/merged_THP1.txt", sep = "\t", quote = F, row.names = F)

# information about session
sessionInfo()