### GSEA analysis
### Laurent Winckers, Maastricht University - BiGCaT
### 14-11-2019

# clean work space
rm(list=ls())
options(stringsAsFactors = F)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# install packages 
library(clusterProfiler)
library(dplyr)

# load in necessary files
data <- read.table("./gene-expression-data/data-output/merged_TiO2.txt", header = T, sep ="\t")

GO <- read.table("./data-output/edges_GO_terms.txt", header = T, sep = "\t")
gene_network <- read.table("./data-output/edges_entrezgene.txt", header = T, sep = "\t")
genesets <- rbind(GO, gene_network)

#GO1 <- read.table("./data-output/ann_GO0006915.txt", header=T, sep="\t")
#GO2 <- read.table("./data-output/ann_GO0006954.txt", header=T, sep="\t")
#GO3 <- read.table("./data-output/ann_GO0006974.txt", header=T, sep="\t")
#GO4 <- read.table("./data-output/ann_GO0034599.txt", header=T, sep="\t")

#GO1 <- GO1[,c("GO.TERM", "entrezgene_id")]
#colnames(GO1) <- c("source", "target")
#GO1$target <- as.numeric(GO1$target)
#GO2 <- GO2[,c("GO.TERM", "entrezgene_id")]
#colnames(GO2) <- c("source", "target")
#GO3 <- GO3[,c("GO.TERM", "entrezgene_id")]
#colnames(GO3) <- c("source", "target")
#GO4 <- GO4[,c("GO.TERM", "entrezgene_id")]
#colnames(GO4) <- c("source", "target")

#pwGO1 <- read.table("./data-output/edges-GO0006915-pw.txt", header=T, sep="\t")
#pwGO2 <- read.table("./data-output/edges-GO0006954-pw.txt", header=T, sep="\t")
#pwGO3 <- read.table("./data-output/edges-GO0006974-pw.txt", header=T, sep="\t")
#pwGO4 <- read.table("./data-output/edges-GO0034599-pw.txt", header=T, sep="\t")

#pwGO1 <- pwGO1[c(1,2)]
#colnames(pwGO1) <- c("source", "target")
#pwGO2 <- pwGO2[c(1,2)]
#colnames(pwGO2) <- c("source", "target")
#pwGO3 <- pwGO3[c(1,2)]
#colnames(pwGO3) <- c("source", "target")
#pwGO4 <- pwGO4[c(1,2)]
#colnames(pwGO4) <- c("source", "target")

#geneset1 <- rbind(GO1, pwGO1)
#geneset2 <- rbind(GO2, pwGO2)
#geneset3 <- rbind(GO3, pwGO3)
#geneset4 <- rbind(GO4, pwGO4)

# prepare geneset
genesets$source <- as.character(genesets$source)
genesets$target <- as.numeric(genesets$target)
colnames(genesets) <- c("ID", "gene")

#geneset1$source <- as.character(geneset1$source)
#geneset1$target <- as.numeric(geneset1$target)
#colnames(geneset1) <- c("ID", "gene")
#geneset2$source <- as.character(geneset2$source)
#geneset2$target <- as.numeric(geneset2$target)
#colnames(geneset2) <- c("ID", "gene")
#geneset3$source <- as.character(geneset3$source)
#geneset3$target <- as.numeric(geneset3$target)
#colnames(geneset3) <- c("ID", "gene")
#geneset4$source <- as.character(geneset4$source)
#geneset4$target <- as.numeric(geneset4$target)
#colnames(geneset4) <- c("ID", "gene")

# save genesets
write.table(genesets, "./data-output/GSEA/geneset.txt", quote = F, row.names = F, sep = "\t")

#write.table(geneset1, "./data-output/GSEA/geneset_GO0006915.txt", quote = F, row.names = F, sep = "\t")
#write.table(geneset2, "./data-output/GSEA/geneset_GO0006954.txt", quote = F, row.names = F, sep = "\t")
#write.table(geneset3, "./data-output/GSEA/geneset_GO0006974.txt", quote = F, row.names = F, sep = "\t")
#write.table(geneset4, "./data-output/GSEA/geneset_GO0034599.txt", quote = F, row.names = F, sep = "\t")

# use specific rank-score calculation 
data_rnk <- data
data_rnk <- data_rnk %>% mutate(caco2_L_rnk = caco2_L_FC * -log10(caco2_L_pval))
data_rnk <- data_rnk %>% mutate(caco2_H_rnk = caco2_H_FC * -log10(caco2_H_pval))
data_rnk <- data_rnk %>% mutate(SAE_L_rnk = SAE_L_FC * -log10(SAE_L_pval))
data_rnk <- data_rnk %>% mutate(SAE_H_rnk = SAE_H_FC * -log10(SAE_H_pval))
data_rnk <- data_rnk %>% mutate(THP1_L_rnk = THP1_L_FC * -log10(THP1_L_pval))
data_rnk <- data_rnk %>% mutate(THP1_H_rnk = THP1_H_FC * -log10(THP1_H_pval))

# clean dataset
data_rnk <- data_rnk[c(1,2,3,(grep("_rnk", names(data_rnk))))]

# extract data 1
GSEAan <- function(GENESET, fileName){
for(i in 4:9) {
  dat <- as.numeric(data_rnk[,i])
  names(dat) <- as.character(data_rnk[,3])
  dat <- sort(dat, decreasing = T)
  # run GSEA
  res <- GSEA(dat, pvalueCutoff = 1,
              minGSSize = 10, maxGSSize = 500,
              pAdjustMethod = "BH", TERM2GENE = GENESET,
              nPerm = 1000)
  loc = paste(getwd(), paste0("/data-output/GSEA/GSEA_",fileName,i-3,".txt"),sep="")
  write.table(res, file = loc, sep="\t", quote = F,
              row.names = F)
}}

GSEAan(GENESET = genesets, fileName = "GSEA")

# test
dat <- as.numeric(data[,4])
names(dat) <- as.character(data[,3])
dat <- sort(dat, decreasing = T)
res <- GSEA(dat, pvalueCutoff = 1,
            minGSSize = 1, maxGSSize = 50000,
            pAdjustMethod = "fdr", TERM2GENE = genesets,
            nPerm = 1000)

# information about the session.
sessionInfo()
