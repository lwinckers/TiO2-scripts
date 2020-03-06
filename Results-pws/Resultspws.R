###
###
###

# set up environment
rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

options(stringsAsFactors = FALSE)

# 
pws <- read.table("C:/Users/Laurent/Desktop/Test/data-output/ann_selected_pws.txt", header = T, sep = "\t")
gsea <- read.table("C:/Users/Laurent/Desktop/Test/data-output/GSEA/testGSEA_GSEA_result1.txt", header = T, sep = "\t")
gsea2 <- read.table("C:/Users/Laurent/Desktop/Test/data-output/GSEA/testGSEA_GSEA_result2.txt", header = T, sep = "\t")
gsea3 <- read.table("C:/Users/Laurent/Desktop/Test/data-output/GSEA/testGSEA_GSEA_result3.txt", header = T, sep = "\t")
gsea4 <- read.table("C:/Users/Laurent/Desktop/Test/data-output/GSEA/testGSEA_GSEA_result4.txt", header = T, sep = "\t")
gsea5 <- read.table("C:/Users/Laurent/Desktop/Test/data-output/GSEA/testGSEA_GSEA_result5.txt", header = T, sep = "\t")
gsea6 <- read.table("C:/Users/Laurent/Desktop/Test/data-output/GSEA/testGSEA_GSEA_result6.txt", header = T, sep = "\t")

#
test <- gsea[gsea$ID %in% pws$ID[pws$GOterm == "GO0006915"],]
test2 <- gsea[gsea$ID %in% pws$ID[pws$GOterm == "GO0006954"],]
test3 <- gsea[gsea$ID %in% pws$ID[pws$GOterm == "GO0006974"],]
test4 <- gsea[gsea$ID %in% pws$ID[pws$GOterm == "GO0034599"],]

sig <- gsea[gsea$pvalue < 0.01,]
sig2 <- gsea2[gsea2$pvalue < 0.01,]
sig3 <- gsea3[gsea3$pvalue < 0.01,]
sig4 <- gsea4[gsea4$pvalue < 0.01,]
sig5 <- gsea5[gsea5$pvalue < 0.01,]
sig6 <- gsea6[gsea6$pvalue < 0.01,]

write.table(sig$ID, "C:/Users/Laurent/Desktop/sig1.txt", row.names = F, sep = "\t", quote = FALSE, col.names = F)
write.table(sig2$ID, "C:/Users/Laurent/Desktop/sig2.txt", row.names = F, sep = "\t", quote = FALSE, col.names = F)
write.table(sig3$ID, "C:/Users/Laurent/Desktop/sig3.txt", row.names = F, sep = "\t", quote = FALSE, col.names = F)
write.table(sig4$ID, "C:/Users/Laurent/Desktop/sig4.txt", row.names = F, sep = "\t", quote = FALSE, col.names = F)
write.table(sig5$ID, "C:/Users/Laurent/Desktop/sig5.txt", row.names = F, sep = "\t", quote = FALSE, col.names = F)
write.table(sig6$ID, "C:/Users/Laurent/Desktop/sig6.txt", row.names = F, sep = "\t", quote = FALSE, col.names = F)

#####################################

setwd("C:/Users/Laurent/Desktop/Test")

library(qusage)
library(plyr)

# load in pathway database gmt files
path = paste0(getwd(), "/data-input")
gmtFile <- list.files(path = path, pattern = ".gmt")
for (i in 1:length(gmtFile)){
  assign(gmtFile[i], read.gmt(paste0(path,"/",gmtFile[i])))
}
gmtFile <- mget(ls(pattern = ".gmt"))
list2env(lapply(gmtFile, function(x){ldply(x,data.frame)}), envir = .GlobalEnv)

databases <- do.call(rbind, lapply(ls(pattern = ".gmt"), get))
colnames(databases) <- c("pathway", "entrezgene")
databases[] <- lapply(databases, gsub, pattern='[[:punct:]]', replacement=' ')


# edge table
edge_table1 <- as.data.frame(databases[databases$pathway %in% sig$ID,])
edge_table2 <- as.data.frame(databases[databases$pathway %in% sig2$ID,])
edge_table3 <- as.data.frame(databases[databases$pathway %in% sig3$ID,])
edge_table4 <- as.data.frame(databases[databases$pathway %in% sig4$ID,])
edge_table5 <- as.data.frame(databases[databases$pathway %in% sig5$ID,])
edge_table6 <- as.data.frame(databases[databases$pathway %in% sig6$ID,])

write.table(edge_table1, "C:/Users/Laurent/Desktop/edge1.txt", row.names = F, sep = "\t", quote = FALSE, col.names = T)
write.table(edge_table2, "C:/Users/Laurent/Desktop/edge2.txt", row.names = F, sep = "\t", quote = FALSE, col.names = T)
write.table(edge_table3, "C:/Users/Laurent/Desktop/edge3.txt", row.names = F, sep = "\t", quote = FALSE, col.names = T)
write.table(edge_table4, "C:/Users/Laurent/Desktop/edge4.txt", row.names = F, sep = "\t", quote = FALSE, col.names = T)
write.table(edge_table5, "C:/Users/Laurent/Desktop/edge5.txt", row.names = F, sep = "\t", quote = FALSE, col.names = T)
write.table(edge_table6, "C:/Users/Laurent/Desktop/edge6.txt", row.names = F, sep = "\t", quote = FALSE, col.names = T)

##############

sig$type <- "Pathway"
sig2$type <- "Pathway"
sig3$type <- "Pathway"
sig4$type <- "Pathway"
sig5$type <- "Pathway"
sig6$type <- "Pathway"

write.table(sig[c(1,12)], "C:/Users/Laurent/Desktop/node1.txt", row.names = F, sep = "\t", quote = FALSE, col.names = T)
write.table(sig2[c(1,12)], "C:/Users/Laurent/Desktop/node2.txt", row.names = F, sep = "\t", quote = FALSE, col.names = T)
write.table(sig3[c(1,12)], "C:/Users/Laurent/Desktop/node3.txt", row.names = F, sep = "\t", quote = FALSE, col.names = T)
write.table(sig4[c(1,12)], "C:/Users/Laurent/Desktop/node4.txt", row.names = F, sep = "\t", quote = FALSE, col.names = T)
write.table(sig5[c(1,12)], "C:/Users/Laurent/Desktop/node5.txt", row.names = F, sep = "\t", quote = FALSE, col.names = T)
write.table(sig6[c(1,12)], "C:/Users/Laurent/Desktop/node6.txt", row.names = F, sep = "\t", quote = FALSE, col.names = T)
