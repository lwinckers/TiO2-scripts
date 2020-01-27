###
###
###

# set up environment
rm(list=ls())

# load libraries
library(dplyr)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# load in GSEA result files
for (i in 1:length(list.files(path="./data-output/GSEA",pattern="GSEA"))){ 
  assign(list.files(path="./data-output/GSEA",pattern="GSEA")[i], 
         read.table(paste0(getwd(), "/data-output/GSEA/", 
                           list.files(path="./data-output/GSEA",pattern="GSEA")[i]), 
                    header = T, sep = "\t"))
}

# clean and merge GSEA result files, select rows only if they contain significant pw in one of the columns
GO0006915 <- mget(ls(pattern = 'GO0006915'))
GO0006915 <- lapply(GO0006915, function(x){x[c(1,5,6)]}) 
for (i in 1:length(GO0006915)){
  colnames(GO0006915[[i]]) <- c("ID", paste0("NES_",i), paste0("pvalue_",i))
}
list2env(GO0006915, envir = .GlobalEnv)
GO0006915 <- Reduce(function(x,y) merge(x = x, y = y, by = "ID"), GO0006915)
write.table(GO0006915, "./data-output/GSEA/grouped_result_GO0006915.txt", quote = F , sep = "\t", row.names = F)
GO0006915 <- GO0006915[apply(GO0006915[grep("pvalue", names(GO0006915), value = T)], 1, function(x) any(x < 0.05)),]
write.table(GO0006915, "./data-output/GSEA/sig_result_GO0006915.txt", quote = F , sep = "\t", row.names = F)

GO0006954 <- mget(ls(pattern = 'GO0006954'))
GO0006954 <- lapply(GO0006954, function(x){x[c(1,5,6)]}) 
for (i in 1:length(GO0006954)){
  colnames(GO0006954[[i]]) <- c("ID", paste0("NES_",i), paste0("pvalue_",i))
}
list2env(GO0006954, envir = .GlobalEnv)
GO0006954 <- Reduce(function(x,y) merge(x = x, y = y, by = "ID"), GO0006954)
write.table(GO0006954, "./data-output/GSEA/grouped_result_GO0006954.txt", quote = F , sep = "\t", row.names = F)
GO0006954 <- GO0006954[apply(GO0006954[grep("pvalue", names(GO0006954), value = T)], 1, function(x) any(x < 0.05)),]
write.table(GO0006954, "./data-output/GSEA/sig_result_GO0006954.txt", quote = F , sep = "\t", row.names = F)

GO0006974 <- mget(ls(pattern = 'GO0006974'))
GO0006974 <- lapply(GO0006974, function(x){x[c(1,5,6)]}) 
for (i in 1:length(GO0006974)){
  colnames(GO0006974[[i]]) <- c("ID", paste0("NES_",i), paste0("pvalue_",i))
}
list2env(GO0006974, envir = .GlobalEnv)
GO0006974 <- Reduce(function(x,y) merge(x = x, y = y, by = "ID"), GO0006974)
write.table(GO0006974, "./data-output/GSEA/grouped_result_GO0006974.txt", quote = F , sep = "\t", row.names = F)
GO0006974 <- GO0006974[apply(GO0006974[grep("pvalue", names(GO0006974), value = T)], 1, function(x) any(x < 0.05)),]
write.table(GO0006974, "./data-output/GSEA/sig_result_GO0006974.txt", quote = F , sep = "\t", row.names = F)

GO0034599 <- mget(ls(pattern = 'GO0034599'))
GO0034599 <- lapply(GO0034599, function(x){x[c(1,5,6)]}) 
for (i in 1:length(GO0034599)){
  colnames(GO0034599[[i]]) <- c("ID", paste0("NES_",i), paste0("pvalue_",i))
}
list2env(GO0034599, envir = .GlobalEnv)
GO0034599 <- Reduce(function(x,y) merge(x = x, y = y, by = "ID"), GO0034599)
write.table(GO0034599, "./data-output/GSEA/grouped_result_GO0034599.txt", quote = F , sep = "\t", row.names = F)
GO0034599 <- GO0034599[apply(GO0034599[grep("pvalue", names(GO0034599), value = T)], 1, function(x) any(x < 0.05)),]
write.table(GO0034599, "./data-output/GSEA/sig_result_GO0034599.txt", quote = F , sep = "\t", row.names = F)

# significant pathways
sigPWs <- rbind(GO0006915, 
              GO0006954, 
              GO0006974, 
              GO0034599)
sigPWs <- as.data.frame(unique(sigPWs$ID))
write.table(sigPWs, "./data-output/GSEA/sig_pws.txt", quote = F , sep = "\t", row.names = F)
