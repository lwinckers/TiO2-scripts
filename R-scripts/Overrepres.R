###
###
###

# set up environment
rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

library(clusterProfiler)
library(qusage)
library(plyr)

# load annotated GO terms
path = paste0(getwd(), "/data-output")
ls_GO <- list.files(path = path, pattern = "^ann_GO")
for (i in 1:length(ls_GO)){assign(ls_GO[i], read.table(paste0(path,"/",ls_GO[i]), header = T, sep = "\t"))}
ls_GO <- mget(ls(pattern = "ann_GO"))

ls_GO <- lapply(ls_GO, function(x){as.numeric(x$entrezgene_id)})

# load gmt files and combine them
path = paste0(getwd(), "/data-input")
gmtFile <- list.files(path = path, pattern = ".gmt")
for (i in 1:length(gmtFile)){assign(gmtFile[i], read.gmt(paste0(path,"/",gmtFile[i])))}

#gmtFile <- mget(ls(pattern = ".gmt"))
#list2env(lapply(gmtFile, function(x){ldply(x,data.frame)}), envir = .GlobalEnv)

databases <- do.call(rbind, lapply(ls(pattern = ".gmt"), get))
colnames(databases) <- c("pathway", "entrezgene")

databases$entrezgene <- as.numeric(databases$entrezgene)

# enricher
for (i in 1:length(ls_GO)){
res <- as.data.frame(enricher(gene = ls_GO[[i]], minGSSize = 400, TERM2GENE = databases, pvalueCutoff = 0.05))
res$GOterm <- names(ls_GO[i])
assign(paste0("result_", names(ls_GO[i])), as.data.frame(res))
write.table(as.data.frame(res), paste0("./data-output/results_", names(ls_GO[i])), quote = F, sep = "\t", row.names = F)
}

# save files
results <- rbind(result_ann_GO0006915.txt, result_ann_GO0006954.txt,
                      result_ann_GO0006974.txt, result_ann_GO0034599.txt)
write.table(results, "./data-output/enricher_Results.txt", quote = F, sep = "\t", row.names = F)

selected_pws <- unique(results[1])
write.table(selected_pws, "./data-output/selected_pws.txt", quote = F, sep = "\t", row.names = F)
