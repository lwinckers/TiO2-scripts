### Combining genes from two GO terms. Annotation of these genes.
### Laurent Winckers, Maastricht University - BiGCaT
### 14-11-2019

### BiomaRt 2.40.5 (through GO-annotation function)

# clean work space
rm(list=ls())
options(stringsAsFactors = F)

# set working directroy
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# link functions 
source("./functions/GO_annotation.R")

##### GENE FILE PER GO TERM #####
path = paste0(getwd(), "/data-input")
ls_GO <- list.files(path = path, pattern = "GO00")
for (i in 1:length(ls_GO)){
  assign(ls_GO[i], read.table(paste0(path,"/",ls_GO[i]), header = T, sep = "\t")[c("GO.TERM", "SYMBOL")])
  }
ls_GO <- mget(ls(pattern = "GO00"))

# map entrez gene IDs for genes
for (i in 1:length(ls_GO)){
GO_annotation(fileName = paste0("ann_", names(ls_GO[i])), data = ls_GO[[i]], values = ls_GO[[i]]$SYMBOL,
           attributes = c('hgnc_symbol', 'entrezgene_id'),
           filter = "hgnc_symbol",
           path = getwd())
}

##### COMBINED GO TERM LIST #####

# load in annotated files
path = paste0(getwd(), "/data-output")
ls_GO <- list.files(path = path, pattern = "ann_GO00")
for (i in 1:length(ls_GO)){assign(ls_GO[i], read.table(paste0(path,"/",ls_GO[i]), header = T, sep = "\t")[c("GO.TERM", "entrezgene_id")])}

# combine all genes from GO terms
GO <- unique(rbind(ls(pattern = "ann_GO00")))
GO <- unique(rbind(ann_GO0006915.txt.txt,ann_GO0006954.txt.txt,ann_GO0006974.txt.txt,ann_GO0034599.txt.txt))
colnames(GO) <- c("hgnc_symbol", "entrezgene")

# save common gene list and mapped gene list
write.table(GO ,"./data-output/ann_genes.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# information about session
sessionInfo()
