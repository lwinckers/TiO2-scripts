### Pathway selection based on common GO genes
### Laurent Winckers, Maastricht University - BiGCaT
### 14-11-2019

# clean work space
rm(list=ls())

# set working directory to where script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# install packages 
source("./functions/autoInstallPackages.R")
source("./functions/pathwaySelection.R")
using("limma", "qusage", "plyr", "dplyr", "tidyr")

# load in genes file
genes <- read.table("./data-output/ann_GO0034599.txt", sep = "\t", header = T)

# load in pathway database gmt files
path = paste0(getwd(), "/data-input")
gmtFile <- list.files(path = path, pattern = ".gmt")
for (i in 1:length(gmtFile)){assign(gmtFile[i], read.gmt(paste0(path,"/",gmtFile[i])))}

gmtFile <- mget(ls(pattern = ".gmt"))
list2env(lapply(gmtFile, function(x){ldply(x,data.frame)}), envir = .GlobalEnv)

databases <- do.call(rbind, lapply(ls(pattern = ".gmt"), get))
colnames(databases) <- c("pathway", "entrezgene")

# use pathwaySelection function
pathwaySelection(pathwayDb = databases, genes = genes$entrezgene_id, geneInfo = TRUE, 
                 number = 10, selected = 1, perc = 30, path = getwd(),
                 fileName = "GO0034599")

# information about the session.
sessionInfo()
