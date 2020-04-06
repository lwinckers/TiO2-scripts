### GSEA analysis function 
### Laurent Winckers, Maastricht University - Department of Bioinformatics BiGCaT
### 2020-03-23

library(enrichplot)

GSEAanalysis <- function(GENESET, fileName, data){
  for(i in 4:9) {
    dat <- as.numeric(data[,i])
    names(dat) <- as.character(data[,3])
    dat <- sort(dat, decreasing = T)
    # run GSEA
    res <- GSEA(dat, pvalueCutoff = 1,
                minGSSize = 10, maxGSSize = 500,
                pAdjustMethod = "BH", TERM2GENE = GENESET,
                nPerm = 1000)
    #p1 <<- gseaplot2(res, geneSetID = 1:3, pvalue_table = TRUE)
    loc = paste(getwd(), paste0("/output/GSEA/GSEA_",fileName,"-",i-3,".txt"),sep="")
    write.table(res, file = loc, sep="\t", quote = F,
                row.names = F)
}}
