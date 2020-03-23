### Function for GO term annotation
### Laurent Winckers
### Maastricht University - Department of Bioinformatics, BiGCaT

GO_annotation <- function(fileName,data,values,attributes,filter,path){

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  library(biomaRt)
  
  ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
  
  genes <<- getBM(
    attributes = attributes, 
    filters = filter,
    values = values,
    mart = ensembl
  )
  
  merged <- merge(data, genes, by.x = "SYMBOL", by.y = "hgnc_symbol")
  
  write.table(merged, paste0(path,fileName,".txt"), sep = "\t", col.names = T,
              row.names = F, quote = F)
}