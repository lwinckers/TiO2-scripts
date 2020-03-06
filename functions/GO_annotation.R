### Function for GO term annotation
### Laurent Winckers
### Maastricht University - Department of Bioinformatics, BiGCaT

GO_annotation <- function(fileName,data,values,attributes,filter,path){

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  library(biomaRt)
    
  if (!dir.exists(file.path(path,"data-output"))){
    dir.create(file.path(path,"data-output"))}
  
  ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
  
  genes <<- getBM(
    attributes = attributes, 
    filters = filter,
    values = values,
    mart = ensembl
  )
  
  merged <- merge(data, genes, by.x = "SYMBOL", by.y = "hgnc_symbol")
  
  write.table(merged, paste0(path,"/data-output/",fileName,".txt"), sep = "\t", col.names = T,
              row.names = F, quote = F)
}