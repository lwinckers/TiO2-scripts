### Funtion pathway selection based on genelist
### Laurent Winckers  
### Maastricht University - Department of Bioinformatics, BiGCaT

pathwaySelection <- function(pathwayDb,genes,geneInfo,number,selected,perc,path,fileName){
  
  library(plyr)
  
  if (!dir.exists(file.path(path, "/data-output"))){
    dir.create(file.path(path, "/data-output"))}
  
  Db <- pathwayDb[pathwayDb[[1]] %in% names(table(pathwayDb[[1]]))[table(pathwayDb[[1]]) >= number],]
  frDb <- plyr::count(Db[[1]])
  Dbpw <- as.data.frame(Db[[1]][Db[[2]] %in% genes])
  Dbpw <- as.data.frame(Dbpw[Dbpw[[1]] %in% names(table(Dbpw[[1]]))[table(Dbpw[[1]]) >= selected],])
  frDbpw <- plyr::count(Dbpw[[1]])
  Dbpw <- as.data.frame(merge(frDbpw, frDb, by = "x"))
  
  colnames(Dbpw)[c(1,2,3)] <- c("pathway", "frequence of selected genes", "frequence of genes")
  Dbpw$percentage <- ((Dbpw[[2]]/Dbpw[[3]]) * 100)
  Dbpw <- Dbpw[Dbpw$percentage >= perc,]
  Dbpw$percentage <- round(Dbpw$percentage, 2)
  write.table(Dbpw, paste0(path,"/data-output/", fileName, "_pathwayStats.txt"), col.names = T, row.names = F, sep = "\t", quote = F)
  
  selpws <- Dbpw[c(-2,-3,-4)]
  write.table(selpws, paste0(path,"./data-output/", fileName, "_selected_pws.txt"), col.names = T, row.names = F, sep = "\t", quote = F)
  
  if (geneInfo == TRUE){
    nGenesTotal <- as.data.frame(Db[[2]][Db[[2]] %in% genes])
    selaDb <- as.data.frame(Db[Db[[1]] %in% Dbpw[[1]],])
    nSelGenesTotal <- as.data.frame(selaDb[[2]][selaDb[[2]] %in% genes])
    
    nGenesTotalF <- plyr::count(nGenesTotal)
    colnames(nGenesTotalF)[2] <- "number of pathways"
    
    nSelGenesTotalF <- plyr::count(nSelGenesTotal)
    colnames(nSelGenesTotalF)[2] <- "number of selected pathways"
    
    write.table(nGenesTotalF, paste0(path,"/data-output/", fileName, "_nGenes.txt"), col.names = T, row.names = F, sep = "\t", quote = F)
    write.table(nSelGenesTotalF, paste0(path,"/data-output/", fileName, "_nSelGenes.txt"), col.names = T, row.names = F, sep = "\t", quote = F)
  }
}