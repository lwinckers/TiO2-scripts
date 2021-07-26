### ORA analysis function 
### Laurent Winckers, Maastricht University - Department of Bioinformatics BiGCaT
### 2021-06-28

enrichment <- function(wpid2gene,wpid2name, data, comparisons, fccutoff, path) {
  for (i in comparisons) {
    buffer <- data %>% dplyr::select(starts_with("entrez") | starts_with(i))
    
    # ORA
     up.genes <- buffer[buffer[,2] > fccutoff & buffer[,4] < 0.05, 1] 
     dn.genes <- buffer[buffer[,2] < -fccutoff & buffer[,4] < 0.05, 1]
     combined <- c(up.genes, dn.genes)
     
      ewp.combined <- clusterProfiler::enricher(
       combined,
       pvalueCutoff = 1,
       qvalueCutoff = 1,
       minGSSize = 10, maxGSSize = 300,
       TERM2GENE = wpid2gene,
       TERM2NAME = wpid2name)
     
     write.table(ewp.combined, file = paste(path, i, ".txt", sep=""), sep="\t", quote = F, row.names = F)
  }
}
