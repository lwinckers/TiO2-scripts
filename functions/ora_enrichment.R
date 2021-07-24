### Overrepresentation analysis
### Laurent Winckers, Maastricht University - Department of Bioinformatics BiGCaT
### 2021-06-28

ora_enrichment <- function(genes, wpid2gene, wpid2name, prefix, path) {
  res <- clusterProfiler::enricher(
    genes,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 10, maxGSSize = 300,
    TERM2GENE = wpid2gene,
    TERM2NAME = wpid2name)
  
  write.table(res, paste0(path, prefix, ".txt"), quote = F, sep = "\t", row.names = F)
}