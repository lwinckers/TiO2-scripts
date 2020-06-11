ora_enrichment <- function(genes, wpid2gene, wpid2name, prefix) {
  res <- clusterProfiler::enricher(
    genes,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 10, maxGSSize = 300,
    TERM2GENE = wpid2gene,
    TERM2NAME = wpid2name)
  
  #dp <- dotplot(res, showCategory= 30, color = "qvalue")
  cp <- cnetplot(res, showCategory=30, node_label="category")  
  #figure <- ggarrange(dp, cp, labels = c("A", "B"), nrow = 2)
  #ggexport(figure, filename = paste("output/",prefix,".png",sep=""), width = 2000, height = 2000)
  ggexport(cp, filename = paste("output/",prefix,".png",sep=""), width = 2000, height = 2000)
  
  
  write.table(res, paste0("output/",prefix,".txt"), quote = F, sep = "\t", row.names = F)
}