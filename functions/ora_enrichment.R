ora_enrichment <- function(genes, wpid2gene, wpid2name, prefix) {
  res <- clusterProfiler::enricher(
    genes,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 10, maxGSSize = 300,
    TERM2GENE = wpid2gene,
    TERM2NAME = wpid2name)
  
#  dp <- dotplot(res, showCategory= 30, color = "qvalue")
  cp <- cnetplot(res, showCategory=30, node_label="category") 
  
  # combined png
#  figure <- ggarrange(dp, cp, labels = c("A", "B"), nrow = 2)
#  ggexport(figure, filename = paste("output/",prefix,".png",sep=""), width = 2000, height = 2000)
  
  # png
#  ggexport(cp, filename = paste("output/",prefix,".png",sep=""), width = 2000, height = 2000)
  
  #svg
  ggsave(file=paste("output/",prefix,".svg",sep=""), plot=cp, width=20, height=20, limitsize = FALSE)

  
  write.table(res, paste0("output/",prefix,".txt"), quote = F, sep = "\t", row.names = F)
}