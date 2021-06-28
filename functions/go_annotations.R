### Selection of genes related to GO-terms based on their evidence 
### Laurent Winckers, Maastricht University - Department of Bioinformatics BiGCaT
### 2021-06-28

go_annotations <- function(go.term, biomart, output) {
  
  res <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id', 'go_linkage_type'), 
                        filters = c("go_parent_term","chromosome_name"), 
                        values = list(go_parent_term=go.term, chromosome_name=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','MT','X','Y')),
                        mart = biomart)
  res <- res[!is.na(res$entrezgene_id),]
  res <- res[res$go_linkage_type != 'ND',]
  res <- res[res$go_linkage_type != 'NAS',]
  res <- res[res$go_linkage_type != 'IEA',]
  res <- res[res$go_linkage_type != 'ISS',]
  res <- res[res$go_linkage_type != 'ISO',]
  res <- res[res$go_linkage_type != 'ISA',]
  res <- res[res$go_linkage_type != 'ISM',]
  res <- res[res$go_linkage_type != 'IGC',]
  res <- res[res$go_linkage_type != 'RCA',]
  res <- res[,-4]
  res <- res[!duplicated(res),]
  
  print(paste(length(res$entrezgene_id)," genes found for ", go.term, sep=""))
  write.table(res, file = output, quote = F, sep = "\t")
  return(res)
}
