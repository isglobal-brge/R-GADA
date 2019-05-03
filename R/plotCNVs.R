plotCNVs <- function(x, range, genome="hg19", drawGenes=FALSE,
                     col.cnvs =c("darkgreen", "darkblue"),
                     group, mosaic = FALSE){
  
  if(missing(range))
    stop("a GenomicRange should be passed from 'range' argument")
  
  chr <- as.character(seqnames(range))
  cnvs.range <- subsetByOverlaps(x, range)
  
  if (missing(group)) {
    if (!mosaic) {
     fill <- ifelse(cnvs.range$State==1, col.cnvs[1], col.cnvs[2]) 
    }
    else{
     fill <- ifelse(cnvs.range$State==1, "orange",
                    ifelse(cnvs.range$State==2, "darkgreen",
                           ifelse(cnvs.range$State==3, "darkblue",
                                  ifelse(cnvs.range$State==4, "tomato",
                                         ifelse(cnvs.range$State==5, "red", "white"))))) 
    }
    
  }
  else {
    gg <- mcols(cnvs.range)[, group]
    gg.l <- levels(as.factor(gg))
    if(length(gg.l) > 2)
      stop("grouping variable must be two levels \n
           send an email to the maintainer for a improved version")
    fill <- ifelse(gg==gg[1], "tomato", "lightblue")
  }
  cnvs.l <- AnnotationTrack(cnvs.range, 
                            fill = fill,
                            name = "Individuals",
                            group = cnvs.range$sample, 
                            cex.group=0.5)
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = genome, 
                          chromosome = chr)
  if (drawGenes) {
    if (genome=="hg19" & requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene")) {
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    }
    else if (genome=="hg18" & requireNamespace("TxDb.Hsapiens.UCSC.hg18.knownGene"))
      txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene
    else {
     warning("Genes are not shown since TxDb database is not installed in you computer")
       drawGenes <- FALSE
    }  
  } 
  
  if(drawGenes) {
    
    allg <- genes(txdb)
    allg.range <- subsetByOverlaps(allg, range)
    allg.range$symbol <- mapIds(Homo.sapiens::Homo.sapiens, 
                                keys=allg.range$gene_id,
                                keytype="ENTREZID",
                                column="SYMBOL")
    
    
    grtrack <- GeneRegionTrack(allg.range, genome = genome,
                               chromosome = chr, showId=TRUE,
                               geneSymbol=TRUE,
                               start = start(range), 
                               end = end(range),
                               name = "Genes")
    
    plotTracks(c(itrack, gtrack, grtrack, cnvs.l), 
               groupAnnotation = "group", from = start(rr),
               to = end(rr))
  }
  else{
    plotTracks(c(itrack, gtrack, cnvs.l), 
               groupAnnotation = "group", from = start(rr),
               to = end(rr))
  }
}
