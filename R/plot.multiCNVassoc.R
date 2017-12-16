


plot.multiCNVassoc<-function(x, map, ...)
 {


plotByChrom <- function(pvals.name, db, threshold="gwas", ymin=0, ymax,
                        cex=0.5, type="p", chrom.use=c(1:22, "X"),
                        main="", xlab="Chromosome", xaxis=TRUE, xaxt="n",
                        col1="green", col2="lightblue", col3="blue", col4="red", ...){

  stopifnot( all(c(pvals.name, "grand.dist", "chromosome", "badQC") %in% colnames(db)) )
  db <- db[ which( db$chromosome %in% chrom.use ), ]
  db <- db[ order(db$grand.dist), ]

  
  pvals  <- db[ , pvals.name]
  db$LOD <- -log10( pvals )

  if (threshold=="gwas")
   {
     threshold <- 0.05/length(pvals)
   }
 

  if(missing(ymax)) ymax <- max( db$LOD, na.rm=TRUE )  
  pt.cex <- cex * (db$LOD + 0.1)

  pt.cols <- ifelse( db$chromosome %in% c(2,4,6,8,10,12,14,16,18,20,22), col2, col3 )  
  pt.cols[ which(pvals < threshold) ] <- col1
  pt.cols[ which(pvals < threshold & db$badQC==1) ] <- col4

  plot( db$grand.dist, db$LOD, type=type, pch=".", xaxt=xaxt, pty="m", bty="L",
        xlab=xlab, ylab=expression( -log[10](P) ), col=pt.cols, cex=pt.cex,
        main=main, ylim=c(ymin, ymax), ...)

  midpts <- tapply( db$grand.dist, db$chromosome, mean )
  if(xaxis) axis(side=1, at=midpts, labels=names(midpts))
}


dataForPlot<-function(map)
 {
   snp.info <- NULL
   snp.info$chromosome <- as.factor(map$chr)
   snp.info$pos <- map$pos.inf
   snp.info$grand.dist <- grand.dist( snp.info$chr, snp.info$pos )
   ans<-data.frame(snp.info)
   ans
 }



grand.dist <- function(chr, coordinate, gap=10000){

  if(!is.factor(chr)) stop("chr must be a factor and levels ordered accordingly")

  n.chr <- length( levels(chr) )

if (F)
  cat("\n", n.chr, "chromosomes detected and arranged as:\n", levels(chr), "\n\n")
  
  chr.first <- tapply( coordinate, chr, min ) ## position of first snp on each chromosome
  chr.last  <- tapply( coordinate, chr, max )
  chr.range <- chr.last - chr.first + gap

  cumulative.start <- c(0, cumsum(chr.range)[-n.chr]) - chr.first
  names(cumulative.start) <- levels(chr)
  
  grand.dist <- cumulative.start[ as.character(chr) ] + coordinate
  ## ## make sure the distances don't overlap by chr
  ## tmp <- do.call(rbind, tapply(grand.dist, chr, range))
  ## tmp[-1,1] - tmp[-n.chr,2]
}


   snp.info <- dataForPlot(map)
   snp.info$badQC<-rep(FALSE,nrow(snp.info))
   snp.info$pval<-unlist(x)
   plotByChrom( "pval", snp.info, ...)


}



