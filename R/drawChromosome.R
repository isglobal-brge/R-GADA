drawChromosome<-function(chr, size=0.05, limits, print.names=TRUE, ...)
 {
  on.exit(par("xpd"=TRUE))
  par(xpd=FALSE)

  if (chr==23)
      chr <- "X"
  if (chr == 24)
      chr <- "Y"


  chromName<-paste("chr",chr,sep="")
  if (missing(limits))
    limits<-c(0,tamanysChroms$size[tamanysChroms$chrom==chromName])

  plot(0, axes = FALSE, xlab = "", type = "n", col = "gray", xlim = limits, ...)

  pp<-par("usr")
  ytop <- pp[3]
  ybottom <- pp[3]+size
  ymig <- (ytop+ybottom)/2


  ### pintem ratlles horitzontals pels cromosomes
  segments(0,ymig,tamanysChroms[tamanysChroms$chrom==chromName,2],ymig,col="gray",cex=2) 
  bandes<-ideogram[ideogram$chrom==chromName,]
  bandes$posGrafic <- bandes$chromStart+(bandes$chromEnd-bandes$chromStart)/2
  numBandes<-dim(bandes)[1]
  gieStain<-levels(bandes$gieStain) 
  colorsBandes<-c("red","gray100","black","gray25","gray50","gray75","gray50","gray50")
  japintat<-0

  for(h in 1:numBandes)
  {
   chromStart<-bandes$chromStart[h]
   chromEnd<-bandes$chromEnd[h]
   giemsaBand<-bandes$gieStain[h]
   colorBanda<-colorsBandes[giemsaBand]
   if(giemsaBand == "acen")
    {
     segments(chromStart,ybottom,chromStart,ytop,col="black")
    } 
   else 
   {
    if(giemsaBand == "gvar")
     {
      rect(chromStart,ybottom,chromEnd,ytop,col=colorBanda,angle=90,density=0.9)
      cylindrect(chromStart,ybottom,chromEnd,ytop,col=colorBanda,gradient="y",nslices=50)
     } 
     else 
     {
      rect(chromStart,ybottom,chromEnd,ytop,col=colorBanda,border="gray50")
      cylindrect(chromStart,ybottom,chromEnd,ytop,col=colorBanda,gradient="y",nslices=50)
     } 
    }
   if(japintat==0 & giemsaBand == "acen")
   {
     japintat <- 1
     points(chromEnd, ymig, pch=19, col = "gray", xpd = TRUE, cex=0.8)
   }
  }

  
# Noms de les bandes
  if (print.names)
   text(bandes$posGrafic,rep(ybottom-1.5*size,numBandes),bandes$name,srt=-90,adj=0,cex=0.6,xpd=TRUE)


 }

