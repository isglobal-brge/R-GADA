`plotRatio.BackwardElimination` <-
function(x, chr, num.points, postscript=FALSE, ...)
{
 
 xx<-attr(x,"data")
 gen.info<-attr(x,"gen.info")
 attr(xx,"gen.info")<-gen.info

 ans<-plotlogRatio(xx, chr=chr, B.allele.freq=FALSE, num.points=num.points, postscript=postscript, ...)

 posOK<-ans$posOK
 
 if (!is.null(gen.info))
  {
   Segments<-summary(x, print=FALSE) 
  }
 else
  {
   Segments<-WextIextToSegments(x)
  }

 if (missing(chr))
  {
   if (is.data.frame(x) | inherits(x,"setupGADA"))
     M<-length(xx$log.ratio)
   else
     M<-length(xx)
   lines(1:M,rep(Segments[,4],Segments[,3]),col='red')
  }

 else
  {
    Segments.selected<-Segments[Segments[,5]==chr,]
    posOK<-NULL
    for (i in 1:nrow(Segments.selected))
      {
       temp<-seq(Segments.selected[i,1],Segments.selected[i,2],length=Segments.selected[i,3])
       posOK<-c(posOK,temp)
      }

    lines(posOK,rep(Segments.selected[,4],Segments.selected[,3]),col='red')
  }

}


