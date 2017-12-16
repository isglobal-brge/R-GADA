`plotlogRatio` <-
function(x, chr, chromStart, chromEnd, B.allele.freq=FALSE, num.points, postscript=FALSE, ylim, ...)
 {

  require(plotrix)
  data(genomicInfo) #only needed for plots.


  if (is.data.frame(x) | inherits(x,"setupGADA"))
    xx<-x$log.ratio
  else
    xx<-x  

  xx.all<-xx

  gen.info<-gen.info.all<-attr(x,"gen.info")

  n.points<-length(xx)
  if (!missing(num.points))
   {
     ss<-seq(1, length(xx), length=num.points)
     xx<-xx[ss]
     gen.info<-gen.info[ss,]
     if (B.allele.freq)
      x$B.allele<-x$B.allele[ss]  

     xxx<-c(1:n.points)[ss]
   }
  else 
   {
    ss<-NULL
    xxx<-c(1:n.points)
   } 

  if (missing(chr))
   {
    if (!postscript)
     {
      X11(width = 15, height = 5)
     }

    if (missing(ylim))
     ylim<-range(xx.all,na.rm=TRUE,finite=TRUE)
 
    posOK<-o<-NULL
    xlab.ok<-ifelse(is.null(gen.info),"Probes","Chromosome")
    plot(xxx,xx,xlab=xlab.ok,ylab="log-ratio",pch=19,cex=0.2,col='darkgreen',axes=FALSE, ylim=ylim)
    axis(2,las=2)
    abline(h=0,lty=2,lwd=2,col="yellow")
  
    if (!is.null(gen.info))
     {
      tt<-cumsum(table(gen.info.all[,2]))
      tt2<-c(0,tt)
      control<-par("usr")
      for (i in 1:length(tt))
       { 
        segments(tt[i],control[3],tt[i],control[4],lty=2,col="blue")
        text(((tt2[i+1]-tt2[i])/2)+tt2[i],control[3],names(tt)[i],xpd=TRUE,cex=0.8)
       }
     } 

    }
  
   else
    {
     old.mar <- par("mar")
     old.mfrow <- par("mfrow")
     on.exit(par(mar = old.mar, mfrow = old.mfrow))

     if (missing(chromStart) & missing(chromEnd))
      {

       o<-gen.info[,2]==chr
       if (any(is.na(o)))
         warning("There are probes with chromosome equal to NA") 
       o[is.na(o)]<-FALSE      

       if (B.allele.freq)
        { 
          m <- matrix(c(1, 2), nrow = 2, ncol = 1, byrow = TRUE)
          layout(m, heights = c(0.3, 0.7))
          par(mar = c(0, old.mar[2], old.mar[3], old.mar[4]))
          plot(x$B.allele[o],pch=19,cex=0.2,axes=FALSE,ylim=c(0,1),col="darkgreen",ylab="B Allele Freq", xlab="")
          axis(2,las=2)
          ll<-seq(0,1,0.2)
          segments(rep(0,length(ll)),ll,rep(length(x$B.allele[o]),length(ll)),ll,lty=2,col="gray80")
         title(paste("Chromosome",chr))
        }

       if (B.allele.freq)
         par(mar = c(old.mar[1], old.mar[2], 1, old.mar[4]))

       posOK<-gen.info[o,3]
       logratioOK<-xx[o]

       if (missing(ylim))
        rr<-range(logratioOK, na.rm=TRUE, finite=TRUE)
       else
        rr<-ylim

       adjusted.size<-(rr[2]-rr[1])/50
       drawChromosome(chr,size=adjusted.size,ylim=rr,ylab="")
       points(posOK,logratioOK,pch=10,cex=0.2,col='darkgreen') 
       axis(2,las=2)
       segments(0,0,posOK[length(posOK)],0,lty=2,lwd=2,col="yellow")   
       if (!B.allele.freq)
         title(paste("Chromosome",chr))
      }
     else
      {

       o<-gen.info[,2]==chr & (gen.info[,3]>=chromStart & gen.info[,3]<=chromEnd)
       if (any(is.na(o)))
         warning("There are probes with chromosome equal to NA")
       o[is.na(o)]<-FALSE  

       if (B.allele.freq)
        { 
         m <- matrix(c(1, 2), nrow = 2, ncol = 1, byrow = TRUE)
         layout(m, heights = c(0.3, 0.7))
         par(mar = c(0, old.mar[2], old.mar[3], old.mar[4]))
         plot(x$B.allele[o],pch=19,cex=0.2,axes=FALSE,ylim=c(0,1),col="darkgreen",ylab="B Allele Freq", xlab="")
         axis(2,las=2)
         ll<-seq(0,1,0.2)
         segments(rep(0,length(ll)),ll,rep(length(x$B.allele[o]),length(ll)),ll,lty=2,col="gray80")
#         title(paste("Chromosome",chr))
        }

       if (B.allele.freq)
         par(mar = c(old.mar[1], old.mar[2], 1, old.mar[4]))

       posOK<-gen.info[o,3]
       logratioOK<-xx[o]
       if (missing(ylim))
        rr<-range(logratioOK, na.rm=TRUE, finite=TRUE)
       else
        rr<-ylim
       adjusted.size<-(rr[2]-rr[1])/50
       drawChromosomeZoom(chr,size=adjusted.size, chromStart=chromStart, chromEnd=chromEnd, ylim=rr,ylab="")
       points(posOK,logratioOK,pch=10,cex=0.2,col='darkgreen') 
       axis(2,las=2)
       segments(chromStart,0,chromEnd,0,lty=2,lwd=2,col="yellow")   
       if (!B.allele.freq) 
         title(paste("Chromosome",chr))

      }
    }
 invisible(list(posOK=posOK,o=o))
 }
